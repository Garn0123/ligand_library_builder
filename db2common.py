"""Shared helpers for the db2 chunk/label/apply pipeline.

db2 files are record-oriented: each molecule starts with an 'M' header line
and ends with an 'E' terminator. Nothing here ever splits a record.

Inputs may be old-style bare `.db2.gz` tranche files or new ZINC22 `.db2.tgz`
archives, each holding thousands of `.db2` members. `iter_sources` reads either
without ever extracting an archive to disk.
"""

import gzip
import os
import re
import sys
import tarfile
import time
import zlib

# Archive vs member extensions, and what find_inputs discovers by default for
# the tar-aware (parallel) tools: old bare gzip tranches + new .db2.tgz archives.
TAR_EXTS = (".tgz", ".tar.gz", ".tar")
DB2_EXTS = (".db2", ".db2.gz")
INPUT_SUFFIXES = (".db2.gz", ".tgz", ".tar.gz", ".tar")


def find_inputs(root, suffixes=(".db2.gz",)):
    """Walk root, return sorted list of files matching any suffix."""
    hits = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            if any(fn.endswith(s) for s in suffixes):
                hits.append(os.path.join(dirpath, fn))
    hits.sort()
    return hits


def dedupe_by_inode(paths):
    """Drop paths resolving to a file already in the list.

    Catches symlinks (os.stat follows them) and hardlinks (same inode,
    different name). Distinct files with identical *content* are not caught.
    Returns (kept, skipped) where skipped is [(dup_path, original_path), ...].
    """
    seen = {}
    kept, skipped = [], []
    for p in paths:
        try:
            st = os.stat(p)
        except OSError:
            kept.append(p)  # let the reader surface the real error
            continue
        key = (st.st_dev, st.st_ino)
        if key in seen:
            skipped.append((p, seen[key]))
        else:
            seen[key] = p
            kept.append(p)
    return kept, skipped


def iter_sources(paths, on_error):
    """Yield (label, binary_fh) for each .db2 stream in `paths`.

    Sources may be ZINC22 `.db2.tgz` archives (read member-by-member in place,
    never extracted -- extracting millions of tiny files is a metadata disaster
    on GPFS/Lustre) or bare `.db2` / `.db2.gz` files. For a tar member the label
    is '<archive>::<member>'; for a bare file it is the path.

    Two error scopes: a container-level failure (a corrupt archive, an
    unopenable file) is reported via on_error(path, msg) and skips that whole
    source; a malformed *record* inside a member is the caller's to catch while
    reading.

    Streaming (`r|gz`): a handle from extractfile() is valid only until the
    generator advances to the next member, so the caller must fully consume it
    before requesting the next. Do not buffer handles or hand them to a pool.
    """
    for path in paths:
        if path.endswith(TAR_EXTS):
            mode = "r|" if path.endswith(".tar") else "r|gz"
            try:
                with tarfile.open(path, mode) as tf:
                    for member in tf:
                        if not member.isfile() or not member.name.endswith(DB2_EXTS):
                            continue
                        fh = tf.extractfile(member)
                        if fh is None:
                            continue
                        label = "{}::{}".format(path, member.name)
                        if member.name.endswith(".gz"):
                            with gzip.GzipFile(fileobj=fh) as gz:
                                yield label, gz
                        else:
                            yield label, fh
            except (tarfile.TarError, OSError, EOFError, zlib.error) as exc:
                on_error(path, "archive unreadable: {}".format(exc))
        else:
            opener = gzip.open if is_gzip(path) else open
            try:
                with opener(path, "rb") as fh:
                    yield path, fh
            except (OSError, EOFError, zlib.error) as exc:
                on_error(path, str(exc))


_ZINC_NAME_RE = re.compile(r"ZINC[A-Za-z0-9]+(?:\.\d+)?")


def id_from_name(name):
    """Full ZINC id (with conformer index) parsed from a filename or label.

    e.g. '<archive>::.../ZINC550000002ttm.0.O.db2' -> 'ZINC550000002ttm.0'.
    Returns 'NO_ID' when no ZINC token is present (e.g. an old tranche path,
    whose per-molecule id lives in the record header instead).
    """
    m = _ZINC_NAME_RE.search(name)
    return m.group(0) if m else "NO_ID"


def _looks_like_id(tok):
    return len(tok) >= 6 and any(ch.isdigit() for ch in tok)


def _looks_like_id_bytes(tok):
    return len(tok) >= 6 and any(48 <= b <= 57 for b in tok)


def id_match(found, expected):
    """True if the id read back from a record matches the expected id.

    Tolerates the ZINC22 truncation: the db2 id field is fixed-width, so a long
    id is left-truncated in the header (record 'NC550000002ttm.0' for manifest
    id 'ZINC550000002ttm.0'). The record id is therefore a suffix of the full
    one.
    """
    if found == expected:
        return True
    return found != "NO_ID" and expected.endswith(found)


def iter_records(fh):
    """Yield (lines, complete) for each record in an open text stream.

    complete=False on a trailing fragment with no 'E' terminator, which
    means the file is truncated.
    """
    buf = []
    for line in fh:
        buf.append(line)
        if line[:1] == "E":
            yield buf, True
            buf = []
    if buf and any(ln.strip() for ln in buf):
        yield buf, False


def iter_records_bytes(fh, bufsize=1 << 22):
    """Fast record splitter operating on raw bytes.

    Yields (record_bytes, complete). Avoids per-line Python iteration, which
    dominates runtime on db2 records with many conformer lines. The separator
    b"\nE\n" guarantees the E is at the start of a line.
    """
    sep = b"\nE\n"
    buf = b""
    first = True
    while True:
        block = fh.read(bufsize)
        if not block:
            break
        buf = buf + block
        if first and buf.startswith(b"E\n"):
            yield b"E\n", True          # degenerate leading terminator
            buf = buf[2:]
        first = False
        start = 0
        while True:
            k = buf.find(sep, start)
            if k < 0:
                break
            yield buf[start:k + len(sep)], True
            start = k + len(sep)
        if start:
            buf = buf[start:]
    if buf.strip():
        yield buf, False


def extract_id_bytes(rec, scan=512):
    """ZINC id from the head of a raw record.

    Old format: a 'ZINC...'-prefixed token on the M header lines. New (ZINC22)
    format: the id is the first token of the first M line and may be truncated
    (e.g. 'NC550000002ttm.0'); accepted only if it looks like an id.
    """
    head_block = rec[:scan]
    for line in head_block.split(b"\n"):
        if not line.startswith(b"M"):
            break
        for tok in line.split():
            if tok.startswith(b"ZINC"):
                return tok.decode("ascii", "replace")
    first = head_block.split(b"\n", 1)[0]
    if first.startswith(b"M"):
        parts = first.split()
        if len(parts) >= 2 and _looks_like_id_bytes(parts[1]):
            return parts[1].decode("ascii", "replace")
    return "NO_ID"


def make_weight_fn_bytes(spec):
    """Weight function for raw-bytes records. See make_weight_fn."""
    if spec == "count":
        return lambda rec: 1
    if spec == "bytes":
        return len
    if spec.startswith("lines:"):
        tag = spec[len("lines:"):].encode()
        if not tag:
            raise SystemExit("--weight lines: needs a line type, e.g. lines:C")
        nl_tag = b"\n" + tag
        return lambda rec: rec.count(nl_tag) + (1 if rec.startswith(tag) else 0)
    raise SystemExit("unknown --weight: {}".format(spec))


def extract_id(lines):
    """Pull the ZINC identifier from a record's M (header) lines.

    Old format: a 'ZINC...'-prefixed token. New (ZINC22) format: the first token
    of the first M line, possibly truncated (e.g. 'NC550000002ttm.0'), accepted
    only if it looks like an id.
    """
    for line in lines:
        if line[:1] != "M":
            break
        for tok in line.split():
            if tok.startswith("ZINC"):
                return tok
    head = lines[0] if lines else ""
    if head[:1] == "M":
        parts = head.split()
        if len(parts) >= 2 and _looks_like_id(parts[1]):
            return parts[1]
    return "NO_ID"


def relabel(lines, orig_id, new_id):
    """Replace orig_id with new_id on the record's M lines only.

    Token-bounded, so ZINC0000003 does not match inside ZINC00000031.
    Whitespace is preserved exactly, but columns to the RIGHT of the id
    shift if the id changes length.
    """
    pat = re.compile(r"(?<![A-Za-z0-9_.])" + re.escape(orig_id) + r"(?![A-Za-z0-9_.])")
    return [pat.sub(new_id, ln) if ln[:1] == "M" else ln for ln in lines]


def is_gzip(path):
    """True if the file starts with the gzip magic bytes 1f 8b.

    The .gz extension is not evidence -- a plaintext file can carry it, and
    vim's gzip.vim plugin hides the difference by decompressing on read.
    """
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def open_gz_text(path, mode="rt", compresslevel=6):
    """Open for reading whether or not the file is actually gzipped.

    Writing always produces real gzip.
    """
    if "r" in mode:
        if is_gzip(path):
            return gzip.open(path, mode, errors="replace")
        return open(path, mode, errors="replace")
    return gzip.open(path, mode, compresslevel=compresslevel)


def read_tsv(path, expect_header):
    """Stream a tsv, yielding field lists. Validates the header."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header != expect_header:
            raise SystemExit(
                "{}: unexpected header\n  got:      {}\n  expected: {}".format(
                    path, header, expect_header))
        for lineno, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != len(expect_header):
                raise SystemExit("{}:{}: expected {} columns, got {}".format(
                    path, lineno, len(expect_header), len(fields)))
            yield fields


class Progress:
    """Time-based progress reporting.

    Emits at most one line every `interval` seconds, so the cadence is
    predictable in wall-clock terms regardless of how the work is
    distributed across files. interval <= 0 disables output entirely.
    """

    def __init__(self, interval=30.0, total=None, unit="files", total_bytes=0):
        self.interval = interval
        self.total = total
        self.unit = unit
        self.total_bytes = total_bytes
        self.start = time.time()
        self.last = self.start

    @staticmethod
    def _size(nbytes):
        for unit in ("B", "KB", "MB", "GB", "TB"):
            if nbytes < 1024 or unit == "TB":
                return "{:.1f} {}".format(nbytes, unit)
            nbytes /= 1024.0

    @staticmethod
    def _hms(seconds):
        seconds = int(seconds)
        return "{:d}:{:02d}:{:02d}".format(
            seconds // 3600, (seconds % 3600) // 60, seconds % 60)

    def tick(self, done, records=None, force=False, done_bytes=None):
        """Report progress. `done_bytes` gives a far better ETA than file
        count when file sizes are uneven, which they are in a tranche tree."""
        if self.interval <= 0:
            return   # 0 means fully silent, including the final summary
        now = time.time()
        if not force and now - self.last < self.interval:
            return
        self.last = now
        elapsed = now - self.start
        parts = ["  " + self._hms(elapsed)]

        if self.total_bytes and done_bytes is not None:
            frac = float(done_bytes) / self.total_bytes
            parts.append("{:.1f}% of {}".format(
                100.0 * frac, self._size(self.total_bytes)))
        elif self.total:
            frac = float(done) / self.total
            parts.append("{}/{} {} ({:.1f}%)".format(
                done, self.total, self.unit, 100.0 * frac))
        else:
            frac = None
            parts.append("{} {}".format(done, self.unit))

        if self.total and self.total_bytes:
            parts.append("{}/{} {}".format(done, self.total, self.unit))
        if records is not None:
            rate = records / elapsed if elapsed > 0 else 0.0
            parts.append("{} molecules".format(records))
            parts.append("{:,.0f} mol/s".format(rate))
        if frac and frac > 0:
            parts.append("ETA " + self._hms(elapsed / frac - elapsed))

        sys.stderr.write("  |  ".join(parts) + "\n")
        sys.stderr.flush()


def make_weight_fn(spec):
    """Build a per-record cost function.

    count     every molecule costs 1 (uniform molecule counts)
    bytes     uncompressed size of the record
    lines:X   number of lines starting with X -- e.g. lines:C weights by
              conformer count, usually the best proxy for docking time.
              Verify which letter your db2 uses before relying on it.
    """
    if spec == "count":
        return lambda lines: 1
    if spec == "bytes":
        return lambda lines: sum(len(ln) for ln in lines)
    if spec.startswith("lines:"):
        tag = spec[len("lines:"):]
        if not tag:
            raise SystemExit("--weight lines: needs a line type, e.g. lines:C")
        n = len(tag)
        return lambda lines: sum(1 for ln in lines if ln[:n] == tag)
    raise SystemExit("unknown --weight: {}".format(spec))


MANIFEST_HEADER = ["chunk", "chunk_idx", "source_file", "source_idx", "original_id"]
LABELS_HEADER = ["chunk", "chunk_idx", "original_id", "new_id"]
