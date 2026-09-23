"""Synthetic db2 record + fixture-tree generators for the test suite.

A db2 record is one or more 'M' header lines, then opaque interior lines
(atoms/bonds/coords/conformers), then a bare 'E' terminator line. Only the M
lines (which carry the id) and the E terminator matter to the pipeline;
everything else is passed through byte-for-byte. The interior generated here is
therefore deliberately minimal but shaped like the real thing.

Every record produced by ``make_record`` is *well formed*: no interior line
begins with 'E', so the text parser (``iter_records``) and the bytes parser
(``iter_records_bytes``) agree on it. Tests that probe the 'E-line' divergence
build their own hand-crafted records.
"""

import gzip
import io
import os
import tarfile
from collections import Counter


def make_record(zinc_id="ZINC00000001", n_conformers=2, n_atoms=3,
                smiles="CCO", extra_m_lines=0):
    """One well-formed db2 record as a str, ending in an 'E' terminator.

    ``n_conformers`` is the number of 'C' lines, which drives the ``lines:C``
    weight metric. ``extra_m_lines`` adds continuation header lines (the id
    still lives on the first M line).
    """
    lines = ["M {} {} ref\n".format(zinc_id, smiles)]
    for i in range(extra_m_lines):
        lines.append("M cont{}\n".format(i))
    for i in range(n_atoms):
        lines.append("A {} C 0.0 0.0 {}\n".format(i, i))
    lines.append("B 0 1 1\n")
    for i in range(n_atoms):
        lines.append("X {} 0.00 0.00 {}.00\n".format(i, i))
    for i in range(n_conformers):
        lines.append("C {} 1.0\n".format(i))
    lines.append("E\n")
    return "".join(lines)


def records_text(specs):
    """Concatenate records for ``specs`` = iterable of (zinc_id, n_conformers)."""
    return "".join(make_record(zid, n_conformers=nc) for zid, nc in specs)


def write_gz(path, text):
    """Write ``text`` as a real gzip file, creating parent dirs."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write(text)


def write_plain(path, text):
    """Write ``text`` uncompressed, even if ``path`` ends in .gz.

    Used to exercise the content-sniffing (``is_gzip``) code path.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def write_corrupt_gz(path, text):
    """Write a gzip file with a valid header but a corrupted deflate body.

    Reading it raises ``zlib.error`` partway through, reproducing a damaged
    download/transfer (the 10-byte gzip header stays intact, so ``is_gzip``
    still reports True and gzip.open proceeds into the bad data). ``text`` must
    be long enough that the compressed body exceeds the header.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as fh:
        fh.write(text.encode())
    data = bytearray(buf.getvalue())
    for i in range(12, min(len(data) - 8, 48)):   # corrupt early deflate bytes
        data[i] ^= 0xFF
    with open(path, "wb") as fh:
        fh.write(bytes(data))


def build_tree(root, files):
    """Write a fixture tree of gzip'd db2 files.

    ``files`` maps a relative path -> list of (zinc_id, n_conformers). Returns a
    Counter of the id multiset written, for lossless-round-trip assertions.
    """
    ids = Counter()
    for rel, specs in files.items():
        write_gz(os.path.join(root, rel), records_text(specs))
        for zid, _nc in specs:
            ids[zid] += 1
    return ids


def make_zinc22_record(short_id="NC550000002ttm.0", n_conformers=2, n_atoms=3):
    """A ZINC22-style record whose header carries the (already left-truncated)
    id as the first token of the first M line, e.g. 'NC550000002ttm.0'.
    """
    lines = ["M {}      fake  {}  {}     23      9      3\n".format(
        short_id, n_atoms, n_atoms - 1)]
    lines.append("M   +1.0000    -43.470     -0.460    -43.930   112.700\n")
    lines.append("M fake\n")
    for i in range(n_atoms):
        lines.append("A {} C 0.0 0.0 {}\n".format(i, i))
    for i in range(n_conformers):
        lines.append("C {} 1.0\n".format(i))
    lines.append("E\n")
    return "".join(lines)


def truncate_zinc_id(full_id):
    """Left-truncate 'ZINC...'+conformer to 16 chars, mirroring the db2 writer
    (e.g. 'ZINC550000002ttm.0' -> 'NC550000002ttm.0')."""
    return full_id[-16:] if len(full_id) > 16 else full_id


def write_db2_tgz(path, members):
    """Write a `.db2.tgz`. `members` maps member-name -> record text.

    A member name ending in `.gz` is gzip-compressed inside the tar; otherwise
    it is stored as a bare `.db2`. Nothing is written to disk uncompressed.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with tarfile.open(path, "w:gz") as tf:
        for name, text in members.items():
            data = text.encode()
            if name.endswith(".gz"):
                buf = io.BytesIO()
                with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
                    gz.write(data)
                data = buf.getvalue()
            info = tarfile.TarInfo(name=name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))


def chunk_ids(path):
    """Ordered list of ids in a finished chunk, via the library's own parser.

    Imported lazily so this module has no import-time dependency on the code
    under test (conftest puts the pipeline dir on sys.path first).
    """
    import db2common
    ids = []
    with db2common.open_gz_text(path, "rt") as fh:
        for lines, complete in db2common.iter_records(fh):
            if not complete:
                break
            ids.append(db2common.extract_id(lines))
    return ids
