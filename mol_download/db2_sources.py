"""Stream .db2 records out of ZINC22 archives without extracting them.

Drop-in for a pipeline loop that currently does::

    for fi, path in enumerate(files):
        opener = gzip.open if is_gzip(path) else open
        with opener(path, "rb") as fh:
            for src_idx, (rec, complete) in enumerate(iter_records_bytes(fh)):
                ...

becomes::

    for fi, (path, fh) in enumerate(iter_sources(files, on_error)):
        for src_idx, (rec, complete) in enumerate(iter_records_bytes(fh)):
            ...

DO NOT extract the tarballs first. Each .db2.tgz holds thousands of tiny
members; across a few thousand archives that is millions of small files on a
shared parallel filesystem. On GPFS you will exhaust the inode allocation long
before the byte quota and the metadata load is antisocial. The archive is also
already the right unit of parallelism.
"""

import gzip
import os
import tarfile
import zlib

TAR_EXTS = (".tgz", ".tar.gz", ".tar")
DB2_EXTS = (".db2", ".db2.gz")

GZIP_MAGIC = b"\x1f\x8b"


def is_gzip(path):
    """True if the file begins with the gzip magic number."""
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == GZIP_MAGIC
    except OSError:
        return False


def iter_sources(paths, on_error=None):
    """Yield ``(label, binary_fh)`` for every .db2 stream in ``paths``.

    ``paths`` may mix tar archives and bare .db2/.db2.gz files. For archive
    members the label is ``"<archive>::<member>"``; for bare files it is just
    the path.

    ``on_error(path, message)`` is called for anything unreadable and iteration
    continues. A corrupt gzip container costs you one archive, not the run.

    IMPORTANT: archives are opened in streaming mode ("r|gz"), so each yielded
    handle is only valid until the generator advances. Consume it fully before
    requesting the next one. Do not collect handles into a list and do not hand
    them to a thread pool.
    """
    if on_error is None:
        def on_error(path, message):
            pass

    for path in paths:
        if path.endswith(TAR_EXTS):
            mode = "r|" if path.endswith(".tar") else "r|gz"
            try:
                with tarfile.open(path, mode) as tf:
                    for member in tf:
                        if not member.isfile():
                            continue
                        if not member.name.endswith(DB2_EXTS):
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
                # zlib.error is NOT an OSError. A .tgz with a valid gzip header
                # but a damaged deflate body raises it from deep inside the
                # read, and would otherwise kill the whole task.
                on_error(path, "archive unreadable: {}".format(exc))
        else:
            opener = gzip.open if is_gzip(path) else open
            try:
                with opener(path, "rb") as fh:
                    yield path, fh
            except (OSError, EOFError, zlib.error) as exc:
                on_error(path, str(exc))


def split_label(label):
    """Split a manifest label back into ``(archive, member)``.

    Bare-file labels return ``(path, None)``.
    """
    if "::" in label:
        archive, member = label.split("::", 1)
        return archive, member
    return label, None


def find_archives(root, pattern="*.tgz"):
    """Sorted list of archives under ``root``. Stable order across shards."""
    import fnmatch
    out = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for name in filenames:
            if fnmatch.fnmatch(name, pattern):
                out.append(os.path.join(dirpath, name))
    out.sort()
    return out


def shard(items, index, total):
    """Deterministic stride-shard. Keeps each archive whole in one shard --
    splitting members of one archive across shards would break per-shard bin
    balancing state (writers / counts / heap)."""
    if total <= 1:
        return list(items)
    return [x for i, x in enumerate(items) if i % total == index]


if __name__ == "__main__":
    import sys

    errs = []
    n = 0
    for label, fh in iter_sources(sys.argv[1:], lambda p, e: errs.append((p, e))):
        n += 1
        fh.read(1)  # prove the handle is live
    print("{} db2 members, {} unreadable archives".format(n, len(errs)))
    for path, msg in errs:
        print("  ERR {}: {}".format(path, msg))
