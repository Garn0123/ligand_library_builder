"""ZINC22 support: read .db2 records from tarballs, new (truncated) id scheme.

New ZINC22 data ships as `.db2.tgz` archives of thousands of small `.db2`
members; we read them in place (never extracting). The molecule id lives in the
member filename (`ZINC550000002ttm.0.O.db2`); the record header carries only a
left-truncated form (`NC550000002ttm.0`) because the db2 id field is fixed-width.
"""

import os
from collections import Counter

import db2common as C
import db2gen


# --- id helpers ----------------------------------------------------------

def test_id_from_name_full_zinc_with_conformer():
    assert C.id_from_name(
        "arc.db2.tgz::H05/H05M000/2t/tm/ZINC550000002ttm.0.O.db2"
    ) == "ZINC550000002ttm.0"
    assert C.id_from_name("/data/H17/H17P200/tranche.db2.gz") == "NO_ID"


def test_extract_id_new_truncated_header():
    rec = db2gen.make_zinc22_record("NC550000002ttm.0")
    assert C.extract_id_bytes(rec.encode()) == "NC550000002ttm.0"
    assert C.extract_id(rec.splitlines(keepends=True)) == "NC550000002ttm.0"


def test_extract_id_old_still_wins():
    rec = db2gen.make_record("ZINC00000042")
    assert C.extract_id_bytes(rec.encode()) == "ZINC00000042"


def test_id_match_tolerates_truncation():
    assert C.id_match("NC550000002ttm.0", "ZINC550000002ttm.0")   # suffix
    assert C.id_match("ZINC00000042", "ZINC00000042")             # exact
    assert not C.id_match("NC999999999x.0", "ZINC550000002ttm.0")  # different
    assert not C.id_match("NO_ID", "ZINC550000002ttm.0")


# --- iter_sources --------------------------------------------------------

def _write_archive(path, specs):
    """specs: list of (member_name, full_zinc_id). Header id is truncated."""
    members = {name: db2gen.make_zinc22_record(db2gen.truncate_zinc_id(zid))
               for name, zid in specs}
    db2gen.write_db2_tgz(path, members)


def test_iter_sources_reads_tar_members(tmp_path):
    arc = str(tmp_path / "H05M000-O-aaaaaa.db2.tgz")
    _write_archive(arc, [
        ("H05/H05M000/2t/tm/ZINC550000002ttm.0.O.db2", "ZINC550000002ttm.0"),
        ("H05/H05M000/2r/Ao/ZINC550000002rAo.2.O.db2", "ZINC550000002rAo.2"),
        ("H05/H05M000/2D/z7/ZINC550000002Dz7.0.O.db2.gz", "ZINC550000002Dz7.0"),
    ])
    errs = []
    ids = []
    for label, fh in C.iter_sources([arc], lambda p, e: errs.append((p, e))):
        assert label.startswith(arc + "::")
        recs = list(C.iter_records_bytes(fh))
        assert len(recs) == 1 and recs[0][1] is True
        ids.append(C.id_from_name(label))
    assert errs == []
    assert set(ids) == {"ZINC550000002ttm.0", "ZINC550000002rAo.2",
                        "ZINC550000002Dz7.0"}


def test_iter_sources_reports_corrupt_archive(tmp_path):
    arc = str(tmp_path / "bad.db2.tgz")
    with open(arc, "wb") as fh:
        fh.write(b"\x1f\x8b\x08" + b"\xff" * 200)   # gzip magic, garbage body
    errs = []
    out = list(C.iter_sources([arc], lambda p, e: errs.append((p, e))))
    assert out == []                 # nothing yielded
    assert len(errs) == 1 and errs[0][0] == arc


def test_iter_sources_bare_files_still_work(tmp_path):
    p = str(tmp_path / "old.db2.gz")
    db2gen.write_gz(p, db2gen.records_text([("ZINC00000001", 2)]))
    got = [(label, list(C.iter_records_bytes(fh)))
           for label, fh in C.iter_sources([p], lambda *_: None)]
    assert len(got) == 1
    assert got[0][0] == p            # bare-file label is the path, no "::"


# --- p1_collect end-to-end over a tarball --------------------------------

def test_parallel_pipeline_over_tarball(tmp_path, run_script):
    src = tmp_path / "src" / "zinc22"
    specs = [
        ("H05/H05M000/2t/tm/ZINC550000002ttm.0.O.db2", "ZINC550000002ttm.0"),
        ("H05/H05M000/2r/Ao/ZINC550000002rAo.2.O.db2", "ZINC550000002rAo.2"),
        ("H05/H05M000/2r/Ao/ZINC550000002rAo.0.O.db2", "ZINC550000002rAo.0"),
        ("H05/H05M000/2D/z7/ZINC550000002Dz7.0.O.db2.gz", "ZINC550000002Dz7.0"),
    ]
    _write_archive(str(src / "H05M000-O-aaaaaa.db2.tgz"), specs)
    # a second archive so sharding has >1 source
    _write_archive(str(src / "H05M000-O-aaaaab.db2.tgz"), [
        ("H05/H05M000/2n/nX/ZINC550000002nnX.0.O.db2", "ZINC550000002nnX.0"),
        ("H05/H05M000/2o/mt/ZINC550000002omt.0.O.db2", "ZINC550000002omt.0"),
    ])
    expected_ids = {zid for _, zid in specs} | {"ZINC550000002nnX.0",
                                                "ZINC550000002omt.0"}

    work = tmp_path / "work"
    out = tmp_path / "out"
    run_script("make_shards.py", "-i", src, "-o", work, "-S", 2)
    for s in range(2):
        run_script("p1_collect.py", "-w", work, "-s", s, "-N", 3)
    for b in range(3):
        run_script("p2_assemble.py", "-w", work, "-b", b, "-o", out)
    run_script("p3_finalize.py", "-w", work, "-o", out)

    rows = list(C.read_tsv(os.path.join(str(out), "manifest.tsv"),
                           C.MANIFEST_HEADER))
    got_ids = Counter(r[4] for r in rows)
    assert set(got_ids) == expected_ids           # full ZINC ids in the manifest
    assert sum(got_ids.values()) == len(expected_ids)
    # provenance: source_file is archive::member for every row
    assert all("::" in r[2] and r[2].endswith(".db2") or r[2].endswith(".db2.gz")
               for r in rows)
    assert all(".db2.tgz::" in r[2] for r in rows)


def test_make_shards_finds_and_counts_archives(tmp_path, run_script):
    src = tmp_path / "src"
    _write_archive(str(src / "a.db2.tgz"),
                   [("p/ZINC55000000000{}.0.O.db2".format(i),
                     "ZINC55000000000{}.0".format(i)) for i in range(6)])
    work = tmp_path / "work"
    res = run_script("make_shards.py", "-i", src, "-o", work, "-S", 1,
                     "--target-per-bin", 3)
    assert "found 1 files" in res.stderr           # the archive itself
    bins = int((work / "bins.txt").read_text().strip())
    assert bins == 2                               # 6 members / 3 per bin


# --- check_inputs on archives --------------------------------------------

def test_check_inputs_classifies_archives(tmp_path, run_script):
    src = tmp_path / "src"
    _write_archive(str(src / "ok.db2.tgz"),
                   [("p/ZINC550000000001.0.O.db2", "ZINC550000000001.0")])
    with open(str(src / "bad.db2.tgz"), "wb") as fh:
        fh.write(b"\x1f\x8b\x08" + b"\xff" * 200)

    good = tmp_path / "good.txt"
    bad = tmp_path / "bad.txt"
    res = run_script("check_inputs.py", "-i", src, "--good-list", good,
                     "--bad-list", bad, "--progress-interval", 0, check=False)
    assert res.returncode != 0
    assert {os.path.basename(p) for p in good.read_text().split()} == {"ok.db2.tgz"}
    assert {os.path.basename(p) for p in bad.read_text().split()} == {"bad.db2.tgz"}


# --- 03_apply relabel with truncated ids ---------------------------------

def test_03_relabels_truncated_ids(tmp_path, run_script):
    # Two records with the same full id (e.g. two state codes of one conformer)
    # -> duplicates -> 02 would suffix them. The chunk header carries the
    # truncated id; 03 must verify and relabel it despite the truncation.
    chunks = tmp_path / "chunks"
    os.makedirs(str(chunks))
    db2gen.write_gz(
        os.path.join(str(chunks), "chunk_00001.db2.gz"),
        db2gen.make_zinc22_record("NC550000002ttm.0")
        + db2gen.make_zinc22_record("NC550000002ttm.0"))

    labels = tmp_path / "labels.tsv"
    with open(str(labels), "w") as fh:
        fh.write("\t".join(C.LABELS_HEADER) + "\n")
        fh.write("chunk_00001.db2.gz\t0\tZINC550000002ttm.0\tZINC550000002ttm.0_1\n")
        fh.write("chunk_00001.db2.gz\t1\tZINC550000002ttm.0\tZINC550000002ttm.0_2\n")

    out = tmp_path / "labelled"
    run_script("03_apply.py", "-c", chunks, "-L", labels, "-o", out,
               "--copy-unedited")

    ids = db2gen.chunk_ids(os.path.join(str(out), "chunk_00001.db2.gz"))
    # the truncated header id gets the same suffix appended
    assert ids == ["NC550000002ttm.0_1", "NC550000002ttm.0_2"]


def test_03_still_rejects_genuinely_stale_labels(tmp_path, run_script):
    chunks = tmp_path / "chunks"
    os.makedirs(str(chunks))
    db2gen.write_gz(os.path.join(str(chunks), "chunk_00001.db2.gz"),
                    db2gen.make_zinc22_record("NC550000002ttm.0"))
    labels = tmp_path / "labels.tsv"
    with open(str(labels), "w") as fh:
        fh.write("\t".join(C.LABELS_HEADER) + "\n")
        fh.write("chunk_00001.db2.gz\t0\tZINC999999999xyz.9\tZINC999999999xyz.9_1\n")

    out = tmp_path / "labelled"
    res = run_script("03_apply.py", "-c", chunks, "-L", labels, "-o", out,
                     check=False)
    assert res.returncode != 0
    assert not os.path.exists(os.path.join(str(out), "chunk_00001.db2.gz"))
