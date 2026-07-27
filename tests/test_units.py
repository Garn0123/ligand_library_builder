"""#5 — pure functions in db2common.

Fast, in-process checks of the primitives every stage relies on: id extraction,
token-bounded relabeling, weight metrics, inode dedup, gzip content-sniffing,
and tsv validation.
"""

import os

import pytest

import db2common as C
import db2gen


# --- id extraction -------------------------------------------------------

def test_extract_id_text_and_bytes_agree():
    rec = db2gen.make_record("ZINC00000042", n_conformers=3)
    lines = rec.splitlines(keepends=True)
    assert C.extract_id(lines) == "ZINC00000042"
    assert C.extract_id_bytes(rec.encode()) == "ZINC00000042"


def test_extract_id_no_zinc_token():
    lines = ["M no zinc here\n", "A 0\n", "E\n"]
    assert C.extract_id(lines) == "NO_ID"
    assert C.extract_id_bytes("".join(lines).encode()) == "NO_ID"


def test_extract_id_scans_continuation_m_lines():
    lines = ["M header only\n", "M ZINC00000007 more\n", "A 0\n", "E\n"]
    assert C.extract_id(lines) == "ZINC00000007"
    assert C.extract_id_bytes("".join(lines).encode()) == "ZINC00000007"


def test_extract_id_stops_at_first_non_header_line():
    # A ZINC-looking token on a non-M line must NOT be picked up as the id.
    lines = ["M header only\n", "A ZINC99999999 not-the-id\n", "E\n"]
    assert C.extract_id(lines) == "NO_ID"
    assert C.extract_id_bytes("".join(lines).encode()) == "NO_ID"


# --- relabel -------------------------------------------------------------

def test_relabel_only_rewrites_m_lines():
    lines = ["M ZINC1 foo\n", "A ZINC1 keep-this\n", "E\n"]
    out = C.relabel(lines, "ZINC1", "ZINC1_1")
    assert out[0] == "M ZINC1_1 foo\n"       # header rewritten
    assert out[1] == "A ZINC1 keep-this\n"   # data line untouched
    assert out[2] == "E\n"


def test_relabel_is_token_bounded():
    # ZINC0000003 must not match inside the longer token ZINC00000031.
    lines = ["M ZINC00000031 x\n", "E\n"]
    out = C.relabel(lines, "ZINC0000003", "REPLACED")
    assert out[0] == "M ZINC00000031 x\n"


def test_relabel_column_shift_is_expected():
    # The documented caveat: a longer id shifts every column to its right.
    lines = ["M ZINC1 a b c\n"]
    out = C.relabel(lines, "ZINC1", "ZINC1_10")
    assert out[0] == "M ZINC1_10 a b c\n"


# --- weight metrics ------------------------------------------------------

def test_weight_count_text_and_bytes():
    rec = db2gen.make_record("ZINC1", n_conformers=5)
    lines = rec.splitlines(keepends=True)
    assert C.make_weight_fn("count")(lines) == 1
    assert C.make_weight_fn_bytes("count")(rec.encode()) == 1


def test_weight_lines_c_text_and_bytes_agree():
    rec = db2gen.make_record("ZINC1", n_conformers=5)
    lines = rec.splitlines(keepends=True)
    assert C.make_weight_fn("lines:C")(lines) == 5
    assert C.make_weight_fn_bytes("lines:C")(rec.encode()) == 5


def test_weight_bytes_metric_agree():
    rec = db2gen.make_record("ZINC1", n_conformers=2)
    lines = rec.splitlines(keepends=True)
    assert C.make_weight_fn("bytes")(lines) == len(rec.encode())
    assert C.make_weight_fn_bytes("bytes")(rec.encode()) == len(rec.encode())


def test_unknown_weight_spec_errors():
    with pytest.raises(SystemExit):
        C.make_weight_fn("nonsense")
    with pytest.raises(SystemExit):
        C.make_weight_fn_bytes("nonsense")


# --- filesystem helpers --------------------------------------------------

def test_dedupe_by_inode_catches_hardlink(tmp_path):
    a = str(tmp_path / "a.db2.gz")
    db2gen.write_gz(a, db2gen.make_record("ZINC1"))
    b = str(tmp_path / "b.db2.gz")
    os.link(a, b)                         # hardlink: same (dev, inode)
    kept, skipped = C.dedupe_by_inode([a, b])
    assert kept == [a]
    assert skipped == [(b, a)]


def test_is_gzip_sniffs_content_not_extension(tmp_path):
    real = str(tmp_path / "real.db2.gz")
    db2gen.write_gz(real, db2gen.make_record("ZINC1"))
    fake = str(tmp_path / "fake.db2.gz")   # plaintext despite the .gz name
    db2gen.write_plain(fake, db2gen.make_record("ZINC1"))
    assert C.is_gzip(real) is True
    assert C.is_gzip(fake) is False


# --- tsv validation ------------------------------------------------------

def test_read_tsv_rejects_wrong_header(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("wrong\theader\n")
    with pytest.raises(SystemExit):
        list(C.read_tsv(str(p), C.MANIFEST_HEADER))


def test_read_tsv_rejects_wrong_column_count(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("\t".join(C.MANIFEST_HEADER) + "\n" + "only\ttwo\n")
    with pytest.raises(SystemExit):
        list(C.read_tsv(str(p), C.MANIFEST_HEADER))


def test_read_tsv_accepts_valid_rows(tmp_path):
    p = tmp_path / "m.tsv"
    rows = [C.MANIFEST_HEADER,
            ["chunk_00001.db2.gz", "0", "src.db2.gz", "0", "ZINC1"],
            ["chunk_00001.db2.gz", "1", "src.db2.gz", "1", "ZINC2"]]
    p.write_text("\n".join("\t".join(r) for r in rows) + "\n")
    got = list(C.read_tsv(str(p), C.MANIFEST_HEADER))
    assert got == rows[1:]
