"""#1 — parser equivalence and record-boundary integrity.

The pipeline's two parsers must agree byte-for-byte on well-formed data: the
serial stages use the line-based ``iter_records``, the parallel collect phase
uses the ~3x-faster bytes-based ``iter_records_bytes``, and the README claims
they are "verified to produce byte-identical records." These tests enforce that
claim (the prerequisite for the planned bytes-parser refactor of 01-04) and pin
the one case where they intentionally diverge.
"""

import io

import db2common as C
import db2gen


def _text_records(data):
    return list(C.iter_records(io.StringIO(data)))


def _byte_records(data, bufsize=1 << 22):
    return list(C.iter_records_bytes(io.BytesIO(data.encode()), bufsize=bufsize))


def test_parsers_byte_identical_on_wellformed():
    specs = [("ZINC{:08d}".format(i), (i % 4) + 1) for i in range(25)]
    data = db2gen.records_text(specs)

    text = _text_records(data)
    byts = _byte_records(data)

    assert len(text) == len(byts) == 25
    for (lines, tc), (rec, bc) in zip(text, byts):
        assert tc is True and bc is True
        assert "".join(lines).encode() == rec   # byte-identical records


def test_bytes_parser_across_buffer_boundary():
    # A tiny read buffer forces records to straddle block reads; the result
    # must not depend on where the buffer happens to cut.
    specs = [("ZINC{:08d}".format(i), 3) for i in range(20)]
    data = db2gen.records_text(specs)

    tiny = _byte_records(data, bufsize=7)
    whole = _byte_records(data, bufsize=1 << 22)

    assert tiny == whole
    assert len(tiny) == 20
    assert all(complete for _rec, complete in tiny)


def test_truncated_tail_is_flagged_incomplete():
    data = db2gen.make_record("ZINC00000001") + "M ZINC00000002 x\nA 0 C 0 0 0\n"

    text = _text_records(data)
    byts = _byte_records(data)

    assert text[0][1] is True and byts[0][1] is True    # first record whole
    assert text[-1][1] is False                         # dangling fragment
    assert byts[-1][1] is False


def test_bytes_parser_handles_leading_terminator():
    # A file that opens with a bare 'E' line: the bytes parser has a special
    # case for this degenerate leading terminator.
    data = "E\n" + db2gen.make_record("ZINC00000001")

    byts = _byte_records(data)

    assert byts[0][0] == b"E\n"
    assert byts[0][1] is True
    assert C.extract_id_bytes(byts[1][0]) == "ZINC00000001"


def test_interior_e_line_divergence_documents_open_question():
    """Pins the README "E-line" open question.

    The text parser breaks a record at ANY line starting with 'E'; the bytes
    parser only at a bare 'E' line (``\\nE\\n``). A data line like 'E7 ...'
    therefore makes them disagree — and the bytes parser is the correct one
    here, keeping the record whole. If the parsers are ever reconciled, this
    test should be updated to match the new, hardened behavior.
    """
    rec = "M ZINC00000001 x\nA 0 C 0 0 0\nE7 not-a-terminator\nC 0 1.0\nE\n"

    text = _text_records(rec)
    byts = _byte_records(rec)

    # Text parser wrongly splits at 'E7' -> two records.
    assert len(text) == 2
    assert text[0][1] is True
    # Bytes parser keeps it whole -> one correct record.
    assert len(byts) == 1
    assert byts[0][0] == rec.encode()
