"""
naming.py -- MolSanitizer's state ids -> NAMING_CONTRACT.md names (DRAP repo).

    C  <12: base id>  <state>  <stereo>  C          exactly 16 characters

What MolSanitizer does with an input id (msani/api.py, moltransform/stereoisomers.py,
0.7.4):

    tautomer x protomer states   ID_1, ID_2, ...        (_01.. when >= 10 states)
    stereoisomers of a state     ID_1.1, ID_1.2, ...    (.01.. when >= 10)
    one state / one isomer       no suffix

and its db2 writer keeps name[-16:] (db2/layout.py, serialize_db2). A 16-character
ZINC22 id with any suffix therefore loses its leading characters inside the db2:
the collapse NAMING_CONTRACT.md invariant 1 exists to prevent. So MolSanitizer
never sees a contract name in the 2D step, and never invents one in the 3D step:

  2D   input id = the parent id (ZINC id verbatim, or LLB + hashed base id);
       its suffixes are PARSED here (parse_state_id) and thrown away.
  name state index  -> base62, renumbered 0.. per parent in MolSanitizer's order;
       stereo index -> base62, renumbered 0.. per state.
  3D   input id = the 16-character contract name, with --no-stereoisomers, so
       MolSanitizer adds no suffix and name[-16:] is the whole name.

The "protomer" field of the contract therefore numbers tautomer x protomer
states, not protomers alone. Nothing downstream reads meaning from it
(NAMING_CONTRACT.md: "Nothing downstream reads meaning from the index").

The parent-id rules are a copy of run_qupkake/llb_ids.py at d84b045 (the one
definition there; tests/test_naming.py fails if the two drift while both exist).

stdlib only, except invertomer_key (RDKit, imported when called).
"""
from __future__ import annotations

import hashlib
import re

BASE62 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
HASH_MARK = "Z"
ZINC_ID = re.compile(r"^ZINC[0-9A-Za-z]{12}$")
ZINC_SHORT = re.compile(r"^ZINC(\d{1,11})$")
LLB_ID = re.compile(r"^LLB(" + HASH_MARK + r"[0-9A-Za-z]{11})$")
CONTRACT = re.compile(r"^C[0-9A-Za-z]{14}C$")
NAME_LEN = 16


# ---- parent ids (copy of run_qupkake/llb_ids.py) ------------------------------

def hashed_base(name: str) -> str:
    """'Z' + 11 base62 characters of sha256(name). Deterministic, shard-independent."""
    n = int.from_bytes(hashlib.sha256(name.encode("utf-8")).digest(), "big")
    out = []
    for _ in range(11):
        n, r = divmod(n, 62)
        out.append(BASE62[r])
    return HASH_MARK + "".join(out)


def parent_id(name: str) -> tuple[str | None, str]:
    """(parent id, kind); parent id is None when the name is rejected.

    kind: 'zinc', 'zinc_padded', 'hashed', or the rejection reason
    'malformed_zinc_id' / 'no_name'.
    """
    if not name:
        return None, "no_name"
    if ZINC_ID.match(name):
        return name, "zinc"
    m = ZINC_SHORT.match(name)
    if m:
        return f"ZINC{int(m.group(1)):012d}", "zinc_padded"
    if name.startswith("ZINC"):
        return None, "malformed_zinc_id"
    return "LLB" + hashed_base(name), "hashed"


def base_id_of(pid: str) -> str:
    """The 12-character base id carried by a parent id."""
    if ZINC_ID.match(pid):
        return pid[4:16]
    m = LLB_ID.match(pid)
    if m:
        return m.group(1)
    raise ValueError(f"not a parent id: {pid!r}")


def is_zinc(pid: str) -> bool:
    return bool(ZINC_ID.match(pid))


# ---- contract names -------------------------------------------------------------

def contract_name(base: str, state: int, stereo: int) -> str:
    if not (0 <= state < len(BASE62) and 0 <= stereo < len(BASE62)):
        raise ValueError(f"index out of base62 range: state {state}, stereo {stereo}")
    return f"C{base}{BASE62[state]}{BASE62[stereo]}C"


def base_id(name: str) -> str:
    return name[1:13]


def parse_state_id(pid: str, msani_id: str) -> tuple[int, int] | None:
    """MolSanitizer's id for one output row -> (state number, stereo number).

    Numbers are MolSanitizer's own (1-based; 0 when the suffix is absent). None
    when the id is not pid plus suffixes -- a format change to stop on, not guess.
    Parent ids are [A-Za-z0-9] only, so '_' and '.' can only be MolSanitizer's.
    """
    m = re.fullmatch(re.escape(pid) + r"(?:_(\d+))?(?:\.(\d+))?", msani_id)
    if not m:
        return None
    return int(m.group(1) or 0), int(m.group(2) or 0)


def assign_indices(keys: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """(state, stereo) numbers as MolSanitizer wrote them -> dense 0-based indices.

    States are renumbered in MolSanitizer's order, isomers within each state the
    same way, so 'ID_2.1' and 'ID_2.2' become (1, 0) and (1, 1) when ID_1 exists.
    """
    states = sorted({s for s, _ in keys})
    state_ix = {s: i for i, s in enumerate(states)}
    stereo_ix: dict[int, dict[int, int]] = {}
    for s in states:
        isomers = sorted({t for ss, t in keys if ss == s})
        stereo_ix[s] = {t: i for i, t in enumerate(isomers)}
    return [(state_ix[s], stereo_ix[s][t]) for s, t in keys]


def invertomer_key(smi: str) -> str:
    """SMILES with the configuration of every N-H+ centre erased.

    Same rule as run_qupkake/assign_names.py. Two isomers sharing this key differ
    only at protonated-amine centres, which interconvert in solution, so the
    matrix keeps the better score of the pair (DRAP stage4/collapse_invertomers.py
    reads library.tsv's invertomer_of). Quaternary N+ (no H) cannot invert and is
    left alone. MolSanitizer says it samples both N configurations inside one
    entry's conformers (paper sec. 2.5); if its stereo step still emits them as
    separate isomers, this is what links them.
    """
    from rdkit import Chem
    m = Chem.MolFromSmiles(smi)
    for a in m.GetAtoms():
        if a.GetSymbol() == "N" and a.GetFormalCharge() == 1 and a.GetTotalNumHs() > 0:
            a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    return Chem.MolToSmiles(m)
