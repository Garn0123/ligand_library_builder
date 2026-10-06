"""
llb_ids.py -- input name -> parent id -> 12-character base id. Defined once.

prepare_parents.py and assign_names.py both import this. A naming rule kept in
two places drifts, and in this project a drifted id rule has already dropped
molecules silently (DRAP stage1/ids.py).

Every input name lands in exactly one of three cases. No flag chooses between
them; the shape of the name does.

  ZINC id      ^ZINC[0-9A-Za-z]{12}$  -> kept verbatim; base id = name[4:16].
               Old unpadded numeric ids (ZINC12345) are zero-padded first:
               ZINC000000012345 is the same molecule's canonical form.
  ZINC-shaped  starts with upper-case "ZINC" but is neither of the above
               (ZINC000012345678.0, ZINC..._1) -> REJECTED, not hashed. It is
               almost certainly a damaged or suffixed ZINC id, and hashing it
               would give a known molecule a second identity next to its real one.
  anything     CHEMBL25, my_benzoic_acid, ... -> parent id "LLB" + base id, where
  else         base id = "Z" + 11 base62 characters of sha256(name).

Why "Z": ZINC22's first id character is the heavy-atom count in base62
(0-9a-zA-Z: H04 -> '4', H22 -> 'm', H28 -> 's'), and the bins stop at H49 ('N');
ZINC20 ids start with a digit. No ZINC base id starts with 'O'..'Z', so a hashed
base id can never equal a ZINC one -- disjoint by construction, not by odds.
Among hashed ids themselves: 62^11 = 5.2e19, so 10 million names collide with
probability ~1e-6, and prepare_parents.py rejects a collision if one happens.

Why foreign parents travel as "LLB" + base id rather than their own name:
QupKake writes processed/{name}.pt and db2_converter runs `rm -r {name}`, so
an arbitrary name is a path and shell hazard. The original name rides along in
parents.tsv and library.tsv (input_name).

stdlib only.
"""
from __future__ import annotations

import hashlib
import re

BASE62 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
HASH_MARK = "Z"
ZINC_ID = re.compile(r"^ZINC[0-9A-Za-z]{12}$")
ZINC_SHORT = re.compile(r"^ZINC(\d{1,11})$")
LLB_ID = re.compile(r"^LLB(" + HASH_MARK + r"[0-9A-Za-z]{11})$")


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
