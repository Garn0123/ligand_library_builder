#!/usr/bin/env python3
"""
build_ligand_relaxed.py -- db2_converter's build_ligand, with its MMFF step
restrained. Same arguments as build_ligand; pass --MMFFopt for the step to run.

    python build_ligand_relaxed.py -i in.smi ... --MMFFopt
    llb build-ligand -i in.smi ... --MMFFopt        (same, in the DB2C env)

-i may be relative: it is made absolute before build_ligand changes into
--workingpath (plain build_ligand would then fail to find it).

db2_converter (63d6656) relaxes conformers with pipeline.mmffopt: plain RDKit
MMFFOptimizeMolecule, dielectric 1, nothing held. On Owen's set (2,735
conformers from conformator) that cleared every H...H clash but moved heavy
atoms 1.83 A and tripled N+...O- contacts (469 -> 1405): charged molecules fold
in vacuum. Restrained (eps = 4r, heavy atoms free within 0.3 A, k = 100 beyond)
left 3 clashes, moved heavy atoms 0.27 A, and added 4% contacts (-> 486). See
llb mmff-test (mmff_variants.py), which this reuses.

How: build_ligand runs its molecules one by one in this process, and
pipeline.gen_conf calls the module-level pipeline.mmffopt, so replacing that one
function before build_ligand.main() runs changes only the relaxation. Nothing
in the db2_converter checkout is modified.

Settings (environment; defaults are what was measured):
    LLB_MMFF_VARIANT     eps4r_restrained   (any mmff_variants.py variant)
    LLB_MMFF_MAX_DISPL   0.3                A, free movement per heavy atom
    LLB_MMFF_K           100                kcal/mol/A^2 beyond it
    LLB_MMFF_MAX_ITERS   200                same as db2_converter's
A conformer RDKit cannot parameterize is kept unrelaxed (db2_converter's own
step would crash on it) and counted in the log.
"""
from __future__ import annotations

import logging
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from rdkit import Chem  # noqa: E402

from db2_converter import pipeline  # noqa: E402
from db2_converter.utils.utils import next_mol2_lines, update_mol2block_from_mol  # noqa: E402
from mmff_variants import relax  # noqa: E402

VARIANT = os.environ.get("LLB_MMFF_VARIANT", "eps4r_restrained")
MAX_DISPL = float(os.environ.get("LLB_MMFF_MAX_DISPL", "0.3"))
K = float(os.environ.get("LLB_MMFF_K", "100"))
MAX_ITERS = int(os.environ.get("LLB_MMFF_MAX_ITERS", "200"))
log = logging.getLogger("llb.mmff")


def mmffopt_restrained(inmol2, outmol2):
    """Drop-in for db2_converter.pipeline.mmffopt(inmol2, outmol2)."""
    blocks = [b for b in next_mol2_lines(inmol2) if b]
    out, unrelaxed = [], 0
    for block in blocks:
        mol = Chem.MolFromMol2Block("".join(block), removeHs=False)
        relaxed = relax(mol, VARIANT, MAX_ITERS, MAX_DISPL, K) if mol is not None else None
        if relaxed is None:
            out.append("".join(block))
            unrelaxed += 1
            continue
        out.append("".join(update_mol2block_from_mol(block, relaxed)))
    with open(outmol2, "w") as fh:
        fh.write("".join(out))
    log.warning(">>> llb restrained MMFF (%s, %.2f A, k=%g): %d conformer(s), %d left "
                "unrelaxed (no MMFF parameters)", VARIANT, MAX_DISPL, K, len(blocks), unrelaxed)


pipeline.mmffopt = mmffopt_restrained


def absolute_input(argv: list[str]) -> list[str]:
    """build_ligand chdirs into --workingpath before reading -i, so a relative
    -i stops resolving; make it absolute against the caller's directory."""
    out = list(argv)
    for k, a in enumerate(out[:-1]):
        if a in ("-i", "--insmi"):
            out[k + 1] = str(Path(out[k + 1]).resolve())
    return out


if __name__ == "__main__":
    from db2_converter.build_ligand import main
    sys.argv = absolute_input(sys.argv)
    sys.exit(main())
