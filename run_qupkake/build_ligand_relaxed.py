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
    LLB_MMFF_VARIANT     eps4r_restrained   (any mmff_variants.py variant, or
                                            db2_converter = leave its step alone)
    LLB_MMFF_MAX_DISPL   0.3                A, free movement per heavy atom
    LLB_MMFF_K           100                kcal/mol/A^2 beyond it
    LLB_MMFF_MAX_ITERS   200                same as db2_converter's
A conformer RDKit cannot parameterize is kept unrelaxed (db2_converter's own
step would crash on it) and counted in the log.
"""
from __future__ import annotations

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
    # print, not logging: db2_converter configures its own loggers, and a
    # message on another logger never reached the shard log (2026-10-14)
    print(f"   INFO: llb restrained MMFF ({VARIANT}, {MAX_DISPL:.2f} A, k={K:g}): "
          f"{len(blocks)} conformer(s), {unrelaxed} left unrelaxed (no MMFF parameters)",
          flush=True)


# "db2_converter" keeps db2_converter's own step (--mmff on, or no --MMFFopt at
# all); anything else is an mmff_variants.py variant swapped in for it.
if VARIANT != "db2_converter":
    pipeline.mmffopt = mmffopt_restrained


# --- per-molecule record --------------------------------------------------------
# build_ligand prints its per-molecule summary (time, conformers, db2 parts,
# failure) once, at the END of a shard; a killed shard loses all of it. This
# appends one row per molecule to OUTPUTPATH/llb_build.tsv as each finishes.
# llb timing summarizes them.
BUILD_COLUMNS = ["name", "started", "elapsed_s", "method", "mmff", "n_in", "n_out",
                 "n_out_roth", "n_ring", "n_db2_part", "db2_bytes", "failed", "host", "smiles"]


def timed(db2_converter_fn):
    import socket
    import time as _time

    def run(line, params, *a, **kw):
        t0 = _time.time()
        started = _time.strftime("%Y-%m-%dT%H:%M:%S")
        result, failed = None, ""
        try:
            result = db2_converter_fn(line=line, params=params, *a, **kw)
            return result
        except Exception as e:                      # recorded, then re-raised
            failed = f"EXCEPTION:{type(e).__name__}"
            raise
        finally:
            elapsed = _time.time() - t0
            lsp = line.split()
            name = lsp[1] if len(lsp) > 1 else "?"
            out = Path(params["outputpath"])
            if result is not None:
                faillist, n_in, n_out, n_roth, n_ring, n_part = result
                failed = "_".join(f[-1] for f in faillist if f) or failed
            else:
                n_in = n_out = n_roth = n_ring = n_part = "NA"
            db2 = out / f"{name}.db2.gz"
            row = [name, started, f"{elapsed:.3f}", params.get("samplopt", ""),
                   (VARIANT if params.get("MMFFopt") else "off"), n_in, n_out, n_roth,
                   n_ring, n_part, db2.stat().st_size if db2.exists() else 0, failed or "-",
                   socket.gethostname(), lsp[0] if lsp else ""]
            path = out / "llb_build.tsv"
            new = not path.exists()
            with open(path, "a") as fh:
                if new:
                    fh.write("\t".join(BUILD_COLUMNS) + "\n")
                fh.write("\t".join(map(str, row)) + "\n")
    return run


def absolute_input(argv: list[str]) -> list[str]:
    """build_ligand chdirs into --workingpath before reading -i, so a relative
    -i stops resolving; make it absolute against the caller's directory."""
    out = list(argv)
    for k, a in enumerate(out[:-1]):
        if a in ("-i", "--insmi"):
            out[k + 1] = str(Path(out[k + 1]).resolve())
    return out


if __name__ == "__main__":
    import db2_converter.build_ligand as bl
    # main() builds partial(db2_converter, ...) from this module global at run
    # time, so wrapping it here times every molecule
    bl.db2_converter = timed(bl.db2_converter)
    sys.argv = absolute_input(sys.argv)
    sys.exit(bl.main())
