#!/usr/bin/env python3
"""
mmff_variants.py -- how much does each kind of MMFF relaxation fix, and collapse?

    llb mmff-test db2_run/out/0003/conformer.NAME.fixed.mol2
    llb mmff-test conformers.sdf --variants plain,eps4r_restrained

db2_converter's --MMFFopt (63d6656, pipeline.mmffopt) is plain RDKit
MMFFOptimizeMolecule: MMFF94, dielectric 1, 200 iterations, nothing held. With
full charges in vacuum that pulls ammonium onto carboxylate: the clashes go,
but the conformers fold up. This runs the same unrelaxed conformers through
several variants and reports, per variant, what changed:

  closest_HH     closest H...H more than two bonds apart (the db2 clash rule
                 uses 1.70 A); median, and how many conformers stay under 1.70
  heavy_rmsd     how far the heavy atoms moved from the input, A
  rg_change      radius of gyration after / before - 1 (negative = shrank)
  salt_bridges   N+ ... O- pairs within 3.2 A, summed over conformers
                 (before -> after)

Variants (all MMFF94s properties, up to --max-iters steps):
  plain                db2_converter's own: dielectric 1, nothing held
  eps4r                distance-dependent dielectric, epsilon = 4r
  noele                electrostatics off
  honly                heavy atoms fixed; only hydrogens move
  restrained           heavy atoms in a flat-bottom well (--max-displ, --k),
                       dielectric 1
  eps4r_restrained     both

Reads mol2 (the way db2_converter does: RDKit MolFromMol2Block, H kept) or SDF.
Writes nothing unless --write VARIANT:OUT.sdf asks for one variant's relaxed
conformers, to look at.
"""
from __future__ import annotations

import argparse
import statistics
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeForceField
from rdkit.Chem import rdForceFieldHelpers
from rdkit.ForceField import rdForceField  # noqa: F401 -- registers MMFFMolProperties

RDLogger.DisableLog("rdApp.*")
VARIANTS = ["plain", "eps4r", "noele", "honly", "restrained", "eps4r_restrained"]


def read(path):
    if path.endswith(".sdf"):
        return [m for m in Chem.SDMolSupplier(path, removeHs=False) if m is not None]
    mols, block = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("@<TRIPOS>MOLECULE") and block:
                mols.append("".join(block))
                block = []
            if line.startswith("#"):
                continue
            block.append(line)
    if block:
        mols.append("".join(block))
    out = []
    for b in mols:
        m = Chem.MolFromMol2Block(b, removeHs=False)
        if m is not None:
            out.append(m)
    return out


def closest_hh(mol):
    conf = mol.GetConformer()
    dm = Chem.GetDistanceMatrix(mol)
    hs = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
    best = float("inf")
    for i in range(len(hs)):
        pi = conf.GetAtomPosition(hs[i])
        for j in range(i + 1, len(hs)):
            if dm[hs[i]][hs[j]] > 2:
                best = min(best, (pi - conf.GetAtomPosition(hs[j])).Length())
    return best


def heavy_rmsd(a, b):
    ca, cb = a.GetConformer(), b.GetConformer()
    idx = [x.GetIdx() for x in a.GetAtoms() if x.GetAtomicNum() > 1]
    return (sum((ca.GetAtomPosition(i) - cb.GetAtomPosition(i)).LengthSq() for i in idx)
            / len(idx)) ** 0.5


def rg(mol):
    conf = mol.GetConformer()
    idx = [x.GetIdx() for x in mol.GetAtoms() if x.GetAtomicNum() > 1]
    c = rdMolTransforms.ComputeCentroid(conf, ignoreHs=True)
    return (sum((conf.GetAtomPosition(i) - c).LengthSq() for i in idx) / len(idx)) ** 0.5


def salt_bridges(mol, cut=3.2):
    conf = mol.GetConformer()
    pos = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "N" and a.GetFormalCharge() > 0]
    neg = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "O" and a.GetFormalCharge() < 0]
    return sum(1 for p in pos for n in neg
               if (conf.GetAtomPosition(p) - conf.GetAtomPosition(n)).Length() < cut)


def relax(mol, variant, max_iters, max_displ, k):
    m = Chem.Mol(mol)
    props = rdForceFieldHelpers.MMFFGetMoleculeProperties(m, mmffVariant="MMFF94s")
    if props is None:
        return None
    if variant == "plain":
        props = rdForceFieldHelpers.MMFFGetMoleculeProperties(m, mmffVariant="MMFF94")
    if variant in ("eps4r", "eps4r_restrained"):
        props.SetMMFFDielectricModel(2)   # RDKit: 1 constant, 2 distance-dependent (eps = D*r)
        props.SetMMFFDielectricConstant(4.0)
    if variant == "noele":
        props.SetMMFFEleTerm(False)
    ff = MMFFGetMoleculeForceField(m, props)
    heavy = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1]
    if variant == "honly":
        for i in heavy:
            ff.AddFixedPoint(i)
    if variant in ("restrained", "eps4r_restrained"):
        for i in heavy:
            ff.MMFFAddPositionConstraint(i, max_displ, k)
    ff.Minimize(maxIts=max_iters)
    return m


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="unrelaxed conformers: mol2 or sdf")
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--max-iters", type=int, default=200, help="db2_converter uses RDKit's 200")
    ap.add_argument("--max-displ", type=float, default=0.3,
                    help="restrained: heavy atoms move freely up to this, A")
    ap.add_argument("--k", type=float, default=100.0,
                    help="restrained: force constant beyond --max-displ, kcal/mol/A^2")
    ap.add_argument("--write", default=None, metavar="VARIANT:OUT.sdf",
                    help="write one variant's relaxed conformers to an SDF")
    args = ap.parse_args(argv)

    mols = [m for p in args.inputs for m in read(p)]
    if not mols:
        raise SystemExit("no conformers read (RDKit could not parse the input)")
    variants = args.variants.split(",")
    write_v, write_path = (args.write.split(":", 1) + [None])[:2] if args.write else (None, None)
    writer = Chem.SDWriter(write_path) if write_path else None

    base_hh = [closest_hh(m) for m in mols]
    base_sb = sum(salt_bridges(m) for m in mols)
    print(f"{len(mols)} conformer(s); input: closest H-H median {statistics.median(base_hh):.2f} A, "
          f"{sum(d < 1.70 for d in base_hh)} under 1.70; salt bridges {base_sb}")
    print(f"{'variant':<18} {'HH median':>9} {'<1.70':>6} {'heavy RMSD':>10} "
          f"{'Rg change':>9} {'salt bridges':>13}")
    for v in variants:
        hh, rmsd, drg, sb, failed = [], [], [], 0, 0
        for m in mols:
            r = relax(m, v, args.max_iters, args.max_displ, args.k)
            if r is None:
                failed += 1
                continue
            hh.append(closest_hh(r))
            rmsd.append(heavy_rmsd(m, r))
            drg.append(rg(r) / rg(m) - 1)
            sb += salt_bridges(r)
            if writer and v == write_v:
                writer.write(r)
        if not hh:
            print(f"{v:<18} no MMFF parameters for this molecule")
            continue
        print(f"{v:<18} {statistics.median(hh):>9.2f} {sum(d < 1.70 for d in hh):>6} "
              f"{statistics.median(rmsd):>10.2f} {100 * statistics.median(drg):>8.1f}% "
              f"{base_sb:>6} -> {sb:<5}" + (f"  ({failed} unparameterized)" if failed else ""))
    if writer:
        writer.close()
        print(f"-> {write_path} ({write_v})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
