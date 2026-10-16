"""Aromatic-N hybrid (QupKake's ring N copies a nearby base's pKa), the tautomer
step in `llb parents`, standardization, and QupKake's embedding pre-check."""
import pytest
from rdkit import Chem

from qupkake_protomers import (AROMATIC_N_VERSION, Site, amino_tautomer, atom_map,
                               hybrid_aromatic_n, molgpka_ref, qupkake_can_embed,
                               standardize)


def test_aromatic_n_version():
    assert AROMATIC_N_VERSION == "2"


def test_tetrazole_beside_a_weak_amine_takes_molgpka():
    # v2 run: QupKake tetrazole N 5.6-5.8 beside an amine both tools put at 5.63
    m = Chem.MolFromSmiles("C#CCN1C[C@@](c2nnnn2CCCF)(N(C)C)CCC1=O")
    qk = [Site(7, "basic", 5.61), Site(8, "basic", 5.81), Site(9, "basic", 5.79),
          Site(15, "basic", 5.63)]
    mg = [Site(7, "basic", 2.69), Site(15, "basic", 5.63)]
    sites, swaps = hybrid_aromatic_n(m, qk, mg)
    assert sorted((s.idx, s.pka) for s in sites) == [(7, 2.69), (15, 5.63)]
    assert {i for i, *_ in swaps} == {7, 8, 9}


def test_aryl_amine_does_not_trigger_the_swap():
    # 4-aminopyridine: QupKake 9.49 is right (measured 9.17), MolGpKa 6.98 is not
    m = Chem.MolFromSmiles("Nc1ccncc1")
    sites, swaps = hybrid_aromatic_n(m, [Site(0, "basic", 9.0), Site(4, "basic", 9.49)],
                                     [Site(4, "basic", 6.98)])
    assert swaps == []
    assert (4, 9.49) in [(s.idx, s.pka) for s in sites]


def test_no_base_off_the_ring_no_swap():
    m = Chem.MolFromSmiles("Cc1cnccn1")
    _, swaps = hybrid_aromatic_n(m, [Site(4, "basic", 2.4)], [Site(4, "basic", 1.8)])
    assert swaps == []


def test_molgpka_sites_map_onto_a_qupkake_mol_with_hydrogens():
    smi = "CN(C)CCc1cocn1"
    parent = Chem.MolFromSmiles(smi)
    qk_mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    assert atom_map(qk_mol, parent) is not None
    assert atom_map(molgpka_ref(smi), parent) is not None


@pytest.mark.parametrize("smi,expect", [
    ("N=c1cccc[nH]1", "Nc1ccccn1"),                                    # 2-iminopyridine
    ("Cc1c(Br)csc1C(=O)N[C@@H](C)c1csc(=N)[nH]1",
     "Cc1c(Br)csc1C(=O)N[C@@H](C)c1csc(N)n1"),                          # stereo kept
    ("CC(=O)N1CCN(c2ccnc(NCCOC(=N)O)n2)CC1",
     "CC(=O)N1CCN(c2ccnc(NCCOC(N)=O)n2)CC1"),                           # imidic acid
    ("NC(=N)c1ccccc1", None),                                          # amidine: no change
    ("N=c1ccn(C)cc1", None),                                           # no amino form
    ("CCC(=N)OCC", None),                                              # imidate ester
    ("CCOC(N)=O", None),
])
def test_amino_tautomer(smi, expect):
    out = amino_tautomer(Chem.MolFromSmiles(smi))
    assert (Chem.MolToSmiles(out) if out is not None else None) == \
        (Chem.CanonSmiles(expect) if expect else None)


@pytest.mark.parametrize("smi", ["CCN=[N+]=[N-]", "CC[N+](=O)[O-]", "c1cc[n+]([O-])cc1",
                                 "C[N+](C)(C)CC(=O)[O-]"])
def test_standardize_keeps_dipoles_and_zwitterions_neutral(smi):
    m = standardize(smi)
    assert Chem.GetFormalCharge(m) == 0


def test_qupkake_embed_precheck():
    assert qupkake_can_embed("CCN1CCN(C)CC1")
    assert not qupkake_can_embed("[C@@H]12C[C@H]1C2")   # impossible bicyclobutane
    # fails ~half of all embeds; passed 3 seeds and crashed a shard of the 50k run twice
    assert not qupkake_can_embed("Cc1noc([C@H](C)NC[C@@]23COCC[C@H]2C3(F)F)n1")
