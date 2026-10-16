"""Site filter (qupkake_protomers.IMPLAUSIBLE_SITES, filter_sites): which predicted
pKa sites are dropped before states are built. Each case is a molecule that
motivated a rule, or one the rule must leave alone. Sources: pka_triage/TRIAGE.md."""
import pytest
from rdkit import Chem

from qupkake_protomers import SITE_FILTER_VERSION, Site, filter_sites, implausible_sites


def reasons(smi, kind, explicit_h=False):
    """{atom symbol+index: rule} for one kind of site on `smi`."""
    m = Chem.MolFromSmiles(smi)
    if explicit_h:
        m = Chem.AddHs(m)
    return {f"{m.GetAtomWithIdx(i).GetSymbol()}{i}": r
            for (i, k), r in implausible_sites(m).items() if k == kind}


def test_version_is_recorded():
    # bump with every rule change; merge refuses mixed versions
    assert SITE_FILTER_VERSION == 6


@pytest.mark.parametrize("explicit_h", [False, True])
@pytest.mark.parametrize("smi,atom,rule", [
    ("CC(=O)NC", "N3", "amide_nh_acid"),          # amide anions: 22 of Owen's 50
    ("CC(=O)Nc1ccccc1", "N3", "amide_nh_acid"),    # anilide
    ("NC(=O)N", "N0", "amide_nh_acid"),            # urea
    ("CCN", "N2", "amine_nh_acid"),
    ("CCO", "O2", "alcohol_oh_acid"),
    ("OC(=O)CC(=O)O", "C3", "carbon_acid"),        # malonate [CH-] trianion (v4)
    ("CCOC(=O)CC(C)=O", "C5", "carbon_acid"),      # ethyl acetoacetate, 10.68
    ("N#CCC#N", "C2", "carbon_acid"),              # malononitrile, 11.81
])
def test_acid_sites_dropped(smi, atom, rule, explicit_h):
    assert reasons(smi, "acidic", explicit_h).get(atom) == rule


@pytest.mark.parametrize("smi,atom", [
    ("O=C1CCC(=O)N1", "N4"),         # imide, genuinely acidic
    ("CC(=O)NS(=O)(=O)C", "N3"),     # acylsulfonamide, kept
    ("CC(=O)CC(C)=O", "C3"),         # acetylacetone 8.95
    ("CC1(C)CC(=O)CC(=O)C1", "C6"),  # dimedone 5.23
    ("CC[N+](=O)[O-]", "C1"),        # nitroethane 8.6
    ("Oc1ccccc1", "O0"),             # phenol
    ("OC(=O)C", "O0"),               # carboxylic acid
])
def test_acid_sites_kept(smi, atom):
    assert atom not in reasons(smi, "acidic")


@pytest.mark.parametrize("smi,atom,rule", [
    ("CCOC", "O2", "oxygen_base"),
    ("CC(=O)NC", "N3", "amide_n_base"),
    ("C1CCSC1", "S3", "sulfur_base"),              # [SH+] seen in v2 (v6)
    ("c1ccsc1", "S3", "sulfur_base"),
    ("CCN=[N+]=[N-]", "N4", "azide_base"),         # N=[N+]=N protonated azides (v6)
    ("CCN=[N+]=[N-]", "N2", "azide_base"),
    ("C=[N+]=[N-]", "C0", "azide_base"),           # diazo
    ("CCNC(N)=NC", "N2", "amidine_amino_n_base"),  # guanidine amino N (v5)
    ("CC(=N)N", "N3", "amidine_amino_n_base"),     # acetamidine amino N
])
def test_base_sites_dropped(smi, atom, rule):
    assert reasons(smi, "basic").get(atom) == rule


@pytest.mark.parametrize("smi,atom", [
    ("CCN", "N2"), ("c1ccncc1", "N3"), ("c1c[nH]cn1", "N4"),
    ("CC(=N)N", "N2"),              # the amidine =N is the base
    ("Nc1ccccn1", "N6"),            # 2-aminopyridine ring N (aromatic C: not an amidine)
    ("Nc1ccccn1", "N0"),            # its NH2 is not touched by the amidine rule
    ("CC[N+](=O)[O-]", "N2"),       # nitro N: no rule (not basic anyway)
])
def test_base_sites_kept(smi, atom):
    assert atom not in reasons(smi, "basic")


def test_guanidine_keeps_one_site_on_the_imine():
    # ZINCll00000k0wI1: QupKake gave =N 11.73 and each NH ~10.8 -> +3 guanidines
    smi = "CCN=C(NCCc1cncc(C)n1)N[C@@H]1C[C@H]1C(F)F"
    m = Chem.MolFromSmiles(smi)
    kept, dropped = filter_sites(m, [Site(2, "basic", 11.73), Site(4, "basic", 10.88),
                                     Site(14, "basic", 10.79)])
    assert [(s.idx, s.pka) for s in kept if s.kind == "basic"] == [(2, 11.73)]
    assert {d.idx for d, r in dropped if r == "amidine_amino_n_base"} == {4, 14}


def test_amidine_site_moves_to_the_imine_when_it_had_none():
    m = Chem.MolFromSmiles("CCN=C(NCCc1cncc(C)n1)N[C@@H]1C[C@H]1C(F)F")
    kept, _ = filter_sites(m, [Site(4, "basic", 10.88), Site(14, "basic", 10.79)])
    assert [(s.idx, s.pka) for s in kept] == [(2, 10.88)]   # best one, on the =N
