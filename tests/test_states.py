"""Protonation-state model: enumerate_states (independent and coupled), dominant
states, and build_protomers' merging of symmetric states."""
import random

import pytest
from rdkit import Chem

from qupkake_protomers import (COUPLING_VERSION, Coupling, Site, build_protomers,
                               dominant_sites, enumerate_states)


class ZeroCoupling:
    def shift(self, a, b):
        return 0.0


def charge_pop(states, charge):
    return sum(p for changes, p in states if sum(d for _, d in changes) == charge)


def test_coupling_version():
    assert COUPLING_VERSION == "1"


def test_zero_penalties_reproduce_independent_states():
    rng = random.Random(1)
    for _ in range(500):
        sites = [Site(i, rng.choice(["acidic", "basic"]), rng.uniform(0, 14))
                 for i in range(rng.randint(0, 6))]
        for ph in (6.4, 7.4, 8.4):
            a = dict(enumerate_states(sites, ph, 1.0, 6))
            b = dict(enumerate_states(sites, ph, 1.0, 6, ZeroCoupling()))
            assert a.keys() == b.keys()
            assert all(abs(a[k] - b[k]) < 1e-9 for k in a)


def test_ethylenediamine_dication_fraction_matches_measured_pka():
    # Williams: 9.98 / 7.52 -> ~43% +2 at pH 7.4; micro pKa 9.68 per N
    m = Chem.MolFromSmiles("NCCN")
    states = list(enumerate_states([Site(0, "basic", 9.68), Site(3, "basic", 9.68)],
                                   7.4, 1.0, 6, Coupling(m)))
    assert charge_pop(states, 2) == pytest.approx(0.37, abs=0.05)


def test_piperazine_uses_the_ring_shift():
    # Williams: 9.82 / 5.68 -> ~16% +2 at pH 6.4
    m = Chem.MolFromSmiles("C1CNCCN1")
    states = list(enumerate_states([Site(2, "basic", 9.52), Site(5, "basic", 9.52)],
                                   6.4, 1.0, 6, Coupling(m)))
    assert charge_pop(states, 2) == pytest.approx(0.17, abs=0.05)


def test_coupling_stops_a_polyamine_at_plus_two():
    smi = "NCCNCCCNCCN"
    m = Chem.MolFromSmiles(smi)
    sites = [Site(a.GetIdx(), "basic", 9.0) for a in m.GetAtoms() if a.GetSymbol() == "N"]
    assert len(dominant_sites(sites, 7.4)) == 4                 # independent: +4
    assert len(dominant_sites(sites, 7.4, Coupling(m))) == 2    # coupled: +2


def test_opposite_charges_are_not_coupled():
    m = Chem.MolFromSmiles("NCC(=O)O")   # glycine
    cp = Coupling(m)
    assert cp.shift(Site(0, "basic", 9.6), Site(4, "acidic", 2.3)) == 0.0


def test_symmetric_states_are_one_protomer():
    m = Chem.MolFromSmiles("C1CNCCN1")
    emitted, *_ = build_protomers("pip", [m], m, [Site(2, "basic", 9.52), Site(5, "basic", 9.52)],
                                  [6.4, 7.4, 8.4], 1.0, 0.01, 8, 6, "carry", coupling=True)
    smiles = [p.smiles for p in emitted]
    assert len(smiles) == len(set(smiles))
    mono = [p for p in emitted if p.charge == 1]
    assert len(mono) == 1 and mono[0].population > 0.9      # both N's, summed
