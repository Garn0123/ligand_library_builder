"""`llb names` (assign_names.py): stereo closure, the 3D embedding check and its
last-isomer guard, invertomer marking, and the name contract."""
import csv

import pytest

import assign_names
from assign_names import closure, embeddable, invertomer_key, needs_embed_check
from qupkake_protomers import PROTOMER_COLUMNS
from rdkit import Chem


def test_cage_isomers_that_cannot_exist_are_not_embeddable():
    # 3chemistrycheck failures on the 10k: conformator built the other isomer
    assert not embeddable("CC(=O)NNC(=O)C[NH2+][C@@H]1[C@@H]2C[C@@H]3[C@H](C2)[C@@]31CO")
    assert not embeddable("C[C@H](Cc1cnc[nH]1)C(=O)N1C[C@H]2C[C@@]1(c1ccccc1)C2")


def test_strained_but_real_isomers_are_embeddable():
    # fast check rejected these; the escalation must rescue them
    assert embeddable("O=C(NC[C@@H]1CCN1Cc1cscc1Br)[C@@H]1CCC[C@H]2C[C@@H]21")
    assert embeddable("CNC(=O)c1ccc(NC[C@]23CCCC[C@@H]2C3)c(F)c1")


def test_cage_closure_keeps_only_buildable_isomers():
    iso = closure("CC(=O)NNC(=O)CNC1C2CC3C(C2)C31CO")
    assert len(iso) == 16
    assert sum(embeddable(s) for s in iso) == 2


def test_only_fused_or_bridged_centres_are_checked():
    assert not needs_embed_check(Chem.MolFromSmiles("C[C@H](N)CC"))
    assert needs_embed_check(Chem.MolFromSmiles("C[C@@H]1[C@H]2CC[C@@H](C2)C1"))


def test_invertomer_key_ignores_only_nh_plus_centres():
    a = "C[N@@H+]1CCN(C)[C@H]2C[C@@H]21"
    b = "C[N@H+]1CCN(C)[C@H]2C[C@@H]21"
    assert invertomer_key(a) == invertomer_key(b)
    assert invertomer_key("C[C@H](N)CC") != invertomer_key("C[C@@H](N)CC")


def write_shard(root, rows):
    d = root / "shard_00000"
    d.mkdir(parents=True)
    (d / "DONE").write_text("shard_index=0\n")
    with open(d / "protomers.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=PROTOMER_COLUMNS + ["site_source"], lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({**{c: "" for c in PROTOMER_COLUMNS}, "site_source": "molgpka",
                        "ph_values": "7.4", "population_estimate": "1.0", **r})


def run_names(tmp_path, rows):
    write_shard(tmp_path / "merged", rows)
    assert assign_names.main(["--shards-dir", str(tmp_path / "merged"), "-o",
                              str(tmp_path / "library"), "--jobs", "1"]) in (0, None)
    with open(tmp_path / "library" / "library.tsv") as fh:
        lib = list(csv.DictReader(fh, delimiter="\t"))
    with open(tmp_path / "library" / "dropped.tsv") as fh:
        dropped = list(csv.DictReader(fh, delimiter="\t"))
    return lib, dropped


def row(pid, smi, p=0, s=0):
    return {"protomer_name": f"{pid}_{p}_{s}", "parent_name": pid, "smiles": smi,
            "protomer_index": str(p), "stereo_index": str(s)}


def test_a_protomer_never_loses_its_last_isomer(tmp_path):
    # bicyclo[1.1.1]pentane: RDKit embeds no labelling, conformator builds it (v1 did)
    lib, dropped = run_names(tmp_path, [
        row("ZINCib00000ABfDp", "NCC(CN)NC(=O)C12CC(CC(F)(F)F)(C1)C2")])
    assert {r["parent_id"] for r in lib} == {"ZINCib00000ABfDp"}
    assert not [d for d in dropped if d["reason"] == "stereo_cannot_embed"]


def test_impossible_isomers_dropped_when_a_buildable_one_exists(tmp_path):
    lib, dropped = run_names(tmp_path, [
        row("ZINCi400000II2VT", "CC(=O)NNC(=O)CNC1C2CC3C(C2)C31CO")])
    assert len(lib) == 2
    assert sum(d["reason"] == "stereo_cannot_embed" for d in dropped) == 14
    assert all(embeddable(r["smiles"]) for r in lib)


def test_names_follow_the_contract_and_mark_invertomers(tmp_path):
    lib, _ = run_names(tmp_path, [
        row("ZINCim000001Q4BG", "C[NH+]1CCN(c2ncc(Cl)cc2Br)[C@H]2C[C@@H]21")])
    names = [r["name"] for r in lib]
    assert all(len(n) == 16 and n[0] == n[-1] == "C" for n in names)
    assert len(names) == 2                                   # the two N-H+ configurations
    inv = [r for r in lib if r["invertomer_of"]]
    assert len(inv) == 1 and inv[0]["invertomer_of"] in names
