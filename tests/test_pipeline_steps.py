"""Steps run end to end on small inputs: parents (tautomers, dedupe), triage routing
(with a fake MolGpKa), the QupKake wrapper (stand-in QupKake), and, with a real
MolGpKa, triage -> QupKake -> reenumerate --molgpka-triage -> merge -> names."""
import csv
from argparse import Namespace
from types import SimpleNamespace

import pytest
from rdkit import Chem

import prepare_parents
import qupkake_protomers
from qupkake_protomers import Site


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def test_parents_fixes_imino_tautomers_and_dedupes(tmp_path):
    (tmp_path / "in.smi").write_text(
        "N=c1[nH]cc(C(=O)N2CCC(O)(CCN)CC2)o1 ZINCi500000cVxcs\n"
        "Nc1ncc(C(=O)N2CCC(O)(CCN)CC2)o1 ZINCi500000zzzzz\n"
        "CCO ZINCi800000bbbbb\n")
    prepare_parents.main([str(tmp_path / "in.smi"), "-o", str(tmp_path / "parents"),
                          "--shard-size", "10"])
    parents = {r["parent_id"]: r for r in read_tsv(tmp_path / "parents" / "parents.tsv")}
    assert parents["ZINCi500000cVxcs"]["tautomer_fix"] == "imino_to_amino"
    assert "=N" not in parents["ZINCi500000cVxcs"]["smiles"]
    assert "ZINCi500000zzzzz" not in parents          # same molecule, amino form
    rejected = read_tsv(tmp_path / "parents" / "rejected.tsv")
    assert [r["reason"] for r in rejected] == ["duplicate_structure"]


class FakeMolGpKa:
    """predict(mol) -> (ref, sites), sites given per canonical SMILES."""
    def __init__(self, sites):
        self.sites = sites

    def predict(self, mol):
        ref = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        return ref, [SimpleNamespace(idx=i, kind=k, pka=p)
                     for i, k, p in self.sites.get(Chem.MolToSmiles(ref), [])]


def triage_args(**kw):
    a = dict(ph=[6.4, 7.4, 8.4], route_window=1.5, margin=1.0, coupling_bonds=3,
             coupled="route", keep_amide_sites=False, min_population=0.01, max_states=8,
             max_ambiguous=6, coupling="on")
    a.update(kw)
    return Namespace(**a)


def test_triage_routes_by_window_and_survives_tied_gaps():
    import triage
    # two different sites equally far from a target pH: min() over (gap, Site) went on
    # to compare the Site objects and crashed the 10k triage (TypeError, 2026-10-14)
    tie = Chem.MolToSmiles(Chem.MolFromSmiles("NCCCCCCN"))
    far = Chem.MolToSmiles(Chem.MolFromSmiles("CC(=O)O"))
    model = FakeMolGpKa({tie: [(0, "basic", 6.9), (7, "basic", 6.9)],
                         far: [(3, "acidic", 2.0)]})
    decisions, _, protomers = triage.triage_shard([(tie, "p_tie"), (far, "p_far")],
                                                  model, triage_args())
    route = {d["parent_id"]: d["route"] for d in decisions}
    assert route == {"p_tie": "qupkake", "p_far": "molgpka"}
    assert any(p.charge == -1 for p, _ in protomers)        # acetate settled by MolGpKa


def test_qupkake_wrapper_leaves_out_unembeddable_and_filters(tmp_path, fake_qupkake):
    (tmp_path / "in.smi").write_text(
        "CCN=C(N)NCCc1cncc(C)n1 guan\n[C@@H]12C[C@H]1C2 cage\nC1CNCCN1 pip\n")
    qupkake_protomers.main([str(tmp_path / "in.smi"), "-o", str(tmp_path / "qk"),
                            "--qupkake-exe", str(fake_qupkake), "--name-style", "index"])
    failed = read_tsv(tmp_path / "qk" / "qupkake_failed.tsv")
    assert [(r["name"], r["reason"]) for r in failed] == [("cage", "qupkake_embed_failed")]
    filt = dict(l.split("=") for l in (tmp_path / "qk" / "site_filter.txt").read_text().split())
    assert filt == {"site_filter_version": "6", "coupling": "1", "aromatic_n": "qupkake"}
    rows = read_csv(tmp_path / "qk" / "protomers.csv")
    pip = [int(r["net_charge"]) for r in rows if r["parent_name"] == "pip"]
    assert len(pip) == len({r["smiles"] for r in rows if r["parent_name"] == "pip"})
    assert max(pip) <= 2 and 1 in pip


def test_qupkake_wrapper_retries_once_after_a_crash(tmp_path, fake_qupkake):
    # QupKake's unseeded embed can miss a molecule the screen let through; the
    # wrapper rescreens with 100 tries and runs once more on a clean folder
    flaky = tmp_path / "flaky_qupkake"
    marker = tmp_path / "crashed_once"
    flaky.write_text(f"""#!/bin/sh
if [ ! -e {marker} ]; then touch {marker}; echo 'ValueError: Bad Conformer Id' >&2; exit 1; fi
exec {fake_qupkake} "$@"
""")
    flaky.chmod(0o755)
    (tmp_path / "in.smi").write_text("CCN guan\nC1CNCCN1 pip\n")
    qupkake_protomers.main([str(tmp_path / "in.smi"), "-o", str(tmp_path / "qk"),
                            "--qupkake-exe", str(flaky), "--name-style", "index"])
    assert marker.exists()
    rows = read_csv(tmp_path / "qk" / "protomers.csv")
    assert {r["parent_name"] for r in rows} == {"guan", "pip"}


@pytest.mark.molgpka
def test_end_to_end_with_aromatic_n_hybrid(tmp_path, fake_qupkake):
    import merge_protomers, reenumerate, triage, assign_names
    (tmp_path / "in.smi").write_text(
        "CCN=C(NCCc1cncc(C)n1)N[C@@H]1C[C@H]1C(F)F ZINCll00000k0wI1\n"
        "Cc1nnnn1CCC(=O)NCC[C@H]1CCCN1C ZINCj500000JdzBJ\n"
        "c1ccncc1 ZINCi100000aaaaa\n")
    P = tmp_path / "parents"
    prepare_parents.main([str(tmp_path / "in.smi"), "-o", str(P), "--shard-size", "10"])
    triage.main([str(P), "-o", str(tmp_path / "tr")])
    # QupKake on every parent (stand-in: every N at 9.0), as for a --validate run
    shard = read_tsv(P / "shards.tsv")[0]
    d = tmp_path / "qk" / "shard_00000"
    qupkake_protomers.main([str(P / shard["path"]), "-o", str(d), "--qupkake-exe",
                            str(fake_qupkake), "--name-style", "index", "--stereo", "carry"])
    (d / "DONE").write_text(f"shard_index=0\nshard_sha256={shard['sha256']}\nph=6.4 7.4 8.4\n"
                            f"margin=1.0\nmin_population=0.01\nmax_states=8\n")
    common = ["--parents", str(P), "--triage", str(tmp_path / "tr"), "--qupkake-parents", str(P)]
    with pytest.raises(SystemExit, match="aromatic_n=qupkake"):
        merge_protomers.main(common + ["--qupkake", str(tmp_path / "qk"), "-o", str(tmp_path / "m0")])
    reenumerate.main([str(P), str(tmp_path / "qk"), str(tmp_path / "qkh"),
                      "--molgpka-triage", str(tmp_path / "tr")])
    merge_protomers.main(common + ["--qupkake", str(tmp_path / "qkh"), "-o", str(tmp_path / "m")])
    assign_names.main(["--shards-dir", str(tmp_path / "m"), "--parents", str(P),
                       "-o", str(tmp_path / "lib"), "--jobs", "1"])
    # the swap acts on QupKake's states whatever the triage routed: check them directly
    qk = {}
    for r in read_csv(tmp_path / "qkh" / "shard_00000" / "protomers.csv"):
        qk.setdefault(r["parent_name"], []).append(r)
    # one guanidinium charge; pyrazine N's take MolGpKa's ~3 instead of the stand-in's 9
    assert {int(r["net_charge"]) for r in qk["ZINCll00000k0wI1"]} <= {0, 1}
    # tetrazole N's swapped out; the tertiary amine keeps its charge
    assert not any("[nH+]" in r["smiles"] for r in qk["ZINCj500000JdzBJ"])
    assert any(int(r["net_charge"]) == 1 for r in qk["ZINCj500000JdzBJ"])
    # pyridine alone has no base off the ring: QupKake's value stands
    assert any("[nH+]" in r["smiles"] for r in qk["ZINCi100000aaaaa"])
    swapped = {r["name"] for r in read_tsv(tmp_path / "qkh" / "shard_00000" / "hybrid_sites.tsv")}
    assert swapped == {"ZINCll00000k0wI1", "ZINCj500000JdzBJ"}
    lib = read_tsv(tmp_path / "lib" / "library.tsv")
    assert {r["parent_id"] for r in lib} == {"ZINCll00000k0wI1", "ZINCj500000JdzBJ",
                                             "ZINCi100000aaaaa"}
