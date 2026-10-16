"""`llb pipeline`'s state machine, with a fake cluster: FakeRunner plays each llb step
by writing the files the real step writes, hands out job ids, and can leave QupKake
shards unfinished or db2 tasks failed. No SLURM, no chemistry."""
import csv
import json

import pytest

import llb_pipeline as P


def write_tsv(path, rows, cols=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = cols or list(rows[0])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


class FakeRunner:
    """Parents p0..p(N-1); p0 and p1 routed to QupKake; each parent gets two names."""

    def __init__(self, run_dir, n=20, rejected=0, qk_missing_rounds=0, db2_failed_rounds=0,
                 build_fail=0, silent_missing=0):
        self.d, self.n, self.rejected = run_dir, n, rejected
        self.qk_missing_rounds, self.db2_failed_rounds = qk_missing_rounds, db2_failed_rounds
        self.build_fail, self.silent_missing = build_fail, silent_missing
        self.calls, self.sbatches, self.next_id = [], [], 100
        self.failed = {}

    def _id(self):
        self.next_id += 1
        return str(self.next_id)

    # -- llb steps ---------------------------------------------------------------
    def llb(self, step, *args):
        self.calls.append((step, *args))
        d = self.d
        if step == "parents":
            pids = [f"p{i}" for i in range(self.n)]
            write_tsv(d / "parents/parents.tsv", [
                {"parent_id": p, "input_name": p, "input_smiles": "CC", "source": f"in.smi:{i}",
                 "heavy_atoms": str(18 + i % 2), "tautomer_fix": ""} for i, p in enumerate(pids)])
            write_tsv(d / "parents/rejected.tsv",
                      [{"input_name": f"bad{i}", "reason": "unparsable", "detail": "x"}
                       for i in range(self.rejected)], ["input_name", "reason", "detail"])
            write_tsv(d / "parents/shards.tsv",
                      [{"index": "0", "path": "shards/shard_00000.smi", "n_parents": str(self.n)}])
        elif step == "triage":
            write_tsv(d / "triage/decisions.tsv", [
                {"parent_id": f"p{i}", "route": "qupkake" if i < 2 else "molgpka",
                 "reason": "near_window" if i < 2 else ""} for i in range(self.n)])
            write_tsv(d / "triage/routed/shards.tsv",
                      [{"index": "0", "path": "shards/shard_00000.smi", "n_parents": "2"}])
        elif step == "reenumerate":
            write_tsv(d / "protomers_qk_h/shard_00000/hybrid_sites.tsv",
                      [{"name": "p0", "idx": "3", "qupkake_pka": "9", "molgpka_pka": "2"}])
        elif step == "merge":
            (d / "merged").mkdir(exist_ok=True)
        elif step in ("validate-names", "compare"):
            pass
        elif step == "check-stereo":
            write_tsv(d / "library/library.stereo_split.tsv",
                      [{"name": "C00000000000001C", "n_isomers": "2", "smiles": "CC"}])
            write_tsv(d / "library/library.stereo_crashed.tsv", [], ["name", "smiles", "error"])
        elif step == "verify-db2":
            names = [r["name"] for r in csv.DictReader(open(d / "library/library.tsv"), delimiter="\t")]
            excluded = {"C00000000000001C"}
            rows, failed, silent = [], 0, 0
            for nm in names:
                pid = "p" + str(int(nm[1:15]) // 10)
                if nm in excluded:
                    st, det = "MISSING", "no db2 and no faillist entry"
                elif failed < self.build_fail:
                    st, det, failed = "BUILD_FAILED", "1generate", failed + 1
                elif silent < self.silent_missing:
                    st, det, silent = "MISSING", "no db2 and no faillist entry", silent + 1
                else:
                    st, det = "OK", ""
                rows.append({"protomer_name": nm, "parent": pid, "expected_charge": "0",
                             "status": st, "detail": det, "db2_file": "", "db2_name": "",
                             "amsol_charge": ""})
            with open(d / "db2_verification.csv", "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
                w.writeheader()
                w.writerows(rows)
        else:
            raise AssertionError(f"unexpected llb step {step}")
        return 0

    def llb_capture(self, step, *args):
        self.calls.append((step, *args))
        if step == "submit-qupkake":
            if self.qk_missing_rounds > 0:
                self.qk_missing_rounds -= 1
            else:
                (self.d / "protomers_qk/shard_00000").mkdir(parents=True, exist_ok=True)
                (self.d / "protomers_qk/shard_00000/DONE").write_text("x")
                write_tsv(self.d / "protomers_qk/shard_00000/qupkake_failed.tsv", [],
                          ["name", "reason", "smiles"])
            return 0, f"Submitted batch job {self._id()}\n"
        if step == "submit-db2":
            jid = self._id()
            if self.db2_failed_rounds > 0:
                self.db2_failed_rounds -= 1
                self.failed[jid] = 3
            return 0, f"Submitted batch job {jid}\n"
        if step == "timing":
            return 0, "== db2 build: fake\n"
        raise AssertionError(step)

    # -- SLURM ---------------------------------------------------------------------
    def sbatch(self, args):
        self.sbatches.append(args)
        jid = self._id()
        wrap = args[args.index("--wrap") + 1]
        if " names " in wrap:                         # the names job: write the library
            lib = []
            for i in range(self.n):
                for s in range(2):
                    lib.append({"name": f"C{i * 10 + s:014d}C", "parent_id": f"p{i}",
                                "protomer_index": "0", "stereo_index": str(s),
                                "site_source": "qupkake" if i < 2 else "molgpka", "smiles": "CC"})
            write_tsv(self.d / "library/library.tsv", lib)
            (self.d / "library/library.smi").write_text(
                "".join(f"CC {x['name']}\n" for x in lib))
            write_tsv(self.d / "library/dropped.tsv", [],
                      ["parent_id", "protomer_index", "reason", "detail"])
            (self.d / "library/stats.json").write_text("{}")
        return jid

    def active(self, ids):
        return False

    def failed_tasks(self, ids):
        return sum(self.failed.get(i, 0) for i in ids)

    def cancel(self, ids):
        self.cancelled = ids


CLUSTER = {"LLB_ACCOUNT": "acct", "LLB_PARTITION": "standard", "LLB_MAIL": "me@x.edu"}


def new_run(tmp_path, **fake):
    (tmp_path / "in.smi").write_text("CC p0\n")
    run_dir = tmp_path / "run"
    settings = json.loads(json.dumps(P.DEFAULTS))
    settings["stop"]["checks"] = 0.05     # the fake library is tiny: 1 excluded of 40 names
    fr = FakeRunner(run_dir, **fake)
    run = P.create_run(run_dir, [str(tmp_path / "in.smi")], settings, CLUSTER,
                       {"llb_commit": "abc"}, fr)
    return run, fr


def advance_until_idle(run, fr, limit=20):
    """Each call is one control job; stop when nothing new was submitted."""
    rcs = []
    for _ in range(limit):
        before = len(fr.sbatches) + len(fr.calls)
        run = P.Run(run.path, fr)
        rcs.append(P.advance(run))
        if rcs[-1] == 3 or run.state.get("finished") or len(fr.sbatches) + len(fr.calls) == before:
            break
    return P.Run(run.path, fr), rcs


def test_happy_path_runs_every_stage_and_writes_the_records(tmp_path):
    run, fr = new_run(tmp_path)
    run, rcs = advance_until_idle(run, fr)
    assert set(rcs) == {0}
    assert all(run.stage(s)["status"] == "done" for s in P.STAGES)
    assert run.f("REPORT.md").exists() and "llb mol2-library" in run.f("REPORT.md").read_text()
    ledger = list(csv.DictReader(open(run.f("ledger.tsv")), delimiter="\t"))
    assert len(ledger) == 20 and {x["fate"] for x in ledger} == {"built"}
    assert ledger[0]["aromatic_n_swapped"] == "1" and ledger[0]["n_excluded_checks"] == "1"
    att = list(csv.DictReader(open(run.f("attrition.tsv")), delimiter="\t"))
    assert att[-1] == {"heavy_atoms": "all", "inputs": "20", "parents": "20", "named": "20",
                       "passed_checks": "20", "built": "20"}
    build = run.f("library/build.smi").read_text()
    assert "C00000000000001C" not in build and build.count("\n") == 39
    controls = [a for a in fr.sbatches if "control" in a[a.index("-J") + 1]]
    assert len(controls) == 3                         # after qupkake, names, db2
    assert "--mail-type=END,FAIL" in controls[-1] and "--mail-type=FAIL" in controls[0]
    assert any(c[0] == "submit-db2" and "--nshards" in c for c in fr.calls)


def test_threshold_stop_then_resume_rechecks_without_rerunning(tmp_path):
    run, fr = new_run(tmp_path, rejected=1)            # 1 of 21 = 4.8% > 1%
    run, rcs = advance_until_idle(run, fr)
    assert rcs[-1] == 3 and run.stage("parents")["status"] == "stopped"
    assert "stop threshold" in run.stage("parents")["message"]

    class A:   # resume with a higher threshold, in this process
        pass
    a = A(); a.run, a.set, a.here = str(run.path), ["stop.parents=0.06"], False
    orig = P.Run.__init__

    def init(self, path, runner=None):
        orig(self, path, fr)
    P.Run.__init__ = init
    try:
        P.cmd_resume(a)                                # submits a control job
    finally:
        P.Run.__init__ = orig
    run, rcs = advance_until_idle(P.Run(run.path, fr), fr)
    assert sum(1 for c in fr.calls if c[0] == "parents") == 1   # not rerun
    assert run.state.get("finished") and run.state["overrides"][0]["key"] == "stop.parents"


def test_unfinished_qupkake_is_retried_once_at_twice_the_time(tmp_path):
    run, fr = new_run(tmp_path, qk_missing_rounds=1)
    run, rcs = advance_until_idle(run, fr)
    subs = [c for c in fr.calls if c[0] == "submit-qupkake"]
    assert [c[-1] for c in subs] == ["--time=00:30:00", "--time=01:00:00"]
    assert run.state.get("finished")


def test_twice_unfinished_qupkake_stops(tmp_path):
    run, fr = new_run(tmp_path, qk_missing_rounds=2)
    run, rcs = advance_until_idle(run, fr)
    assert rcs[-1] == 3
    assert run.stage("qupkake")["phase"] == "collect"


def test_failed_db2_tasks_are_retried_and_old_failures_forgotten(tmp_path):
    run, fr = new_run(tmp_path, db2_failed_rounds=1)
    run, rcs = advance_until_idle(run, fr)
    assert [c[-1] for c in fr.calls if c[0] == "submit-db2"] == ["01:00:00", "02:00:00"]
    assert run.state.get("finished")


def test_a_silent_db2_loss_stops_the_run(tmp_path):
    run, fr = new_run(tmp_path, silent_missing=1)
    run, rcs = advance_until_idle(run, fr)
    assert rcs[-1] == 3 and run.stage("verify")["status"] == "stopped"
    assert "missing with no reason" in run.stage("verify")["message"]


def test_parent_loss_at_db2_above_threshold_stops(tmp_path):
    run, fr = new_run(tmp_path, build_fail=4)          # p0 (other name excluded) and p1 lose all -> 2/20 = 10%
    run, rcs = advance_until_idle(run, fr)
    assert rcs[-1] == 3 and "named parents with nothing built" in run.stage("verify")["message"]


def test_settings_and_times():
    assert P.double_time("00:30:00") == "01:00:00"
    assert P.double_time("18:00:00") == "1-12:00:00"
    s = json.loads(json.dumps(P.DEFAULTS))
    P.set_setting(s, "stop.db2", "0.08")
    P.set_setting(s, "db2_names_per_task", "30")
    P.set_setting(s, "validate", "yes")
    assert s["stop"]["db2"] == 0.08 and s["db2_names_per_task"] == 30 and s["validate"] is True
    with pytest.raises(SystemExit):
        P.set_setting(s, "stop.nonsense", "1")


def test_status_and_cancel_commands(tmp_path, capsys):
    run, fr = new_run(tmp_path)
    advance_until_idle(run, fr)

    class A:
        pass
    a = A(); a.run = str(run.path)
    orig = P.Run.__init__
    P.Run.__init__ = lambda self, path, runner=None: orig(self, path, fr)
    try:
        P.cmd_status(a)
        P.cmd_cancel(a)
    finally:
        P.Run.__init__ = orig
    out = capsys.readouterr().out
    assert "verify" in out and "finished" in out and "llb mol2-library" in out
    assert fr.cancelled                      # every job and control job id
