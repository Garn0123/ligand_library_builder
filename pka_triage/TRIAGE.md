# pKa triage: MolGpKa first, QupKake only where it matters

Branch `feat/pka-triage`. Status 2026-10-08: environment + MolGpKa port done and
checked; triage step, merge and validation harness not written yet.

## Why

QupKake on one core measured about 3.6 min/molecule on Owen's 10 largest
(2026-10-08), almost all of it xtb. A million molecules would take on the order
of 60,000 CPU-hours. MolGpKa (2D GCN) runs at **4.5 ms/molecule on one thread**
(port below, 960 drug-like molecules, laptop core).

## What exists on this branch

| file | what |
|---|---|
| `environments/llb_pka.yml` | new CPU-only env (python 3.11, rdkit 2024.03.5 = db2_converter's, torch 2.4 `cpu*` build, PyG only for the check). Solved for osx-arm64 (installed and run) and linux-64 (dry run: `cpu_mkl` torch, no CUDA). |
| `pka_triage/molgpka.py` | MolGpKa network + featurizer in plain torch. Weights and SMARTS read from a MolGpKa clone (`MOLGPKA_DIR`), sha256-checked. `--check` runs the upstream code beside it: max difference **3.8e-6 pKa** over 16 molecules / 47 sites. |
| `llb molgpka`, `PKA_ENV`, `MOLGPKA_DIR` | dispatcher step and config keys. |

Upstream problems the port fixes or records:
- upstream reloads both weight files for every molecule, and finds its SMARTS
  table relative to the current directory (works only from `src/`);
- its H-bond donor/acceptor feature columns are **always 0** (an `int in
  list-of-tuples` test). The weights were trained that way, so the port keeps
  them 0; "fixing" them moves predictions by up to 3 units;
- acid sites are reported by H index, once per H (NH2 twice). The port maps
  them to the heavy atom and collapses duplicates;
- `torch_scatter` and the removed `torch_geometric.data.DataLoader` are not needed.

## Where the swap goes

`run_qupkake/qupkake_protomers.py` already separates the two halves:

1. **site source**: QupKake SDF -> `{name: (mol, [Site(idx, kind, pka)])}`
   (`parse_qupkake_sdf`);
2. **state model**: Henderson-Hasselbalch per site, `--margin`, population
   products, `--max-states`, and `apply_changes`, which does the H/charge
   bookkeeping on the heavy atom at `idx`.

MolGpKa produces exactly half 1's shape (heavy-atom index, acidic/basic, pKa), so
the triage step is a second site source feeding the same state model, plus a
routing decision. **Dimorphite-DL is not needed to build protomers**:
`apply_changes` already does the work, keyed by the atom index MolGpKa returns.
Dimorphite matches its own SMARTS, so its sites would have to be mapped back to
MolGpKa's atoms. It stays in the env as an independent third opinion for the
harness.

## Things that assume QupKake, or are wrong for both tools

1. **Neither tool, as we run it, models site-site coupling.** QupKake predicts
   every site on the input molecule; each site's conjugate is built from that one
   state. Our state model then multiplies independent per-site fractions. MolGpKa
   is the same. So routing rule (b), "coupled sites -> QupKake", does not buy
   coupling. Measured example: piperazine, MolGpKa **10.46 on both N**, so the
   state model locks the **dication** at every target pH. Experimentally the
   second pKa is about 5.6, so the monocation dominates at 7.4. Any fix is new
   work in either path:
   - sequential re-prediction on the singly-protonated state. Cheap with
     MolGpKa, but upstream un-charges its input and was trained on neutral
     molecules, so whether it can predict on a cation is an open question;
   - or a rule-based penalty for sites within N bonds.
   The harness should report coupled pairs separately so we see the damage.
2. **Two thresholds, not one.** `--margin` (default 1.0) decides which sites are
   enumerated both ways. The proposed 1.5 decides which *molecules* go to
   QupKake. Keep them as separate config values.
3. **Atom-index frame.** QupKake's `idx` refers to its own re-parse of the
   input SMILES; MolGpKa's to the un-charged canonical re-parse. Inputs are
   already neutral and canonical after `standardize`, so the two should
   coincide, but `skeleton_matches` must be checked per molecule as it is now,
   not assumed.
4. **Different site sets.** MolGpKa finds sites with 143 SMARTS. QupKake's site
   detector is its own. Nifedipine has no MolGpKa sites. "No sites" has to be a
   recorded outcome (`note=no_sites_predicted`, `site_source=molgpka`), and the
   harness compares site sets, not only pKa values.
5. **Shard outputs.** `assign_names.py` reads `protomers.csv`,
   `qupkake_failed.tsv` and `DONE` per shard, and groups by `parent_name`. A
   parent present in both a MolGpKa shard and a QupKake shard would collide. The
   merge (QupKake wins for routed parents) must happen before naming, keeping
   the same parent ids. `library.tsv` should gain a `site_source` column.
6. **`audit_qupkake`, the xtb failure modes, and the `XTBPATH`/exe provenance in
   `DONE`** are QupKake-only. The MolGpKa path records the weight sha256 instead.
7. **What agreement means.** MolGpKa was trained, per its paper, on ChEMBL pKa
   values computed by another program, not measured ones. MolGpKa-vs-QupKake
   agreement therefore measures model against model. Use it to tune the
   *routing* threshold, and keep a small experimental set to see which side is
   right. Spotted already: tyrosine phenol at 6.45 from MolGpKa (real about 10).
   It falls inside the window, so it would be routed.
8. **Config location.** The pH list and thresholds belong in this repo's run
   configuration (recorded in each shard's `DONE`), not in db2_converter. That
   repo is pinned upstream code (hnlab @ 63d6656) and only builds conformers.
9. **Caching by ionizable fragment** would change answers. MolGpKa is a
   whole-graph model: 5 GCN layers plus attention over every atom, so a site's
   pKa depends on the whole molecule. At 4.5 ms/molecule it is not needed either
   (1M molecules is about 1.3 CPU-hours). Exact-structure dedupe already happens
   in `prepare_parents.py`.

## Next steps

1. `pka_triage/triage.py`: per parent, MolGpKa sites -> route decision
   (`decided` / `near_window` / `coupled`, with the nearest |pKa - pH|), and for
   decided parents the protomers written through the existing state model, in
   the shard schema `assign_names.py` reads. `routed.tsv` goes to the QupKake
   array.
2. A merge step so `llb names` sees one result per parent.
3. Harness on ~10k parents run both ways: site-set agreement, pKa scatter by
   site class, protomer-set agreement at each pH, and the fraction routed as a
   function of the routing threshold (curve, not one number).
