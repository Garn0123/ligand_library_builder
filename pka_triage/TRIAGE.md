# pKa triage: MolGpKa first, QupKake only where it matters

Branch `feat/pka-triage`. Status 2026-10-08: environment, MolGpKa port, triage
and merge written and tested end to end with MOCK QupKake output (no real
QupKake or Great Lakes run yet). Validation harness not written yet.

## Running it

```
llb parents my.smi -o parents                         # as before
llb triage parents -o triage                          # MolGpKa, ~5 ms/parent, one core
llb submit-qupkake triage/routed protomers_qk --time=01:30:00
llb merge --parents parents --triage triage --qupkake protomers_qk -o merged
llb names --shards-dir merged --parents parents -o library
```

`triage/routed/` is an ordinary parents directory (shards of
`--qupkake-shard-size`, default 10, plus `preflight.smi` for timing). For a big
library, run slices `--part K/N` as separate jobs, then `--finish`.
`triage/stats.json` gives the routed fraction and what it would have been at
other `--route-window` values. `library.tsv` has a `site_source` column:
molgpka, qupkake or molgpka_fallback. The fallback is used when QupKake failed
on a routed parent; `--on-qupkake-fail drop` drops those parents instead.

Settings: `--ph` (6.4 7.4 8.4), `--route-window` (1.5), `--margin` (1.0, same
as the QupKake path), `--coupling-bonds` (3), `--coupled route|ignore`,
`--min-population`, `--max-states`, `--max-ambiguous`, `--coupling on|off`
(on; see below). They are recorded in
`triage/settings.json`; a second run with different settings into the same
`-o` is refused. Merge refuses QupKake shards run with a different pH list,
margin, min population, max states, site filter or coupling.

Also filtered, in BOTH paths (`qupkake_protomers.IMPLAUSIBLE_SITES`,
`--keep-amide-sites` turns it off):
- amide-type N as a base (MolGpKa: anilides 3.9-6.0, benzamide 3.4; real ~ -1);
- amide-type N-H as an acid (amides, anilides, ureas, carbamates, hydrazides,
  thioamides; real ~13-17). On Owen's 50 QupKake put 44 of these at 2.4-6.7,
  so the all-QupKake library had an amide ANION as the dominant state at
  pH 7.4 in 22 of 50 parents; all 5 "settled but disagree" parents were this.
  Imide N-H and N-sulfonyl N-H are kept (genuinely acidic);
- the amino N's of an amidine or guanidine as bases (`amidine_amino_n_base`,
  site filter v5). One proton per group, on the =N, charge shared: guanidine
  13.71, acetamidine 12.52, benzamidine 11.6, phenylguanidine 10.88, one pKa
  each; biguanide 13.25 then 3.07 (Williams p. 29). QupKake gave each amino N
  its own site at 10.7-11.1 beside the =N, building +2/+3 guanidines (and the
  =NH protomers that crash RDKit's stereo enumeration). If the =N got no
  site, the best amino-N site moves onto it (`filter_sites`), so the
  protonation is never lost; triage `sites.tsv` marks it
  `moved_from_amidine_amino_n`. Aromatic C is not matched (2-aminopyridine:
  the ring N is the base);
- carbon as an acid (`carbon_acid`, site filter v4), except a CH between two
  ketone/aldehyde carbonyls or a CH carrying a nitro group. MolGpKa put
  malonic acid's CH2 at 7.5 and built a [CH-] trianion.
Sources for the amide numbers (checked 2026-10-08):

| value | compound | solvent | source |
|---|---|---|---|
| 15.1 (the table marks water values > 14 as extrapolated) | acetamide N-H | water | Ripin & Evans, pKa table (Harvard Chem 206), "Amides & carbamates" |
| 25.5 / 25.9 / 23.3 / 21.5 | acetamide / N-methylacetamide / benzamide / acetanilide N-H | DMSO | Reich, *Equilibrium pKa Table (DMSO)*, ACS Division of Organic Chemistry, updated 2024-07-26, "Amides (NH)" (Bordwell data) |
| 26.9 / 24.2 / 18.5 | urea / ethyl carbamate / thioacetamide N-H | DMSO | same table, "Amides (NH)", "Carbamates (NH)" |
| 14.7 | succinimide N-H | DMSO | same table, "Imides (NH)" |
| 8.30 | phthalimide N-H | water | Ripin & Evans, "Imides" |
| 38 (extrapolated) / 41 | ammonia N-H | water / DMSO | Ripin & Evans, "Amines" |
| 44 / 30.6 / 25.0 | pyrrolidine / aniline / diphenylamine N-H | DMSO | Reich DMSO table, "Amines (NH)" |
| -12.4 / -7.8 / -3.8 / -2.2 | protonated nitrobenzene / benzoic acid / dimethyl ether / methanol | water | Reich, *Equilibrium pKa Table (H2O)*, "Oxygen Acids" |
| 16 / 20 (approximate) / 10.0 | methanol / tert-butanol / phenol O-H | water | same page |
| -0.62 +/- 0.07, -0.42, -0.28 | protonated acetamide, N-methylacetamide, N,N-dimethylacetamide | water, 298 K | Grant, McTigue & Ward, *Aust. J. Chem.* 1983, 36, 2211-2218, doi:10.1071/CH9832211 |

Carbon acids, water, from R. Williams, *pKa Data Compiled by R. Williams*
(ACS Organic Division, updated 4/7/2022), "Carbon Acids", pp. 31-32, and
dimedone p. 29 (checked 2026-10-15):

| value | compound | kept? |
|---|---|---|
| 8.95 / 9 | acetylacetone | kept (1,3-diketone) |
| 8.23 (enol), 9.4 | benzoylacetone | kept |
| 5.23 / 5.85 / 5 | dimedone / triacetylmethane / malondialdehyde | kept |
| 4.7 / 6.82 | CF3COCH2COCH3 / C6H5COCH2COCF3 | kept |
| 7.74 / 8.6 / 10.29 | 2-nitropropane / nitroethane / nitromethane | kept (nitro) |
| 10.68 / 10 | ethyl / methyl acetoacetate | dropped |
| 9 / 11.81 | ethyl cyanoacetate / malononitrile | dropped |
| 13.3 | diethyl malonate | dropped |
| ~20 / ~25 / ~23 | acetone / acetonitrile / dimethyl sulfone | dropped |

The `carbon_acid` row: only the 1,3-dicarbonyl (ketone/aldehyde) and nitro
classes reach the 6.4-8.4 window; ester-, nitrile- and sulfone-activated CH
sit at 9 or above. Known miss: Meldrum's acid (cyclic malonate diester, water
pKa about 4.8; not in this compilation) is dropped by the rule, as is any
other doubly ester-activated CH that is unusually acidic. Enol forms are O-H
sites and are not affected.

The oxygen rows back `oxygen_base` and `alcohol_oh_acid` (site filter v3,
2026-10-14): Owen's rebuild had C(=O)[OH2+], a protonated nitro group and an
alkoxide among QupKake's protomers. Phenols and enols stay.

The amine rows back `amine_nh_acid` (site filter v2, 2026-10-08): Owen's set
produced [NH-]/[N-] amine-anion protomers. (The three that then failed
db2_converter with `1generate` belong to a parent that fails the same way with
sensible protomers too, so that failure is the scaffold's, not the filter's.)

What these support: no listed value puts an ordinary amide, anilide, urea or
carbamate N-H anywhere near 2-7 in water. The one water value is 15.1
(acetamide, extrapolated), and in DMSO every one of them sits 4-12 units above
an imide, whose water pKa is about 8 (phthalimide). That supports dropping
them. It does NOT give a water pKa for anilides or ureas specifically, and the
acylsulfonamide "keep" is chemistry knowledge, not a number checked here.

Dropped sites are listed per shard in `dropped_sites.tsv` (QupKake path) and
`decisions.tsv` / `sites.tsv` kept=0 (triage). QupKake runs enumerated before
this filter are rebuilt from their SDFs with `llb reenumerate` (no xtb);
`llb merge` refuses them until then.

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

1. **Neither predictor models site-site coupling; the state model now does
   (`--coupling on`, the default since 2026-10-15).** QupKake and MolGpKa each
   predict every site with the others neutral, so the independent
   Henderson-Hasselbalch product charges every amine whose pKa clears the pH.
   On the drug-like 10k (before the embed fix), about 95 of the 165
   settled-but-disagree parents were MolGpKa +1 against QupKake +2 to +5, at
   every pH. MolGpKa's own piperazine
   (10.46 on both N) locked the dication the same way; experimentally the
   second pKa is 5.68 and the monocation dominates at 7.4.

   The fix (`qupkake_protomers.Coupling`) scores every combination of sites
   with a penalty for each two same-kind sites charged together, by bonds
   apart: micro-pKa shifts from symmetric diamines and diacids, (pKa1 - pKa2)
   - log10(4), all from the Williams compilation (ACS Organic Division, rev.
   4/7/2022; values and pages in the code comment):

   | bonds apart | bases | acids | from |
   |---|---|---|---|
   | <= 3, same ring | 3.5 | 3.5 | piperazine 9.82/5.68, DABCO 8.19/4.18 |
   | <= 3 | 2.2 | 2.4 | ethylenediamine 9.98/7.52, 1,2-diaminopropane, 2,3-diaminobutane; oxalic 1.23/4.19 |
   | 4 | 1.1 | 2.3 | 1,3-diaminopropanes 1.4 and 0.85; malonic 2.83/5.69 |
   | 5 | 0.9 | 0.7 | putrescine 10.84/9.31; succinic 4.19/5.48 |
   | 6 | 0.7 | 0.5 | cadaverine 11.05/9.74; glutaric 4.34/5.42 |
   | 7 | 0.4 | 0.4 | adipic 4.42/5.41 (bases borrow it) |
   | 8-10 | 0.3 | 0.3 | extrapolated |

   Checked against the same tables: ethylenediamine at pH 7.4 is 37% dication
   in the model, 43% from its pKas; piperazine at 6.4 is 17% against 16%.
   Limits: opposite charges are not coupled (the correction only removes
   charge); conjugated paths (pyrazine, pyrimidine, amidines) couple more
   strongly than these aliphatic references, so they are under-corrected; the
   predicted pKas themselves are unchanged (MolGpKa's 10.46 for piperazine
   still makes the dication 65% at 6.4). `--coupling off` restores the
   independent product exactly. Merge refuses halves built with different
   coupling; `llb reenumerate` rebuilds QupKake shards without xtb.

   States are also merged by SMILES now: in a symmetric molecule protonating
   either N is one species, and was written twice before (both modes).

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

1. First real run: Owen's 50, `llb triage` then QupKake on `triage/routed`,
   to check that QupKake's real shard outputs merge cleanly.
2. `llb compare` (written 2026-10-08) is the harness: per-site pKa by site
   class, dominant microstate agreement per pH, and the routed fraction and
   miss rate at every window. First real use: Owen's 50. Then ~10k parents run both ways (`--route-window 99` routes
   everything, so triage's protomers.csv and sites.tsv sit next to QupKake's
   for every parent): site-set agreement, pKa scatter by site class,
   protomer-set agreement per pH, routed fraction against the window.
3. A SLURM wrapper for `triage --part` once a library is big enough to need it
   (1M parents is about 1.5 CPU-hours).
