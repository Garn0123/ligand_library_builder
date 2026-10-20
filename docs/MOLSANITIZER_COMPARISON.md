# LLB vs MolSanitizer: what differs, what we checked, what to test next

2026-10-20. Against LLB `feat/pipeline` at 8c0cff3 (site filter v6, coupling v1,
aromatic-N v2) and MolSanitizer 0.7.4 (PyPI). Paper: Lam, Pach, Ullmann, …, Irwin,
Carlsson, Cabeza de Vaca, *MolSanitizer: an open-source pipeline to prepare
small-molecule databases for large-scale virtual screening*, ChemRxiv preprint,
9 Oct 2026, doi:10.26434/chemrxiv.15010189/v1 (not peer reviewed). Page numbers
below are the manuscript's own.

**Why this matters for us.** The docking scores train a generative model, so a
wrong protonation state or tautomer becomes a wrong training label. MolSanitizer
is the closest published, benchmarked pipeline to LLB (same DOCK db2 target, same
AMSOL), so a reviewer will ask how LLB compares to it.

## Bottom line

1. **The two pipelines get protonation states in different ways.** LLB predicts a
   pKa for every site with ML (MolGpKa for every parent; QupKake for the ~36% near
   the pH window). It then builds populations at pH 6.4, 7.4 and 8.4. MolSanitizer
   applies 144 hand-curated SMARTS rules (114 acid, 30 base), each with an
   *integer* pKa, and keeps one state per pH. Its CLI takes only integer pH.
   On our probe sets neither tool wins outright (table below). The paper's
   measured-pKa benchmarks (pp. 28–31) would place LLB on the same chart as
   MolSanitizer, Epik and ChemAxon. That is test **C**, the most useful test
   still to run.
2. **Tautomers are LLB's one real gap.** LLB moves protons only for exocyclic
   `=NH` and imidic acids (DECISIONS 2.2/2.3) and otherwise keeps ZINC's drawn
   tautomer. MolSanitizer has a two-layer tautomerizer and leads its TautoBase
   benchmark with 70% correct (p. 24). On 13 textbook probes it got 13/13 and
   LLB 5/13. What we don't know is **how often ZINC22 SMILES actually arrive
   in a minor tautomer**. Test **B** measures that on the drug-like 10k in
   minutes.
3. **The conformer and db2 stage is a possible future swap, not a defect.**
   MolSanitizer's own generator (ETKDG plus a TorLib torsion library) matched
   Conformator on bioactive recovery (Platinum set, 57% within 0.5 Å for both).
   It took a quarter of Conformator's time (1.51 h vs 6.02 h), and it also writes
   db2 with a C++ port of AMSOL 7.1. A laptop run built db2s from LLB-style
   names. It is not a drop-in for LLB's naming contract, though: its stereo
   suffix `NAME.1` is truncated to the last 16 characters in the db2 (shown
   below). Feeding it LLB's stereo-closed, contract-named `build.smi` with
   `--no-neutralize --no-stereoisomers` keeps the names intact.
4. **Fixes for LLB that came out of the probes** (each small; details under Test A):
   - (a) biguanides are built as +2;
   - (b) Meldrum's-type CH acids are dropped (a miss DECISIONS already knows about);
   - (c) piperazine keeps a minor +2 state at 7.4;
   - (d) DECISIONS 5.1 says ethylenediamine is 43% dication at pH 7.4 "from its
     pKas". 9.98/7.52 give **57%**.

## 1. What MolSanitizer does

Built on RDKit 2025.09.01 (p. 8). Steps run in this order: filters → tautomerizer →
ionizer → stereoisomers → conformers → output (Fig. 1, p. 8). Each step is a CLI flag
and also a Python class.

| step | what it does | defaults and numbers |
| --- | --- | --- |
| **Salt removal / standardize** | RDKit SaltRemover with a curated list of 65 counterions, then largest fragment by atom count (pp. 9–10). Neutralization runs before tautomers and protonation. `--standardize` is a separate ML-prep mode (`rdMolStandardize.Cleanup`, metal disconnection, neutralize, organometallics removed) that **ignores all other flags** (p. 21). Metal complexes: `--metal strict/soft/off` | neutralize on whenever tautomers or protonation are requested |
| **Filters** | 7 descriptors (HBA, HBD, MW, logP, heavy atoms, Fsp3, TPSA); PAINS (RDKit's 3 catalogs); 73 "unwanted" SMARTS in regular/optional/special tiers; custom SMARTS (pp. 10–11) | off unless asked |
| **Tautomers** (`--tautomers`) | Layer 1: RDKit `TautomerEnumerator`, then drop candidates that change "non-labile" groups (amide-like and ketone-like C, carboxylic acids, isolated C=C, methyl; phenol/alcohol O; sulfoximine S) or specified stereo. Rank with a modified Sitzmann score plus literature preferences (Table S1). **Keep the input unless something scores higher**; otherwise keep every candidate within 4 score units of the best, and avoid moving isolated C=C. Layer 2: literature-derived correction rules (File S3) that assign the preferred form and may emit several tautomers (pp. 11–13, Fig. 2). `--extended-tautomers`: up to 20 forms | `tautomers_v3.txt`, 357 lines |
| **Protonation** (`--protonation`) | SMARTS reaction rules, each with a fixed pKa and type (File S4; `ionizations_v3.txt`). An acid is deprotonated when pH > pKa, a base protonated when pH < pKa, both forms kept when pH = pKa. `--pH_range r` evaluates pH−r, pH, pH+r (each ±0.1) and keeps the union (pp. 13–14, Fig. 3). Rules encode context, e.g. no doubly protonated piperazine, no protonated amide, no phenolate without EWG; some moieties enumerate protomers (asymmetric piperazines). Boronates get OH⁻ | **pH 7, range 0**. The CLI rejects non-integer pH (`invalid int value: '7.4'`); the Python `Ionizer(pH=7.4)` accepts it. All rule pKas are integers |
| **Stereo** (`--stereoisomers`) | RDKit enumeration of unspecified centres and double bonds; isomers that fail 3D embedding are discarded. A protonated tertiary amine N is a *temporary pseudo-stereocentre*: both invertomers are sampled inside one entry's conformer set, not named separately (p. 15) | **max 8 isomers**; 60 s timeout; on by default with `--gen3d` |
| **Conformers** (`--gen3d`, C++ via pybind11) | srETKDGv3 embedding (50 tries for non-aromatic rings, 1 otherwise), MMFF94s with electrostatics off, ring-conformation library (38 ring systems, File S5), then the lowest-energy structure. Stochastic torsion sampling from a modified TorLib v3 (30° bins, "tolerance 2"). Bonds near the centre are picked more often. Clash check at 0.6×vdW, 25 kcal/mol window, 0.5 Å RMSD dedupe, then polar-H rotation (3 positions aliphatic, 2 aromatic). Fallback: Open Babel after a 2 min ETKDG timeout; CORINA optional (pp. 15–20) | **2,000 conformers**, `--nringconfs 1`, energy window 25, RMSD 0.5 |
| **Output** | DB2 (AMSOL 7.1 ported to C++, `db2writer`), PDBQT (Meeko), MOL2, SDF, SMILES; `msani_batch` for SLURM (pp. 9, 20, 23) | `db2.tgz` |

**Claims and benchmarks** (pp. 24–37):

| benchmark | set | metric | MolSanitizer | others |
| --- | --- | --- | --- | --- |
| tautomers | TautoBase curated subset, 873 aqueous equilibria; input = the *minor* tautomer | correct / incomplete / exceed / incorrect | **70% correct**, 20% incorrect | ChemAxon 48%, Epik Classic 23%, Unicon 12%, Scrubber 10%, Gypsum-DL 7% (Fig. 5) |
| protonation, simple | monoprotic, measured pKa: Reaxys + IUPAC digitized (Zheng 2026) + Lee 2008 | same, at pH 5–9 | **86%** | ChemAxon 86, Epik 85, Unicon 73, Scrubber 71, Gypsum 49 (Fig. 7A) |
| protonation, drug-like | AstraZeneca (ChEMBL3301361), Novartis (Liao & Nicklaus 2009), Vertex (Settimo 2014) | same | **77 / 78 / 69%** | ChemAxon 78 / 80 / 64; Epik 60 / 73 / 69 (Fig. 7B–D) |
| conformers | Platinum Diverse 2017, 2,859 bound ligands, ≤600 conformers | best heavy-atom RMSD | ≤0.5 Å **57%**, ≤2 Å **99.5%**, median 115 conformers, **1.51 h** total | Conformator 57% / 99.3% / 426 / 6.02 h; Conforge 60% / 98.7% / 72 / 2.13 h; RDKit 45% / 99.7% / 84 / 30.98 h; ZINC22 reference (CORINA + OMEGA) 48% / median 43 (Fig. 9) |
| docking | DUDE-Z, 43 targets, 2,312 actives, 69,904 decoys, DOCK3.8. **Conformers only**: their tautomers and protonation were *not* applied | adjusted logAUC | median 17, mean 22; better than reference on 33/43 targets | ZINC reference db2s: median 15, mean 18 (Fig. 10) |

Failures they report: clofazimine (missing iminophenazine rule); long-chain amines
put at pKa > 7 when measured < 6; phenols with uncounted substituent effects
(pp. 30–32); a layer-1 bias toward aromatic tautomers; macrocycles (≥9-membered
rings) not evaluated (p. 40).

**Getting it.** Apache 2.0. Code: github.com/carlssonlab/MolSanitizer. Docs:
msani.readthedocs.io. `pip install molsanitizer` installs prebuilt wheels:
macOS arm64/x86_64, manylinux_2_28 x86_64/aarch64, Windows; Python 3.10–3.14.
Dependencies: rdkit ≥ 2025.3.1, pandas, numpy, pyaml, platformdirs. The `pdbqt`
extra adds Meeko, scipy and gemmi; `oe` adds OpenEye. Open Babel is needed only
for its embedding fallback. There is a web server, carlssonlabtools.icm.uu.se/molsani,
with a 10,000-molecule limit. The benchmark inputs and every program's outputs
are Supplementary Files S7 (tautomers), S8 and S9 (protonation) on ChemRxiv;
they are not in the GitHub repo.

## 2. Stage-by-stage mapping

Legend: **=** the two overlap; **LLB+** only LLB does it; **MS+** only MolSanitizer
does it; **≠** both do it with different defaults.

| LLB stage | MolSanitizer counterpart | relation | the difference that matters |
| --- | --- | --- | --- |
| `parents`: largest fragment, `Cleanup`, `Uncharger`, dedupe by structure; ZINC-id validation and hashed ids | salt list (65) + largest fragment; Neutralizer | **≠** | LLB dedupes and validates ids (needed for the naming contract). MS strips named counterions first. On ZINC22 2D inputs (already single-component) these should agree; Test B checks |
| `parents`: imino → amino, imidic acid → amide (DECISIONS 2.2, 2.3) | Tautomerizer layers 1 + 2 | **≠ / MS+** | both fix `N=c1[nH]…` and `RO-C(=NH)OH` (probes: 4/4 each). MS also handles lactam/lactim, thione/thiol, keto/enol, guanine N-H position, dearomatized forms. LLB keeps ZINC's tautomer for all of these by design (DECISIONS 2.2 rejected a general canonical tautomer). MS can return *several* tautomers per parent |
| — | Filters (PAINS, unwanted, descriptors) | **MS+** | LLB filters nothing beyond ZINC's own tranche selection. For a training set that is a choice, not a gap: PAINS are not docking-wrong |
| `triage`: MolGpKa per-site pKa for every parent; route to QupKake if a site is within 1.5 of a target pH or same-kind sites ≤3 bonds apart | Ionizer: rule pKa per functional group | **≠** | continuous ML pKa per site versus integer pKa per SMARTS class. LLB is whole-molecule aware (a GCN over the molecule); MS knows only the environment its SMARTS spells out. MS has nothing like routing (it doesn't need it) |
| QupKake on routed parents (xtb + GNN) | — | **LLB+** | about 32 s/parent and 320 CPU-h per 100k parents. MS has no QM step |
| site filter v6 (amide N-H, amine N-H, O base, alcohol O-H, C acid except 1,3-dicarbonyl/nitro, one site per amidine/guanidine, S/azide/diazo) | absent by construction: a rule only exists for groups meant to ionize | **=** in outcome | both leave amides, sulfides and azides alone (probes: all agree). LLB needs the filter because ML predictors propose these sites; MS never proposes them. Exception: Meldrum's acid, which MS ionizes (rule) and LLB's `carbon_acid` drops |
| state model: Henderson–Hasselbalch per site, margin 1.0, min population 0.01, ≤8 states, pH 6.4/7.4/8.4 | apply rules at pH ± range | **≠** | LLB keeps every state above 1% within ±1 unit of pH (median ~1.9 names/parent including stereo). MS keeps one state per pH value (both forms only for rules flagged `Enumerate`, or when pH equals a rule's pKa): compact but "incomplete" when a site sits near the pH (imidazole, 2-aminopyridine, phthalimide). Default pH 7 vs LLB's 7.4 |
| coupling model v1 (bonds-apart penalty from diamine/diacid micro-pKa shifts) | encoded inside individual SMARTS (no doubly protonated piperazine; biguanide singly protonated) | **≠** | LLB is general but approximate. Its known limit, conjugated paths under-corrected, shows up as metformin +2. MS is exact where a rule exists and blind where none does |
| aromatic-N hybrid v2 (MolGpKa's ring-N value in routed parents with an off-ring base) | pyridine/imidazole rules with fixed pKas (imidazole 6, 2-aminopyridine 6, 4-aminopyridine 9, pyridine-3-amine 6) | **LLB+** in mechanism | the hybrid fixes a QupKake artefact that MS doesn't have. On the tethered reference molecules both keep the ring N neutral beside a charged base (all 24 agree) |
| `merge` (QupKake wins for routed; fallback to MolGpKa) | — | **LLB+** | bookkeeping for two predictors |
| `names`: stereo closure under db2_converter's options, ≤32 isomers, 16-char contract names, invertomers named separately (`invertomer_of`) | RDKit enumeration, ≤8 isomers, `NAME.1/.2` suffixes; invertomers sampled inside one entry | **≠** | (1) MS's suffix breaks 16-char names: `Ch300000BlLjl00C.1` → `M 300000BlLjl00C.1` in the db2 (run 2026-10-20). (2) MS's 8-isomer cap is below LLB's 32. (3) MS's in-entry invertomer sampling is the "collapse before building" alternative DECISIONS 10.1 deferred: one name, both N configurations, the docking minimum taken automatically |
| 3D stereo check: embed isomers with a centre on a multi-ring atom, 6-try escalation, last-isomer guard | drop stereoisomers RDKit cannot embed | **≠** | same idea. LLB escalates and never drops a protomer's last isomer (DECISIONS 7.3–7.4: bicyclo[1.1.1]pentanes); MS documents neither safeguard |
| db2: db2_converter (hnlab 63d6656), Conformator, 600 conformers, restrained MMFF (ε = 4r, 0.3 Å), `--reseth --rotateh`, AMSOL 7.1, mol2db2 | confgen (ETKDG + ring library + TorLib), MMFF94s without electrostatics, 2,000 conformers, AMSOL 7.1 C++ port, db2writer | **≠** | different generators with similar recovery in the paper (Conformator 57%, MS 57% at 0.5 Å). MS is about 4× cheaper on Platinum, and LLB adds restrained MMFF on top of Conformator. LLB's measured cost: median 45.8 s/name. MS on the laptop: 3 molecules → 4 db2 in 12 s wall on 3 cores (not comparable hardware) |
| `verify-db2`: AMSOL charge vs requested formal charge, name and rename checks, faillist accounting | none described | **LLB+** | LLB's 0-charge-mismatch guarantee has no MS counterpart |
| `mol2` / `mol2-library` from the db2 | MOL2/SDF/PDBQT writers | **≠** | MS also writes PDBQT |

## 3. What already ran (laptop, 2026-10-20, a few seconds of CPU)

Setup: MolSanitizer 0.7.4 in a scratch venv (Python 3.14.6, RDKit 2026.9.1,
newer than the paper's tested 2026.3.6). LLB's `parents` + `triage` at 8c0cff3 in a
scratch `llb_pka` env (RDKit 2024.03.5, torch 2.4.1, MolGpKa 4dc8352). **QupKake
was not run** (it needs xtb), so for routed parents the LLB column is its
MolGpKa-path state, which production replaces with QupKake's. Exact versions:
`docs/molsanitizer_probes/results/versions.txt`.

Scoring uses the paper's classes over the set of net charges at pH 7.4.
**correct** means the same set. **incomplete** means a strict subset.
**exceed** means a strict superset. **incorrect** means the expected charge is
missing. "Expected" holds every charge with ≥10% population from the measured
pKa. These sets were built from LLB's own decision list, so they are a
rule-by-rule check, not a benchmark. Both tools were plausibly tuned on these
textbook classes.

### Test A1: LLB rule probes (35 molecules, `llb_rule_probes.tsv`)

| | correct | incomplete | exceed | incorrect |
| --- | --- | --- | --- | --- |
| LLB, MolGpKa path | 29 | 2 | 3 | 1 |
| MolSanitizer pH 7 (CLI default) | 30 | 5 | 0 | 0 |
| MolSanitizer pH 7.4 ± 1 (API) | 31 | 4 | 0 | 0 |

LLB misses, with what production would do:
- **metformin**: LLB gives +1,+2 at 7.4, with the +2 dominant (MolGpKa 9.57 and 9.18
  on the two =N; the 4-bond coupling penalty of 1.1 is too small). Measured
  13.25/3.07 → +1. The parent is routed (near_window), so production uses
  QupKake's sites, which are unchecked. MS: +1 (its `biguanides` rule).
  **Candidate site filter v7: one basic site per biguanide**, the same reasoning as v5.
- **piperazine**: +1,+2 (coupled MolGpKa model: 27% dication at 7.4; measured
  2%). The extra +2 is docked as a separate state. MS: +1.
- **4-aminopyridine**: 0,+1 (MolGpKa 6.98; measured 9.17). Routed; QupKake has 9.49,
  so production is right. Residual risk: the aromatic-N hybrid would hand a
  4-aminopyridine ring N MolGpKa's 6.98 when the molecule also has an off-ring
  base. MS: +1.
- **Meldrum's acid**: 0, because `carbon_acid` drops the site (measured 4.97;
  DECISIONS 4.4 lists this as a known miss). MS: −1.
- **ethylenediamine**: +2 only (MolGpKa 10.63 on both N; real 57/43% +2/+1).
  MS: +1 only. Both incomplete, from opposite sides.
- **phthalimide**: neutral (MolGpKa 10.10; measured 8.30, 11% anion). MS: neutral.

MS misses: DABCO (keeps +1, measured 14% neutral), 2-aminopyridine (rule pKa 6 vs
6.86 measured; misses the 22% cation), p-nitrophenol at pH 7 (rule pKa 8 vs 7.15
measured, so the dominant anion is missing unless range 1 is used), and
ethylenediamine and phthalimide, as above.

Representation note: deprotonated dimedone comes out as `[CH-]` from LLB and as
the enolate `C=C([O-])` from MS. RDKit perceives LLB's carbanion as sp2, and MMFF
keeps it planar (checked). Conformator's handling is unchecked.

### Test A2: aromatic-N reference set (`pka_triage/reference/aromatic_n_bases`, 53 molecules)

20 single-ring molecules have a measured pKa:

| | correct | incomplete | exceed |
| --- | --- | --- | --- |
| LLB, MolGpKa path | 18 | 0 | 2 (4-aminopyridine; benzotriazole −1,0) |
| LLB, production estimate (QupKake's ring-N value for routed rows, LLB margin rule)¹ | 20 | 0 | 0 |
| MolSanitizer pH 7.4, range 0 or 1 | 16 | 4 | 0 |

¹ The estimate models only the ring-N base. For benzotriazole, production would also
carry QupKake's N-H acid site, so it probably repeats the MolGpKa path's −1,0.

MS's four misses are imidazole (6.95), 2-methylimidazole (7.75), 4-methylimidazole
(7.45) and 2-aminopyridine (6.86). It leaves all four neutral at 7.4 ± 1 because
its rules put both classes at pKa 6. For 2-methylimidazole that loses the
*dominant* form (69% cation). With the integer CLI at pH 7 ± 1 (pH 6, 7, 8),
pH 6 equals the rule pKa, so both forms come back. The benzotriazole "exceed" is
LLB keeping the N-H anion; benzotriazole's N-H acidity (~8.4, not in the
reference set) puts about 10% anion at 7.4, so it is defensible. For the 24
tethered molecules (ring + guanidine or dimethylamine), MS and LLB's
production-path states agree: the base is charged and the ring N neutral.

### Test A3: tautomer probes (13 molecules, `tautomer_probes.tsv`; input = the minor form)

| | correct | incorrect |
| --- | --- | --- |
| LLB `parents` (standardize + `amino_tautomer`) | 5 | 8 |
| MolSanitizer `Tautomerizer()` | 13 | 0 |

LLB gets exactly its four `=NH` and imidic-acid cases plus cytosine, and keeps
the input everywhere else: 2-/4-hydroxypyridine, 2-quinolinol, uracil dienol,
2-mercaptopyridine, acyclovir 3H, and the cyclohexadienone and catechol keto
forms. That is the expected behaviour by design. The open question is the
frequency of such inputs in ZINC22, which Test B answers.

## 4. Test plan, ranked by value per cost

| rank | test | question | input | metric → what would change LLB | cost |
| --- | --- | --- | --- | --- | --- |
| — | **A** (done) | do LLB's rule-level decisions match an independent curated rule set? | the probe sets above | see section 3 → v7 biguanide filter; reconsider Meldrum; fix DECISIONS 5.1 (57%) | seconds, laptop |
| 1 | **B. 10k diff** | how often do LLB and MS disagree on real ZINC22 parents, in which classes, and how often does MS change the tautomer? | drug-like 10k: `druglike_10k_v2` `parents/parents.tsv` + `library_v6b/library.tsv` (on scratch, purged ~60 days after last use) | `agree_dominant` (MS's state is among LLB's at 7.4), `agree_set_r1`, `ms_tautomer_changed`, per class; then a person judges the seeded 100 in `review.tsv` against Williams/Reich values. **Change LLB if** a class with ≥0.5% of parents disagrees and review sides with MS ≥2:1 (add a site-filter rule or a fixed-pKa override). **Add a tautomer layer** if MS changes ≥1% of parents and review sides with MS ≥2:1 on the non-trivial changes (ignore amidine/guanidine double-bond placement, which is the same species) | ~0.03 s/parent on 1 core → **<0.2 CPU-h**; laptop once the two TSVs are copied, or one login-node or `srun` core |
| 2 | **C. measured pKa** | where does LLB sit on the paper's Fig. 7 chart (monoprotic; AZ/NV/VX drug-like)? | ChemRxiv SI Files **S8, S9** (inputs + every program's outputs, so LLB is scored against identical reference states). Fallback without the SI: AZ ChEMBL3301361, Liao & Nicklaus 2009 SI, Settimo 2014 SI, Zheng 2026 IUPAC set, Lee 2008 SI (Reaxys part not public) | the paper's 4-class outcome at its pH convention. **Change LLB if** it is below MS by ≥5 points on a drug-like set, *and* the misses cluster in a class we can fix. Also report LLB's "exceed" rate, which is the price of its population model | per 1,000 molecules: MolGpKa ~10 s, MS ~30 s, QupKake on the ~36% routed ≈ **3.2 CPU-h** → Great Lakes array (1 core/task), everything else on the laptop |
| 3 | **D. TautoBase** | LLB's tautomer accuracy on the paper's 873 entries | SI File **S7** (or TautoBase, github.com/WahlOliver/tautobase, with the curation in their Supplementary Methods) | 4-class outcome. The expected result is "low, because LLB keeps the input", so the value is the number for a methods section and the list of classes LLB misses. Feeds the decision from B | seconds, laptop |
| 4 | **E. docking sensitivity** | when the states or tautomers differ, do the **docking scores** differ enough to matter for the training data? | the B disagreements (or a seeded 1,000 of them), each built twice through *LLB's* db2 path (LLB states vs MS states/tautomers), so only the chemistry differs; 30 pockets, one per family, as in the docking pilot | per parent and pocket, Δ(min score over the parent's names); Spearman across pockets. **Change LLB if** a class shows median \|Δ\| > 1 kcal/mol or ρ < 0.9 *and* B/C say MS is right for that class | per 1,000 parents: db2 ≈ 2 × 1.9 names × 52 s ≈ 55 CPU-h + docking ≈ 2 × 1,900 × 30 × 1.28 s ≈ 40 CPU-h → **~100 CPU-h**, Great Lakes. Needs a ~50-line adapter that writes MS states as a merge-style shard so `llb names` closes stereo and assigns names |
| 5 | **F. MS db2 at scale** | could MS's confgen+db2 replace db2_converter for H23–30 and later: build rate, cost, name and charge integrity? | `library/build.smi` of the drug-like 10k (19,102 names), `msani -i build.smi --no-neutralize --no-stereoisomers --gen3d -f db2 --numconfs 600` | built fraction vs LLB's 18,620/19,102; CPU-h vs LLB's 268; `llb verify-db2` on MS output (charges, names); DOCK6 reads it (`db2tool audit`). **Switch** only if ≥ LLB's build rate, 0 charge mismatches, names intact, and E-style docking on 30 pockets shows no score shift | MS maybe ~3–5 s/name (laptop figure, unverified on GL) → **~20–30 CPU-h**; `msani_batch` or an LLB array |
| 6 | **G. Platinum conformers** | does LLB's actual db2 conformer set (Conformator + restrained MMFF + mol2db2 hierarchy) recover bioactive poses as well as raw Conformator and MS? | Platinum Diverse 2017, 2,859 ligands (public) | best heavy-atom RMSD at ≤0.5/1/2 Å; conformer count; CPU-h. **Revisit restrained MMFF** if LLB's ≤0.5 Å recovery drops >5 points below the paper's raw Conformator figure (57%) | LLB ≈ 2,859 × 46 s ≈ **37 CPU-h** (GL); MS ≈ 1.5 h on their hardware |
| 7 | **H. DUDE-Z enrichment** | end to end, LLB vs MS vs the ZINC reference db2s on 43 protein targets | DUDE-Z (dudez.docking.org): SMILES, grids, reference db2s | adjusted logAUC per target. Low relevance: our targets are RNA pockets and the use is training labels, not enrichment | ~72k molecules × ~1.9 names × 52 s ≈ **~2,000 CPU-h** for LLB + DOCK3.8 runs; run on 5–8 targets first if at all |

Order of work: **B now** (copy the two TSVs before scratch purges). C and D as
soon as the SI files are downloaded. E only for the classes B and C flag. F
when planning the next library after H23–30. G and H only if a reviewer asks.

## 5. Commands

### Install MolSanitizer (no system-wide changes)

Laptop (what was used here):

```bash
python3 -m venv ~/envs/msani && ~/envs/msani/bin/pip install molsanitizer==0.7.4
~/envs/msani/bin/msani --version
```

Great Lakes. The wheels are manylinux_2_28, so check `ldd --version` ≥ 2.28 first:

```bash
conda create -y -p /nfs/turbo/umms-maom/johnbick/opt/envs/msani -c conda-forge python=3.12
/nfs/turbo/umms-maom/johnbick/opt/envs/msani/bin/pip install molsanitizer==0.7.4
```

### Test A (reproduce section 3)

From the repo root, with an `llb_pka` env (`PKA_ENV`) and a MolGpKa checkout at
4dc8352 (`MOLGPKA_DIR`):

```bash
P=docs/molsanitizer_probes; OUT=/tmp/msani_probes; mkdir -p $OUT
tail -n +2 $P/llb_rule_probes.tsv | awk -F'\t' '{print $2" "$1}' > $OUT/probes.smi
llb parents $OUT/probes.smi -o $OUT/probes.parents
llb triage $OUT/probes.parents -o $OUT/probes.triage
~/envs/msani/bin/python $P/compare_probes.py $P/llb_rule_probes.tsv $OUT/probes.parents $OUT/probes.triage -o $OUT/probes.compare.tsv
```

```bash
llb parents pka_triage/reference/aromatic_n_bases.smi -o $OUT/aromn.parents
llb triage $OUT/aromn.parents -o $OUT/aromn.triage
git show data/validation-2026-10:validation/2026-10_druglike10k/druglike_10k/ref_qk/reference_check.tsv > $OUT/reference_check.tsv
~/envs/msani/bin/python $P/compare_aromatic_n.py $OUT/reference_check.tsv $OUT/aromn.parents $OUT/aromn.triage -o $OUT/aromn.compare.tsv
```

```bash
llb env PKA python $P/compare_tautomers.py llb $P/tautomer_probes.tsv > $OUT/t_llb.tsv
~/envs/msani/bin/python $P/compare_tautomers.py msani $P/tautomer_probes.tsv > $OUT/t_ms.tsv
~/envs/msani/bin/python $P/compare_tautomers.py score $P/tautomer_probes.tsv $OUT/t_llb.tsv $OUT/t_ms.tsv
```

The outputs from 2026-10-20 are in `docs/molsanitizer_probes/results/`.

### Test B (the 10k diff; <15 min)

First copy the inputs off scratch (they purge). Check the exact subdirectory names
in the run directory first:

```bash
R=/scratch/maom_root/maom99/johnbick/2026_10_08_build_small_mols/druglike_10k_v2
mkdir -p ~/msani_b && cp $R/parents/parents.tsv ~/msani_b/parents.tsv && cp $R/library_v6b/library.tsv ~/msani_b/library.tsv
```

Then, anywhere with the msani env:

```bash
python docs/molsanitizer_probes/compare_library.py ~/msani_b/parents.tsv ~/msani_b/library.tsv -o ~/msani_b/out -j 4
```

It writes `out/summary.json` (overall and per-class agreement, tautomer-change rate,
seconds), `out/parents.tsv`, and `out/review.tsv`: 100 seeded disagreements with
empty `verdict` / `evidence` columns for a person to fill. The same command on the
H18–22 production library (`parents/parents.tsv` + `library/library.tsv` of that
run) gives the 50k version.

### Test C and D inputs

Download the Supplementary Information from the ChemRxiv page
(doi:10.26434/chemrxiv.15010189/v1, "Supplementary materials"). S7 holds the
TautoBase inputs and outputs; S8 and S9 the protonation sets. Then:
- run `llb parents` + `llb triage` on each input set;
- `llb submit-qupkake triage/routed …` on Great Lakes;
- `llb merge`, then score LLB's 7.4 states with the paper's classes against the
  reference states in S8/S9.
`compare_probes.py` shows the scoring; it needs a reader for the SI's column
layout, which is unknown until the files are in hand.

## 6. Smaller findings to file

- DECISIONS 5.1 says ethylenediamine is "43% from its pKas" at pH 7.4. From
  9.98/7.52 the dication fraction is 1/(1 + 10^(7.4 − 7.52)) = **57%**; the line
  looks swapped. The model's 37% is then 20 points low, not 6.
- MolSanitizer's CLI rejects `--pH 7.4` (integers only); 7.4 needs the Python API
  (`Ionizer(pH=7.4)`). Record which was used in any comparison: pH 7 and pH 7.4
  give different states for rules whose pKa is 7 or 8.
- If MS's db2 path is ever used, run it on stereo-closed, contract-named input
  with `--no-neutralize --no-stereoisomers`. Otherwise the `.N` suffix is
  truncated (`M 300000BlLjl00C.1`), as db2_converter's `--checkstereo` was
  before the closure step.
- MS's in-entry invertomer sampling gives DECISIONS 10.1's deferred
  alternative, "collapse before building", for free. If the docking run's
  `non_first_won_rate` comes out small, that is the design to copy: one name,
  both N configurations in one db2.
