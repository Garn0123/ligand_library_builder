# Decision record: how ligand_library_builder chooses what to build

> **Code at tag `llb-full-2026-10`.** The construction code this record cites
> (`pka_triage/`, `run_qupkake/`, `pipeline/`, `environments/`) was removed on
> 2026-10-26, when ligand construction moved to MolSanitizer. Paths and commits
> here refer to that tag (`git checkout llb-full-2026-10`). Kept in the tree
> because `lib_h18_22` was built and docked with that code.

Every choice that decides which molecules, protonation states, tautomers and
stereoisomers end up in a library, with the setting it governs, the commit that
made it, the evidence behind it, and what was tried or considered instead.
Written for a paper's methods section and for answering a reviewer.

This file is canonical. `pka_triage/TRIAGE.md` keeps the longer tables and the
operational detail; code comments repeat the short form next to each constant.
If they disagree, this file and the code are what to check, and the
disagreement is a bug to fix.

Draft of 2026-10-16, against `main` at `915b93c` (PR #1 merge); evidence
files from branch `data/validation-2026-10` (`5b64a1f`).

## How to read an entry

| field | meaning |
| --- | --- |
| **Decision** | what the pipeline does |
| **Governs** | the setting, constant or version string it sets, and where it lives |
| **Commits** | where it came in or changed (`git show <hash>`) |
| **Evidence** | a citation (author, title, table, page) or a measured number with the run it came from |
| **Alternatives** | what else was considered, and why not |
| **Status** | `sourced` (rests on a cited value), `measured` (rests on a number from our runs), `judgment` (rests on chemistry knowledge or a design argument, not on a checked number), or a combination |

A `judgment` status is not a weakness to hide. It marks where a reviewer can
fairly ask for more, and where a future measurement would replace the argument.

## Sources cited

| short name | full reference |
| --- | --- |
| **Williams** | R. Williams, *pKa Data Compiled by R. Williams*, ACS Division of Organic Chemistry, updated 4/7/2022. organicchemistrydata.org/hansreich/resources/pka/pka_data/pka-compilation-williams.pdf. Water, 25 °C unless marked. |
| **Reich DMSO** | H. J. Reich, *Equilibrium pKa Table (DMSO Solvent and Reference)*, ACS Division of Organic Chemistry, updated 2024-07-26 (Bordwell data). organicchemistrydata.org/hansreich/resources/pka/pka_data/pka-compilation-reich-bordwell.pdf |
| **Reich H2O** | H. J. Reich, *Equilibrium pKa Table (H2O)*, p. 9 of the Reich–Bordwell PDF above ("Oxygen Acids" and protonated species; several entries marked approximate) |
| **Ripin & Evans** | D. H. Ripin and D. A. Evans, pKa table (Harvard Chem 206), updated 4/16/2024, evans_pKa_table.pdf (scanned; read as page images) |
| **Grant 1983** | H. M. Grant, P. McTigue and D. G. Ward, *Aust. J. Chem.* 1983, 36, 2211–2218, doi:10.1071/CH9832211 |
| **MolGpKa** | X. Pan, H. Wang, C. Li, J. Z. H. Zhang and C. Ji, *J. Chem. Inf. Model.* 2021, 61, 3159–3165; github.com/Xundrug/MolGpKa |
| **QupKake** | O. Abarbanel and G. R. Hutchison, *J. Chem. Theory Comput.* 2024, 20, 6946 |
| **Angyal 1952** | S. J. Angyal and C. L. Angyal, *J. Chem. Soc.* 1952, 1461; S. J. Angyal and W. K. Werner, *J. Chem. Soc.* 1952, 2911 |
| **Primary sources from the sourcing pass** | Lassalas et al. 2016; Pihlaja & Seilo 1969; Pinching & Bates 1948, 1950; Katritzky, Karelson & Harris 1991; Forlani & De Maria 1982; Bonvicini et al. 1972; McGarrity & Smyth 1977: full references in section 0 |
| **Li 1994** | K.-H. Li, *ACM Trans. Math. Softw.* 1994, 20, 481–493 (reservoir sampling, Algorithm L) |

Section 0 lists the sourcing pass run for this record (2026-10-16) and what it
found for the decisions that rested on chemistry knowledge.

## Runs cited

Measured numbers name the run they came from. The small summary files from the
two validation runs are in this repo under `validation/2026-10_druglike10k/` (branch
`data/validation-2026-10`, commit `5b64a1f`; its README gives the run and code
state behind each file). Paths below are relative to the repo root. Numbers
marked *read off the terminal* were read from a run's console output at the
time and are not in those files; they can be recomputed from the run
directories (`/scratch/maom_root/maom99/johnbick/2026_10_08_build_small_mols/druglike_10k{,_v2}`
on Great Lakes) until scratch purges them.

| run | what it is | evidence files |
| --- | --- | --- |
| **Owen's 50** | 50 large, flexible, peptide-like molecules from a collaborator (up to 39 heavy atoms), non-ZINC names; first real QupKake + triage + db2 run, 2026-10-08 to 10-14 | none; *read off the terminal* |
| **druglike_10k (v1)** | 10,000 ZINC22 2D parents, H18–H22, one pooled uniform draw (seed 20260923); QupKake on all 10k; site filter v3 then v4, coupling off then on; 2026-10-15 | `validation/2026-10_druglike10k/druglike_10k/triage/stats.json`, `validation/2026-10_druglike10k/druglike_10k/parents/rejected.tsv`, `validation/2026-10_druglike10k/druglike_10k/library/stats.json`, `validation/2026-10_druglike10k/druglike_10k/compare/` (summary.json, sites.tsv, parents.tsv, window.png, scatter.png), `validation/2026-10_druglike10k/druglike_10k/db2_verification.csv`, `validation/2026-10_druglike10k/druglike_10k/ref_qk/reference_check.tsv` |
| **druglike_10k_v2 (v6, v6b)** | same 10,000 parents rebuilt with every fix through site filter v6, aromatic-N v2, the escalated 3D check; QupKake on the routed 35.8% only; `library_v6` before and `library_v6b` after the last-isomer guard; 2026-10-16 | `validation/2026-10_druglike10k/druglike_10k_v2/triage_v6/stats.json`, `validation/2026-10_druglike10k/druglike_10k_v2/library_v6b/stats.json`, `validation/2026-10_druglike10k/druglike_10k_v2/library_v6b/dropped.tsv`, `validation/2026-10_druglike10k/druglike_10k_v2/db2_verification_v6b.csv`, `validation/2026-10_druglike10k/druglike_10k_v2/timing_v6.txt`, `validation/2026-10_druglike10k/druglike_10k_v2/regressed.txt` |
| **aromatic-N reference** | 53 molecules in `pka_triage/reference/aromatic_n_bases.{smi,tsv}` through QupKake and MolGpKa (`llb pka-reference`), run inside druglike_10k | `validation/2026-10_druglike10k/druglike_10k/ref_qk/reference_check.tsv` |

Headline result of the pipeline as it stands (druglike_10k_v2, v6b): 9,797 of
10,000 parents built (98.0%); 18,620 of 19,102 contract names built; 0 charge
mismatches between library and db2; 135 parents lost because no stereoisomer
as ZINC specifies it can be built, 65 to conformator/chemistry-check failures,
3 to stereo-enumeration crashes or skipped names (`validation/2026-10_druglike10k/druglike_10k_v2/db2_verification_v6b.csv`,
`validation/2026-10_druglike10k/druglike_10k_v2/library_v6b/dropped.tsv`; per-parent split *read off the terminal*). Cost per 100k
parents at H18–H22: ~320 CPU-h QupKake (35.8% routed, 31.8 s/parent) +
~2,700 CPU-h db2 (1.91 names/parent, mean ~52 s/name), ~90 GB of db2
(`validation/2026-10_druglike10k/druglike_10k_v2/timing_v6.txt`).

## Versions in force

| constant | value | file |
| --- | --- | --- |
| `SITE_FILTER_VERSION` | 6 | `run_qupkake/qupkake_protomers.py` |
| `COUPLING_VERSION` | "1" | same |
| `AROMATIC_N_VERSION` | "2" | same |
| `--ph` | 6.4 7.4 8.4 | `qupkake_protomers.py`, `pka_triage/triage.py` |
| `--route-window` | 1.5 | `triage.py` |
| `--margin` | 1.0 | both |
| `--min-population` | 0.01 | both |
| `--max-states` | 8 | both |
| `--max-ambiguous` | 6 | both |
| `--coupling-bonds` (routing) | 3 | `triage.py` |
| `MAX_EXACT_SITES` | 12 | `qupkake_protomers.py` |
| `EMBED_TRIES` | 6 tries, below | `run_qupkake/assign_names.py` |
| `--max-stereo` | 32 | `assign_names.py`, `qupkake_protomers.py` (`DEFAULT_STEREO_CAP`) |
| `LLB_MMFF_VARIANT` / `_MAX_DISPL` / `_K` / `_MAX_ITERS` | eps4r_restrained / 0.3 Å / 100 kcal mol⁻¹ Å⁻² / 200 | `run_qupkake/build_ligand_relaxed.py` |
| RDKit | 2024.03.5 (db2_converter's pin; also `environments/llb_pka.yml`) | |
| db2_converter | hnlab @ 63d6656, unmodified | `config/hpc.env` |

Every shard records the versions it was built under, and `llb merge` refuses
to combine halves built under different ones (site filter, coupling,
aromatic-N, pH list, margin, min population, max states).

---

## 0. Sourcing pass (2026-10-16)

For each decision that the first draft labelled "chemistry, not a checked
number", or that rested on a secondary table, a primary source was looked for
and read (full text unless marked *abstract only*; publisher bot checks on some
ACS/RSC pages were not bypassed). Notes and downloaded copies:
`side-session-notes/.../llb-decisions-sourcing.md` (outside the repo).

| claim | result | source | effect on the record |
| --- | --- | --- | --- |
| acylsulfonamide / sulfonylurea N-H are acids near pH 7 | **sourced** | Lassalas et al., *J. Med. Chem.* 2016, 59, 3183–3203, doi:10.1021/acs.jmedchem.5b01963, Table 1 (capillary electrophoresis): acylsulfonamide PhCH2CH2CONHSO2Me 4.94; sulfonylurea BnNHCONHSO2Me 5.04; acyl dimethylsulfamide 5.86; 3-phenylpropionic acid 4.64; plain sulfonamides 10.04 and > 12 | 4.1b becomes sourced. The dimethylsulfamide (5.86) is only ~70–90% ionized at 6.4, which the state model handles. |
| sulfur is never a base | **sourced** (thioethers) | Reich H2O (Reich–Bordwell PDF p. 9, marked approximate): protonated Me2S −5.4, tetrahydrothiophene −4.4 (vs Me2O −3.8, THF −2.1). Bonvicini, Levi, Lucchini & Scorrano, *J. Chem. Soc., Perkin Trans. 2* 1972, 2267, doi:10.1039/P29720002267 (*abstract only*): alkyl methyl sulfides pK_BH+ −6.7 to −6.9 in aqueous H2SO4. Perdoncin & Scorrano, *J. Am. Chem. Soc.* 1977, 99, 6983 (title verified, values not read: paywall). | 4.6 sulfur part becomes sourced for sulfides. Thiophene, thiocarbonyl and sulfone S: no value found, judgment. |
| azide atoms are never bases | **not found** | no measured basicity of an organic azide found (HN3's 4.7 is the acidity of HN3 itself, not protonation of R-N3) | stays judgment |
| diazo atoms are never bases | **contradicted as worded** | Reich H2O p. 9: methyldiazonium CH3N2+ is a carbon acid of pKa 10.0, so diazomethane is thermodynamically C-protonated below pH ~10. The cation does not persist (water attacks it, N2 is lost; protonation by H3O+ ~10⁸ M⁻¹ s⁻¹, McGarrity & Smyth, *J. Chem. Soc., Chem. Commun.* 1977, 347, *abstract only*). No value found for diazocarbonyls. | 4.6 diazo part re-justified: the protonated form is a reactive intermediate, not a species to dock; the rule stays, the reason changes. |
| amino tautomer dominates (2-aminopyridine, -thiazole, -oxazole) | **partial** | Katritzky, Karelson & Harris, *Heterocycles* 1991, 32, 329, doi:10.3987/REV-90-421 (review), Scheme 24, p. 344: pK_T = 8.7 in water for 4-aminopyridine. Forlani & De Maria, *J. Chem. Soc., Perkin Trans. 2* 1982, 535, doi:10.1039/P29820000535 (*abstract only*): 2-aminothiazoles mainly amino in water, 25 °C; but 2-p-tosylaminothiazole mainly **imino**. 2-aminopyridine: no primary K_T found (a secondary overview says 1000:1, conflicting with the often-quoted ~10⁵). 2-aminooxazole: nothing. Angyal 1952 not readable (blocked). | 2.2: numeric only for 4-aminopyridine; 2-amino cases qualitative. N-sulfonyl aminothiazoles are a known exception, but the rule does not touch them (their =N carries the sulfonyl, not H). |
| Meldrum's acid water pKa ~4.8 | **sourced, value 4.97** | Pihlaja & Seilo, *Acta Chem. Scand.* 1969, 23, 3003–3010, doi:10.3891/acta.chem.scand.23-3003, Table 5, p. 3007: pKa° 4.97 (UV, zero ionic strength); the often-quoted 4.83 is probably Arnett & Harrelson, *J. Am. Chem. Soc.* 1987, 109, 809 (value not read). DMSO 7.3 (Reich DMSO p. 1). | 4.4 uses 4.97. |
| oxalic acid pKa2 (Williams 4.19 vs CRC 3.81) | **resolved toward Williams** | Pinching & Bates, *J. Res. Natl. Bur. Stand.* 1948, 40, 405, doi:10.6028/jres.040.031: thermodynamic pK2 4.266 at 25 °C, zero ionic strength. Williams itself lists 1.23/4.19 (p. 6) and 1.25/4.14 (p. 8). CRC 3.81 not read (paywall); a value that low fits a constant at finite ionic strength (inference). | 5.1: the 3-bond acid shift (2.4) stands; with 3.81 it would be 2.0. |
| succinic acid pKa2 | **Williams low by 0.16** | Pinching & Bates, *J. Res. Natl. Bur. Stand.* 1950, 45, 322, doi:10.6028/jres.045.035: pK2 5.638, pK1 ~4.21, 25 °C | 5.1: the 5-bond acid shift would be 0.8, not 0.7; within the table's precision, recorded, not changed. |
| protonated ethers −2 to −3.8 | **sourced** | Ripin & Evans (updated 4/16/2024), p. 1, "Protonated species": Me2O −3.8, THF −2.05, MeOH2+ −2.2 (negative values marked extrapolated) | as cited |
| Williams values used in the code | **all match; most page numbers do not** | checked against the PDF (rev. 4/7/2022): diacids are on p. 6 (code: 7); ethylenediamine, putrescine, cadaverine p. 14 (15); other diamines, piperazine, DABCO p. 19 (20); guanidines, amidines, biguanide, dimedone p. 28 (29); 2-aminopyridine p. 24 (27); 2-aminothiazole p. 27 (28); carbon acids pp. 30–31 (31–32). Williams is internally inconsistent in places (DABCO also 8.8 "?" on p. 14; nitromethane "strong" on p. 30, 10.29 on p. 31). The Reich H2O table differs slightly: diethyl malonate 12.9 (13.3), malononitrile 11.2 (11.81), benzamidine 11.2 (11.6). | This record cites the correct pages; the code comments, TRIAGE.md and the reference set's `source` column still carry the old ones (see "Inconsistencies found"). |

---

## 1. Sampling

### 1.1 Sample uniformly over whole heavy-atom bins, from files on disk

- **Decision.** A sample from a ZINC22 heavy-atom bin is a uniform draw over
  every molecule in every tranche file of that bin, read from a full local
  download: one reservoir pass (Algorithm L) across the files.
- **Governs.** `llb sample --local ROOT --urls LIST` (`mol_download/sample_2d.py`,
  `reservoir`, `sample_local`); `--urls` makes the sampler refuse a bin with any
  listed file missing, so a half-finished download cannot skew the logP mix.
  Seed 20260923.
- **Commits.** `5e728fd` (sampler, local mode, the bias measurement), `690f304`
  (refactor onto shared `bin_files`/`reservoir`/`write_bin`; per-bin output
  byte-identical to before).
- **Evidence.** Measured 2026-09-23 on H17M100, 3,600 molecules: the front of a
  tranche file is one registration batch. Prefix sample vs two uniform draws:
  unique Murcko scaffolds 1,479 vs 1,654 / 1,591; mean nearest-neighbour
  Tanimoto 0.469 vs 0.432 / 0.435; rings per molecule 1.378 vs 1.463 (KS
  p = 2e-5, against p = 0.86 between the two uniform draws); amides per
  molecule 1.502 vs 1.423 (KS p = 0.001 vs 1.0). String features (length,
  stereo, charge) are flat across file deciles and do not reveal the bias
  (`sample_2d.py` docstring). Algorithm L: Li 1994.
- **Alternatives.** *Streaming the front of each remote file* (the default mode,
  a few MB per bin): rejected for any real sample because it is clustered and
  less diverse, the wrong error for a chemical-space study; kept only for
  plumbing runs, and it prints a warning. *Random byte offsets into remote
  files*: not built; the server's per-connection rate (~14 MB/s measured) makes
  the full download (~5–19 h for H14–H28) affordable, and a local pass reads
  ~16M lines/s.
- **Status.** measured.

### 1.2 Proportional to bin size for the validation set; per-bin as an option

- **Decision.** The druglike_10k validation set is one uniform draw of 10,000
  over the pooled H18–H22 bins, so each bin is represented in proportion to its
  molecule count. Equal counts per bin (`--per-bin N`) remain available.
- **Governs.** `llb sample --total N` (`sample_local_pooled`; writes
  `POOLED.done` with each bin's population and draw).
- **Commits.** `690f304`.
- **Evidence.** Draw recorded in `POOLED.done` (*read off the terminal*): population 3,797,772,950; H18 529, H19 865, H20 1,576, H21
  2,455, H22 4,575. Tested before use on a synthetic tree of bins of 1,000 /
  3,000 / 6,000: a draw of 1,000 gave 91 / 296 / 613, mean over 20 seeds
  96 / 303 / 600.
- **Alternatives.** *Equal per bin* (2,000 × 5 was the first plan): rejected by
  John for the validation run, because the validation numbers (routed fraction,
  attrition, cost) should describe the library as it will be drawn. Per-bin is
  the better choice when a model should see every size range evenly; that is a
  per-use choice, not fixed here. Bin range H18–H22: chosen as drug-like
  (judgment); larger bins are untested by the triage (see 3.2).
- **Status.** judgment (the choice), measured (the draw).

---

## 2. Standardization and tautomers

### 2.1 Standardize: largest fragment, cleanup, neutralize

- **Decision.** Every input SMILES is reduced to its largest fragment, cleaned
  and neutralized (RDKit `rdMolStandardize`: `FragmentParent`, `Cleanup`,
  `Uncharger`) before anything else, and deduplicated on that structure. Both
  pKa predictors are given this neutral parent.
- **Governs.** `qupkake_protomers.standardize` (imported by
  `run_qupkake/prepare_parents.py`, so the dedupe key is exactly the structure
  QupKake sees).
- **Commits.** `5e728fd`.
- **Evidence.** QupKake expects neutral, single-fragment input (each site's conjugate is
  made from the input state); MolGpKa uncharges its input
  itself (MolGpKa `predict()`), so both tools see the same molecule.
- **Alternatives.** Keeping ZINC's supplied charge state: rejected, because the
  pipeline's job is to assign the states itself at three pH values.
- **Status.** judgment.

### 2.2 Imino → amino tautomer for unsubstituted exocyclic =NH

- **Decision.** A parent with an unsubstituted exocyclic `=NH`, ring
  (`N=c1[nH]…`: the imino form of 2-/4-aminopyridine, aminothiazole,
  aminooxazole, aminotriazole) or acyclic, is moved to the tautomer without it.
  Only protons on N move; stereo elsewhere is kept. Molecules with no such
  tautomer (N-alkyl iminopyridines, plain amidines `C(=NH)NH2`) are left as
  given.
- **Governs.** `qupkake_protomers.amino_tautomer`, SMARTS `[NX2H1;!R]=[#6]`;
  `parents.tsv` column `tautomer_fix`; `--keep-imino` turns it off.
- **Commits.** `7676adf` (ring imines), `ba0b54d` (acyclic, see 2.3).
- **Evidence.** *Chemistry:* the amino form dominates in water for these
  heterocycles: qualitative for the 2-amino cases (Angyal 1952; for
  2-aminothiazole, Forlani & De Maria 1982), numeric only for 4-aminopyridine,
  pK_T 8.7 (Katritzky, Karelson & Harris 1991, Scheme 24; section 0), and the
  measured base pKa is the amino form's ring N: 2-aminopyridine 6.86,
  2-aminothiazole 5.39 (Williams pp. 24, 27; also the aromatic-N reference set,
  where both tools reproduce these to within 0.5 alone). *Build:* the `=NH`
  bond crashes RDKit 2024.03.5's `EnumerateStereoisomers` ("Range Error idx2":
  no atom to define E/Z), and db2_converter then skips the molecule. On
  druglike_10k v1, 247 names crashed check-stereo and had no db2 (*read off the
  terminal*: `library.stereo_crashed.tsv` was not copied; the 248
  `MISSING` in `validation/2026-10_druglike10k/druglike_10k/db2_verification.csv` are these plus one
  hand-excluded name). In v2, 65 of 10,000 parents were converted, none
  duplicating another parent; stereo-crash names fell from 247 to 5 (with
  site filter v5/v6; *read off the terminal*).
- **Alternatives.** *QupKake's own tautomer search* (`-t`, xtb): rejected as much
  slower, and it rewrites the skeleton, breaking the atom-order correspondence
  the stereo carry-over relies on. *A general canonical-tautomer step* (RDKit's
  `TautomerEnumerator.Canonicalize`): rejected because it moves protons on
  molecules that need no change and can drop stereo; the rule touches only the
  `=NH` case that is both chemically wrong and a build crash. *Keep and skip*
  (accept the crash): rejected, it silently loses ~1% of names.
- **Status.** sourced (pKa values; dominance qualitative unless section 0 says
  otherwise) + measured (crash counts).

### 2.3 Imidic acid → amide; acyclic =NH

- **Decision.** The same step also covers acyclic `=NH` and lets the proton come
  from O as well as N: a primary carbamate drawn as its imidic acid
  `RO-C(=NH)-OH` becomes `RO-C(=O)-NH2`; an amidoxime ether with the `=NH` on
  the wrong N is moved. Plain amidines/guanidines, amides, carbamates, imidate
  esters and N-alkyl iminopyridines are untouched.
- **Governs.** `amino_tautomer` (changed atoms must be N or O and must include
  the exocyclic `=NH`; ties broken by fewest changed atoms, then RDKit's
  tautomer score).
- **Commits.** `ba0b54d`.
- **Evidence.** druglike_10k_v2 (v6): the 5 remaining stereo-crash names came
  from 2 parents with exactly these two patterns (*read off the terminal*;
  `library_v6/library.stereo_crashed.tsv`). Amide over imidic acid in water:
  judgment (standard; no number checked here).
- **Alternatives.** As 2.2.
- **Status.** measured (the trigger) + judgment (amide dominance).
- **Note.** Takes effect only from a fresh `llb parents`; the v6/v6b library was
  built before it, which is where its last few crash names come from.

---

## 3. pKa prediction and triage

### 3.1 Two predictors: MolGpKa for every parent, QupKake for the ambiguous ones

- **Decision.** MolGpKa predicts per-site pKa for every parent; parents it
  cannot settle are routed to QupKake. Protonation states are built by one
  shared state model from whichever tool's sites the parent ends up with
  (`site_source` = molgpka, qupkake, or molgpka_fallback when QupKake failed on
  a routed parent).
- **Governs.** `llb triage` (`pka_triage/triage.py`), `llb merge`
  (`merge_protomers.py`, `--on-qupkake-fail molgpka` default),
  `qupkake_protomers.build_protomers` (shared). MolGpKa weights sha256-checked
  from `MOLGPKA_DIR`.
- **Commits.** `a0b876e` (MolGpKa port), `543e962` (triage, merge, shared state
  model), `a83d287` (`llb compare`).
- **Evidence.** *Cost:* QupKake ~3.6 min/molecule on one core on Owen's 10
  largest, almost all xtb, i.e. ~60,000 CPU-h per million (*read off the terminal*);
  28.0 s/parent on druglike_10k (77.83 CPU-h for 10,000; `llb timing`,
  *read off the terminal*); 31.8 s/parent on the routed v2 parents
  (`validation/2026-10_druglike10k/druglike_10k_v2/timing_v6.txt`). MolGpKa port: 4.5 ms/molecule on one laptop thread;
  triage end to end 9.7–13 ms/parent on Great Lakes. *Port fidelity:* max
  difference from upstream 3.8e-6 pKa (16 molecules / 47 sites, local) and
  6.7e-6 (Owen's 50, 242 sites, Great Lakes). The port keeps upstream's
  always-zero H-bond donor/acceptor feature columns, because the weights were
  trained that way; "fixing" them moves predictions by up to 3 units.
  *Agreement on druglike_10k* (QupKake on all 10k, 9,772 compared; filter v3,
  coupling off; `validation/2026-10_druglike10k/druglike_10k/compare/summary.json`): per-site MAE 1.25 over
  10,865 shared sites; aliphatic amine 0.75 (n = 3,037), carboxylic acid 0.49
  (293), pyridine-like N 1.15 (1,920, MolGpKa lower by 0.86 on average),
  imidazole-like N 1.16 (2,823, lower by 0.88); dominant microstate agrees
  83% / 85% / 86% at pH 6.4 / 7.4 / 8.4.
- **Alternatives.** *QupKake on everything:* ~780 CPU-h per 100k parents at
  28 s; affordable at 10k, not at the millions the project needs, and the
  reference set (6.1) shows it is not uniformly better. *MolGpKa alone:* no
  second opinion where a site sits in the window; and MolGpKa was trained on
  computed, not measured, ChEMBL pKa values (MolGpKa paper), so its errors
  near the window matter most. *Dimorphite-DL for the states:* not needed;
  `apply_changes` builds states from the predictor's own atom indices, and
  Dimorphite's SMARTS would have to be mapped back onto them; kept in the env
  as an independent check only. *Caching by ionizable fragment:* rejected,
  MolGpKa is whole-graph (5 GCN layers + attention), so a site's pKa depends on
  the whole molecule.
- **Status.** measured.
- **Caveat a reviewer will raise.** MolGpKa–QupKake agreement is model against
  model. Where one tool was right was settled only on the 53-molecule
  aromatic-N reference set (6.1) and on chemistry with cited values (section
  4). There is no broad experimental pKa benchmark in this record.

### 3.2 Route window 1.5 pKa units

- **Decision.** A parent goes to QupKake if any MolGpKa site (after the site
  filter) lies within 1.5 units of any target pH.
- **Governs.** `--route-window 1.5` (`triage.py`); `triage/stats.json` records the
  routed fraction at windows 0–4.
- **Commits.** `543e962`.
- **Evidence.** The 1.5 was John's proposal before any data (|pKa − pH| > 1.5
  puts a site > 97% in one form); it was then measured on druglike_10k, QupKake
  on all 10k, filter v4, coupling on (`compare_c`, *read off the terminal*; the
  filter-v3 `validation/2026-10_druglike10k/druglike_10k/compare/` gives the same curve within
  1 point):

  | window | routed | settled | settled but disagree with QupKake |
  | --- | --- | --- | --- |
  | 0.5 | 24% | 7,432 | 679 (9.1%) |
  | 1.0 | 30% | 6,823 | 378 (5.5%) |
  | **1.5** | **36%** | **6,255** | **207 (3.3%)** |
  | 2.0 | 41% | 5,748 | 129 (2.2%) |
  | 3.0 | 47% | 5,154 | 72 (1.4%) |

  Reading the 165 disagreements of the first compare by charge pattern
  (*read off the terminal*): ~95 were MolGpKa +1 against QupKake +2 to +5 at every
  pH, i.e. QupKake overcharging (later traced to guanidine splitting, 4.5, and
  the aromatic-N artefact, 6.1); ~25 MolGpKa −1 vs 0; ~30 at pH 8.4 only. So
  most "misses" at 1.5 were QupKake errors, and widening the window would
  route more of them in. Cost: moving 1.5 → 2.0 adds ~40 CPU-h per 100k
  parents, small next to db2, so the window is an accuracy choice, not a cost
  one. v2 routed 35.8%.
- **Alternatives.** *2.0 or wider:* rejected, since most added disagreements are
  QupKake errors. *1.0:* defensible (5.5% disagree), not chosen.
- **Status.** measured (H18–H22 only; untested on larger molecules; the
  pipeline plan runs a `--validate` comparison on the H18–H30 set).

### 3.3 Target pH 6.4, 7.4, 8.4

- **Decision.** States are enumerated at pH 6.4, 7.4 and 8.4; a protomer is kept
  if it passes the population filter at any of them.
- **Governs.** `--ph 6.4 7.4 8.4` (both paths; recorded in every shard; merge
  refuses a mismatch).
- **Commits.** Default present since `5e728fd` (2026-09-23, when `run_qupkake`
  joined the repo); predates this repo's history.
- **Evidence.** None checked. 7.4 is physiological pH; ±1 brackets the range of
  binding-site and assay pH the library should cover. These three values were
  already in John's QupKake setup when this work started.
- **Alternatives.** *7.4 only:* rejected as too narrow for targets whose pocket
  pH or local environment shifts protonation. *A continuous range* (e.g. ZINC's
  "mid" 6–8 style): equivalent in spirit; three points are what the state model
  evaluates.
- **Status.** judgment. A reviewer can ask for a citation here; we have none in
  hand.

### 3.4 Margin 1.0, minimum population 0.01, at most 8 states, at most 6 ambiguous sites

- **Decision.** A site within 1.0 unit of the pH (after coupling, its
  population-weighted log ratio) is enumerated both ways; outside it is locked
  to its dominant form. A state is dropped below 1% estimated population; at
  most 8 states per parent per pH are kept (highest population first); at most
  6 sites are enumerated simultaneously (2⁶ combinations), the rest locked.
- **Governs.** `--margin 1.0`, `--min-population 0.01`, `--max-states 8`,
  `--max-ambiguous 6` (same meaning in both paths).
- **Commits.** Defaults since `5e728fd`; margin kept separate from the route
  window by design in `543e962`.
- **Evidence.** At |pKa − pH| = 1.0 the minor form is ~9%. No measurement chose
  these values.
- **Alternatives.** Using one threshold for routing and enumeration: rejected;
  the route window decides which *molecules* get the expensive predictor, the
  margin which *sites* get both forms, and they answer different questions.
- **Status.** judgment.

### 3.5 Route same-kind charged pairs within 3 bonds

- **Decision.** A parent is also routed if two same-kind MolGpKa sites (two
  bases or two acids) are both on the charged side at some pH within 3 bonds.
- **Governs.** `--coupling-bonds 3`, `--coupled route` (`triage.py`).
- **Commits.** `543e962`.
- **Evidence.** druglike_10k: 260 parents routed for coupling (251 also near the
  window, 9 coupling alone; `validation/2026-10_druglike10k/druglike_10k/triage/stats.json`). Since
  `85cdb85` the state model itself corrects coupling (section 5), so this route
  is a review flag rather than the fix.
- **Alternatives.** Not routing them: possible now that coupling is modelled;
  kept because it costs 9 parents in 10,000.
- **Status.** judgment.

### 3.6 Molecules QupKake cannot embed, or fails on, get MolGpKa's states

- **Decision.** Before calling QupKake, each routed molecule is put through the
  same RDKit calls QupKake's input step makes (`AddHs`, `EmbedMolecule`,
  `MMFFOptimizeMolecule`) under three seeds; one that fails any is left out
  (`qupkake_embed_failed`) and the merge gives it MolGpKa's states.
- **Governs.** `qupkake_protomers.qupkake_can_embed(tries=3)`;
  `--on-qupkake-fail molgpka` (default; `drop` is the alternative).
- **Commits.** `e547c3b`.
- **Evidence.** druglike_10k: 189 of 1,000 ten-molecule QupKake shards died with
  "Bad Conformer Id" because one molecule in each did not embed (QupKake's CLI
  embeds every row without checking, then exits). After the fix, 206 molecules
  were left out and 22 failed inside QupKake; 73 routed parents in the v1 and
  v2 merges are `molgpka_fallback` (*read off the terminal*).
- **Alternatives.** Dropping those parents (`--on-qupkake-fail drop`): rejected
  as default, since MolGpKa's states are a reasonable answer and the parent
  stays in the library, marked.
- **Status.** measured.

---

## 4. Site filter (implausible sites, v1–v6)

Sites no predictor should be trusted on are dropped before any state is built,
in both paths. Each rule is a SMARTS in `qupkake_protomers.IMPLAUSIBLE_SITES`;
dropped sites are listed per shard (`dropped_sites.tsv`, triage `sites.tsv`
`kept=0`, `dropped_reason`). `SITE_FILTER_VERSION` is written into every output
and merge refuses mixed versions; old QupKake runs are rebuilt with `llb
reenumerate` (no xtb). `--keep-amide-sites` turns the **whole** filter off,
not only the amide rules.

MolGpKa sites dropped on druglike_10k_v2 (`validation/2026-10_druglike10k/druglike_10k_v2/triage_v6/stats.json`): amide_nh_acid
7,568; alcohol_oh_acid 2,160; amide_n_base 1,440; amine_nh_acid 451;
carbon_acid 54; sulfur_base 17; oxygen_base 9; amidine_amino_n_base 5. Most of
these were far from the window and would never have mattered; the rules exist
for the minority that were not.

### 4.1 v1: amide-type N-H is not an acid; amide-type N is not a base

- **Decision.** Drop as acids: N-H of amides, anilides, ureas, carbamates,
  hydrazides, thioamides (`amide_nh_acid`). Drop as bases: any N on C=O, C=S
  or S(=O)=O (`amide_n_base`). Keep as acids: imide N-H (between two C=O) and
  N-sulfonyl N-H (acylsulfonamides, sulfonylureas; see 4.1b).
- **Governs.** `amide_nh_acid`, `amide_n_base`; `SITE_FILTER_VERSION` 1.
- **Commits.** `543e962` (amide N bases, MolGpKa side), `5cb1450` (both rules,
  both paths), `1a1fb8b` (sources).
- **Evidence.** *Measured:* Owen's 50: QupKake put 44 amide-type N-H at pKa
  2.4–6.7, making an amide anion the dominant state at pH 7.4 in 22 of 50
  parents; all 5 parents the triage settled "against" QupKake were this
  (*read off the terminal*). druglike_10k: 255 parents had a QupKake amide/imide anion
  at mid pH before the filter (`validation/2026-10_druglike10k/druglike_10k/compare/summary.json`). MolGpKa put
  anilides at 3.9–6.0 as bases and some anilides at 6–7 as acids. *Sourced:*
  acetamide N-H 15.1 in water (Ripin & Evans, "Amides & carbamates"; marked
  extrapolated); DMSO: acetamide 25.5, N-methylacetamide 25.9, benzamide 23.3,
  acetanilide 21.5, urea 26.9, ethyl carbamate 24.2, thioacetamide 18.5, against
  succinimide 14.7 (Reich DMSO, "Amides (NH)", "Carbamates (NH)", "Imides
  (NH)"); phthalimide 8.30 in water (Ripin & Evans, "Imides"). Amide
  protonation is on O and weak: protonated acetamide −0.62 ± 0.07,
  N-methylacetamide −0.42, N,N-dimethylacetamide −0.28 (Grant 1983, water,
  298 K).
- **What the sources do not give.** A water pKa for anilides or ureas
  specifically; that part rests on the DMSO ordering (every amide-type N-H 4–12
  units above an imide, whose water value is ~8).
- **Alternatives.** Trusting the predictors: rejected on the numbers above.
  Dropping imides too: rejected, phthalimide is 8.30.
- **Status.** sourced + measured.

### 4.1b Acylsulfonamide and sulfonylurea N-H kept as acids

- **Decision.** N-H between C=O and S(=O)=O stays an acid site.
- **Governs.** the `!$([#7]-[#16](=[#8])=[#8])` exclusion in `amide_nh_acid`.
- **Commits.** `5cb1450`.
- **Evidence.** Lassalas et al. 2016, Table 1 (measured, capillary
  electrophoresis): acylsulfonamide PhCH2CH2CONHSO2Me 4.94, sulfonylurea
  BnNHCONHSO2Me 5.04, acyl dimethylsulfamide 5.86, against 4.64 for the parent
  carboxylic acid. Ionized at pH 7.4.
- **Alternatives.** Dropping them with the amides: rejected; they are carboxylic
  acid isosteres, ionized at pH 7.4.
- **Status.** sourced (sourcing pass, section 0).

### 4.2 v2: amine N-H is not an acid

- **Decision.** Drop N-H acid sites on alkyl, cyclic, aryl, diaryl amines and
  enamines (`[N-]`/`[NH-]` anions). N bonded to anything but C/H, or to C=N/O/S,
  is left to other rules.
- **Governs.** `amine_nh_acid`; `SITE_FILTER_VERSION` 2 (first version written
  into shards and checked by merge).
- **Commits.** `e9424c1`.
- **Evidence.** Owen's first db2 build had `[NH-]`/`[N-]` amine-anion protomers
  (*read off the terminal*). NH3 38 in water (extrapolated) and 41 in DMSO (Ripin &
  Evans, "Amines"); DMSO: pyrrolidine 44, aniline 30.6, diphenylamine 25.0
  (Reich DMSO, "Amines (NH)").
- **Correction recorded.** The commit first blamed three `1generate` build
  failures on these anions; Owen's rebuild under v2 showed the parent fails the
  same way with sane protomers (`ac68df8`), so the failure was the scaffold's.
- **Status.** sourced.

### 4.3 v3: oxygen is never a base; alcohol O-H is never an acid

- **Decision.** Drop every O as a base (`C(=O)[OH2+]`, protonated nitro,
  oxonium) and O-H on an sp3 carbon as an acid (alkoxide). Phenols, enols,
  carboxylic acids, oximes and hydroxamic acids keep their acid sites.
- **Governs.** `oxygen_base` (`[#8]`), `alcohol_oh_acid` (`[#8X2;H1][CX4]`);
  `SITE_FILTER_VERSION` 3.
- **Commits.** `121e4ad`.
- **Evidence.** Owen's rebuild (2026-10-14) had `C(=O)[OH2+]`, a protonated nitro
  and an alkoxide among QupKake's protomers (*read off the terminal*). Water (Reich
  H2O, "Oxygen Acids"): protonated nitrobenzene −12.4, protonated benzoic acid
  −7.8, Me2OH+ −3.8, MeOH2+ −2.2; methanol 16, tert-butanol ~20, phenol 10.0.
- **Status.** sourced.

### 4.4 v4: carbon is not an acid, except 1,3-dicarbonyl CH and nitro CH

- **Decision.** Drop C-H acid sites unless the carbon sits between two
  ketone/aldehyde carbonyls or carries a nitro group.
- **Governs.** `carbon_acid`; `SITE_FILTER_VERSION` 4.
- **Commits.** `25c3b01`.
- **Evidence.** MolGpKa put malonic acid's CH2 at 7.46 and, with coupling off,
  the state model built a `[CH-]` trianion at 90% (local test, *read off the
  terminal*). Water (Williams "Carbon Acids", pp. 30–31; dimedone p. 28): kept
  classes reach the window — acetylacetone 8.95, benzoylacetone 8.23 (enol)/9.4,
  dimedone 5.23, triacetylmethane 5.85, malondialdehyde 5, 2-nitropropane 7.74,
  nitroethane 8.6, nitromethane 10.29; dropped classes do not — ethyl/methyl
  acetoacetate 10.68/10, ethyl cyanoacetate 9, malononitrile 11.81, diethyl
  malonate 13.3, acetone ~20, acetonitrile ~25, dimethyl sulfone ~23. 54
  MolGpKa sites dropped on 10k parents.
- **Known miss.** Meldrum's acid (cyclic malonate diester; water pKa°
  4.97, Pihlaja & Seilo 1969, Table 5) is dropped by the rule, as is any other
  doubly ester-activated CH that is unusually acidic. Enol forms are O-H sites
  and unaffected. Accepted as rare in drug-like ZINC; not counted.
- **Alternatives.** Keeping ester-activated CH: rejected; they sit at 9 and above
  except the Meldrum's-type cyclic case.
- **Status.** sourced; the Meldrum's exception is a recorded, accepted miss.

### 4.5 v5: one basic site per amidine or guanidine, on the =N

- **Decision.** Drop basic sites on the amino N's of an amidine or guanidine
  (non-aromatic `N-C=N`). If the predictor gave the =N no site, the
  highest-pKa amino-N site moves onto the =N (`moved_from_amidine_amino_n`), so
  the protonation is never lost.
- **Governs.** `amidine_amino_n_base` (`[#7;!a;$([#7]-[#6;!a]=[#7;!a])]`) and
  `filter_sites`; `SITE_FILTER_VERSION` 5. Aromatic C is not matched
  (2-aminopyridine: the ring N is the base).
- **Commits.** `19dedb6`.
- **Evidence.** druglike_10k: QupKake gave each amino N its own site at
  10.7–11.1 beside the =N at 11.5–11.8, building +2/+3 guanidines, and those
  `NC(=N)[NH3+]`-type protomers crashed RDKit's stereo enumeration (*read off the
  terminal*). One pKa each: guanidine 13.71, acetamidine 12.52, benzamidine 11.6,
  phenylguanidine 10.88; biguanide 13.25 then 3.07, showing how far a second
  proton on one conjugated system sits (Williams p. 28).
- **Alternatives.** Leaving it to the coupling penalty: rejected; a 2-bond
  penalty of 2.2 does not bring sites 3–4 units above pH back below it.
- **Status.** sourced + measured.

### 4.6 v6: sulfur is never a base; azide and diazo atoms are never bases

- **Decision.** Drop every S as a base (`sulfur_base`) and every atom of an azide
  or diazo group as a base (`azide_base`).
- **Governs.** `sulfur_base` (`[#16]`), `azide_base`; `SITE_FILTER_VERSION` 6.
- **Commits.** `2fac899`.
- **Evidence.** druglike_10k_v2 (v2 merge): a sulfonium `[SH+]` among QupKake's
  +2 protomers, and azides protonated on the terminal N (`N=[N+]=N`) in every
  protomer of 38 names, which also crash RDKit's stereo enumeration
  (*read off the terminal*). 17 MolGpKa sulfur sites dropped on v6. *Sulfur:* protonated
  sulfides are stronger acids than protonated ethers: Me2SH+ −5.4,
  tetrahydrothiophene-H+ −4.4 (Reich H2O p. 9, approximate); alkyl methyl
  sulfides −6.7 to −6.9 (Bonvicini et al. 1972, abstract); protonated ethers
  −2.05 to −3.8 (Ripin & Evans p. 1). No value found for thiophene, C=S or
  sulfone S. *Azides:* no measured basicity found; not protonated at pH 6–8 is
  judgment. *Diazo:* methyldiazonium is a carbon acid of pKa 10.0 (Reich H2O
  p. 9), so diazomethane is thermodynamically protonated at pH 7; the rule is
  kept because CH3N2+ is a reactive intermediate (it loses N2 in water), not a
  species to dock, and no value exists for the diazocarbonyls that make up
  drug-like diazo compounds.
- **Status.** sourced (sulfides) + judgment (thiophene, C=S, sulfone S,
  azides, diazo; the diazo rule rests on reactivity, not basicity).

---

## 5. Coupling model

### 5.1 Penalize same-kind charges by bonds apart, from diamine/diacid micro-pKa shifts

- **Decision.** Every combination of sites is scored together. A state's log10
  weight is the sum over charged sites of (pKa − pH) for bases or (pH − pKa)
  for acids, minus a penalty SHIFT for each pair of same-kind charged sites:

  | bonds apart | bases | acids | from (Williams, water, 25 °C) |
  | --- | --- | --- | --- |
  | ≤ 3, same ring | 3.5 | 3.5 | piperazine 9.82/5.68 → 3.5; DABCO 8.19/4.18 → 3.4 (p. 19) |
  | ≤ 3 | 2.2 | 2.4 | ethylenediamine 9.98/7.52 → 1.9, 1,2-diaminopropane 10.00/7.13 → 2.3, 2,3-diaminobutane 10.00/6.91 → 2.5 (pp. 14, 19); oxalic 1.23/4.19 → 2.4 (p. 6) |
  | 4 | 1.1 | 2.3 | 1,3-diamino-2,2-dimethylpropane 10.22/8.18 → 1.4, 1,3-diaminopropan-2-ol 9.68/8.23 → 0.85 (p. 19); malonic 2.83/5.69 → 2.3 (p. 6) |
  | 5 | 0.9 | 0.7 | putrescine 10.84/9.31 (p. 14); succinic 4.19/5.48 (p. 6) |
  | 6 | 0.7 | 0.5 | cadaverine 11.05/9.74 (p. 14); glutaric 4.34/5.42 (p. 6) |
  | 7 | 0.4 | 0.4 | adipic 4.42/5.41 (p. 6); bases borrow the acid value |
  | 8–10 | 0.3 | 0.3 | extrapolated, not data |
  | > 10 | 0 | 0 | |

  Each SHIFT is (pKa2 − pKa1) corrected by log10(4), the statistical part of
  the gap for two equivalent sites. Sites more than `margin` on the uncharged
  side are held uncharged (the penalty only pushes them further); up to 12
  variable sites are scored exactly (2¹²), beyond that the most strongly
  charged are held charged with their penalties still counted. A site is
  enumerated both ways when its population-weighted log ratio is within the
  margin. States are merged by SMILES, populations summed (piperazine's two
  monocations are one species).
- **Governs.** `qupkake_protomers.Coupling`, `BASE_SHIFT`, `ACID_SHIFT`,
  `RING_SHIFT`, `FAR_SHIFT`, `FAR_BONDS`, `MAX_EXACT_SITES = 12`;
  `COUPLING_VERSION = "1"`; `--coupling on` (default in QupKake path, triage,
  compare, reenumerate).
- **Commits.** `85cdb85`.
- **Evidence.** *Why:* both predictors predict each site with the others
  neutral, so the independent Henderson–Hasselbalch product charges every amine
  whose pKa clears the pH. MolGpKa's own piperazine (10.46 on both N) locked
  the dication; experimentally the second pKa is 5.68 and the monocation
  dominates at 7.4. *Checks* (local, *read off the terminal*): with all penalties
  zero the model reproduces the independent states exactly (3,000 random site
  sets, 0 differences); ethylenediamine at pH 7.4: 37% dication in the model vs
  43% from its measured pKas; piperazine at 6.4: 17% vs 16%; a
  piperazine–propyl–piperazine went from +4 to +2. *On druglike_10k:*
  settled-but-disagree stayed at 207 with coupling on (it was not the main
  cause; see 4.5 and 6.1); dominant-state agreement moved from 83/85/86% to
  84/86/87%. *Source check* (section 0):
  oxalic pK2 is 4.266 thermodynamic (Pinching & Bates 1948), supporting
  Williams' 4.19 over the CRC's 3.81; succinic pK2 is 5.638 (Pinching & Bates
  1950), so the 5-bond acid shift would be 0.8 rather than 0.7, left as is.
- **Limits, stated in the code.** Opposite charges are not coupled (the
  correction only removes charge; zwitterions are scored independently).
  Through conjugated paths (pyrazine, pyrimidine, amidines) real coupling is
  stronger than these aliphatic references, so it is under-corrected there. The
  predicted pKas themselves are unchanged: MolGpKa's 10.46 piperazine is still
  65% dication at 6.4.
- **Alternatives.** *Cap total charge by counting basic centres:* rejected as
  cruder than a measured shift. *A fixed 1–2 unit shift within 3 bonds:* the
  first proposal; replaced by the table so every number has a source. *Rely on
  QupKake's micro-pKas:* not possible; QupKake as run also predicts each site
  from the neutral input state.
- **Status.** sourced (shifts ≤ 7 bonds) + judgment (8–10-bond extrapolation,
  ≤ 2-bond use of the 3-bond value, bases at 7 bonds) + measured (checks).

---

## 6. Aromatic-N hybrid

### 6.1 In routed parents with a base off the ring, aromatic-N basic sites come from MolGpKa

- **Decision.** In a routed parent that has any basic site (either tool's, any
  pKa) on an atom that is neither aromatic nor bonded to an aromatic atom, the
  basic sites on aromatic N are MolGpKa's; QupKake keeps everything else.
  Applied after the QupKake array, without xtb (`llb reenumerate
  --molgpka-triage`); every swap goes to `hybrid_sites.tsv` and the protomers
  carry `aromatic_n_from_molgpka`.
- **Governs.** `qupkake_protomers.hybrid_aromatic_n`, `AROMATIC_N_VERSION = "2"`;
  merge refuses QupKake shards without it unless `--qupkake-aromatic-n`.
- **Commits.** `3c5202c` (reference set, `llb pka-reference`), `11f4a56` (v1),
  `586401d` (v2).
- **Evidence.** druglike_10k: QupKake put pyrazine, oxazole, 1,2,3-triazole and
  tetrazole ring N at 9–10 where MolGpKa gave 2.7–4.2 (pyrazine measured 0.6);
  186 of 4,743 shared aromatic-N sites had QupKake ≥ 6.4 with MolGpKa ≤ 4
  (*read off the terminal*). The reference set separated the cases
  (`validation/2026-10_druglike10k/druglike_10k/ref_qk/reference_check.tsv`; measured values from Williams,
  cited per row in `pka_triage/reference/aromatic_n_bases.tsv`):

  | | MolGpKa | QupKake |
  | --- | --- | --- |
  | measured rings, alone | MAE 0.51 (19 rings), worst −2.19 | MAE 0.66 (20 rings), worst +1.97 |
  | ring two carbons from an N-ethylguanidine (~12) | 3.1–7.5 | 9.8–10.9, every ring |
  | ring two carbons from a dimethylamine (~8.5) | 3.1–6.7 | 7.7–9.8 |

  QupKake alone → tethered (guanidine / dimethylamine): oxazole 2.6 → 10.3 /
  8.6; pyrazine 2.9 → 10.5 / 8.9; pyrimidine 2.2 → 10.6 / 9.3; tetrazole 3.0 →
  9.8 / 7.7. The ring N takes the strong base's value whatever the ring, so it
  is an artefact, not chemistry. Mechanism: unverified guess that xtb moves the
  proton to the strong base while optimizing the ring-protonated structure.
  Alone, QupKake is the better tool on 4-aminopyridine (9.49 vs measured 9.17;
  MolGpKa 6.98), which is why the rule is conditional. v2 run: the swap applied
  to 735 parents (21% of routed); QupKake protomers at +3/+4 fell from 269 to
  76 together with v5, and from 3.05 to 2.53 protomers per routed parent
  (*read off the terminal*).
- **v1 → v2.** v1 (`11f4a56`) triggered only on an off-ring base at pKa ≥ 7. On
  the v2 run a routed tetrazole kept QupKake's ring-N values (5.61–5.81)
  beside a dimethylamine both tools put at 5.63 (lowered by the adjacent
  tetrazole), and built a tetrazolium (MolGpKa: 2.69). Same artefact, weaker
  base; the threshold was removed (`586401d`).
- **Exclusion.** Aryl amines (aniline 4.6, the NH2 of an aminopyridine) do not
  trigger it: they are weak bases whatever a tool predicts, and not the kind of
  group that pulled QupKake's ring values up.
- **Alternatives.** *Always MolGpKa for aromatic N:* slightly better alone overall
  but loses QupKake's 4-aminopyridine-type values; not chosen. *Stop routing
  for aromatic-N basicity:* rejected as it would also drop the cases QupKake
  gets right.
- **Open.** MolGpKa raises a ring N by 2–4 units next to a neutral amine
  (pyrazine 2.1 → 6.1); no measurement exists for those molecules. The coupling
  penalty (the amine is charged first, 0.9–1.1 at 4–5 bonds) is what keeps it
  in check.
- **Status.** measured (reference set) + sourced (the measured ring pKas).

---

## 7. Stereo closure, 3D embedding check, escalation, last-isomer guard

### 7.1 Close every protomer under db2_converter's stereo enumeration

- **Decision.** Every protomer is expanded with exactly the options db2_converter
  uses (`StereoEnumerationOptions(tryEmbedding=False)`), each isomer gets its
  own stereo index, and the result must be a fixed point: re-enumerating any
  output SMILES gives one isomer, so `--checkstereo` renames nothing.
  `check-stereo` re-runs that test with db2_converter's own function in its
  own environment before any build.
- **Governs.** `assign_names.closure`, `DB2C_OPTS`; `--max-stereo 32`
  (db2_converter skips > 2⁵ as `0enumerate_many`); `llb check-stereo`.
- **Commits.** `5e728fd`; `d140a24`, `dd4a03d` (check-stereo reports RDKit
  crashes instead of stopping).
- **Evidence.** Protonation creates stereocentres the neutral parent lacks:
  `CCN(C)Cc1ccccc1` 1 isomer, `CC[NH+](C)Cc1ccccc1` 2, `C[C@H]1CC[NH+](C)CC1` 2.
  db2_converter would rename these `NAME.0/NAME.1`, and mol2db2 keeps the last
  16 characters, breaking the contract inside the db2. druglike_10k v1: closure
  added 5,499 isomers (1,828 of them N-H+ invertomers;
  `validation/2026-10_druglike10k/druglike_10k/library/stats.json`); v6b: 19,102 names,
  0 renamed in the db2 (`validation/2026-10_druglike10k/druglike_10k_v2/db2_verification_v6b.csv`). One name per run still
  differs in perception between the naming and db2_converter environments
  (`Cm50000037lNj10C`) and is excluded by hand (see open items).
- **Alternatives.** Stereo enumerated on the neutral parent only: rejected, it
  misses protonation-created centres (above).
- **Status.** measured.

### 7.2 Drop stereoisomers that cannot exist in 3D

- **Decision.** Isomers with a specified tetrahedral centre on an atom in two or
  more rings (fused, bridged, spiro, cage) are embedded first; those that
  cannot be embedded are dropped (`stereo_cannot_embed` in `dropped.tsv`).
- **Governs.** `assign_names.needs_embed_check`, `embed_verdict`;
  `--no-embed-check`.
- **Commits.** `aa060ef`.
- **Evidence.** The closure enumerates each centre independently, as
  db2_converter does, so on cages it names configurations with no 3D
  structure. druglike_10k v1: 11 of 14 inspected `3chemistrycheck` failures
  fail to embed (conformator builds the possible isomer, which then fails the
  stereo check), and 285 of 314 `1generate` names cannot embed (296 have a
  centre on a fused/bridged atom); 64 parents lost every name this way though
  each had buildable isomers (*read off the terminal*; the failure lists are in
  `validation/2026-10_druglike10k/druglike_10k/db2_verification.csv`). The RDKit db2 fallback built 36 of the
  314 (see 9.2), consistent with these being impossible.
- **Status.** measured.

### 7.3 Escalate before calling an isomer impossible

- **Decision.** An isomer is rejected only after all six tries fail: ETKDGv3
  with 10 iterations at seeds 0xC0FFEE and 0xBEEF; random coordinates with 200
  iterations (seed 7); random coordinates with 10 iterations at seeds 1, 2, 3.
  Verdicts are computed up front in a fork pool.
- **Governs.** `EMBED_TRIES`; `--jobs` (default `SLURM_CPUS_PER_TASK`).
- **Commits.** `1d9d8af`.
- **Evidence.** v6 run: of 300 rejected isomers from the 245 parents lost at
  names, 9 embedded with more effort (fused cyclopropanes, bicyclo[2.1.1]);
  the escalation rescues 8 of the 9, and still rejects the impossible cage
  isomers (40-isomer cage set keeps the same 5) (*read off the terminal*). Uncapped
  ETKDG was left out: it rescued fewer and spends ~20 s on each impossible
  isomer. Cost ~8 s CPU per rejection, ~800 rejections per 10k parents.
- **Status.** measured.

### 7.4 Never drop a protomer's last isomer on the 3D check

- **Decision.** If every isomer of a protomer fails the check, all are kept and
  db2_converter decides. The check only removes impossible isomers where a
  buildable alternative exists.
- **Governs.** guard in `assign_names.main`; counted as
  `protomers_none_embed_kept`.
- **Commits.** `29a0230`.
- **Evidence.** v6 built 9,698 parents against v1's 9,771; all 100 parents
  that built in v1 and not in v6 were lost to `stereo_cannot_embed`
  (`validation/2026-10_druglike10k/druglike_10k_v2/regressed.txt`), mostly bicyclo[1.1.1]pentanes (whose
  bridgehead labels are not real stereocentres; RDKit rejects every
  labelling) and fused cyclopropanes that conformator builds. With the guard
  (v6b; 357 impossible isomers dropped, 310 protomers kept whole,
  `validation/2026-10_druglike10k/druglike_10k_v2/library_v6b/stats.json`): 9,797 parents built, 1 parent built in v1 and not in v6b
  (`validation/2026-10_druglike10k/druglike_10k_v2/db2_verification_v6b.csv`).
- **Alternatives.** More embedding effort: rejected; these fail at any effort in
  RDKit while conformator builds them.
- **Status.** measured.

---

## 8. Naming contract

### 8.1 One 16-character name per protomer and stereoisomer

- **Decision.** Every built molecule is named `C<12 base id><protomer><stereo>C`,
  exactly 16 characters, base62 fields. The base id is the ZINC id's last 12
  characters verbatim (unpadded ZINC ids zero-padded first), or `Z` + 11 base62
  of the sha256 of any other name; ZINC-shaped but malformed ids are rejected.
  One structure keeps one name (first by ZINC id); a second parent reaching it
  is dropped as `duplicate_structure`.
- **Governs.** `run_qupkake/llb_ids.py`, `assign_names.py` (asserts invariants
  1–3 before writing), `llb validate-names` (DRAP's `stage4/validate_names.py`).
- **Commits.** `5e728fd`, `3f7e5ff` (hashed base ids for non-ZINC names).
- **Evidence.** The contract and the failure each invariant prevents are
  specified in DRAP's
  [NAMING_CONTRACT.md](https://github.com/Garn0123/drap/blob/main/NAMING_CONTRACT.md)
  (top level of the DRAP repo). db2's M line stores `name[-16:]`
  (mol2db2/hierarchy.py); `--mergeiso` truncates at the first dot; db2_converter
  uses the name unquoted in `rm -r` with `shell=True`. `Z` cannot start a ZINC
  base id (ZINC22's first character is the heavy-atom count in base62, ≤ `N`;
  ZINC20's is a digit). v6b: all invariants hold on 19,102 names (v6: 9,771
  base ids, at most 64 names per base id).
- **Alternatives.** Readable suffixes (`_p1s0`, `_pH7p4`): kept as non-default
  name styles for debugging; rejected for libraries because db2 truncation
  loses the front of the name.
- **Status.** judgment (design), enforced mechanically.

---

## 9. db2 build

### 9.1 Restrained MMFF (ε = 4r, 0.3 Å flat-bottom, k = 100) replaces db2_converter's plain MMFF

- **Decision.** Conformers from conformator are relaxed with MMFF94s under a
  distance-dependent dielectric (ε = 4r; RDKit dielectric model 2) with each
  heavy atom free within 0.3 Å of its input position and a 100 kcal mol⁻¹ Å⁻²
  spring beyond, up to 200 iterations. A conformer RDKit cannot parameterize is
  kept unrelaxed and counted. db2_converter is not modified: the launcher
  replaces `pipeline.mmffopt` in-process before `build_ligand.main()` runs.
- **Governs.** `run_qupkake/build_ligand_relaxed.py`; `LLB_MMFF_VARIANT=eps4r_restrained`,
  `LLB_MMFF_MAX_DISPL=0.3`, `LLB_MMFF_K=100`, `LLB_MMFF_MAX_ITERS=200`;
  `submit_db2.sh --mmff restrained` (default; `on` = db2_converter's step,
  `off` = none).
- **Commits.** `81e3076` (`llb mmff-test`), `51a18b4` (default), `6fb4496`
  (logging).
- **Evidence.** Owen's set with no relaxation: conformator left starting
  geometry in some molecules (an N-H 1.149 Å from a C-H in every conformer),
  so every db2 set was flagged broken (98 of 429 hierarchies all-broken). On
  2,735 conformers from three such families (`llb mmff-test`, *read off the
  terminal*; input 2,011 under the 1.70 Å H…H rule, 469 N+…O− contacts < 3.2 Å):

  | variant | still < 1.70 Å | heavy-atom RMSD | Rg change | N+…O− contacts |
  | --- | --- | --- | --- | --- |
  | plain (db2_converter) | 0 | 1.83 Å | −0.5% | 1,405 (3×) |
  | ε = 4r | 0 | 1.26 | +3.6% | 750 |
  | electrostatics off | 0 | 1.24 | +3.5% | 615 |
  | hydrogens only | 373 | 0 | 0 | 469 |
  | restrained 0.3 Å | 9 | 0.28 | +2.6% | 537 |
  | **ε = 4r + restrained** | **3** | **0.27** | **+2.7%** | **486 (+4%)** |

  Plain MMFF clears the clashes but folds charged molecules in vacuum
  (ammonium onto carboxylate). After the switch, Owen's rebuild had 2 of 554
  hierarchies all-broken (from 98 of 429), and the bad conformer's closest
  contact went 1.15 → 1.75–1.87 Å; checked through the pinned db2_converter
  on the same conformer, 1.15 → 1.80 Å, bond table unchanged. druglike_10k:
  0 charge mismatches in 21,574 (v1) and 18,620 (v6b) builds; 1.6–1.7% of
  names under 25 conformers.
- **Alternatives.** As in the table: plain (folds), hydrogens only (leaves 373
  clashes, heavy atoms crowded too), ε = 4r or no electrostatics alone (still
  1.2 Å of movement and 30–60% more contacts), restraint alone (9 clashes, 15%
  more contacts). The 0.3 Å / k = 100 values were not scanned; `--max-displ
  0.5` was the fallback if 0.3 had not cleared the clashes.
- **Status.** measured (Owen's set; not re-measured on drug-like molecules,
  though their build statistics are clean).

### 9.2 No RDKit conformer fallback

- **Decision.** Names conformator cannot build are not retried with RDKit
  sampling.
- **Governs.** nothing in the default path; `verify-db2` still accepts
  `--outputpath` more than once (`ad4ed1e`) for a hand-run second build.
- **Commits.** `ad4ed1e` (multi-path verify); fallback steps removed from the
  runbooks.
- **Evidence.** RDKit did rescue one of Owen's molecules (owen_13: 587 of 600
  conformers, 544 s, against conformator's 0). On druglike_10k v1, the 314
  `1generate` names rebuilt with `--method rdkit` (*read off the terminal*;
  `llb timing db2_fallback`): 36 built (11.5%), 277 failed
  `3chemistrycheck`; median 133 s per name, max 4,950 s; 19% under 25
  conformers; parents lost to `1generate` moved columns (135 → 1) but parents
  lost overall did not fall (186 `3chemistrycheck` after). 285 of the 314 were
  stereoisomers that cannot exist (7.2): conformator correctly refused them,
  RDKit built the possible shape, which then failed the stereo check. The
  embedding check removes them upstream instead.
- **Alternatives.** RDKit for everything: ~100× slower (9–13 min per protomer on
  Owen's set). RDKit only for low-count builds (< 25 conformers): the planned
  threshold; dropped once the failures were shown to be stereo, not sampling.
- **Status.** measured.

---

## 10. Invertomers

### 10.1 Build and dock both N-H+ invertomers; collapse on scores

- **Decision.** A protonated amine with three different substituents is a
  stereocentre in a db2 (conformers never invert it), so both configurations
  are named, built and docked. `library.tsv` links each to its group's first
  member (`invertomer_of`). After docking, scores are collapsed to one row per
  species: per target, the best score in the group (lower is better for DOCK
  Grid_Score), carried under the first member's name, with a provenance table
  of winner and margin.
- **Governs.** `assign_names.invertomer_key` (configuration of every N-H+ erased;
  quaternary N+ left alone, it cannot invert); `llb collapse-invertomers`
  (`run_qupkake/collapse_invertomers.py`); DRAP `stage4/collapse_invertomers.py`
  (streaming port with a parity test against LLB's script), run first in stage
  4 so the row axis is fixed after collapse; `join_library`/`build_matrix`
  refuse an uncollapsed table.
- **Commits.** `5e728fd` (`invertomer_of`), `f6578cf` (LLB collapse); DRAP
  `1e2ea63`.
- **Evidence.** v1: 1,828 of 22,613 names (8.1%) are invertomers of another row
  (`validation/2026-10_druglike10k/druglike_10k/library/stats.json`); v6b:
  1,471 of 19,102 (7.7%; `validation/2026-10_druglike10k/druglike_10k_v2/library_v6b/stats.json`). The two are different 3D shapes for
  docking (a piperidinium substituent axial or equatorial; the N-H donor
  pointing one way or the other); in water they interconvert (deprotonate,
  invert, reprotonate), so they are one species, and a binding site selects
  one. v6: 24 of 179 `3chemistrycheck` failures were invertomers where
  conformator built the other one, which shows the geometries are not
  interchangeable in the build (*read off the terminal*).
- **Alternatives.** *Collapse before building* (keep one invertomer): saves
  ~8% of names, ~200 CPU-h of db2 per 100k parents and 8% of all docking; loses
  the pose the protein may prefer. Rejected for now: the saving is fixed and
  modest, the coverage loss unknown. `collapse-invertomers` reports how often a
  non-first invertomer wins and by how much; if it rarely matters, collapsing
  upstream becomes the cheaper option. *Treat invertomers as separate rows in
  the matrix:* rejected, near-duplicate rows for one species.
- **Status.** judgment (pending the first docking run's non-first-win rate).

---

## Open items a reviewer may ask about

- **No broad experimental pKa benchmark.** Triage agreement is model vs model
  (3.1); the only measured-value check is the 53-molecule aromatic-N set.
- **pH list, margin, min population, max states** (3.3, 3.4) are judgment.
- **Route window measured on H18–H22 only** (3.2).
- **Coupling under-corrects conjugated paths and ignores opposite charges** (5.1).
- **One name per run differs in stereo perception** between the naming RDKit and
  db2_converter's RDKit (`Cm50000037lNj10C`, a protonated bicyclic amine); it
  is excluded by hand before the build.
- **Restrained MMFF parameters** were chosen on Owen's set, not scanned (9.1).
- **Evidence gaps.** Numbers marked *read off the terminal* have no file in
  `validation/2026-10_druglike10k/`; the missing ones that matter most are `compare_c/summary.json`
  (the coupled window curve), v1's `llb timing`
  output, the `db2_fallback` timing, `POOLED.done`,
  `library_v6/library.stereo_crashed.tsv`, the merged charge spreads, the 735
  aromatic-N swaps and the per-parent loss breakdowns.

## Inconsistencies found while writing this record (not fixed here)

Code and docs were left unchanged; each is a follow-up.

1. **Williams page numbers** in `qupkake_protomers.py` (coupling and
   `IMPLAUSIBLE_SITES` comments), `pka_triage/TRIAGE.md` and the `source` column
   of `pka_triage/reference/aromatic_n_bases.tsv` are mostly wrong by one or
   more pages (the values are right); correct pages in section 0.
2. **Meldrum's acid** is given as "about 4.8" in `TRIAGE.md`; the primary water
   value is 4.97.
3. **Diazo justification**: the `azide_base` comment and `TRIAGE.md` say diazo
   groups are not protonated at pH 6–8; methyldiazonium's pKa is 10.0, so the
   reason is reactivity, not basicity (4.6).
4. **Stale "independent-site" descriptions** after coupling became the default:
   the `qupkake_protomers.py` module docstring ("This ignores site-site
   coupling") and the `triage.py` docstring ("QupKake as we run it is ALSO
   independent-site, so routing these is a flag for review rather than a
   fix").
5. **Stale branch/status text**: `TRIAGE.md` opens with "Branch
   `feat/pka-triage`. Status 2026-10-08 ... no real QupKake or Great Lakes run
   yet", and its "Next steps" are done; `README.md` still says pka_triage is on
   branch `feat/pka-triage`.
6. **`--keep-amide-sites`** turns off the whole site filter (all v1–v6 rules), not
   only the amide rules its name and help text describe; the
   `IMPLAUSIBLE_SITES` comment says "turns both off".
7. **Comment order in `IMPLAUSIBLE_SITES`**: the carbon-acid comment sits above
   the amidine entry, and the carbon-acid entry is last.
8. **Dates**: code comments and `TRIAGE.md` date v5, v6, aromatic-N and the 3D
   check fixes 2026-10-16, while their commits are dated 2026-10-15 (UTC vs
   local time); harmless but confusing in a paper trail.
9. **Two citations for Me2O −3.8**: the `sulfur_base` comment cites Ripin &
   Evans, `TRIAGE.md` cites Reich H2O; both tables give it.
