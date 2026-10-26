# Drug-like 10k validation (October 2026)

> **Code at tag `llb-full-2026-10`.** The construction code this record cites
> (`pka_triage/`, `run_qupkake/`, `pipeline/`, `environments/`) was removed on
> 2026-10-26, when ligand construction moved to MolSanitizer. Paths and commits
> here refer to that tag (`git checkout llb-full-2026-10`). Kept in the tree
> because `lib_h18_22` was built and docked with that code.

Summary files from the two validation runs behind `DECISIONS.md` and
`pka_triage/TRIAGE.md`. Same 10,000 parents both times: a uniform draw from ZINC22
2D H18-H22, proportional to bin size (`llb sample --heavy 18-22 --total 10000 --local
...`). Run directories on Great Lakes:
`/scratch/maom_root/maom99/johnbick/2026_10_08_build_small_mols/druglike_10k{,_v2}`
(scratch: purged after 60 days unused; these copies are the record).

## druglike_10k/ (v1: before the October 14-16 fixes)

Code: `feat/pka-triage` at site filter v3, independent sites (no coupling), no
aromatic-N hybrid, no tautomer step, no 3D stereo check. Restrained MMFF.

| file | what it shows |
| --- | --- |
| `triage/stats.json` | triage routed 36.1% (route window 1.5, site filter v3) |
| `compare/` | MolGpKa vs QupKake on 9,772 parents after the QupKake embed fix (e547c3b): dominant state agrees 83.4 / 85.3 / 86.0% at pH 6.4 / 7.4 / 8.4; 207 settled-but-disagree at window 1.5; `window.png` the routed fraction and disagreement per window, `sites.tsv` every site pair with class and context |
| `parents/rejected.tsv` | no inputs rejected |
| `library/stats.json` | 9,997 parents -> 22,613 names (2.26 per parent) |
| `db2_verification.csv` | 21,574 of 22,613 names built; 477 `3chemistrycheck`, 314 `1generate`, 248 missing (stereo-enumeration crashes + 1 left out); 0 charge mismatches |
| `ref_qk/reference_check.tsv` | aromatic-N reference set (`llb pka-reference`): QupKake vs MolGpKa vs measured, alone and tethered to a guanidine or dimethylamine (commit 3c5202c) |

## druglike_10k_v2/ (v6b: with every fix through 29a0230)

Code: site filter v6, coupling v1, aromatic-N hybrid v2, imino -> amino tautomers,
3D stereo check with escalation and the last-isomer guard. QupKake ran on the routed
parents only.

| file | what it shows |
| --- | --- |
| `triage_v6/stats.json` | triage routed 35.8% (site filter 6, coupling 1) |
| `library_v6b/stats.json` | 9,997 parents -> 19,102 names (1.91 per parent); 357 impossible stereoisomers dropped; 310 protomers kept whole (no isomer embedded) |
| `library_v6b/dropped.tsv` | every dropped isomer and why |
| `db2_verification_v6b.csv` | 18,620 of 19,102 names built (verify over the reused v6 db2s + the v6b build); 9,797 of 10,000 parents built; 0 charge mismatches |
| `timing_v6.txt` | cost: triage 10.7 ms/parent; QupKake 31.8 s/parent; db2 median 45.8 s/name, 268 CPU-h for 18,677 names |
| `regressed.txt` | the 100 parents v1 built that v6 (before the last-isomer guard, 29a0230) lost: all `stereo_cannot_embed`, mostly bicyclo[1.1.1]pentanes and fused cyclopropanes; the reason for the guard |

Numbers quoted elsewhere that are not in these files (they were read off the
terminal at the time): the QupKake +2/+3/+4 protomer counts before and after the
fixes (v1: 1,824 / 257 / 12; v2 merge at site filter v5 and aromatic-N v1: 1,050 / 75 / 1), the 735 parents with an aromatic-N
swap, and the per-parent loss breakdowns (v1: 135 `1generate`, 64
`3chemistrycheck`, 27 missing; v6b: 135, 65). They can be recomputed from the run
directories while those exist.
