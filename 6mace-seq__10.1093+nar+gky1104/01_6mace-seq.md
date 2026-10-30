# 6mACE-seq (DNA)

Defining DOI: https://doi.org/10.1093/nar/gky1104.

## Sources inspected

🟢 Published main Methods82–84; published TableS1 rows3–7; Supplementary Methods/legends relevant headings and FigS1–S2. Original cache: `/husky/henriksson/NGS101_papers/sources/6mace-seq__10.1093+nar+gky1104/nar_gky1104_PMC6294517.xml.txt`. Mandatory whole-bundle oligo scraper run; source hits checked against original text/table, not treated as independent evidence.

## Molecular scope

DNA-specific protection branch. Antibody protection is incomplete/selective in practice; no perfect specificity, complete digestion or complete protein removal is asserted. Exact assembled end transitions and full TruSeq PCR/read geometry are deferred. Lambda and RecJf are separate reactions with washes; RNA XRN1 is not substituted.

1. 🟢 Covaris near150bp; RNaseA/T1 and AMPure source workup. Molecular representation: fragmented RNA-depleted genomic DNA; exact residual RNA unquantified.

2. 🟢 Anti6mA incubation then254nmUV6x0.15J/cm2. Molecular representation: genomic DNA population with UV antibody adducts at selected sites and uncrosslinked background.

3. 🟢 ProteinA bead binding and source sequential buffer washes. Molecular representation: bead-associated DNA-antibody complexes and possible nonspecific DNA.

4. 🟢 T4DNA polymerase12C30min then source bead washes. Molecular representation: end-repaired captured DNA; exact strand architecture unassembled.

5. 🟢 T4PNK37C30min then source bead washes. Molecular representation: phosphorylated captured DNA class.

6. 🟢 PG160/PG161 reverse adapter ligation25C2h and washes. Molecular representation: reverse-adapter-associated DNA; exact sealed junction geometry unassembled.

7. 🟢 Phi29 fill and source bead washes. Molecular representation: Phi29 fill products on captured DNA; no exact end product inferred.

8. 🟢 Second T4PNK then source bead washes. Molecular representation: rephosphorylated captured DNA population.

9. 🟢 Lambda37C30min. Molecular representation: antibody-protected DNA and Lambda-generated digestion debris, all still in tube.

10. 🟢 Actual source bead washes after Lambda. Molecular representation: washed protected DNA; retained/carryover debris unquantified.

11. 🟢 RecJf37C30min. Molecular representation: protected DNA and RecJf-generated debris, retained until source washes.

12. 🟢 Actual source bead washes after RecJf. Molecular representation: washed protected DNA fraction; no perfect depletion assumed.

13. 🟢 SDS elution and proteinaseK overnight65C followed ethanol precipitation. Molecular representation: recovered protected ssDNA with possible residual peptide/adduct material; no perfect chemical reversal asserted.

## Readout and remaining work

Source MiSeq/NextSeq single-end; read-start pileups report protected6mA positions computationally, not a chemical conversion of adenine.

Full molecular end-transition assembly and source TruSeq PCR primer/run geometry remain unreviewed. No sequence-absence claim.

This capped pass is incomplete and makes no accepted missing-information claim. Local preview has no finished-publication marker.

## Ordering oligos and post-exonuclease route

🟢 Published TableS1 rows3–7 (5′→3′):

| Name | Verbatim sequence | Use |
|---|---|---|
| PG160 | TGGAATTCTCGGGTGCC | Reverse-adapter strand |
| PG161 | /5Phos/CCTTGGCACCCGAGAATTCCA | Reverse-adapter partner |
| PG162 | GCCTTGGCACCCGAGAATTCCA | Post-release Phi29 primer |
| PG163 | GTTCAGAGTTCTACAGTCCGACGATC | Forward-adapter strand |
| PG164 | GATCGTCGGACTGTAGAACTCTGAACCACCA | Forward-adapter partner |

🟢 After ethanol precipitation, source Methods84 anneals PG162 to the reverse adapter and extends with Phi29 at30°C30min; forward adapter is ligated with T4ligase25°C2h, followed by1.7×AMPure, Phusion/TruSeq PCR (human18cycles, bacterial13cycles), then0.8×AMPure dimer removal. The graph now includes these source-known operations with explicitly unassembled end and PCR geometry; no exact library or primer-binding claim.

The bait/control oligos elsewhere in TableS1 belong to DNA pulldown, demethylation, anisotropy or validation assays and are not substituted for library adapters. The source scraper can misparse embedded `6mA` as an ordinary-base substring; those hits do not establish adapter identity.

## Component ancestry

Explicit DNA6mACE→RNA m6ACE development: published RNA m6ACE Results39 / Figure1a / reference31, associated source `s41467-019-13561-z_PMC6904561.xml.txt`. The authors state they capitalised on the recent DNA6mA method. Scope: antibody photo-crosslink and exonuclease-protected read-start strategy. DNA Lambda/RecJf/Phi29/adapters are not inherited by the RNA/XRN1 preparation.
