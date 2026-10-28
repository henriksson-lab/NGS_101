# eTAM-seq: enzyme-assisted adenosine deamination

Xiao et al., *Nature Biotechnology* (2023), [doi:10.1038/s41587-022-01587-6](https://doi.org/10.1038/s41587-022-01587-6).

The current scope is the transcriptome-wide NGS library, with FTO−/FTO+/IVT controls. Site-specific Sanger assays are not counted as another NGS protocol. This is a bounded source pass and partial model; final published procedural identity, complete copying/library geometry and calibration remain unresolved.

## Sources and actual reading

- `s41587-022-01587-6_PMC10625715.xml.txt` and corresponding PDF are an accepted author manuscript (manifest author-manuscript true). Inspected spike-ins 173–175; IVT preparation 177–181; enzyme assays 183–195; NGS library 203–211; distinct site-specific assay 213–217; Fig. 1 legend 843–845. Reference R14 was resolved in the same original.
- Publisher `s41587-022-01587-6_41587_2022_1587_MOESM1_ESM.pdf.txt`: identity/front matter; actual Supplementary Notes 1–2, lines 90–128, describing the temperature schedule and IVT accessibility model. Other statistical notes, full tables and image panels were not systematically inspected. The NIHMS counterpart is cached separately; complete version equivalence was not presumed.
- The complete-bundle oligo scraper was run; immediate library oligos were verified at the manuscript Methods passages. Large concatenated hits from synthetic probe tables are not a verified transcriptome-library construct.

## Covalent RNA, complementary DNA and reported calls

🟢 TadA8.20 deaminates susceptible RNA A to **inosine**, while m6A is relatively resistant. RNA I templates complementary first-strand DNA C; reference-oriented sequence reporting gives G. The local `SitePath` preserves all three identities separately. It does not rewrite RNA as DNA G, remove m6A chemically during RT, or treat every persistent A as a methylation call without controls.

The constructor models selected reacted A→I and protected m6A paths. Unconverted A, inaccessible sites, background and FTO residual fractions remain unassigned. The publisher SI explicitly frames IVT accessibility matching as an assumption for statistical correction, not a known molecular identity. No methylation stoichiometry is estimated by the current code.

## Source-defined NGS preparation

1. 🟢 Poly(A)+ target RNA and separately prepared IVT control RNA undergo oligo(dT)/RNase H treatment, zinc fragmentation at 70 °C for 5 min, EDTA quench, PNK treatment and cleanup. Spike-in probes are added at the source's 2% mass fraction.
2. 🟢 Attach the preadenylated 3′-biotin DNA adapter with truncated KQ Rnl2; digest excess adapters with deadenylase and RecJf. The target preparation is divided into FTO− and FTO+ halves; IVT is a separate control. Bind the three preparations to C1 streptavidin beads.
3. 🟢 Treat the FTO+ group with FTO at 37 °C for 1 h in the reported buffer, then wash. No complete erasure or identical accessibility is assumed.
4. 🟢 Perform TadA8.20 deamination at 53 °C for 1 h, then two fresh 44 °C 1-h rounds. Source supernatants are drained and beads washed. The **site-specific assay** uses a different 37 °C schedule and is not substituted here.
5. 🟢 Maxima H− RT acts at 50 °C for 1 h. Heat/SDS elution and cleanup precede cDNA 3′ adapter ligation, another cleanup and NEBNext Ultra II Q5/E6440S unique-dual-index amplification (typically 10–11 cycles), followed by AMPure purification and sequencing.

The graph begins at prepared bead-bound target RNA and stops after the selected-site RNA–cDNA hybrid. It retains the RNA template at RT and marks soluble material as waste only at documented physical washes. Upstream fragment/adapter cleanup, IVT synthesis, calibration probes and full downstream PCR geometry are not claimed complete.

## Verified immediate library oligos

| Component | Accepted-manuscript sequence/modifications |
|---|---|
| RNA 3′ adapter | `/5rApp/AGATCGGAAGAGCGTCGTG/3Bio/` |
| RT primer | `ACACGACGCTCTTCCGATCT` |
| cDNA adapter | `/5Phos/NNNNNNAGATCGGAAGAGCACACGTCTG/3SpC3/` |

The C3 block is non-nucleotide; it is not deleted to create an ordinary DNA string. The local cDNA-adapter inventory carries a structured random-identifier feature; UMI semantics remain an explicit inference pending full processing review. Actual E6440S index members, sample-sheet values, acquired read cycles and exact run-primer geometry are unassigned, not declared absent.

The IVT target is modification-free RNA produced through oligo(dT), a biotin/T7 template-switch oligo, second-strand synthesis and T7 IVT (Methods 179–181). This does not imply that subsequently added methylated calibration probes are modification-free. The local IVT-site paths describe the target component; calibration material is separate and not yet modelled.

## Component ancestry and remaining work

**Explicit library-strategy predecessor:** `m6a-sac-seq__10.1038+s41587-022-01243-z` → this eTAM-seq directory. Methods 205 and 211 cite ref. 14, which resolves to Hu 2022 m6A-SAC-seq. Scope: poly(A) removal/end repair/RNA-adapter and cDNA-adapter library subroute, with local changes. It does not establish inherited allylation/cyclization chemistry or equal adapter sequences, enzyme conditions, controls or readout. Each remains declared locally; no recipe constructor is imported.

Proposed metadata: RNA; bulk poly(A)+ input; enzyme-assisted deamination; affinity-tethered treatment; RNA/cDNA adapter ligation; RT/PCR; IVT/FTO/calibration controls; dual-index Illumina library; published instance. Keep `notes`/`pending_source_review`. Deferred: final VOR comparison, remaining supplements, full IVT/spike-in states, exact cDNA adapter/PCR products/index assignments/acquisition, and statistical correction. No confirmed missing-information claim.
