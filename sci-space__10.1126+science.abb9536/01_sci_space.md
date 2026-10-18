# sci-Space — research notes

## Active native migration, 2026-10-18

Source scope: Published defining JATS science.abb9536_PMC9118175.xml.txt overview and linked supplemental-method identities; published Science Supplement caltech_abb9536-Srivatsan-SM.pdf.txt pp1–5 (Overview24–36; grid37–62; transfer117–144; native library145–232; computational filtering233–249). Exact cited gridguide protocolsio_creating-sci-space-grids-for-spatial-barcoding-5qpvoy58bg4o.pdf.txt first-page/printedidentity DOI10.17504/protocols.io.bm64k9gw,Feb13,2021,lastmodifiedFeb13,2021,ID42940; Step7.1 119–125 hash scaffold/order-file links. Normalized cache/URL is not revisionidentity; printed original matches cited DOI.

Own17step graph independently retains transcriptRNA,spotDNA,sectorDNA,sectionDNA and fourcopiedpopulations. SourcegridDNAprefix+barcode10+B(C/G/T)+A32 is verified; copied Bcomplement remains A/C/Grole, not fabricated literalbase. SelectedSlides1–6 sort200nuclei/PCRwell follows pooledindexedRT; independentSlides7–14dilution250cohort is descriptive only, not graphsplit. Exact RTprimer order/placement unassigned; typed UMI8/RTbarcode10 window layout is explicitly inferred from sourceR1content, no defaultfullprimer. Actual ssRTcopies→native dsDNA→customn7tagged/gapped states→75°Cextension(with actualPCRprimers/mix) enter PCR; originalRNA/DNAhashes/RTcopies/nonamplifying products retain unresolved molecularfates. PhysicalPCRpool/Zymo/0.8×SPRI explicit. SourcePCRindex10+10 and run18/10/10/52 are independently stated; unknownprimers/nativeends not canonical shells. CitedgridDOIextra_sources locator was corrected by parent to exactcachedPDFURL after printedidentity verification; no neworiginal/fetch. Historical guide-versus-study condition notices remain scoped history, not newly reread sources.

Mandatory scraper: 3hits2files;22MBDataFileS3 skipped (analytical gene-model data, not inspected as chemistry). Cited gridscaffold and native SMwidths checked in originals.

Deferred: Cited native two-level sci-RNAseq guide bquamwse fullRT/customn7/orderedindexedPCR species and revision applicability; GridguideStep7.1linked oligo-order files; member choices/attachment identity not inspected aftercap; Full section-hash sequence/tail, native RTwindow ordering/extent, exact copied/taggedend geometry; IndependentSlides7–14dilution cohort and other validation experiments outside selected branch.

This is a capped implementation/source pass, not a complete chemistry review or accepted missing-information finding. Historical source passes below retain their original scope; sources were frozen before implementation and no source expansion followed the cap.


Srivatsan et al., *Science* (2021),
[doi:10.1126/science.abb9536](https://doi.org/10.1126/science.abb9536).

## Sources read

- 🟢 Defining article (PMCID: PMC9118175, author manuscript JATS).
- 🟢 Published Science Supplementary Materials (`abb9536-Srivatsan-SM.pdf`, "Published 2 July
  2021, Science 373, 111"), obtained from the CaltechAUTHORS repository record
  [zd6v5-5nb93](https://authors.library.caltech.edu/records/zd6v5-5nb93) because the PMC
  author-manuscript copy is behind PMC's download gate. Materials and Methods pp. 1–5 read
  (extracted lines 25–245). Data Files S1 (cluster marker genes) and S3 (trajectory
  model terms) were also fetched and are gene tables, not chemistry.
- 🟢 Author-linked protocols.io guides for [grid preparation](https://doi.org/10.17504/protocols.io.bm64k9gw), [spatial transfer and imaging](https://doi.org/10.17504/protocols.io.btw5npg6), and [two-level sci-RNA-seq library preparation](https://doi.org/10.17504/protocols.io.bquamwse). Their fetched PDFs were checked with the oligo scraper and read at the cited steps.
- Not read: Cao et al. 2017 two-level sci-RNA-seq chemistry,
  [doi:10.1126/science.aam8940](https://doi.org/10.1126/science.aam8940), cited by the
  Supplementary Materials (ref. 39) for the indexed oligo(dT) and indexed P5/P7 primers.

## What is spatially encoded

🟢 The agarose-coated slides carry a printed 18 mm × 18 mm grid of 7,056 uniquely
barcoded spots. Hash oligos are single-stranded DNA with poly(A) tails. About 5% of the
spots also contain SYBR Green, forming fluorescent fiducials used to register the ideal
grid to the tissue image. The reported mean spot radius was 73.2 µm and center spacing
222 µm.

🟢 A separate hash identifies the tissue section. It permits sections to be pooled,
but is not itself the x–y coordinate. 🟡 The linked transfer guide calls for a unique
slide-specific oligo in lysis buffer (step 1.3) without printing its sequence or tail;
its poly(A) structure cannot yet be verified from the readable sources.

🟢 The grid-preparation guide, step 7.1, prints the spot/sector hash ordering structure
verbatim as `5′-GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG[10bp-barcode sequence]-BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA-3′`,
where `B` is C, G or T. The 10-base barcode differs among ordered oligos. The guide
points to `all_plates_for_ordering.xlsx` for the complete panel, but that attachment is
not included in the fetched PDF. The printed structure is also shown on the public page;
the individual plate sequences remain an open source item, not a confirmed absence.

## Construct-changing steps

1. 🟢 Cryosection fresh-frozen tissue, permeabilize it with the section-hash oligo and
   juxtapose it with the printed agarose grid. Grid hashes diffuse into nearby nuclei.
2. 🟢 Image the sandwich during transfer so fluorescent fiducials preserve the transform
   between barcode coordinates and tissue anatomy.
3. 🟢 Extract nuclei, fix them and perform sci-RNA-seq. The barcoded oligo(dT) RT primer
   captures endogenous poly(A) RNA and the printed hashes in the same nucleus. The
   section hash is also recovered with the nucleus, but its exact capture-ready tail
   structure is not established by the readable guide.
4. 🟢 Subsequent pooling, redistribution, second-strand synthesis, Tn5 fragmentation and
   indexed PCR create separate transcript and hash molecules carrying the same nucleus
   identity.

🟢 The author-linked two-level sci-RNA-seq guide prints the terminal amplification
primers as `5′-AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′`
and `5′-CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG-3′` (step 51, extracted
lines 592–601; the P5 primer's final `TCT-3′` continues after the page 13/14 break,
which an earlier version of this note dropped). It describes an n7-loaded Tn5 or dual-loaded Illumina Tn5
before PCR (steps 37–39). These are linked predecessor-guide statements; the main
sci-Space paper does not itself spell out every library-primer sequence.

## Library boundary

🟢 A spatial hash is not covalently fused to a cell’s transcripts. Association is by the
shared sci-RNA-seq combinatorial indices: transcript reads define expression and hash
reads define the most enriched printed coordinate and section. The schematic therefore
shows a representative hash molecule and separately states the transcript branch rather
than drawing a chemically false combined molecule.

🟡 The public diagram uses the established two-level sci-RNA-seq outer primer geometry.
Individual barcode identities are represented as semantic regions because the experiment
uses large panels; they are not synthesized as strings of repeated placeholder bases.

## Bounded pass, 2026-10-17 (historical; superseded where noted below)

Read: PMC9118175 JATS main text lines 93–97 (four-step protocol, grid, slide oligo)
and references 37–38; the two-level sci-RNA-seq guide (doi:10.17504/protocols.io.bquamwse,
fetched PDF dated 13 Feb 2021, printed DOI matches the article's cited v1; no version
label printed) steps 15–20, 37–59; the transfer guide steps 1.3 and 8–10; grid guide
step 7.1 line 122 re-confirmed for the hash structure. The Supplemental Materials and
Methods PDF and Files S1–S3 remain behind the PMC download gate
(<https://pmc.ncbi.nlm.nih.gov/articles/instance/9118175/bin/NIHMS1801003-supplement-Supplemental_Material_and_Methods.pdf>).

Model corrections:

- 🟢 Second strand uses the NEBNext Ultra II Non-Directional RNA Second Strand
  Synthesis module (16 °C, 150 min); the earlier separate "second-strand primer"
  input was a model assumption and is removed.
- 🟢 Tn5 is in-house N7-loaded or dual-loaded Illumina Tn5 (55 °C, 5 min), followed by
  Zymo binding buffer and a per-well SPRI cleanup; PCR uses NEBNext High-Fidelity 2×
  Master Mix (72 °C 5 min, then 18–22 cycles) and pooled 1× AMPure cleanup.
- 🟢 The slide oligo enters with indexed cold lysis buffer (permeabilization), and grid
  hashes diffuse into nuclei; neither is modelled as hybridization. Scraping into 5%
  PFA with sonication is now a fixation step.
- 🟡 The final hash molecule now includes the oligo(dT) copy of the hash poly(A):
  the step-7.1 hash carries S7 and the mosaic end, so its RT copy already ends in
  the P7-primer site (derived; not stated in the read sources).
- Tagmentation by-products are not physically separated, so their fate is unresolved.
- Open: the guide recommends Read 1 18, Index 1 10, Index 2 10, Read 2 52 cycles;
  the sci-Space run cycles are not stated in the read sources and the gated
  supplementary methods were not read, so no acquired cycle count is assigned. The RT-primer sequence, index members,
  section-oligo sequence and the ordering workbook remain unread.

## Deferred follow-up pass, 2026-10-17 (supplementary methods read)

Read: Supplementary Materials pp. 1–5 (extracted lines 25–245: grid printing, hash
tests, transfer, library preparation, sequencing); scraper found no sequences in that
PDF. Corrections relative to the guide-based model above:

- 🟢 Printed wells hold 15 µM spot oligo and 2.5 µM sector oligo with 0.5% glycerol
  (± 5× SYBR Green) (p. 1). Slide hash: 5 µL of 10 µM in 500 µL permeabilization buffer
  (10 mM Tris-HCl pH 7.4, 10 mM NaCl, 3 mM MgCl2, 1% SUPERase-In, 0.1% IGEPAL) with DAPI (p. 3).
- 🟢 Embryo sections were scraped into 4% PFA, slides 7–14 sonicated, 15 min on ice,
  pooled in NSB and frozen (p. 3). The transfer guide and the SM hash-test section (p. 2)
  use 5% PFA; the model follows the embryo-experiment value and the page shows the difference.
- 🟢 Thawed nuclei were permeabilized (NSB + 0.25% Triton X-100, 3 min), sonicated 12 s and
  distributed at 5000 per RT well; 1 µL of 25 µM indexed oligo(dT) (guide: 10 µM), 55 °C 5 min,
  SuperScript IV ramp 10/20/30/40/50 °C × 2 min then 55 °C 15 min (guide ramp starts at 4 °C) (p. 4).
- 🟢 Second strand: NEB mRNA second strand synthesis buffer and enzyme, 16 °C 150 min, after
  sorting 200 nuclei/well (slides 1–6) or diluting to 250/well (slides 7–14) (p. 4).
- 🟢 Tagmentation with 0.01 µL custom n7-loaded Tn5 in 2× Nextera TD buffer, 55 °C 5 min;
  12 µL Zymo DNA binding buffer; 36 µL AMPure XP; 17 µL EB elution (p. 4).
- 🟢 PCR: NEBNext HF 2× master mix, 75 °C 3 min, 98 °C 30 s, 18 cycles (98 °C 10 s, 66 °C 30 s,
  72 °C 1 min), 72 °C 5 min; pooled, Zymo DNA Clean & Concentrator, 0.8× AMPure XP (p. 5).
- 🟢 Acquired run: NextSeq 500 high-output 75-cycle kit, onboard primers, Read 1 18 (8-nt
  UMI + 10-nt RT barcode), Index 1 10, Index 2 10, Read 2 52 (transcript or hash) (p. 5).
  These are now the page's acquired read lengths.

Still open: the RT-primer, Tn5 adapter and P5/P7 index member sequences (SM cites
Cao 2017, ref. 39, which was not read in this pass); the section-hash sequence/tail;
the `all_plates_for_ordering.xlsx` hash panel. These are open questions, not confirmed absences.
