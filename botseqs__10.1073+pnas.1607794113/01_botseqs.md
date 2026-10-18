# BotSeqS

**Source.** Hoang et al., *PNAS* (2016), [doi:10.1073/pnas.1607794113](https://doi.org/10.1073/pnas.1607794113), PMCID PMC5024639.

- 🟢 Libraries used the Illumina TruSeq DNA PCR-Free kit.
- 🟢 The bottleneck dilution is made after adapter ligation and before PCR, at about one haploid-genome equivalent.
- 🟢 For 500 ng input, the typical post-ligation dilution was 100,000-fold. The paper used 34 ng–1 µg input, 350-bp insert selection and 2 × 100-base HiSeq 2500 sequencing.
- 🟢 Read pairs sharing both genomic coordinates and orientation form a family; P5-derived Read 1 and P7-derived Read 2 distinguish Watson and Crick strands.
- 🟢 Native fragment endpoints are endogenous identifiers; BotSeqS adds no random UMI.
- 🟡 The paper says standard Illumina Y adapters but does not specify a dual-index arrangement. The schematic therefore shows the canonical single-index TruSeq molecule.

Dilution restricts the adapter-ligated originals entering PCR. Deep paired-end sequencing then supplies redundant Watson and Crick families for consensus.

## Bounded source pass 2026-10-17 (round breadth-first-07)

Read: PMC5024639 article HTML text lines 1243 (dilution rationale), 1246–1256 (library
generation, Watson/Crick definition, HiSeq 2500 rapid run 2 × 100, Fig. 1 legend) and
1460–1463 (Methods: DNA purification, shearing, Y-adapter ligation, dilution). The two
HTML copies on disk differ only in page boilerplate. The SI Appendix PDF and Dataset S1
XLSX were not obtainable: PMC returns a download gate, PNAS returns 403, and PMC Cloud
has no copy (404). Scraper: 0 oligo hits; the article prints no oligo sequences.

- 🟢 DNA came from Qiagen AllPrep, 34 ng–1 µg in 55 µL TE. It was sheared on a BioRuptor at
  high intensity, 15 s on / 90 s off, seven cycles at 3 °C. These settings are now in the graph.
- 🟢 Y-adapters were ligated with the TruSeq DNA PCR-Free kit using the standard low-input
  protocol, with selection for 350-bp inserts. The graph now has a size-selection step. 🟡 Its
  placement after end repair follows the October 2017 TruSeq DNA PCR-Free guide
  (1000000039279 v00), which postdates the paper; the 2016 kit revision was not inspected.
- 🟢 The dilution factor was chosen from a 10-fold dilution series assessed on MiSeq
  (SI Appendix Fig. S1). The typical dilution was 10^5-fold for 500 ng input.
- Model correction: the post-dilution PCR primers had been labelled "TruSeq … indexed
  reverse primer". The article names no PCR primers, polymerase or cycle number, and the
  index is carried by the ligated PCR-free adapter. The graph now marks these as
  identified only in the uninspected SI Appendix. The acquired reads are limited to the
  stated Read 1/Read 2 × 100; whether an index read was acquired is not stated.

Open: the SI Appendix (dilution, PCR, sequencing details) and Dataset S1; the 2016 TruSeq
DNA PCR-Free kit revision (adapter index length/arrangement, size-selection placement).
The drawn canonical library leaves the i7 window at an unassigned length
(`batch_ngs.UNASSIGNED_LENGTH`), because the 2016 adapter index length is not stated.

## Deferred follow-up 2026-10-17 (round deferred-followup-02)

Routes for the PNAS SI Appendix and Dataset S1, all unsuccessful: `get_sources.py --slugs`
(PMC bin URLs return a proof-of-work "Preparing to download" page, not solved; Europe PMC
fullTextXML HTTP 500 and supplementaryFiles "not open access"; PMC Cloud `PMC5024639.*`
404); pnas.org `doi/suppl/.../pnas.1607794113.sapp.pdf` and europepmc.org bin URL are
Cloudflare 403; Unpaywall/OpenAlex list only PMC and the publisher; a lab page found by
search (labs.pathology.jhu.edu) returned 403. Manual URL:
https://pmc.ncbi.nlm.nih.gov/articles/instance/5024639/bin/pnas.1607794113.sapp.pdf.
Illumina support URLs guessed for the 2016 TruSeq DNA PCR-Free guide 15036187 returned 404.

Read instead, as a **labelled lead**: US 12,359,250 B2 "Bottleneck sequencing" (Johns
Hopkins; PCT/US2017/015229; provisional 62/288,869, 29 Jan 2016), USPTO image PDF, OCR
twin. Example 1 ("Methods", col. 12–13) uses the same sample IDs and states that its
supplementary tables are those of Hoang et al. PNAS 113(35), so it describes this study.
The oligo bases were checked against the rendered page image, not only the OCR.

- 🟡 (patent) Shearing, 34 ng–1 µg input, TruSeq PCR-Free low-input ligation with 350-bp
  selection and 20 µL ligated volume agree with the PNAS Methods.
- 🟡 (patent) Dilution: five serial 10-fold dilutions, 2 µL into 18 µL TE with 1 ng/µL
  pBlueScript carrier; 2 µL of each used as PCR template; final factors 10^3–10^7.
- 🟡 (patent) PCR: HPLC-purified IDT primers with one 3′ phosphorothioate bond,
  TS-PCR Oligo1 5′-AATGATACGGCGACCACCGAG*A-3′ and TS-PCR Oligo2
  5′-CAAGCAGAAGACGGCATACGA*G-3′, 0.5 µM each; Q5 Hot Start 2x master mix, 50 µL;
  98 °C 30 s, 18 × (98 °C 10 s, 69 °C 30 s, 72 °C 30 s), 72 °C 2 min; 1.0x AMPure XP.
  Oligo1 is a 5′ prefix of P5 and Oligo2 of P7 (computed by the scraper and the PCR
  geometry), so the
  graph now computes both inward sites on the canonical TruSeq termini. The AMPure
  cleanup is not modelled.
- 🟡 (patent) MiSeq dilution test: ~5 M clusters, 2 × 50, three well-barcoded samples per
  lane; optimum 5–10 duplicates per molecule scaled to ½ HiSeq rapid lane.
- 🟡 (patent) A WGS library was amplified from 10 µL of the same ligation with TruSeq PCR
  Primer Cocktail; that branch is outside the BotSeqS schematic.

These values replace the "unidentified" PCR inputs as inferred, patent-sourced inputs;
they are not verified against the PNAS SI, and the page caveat says so.

Open: PNAS SI Appendix and Dataset S1 (manual); 2016 TruSeq DNA PCR-Free kit revision
(adapter index, size-selection placement); index read acquisition.

## Deferred follow-up 2026-10-17 (round deferred-followup-03): PNAS SI Appendix read

The maintainer hand-placed the PNAS SI Appendix (`pnas.1607794113.sapp.pdf`, first page
"Supporting Information Appendix / Hoang et al., PNAS 2016", running footer "SI Appendix
for BotSeqS") and Dataset S1 (`pnas.1607794113.sd01.xlsx`, sheets Table S1–S9; Table S1
"Human samples, bottleneck sequencing library preparation, and sequencing information").
Both match the defining article. Read: SI Materials and Methods lines 1–61 (dilution, PCR,
MiSeq run, WGS), Fig. S1 legend lines 368–380, and Dataset S1 Table S1 (44 BotSeqS
libraries, lines 1–53 of the `.txt` twin). Scraper: the SI prints the two TS-PCR primers
(line 19); the only other SI hit (line 1955) is a genomic sequence in an analysis listing.

- 🟢 Dilution (SI lines 4–17): the 20 µL of TruSeq adapter-ligated molecules is serially
  diluted, 2 µL into 18 µL TE with 1 ng/µL pBlueScript, five times; 2 µL of each
  dilution is the PCR template, giving final factors 10^3–10^7. The minimum input per
  library was 3.4 ng.
- 🟢 PCR (SI lines 18–28): custom HPLC-purified IDT primers with one 3′ phosphorothioate
  bond, TS-PCR Oligo1 5′-AATGATACGGCGACCACCGAG*A and TS-PCR Oligo2
  5′-CAAGCAGAAGACGGCATACGA*G, 0.5 µM each, 1x Q5 Hot Start High-Fidelity master mix,
  50 µL; 98 °C 30 s, 18 or 20 × (98 °C 10 s, 69 °C 30 s, 72 °C 30 s), 72 °C 2 min;
  1.0x AMPure XP. These match the patent values read earlier. The model now cites the SI,
  drops the inferred marking from the primers and polymerase, and adds the 1.0x AMPure
  cleanup with its discarded supernatant.
- 🟢 Dataset S1 Table S1: 36 of the 44 libraries used a 10^5-fold final dilution and
  18 cycles; six kidney libraries (Bot01–Bot06) used 10^7-fold and 20 cycles (footnote
  "100-fold difference in dilution in replicate samples"); the two low-input PMS2-/-
  libraries (Bot15, Bot16; 129 and 34 ng) are listed at 10-fold and 20 cycles. All 44
  were sequenced PE 2 × 100 on HiSeq 2500 Rapid Run. The SI text gives the serial scheme
  as 10^3–10^7, so the 10-fold entries presumably mean template taken without serial
  dilution; the SI does not say so, and the model keeps the typical 10^5 value.
- 🟢 MiSeq test (SI lines 29–47, Fig. S1): ~5 M clusters per library, 2 × 50; three
  well-barcoded samples per lane; optimum 5–10 duplicates per molecule when scaled to
  ½ HiSeq 2500 rapid lane (~70 M clusters).
- 🟢 WGS branch (SI lines 51–61): 10 µL of the 20 µL ligation amplified with TruSeq PCR
  Primer Cocktail and Master Mix, 1.0x AMPure, PE 2 × 100 on HiSeq 2500 at >30x. Used for
  germline and clonal filtering; kept outside the BotSeqS schematic.

Still open: the 2016 TruSeq DNA PCR-Free guide revision (adapter index set and length,
size-selection placement) and whether an index read was acquired. The SI says the MiSeq
samples were "well-barcoded" but names no index set.

## Native graph migration, 2026-10-18

Read three published originals: article PMC5024639 (1246–1256, 1460–1463), SI Appendix Materials and Methods (4–61), and Dataset S1 Table S1 (Bot01–Bot44). The bundle oligo scraper was run; both modified TS-PCR primers were checked in the SI. The patent and later vendor guide were not used to supply new sequence facts in this pass.

The graph owns the native-end family key, modified primers, Q5 PCR, 1.0X AMPure cleanup and PE100 readout locally. It separates physical material from computational Watson/Crick consensus families. The post-ligation stock also has a source-supported optional matched-WGS branch: the SI uses 10 µL of the original 20 µL stock before dilution, with TruSeq PCR Primer Cocktail/Master Mix, 1.0X cleanup and PE100 WGS. This describes the source's 32 WGS samples; two other individuals used Sanger confirmation. The stock outside the sampled BotSeqS PCR aliquot is unsampled material, not an asserted physical discard. WGS stock is retained separately.

The canonical single-index library remains an explicitly inferred compatibility model: exact 2016 PCR-free adapter/index revision and native Y-end geometry are unassigned. Adapter ligation and both PCR declarations therefore retain actual incoming substrates with unavailable native geometry; a finished compatible duplex is not used as proof of the incoming Y-adapter geometry. The older post-end-repair size-selection placement remains a labelled inference from the 2017 vendor guide, not a verified 2016 kit instruction. Pilot MiSeq molecular barcodes are not assigned to the main native-end BotSeqS key.

Recounting every Table S1 row gives 44 libraries: 36 at 10^5-fold/18 cycles, six at 10^7-fold/20 cycles and two at 10-fold/20 cycles, all PE100. Earlier note counts of 43/35 were repository counting errors, not publication discrepancies. The exact meaning of the two 10-fold entries remains unassigned; no extra dilution recipe is inferred. This capped pass makes no absence claim.

Implementation review: the adapter-ligated stock, retained WGS aliquot and diluted BotSeqS aliquot now share an explicit unresolved precursor Scene. Both unavailable-PCR views preserve that state without drawing a complementary strand from the inferred completed library. The canonical compatibility construct is used only for post-PCR products. This is a model correction using the source scope above; no additional originals were read.
