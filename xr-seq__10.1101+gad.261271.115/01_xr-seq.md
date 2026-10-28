# XR-seq: sequence the excised repair oligomer

Hu et al. (2015), [Genome-wide analysis of human global and transcription-coupled excision repair of UV damage at single-nucleotide resolution](https://doi.org/10.1101/gad.261271.115).

## Capped sources

One reader,28October2026; published article PMC4421983 (Results and complete XR-seq/Oligonucleotides Methods); published Supplemental Material (workflow/length and context figures/captions); published Supplemental Table5 workbook (sample identity, PCR cycles, input dishes, alignment counts). Mandatory full-bundle oligo scraper run; every adapter, spacer, phosphate and PCR-polymer hit checked against the article. No preprint or later optimized protocol substituted. This is a capped first pass, not a completed source audit.

## Preparation and physical selection

1. UV-treated cells generate excised nominal30mers during nucleotide excision repair. The represented NHF1/CS-B route uses254nm,20J/m²; recovery time is experiment-specific. Dounce lysis and centrifugation separate chromatin from soluble material. **TFIIH-associated excision products** are immunoprecipitated with anti-p89/XPB and anti-p62, with RNaseA. This is not initial anti-lesion capture of whole genomic DNA.
2. Wash TFIIH immunoprecipitates, elute65°C in SDS, PCI extract and ethanol precipitate. Treat with RNaseA/T1; PCI/G50/ethanol recover DNA. RNA digestion products stay in the preparation until physical purification.
3. Anneal the two splinted adapters and ligate with T4DNA ligase overnight16°C, then PCI/ethanol purify. Each splint has five random positions for end pairing; they are **not UMIs or sample barcodes**. Source5′ and3′ adapters remain distinct; SpacerC3 blocks unwanted3′ ends.
4. For **NHF1 and CS-B**, anti-CPD capture precedes anti-(6-4)PP capture of its **flowthrough**. XP-C reverses this order; it is a separate preparation variant. The graph follows the first order and does not discard the flowthrough before the second capture. Detailed lesion-IP buffers are source-linked to Hu2013 and deferred here; no invented recipe fills that boundary.
5. Repair the purified selected oligos with the appropriate photolyase:20pmol Drosophila CPD photolyase or1.25pmol E.coli(6-4)PP photolyase,2h on ice under366nm light. The physical photoproduct is removed before PCR. The source reserves1% unrepaired or repaired DNA for QC; this is not additional library input.
6. PCI/ethanol purify, KAPA HiFi HotStart PCR with RP1/RPIn, then PCI/ethanol and10%nativePAGE. Excise130–155nt bands, elute in0.3MNaCl6h, ethanol precipitate, quantify and pool. HiSeq2000,1×50. Table5 PCR cycles range6–15 depending on sample; no fixed universal cycle count is assigned.

## Verbatim5′→3′ oligos

| Name | Printed sequence/modification |
|---|---|
| A5F | GTTCAGAGTTCTACAGTCCGACGATC |
| A5R | NNNNNGATCGTCGGACTGTAGAACTCTGAAC-SpC3 |
| A3F | phos-TGGAATTCTCGGGTGCCAAGG-SpC3 |
| A3R | CCTTGGCACCCGAGAATTCCANNNNN-SpC3 |
| RP1 | AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA |
| RPIn | CAAGCAGAAGACGGCATACGAGATXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA |

Six X positions specify a selected TruSeq small-RNA-compatible index. Exact sample assignments remain unassigned. SpacerC3 contributes no nucleotide; the splint strands are not covalently inserted into the product. The conditional PCR-reference product is derived from the supplied primer polymers and their inward binding sites, retaining the complete adapter-to-insert junction rather than assuming full primer overlap.

## Meaning and limits

XR-seq sequences the repair-excised strand. It reports repair activity, not all genomic lesions, and does not itself separate damage formation from repair rate. Most observed inserts are20–30nt;26nt is an average, not a required molecule length. The lesion's enrichment near six bases from the3′ end is a population observation, not a coordinate hard-coded into every molecule. No adapter-degenerate window is extracted as a UMI.

Exact run-primer acquisition/index mask, source-linked Hu2013 IP details and photolyase-production originals, site-specific modification/strand-preservation validation and an independent reader remain deferred. No missing-information claim is accepted. Source citation to Hu2013 supports TFIIH/excision-product isolation ancestry; source photolyase citations support enzyme production, not inherited NGS adapters. Later XR-seq protocols are follow-ups, not substitutes for this original.
