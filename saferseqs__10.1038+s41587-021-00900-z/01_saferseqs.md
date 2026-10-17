# SaferSeqS — research notes

**Defining paper:** Cohen *et al.*, “Detection of low-frequency DNA variants by targeted sequencing of the Watson and Crick strands,” *Nature Biotechnology* (2021), DOI [10.1038/s41587-021-00900-z](https://doi.org/10.1038/s41587-021-00900-z).

## Sources read

- 🟢 Defining open author manuscript, Methods and Fig. 1 (PMCID PMC8627329).
- 🟢 Publisher workbook `41587_2021_900_MOESM3_ESM.xlsx`, especially Supplementary Tables 5 and 7. Table 7 is the authority for adapter, nested-PCR and custom sequencing-primer oligos; Table 5 supplies target-specific GSPs.

## Construct-changing steps

Sections read in the 2026-10-17 bounded pass: JATS Methods "Library preparation",
"Anchored hemi-nested PCR", "Sequencing" and "Mutation calling" (`PMC8627329_PMC8627329.xml.txt`
lines 123–141); MOESM1 Supplementary Note "Library construction" (pdf.txt lines 180–209)
and bioinformatic pipeline (lines 314–333); MOESM3 Supplementary Tables 5 and 7.

1. 🟢 DNA is pretreated with 3 U USER enzyme (15 min, 37 °C), then taken through the Accel-NGS 2S DNA Library Kit (Swift 21024) with custom adapters; end repair 1 and 2 dephosphorylate and repair the ends (SPRI 2.0× and 1.8×).
2. 🟢 Ligation 1: a custom 50 µM 3′ adapter (replaces kit Y2) containing a 14-nt UID is ligated to the 3′ fragment ends (SPRI 1.2×).
3. 🟢 Ligation 2: a custom 42 µM 5′ adapter (replaces kit B2) is attached by a concerted polymerase / cohesive-end ligase / uracil DNA glycosylase reaction that synthesizes the UID complement, degrades the blocking Oligo #2 and ligates the extended adapter to the 5′ fragment end (Supplementary Note; SPRI 1.05×).
4. 🟢 Library PCR with library amplification primers #1 and #2 (2 µM each, NEBNext Ultra II Q5; 5, 7 or 11 cycles, 98 °C 10 s / 65 °C 75 s; eight 50-µl aliquots with extra Q5/dNTP/MgCl2 for 11 cycles; 1.8× SPRI) creates redundant copies of both original strands.
5. 🟢 Separate Watson and Crick first-round anchored PCRs use the same GSP1 with different short anchor primers (19 cycles); 1% of each is the template for a second round with the GSP2 primer and an anchor/indexing primer (17 cycles). Products are pooled and purified with 1.8× SPRI.
6. 🟢 Paired-end sequencing with eight-base dual indexing on MiSeq, HiSeq 2500 or HiSeq 4000 (read lengths not stated). Custom Read 1, index and Read 2 primers are combined with the standard Illumina primers at 1 µM. Reads must carry the constant tag GCCGTCGTTTTAT immediately after the UID (Supplementary Note pipeline).

## Watson/Crick labelling discrepancy (published)

🟢 Supplementary Table 7 assigns the P7-side anchor (`GACTGGAGTTCAGACGT*G*T`) and the P7
anchor/indexing primer to **Watson**, and the P5-side anchor and P5 anchor/indexing primer
to **Crick**. The Methods first-round passage agrees (Watson = "P7 short anchor"), but the
Methods second-round passage gives Watson the "anchor P5 indexing primer" and Crick the
"anchor P7 indexing primer", and the analysis passage says the Watson UID is the first 14
bases of Read 1 and the Crick UID the first 14 bases of Read 2.

🟡 Computed from the oligos: 3′ N14 Adapter Oligo #1 contains the reverse complement of
library amplification primer #2 / the Watson anchor immediately 3′ of its UID, so a Table 7
Watson product carries the UID beside P7 (read by Read 2) and a Table 7 Crick product
carries the copied UID beside P5 (read by Read 1, followed by GCCGTCGTTTTAT). The model
uses the Table 7 names and draws both products; which label set the authors intended is
unresolved. A P5 anchor cannot extend a first-round product made with the P7-side anchor,
so the Methods second-round names are not self-consistent with the first round.

🔴 Supplementary Table 7 lists no custom Index 2 primer for the Watson-primer-set product,
where the standard reverse-complement i5 primer has no site (the Watson GSP2 handle
follows i5). How that i5 was read is not stated in the sections inspected.

## Recovered oligos

🟢 Supplementary Table 7 gives both strands of the 3′ adapter, including the 14-nt
random UID, `/5Phos/`, `/3SpC3/`, deoxyuridines and `/33dT/`; the 5′ adapter; both
library-amplification primers; the Watson and Crick nested-PCR shells; and the custom
Read 1, index and Read 2 sequencing primers. These are now transcribed verbatim on the
schematic page. The only variable bases are the stated eight-base indices and the
target-specific GSP portions listed separately in Supplementary Table 5. The page instantiates
amplicon TP53_342 (GSP1 `GGGCATCCTTGAGTTCC*A*A`, GSP2 `CTTGAGTTCCAAGGCCTC*A*T`).

## Interpretation

🟡 “Identical” barcode means that the UID installed at one fragment end is copied into both complementary strands before amplification. Watson and Crick families are separated by oriented target PCR and rejoined computationally through that shared molecular identity.
