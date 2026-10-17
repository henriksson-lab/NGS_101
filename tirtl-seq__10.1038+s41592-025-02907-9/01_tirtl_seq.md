# TIRTL-seq

**Defining paper:** Pogorelyy et al., *Nature Methods* (2026), doi:10.1038/s41592-025-02907-9. Read from the open PMC copy (PMC12791010, JATS `s41592-025-02907-9_PMC12791010.xml.txt`), Methods "TIRTL-seq library preparation" lines 214–226, plus Supplementary Data 1 (`41592_2025_2907_MOESM3_ESM.xls.txt`). The protocols.io "TIRTL-seq 96-well" V.2 (dx.doi.org/10.17504/protocols.io.kqdg3k6xpv25/v2) was consulted only for the manual 96-well variant.

## Molecular boundary

🟢 Cells are distributed over replicate wells and lysed during miniaturized RT. A first multiplex PCR uses V- and C-segment primers carrying partial Nextera handles; C-segment reverse primers also carry a plate barcode.

🟢 A second PCR adds full Illumina adapters and one well-specific UDI pair. Alpha/beta pairing is inferred from clonotypes' shared occurrence patterns across replicate wells, not from two chains joined in one molecule.

## Reaction order (2026-10-17 bounded pass)

🟢 Methods line 218: 0.4 µl lysis/RT mix per well under 3 µl Vapor-Lock — 0.1% Triton X-100, 0.5 mM dNTP, 2.5 ng/µl random hexamer, 2.5 µM oligo(dT) (both Thermo), Maxima H RT buffer, 5 mM DTT, RNasin, 4 U/µl Maxima H Minus RT — then 0.25 µl cells; 42 °C 5 min, 25 °C 10 min, 50 °C 60 min, 94 °C 5 min.

🟢 Line 220: PCR I in the same well with KAPA2G Fast Multiplex Mix, 88 V-segment forward primers "containing Illumina Nextera reverse adapter sequence" and TRAC/TRBC reverse primers with "6- or 8-nt barcodes" and the "Nextera forward adapter sequence"; 95 °C 3 min, 20 × (95 °C 15 s, 59 °C 30 s, 72 °C 1 min), 72 °C 5 min.

🟢 Line 222: PCR-I product diluted ~1:10; 1 µl into PCR II with Q5 Hot Start and 1 µl well-specific i5/i7 UDI primer mix; 98 °C 30 s, 15 × (98 °C 10 s, 58 °C 10 s, 72 °C 50 s), 72 °C 2 min. Wells are pooled by inverted centrifugation, purified with 1X AMPure XP and sequenced on NovaSeq at 300,000 reads per well. Read lengths and index cycles are not printed in the article; the protocols.io 96-well V.2 (step 50) prefers 2×150 nt and gives a ~340 nt final library (step 49) for the manual variant only.

🟡 Orientation, derived from Supplementary Data 1: the C-segment primers begin with the Nextera Read 1 arm `TCGTCGGCAGCGTC…` and the i5 primers end in that arm, while V-segment primers and i7 primers carry the Read 2 arm `GTCTCGTGGGCTCGG…`. The final molecule is therefore P5 – i5 – Read 1 arm – plate barcode – C primer site – antisense CDR3/V(D)J – V primer site – Read 2 arm – i7 – P7, so Read 1 begins at the plate barcode. An earlier model placed the plate barcode beside the P7 end and drew a single UDI segment; this was a repository modelling error, corrected in the declaration.

🟡 Plate-barcode window: in both the alpha and beta C-primer sheets the 12 ID rows differ only at positions 34–40 after the mosaic end. The Methods say 6- or 8-nt barcodes; the supplement does not mark the boundary, so the page marks the variable window, not a declared barcode length.

## Evidence boundary

🟢 The published Supplementary Data 1 (`41592_2025_2907_MOESM3_ESM.xls.txt`) gives the full ordered panel in four sheets: 88 V-segment PCR-I forward primers, 24 plate-barcoded C-segment PCR-I reverse primers (12 alpha and 12 beta), 384 well-specific i5 PCR-II forward primers and 384 matching i7 PCR-II reverse primers. The V primers share the printed `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG` partial Nextera arm; the C primers begin `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG` before a 6- or 8-base plate identifier and the C-specific arm. The i5 and i7 sheets print complete P5/P7-tailed UDI primer sequences. There are no modification annotations in these four sheets.

🟢 The main article methods (`s41592-025-02907-9_PMC12791010.xml.txt`, line 218) also specify commercial random hexamer and oligo(dT) RT primers. Their synthesized sequences or oligo(dT) length are not printed and are not members of the 880-oligo PCR panel.

🟡 The schematic currently depicts the primer, plate barcode and well UDI **roles**; it does not encode all 880 panel members as segments or derive primer-binding and read-cycle spans from a selected panel row. The panel is exact for the published implementation, while target choice and index combinations can vary. No bases were inferred from raw reads because Supplementary Data 1 supplies the sequences directly.

## One cell per well

🟢 The defining paper's Methods describe FACS sorting one live CD3-positive T cell into
each 384-well lysis/RT reaction (original JATS line 204), followed by TIRTL-seq.
🟡 For that experiment, chain co-occurrence within a well gives pairing directly;
the pooled-cell implementation infers pairing from clonotype overlap across
replicate wells. The one-cell experiment changes input and interpretation while
using the published library chemistry, so it shares this chemistry page.

🟡 In the published reviewer response (MOESM2 text lines 709–722), the authors
suggest adding gene-specific primers during PCR I to phenotype or genotype a
single sorted cell. They do not give a tested primer panel or construct for that
extension. It is a proposed extension, not part of the one-cell TCR assay modeled
here.
