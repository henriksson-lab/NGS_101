# scID-seq — plate-barcoded antibody-DNA detection

van Eijl et al. (2018), *Single-Cell ID-seq Reveals Dynamic BMP Pathway Activation Upstream of the MAF/MAFB-Program in Epidermal Differentiation*, [10.1016/j.isci.2018.11.009](https://doi.org/10.1016/j.isci.2018.11.009).

## Sources inspected

Three published originals: defining article XML (Results scID-seq workflow),
`j.isci.2018.11.009_supp_mmc1.pdf` Transparent Methods, antibody conjugation,
staining/sorting and library preparation (135–185), and published
`j.isci.2018.11.009_supp_mmc2.xlsx` Table S1, antibody/barcode overview.
The full source-bundle scraper was run; its scID primer hits were checked against
Transparent Methods. The separate bulk CEL-seq2 and ChIP-seq validation assays
are not the single-cell protein library and their primers are not substituted.

## Source workflow

🟢 dsDNA/antibody functionalization uses NHS-S-S-PEG4-tetrazine and
DBCO-PEG12-TCO, followed by conjugation and quenching. Cells are PFA-fixed,
glycine-quenched, permeabilized/stained, washed and single-cell sorted into
96-well plates with DTT release buffer under Vapor-lock. The source uses a
15-base UMI and antibody identities from its barcode dictionary.

🟢 PCR1 preamplifies reporters for 15 cycles **before** cell-specific barcoding.
PCR2 adds a well-specific barcode in five cycles. Material is pooled per plate,
oil removed, QIAquick-purified, treated with Exonuclease I and purified again.
PCR3 uses a common long forward and plate-specific reverse primer for another
five cycles, followed by cleanup. The expected library is approximately 185 bp
and is sequenced on NextSeq500. Exact well/plate index sequences and cycle spans
remain unresolved without the cited predecessor material.

## Verbatim scID primers and identity fragments

🟢 Forward amplification: `CACGACGCTCTTCCGATCT`.
Reverse amplification: `TCGCTTATCTGTTGACTGAT`.
Fw_long: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`.
The model imports canonical Illumina regions where they match these source
sequences. Table S1 assigns Actin barcode 1 `GTATCGTCGT`; this fragment and the
15-base UMI are shown separately, not joined into an invented full reporter.

The local graph begins at source-prepared, sorted wells and represents two
symbolic well classes to make plate pooling explicit. It does not invent
experimental well IDs or barcode bases. All full PCR product geometries remain
unassigned. Plate/well identity semantics are distinguished from antibody
identity and the UMI; the source's bulk RNA reads are not a second scID output.

## Scope and lineage

The Transparent Methods explicitly inherits downstream barcode sequences from
**ID-seq (van Buggenum et al., 2018)** and conjugation from van Buggenum et al.
(2016). These are component-level ancestry candidates. Neither cited original
was read within this cap; no hidden predecessor constructor supplies chemistry.
Deferred: both originals and their exact oligo lists; complete dsDNA reporter,
well/plate primer designs and source-specific read layout; remaining original
sections and independent review. Status remains notes; no missing component is
confirmed absent.


## Public preparation view, 2026-10-30

The public representation is explicitly **illustrative_example**: two invented symbolic well classes illustrate the source per-well procedure and all-material per-plate pooling. They are not two experimental cell identities or a claim of96occupied wells. The source96wellplate capacity is stated separately. Input cells are already fixed, immunostained and sorted; manufacture/staining of the antibody reporters remains an explicit preparation boundary. Source Actinbarcode1 and UMI15 are separate identity fragments, not an assembled full reporter. DTT is a molecular contribution in the actual1uL releasebuffer under7uL Vaporlock. PCR1adds4uL enzyme/dNTP/buffer/primer mix; PCR2adds1uL buffer/forward/wellreverse mix and uses enzyme/nucleotide carryover. These component ports do not mean fresh independent additions at PCR2. PCR3after two purifications receives a fresh supply interface with enzyme identity/formulation unassigned in the inspected20uL recipe, rather than presumed PCR1enzyme carryover. Final agarose/BioAnalyzer material, retained library, loading aliquot and untransferred stock remain distinct.

NextSeq500 alone does not assign a readlength, acquired indexmask, sequencing-primer landing or complete well/plate barcode position. Full ID-seq predecessor strands and exact well/plate primer sequences remain unassembled; the known local amplification primers remain source-exact. No bulk RNA/ChIP validation library is used as an extra scID output.



Three selected published originals: mainPMC6249387 ResultsUMI/workflow; TransparentMethodsmmc1 text135–185 with actualPDFp12; TableS1mmc2 originalZIP/XML headers/Actinrow3. No predecessor original or bulkRNA/ChIP recipe adopted.

Source investigation closed before implementation; mandatory scrapers ran during the
read-only preflight. Literal curator source inventories include the stated technical
slot/pair counts, not graph-enumerated automatic coverage. Normal gates, source
transcription, physical-fate checks and full workflow rendering pass with unresolved
chemistry visible. No sequence-exact final product, assay performance, full independent
source audit or accepted missing-information claim is asserted.
