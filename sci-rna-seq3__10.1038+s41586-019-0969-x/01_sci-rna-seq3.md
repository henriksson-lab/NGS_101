# sci-RNA-seq3 — research notes

sci-RNA-seq3 is documented from Cao *et al.*, 2019,
doi:[10.1038/s41586-019-0969-x](https://doi.org/10.1038/s41586-019-0969-x).

The [2017 predecessor research note](../sci-rna-seq__10.1126+science.aam8940/01_sci-rna-seq.md)
describes the earlier two-level implementation. The 2019 source checks below and
the public schematic use the sci-RNA-seq3 paper and Supplementary Table S11 as
their molecular sources; precursor and successor chemistry have separate evidence
boundaries.

## Linked author sources checked on 2026-10-16

The 2019 paper's reporting summary links the
[author pipeline repository](https://github.com/JunyueC/sci-RNA-seq3_pipeline).
Its `sci3_primer_sequences_plate.xls` gives plate positions and sequences for 384 RT
primers, 384 ligation hairpins, 96 P5 primers and 96 P7 primers. The representative
members used in the model occur in both this workbook and Supplementary Table S11. Its
`RT_barcode_embryo_information.csv` has 382 distinct RT-group barcodes assigned to 61
embryos; all 382 match the 384 RT-barcode sequences in the primer workbook. Thus the
embryo-to-RT-barcode mapping is published in linked author material. The CSV does not
give each nucleus's subsequent ligation well or P5/P7 pair.

The reporting summary also links [GEO GSE119945](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119945).
Its cell annotation `sample` field records an observed PCR-group label followed by
the ligation and RT barcode sequences for each called cell. For example,
`sci3-me-001.GTCGGAGTTTGAGGTAGAA` contains a nine-base ligation read-side barcode
and the ten-base RT barcode `TGAGGTAGAA`, both present in the author workbook. The GEO
sample metadata says cell IDs attach RT and ligation barcodes; the author pipeline
calls `sci3-me-001` a PCR group. The remaining question is which ordered P5/P7 primer
pair generated each numeric PCR-group label. These observed barcode combinations
must not be described as unpublished. A second independent read of GEO is pending.

The [later author protocol](https://doi.org/10.17504/protocols.io.9yih7ue) prints the
RT, ligation, P5 and P7 architectures and the operational workflow. Its sequencing
example uses Read 1 34, Read 2 100 and two 10-cycle index reads; the defining 2019
article specifies Read 1 34, Read 2 52 and two 10-cycle index reads. The schematic
uses the defining article's 2019 run. This difference may reflect a later run setup;
the available sources do not establish that either printed recipe is erroneous.

There is a separate source discrepancy: the 2019 Methods print `/dideoxyU/` for the
internal hairpin modification, whereas Supplementary Table S11 and the author plate
workbook print `/ideoxyU/`. The model follows the two explicit oligo lists and the
subsequent USER cleavage step. The intended chemical specification behind the Methods
wording remains uncertain; the public schematic identifies both source values.

## Bounded Methods pass, 2026-10-17

The 2019 Methods section "sci-RNA-seq3 library preparation and sequencing"
(PMC6434952 JATS lines 176–190) was reread against the graph. The model now names
SuperScript IV RT (4–55 °C gradient), T4 DNA ligase with betaine and PEG 8000
(16 °C, 3 h), the NEB mRNA Second Strand Synthesis buffer and enzyme (16 °C,
180 min), i7-only TDE1 (55 °C, 5 min), the Zymo binding-buffer stop and 1.5× AMPure
cleanup, on-bead USER elution (37 °C, 15 min), NEBNext High-Fidelity 2× PCR Master
Mix (72 °C 5 min, then 12–14 cycles) and the pooled 0.8× AMPure cleanup. The
earlier "random second-strand primer" input was a model assumption: the Methods
name only the NEB module, so no primer is declared. Fragments lacking the hairpin
end and the USER-released hairpin arm are not physically removed (only the beads
are removed after USER), so their fates are unresolved rather than discarded.

The Methods paragraph opening "After ligation reaction" (line 180) describes the
ligation itself, and the one opening "After RT reaction" (line 182) describes the
step after ligation. The openings appear transposed; the model keeps RT →
ligation → second strand, and the schematic shows both locations.

The source checklist and two-reader status are in
`catalogue/source_reviews/sci-rna-seq3__10.1038+s41586-019-0969-x.json`. The cited
2017 predecessor's supplement is still inaccessible from the current source routes,
so inherited sequence claims and the PCR-group-to-P5/P7 pair mapping are blocked
from being called confirmed absent.
