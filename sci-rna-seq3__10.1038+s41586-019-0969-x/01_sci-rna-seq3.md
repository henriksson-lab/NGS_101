# sci-RNA-seq3 — research-note pointer

sci-RNA-seq3 is documented from Cao *et al.*, 2019,
doi:[10.1038/s41586-019-0969-x](https://doi.org/10.1038/s41586-019-0969-x).

The detailed source record was originally written jointly for both versions and remains
available as provenance in the
[combined sci-RNA-seq / sci-RNA-seq3 research note](../sci-rna-seq__10.1126+science.aam8940/01_sci-rna-seq.md#sci-rna-seq3-nature-methods-unless-marked).
The public schematic uses the sci-RNA-seq3 paper and Supplementary Table S11 as its
authoritative molecular sources.

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

The source checklist and two-reader status are in
`catalogue/source_reviews/sci-rna-seq3__10.1038+s41586-019-0969-x.json`. The cited
2017 predecessor's supplement is still inaccessible from the current source routes,
so inherited sequence claims and the PCR-group-to-P5/P7 pair mapping are blocked
from being called confirmed absent.
