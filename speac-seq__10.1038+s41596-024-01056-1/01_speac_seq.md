# SPEAC-seq

**Defining source:** doi:10.1038/s41596-024-01056-1.

## Molecular boundary

🟢 One CRISPR-perturbed cell and one fluorescent reporter cell are co-cultured per droplet; responsive droplets are sorted and their integrated guide cassettes PCR-amplified.

## Evidence boundary

🟢 This is a pooled functional screen. Sequencing measures guide enrichment by sorted fraction, not a cell-resolved transcriptome.

## Oligo source boundary

🟢 The original SPEAC-seq study, Wheeler et al., *Science* 2023
(doi:10.1126/science.abq4822; PMC10066924), confirms that genomic DNA
from sorted droplets is PCR-amplified to recover integrated sgRNAs for
Illumina sequencing. The Nature Protocols method describes the same
functional assay and links the original study as its key reference.

🔴 The exact guide-flanking genomic PCR primers, indexing primers and
sequencing-primer site cannot be transcribed from the publicly retrievable
full text. The Nature Protocols article is subscription-gated, and its
downloadable attachments are device CAD files rather than an oligo table.
The original *Science* supplementary PDF and data ZIP are exposed at
PMC10066924 but the source fetcher records a PMC download gate for both.
Those files must be reviewed before interpreting their absence as
unpublished. GSE200457 is a candidate raw-read source for a future
**inference** of invariant adapter sequence, but it cannot identify
the chemical identity of a synthesis oligo or prove primer names.
