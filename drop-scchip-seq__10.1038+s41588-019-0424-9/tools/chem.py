"""Drop-scChIP-seq chromatin fragmentation, immunocapture and droplet indexing."""
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import indexed_molecule
TITLE="Drop-scChIP-seq — droplet single-cell ChIP-seq"
NOTES="01_drop-scchip-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41588-019-0424-9">Rotem et al., <i>Nature Genetics</i> (2019)</a>.'
SUMMARY="MNase-released, antibody-selected chromatin fragments are compartmentalized with barcode beads so each retained fragment acquires a cell identity before Illumina amplification."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("Drop-scChIP-seq","immunoprecipitated chromatin",cell=12,umi=8)
FINAL_CAPTION="Cell barcode and UMI identify immunoprecipitated chromatin fragments."
SEQUENCING_INTRO="Declared Illumina primers bind the finished Drop-scChIP-seq library."
def sections():
    return [("Digest chromatin with MNase",[Row(chunks=[("nucleosome-protected DNA — MNase-cut ends", None, False)])],"MNase releases chromatin fragments while preserving histone-bound material."),("Immunoprecipitate the chosen chromatin mark",[Row(chunks=[("fragment — marked nucleosome   ↑ antibody capture", "w1", False)])],"The antibody defines the assayed target; it is not encoded as a sequence barcode."),("Transfer droplet barcode and UMI",indexed_molecule(payload="ChIP fragment",barcode_parts=(("cell barcode",12),),umi=8).rows(),"Droplet co-encapsulation assigns cellular and molecular identity before PCR."),]
