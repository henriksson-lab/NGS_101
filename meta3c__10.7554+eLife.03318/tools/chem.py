"""meta3C restriction-ligation library from mixed microbial populations."""
from multibranch import illumina_branch
from protocol_extensions import proximity_product
from chemdraw import Row
TITLE="meta3C — metagenomic chromosome-conformation capture"
NOTES="01_meta3c.html"
SOURCE='Defining source: <a href="https://doi.org/10.7554/eLife.03318">Marbouty et al., <i>eLife</i> (2014)</a>.'
SUMMARY="Crosslinked cells from an intact microbial mixture are restriction-digested and diluted-ligated; ordinary shotgun library construction preserves both contact junctions and unligated genomic fragments without biotin enrichment."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("meta3C","microbial 3C fragment")
FINAL_CAPTION="Paired-end reads include ordinary genomic fragments and restriction-junction contact products."
SEQUENCING_INTRO="Declared Illumina primers bind the finished meta3C library."
def sections():
    return [("Crosslink the mixed microbial community",[Row(chunks=[("species A cells + species B cells + … → formaldehyde-fixed chromosomes", None, False)])],"Crosslinking occurs before lysis so contacts retain cellular provenance."),("Digest fixed chromatin",[Row(chunks=[("bacteria: HpaII C^CGG; yeast: DpnII G^ATC", "me", False)])],"The enzyme is chosen for the organisms and crosslinking conditions."),("Dilute and ligate proximal restriction ends",proximity_product(left="microbial locus A",right="microbial locus B",bridge="restriction junction").rows(),"Ligation produces within-cell contact chimeras; no biotin fill-in or junction pull-down is used."),("Reverse crosslinks, shear and make a paired-end library",[Row(chunks=[("3C pool → Covaris shear → custom PE adapters → size select → PCR", "r1", False)])],"Keeping non-junction fragments also supplies shotgun coverage for metagenome assembly."),]
