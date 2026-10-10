"""Spear-ATAC accessibility and guide-capture branches."""
from chemdraw import Row
from multibranch import illumina_branch
TITLE="Spear-ATAC — single-cell ATAC CRISPR screening"
NOTES="01_spear-atac.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-021-23213-w">Pierce et al., <i>Nature Communications</i> (2021)</a>.'
SUMMARY="A 10x scATAC-derived workflow captures accessible chromatin and enriches guide identity from the same barcoded nuclei."
CAVEAT="The paper defines the shared cell barcode and the two enrichment branches, but does not print each completed library as one contiguous molecule. Outer sequencing shells are marked inferred."
ATAC,AP=illumina_branch("Spear-ATAC chromatin library","accessible chromatin fragment",architecture="nextera",cell=16,inferred_adapters=True)
GUIDE,GP=illumina_branch("Spear-ATAC guide library","guide-capture amplicon",cell=16,guide=20,inferred_adapters=True)
FINAL_LIBRARIES=(("Chromatin-accessibility library",ATAC,AP,"Tagmented chromatin fragments carry the droplet cell barcode.","Nextera-compatible primers read the ATAC branch."),("Guide-identity library",GUIDE,GP,"Guide enrichment preserves the same droplet cell barcode.","Illumina primers read the guide branch."))
def sections(): return [
 ("Tagment nuclei before partitioning",[Row(chunks=[("accessible chromatin ** Tn5 mosaic ends","me",False)])],"Bulk tagmentation marks open chromatin."),
 ("Barcode nuclei in droplets",[Row(chunks=[("ATAC fragment + gel-bead cell barcode","cbc",False)])],"Gap fill transfers one droplet identity to accessible fragments."),
 ("Amplify the shared barcoded material",[Row(chunks=[("barcoded ATAC products → PCR pool","r1",False)])],"The pool contains genomic ATAC molecules and guide-associated material."),
 ("INFERRED — complete the two sequencing shells",[Row(chunks=[("ATAC library | guide-specific linear PCR → guide library","r2",False)])],"Guide enrichment is source-established; the fully drawn outer adapter completion is inferred."),]
