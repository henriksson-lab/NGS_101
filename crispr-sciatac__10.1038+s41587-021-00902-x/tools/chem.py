"""CRISPR-sciATAC guide and accessibility readouts."""
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import indexed_molecule
TITLE="CRISPR-sciATAC — pooled perturbations with accessibility readout"
NOTES="01_crispr-sciatac.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41587-021-00902-x">Liscovitch-Brauer et al., <i>Nature Biotechnology</i> (2021)</a>.'
SUMMARY="A lentiviral guide cassette and accessible chromatin are recovered from the same combinatorially indexed nucleus, linking each perturbation to its ATAC profile."
CELL_PARTS=(("transposition index",8),("PCR-well index",12))
ATAC,AP=illumina_branch("CRISPR-sciATAC chromatin","accessible genomic DNA",architecture="nextera",cell_parts=CELL_PARTS)
GUIDE,GP=illumina_branch("CRISPR-sciATAC guide","integrated guide cassette",cell_parts=CELL_PARTS,guide=20)
FINAL_LIBRARIES=(("Accessibility library",ATAC,AP,"Indexed ATAC fragments report open chromatin.","Nextera primers on the accessibility branch."),("Guide-identity library",GUIDE,GP,"Targeted amplification recovers the perturbation guide with the same cell-index combination.","Illumina primers on the guide branch."))
def sections():
    return [("Introduce the lentiviral guide cassette",[Row(chunks=[("provirus — U6 promoter — guide spacer — scaffold", "cbc", False)])],"The integrated cassette is the recoverable perturbation identity."),("Index accessible chromatin in successive pools",indexed_molecule(payload="tagmented accessible DNA",barcode_parts=(("transposition index",8),("PCR-well index",12),)).rows(),"The barcode combination identifies a nucleus."),("Recover guide identity with the same index combination",indexed_molecule(payload="guide-cassette amplicon",barcode_parts=(("transposition index",8),("PCR-well index",12),),target_barcode=20).rows(),"Targeted semi-suppressive PCR links the integrated guide to the cell identity."),]
