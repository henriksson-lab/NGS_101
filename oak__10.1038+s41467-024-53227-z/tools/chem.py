"""OAK overloading-and-unpacking combinatorial indexing."""
from multibranch import illumina_branch
from protocol_extensions import analyte_split, overload_then_unpack
from chemdraw import Row
from multimodal_spatial import (finish_pcr_branch, overload_carriers,
                                recover_and_redistribute_carriers)
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
TITLE="OAK — overloading and unpacking droplet barcodes"
NOTES="01_oak.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-024-53227-z">Darmanis et al., <i>Nature Communications</i> (2024)</a>.'
SUMMARY="Fixed cells or nuclei are overloaded into Chromium droplets for in situ first-round barcoding, recovered intact and redistributed into aliquots for a second index; RNA, ATAC and antibody-tag branches can share the resulting cell identity."
CAVEAT="The source prints example RNA aliquot and library PCR primers, listed below. The proprietary Chromium bead oligo, complete index families and Multiome kit oligos remain unavailable; the graph therefore leaves their reaction geometry unresolved."
RNA_TSO_PCR="AAGCAGTGGTATCAACGCAGAGT"
RNA_ALIQUOT_I5="AATGATACGGCGACCACCGAGATCTACACAACGTGATACACTCTTTCCCTACACGACGCTCTTCCGATCT"
RNA_PARTIAL_P5="AATGATACGGCGACCACCGAGA"
RNA_LIBRARY_I7="CAAGCAGAAGACGGCATACGAGATCGCATGTTACGTGACTGGAGTTCAGACGTGT"
def oligos():
    return [f"<b>{name}</b> (published RNA route, 5′→3′): <code>{sequence}</code>"
            for name, sequence in (("TSO recognition PCR primer", RNA_TSO_PCR),
                                   ("example aliquot i5 primer", RNA_ALIQUOT_I5),
                                   ("partial P5 library primer", RNA_PARTIAL_P5),
                                   ("example library i7 primer", RNA_LIBRARY_I7))]
CELL_PARTS=(("droplet barcode",16),("aliquot index",8))
RNA,RP=illumina_branch("OAK RNA","3-prime cDNA",cell_parts=CELL_PARTS,umi=10)
ATAC,AP=illumina_branch("OAK ATAC","accessible genomic DNA",architecture="nextera",cell_parts=CELL_PARTS)
FINAL_LIBRARIES=(("RNA library",RNA,RP,"Droplet and aliquot indexes jointly identify a cell's transcripts.","Illumina primers on the OAK RNA branch."),("ATAC library",ATAC,AP,"The matching two-round identity accompanies Multiome ATAC fragments.","Nextera primers on the OAK ATAC branch."))
def sections():
    return [("Fix cells or nuclei so they survive droplet recovery",[Row(chunks=[("fixed intact cell/nucleus containing retained indexed molecules", None, False)])],"The cell itself becomes the carrier between indexing rounds."),("Overload Chromium droplets and perform first-round barcoding",overload_then_unpack(droplet_barcode=16,well_barcode=8)[:1],"Many cells use each droplet's bead barcode."),("Break droplets, mix cells and redistribute into aliquots",overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:],"Aliquot-specific PCR supplies the orthogonal secondary index."),("Separate configured modalities",analyte_split("RNA / antibody-tag library","ATAC library when Multiome is used"),"The demonstration supports RNA alone and joint RNA–ATAC; the cell identity is shared."),]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    sample=g.input("sample","cells or nuclei for OAK",MaterialKind.BIOLOGICAL,primary=True)
    fixative=g.input("fixative","carrier-preserving fixation condition",MaterialKind.CONDITION)
    fixed=g.product("fixed","fixed cells or nuclei",MaterialKind.MOLECULAR,state=sections()[0][1])
    g.react("fix","Fix cells or nuclei so they survive droplet recovery",StepKind.CROSS_LINKING,((sample,InputRole.TEMPLATE),(fixative,InputRole.CONDITION)),(fixed,))

    droplets,bead,overload_reviewed=overload_carriers(
        g,fixed,prefix="droplet",
        state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[:1],
        bead_label="proprietary Chromium barcode bead oligos",inferred=True,
        source=SOURCE+" (proprietary bead sequence not disclosed)")

    rt=g.input("reverse_transcriptase","in-droplet reverse transcriptase",MaterialKind.CATALYST)
    copied=g.product("copied","fixed carriers containing droplet-barcoded cDNA and chromatin",MaterialKind.MOLECULAR,state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[:1],inferred=True,source=SOURCE+" (proprietary bead oligo structure not disclosed)")
    g.react("reverse_transcribe","Transfer the droplet barcode during in-situ reverse transcription",StepKind.REVERSE_TRANSCRIPTION,((droplets,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))

    transposome=g.input("transposome","Multiome ATAC transposome carrying the droplet identity route",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (proprietary loaded-adapter structure not disclosed)")
    tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=g.product("tagged","fixed carriers containing droplet-barcoded cDNA and tagged chromatin",MaterialKind.MOLECULAR,state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[:1],inferred=True,source=SOURCE+" (proprietary Multiome adapter structure not disclosed)")
    tagment_waste=g.product("tagment_waste","non-library chromatin tagmentation products",MaterialKind.WASTE)
    g.react("tagment","Tag chromatin for the configured Multiome ATAC branch",StepKind.TAGMENTATION,((copied,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,tagment_waste))
    g.discard(tagment_waste,"does not retain a library-eligible heterologous ATAC fragment")

    aliquoted,recovery_reviewed=recover_and_redistribute_carriers(
        g,tagged,prefix="carrier",
        recovered_state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:],
        aliquoted_state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:])

    rna=g.product("rna_template","OAK RNA templates carrying the shared droplet identity",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE+" (terminal primer sequences not transcribed)")
    atac=g.product("atac_template","OAK ATAC templates carrying the shared droplet identity",MaterialKind.MOLECULAR,state=ATAC,inferred=True,source=SOURCE+" (terminal primer sequences not transcribed)")
    g.react("separate_modalities","Separate the configured RNA and ATAC molecular branches",StepKind.SPLIT,((aliquoted,InputRole.TEMPLATE),),(rna,atac),changes_molecule=False)

    reviewed=["sample","fixative","fix",*overload_reviewed,"reverse_transcriptase","reverse_transcribe","transposome","tn5","tagment",*recovery_reviewed,"separate_modalities"]
    finish_pcr_branch(g,"rna",rna,RNA,RP,reviewed,label="OAK RNA",inferred_primers=True,limitation="The proprietary Chromium and aliquot-index primer sequences are not disclosed on this page.")
    finish_pcr_branch(g,"atac",atac,ATAC,AP,reviewed,label="OAK ATAC",inferred_primers=True,limitation="The proprietary Multiome and aliquot-index primer sequences are not disclosed on this page.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
