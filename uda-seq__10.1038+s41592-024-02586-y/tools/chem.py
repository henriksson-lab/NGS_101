"""UDA-seq universal droplet-first combinatorial indexing."""
from multibranch import illumina_branch
from multimodal_spatial import (finish_pcr_branch, overload_carriers,
                                recover_and_redistribute_carriers)
from protocol_extensions import analyte_split, overload_then_unpack
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
TITLE="UDA-seq — droplet-first combinatorial indexing"
NOTES="01_uda-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41592-024-02586-y">Wang et al., <i>Nature Methods</i> (2025)</a>.'
SUMMARY="Fixed cells or nuclei are intentionally overloaded into droplets for the first barcode, recovered intact, redistributed to PCR wells and given a second index; the barcode pair resolves cells for several droplet-compatible modalities."
CAVEAT="UDA-seq wraps several assay chemistries. Supplementary Table 7 prints UDA indexing primers and the RNA-route Tn5 cassette, listed below. The complete proprietary droplet bead oligo and branch-specific kit primer geometry remain unresolved in this compact graph."
I5_FIRST="AATGATACGGCGACCACCGAGATCT-GGTTTACT-ACACTCTTTCCCTACACGACGCTCTTCCGATC-s-T"
I7_FIRST="CAAGCAGAAGACGGCATACGAGAT-CGCCGTTT-GTCTCGTGGGCTCGG"
TN5_TOP="[5Phos]CTGTCTCTTATACACATCT[3ddC]"
TN5_BOTTOM="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
def oligos():
    return [f"<b>{name}</b> (Supplementary Table 7, 5′→3′): <code>{sequence}</code>"
            for name, sequence in (("first i5-family primer", I5_FIRST),
                                   ("first i7-family primer", I7_FIRST),
                                   ("RNA-route Tn5-top_ME", TN5_TOP),
                                   ("RNA-route Tn5-bottom_Read2N", TN5_BOTTOM))]
CELL_PARTS=(("droplet barcode",16),("well barcode",8))
RNA,RP=illumina_branch("UDA-seq RNA","RNA-derived cDNA",cell_parts=CELL_PARTS,umi=10)
ATAC,AP=illumina_branch("UDA-seq ATAC","accessible genomic DNA",architecture="nextera",cell_parts=CELL_PARTS)
FINAL_LIBRARIES=(("Representative RNA library",RNA,RP,"Droplet barcode × well barcode identifies each recovered cell.","Illumina primers on an RNA UDA-seq branch."),("Representative ATAC library",ATAC,AP,"The same indexing logic can wrap a Multiome ATAC branch.","Nextera primers on an ATAC UDA-seq branch."))
def sections():
    return [("Fix cells or nuclei and prepare the chosen assay",overload_then_unpack(droplet_barcode=16,well_barcode=8),"The molecular target may be RNA, ATAC, VDJ or guide capture; UDA changes the indexing wrapper."),("Overload the microfluidic droplets",overload_then_unpack(droplet_barcode=16,well_barcode=8)[:1],"Multiple intact cells share one droplet barcode by design."),("Break the emulsion and recover intact cells",overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:],"Recovered cells are randomly redistributed among 96 or 384 wells for index PCR."),("Build modality-specific sublibraries",analyte_split("RNA / VDJ / guide branch as configured","DNA / ATAC branch as configured"),"The final construct depends on the wrapped assay; the two representative branches show the published Multiome use."),]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    sample=g.input("sample","cells or nuclei for UDA-seq",MaterialKind.BIOLOGICAL,primary=True)
    fixative=g.input("fixative","carrier-preserving fixation condition",MaterialKind.CONDITION)
    fixed=g.product("fixed","fixed cells or nuclei",MaterialKind.MOLECULAR,state=sections()[0][1])
    g.react("fix","Fix cells or nuclei for overloaded-droplet recovery",StepKind.CROSS_LINKING,((sample,InputRole.TEMPLATE),(fixative,InputRole.CONDITION)),(fixed,))

    droplets,bead,overload_reviewed=overload_carriers(
        g,fixed,prefix="droplet",
        state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[:1],
        bead_label="proprietary droplet-barcode bead oligos",inferred=True,
        source=SOURCE+" (proprietary bead sequence not disclosed)")
    rt=g.input("reverse_transcriptase","in-droplet reverse transcriptase",MaterialKind.CATALYST)
    copied=g.product("copied","fixed carriers containing droplet-barcoded cDNA and chromatin",MaterialKind.MOLECULAR,state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[:1],inferred=True,source=SOURCE+" (proprietary bead oligo structure not disclosed)")
    g.react("reverse_transcribe","Transfer the droplet barcode during RNA reverse transcription",StepKind.REVERSE_TRANSCRIPTION,((droplets,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))

    transposome=g.input("transposome","Multiome ATAC transposome",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (assay-dependent loaded-adapter structure not transcribed)")
    tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=g.product("tagged","fixed carriers containing barcoded cDNA and tagged chromatin",MaterialKind.MOLECULAR,state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[:1],inferred=True,source=SOURCE+" (assay-dependent loaded-adapter structure not transcribed)")
    tagment_waste=g.product("tagment_waste","non-library chromatin tagmentation products",MaterialKind.WASTE)
    g.react("tagment","Tag chromatin for the representative Multiome ATAC branch",StepKind.TAGMENTATION,((copied,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,tagment_waste))
    g.discard(tagment_waste,"does not retain a library-eligible heterologous ATAC fragment")

    aliquoted,recovery_reviewed=recover_and_redistribute_carriers(
        g,tagged,prefix="carrier",
        recovered_state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:],
        aliquoted_state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:])
    rna=g.product("rna_template","UDA-seq RNA templates carrying droplet identity",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE+" (terminal assay primers not transcribed)")
    atac=g.product("atac_template","UDA-seq ATAC templates carrying droplet identity",MaterialKind.MOLECULAR,state=ATAC,inferred=True,source=SOURCE+" (terminal assay primers not transcribed)")
    g.react("separate_modalities","Separate representative RNA and ATAC molecular branches",StepKind.SPLIT,((aliquoted,InputRole.TEMPLATE),),(rna,atac),changes_molecule=False)

    reviewed=["sample","fixative","fix",*overload_reviewed,"reverse_transcriptase","reverse_transcribe","transposome","tn5","tagment",*recovery_reviewed,"separate_modalities"]
    finish_pcr_branch(g,"rna",rna,RNA,RP,reviewed,label="UDA-seq RNA",inferred_primers=True,limitation="The assay-dependent RNA terminal primer pair is not disclosed in this compact module.")
    finish_pcr_branch(g,"atac",atac,ATAC,AP,reviewed,label="UDA-seq ATAC",inferred_primers=True,limitation="The assay-dependent ATAC terminal primer pair is not disclosed in this compact module.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
