"""scRICA-seq circular cDNA isoform reconstruction plus ATAC."""
from multibranch import illumina_branch
from multimodal_spatial import finish_pcr_branch
from protocol_extensions import analyte_split, indexed_molecule
from chemdraw import Row
from circular import CircularizationGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
TITLE="scRICA-seq — RNA isoforms and chromatin accessibility"
NOTES="01_scrica-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-024-52335-0">Jiang et al., <i>Nature Communications</i> (2024)</a>.'
SUMMARY="Full-length cDNA receives end UMIs, is circularized and amplified so short Tn5 fragments can be grouped back into isoforms; a matched ATAC branch reports chromatin accessibility from the same cell."
CAVEAT="The source establishes end labelling, circularization, rolling-circle amplification and matched ATAC processing. This compact module has not transcribed the complete end-labelling and terminal index-primer oligos, so those structural inputs remain explicitly inferred."
RNA,RP=illumina_branch("scRICA-seq isoform","RCA-derived full-length-cDNA fragment",architecture="nextera",cell=16,umi=12)
ATAC,AP=illumina_branch("scRICA-seq ATAC","accessible genomic DNA",architecture="nextera",cell=16)
FINAL_LIBRARIES=(("Isoform library",RNA,RP,"Tagmented RCA copies share the UMI of one original full-length cDNA.","Nextera primers on the isoform branch."),("Accessibility library",ATAC,AP,"ATAC fragments carry the matching cell barcode.","Nextera primers on the accessibility branch."))
def sections():
    return [("Capture RNA and accessible DNA from the same cell",analyte_split("full-length RNA/cDNA branch","ATAC branch"),"Low-throughput separation or the high-throughput microfluidic implementation preserves cellular pairing."),("Label both cDNA ends with molecule identity",indexed_molecule(payload="full-length cDNA",barcode_parts=(("cell barcode",16),),umi=12).rows(),"End labels let fragments from one original transcript be regrouped."),("Circularize and amplify full-length cDNA",[Row(chunks=[("end-labelled cDNA ** end-to-end ligation → circle → rolling-circle copies", "w1", False)])],"Circular amplification makes repeated copies of the same full-length molecule."),("Tagment the amplified circle and build both libraries",[Row(chunks=[("RCA concatemer → random Tn5 fragments → isoform library", "me", False)]),Row(chunks=[("accessible genomic DNA → ATAC library", "me", False)])],"Short isoform fragments sharing a UMI are integrated computationally across TSS, TES and exons."),]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cell=g.input("cell","one cell containing RNA and chromatin",MaterialKind.BIOLOGICAL,primary=True)
    rna=g.product("rna","full-length RNA branch",MaterialKind.MOLECULAR,state=analyte_split("full-length RNA/cDNA branch","ATAC branch")[:1])
    chromatin=g.product("chromatin","matched accessible-chromatin branch",MaterialKind.MOLECULAR,state=analyte_split("full-length RNA/cDNA branch","ATAC branch")[1:])
    g.react("separate","Capture RNA and accessible DNA from the same cell",StepKind.SPLIT,((cell,InputRole.TEMPLATE),),(rna,chromatin),changes_molecule=False)

    rt_primer=g.input("rt_primer","full-length RNA reverse-transcription primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete oligo not transcribed in compact module)")
    rt=g.input("reverse_transcriptase","full-length reverse transcriptase",MaterialKind.CATALYST)
    cdna=g.product("cdna","full-length cDNA",MaterialKind.MOLECULAR,inferred=True,source=SOURCE+" (complete oligo not transcribed in compact module)")
    g.react("reverse_transcribe","Reverse-transcribe the full-length RNA branch",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
    end_tags=g.input("end_tags","cell-barcode and molecule-UMI end-tag oligos",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (complete end-tag oligos not transcribed in compact module)")
    labelled=g.product("labelled","full-length cDNA labelled at both ends",MaterialKind.MOLECULAR,state=indexed_molecule(payload="full-length cDNA",barcode_parts=(("cell barcode",16),),umi=12),inferred=True,source=SOURCE+" (complete end-tag oligos not transcribed in compact module)")
    g.react("label_ends","Label both cDNA ends with cell and molecule identity",StepKind.ADAPTER_ATTACHMENT,((cdna,InputRole.TEMPLATE),(end_tags,InputRole.STRUCTURAL)),(labelled,))
    circle_ligase=g.input("circle_ligase","cDNA circularization ligase",MaterialKind.CATALYST)
    circle=g.product("circle","end-labelled circular full-length cDNA",MaterialKind.MOLECULAR,state=sections()[2][1])
    g.react("circularize","Circularize the end-labelled full-length cDNA",StepKind.CIRCULARIZATION,((labelled,InputRole.TEMPLATE),(circle_ligase,InputRole.CATALYST)),(circle,),geometry=CircularizationGeometry.unavailable("The compact source model establishes end-to-end cDNA circularization but has not transcribed the ligating end geometry."))
    rca_polymerase=g.input("rca_polymerase","rolling-circle amplification polymerase",MaterialKind.CATALYST)
    concatemer=g.product("concatemer","rolling-circle copies of full-length cDNA",MaterialKind.MOLECULAR,state=sections()[2][1])
    g.react("rolling_circle","Amplify the circular cDNA by rolling-circle extension",StepKind.EXTENSION,((circle,InputRole.TEMPLATE),(rca_polymerase,InputRole.CATALYST)),(concatemer,))

    rna_transposome=g.input("rna_transposome","isoform-branch loaded Tn5 adapters",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (complete loaded-adapter oligos not transcribed)")
    rna_tn5=g.input("rna_tn5","isoform-branch Tn5 transposase",MaterialKind.CATALYST)
    rna_template=g.product("rna_template","tagmented rolling-circle cDNA copies",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE+" (terminal oligos not transcribed)")
    rna_waste=g.product("rna_tagment_waste","non-library rolling-circle tagmentation products",MaterialKind.WASTE)
    g.react("tagment_rna","Tagment the rolling-circle cDNA copies",StepKind.TAGMENTATION,((concatemer,InputRole.TEMPLATE),(rna_transposome,InputRole.STRUCTURAL),(rna_tn5,InputRole.CATALYST)),(rna_template,rna_waste));g.discard(rna_waste,"does not retain a library-eligible heterologous isoform fragment")

    atac_transposome=g.input("atac_transposome","ATAC loaded Tn5 adapters",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (complete loaded-adapter oligos not transcribed)")
    atac_tn5=g.input("atac_tn5","ATAC Tn5 transposase",MaterialKind.CATALYST)
    atac_template=g.product("atac_template","tagmented accessible genomic DNA",MaterialKind.MOLECULAR,state=ATAC,inferred=True,source=SOURCE+" (terminal oligos not transcribed)")
    atac_waste=g.product("atac_tagment_waste","non-library chromatin tagmentation products",MaterialKind.WASTE)
    g.react("tagment_atac","Tagment accessible chromatin",StepKind.TAGMENTATION,((chromatin,InputRole.TEMPLATE),(atac_transposome,InputRole.STRUCTURAL),(atac_tn5,InputRole.CATALYST)),(atac_template,atac_waste));g.discard(atac_waste,"does not retain a library-eligible heterologous ATAC fragment")

    reviewed=["cell","separate","rt_primer","reverse_transcriptase","reverse_transcribe","end_tags","label_ends","circle_ligase","circularize","rca_polymerase","rolling_circle","rna_transposome","rna_tn5","tagment_rna","atac_transposome","atac_tn5","tagment_atac"]
    finish_pcr_branch(g,"rna",rna_template,RNA,RP,reviewed,label="scRICA-seq isoform",inferred_primers=True,limitation="The complete isoform-branch index-primer pair has not been transcribed in this compact module.")
    finish_pcr_branch(g,"atac",atac_template,ATAC,AP,reviewed,label="scRICA-seq accessibility",inferred_primers=True,limitation="The complete ATAC-branch index-primer pair has not been transcribed in this compact module.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
