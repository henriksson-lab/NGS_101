"""scRICA-seq circular cDNA isoform reconstruction plus ATAC."""
import illumina as il
import nextera as nx
from multibranch import illumina_branch
from multimodal_spatial import finish_pcr_branch
from protocol_extensions import analyte_split, indexed_molecule
from chemdraw import Construct, Row, Segment, feature, oligo, revcomp
from circular import CircularizationGeometry, ClosedCircle, rolling_circle
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
TITLE="scRICA-seq — RNA isoforms and chromatin accessibility"
NOTES="01_scrica-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-024-52335-0">Jiang et al., <i>Nature Communications</i> (2024)</a>.'
SUMMARY="Full-length cDNA receives end UMIs, is circularized and amplified so short Tn5 fragments can be grouped back into isoforms; a matched ATAC branch reports chromatin accessibility from the same cell."
CAVEAT="Supplementary Data 2 prints the complete low-throughput end-label, overlap-assembly, RCA, TSS/TES enrichment and terminal index-primer oligos. The page models that sequence-explicit branch; the separate high-throughput variant starts from proprietary 10x Multiome cDNA and is not substituted for it."
RNA,RP=illumina_branch("scRICA-seq isoform","RCA-derived full-length-cDNA fragment",architecture="nextera",cell=16,umi=12)
ATAC,AP=illumina_branch("scRICA-seq ATAC","accessible genomic DNA",architecture="nextera",cell=16)

DT_HANDLE="AAGCAGTGGTATCAACGCAGAGT"
TSO="GCAATGAAGTCGCAGGGTTG"
TAG1="AGTTAACGCTATACGGCACGCG"
LINKER=revcomp(DT_HANDLE)+TSO
TSO_ISPCR=TSO
OLIGO_PRIMER=DT_HANDLE
PHI29_PRIMER=TSO
TN5_I5=nx.S5+nx.ME[:6]
TN5_I7=nx.S7+nx.ME[:5]
ANTI_TAG1=revcomp(TAG1)
ANTI_TSO=revcomp(TSO)
I7_ANTI_TAG1=nx.ADAPTOR_S7+ANTI_TAG1
I7_ANTI_TSO=nx.ADAPTOR_S7+ANTI_TSO
INDEX_N7=il.P7+"N"*8+nx.ADAPTOR_S7
INDEX_N5=il.P5+"N"*8+nx.ADAPTOR_S5
POLYT="T"*30+"VN"
CELL=feature("cell_barcode","cell_barcode","random")
UMI_F=feature("umi","umi","random")
CIRCLE_TEMPLATE=Construct([
    Segment("Phi29 primer site",revcomp(PHI29_PRIMER)),
    Segment("full-length cDNA","X"*30,placeholder=True),
    Segment("oligo(dT) handle",DT_HANDLE),
    Segment("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL),
    Segment("UMI","U"*12,"umi",placeholder=True,feature=UMI_F),
    Segment("tag1",TAG1),
],name="one sequence-explicit strand of end-labelled full-length cDNA")
CIRCLE=ClosedCircle(CIRCLE_TEMPLATE,"tag1","Phi29 primer site")
RCA=rolling_circle(CIRCLE,PHI29_PRIMER,copies=3)
N5_PCR=PCRPrimer("Index 2 primer (N5 pool)",INDEX_N5,anneal_nt=len(nx.ADAPTOR_S5))
N7_PCR=PCRPrimer("Index 1 primer (N7 pool)",INDEX_N7,anneal_nt=len(nx.ADAPTOR_S7))
RNA_GEOMETRY=PCRGeometry.exact(RNA,N5_PCR,N7_PCR,statement="Supplementary Data 2 prints both complete N8 index-primer pools; their 3-prime Nextera arms land exactly on the endpoint library.")
ATAC_GEOMETRY=PCRGeometry.exact(ATAC,N5_PCR,N7_PCR,statement="The matched ATAC endpoint uses the same source-printed Nextera index-primer geometry.")
FINAL_LIBRARIES=(("Isoform library",RNA,RP,"Tagmented RCA copies share the UMI of one original full-length cDNA.","Nextera primers on the isoform branch."),("Accessibility library",ATAC,AP,"ATAC fragments carry the matching cell barcode.","Nextera primers on the accessibility branch."))

def oligos():
    yield oligo("Oligo(dT) primer",[Segment("SMART handle",DT_HANDLE),Segment("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL),Segment("UMI","U"*12,"umi",placeholder=True,feature=UMI_F),Segment("tag1",TAG1),Segment("oligo(dT)VN",POLYT,placeholder=True)])
    yield oligo("TSO primer",[Segment("TSO",TSO)],mods="/rG//rG//rG/ at 3-prime end")
    for name,seq in (("TSO ISPCR primer",TSO_ISPCR),("linker primer",LINKER),("oligo primer",OLIGO_PRIMER),("Phi29 primer",PHI29_PRIMER),("Tn5-i5 primer",TN5_I5),("Tn5-i7 primer",TN5_I7),("anti-tag1 primer",ANTI_TAG1),("anti-TSO primer",ANTI_TSO),("i7-anti-tag1 primer",I7_ANTI_TAG1),("i7-anti-TSO primer",I7_ANTI_TSO),("Index 1 primer (N7)",INDEX_N7),("Index 2 primer (N5)",INDEX_N5)):
        yield oligo(name,[Segment(name,seq,placeholder="N" in seq)])
def sections():
    return [("Capture RNA and accessible DNA from the same cell",analyte_split("full-length RNA/cDNA branch","ATAC branch"),"Low-throughput separation or the high-throughput microfluidic implementation preserves cellular pairing."),("Label both cDNA ends with molecule identity",indexed_molecule(payload="full-length cDNA",barcode_parts=(("cell barcode",16),),umi=12).rows(),"The printed oligo(dT) primer carries a 16-base cell barcode, 12-base UMI and tag1 after its SMART handle."),("Circularize full-length cDNA",CIRCLE.rows("NEBuilder overlap assembly"),"The linker is exactly reverse-complementary to the oligo(dT) handle followed by the TSO, establishing the end overlap used for circle assembly."),("Rolling-circle amplify the closed cDNA",RCA.rows("full-length-cDNA concatemer"),"The printed Phi29 primer binds the complementary TSO-derived site; the concatemer is computed from the circle."),("Tagment the amplified circle and build both libraries",[Row(chunks=[("RCA concatemer → random Tn5 fragments → TSS and TES libraries", "me", False)]),Row(chunks=[("accessible genomic DNA → ATAC library", "me", False)])],"Printed end-specific primers enrich TSS or TES fragments before the source-printed N5/N7 terminal PCR."),]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cell=g.input("cell","one cell containing RNA and chromatin",MaterialKind.BIOLOGICAL,primary=True)
    rna=g.product("rna","full-length RNA branch",MaterialKind.MOLECULAR,state=analyte_split("full-length RNA/cDNA branch","ATAC branch")[:1])
    chromatin=g.product("chromatin","matched accessible-chromatin branch",MaterialKind.MOLECULAR,state=analyte_split("full-length RNA/cDNA branch","ATAC branch")[1:])
    g.react("separate","Capture RNA and accessible DNA from the same cell",StepKind.SPLIT,((cell,InputRole.TEMPLATE),),(rna,chromatin),changes_molecule=False)

    rt_primer=g.input("rt_primer","printed barcoded UMI oligo(dT) reverse-transcription primer",MaterialKind.OLIGO)
    rt=g.input("reverse_transcriptase","full-length reverse transcriptase",MaterialKind.CATALYST)
    cdna=g.product("cdna","full-length cDNA",MaterialKind.MOLECULAR)
    g.react("reverse_transcribe","Reverse-transcribe the full-length RNA branch",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
    end_tags=g.input("end_tags","printed TSO, oligo(dT), ISPCR and linker oligos",MaterialKind.ADAPTER)
    labelled=g.product("labelled","full-length cDNA labelled at both ends",MaterialKind.MOLECULAR,state=indexed_molecule(payload="full-length cDNA",barcode_parts=(("cell barcode",16),),umi=12))
    g.react("label_ends","Label both cDNA ends with cell and molecule identity",StepKind.ADAPTER_ATTACHMENT,((cdna,InputRole.TEMPLATE),(end_tags,InputRole.STRUCTURAL)),(labelled,))
    circle_ligase=g.input("circle_ligase","cDNA circularization ligase",MaterialKind.CATALYST)
    circle=g.product("circle","end-labelled circular full-length cDNA",MaterialKind.MOLECULAR,state=CIRCLE)
    g.react("circularize","Circularize the end-labelled full-length cDNA",StepKind.CIRCULARIZATION,((labelled,InputRole.TEMPLATE),(circle_ligase,InputRole.CATALYST)),(circle,),geometry=CircularizationGeometry.known("overlap assembly",statement="The printed linker is reverse-complementary to the oligo(dT) handle followed by the TSO; the Methods specify NEBuilder HiFi assembly."))
    rca_polymerase=g.input("rca_polymerase","rolling-circle amplification polymerase",MaterialKind.CATALYST)
    rca_primer=g.input("rca_primer","printed Phi29 primer",MaterialKind.OLIGO)
    concatemer=g.product("concatemer","rolling-circle copies of full-length cDNA",MaterialKind.MOLECULAR,state=RCA)
    g.react("rolling_circle","Amplify the circular cDNA by rolling-circle extension",StepKind.EXTENSION,((circle,InputRole.TEMPLATE),(rca_primer,InputRole.STRUCTURAL),(rca_polymerase,InputRole.CATALYST)),(concatemer,))

    rna_transposome=g.input("rna_transposome","source-printed isoform-branch Tn5 adapter primers",MaterialKind.ADAPTER)
    rna_tn5=g.input("rna_tn5","isoform-branch Tn5 transposase",MaterialKind.CATALYST)
    rna_template=g.product("rna_template","tagmented and end-enriched rolling-circle cDNA copies",MaterialKind.MOLECULAR,state=RNA)
    rna_waste=g.product("rna_tagment_waste","non-library rolling-circle tagmentation products",MaterialKind.WASTE)
    g.react("tagment_rna","Tagment the rolling-circle cDNA copies",StepKind.TAGMENTATION,((concatemer,InputRole.TEMPLATE),(rna_transposome,InputRole.STRUCTURAL),(rna_tn5,InputRole.CATALYST)),(rna_template,rna_waste));g.discard(rna_waste,"does not retain a library-eligible heterologous isoform fragment")

    atac_transposome=g.input("atac_transposome","ATAC loaded Tn5 adapters",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (complete loaded-adapter oligos not transcribed)")
    atac_tn5=g.input("atac_tn5","ATAC Tn5 transposase",MaterialKind.CATALYST)
    atac_template=g.product("atac_template","tagmented accessible genomic DNA",MaterialKind.MOLECULAR,state=ATAC,inferred=True,source=SOURCE+" (terminal oligos not transcribed)")
    atac_waste=g.product("atac_tagment_waste","non-library chromatin tagmentation products",MaterialKind.WASTE)
    g.react("tagment_atac","Tagment accessible chromatin",StepKind.TAGMENTATION,((chromatin,InputRole.TEMPLATE),(atac_transposome,InputRole.STRUCTURAL),(atac_tn5,InputRole.CATALYST)),(atac_template,atac_waste));g.discard(atac_waste,"does not retain a library-eligible heterologous ATAC fragment")

    reviewed=["cell","separate","rt_primer","reverse_transcriptase","reverse_transcribe","end_tags","label_ends","circle_ligase","circularize","rca_primer","rca_polymerase","rolling_circle","rna_transposome","rna_tn5","tagment_rna","atac_transposome","atac_tn5","tagment_atac"]
    finish_pcr_branch(g,"rna",rna_template,RNA,RP,reviewed,label="scRICA-seq isoform",geometry=RNA_GEOMETRY)
    finish_pcr_branch(g,"atac",atac_template,ATAC,AP,reviewed,label="scRICA-seq accessibility",geometry=ATAC_GEOMETRY)
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
