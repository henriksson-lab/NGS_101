"""scRICA-seq circular cDNA isoform reconstruction plus ATAC."""
import illumina as il
import nextera as nx
from multibranch import illumina_branch
from multimodal_spatial import finish_pcr_branch
from protocol_extensions import analyte_split, indexed_molecule
from chemdraw import Construct, Row, Segment, feature, oligo, revcomp
from circular import CircularizationGeometry, ClosedCircle, rolling_circle
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (AnnotationSlot, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)
TITLE="scRICA-seq — RNA isoforms and chromatin accessibility"
NOTES="01_scrica-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-024-52335-0">Jiang et al., <i>Nature Communications</i> (2024)</a>.'
SUMMARY="Full-length cDNA receives end UMIs, is circularized and amplified so short Tn5 fragments can be grouped back into isoforms; a matched ATAC branch reports chromatin accessibility from the same cell."
CAVEAT="Supplementary Data 2 prints the low-throughput N6 cell barcode, end-label, overlap-assembly, RCA, TSS/TES enrichment and terminal index-primer oligos. This page models that branch, with TSS/TES PCR1/PCR2 intermediate ends still unresolved in the molecular model. The separate high-throughput 10x Multiome variant uses a 16-base cell barcode and is not substituted for it."
TSS_LONG,TSS_PRIMERS=illumina_branch("scRICA-seq TSS long","RCA-derived long TSS fragment",architecture="nextera",cell=6,umi=12)
TSS_SHORT,_=illumina_branch("scRICA-seq TSS short","RCA-derived short TSS fragment",architecture="nextera",cell=6,umi=12)
TES_LONG,TES_PRIMERS=illumina_branch("scRICA-seq TES long","RCA-derived long TES fragment",architecture="nextera",cell=6,umi=12)
TES_SHORT,_=illumina_branch("scRICA-seq TES short","RCA-derived short TES fragment",architecture="nextera",cell=6,umi=12)
ATAC,AP=illumina_branch("scRICA-seq ATAC","accessible genomic DNA",architecture="nextera")

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
    Segment("cell barcode","B"*6,"cbc",placeholder=True,feature=CELL),
    Segment("UMI","U"*12,"umi",placeholder=True,feature=UMI_F),
    Segment("tag1",TAG1),
],name="one sequence-explicit strand of end-labelled full-length cDNA")
CIRCLE=ClosedCircle(CIRCLE_TEMPLATE,"tag1","Phi29 primer site")
RCA=rolling_circle(CIRCLE,PHI29_PRIMER,copies=3)
N5_PCR=PCRPrimer("Index 2 primer (N5 pool)",INDEX_N5,anneal_nt=len(nx.ADAPTOR_S5))
N7_PCR=PCRPrimer("Index 1 primer (N7 pool)",INDEX_N7,anneal_nt=len(nx.ADAPTOR_S7))
TSS_GEOMETRY=PCRGeometry.exact(TSS_LONG,N5_PCR,N7_PCR,statement="Supplementary Data 2 prints both complete N8 index-primer pools; their 3-prime Nextera arms land exactly on the TSS endpoint library.")
TES_GEOMETRY=PCRGeometry.exact(TES_LONG,N5_PCR,N7_PCR,statement="Supplementary Data 2 prints both complete N8 index-primer pools; their 3-prime Nextera arms land exactly on the TES endpoint library.")
TSS_ENRICH_GEOMETRY=PCRGeometry.unavailable(TSS_LONG,reason="The published enrichment primers are known, but the two intermediate TSS fragment ends are not represented as exact states.")
TES_ENRICH_GEOMETRY=PCRGeometry.unavailable(TES_LONG,reason="The published enrichment primers are known, but the two intermediate TES fragment ends are not represented as exact states.")
ATAC_GEOMETRY=PCRGeometry.exact(ATAC,N5_PCR,N7_PCR,statement="The matched ATAC endpoint uses the same source-printed Nextera index-primer geometry.")
FINAL_LIBRARIES=(
    ("TSS long-fragment library",TSS_LONG,TSS_PRIMERS,"Long TSS fragments from RCA copies retain the source-cell barcode and molecule UMI.","Nextera primers on the TSS long-fragment branch."),
    ("TSS short-fragment library",TSS_SHORT,TSS_PRIMERS,"Short TSS fragments are a separately retained size fraction of the same indexed TSS pool.","Nextera primers on the TSS short-fragment branch."),
    ("TES long-fragment library",TES_LONG,TES_PRIMERS,"Long TES fragments from RCA copies retain the source-cell barcode and molecule UMI.","Nextera primers on the TES long-fragment branch."),
    ("TES short-fragment library",TES_SHORT,TES_PRIMERS,"Short TES fragments are a separately retained size fraction of the same indexed TES pool.","Nextera primers on the TES short-fragment branch."),
    ("Accessibility library",ATAC,AP,"The low-throughput ATAC fraction is paired to its RNA fraction by the source cell.","Nextera primers on the accessibility branch."))

def oligos():
    yield oligo("Oligo(dT) primer",[Segment("SMART handle",DT_HANDLE),Segment("cell barcode","B"*6,"cbc",placeholder=True,feature=CELL),Segment("UMI","U"*12,"umi",placeholder=True,feature=UMI_F),Segment("tag1",TAG1),Segment("oligo(dT)VN",POLYT,placeholder=True)])
    yield oligo("TSO primer",[Segment("TSO",TSO)],mods="/rG//rG//rG/ at 3-prime end")
    for name,seq in (("TSO ISPCR primer",TSO_ISPCR),("linker primer",LINKER),("oligo primer",OLIGO_PRIMER),("Phi29 primer",PHI29_PRIMER),("Tn5-i5 primer",TN5_I5),("Tn5-i7 primer",TN5_I7),("anti-tag1 primer",ANTI_TAG1),("anti-TSO primer",ANTI_TSO),("i7-anti-tag1 primer",I7_ANTI_TAG1),("i7-anti-TSO primer",I7_ANTI_TSO),("Index 1 primer (N7)",INDEX_N7),("Index 2 primer (N5)",INDEX_N5)):
        yield oligo(name,[Segment(name,seq,placeholder="N" in seq)])
def sections():
    return [("Capture RNA and accessible DNA from the same cell",analyte_split("full-length RNA/cDNA branch","ATAC branch"),"Low-throughput separation keeps the two fractions paired by the source cell."),("Label both cDNA ends with molecule identity",indexed_molecule(payload="full-length cDNA",barcode_parts=(("cell barcode",6),),umi=12).rows(),"The printed low-throughput oligo(dT) primer carries a six-base cell barcode, 12-base UMI and tag1 after its SMART handle."),("Circularize full-length cDNA",CIRCLE.rows("NEBuilder overlap assembly"),"The linker is exactly reverse-complementary to the oligo(dT) handle followed by the TSO, establishing the end overlap used for circle assembly."),("Rolling-circle amplify the closed cDNA",RCA.rows("full-length-cDNA concatemer"),"The printed Phi29 primer binds the complementary TSO-derived site; the concatemer is computed from the circle."),("Tagment the amplified circle and build both libraries",[Row(chunks=[("RCA concatemer → random Tn5 fragments → TSS and TES libraries", "me", False)]),Row(chunks=[("accessible genomic DNA → ATAC library", "me", False)])],"Printed end-specific primers enrich TSS or TES fragments before the source-printed N5/N7 terminal PCR."),]


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
    labelled=g.product("labelled","full-length cDNA labelled at both ends",MaterialKind.MOLECULAR,state=indexed_molecule(payload="full-length cDNA",barcode_parts=(("cell barcode",6),),umi=12))
    g.react("label_ends","Label both cDNA ends with cell and molecule identity",StepKind.ADAPTER_ATTACHMENT,((cdna,InputRole.TEMPLATE),(end_tags,InputRole.STRUCTURAL)),(labelled,))
    assembly_mix=g.input("assembly_mix","NEBuilder HiFi DNA Assembly Master Mix",MaterialKind.CATALYST)
    circle=g.product("circle","end-labelled circular full-length cDNA",MaterialKind.MOLECULAR,state=CIRCLE)
    g.react("circularize","Circularize the end-labelled full-length cDNA",StepKind.CIRCULARIZATION,((labelled,InputRole.TEMPLATE),(assembly_mix,InputRole.CATALYST)),(circle,),geometry=CircularizationGeometry.known("overlap assembly",statement="The printed linker is reverse-complementary to the oligo(dT) handle followed by the TSO; the Methods specify NEBuilder HiFi assembly."))
    rca_polymerase=g.input("rca_polymerase","rolling-circle amplification polymerase",MaterialKind.CATALYST)
    rca_primer=g.input("rca_primer","printed Phi29 primer",MaterialKind.OLIGO)
    concatemer=g.product("concatemer","rolling-circle copies of full-length cDNA",MaterialKind.MOLECULAR,state=RCA)
    g.react("rolling_circle","Amplify the circular cDNA by rolling-circle extension",StepKind.EXTENSION,((circle,InputRole.TEMPLATE),(rca_primer,InputRole.STRUCTURAL),(rca_polymerase,InputRole.CATALYST)),(concatemer,))

    rna_transposome=g.input("rna_transposome","source-printed isoform-branch Tn5 adapter primers",MaterialKind.ADAPTER)
    rna_tn5=g.input("rna_tn5","isoform-branch Tn5 transposase",MaterialKind.CATALYST)
    rna_template=g.product("rna_template","tagmented rolling-circle cDNA fragments",MaterialKind.MOLECULAR)
    rna_waste=g.product("rna_tagment_waste","non-library rolling-circle tagmentation products",MaterialKind.WASTE)
    g.react("tagment_rna","Tagment the rolling-circle cDNA copies",StepKind.TAGMENTATION,((concatemer,InputRole.TEMPLATE),(rna_transposome,InputRole.STRUCTURAL),(rna_tn5,InputRole.CATALYST)),(rna_template,rna_waste));g.discard(rna_waste,"does not retain a library-eligible heterologous isoform fragment")
    tss_fragments=g.product("tss_fragments","TSS-directed aliquot of tagmented cDNA",MaterialKind.MOLECULAR)
    tes_fragments=g.product("tes_fragments","TES-directed aliquot of tagmented cDNA",MaterialKind.MOLECULAR)
    g.react("split_ends","Split tagmented cDNA for separate TSS and TES library construction",StepKind.SPLIT,((rna_template,InputRole.TEMPLATE),),(tss_fragments,tes_fragments),changes_molecule=False)
    tn5_i5=g.input("tn5_i5","published Tn5-i5 enrichment primer",MaterialKind.OLIGO)
    anti_tag1=g.input("anti_tag1","published anti-tag1 TSS primer",MaterialKind.OLIGO)
    anti_tso=g.input("anti_tso","published anti-TSO TES primer",MaterialKind.OLIGO)
    i7_anti_tag1=g.input("i7_anti_tag1","published i7-anti-tag1 TSS primer",MaterialKind.OLIGO)
    i7_anti_tso=g.input("i7_anti_tso","published i7-anti-TSO TES primer",MaterialKind.OLIGO)
    enrich_polymerase=g.input("enrich_polymerase","KAPA HiFi HotStart ReadyMix",MaterialKind.CATALYST)
    tss_first=g.product("tss_first","TSS-enriched PCR1 fragments",MaterialKind.MOLECULAR)
    tes_first=g.product("tes_first","TES-enriched PCR1 fragments",MaterialKind.MOLECULAR)
    tss_second=g.product("tss_second","TSS-enriched PCR2 fragments",MaterialKind.MOLECULAR)
    tes_second=g.product("tes_second","TES-enriched PCR2 fragments",MaterialKind.MOLECULAR)
    g.react("tss_pcr1","Enrich TSS fragments with Tn5-i5 and anti-tag1",StepKind.PCR,((tss_fragments,InputRole.TEMPLATE),(tn5_i5,InputRole.STRUCTURAL),(anti_tag1,InputRole.STRUCTURAL),(enrich_polymerase,InputRole.CATALYST)),(tss_first,),geometry=TSS_ENRICH_GEOMETRY)
    g.react("tes_pcr1","Enrich TES fragments with Tn5-i5 and anti-TSO",StepKind.PCR,((tes_fragments,InputRole.TEMPLATE),(tn5_i5,InputRole.STRUCTURAL),(anti_tso,InputRole.STRUCTURAL),(enrich_polymerase,InputRole.CATALYST)),(tes_first,),geometry=TES_ENRICH_GEOMETRY)
    g.react("tss_pcr2","Add the TSS i7-anti-tag1 arm",StepKind.PCR,((tss_first,InputRole.TEMPLATE),(tn5_i5,InputRole.STRUCTURAL),(i7_anti_tag1,InputRole.STRUCTURAL),(enrich_polymerase,InputRole.CATALYST)),(tss_second,),geometry=TSS_ENRICH_GEOMETRY)
    g.react("tes_pcr2","Add the TES i7-anti-TSO arm",StepKind.PCR,((tes_first,InputRole.TEMPLATE),(tn5_i5,InputRole.STRUCTURAL),(i7_anti_tso,InputRole.STRUCTURAL),(enrich_polymerase,InputRole.CATALYST)),(tes_second,),geometry=TES_ENRICH_GEOMETRY)
    n5_index=g.input("n5_index","published N5 index-primer pool",MaterialKind.OLIGO)
    n7_index=g.input("n7_index","published N7 index-primer pool",MaterialKind.OLIGO)
    index_polymerase=g.input("index_polymerase","KAPA HiFi index-PCR polymerase",MaterialKind.CATALYST)
    tss_indexed=g.product("tss_indexed","indexed TSS fragments before size separation",MaterialKind.MOLECULAR,state=TSS_LONG)
    tes_indexed=g.product("tes_indexed","indexed TES fragments before size separation",MaterialKind.MOLECULAR,state=TES_LONG)
    g.react("tss_index_pcr","Index the TSS library",StepKind.PCR,((tss_second,InputRole.TEMPLATE),(n5_index,InputRole.STRUCTURAL),(n7_index,InputRole.STRUCTURAL),(index_polymerase,InputRole.CATALYST)),(tss_indexed,),geometry=TSS_GEOMETRY)
    g.react("tes_index_pcr","Index the TES library",StepKind.PCR,((tes_second,InputRole.TEMPLATE),(n5_index,InputRole.STRUCTURAL),(n7_index,InputRole.STRUCTURAL),(index_polymerase,InputRole.CATALYST)),(tes_indexed,),geometry=TES_GEOMETRY)
    size_beads=g.input("size_beads","AMPure XP two-cut size-selection beads",MaterialKind.SELECTION)
    tss_long=g.product("tss_long","long-fragment TSS library",MaterialKind.LIBRARY,state=TSS_LONG)
    tss_short=g.product("tss_short","short-fragment TSS library",MaterialKind.LIBRARY,state=TSS_SHORT)
    tes_long=g.product("tes_long","long-fragment TES library",MaterialKind.LIBRARY,state=TES_LONG)
    tes_short=g.product("tes_short","short-fragment TES library",MaterialKind.LIBRARY,state=TES_SHORT)
    g.react("select_tss_sizes","Retain long and short TSS size fractions",StepKind.SELECTION,((tss_indexed,InputRole.TEMPLATE),(size_beads,InputRole.SELECTION)),(tss_long,tss_short),changes_molecule=False)
    g.react("select_tes_sizes","Retain long and short TES size fractions",StepKind.SELECTION,((tes_indexed,InputRole.TEMPLATE),(size_beads,InputRole.SELECTION)),(tes_long,tes_short),changes_molecule=False)
    for material,construct,primers in ((tss_long,TSS_LONG,TSS_PRIMERS),(tss_short,TSS_SHORT,TSS_PRIMERS),(tes_long,TES_LONG,TES_PRIMERS),(tes_short,TES_SHORT,TES_PRIMERS)):
        g.sequence(material,SequencingPlan(library=construct,primers=tuple(primers),acquired_roles=tuple(p.role for p in primers)))

    atac_transposome=g.input("atac_transposome","ATAC loaded Tn5 adapters",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (complete loaded-adapter oligos not transcribed)")
    atac_tn5=g.input("atac_tn5","ATAC Tn5 transposase",MaterialKind.CATALYST)
    atac_template=g.product("atac_template","tagmented accessible genomic DNA",MaterialKind.MOLECULAR,state=ATAC,inferred=True,source=SOURCE+" (terminal oligos not transcribed)")
    atac_waste=g.product("atac_tagment_waste","non-library chromatin tagmentation products",MaterialKind.WASTE)
    g.react("tagment_atac","Tagment accessible chromatin",StepKind.TAGMENTATION,((chromatin,InputRole.TEMPLATE),(atac_transposome,InputRole.STRUCTURAL),(atac_tn5,InputRole.CATALYST)),(atac_template,atac_waste));g.discard(atac_waste,"does not retain a library-eligible heterologous ATAC fragment")

    reviewed=["cell","separate","rt_primer","reverse_transcriptase","reverse_transcribe","end_tags","label_ends","assembly_mix","circularize","rca_primer","rca_polymerase","rolling_circle","rna_transposome","rna_tn5","tagment_rna","split_ends","tn5_i5","anti_tag1","anti_tso","i7_anti_tag1","i7_anti_tso","enrich_polymerase","tss_pcr1","tes_pcr1","tss_pcr2","tes_pcr2","n5_index","n7_index","index_polymerase","tss_index_pcr","tes_index_pcr","size_beads","select_tss_sizes","select_tes_sizes","tss_long","tss_short","tes_long","tes_short","atac_transposome","atac_tn5","tagment_atac"]
    finish_pcr_branch(g,"atac",atac_template,ATAC,AP,reviewed,label="scRICA-seq accessibility",geometry=ATAC_GEOMETRY)
    g.annotate(SourceDiscrepancy(
        target="rt_primer", title="Low-throughput barcode length in Supplementary Data 2",
        claims=(
            SourceClaim("Jiang et al. 2024, Supplementary Data 2, Oligo(dT) primer row",
                        "The oligo is printed as (N)6(N)12 before tag1."),
            SourceClaim("Jiang et al. 2024, Supplementary Data 2, table footnote",
                        "The footnote calls (N)6 a 16bp cell barcode."),
            SourceClaim("Jiang et al. 2024, Supplementary Fig. 3, low-throughput chemistry",
                        "The low-throughput drawing labels N6 as a 6bp cell barcode; Supplementary Fig. 2 separately labels high-throughput N16 as 16bp.")),
        model_choice="The low-throughput oligo and isoform library use a six-base cell barcode; the paired low-throughput ATAC library has no molecular cell-barcode segment.",
        evidence="The printed oligo length and low-throughput molecular drawing agree on N6. Cells are paired by physical separation in this route.",
        uncertainty="The supplementary table footnote appears to carry over the high-throughput N16 length; the authors did not issue a correction in the checked sources.",
        slot=AnnotationSlot.INTRO))
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
