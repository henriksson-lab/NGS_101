from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import joined_scene,seg,truseq_library
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState,Row,Scene,Segment,Workflow,feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE="GUIDE-seq — cellular capture of nuclease-cut sites"
NOTES="01_guide-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nbt.3117">Tsai et al., <i>Nature Biotechnology</i> (2015)</a>.'
SUMMARY="A protected, phosphorylated blunt dsODN integrates into nuclease-created breaks in living cells; single-tail adapter/tag PCR then recovers genomic DNA immediately beside either dsODN end."
CAVEAT="The two STAT-PCR orientations are parallel libraries from the two dsODN strands. The schematic shows one orientation and names the mirrored reaction explicitly."
UMI=feature("guide_seq_umi","umi","random")
FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([seg("8-nt molecular barcode","U"*8,"umi",placeholder=True,feature=UMI),seg("integrated dsODN end","D"*18,placeholder=True),seg("adjacent genomic DNA","X"*42,placeholder=True)],"GUIDE-seq STAT-PCR library")
FINAL_CAPTION="STAT-PCR joins one dsODN-specific end to a single-tail adapter, while an 8-nt molecular barcode distinguishes captured molecules."
SEQUENCING_INTRO="Read 1 enters the dsODN/genome junction; the opposite read and index sites are supplied by nested library PCR."

def integration_scene():
    left=[Segment("genome left","X"*24,placeholder=True),Segment("dsODN","D"*34,"cbc",placeholder=True),Segment("genome right","X"*24,placeholder=True)]
    sc=Scene.duplex(left,label="NHEJ-tagged break")
    sc.junction("top","genome left","dsODN","NHEJ capture")
    sc.junction("top","dsODN","genome right","NHEJ capture")
    sc.labels("top"); return sc

def sections(): return [
 ("Integrate the protected blunt dsODN at cellular DSBs",integration_scene().rows(),"The 34-bp dsODN is 5-prime phosphorylated and terminally phosphorothioate-protected; NHEJ captures it at on- and off-target breaks."),
 ("Shear genomic DNA and ligate a single-tail adapter",joined_scene([seg("single-tail adapter","X"*20,"r1",placeholder=True),seg("random genomic shear fragment","X"*34,placeholder=True),seg("integrated dsODN","D"*18,"cbc",placeholder=True)],(("single-tail adapter","random genomic shear fragment","adapter ligation"),),label="single-tail capture molecule").rows(),"The adapter has only one amplifiable tail, suppressing background from fragments lacking the dsODN tag."),
 ("Run dsODN-to-adapter STAT-PCR in both orientations",[Row(chunks=[("dsODN strand A primer ---> flank A    |    dsODN strand B primer ---> flank B", "r1",False)])],"Separate reactions recover genomic sequence from both sides of the integrated tag."),
 ("Nested PCR adds UMI, index and complete Illumina arms",[Row(chunks=[("P5 -- UMI -- dsODN -- genomic flank -- P7",None,False)])],"The random 8-bp molecular barcode supports PCR-bias correction before sequencing."),
]

def workflow():
 s=sections();w=Workflow(MolecularState(s[0][0],tuple(s[0][1])))
 g=AdapterLigationGeometry.unavailable("The GUIDE-seq paper and supplement identify a single-tail adapter, but its cut-facing strands and prepared shear-end geometry have not yet been verified for this schematic.")
 w.ligate_adapter(s[1][0],g,s[1][1],note=s[1][2]);w.react(s[2][0],s[2][1],note=s[2][2]);w.react(s[3][0],s[3][1],note=s[3][2]);return w


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    cells=graph.input("cells","living cells with genomic targets",MaterialKind.BIOLOGICAL,
                      primary=True)
    guide=graph.input("guide_rna","targeting guide RNA",MaterialKind.OLIGO)
    nuclease=graph.input("programmable_nuclease","programmable nuclease",
                         MaterialKind.CATALYST)
    cut=graph.product("cut_genome","cells containing nuclease-created DSBs",
                      MaterialKind.MOLECULAR)
    graph.react("cut","Create on- and off-target genomic DSBs",StepKind.GENOME_EDITING,
                ((cells,InputRole.TEMPLATE),(guide,InputRole.STRUCTURAL),
                 (nuclease,InputRole.CATALYST)),(cut,))
    dsodn=graph.input("dsodn","protected phosphorylated blunt dsODN tag",
                      MaterialKind.ADAPTER)
    nhej=graph.input("nhej","cellular NHEJ machinery",MaterialKind.CATALYST)
    tagged=graph.product("tagged_genome","dsODN-tagged genomic break junctions",
                         MaterialKind.MOLECULAR,state=integration_scene())
    graph.react("integrate","Capture the dsODN at cellular DSBs",StepKind.INTEGRATION,
                ((cut,InputRole.TEMPLATE),(dsodn,InputRole.STRUCTURAL),
                 (nhej,InputRole.CATALYST)),(tagged,))
    shear=graph.input("shear_condition","genomic DNA shearing condition",
                      MaterialKind.CONDITION)
    fragments=graph.product("fragments","sheared dsODN-tagged genomic fragments",
                            MaterialKind.MOLECULAR)
    graph.react("shear","Shear extracted genomic DNA",StepKind.FRAGMENTATION,
                ((tagged,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),(fragments,))
    adapter=graph.input("single_tail_adapter","single-tail capture adapter",
                        MaterialKind.ADAPTER,inferred=True,
                        source=SOURCE+" (cut-facing strands not verified)")
    ligase=graph.input("adapter_ligase","single-tail adapter ligase",
                       MaterialKind.CATALYST)
    adapted=graph.product("adapted","single-tail-adapted dsODN junction fragments",
                          MaterialKind.MOLECULAR)
    graph.react("ligate_adapter","Ligate the single-tail adapter",StepKind.LIGATION,
                ((fragments,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The source identifies the adapter but its cut-facing strands and prepared shear-end geometry are not verified."))
    flank_a=graph.product("flank_a","aliquot for dsODN strand-A STAT-PCR",
                          MaterialKind.MOLECULAR)
    flank_b=graph.product("flank_b","aliquot for dsODN strand-B STAT-PCR",
                          MaterialKind.MOLECULAR)
    graph.react("split_orientations","Split into the two junction orientations",
                StepKind.SPLIT,((adapted,InputRole.TEMPLATE),),(flank_a,flank_b),
                changes_molecule=False)
    adapter_primer=graph.input("adapter_primer","single-tail adapter PCR primer",
                               MaterialKind.OLIGO)
    ds_a=graph.input("dsodn_a_primer","dsODN strand-A primer",MaterialKind.OLIGO)
    ds_b=graph.input("dsodn_b_primer","dsODN strand-B primer",MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","STAT-PCR polymerase",MaterialKind.CATALYST)
    amp_a=graph.product("amp_a","strand-A genomic-flank amplicons",MaterialKind.MOLECULAR)
    amp_b=graph.product("amp_b","strand-B genomic-flank amplicons",MaterialKind.MOLECULAR)
    graph.react("stat_a","Amplify genomic flank A",StepKind.PCR,
                ((flank_a,InputRole.TEMPLATE),(adapter_primer,InputRole.STRUCTURAL),
                 (ds_a,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (amp_a,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The variable genomic shear fragment prevents one fixed primer-site Construct."))
    graph.react("stat_b","Amplify genomic flank B",StepKind.PCR,
                ((flank_b,InputRole.TEMPLATE),(adapter_primer,InputRole.STRUCTURAL),
                 (ds_b,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (amp_b,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The variable genomic shear fragment prevents one fixed primer-site Construct."))
    pooled=graph.product("pooled","pooled two-orientation STAT-PCR products",
                         MaterialKind.MOLECULAR)
    graph.react("pool","Pool both captured break orientations",StepKind.POOL,
                ((amp_a,InputRole.TEMPLATE),(amp_b,InputRole.TEMPLATE)),(pooled,),
                changes_molecule=False)
    nested_f=graph.input("nested_forward","nested UMI-bearing forward primer",
                         MaterialKind.OLIGO)
    nested_r=graph.input("nested_reverse","nested indexed reverse primer",
                         MaterialKind.OLIGO)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    graph.react("nested_pcr","Add UMI, index and complete Illumina arms",
                StepKind.PCR,
                ((pooled,InputRole.TEMPLATE),(nested_f,InputRole.STRUCTURAL),
                 (nested_r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The nested endpoint is exact, but the two variable STAT-PCR precursor orientations are not one fixed duplex."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "cells","guide_rna","programmable_nuclease","cut","dsodn","nhej",
        "integrate","shear_condition","shear","single_tail_adapter",
        "adapter_ligase","ligate_adapter","split_orientations","adapter_primer",
        "dsodn_a_primer","dsodn_b_primer","pcr_polymerase","stat_a","stat_b",
        "pool","nested_forward","nested_reverse","nested_pcr","library"),SOURCE)
    return graph.require_publishable()
