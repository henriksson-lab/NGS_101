"""Molecular construct model for Plate_scATAC-seq (Chen et al., 2018)."""
from __future__ import annotations
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,complement_segments,feature,revcomp
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, AnnotationSlot, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)

def seg(name,top,tag=None,**kw): return Segment(name=name,top=top,tag=tag,**kw)

I7_N701="TCGCCTTA"; I5_S502="CTCTCTAT"
N701=il.P7+I7_N701+nx.S7
S502=il.P5+I5_S502+nx.S5

def n701(): return [seg("P7",il.P7,"p7"),seg("i7",I7_N701,"cbc"),seg("s7",nx.S7,"s7")]
def s502(): return [seg("P5",il.P5,"p5"),seg("i5",I5_S502,"cbc"),seg("s5",nx.S5,"s5")]

def tagmented_scene():
    gap=nx.TAGMENTATION_GAP; core=seg("genomic core","X"*24,placeholder=True)
    top=[seg("s5",nx.S5,"s5"),seg("left ME",nx.ME,"me"),seg("left gap","X"*gap,placeholder=True),core]
    bot=[seg("s7",nx.S7,"s7"),seg("right ME",nx.ME,"me"),seg("right gap","X"*gap,placeholder=True),*complement_segments([core])]
    sc=Scene(); sc.strand("top",top,label="nuclear DNA")
    sc.anneal("bottom",bot,to="top",pair=("genomic core'","genomic core"),label="nuclear DNA")
    sc.anneal("left ME bottom",[seg("left ME'",nx.ME_RC,"me")],to="top",pair=("left ME'","left ME"),label="",mod5="p")
    sc.anneal("right ME bottom",[seg("right ME'",nx.ME_RC,"me")],to="bottom",pair=("right ME'","right ME"),label="",mod5="p",above=True)
    sc.mark("top","left gap","9-nt gap"); sc.mark("bottom","right gap","9-nt gap"); return sc

SEQ_PRIMERS=tuple(sp.NEXTERA[k] for k in ("R1","I1","I2","R2"))
def final_library(insert_nt=34):
    lib=Construct([seg("P5",il.P5,"p5"),seg("i5",I5_S502,"cbc",feature=feature("sample_i5","sample_index","whitelist",whitelist="Nextera index set")),seg("s5",nx.S5,"s5"),seg("ME",nx.ME,"me"),
      seg("accessible genomic DNA","X"*insert_nt,placeholder=True),seg("ME reverse complement",nx.ME_RC,"me"),seg("s7 reverse complement",nx.S7_RC,"s7"),
      seg("i7 reverse complement",revcomp(I7_N701),"cbc",feature=feature("sample_i7","sample_index","whitelist",whitelist="Nextera index set")),seg("P7 reverse complement",il.P7_RC,"p7")],name="Plate_scATAC-seq library")
    problems=sp.verify(lib,SEQ_PRIMERS)
    if problems: raise ValueError("invalid Plate_scATAC-seq library: "+"; ".join(problems))
    return lib
def primer_landings():
    lib=final_library(); return [(p,sp.locate(lib,p)) for p in SEQ_PRIMERS]

def protocol_graph():
    source="Defining source: https://doi.org/10.1038/s41467-018-07771-0"
    graph=ProtocolGraph("Plate_scATAC-seq",authority=source)
    nuclei=graph.input("nuclei","native nuclei",MaterialKind.BIOLOGICAL,primary=True)
    transposome=graph.input("transposome","Nextera-loaded Tn5 transposome",MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagmented_nuclei","intact nuclei containing tagmented chromatin",
                         MaterialKind.MOLECULAR,state=tagmented_scene())
    waste=graph.product("tagmentation_waste","candidate non-library tagmentation products",MaterialKind.MOLECULAR)
    graph.react("tagment","Tagment native nuclei in bulk: 37°C for30min",StepKind.TAGMENTATION,
                ((nuclei,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,waste))
    graph.mark_terminal(Terminal(waste.id,Disposition.UNRESOLVED,"Amplification eligibility is not evidence of physical removal of candidate insertion products."))
    sorter=graph.input("single_nucleus_sort","one-nucleus-per-well sorting rule",MaterialKind.SELECTION)
    wells=graph.product("single_nucleus_wells","single tagmented nuclei in indexed wells",MaterialKind.MOLECULAR)
    excluded=graph.product("excluded_events","empty wells and multiplet sorting events",MaterialKind.WASTE)
    graph.react("sort","Sort one tagmented nucleus into each pre-indexed well",StepKind.SELECTION,
                ((tagged,InputRole.TEMPLATE),(sorter,InputRole.SELECTION)),(wells,excluded))
    graph.discard(excluded,"empty and multiplet events are excluded from single-cell libraries")
    release=graph.input("release_mix","SDS, proteinase K and Tween-20 release conditions",MaterialKind.CONDITION)
    released=graph.product("released_fragments","Tn5-released tagmented fragments",MaterialKind.MOLECULAR)
    debris=graph.product("protein_debris","solubilized Tn5 and nuclear components",MaterialKind.MOLECULAR)
    graph.react("release","Release Tn5 and nuclear proteins in each well",StepKind.SOLUBILIZATION,
                ((wells,InputRole.TEMPLATE),(release,InputRole.CONDITION)),(released,debris))
    graph.mark_terminal(Terminal(debris.id,Disposition.UNRESOLVED,"Solubilization does not itself establish a physical discarded fraction; detailed recovery remains deferred."))
    polymerase=graph.input("gap_fill_polymerase","gap-fill polymerase",MaterialKind.CATALYST)
    filled=graph.product("gap_filled_fragments","gap-filled tagmented fragments",MaterialKind.MOLECULAR)
    graph.react("gap_fill","Fill the transposition gaps at 72 degrees C",StepKind.EXTENSION,
                ((released,InputRole.TEMPLATE),(polymerase,InputRole.CATALYST)),(filled,))
    s5=graph.input("s5_index_primer","well-specific S5xx primer",MaterialKind.OLIGO,state=tuple(s502()))
    n7=graph.input("n7_index_primer","well-specific N7xx primer",MaterialKind.OLIGO,state=tuple(n701()))
    library=graph.product("library",final_library().name,MaterialKind.LIBRARY,state=final_library())
    graph.react("indexed_pcr","Amplify with well-specific S5xx/N7xx: 18 total cycles in this study",StepKind.PCR,
                ((filled,InputRole.TEMPLATE),(s5,InputRole.STRUCTURAL),
                 (n7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(final_library(),
                    reason="The gap-filled pre-PCR construct is drawn as a molecular scene; exact final Nextera arms are verified separately."))
    graph.sequence(library,SequencingPlan(
        library=final_library(),primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_lengths={"Index 1 (i7)":8,"Index 2 (i5)":8},
        read_uses={"Read 1":"include","Index 1 (i7)":"demultiplex",
                   "Index 2 (i5)":"demultiplex","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "nuclei","transposome","tn5","tagment","single_nucleus_sort","sort",
        "release_mix","release","gap_fill_polymerase","gap_fill","s5_index_primer",
        "n7_index_primer","indexed_pcr","library"),source)
    graph.annotate(SourceDiscrepancy(
        target="gap_fill", title="Published gap-fill timing difference",
        claims=(
            SourceClaim("Chen et al. main text, qPCR optimization", "Gap fill at 72 °C for 5 min."),
            SourceClaim("Chen et al. Supplementary Methods, library preparation", "Gap fill at 72 °C for 10 min.")),
        model_choice="The page follows the 10 min library-preparation instruction.",
        evidence="The supplementary step is the production library recipe represented by this graph.",
        uncertainty="The sources do not establish whether the qPCR and production timings differ intentionally.",
        slot=AnnotationSlot.BEFORE_WORKFLOW))
    return graph.require_publishable()
