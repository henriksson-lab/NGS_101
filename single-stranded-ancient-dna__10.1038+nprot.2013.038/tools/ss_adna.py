"""Gansauge & Meyer single-stranded ancient-DNA library model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Row, Scene, Workflow, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "Single-stranded library preparation for ancient or damaged DNA"
NOTES = "01_single-stranded-library.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nprot.2013.038">Gansauge & Meyer, <i>Nature Protocols</i> (2013)</a>.'
SUMMARY = "Heat-denatured individual DNA strands receive a biotinylated 3′ adapter, are captured on streptavidin beads, copied into duplexes, blunt-ligated to a second adapter and completed by indexed PCR."
CAVEAT = "Modified oligos are shown with their published modifications. The damaged insert remains a placeholder because every recovered strand can have a different length and lesion pattern."

CL72 = "ACACTCTTTCCCTACACGACGCTCTTCC"
GES_INDEX2 = "GGAAGAGCGTCGTGTAGGGAAAGAGTGT"
FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"), seg("i5", "J" * 7, "cbc", placeholder=True,
                                  feature=feature("sample_i5", "sample_index", "unknown")),
    seg("CL72 custom Read 1 site", CL72, "r1"),
    seg("ancient DNA insert", "X" * 34, placeholder=True), seg("dA junction", "A"),
    seg("CL9 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
    seg("i7 reverse complement", "I" * 7, "cbc", placeholder=True,
        feature=feature("sample_i7", "sample_index", "unknown")),
    seg("P7 reverse complement", il.P7_RC, "p7")], name="single-stranded ancient-DNA library")
SEQ_PRIMERS = (sp.custom("Read 1", "CL72", CL72, "Gansauge & Meyer 2013"),
               sp.TRUSEQ["I1"],
               sp.custom("Index 2 (i5)", "Gesaffelstein index 2", GES_INDEX2,
                         "Gansauge single-strand library sequencing setup"),
               sp.TRUSEQ["R2"])
problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                     required_roles=tuple(p.role for p in SEQ_PRIMERS))
if problems: raise ValueError("ancient-DNA library: " + "; ".join(problems))
FINAL_CAPTION = "Indexed PCR supplies full flow-cell arms around the recovered ancient strand. The P5-side primer site is deliberately truncated and uses the published CL72 custom Read 1 primer."
SEQUENCING_INTRO = "CL72 replaces the standard Read 1 primer because the P5-side library arm is truncated. The read begins directly in the captured ancient DNA strand."

def _cl78_geometry():
    sc=Scene(); sc.strand("CL78",[seg("5′ phosphorylated CL78 ligating end","AGATCGGAAG"),seg("C3 spacers + TEG-biotin","[C3 spacer]10-[TEG-biotin]",placeholder=True)],label="published CL78")
    sc.ligation_end("CL78","5′ phosphorylated CL78 ligating end",kind="single-strand",polarity="5-prime",owner="adapter",text="5′-p single-strand ligating end")
    return AdapterLigationGeometry.known("single-strand",sc,statement="The protocol specifies CircLigase joining the input 3′-OH to the phosphorylated 5′ end of CL78; this is not a duplex overhang.")

def _cl53_cl73_geometry():
    sc=Scene.duplex([seg("CL53/CL73 paired end","X"*20,placeholder=True)],label="CL53/CL73 adapter")
    sc.ligation_end("top","CL53/CL73 paired end",kind="blunt",owner="adapter",text="blunt adapter end")
    return AdapterLigationGeometry.known("blunt",sc,statement="The protocol explicitly blunt-ligates the CL53/CL73 duplex to the copied ancient-DNA duplex.")

def workflow():
    ss=Scene(); ss.strand("damaged strand",[seg("damaged DNA","X"*30,placeholder=True)],label="heat-denatured input")
    w=Workflow(MolecularState("heat-denatured ancient DNA",tuple(ss.rows())))
    joined=Scene(); joined.strand("captured strand",[seg("damaged DNA","X"*30,placeholder=True),seg("CL78","AGATCGGAAG"),seg("C3 spacers + TEG-biotin","[C3 spacer]10-[TEG-biotin]",placeholder=True)],label="CL78-ligated strand"); joined.junction("captured strand","damaged DNA","CL78","CircLigase")
    w.ligate_adapter("Ligate CL78 to individual strands",_cl78_geometry(),joined.rows(),name="CL78-ligated ancient DNA",note="The terminal biotin enables capture.")
    w.react("Immobilize and copy the captured strand",[Row(chunks=[("streptavidin bead—biotin—CL78—ancient-DNA duplex","cbc",False)])],name="copied ancient-DNA duplex",note="The bead retains the adapter-bearing template through washes.")
    product=Scene.duplex([seg("CL53/CL73 adapter","X"*18,placeholder=True),seg("copied ancient DNA","X"*30,placeholder=True),seg("first-adapter end","X"*16,placeholder=True)],label="dual-adapter molecule"); product.junction("top","CL53/CL73 adapter","copied ancient DNA","blunt ligation")
    w.ligate_adapter("Blunt-ligate the second adapter",_cl53_cl73_geometry(),product.rows(),name="dual-adapter ancient-DNA duplex",note="Indexed PCR releases and completes the sequenceable molecule.")
    return w


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("damaged_dna","purified damaged or ancient DNA",
                    MaterialKind.MOLECULAR,primary=True)
    heat=graph.input("denaturation_heat","heat-denaturation condition",
                     MaterialKind.CONDITION)
    strand_a=graph.product("strand_a","individual denatured ancient-DNA strand A",
                           MaterialKind.MOLECULAR)
    strand_b=graph.product("strand_b","individual denatured ancient-DNA strand B",
                           MaterialKind.MOLECULAR)
    graph.react("denature","Heat-denature damaged DNA into individual strands",
                StepKind.STRAND_SEPARATION,
                ((dna,InputRole.TEMPLATE),(heat,InputRole.CONDITION)),
                (strand_a,strand_b))
    strands=graph.product("single_strands","pool of independently recoverable DNA strands",
                          MaterialKind.MOLECULAR)
    graph.react("collect_strands","Collect both denatured-strand populations",
                StepKind.POOL,
                ((strand_a,InputRole.TEMPLATE),(strand_b,InputRole.TEMPLATE)),
                (strands,),changes_molecule=False)
    cl78=graph.input("cl78","phosphorylated biotinylated CL78 adapter",
                     MaterialKind.ADAPTER)
    circligase=graph.input("circligase","CircLigase",MaterialKind.CATALYST)
    adapted=graph.product("cl78_ligated","CL78-ligated ancient-DNA strands",
                          MaterialKind.MOLECULAR)
    graph.react("ligate_cl78","Ligate CL78 to individual strand 3-prime ends",
                StepKind.LIGATION,
                ((strands,InputRole.TEMPLATE),(cl78,InputRole.STRUCTURAL),
                 (circligase,InputRole.CATALYST)),(adapted,),geometry=_cl78_geometry())
    streptavidin=graph.input("streptavidin","streptavidin capture beads",
                            MaterialKind.SELECTION)
    captured=graph.product("captured","bead-captured CL78-bearing strands",
                           MaterialKind.MOLECULAR)
    uncaptured=graph.product("uncaptured","unligated or unbiotinylated DNA",MaterialKind.WASTE)
    graph.react("capture","Capture biotinylated CL78-bearing strands",
                StepKind.AFFINITY_CAPTURE,
                ((adapted,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),
                (captured,uncaptured))
    graph.discard(uncaptured,"does not carry the terminal biotinylated CL78 adapter")
    cl9=graph.input("cl9","CL9 extension primer",MaterialKind.OLIGO)
    copy_polymerase=graph.input("copy_polymerase","strand-copying DNA polymerase",
                                MaterialKind.CATALYST)
    duplex=graph.product("copied_duplex","copied ancient-DNA duplex",
                         MaterialKind.MOLECULAR)
    graph.react("copy_strand","Copy the captured ancient-DNA strand into a duplex",
                StepKind.EXTENSION,
                ((captured,InputRole.TEMPLATE),(cl9,InputRole.STRUCTURAL),
                 (copy_polymerase,InputRole.CATALYST)),(duplex,))
    second_adapter=graph.input("second_adapter","CL53/CL73 blunt duplex adapter",
                               MaterialKind.ADAPTER)
    ligase=graph.input("blunt_ligase","blunt-end DNA ligase",MaterialKind.CATALYST)
    dual=graph.product("dual_adapter","dual-adapter ancient-DNA duplex",
                       MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    graph.react("ligate_second","Blunt-ligate the second adapter",
                StepKind.LIGATION,
                ((duplex,InputRole.TEMPLATE),(second_adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(dual,),geometry=_cl53_cl73_geometry())
    forward=graph.input("pcr_forward","library PCR forward primer",MaterialKind.OLIGO)
    reverse=graph.input("pcr_reverse","library PCR reverse/index primer",MaterialKind.OLIGO)
    pcr_polymerase=graph.input("pcr_polymerase","library PCR polymerase",
                               MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    graph.react("index_pcr","Release, amplify and index the captured library",
                StepKind.PCR,
                ((dual,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(pcr_polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The page models the exact sequencing endpoint, but does not assign one index-primer member as a universal PCR reagent."))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "damaged_dna","denaturation_heat","denature","collect_strands","cl78",
        "circligase","ligate_cl78","streptavidin","capture","cl9",
        "copy_polymerase","copy_strand","second_adapter","blunt_ligase",
        "ligate_second","pcr_forward","pcr_reverse","pcr_polymerase",
        "index_pcr","library"),SOURCE)
    return graph.require_publishable()
