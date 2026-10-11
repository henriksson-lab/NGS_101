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

TITLE = "Single-stranded library preparation for ancient or damaged DNA"
NOTES = "01_single-stranded-library.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nprot.2013.038">Gansauge & Meyer, <i>Nature Protocols</i> (2013)</a>.'
SUMMARY = "Heat-denatured individual DNA strands receive a biotinylated 3′ adapter, are captured on streptavidin beads, copied into duplexes, blunt-ligated to a second adapter and completed by indexed PCR."
CAVEAT = "Modified oligos are shown with their published modifications. The damaged insert remains a placeholder because every recovered strand can have a different length and lesion pattern."

CL72 = "ACACTCTTTCCCTACACGACGCTCTTCC"
GES_INDEX2 = "GGAAGAGCGTCGTGTAGGGAAAGAGTGT"
FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"), seg("i5", "J" * 8, "cbc", placeholder=True,
                                  feature=feature("sample_i5", "sample_index", "unknown")),
    seg("CL72 custom Read 1 site", CL72, "r1"),
    seg("ancient DNA insert", "X" * 34, placeholder=True), seg("dA junction", "A"),
    seg("CL9 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
    seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
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
    sc=Scene(); sc.strand("CL78",[seg("5′ phosphorylated CL78 ligating end","X"*12,placeholder=True),seg("C3 spacers + TEG-biotin","[C3 spacer]10-[TEG-biotin]",placeholder=True)],label="published CL78")
    sc.ligation_end("CL78","5′ phosphorylated CL78 ligating end",kind="single-strand",polarity="5-prime",owner="adapter",text="5′-p single-strand ligating end")
    return AdapterLigationGeometry.known("single-strand",sc,statement="The protocol specifies CircLigase joining the input 3′-OH to the phosphorylated 5′ end of CL78; this is not a duplex overhang.")

def _cl53_cl73_geometry():
    sc=Scene.duplex([seg("CL53/CL73 paired end","X"*20,placeholder=True)],label="CL53/CL73 adapter")
    sc.ligation_end("top","CL53/CL73 paired end",kind="blunt",owner="adapter",text="blunt adapter end")
    return AdapterLigationGeometry.known("blunt",sc,statement="The protocol explicitly blunt-ligates the CL53/CL73 duplex to the copied ancient-DNA duplex.")

def workflow():
    ss=Scene(); ss.strand("damaged strand",[seg("damaged DNA","X"*30,placeholder=True)],label="heat-denatured input")
    w=Workflow(MolecularState("heat-denatured ancient DNA",tuple(ss.rows())))
    joined=Scene(); joined.strand("captured strand",[seg("damaged DNA","X"*30,placeholder=True),seg("CL78","X"*16,placeholder=True)],label="CL78-ligated strand"); joined.junction("captured strand","damaged DNA","CL78","CircLigase")
    w.ligate_adapter("Ligate CL78 to individual strands",_cl78_geometry(),joined.rows(),name="CL78-ligated ancient DNA",note="The terminal biotin enables capture.")
    w.react("Immobilize and copy the captured strand",[Row(chunks=[("streptavidin bead—biotin—CL78—ancient-DNA duplex","cbc",False)])],name="copied ancient-DNA duplex",note="The bead retains the adapter-bearing template through washes.")
    product=Scene.duplex([seg("CL53/CL73 adapter","X"*18,placeholder=True),seg("copied ancient DNA","X"*30,placeholder=True),seg("first-adapter end","X"*16,placeholder=True)],label="dual-adapter molecule"); product.junction("top","CL53/CL73 adapter","copied ancient DNA","blunt ligation")
    w.ligate_adapter("Blunt-ligate the second adapter",_cl53_cl73_geometry(),product.rows(),name="dual-adapter ancient-DNA duplex",note="Indexed PCR releases and completes the sequenceable molecule.")
    return w
