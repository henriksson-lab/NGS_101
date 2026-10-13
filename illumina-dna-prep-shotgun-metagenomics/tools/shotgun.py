"""Illumina DNA Prep shotgun-metagenomics model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from targeted_ngs import shotgun_library
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

TITLE="Illumina DNA Prep shotgun metagenomics"
NOTES="01_shotgun-metagenomics.html"
SOURCE='Illumina <a href="https://support.illumina.com/downloads/illumina-dna-prep-reference-guide-1000000025416.html">DNA Prep Reference Guide 1000000025416 v12</a> and <a href="https://www.illumina.com/content/dam/illumina/gcs/assembled-assets/marketing-literature/novaseq-shotgun-metagenomics-app-note-m-na-00010/novaseq-shotgun-metagenomics-app-note-m-na-00010.pdf">shotgun-metagenomics application note M-NA-00010 v2.0</a>.'
SUMMARY="Total community DNA is tagmented by bead-linked transposomes and index-PCR-completed into a dual-index Nextera library without locus selection."
CAVEAT="Bead-linked transposome reagent composition is proprietary; its public molecular consequence—Nextera adapter ends on fragmented DNA—is modeled exactly."
FINAL_LIBRARY,SEQ_PRIMERS=shotgun_library()
FINAL_CAPTION="A representative community-DNA fragment after bead-linked tagmentation and dual-index PCR."
SEQUENCING_INTRO="Standard Nextera primers read a random community fragment and both sample indexes."

def sections():
 return [
  ("Extract total community DNA",[
   Row(chunks=[("bacteria + archaea + fungi + viruses → pooled genomic DNA",None,False)]),
   Row(chunks=[("no marker-gene primer; no locus selection",None,False)]),
  ],"Extraction determines which community molecules reach library preparation."),
  ("Tagment with bead-linked transposomes",[
   Row(chunks=[("BLT bead ●—Tn5 + community dsDNA → bead-bound tagged fragments", "me", False)]),
   Row(chunks=[("wash bead → remove unbound DNA and reaction components",None,False)]),
  ],"Tn5 fragmentation and adapter transfer occur together; bead immobilization permits washing and saturation-based normalization."),
  ("Index PCR and pool",[
   Row(chunks=[("tagged fragment → reduced-cycle PCR → P5—i5 … insert … i7—P7", "cbc", False)]),
  ],"PCR completes the platform arms and labels the sample, not the organism or molecule."),
 ]


def protocol_graph():
 graph=ProtocolGraph(TITLE,authority=SOURCE)
 dna=graph.input("community_dna","extracted total community DNA",MaterialKind.BIOLOGICAL,state=molecule("community DNA",(role("pooled microbial genomes"),)),primary=True)
 adapter=graph.input("blt_adapter","bead-linked Nextera adapter transposomes",MaterialKind.ADAPTER)
 tn5=graph.input("tn5","bead-linked Tn5 transposase",MaterialKind.CATALYST)
 tagged=graph.product("tagged","bead-bound tagged community-DNA fragments",MaterialKind.MOLECULAR,state=dna.state)
 waste=graph.product("tagmentation_waste","unbound and unusable tagmentation material",MaterialKind.WASTE)
 graph.react("tagment","Tagment community DNA with bead-linked transposomes",StepKind.TAGMENTATION,((dna,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,waste))
 graph.discard(waste,"does not form a retained bead-bound tagged fragment")
 wash=graph.input("wash_beads","BLT wash and retention matrix",MaterialKind.SELECTION)
 retained=graph.product("retained","washed bead-bound tagged fragments",MaterialKind.MOLECULAR,state=tagged.state)
 washed=graph.product("wash_waste","removed unbound DNA and reaction components",MaterialKind.WASTE)
 graph.react("wash","Wash bead-bound tagged fragments",StepKind.CLEANUP,((tagged,InputRole.TEMPLATE),(wash,InputRole.SELECTION)),(retained,washed),changes_molecule=False)
 graph.discard(washed,"is removed by the post-tagmentation bead wash")
 i5=graph.input("i5_primer","DNA Prep i5 index primer",MaterialKind.OLIGO)
 i7=graph.input("i7_primer","DNA Prep i7 index primer",MaterialKind.OLIGO)
 polymerase=graph.input("polymerase","reduced-cycle index-PCR polymerase",MaterialKind.CATALYST)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("index_pcr","Complete adapters and index the sample",StepKind.PCR,((retained,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),(i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason="The proprietary bead-bound pre-PCR population is not an exact Construct; the public final Nextera sites are checked."))
 graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS)))
 graph.declare_reviewed_scope(["community_dna","blt_adapter","tn5","tagment","wash_beads","wash","i5_primer","i7_primer","polymerase","index_pcr","library"],SOURCE)
 return graph.require_publishable()
