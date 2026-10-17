"""Illumina DNA Prep shotgun-metagenomics model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import nextera as nx
from chemdraw import Scene, Segment
from batch_ngs import nextera_library, seg
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, StepKind)
from protocol_model import molecule, role

def shotgun_library():
    return nextera_library([seg("metagenomic DNA fragment", "N" * 36,
                                placeholder=True)],
                           "Illumina DNA Prep shotgun library", i5_nt=10, i7_nt=10)


TITLE="Illumina DNA Prep shotgun metagenomics"
NOTES="01_shotgun-metagenomics.html"
SOURCE='Illumina <a href="https://support.illumina.com/downloads/illumina-dna-prep-reference-guide-1000000025416.html">DNA Prep Reference Guide 1000000025416 v12</a> and <a href="https://www.illumina.com/content/dam/illumina/gcs/assembled-assets/marketing-literature/novaseq-shotgun-metagenomics-app-note-m-na-00010/novaseq-shotgun-metagenomics-app-note-m-na-00010.pdf">shotgun-metagenomics application note M-NA-00010 v2.0</a>.'
SUMMARY="Total community DNA (10 ng in the application note) is tagmented by bead-linked transposomes, stopped and washed on the bead, index-PCR-completed with IDT for Illumina UD Set A adapters and bead-purified into a dual-index Nextera library without locus selection."
CAVEAT="Bead-linked transposome reagent composition is proprietary; its public molecular consequence—Nextera adapter ends on fragmented DNA—is modeled exactly."
FINAL_LIBRARY,SEQ_PRIMERS=shotgun_library()
FINAL_CAPTION="A representative community-DNA fragment after bead-linked tagmentation and dual-index PCR."
SEQUENCING_INTRO="Standard Nextera primers read a random community fragment and both sample indexes."

# Application note M-NA-00010 v2.0 lines 215-221: 10 ng community DNA in 30 µl, Illumina
# DNA Prep (M) Tagmentation, IDT for Illumina DNA/RNA UD Indexes Set A; lines 239-242:
# 2 x 251 bp on HiSeq 2500 Rapid Run and NovaSeq 6000 SP.  Guide p.26: UD index codes
# are 10 bp and the run must be configured for 10 bp index reads.  Index read cycles are
# taken from the guide, not printed in the note.
INPUT_NG = 10
# Guide 1000000025416 v12 p.27 BLT PCR table row "10-24 ng -> 8 cycles"; the note does not
# print its cycle number, so this is a derived assignment.
PCR_CYCLES_10_TO_24_NG = 8
READ_LENGTHS = {"Read 1": 251, "Index 1 (i7)": 10, "Index 2 (i5)": 10, "Read 2": 251}


def _tagged_state():
 """Public Nextera transposon ends on a community fragment (BLT composition proprietary)."""
 parts=[Segment("S5",nx.S5,"s5"),Segment("left mosaic end",nx.ME,"me"),
        Segment("community DNA fragment","N"*36,placeholder=True),
        Segment("right mosaic end reverse complement",nx.ME_RC,"me"),
        Segment("S7 reverse complement",nx.S7_RC,"s7")]
 sc=Scene.duplex(parts,label="tagmented community DNA")
 sc.junction("top","left mosaic end","community DNA fragment","transposition boundary")
 sc.junction("top","community DNA fragment","right mosaic end reverse complement","transposition boundary")
 return sc


def protocol_graph():
 graph=ProtocolGraph(TITLE,authority=SOURCE)
 dna=graph.input("community_dna","extracted total community DNA (10 ng in 30 µl)",MaterialKind.BIOLOGICAL,state=molecule("community DNA",(role("pooled microbial genomes"),)),primary=True)
 adapter=graph.input("blt_adapter","bead-linked Nextera adapter transposomes (BLT)",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" Public endpoint architecture is used; BLT reagent composition is proprietary.")
 tn5=graph.input("tn5","bead-linked Tn5 transposase",MaterialKind.CATALYST)
 tag=graph.input("tag_program","TB1; TAG program 55°C 15 min",MaterialKind.CONDITION)
 tagged=graph.product("tagged","bead-bound tagged community-DNA fragments",MaterialKind.MOLECULAR,state=_tagged_state())
 waste=graph.product("tagmentation_waste","non-library or same-end tagmentation products",MaterialKind.WASTE)
 graph.react("tagment","Tagment community DNA with bead-linked transposomes",StepKind.TAGMENTATION,((dna,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST),(tag,InputRole.CONDITION)),(tagged,waste))
 graph.discard(waste,"does not form a heterologous amplifiable library fragment")
 tsb=graph.input("tsb","Tagment Stop Buffer (TSB); PTC program 37°C 15 min",MaterialKind.CONDITION)
 stopped=graph.product("stopped","stopped bead-bound tagged fragments",MaterialKind.MOLECULAR,state=tagged.state)
 graph.react("stop_tagmentation","Stop tagmentation with TSB",StepKind.QUENCHING,((tagged,InputRole.TEMPLATE),(tsb,InputRole.CONDITION)),(stopped,),changes_molecule=False)
 wash=graph.input("wash_beads","magnet: discard supernatant; TWB 100 µl × 3, last TWB kept until PCR",MaterialKind.CONDITION)
 retained=graph.product("retained","washed bead-bound tagged fragments",MaterialKind.MOLECULAR,state=tagged.state)
 washed=graph.product("wash_waste","supernatant and TWB: unbound DNA and reaction components",MaterialKind.WASTE)
 graph.react("wash","Wash bead-bound tagged fragments",StepKind.CLEANUP,((stopped,InputRole.TEMPLATE),(wash,InputRole.CONDITION)),(retained,washed),changes_molecule=False)
 graph.discard(washed,"is removed with the post-tagmentation supernatant and TWB washes")
 i5=graph.input("i5_primer","i5 index adapter (IDT for Illumina DNA/RNA UD Indexes Set A well)",MaterialKind.OLIGO)
 i7=graph.input("i7_primer","i7 index adapter (same pre-paired Set A well)",MaterialKind.OLIGO)
 polymerase=graph.input("polymerase","Enhanced PCR Mix (EPM)",MaterialKind.CATALYST)
 program=graph.input("pcr_program",f"BLT PCR: 68°C 3 min, 98°C 3 min, {PCR_CYCLES_10_TO_24_NG} × (98°C 45 s, 62°C 30 s, 68°C 2 min), 68°C 1 min (cycle number derived from 10 ng input)",MaterialKind.CONDITION)
 amplified=graph.product("amplified","amplified dual-index library (unpurified)",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
 graph.react("index_pcr","Complete adapters and index the sample",StepKind.PCR,((retained,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),(i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST),(program,InputRole.CONDITION)),(amplified,),geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason="The proprietary bead-bound pre-PCR population is not an exact Construct and the Set A well is unassigned; the public final Nextera sites are checked."))
 ipb=graph.input("ipb","Illumina Purification Beads (IPB), double-sided cleanup",MaterialKind.SELECTION)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 large=graph.product("large_fraction","first-bead fraction (larger fragments)",MaterialKind.WASTE)
 small=graph.product("small_fraction","second supernatant (short products and reaction components)",MaterialKind.WASTE)
 graph.react("ipb_cleanup","Double-sided IPB cleanup, ethanol washes, elute in RSB",StepKind.CLEANUP,((amplified,InputRole.TEMPLATE),(ipb,InputRole.SELECTION)),(library,large,small),changes_molecule=False)
 graph.discard(large,"is discarded with the first bead plate")
 graph.discard(small,"is discarded with the second-bind supernatant")
 graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),read_lengths=dict(READ_LENGTHS)))
 graph.annotate(Highlight("pcr_program","PCR cycle number is derived",
  "The application note (M-NA-00010 v2.0) states 10 ng input but not the cycle number; the guide's BLT PCR table assigns 8 cycles to 10–24 ng. The note's Table 2 lists a 350 bp insert for these libraries; the guide itself does not set an insert target."))
 graph.annotate(Highlight("library","Read configuration",
  "Application note: 2 × 251 bp on HiSeq 2500 Rapid Run and NovaSeq 6000 SP. The 10-cycle index reads follow the guide's instruction for 10 bp UD index codes; the note does not print index cycles.",slot=AnnotationSlot.SEQUENCING))
 graph.declare_reviewed_scope(["community_dna","blt_adapter","tn5","tag_program","tagment","tsb","stop_tagmentation","wash_beads","wash","i5_primer","i7_primer","polymerase","pcr_program","index_pcr","ipb","ipb_cleanup","library"],SOURCE)
 return graph.require_publishable()
