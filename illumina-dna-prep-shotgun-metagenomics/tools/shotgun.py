"""Illumina DNA Prep shotgun-metagenomics model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import nextera as nx
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, StepKind, Terminal, Disposition)
from protocol_model import molecule, role

# Guide v12 defines 10-base UD codes. Official adapter table specifies these arms;
# the application note selects Set A, but leaves the actual well unassigned.
I5 = feature("dna_prep_i5", "sample_index", "unknown", note="Set A well unassigned")
I7 = feature("dna_prep_i7", "sample_index", "unknown", note="Set A well unassigned")

def index_primer(side):
 if side == "i5":
  return Construct([Segment("P5",il.P5,"p5"),role("i5",length_bp=10,molecular_feature=I5),Segment("i5 annealing arm",nx.S5,"s5")],name="Set A i5 primer")
 return Construct([Segment("P7",il.P7,"p7"),role("i7",length_bp=10,molecular_feature=I7),Segment("i7 annealing arm",nx.S7,"s7")],name="Set A i7 primer")

def shotgun_library():
 lib=Construct([Segment("P5",il.P5,"p5"),role("i5",length_bp=10,molecular_feature=I5),
  Segment("S5",nx.S5,"s5"),Segment("left mosaic end",nx.ME,"me"),role("metagenomic DNA fragment; variable length"),
  Segment("right mosaic end reverse complement",nx.ME_RC,"me"),Segment("S7 reverse complement",nx.S7_RC,"s7"),
  role("i7 reverse complement",length_bp=10,molecular_feature=I7),Segment("P7 reverse complement",il.P7_RC,"p7")],name="Illumina DNA Prep shotgun library")
 primers=tuple(sp.NEXTERA[x] for x in ("R1","I1","I2","R2"))
 problems=sp.verify(lib,primers,required_roles=tuple(x.role for x in primers))
 if problems:raise ValueError("; ".join(problems))
 return lib,primers


TITLE="Illumina DNA Prep shotgun metagenomics"
NOTES="01_shotgun-metagenomics.html"
SOURCE='Illumina <a href="https://support.illumina.com/downloads/illumina-dna-prep-reference-guide-1000000025416.html">DNA Prep Reference Guide 1000000025416 v12</a> and <a href="https://www.illumina.com/content/dam/illumina/gcs/assembled-assets/marketing-literature/novaseq-shotgun-metagenomics-app-note-m-na-00010/novaseq-shotgun-metagenomics-app-note-m-na-00010.pdf">shotgun-metagenomics application note M-NA-00010 v2.0</a>.'
SUMMARY="Selected reference implementation: DNA Prep guide1000000025416 v12 (April2026) at10ng community DNA, with on-bead tagmentation, wash, UD-index PCR and IPB recovery. The 2024 application note supplies the metagenomics context and Set A example; its historical executed recipe is not equated with the later guide."
CAVEAT="Public transposase adapter and index-primer arms are source-known. Native BLT strand gaps, bead attachment and pre-PCR end occupancy remain unassigned; the completed final library is distinct from that precursor. The selected Set A well remains unassigned. The guide is April2026; the ©2024 application note reports2×251 sequencing but leaves index cycles unassigned. The drawn guide-recommended2×151/10/10 run is a separate reference configuration."
FINAL_LIBRARY,SEQ_PRIMERS=shotgun_library()
FINAL_CAPTION="A representative community-DNA fragment after bead-linked tagmentation and dual-index PCR."
SEQUENCING_INTRO="Selected v12 guide recommendation:2×151 biological reads plus10 cycles per index. The ©2024 metagenomics application note reports2×251 on HiSeq2500 Rapid Run/NovaSeq6000 SP, with index acquisition unassigned in that historical run."

# Application note M-NA-00010 v2.0 lines 215-221: 10 ng community DNA in 30 µl, Illumina
# DNA Prep (M) Tagmentation, IDT for Illumina DNA/RNA UD Indexes Set A; lines 239-242:
# 2 x 251 bp on HiSeq 2500 Rapid Run and NovaSeq 6000 SP.  Guide p.26: UD index codes
# are 10 bp and the run must be configured for 10 bp index reads.  Index read cycles are
# taken from the guide, not printed in the note.
INPUT_NG = 10
# Guide 1000000025416 v12 p.27 BLT PCR table row "10-24 ng -> 8 cycles"; the note does not
# print its cycle number, so this is a derived assignment.
PCR_CYCLES_10_TO_24_NG = 8
READ_LENGTHS = {"Read 1": 151, "Index 1 (i7)": 10, "Index 2 (i5)": 10, "Read 2": 151}


def blt_adapter_inventory():
 sc=Scene()
 for name,arm in (("S5",nx.S5),("S7",nx.S7)):
  sc.strand(name,[Segment(name+" entry arm",arm),Segment(name+" mosaic end",nx.ME)],label="public transposase adapter family; native bead linkage/partner unassigned")
 return sc

def _tagged_state():
 """Known supplied arms retained without assuming a gap-filled native duplex."""
 sc=Scene()
 for name,arm in (("S5",nx.S5),("S7",nx.S7)):
  sc.strand(name,[Segment(name+" entry arm",arm),Segment(name+" mosaic end",nx.ME),role(name+"-tagged genomic strand; variable extent")],label="opposing tagged-strand class; native gaps/pairing and bead occupancy unassigned")
 return sc

def oligos():
 from chemdraw import oligo
 return tuple(oligo(x.name,list(x)) for x in (index_primer("i5"),index_primer("i7")))


def protocol_graph():
 graph=ProtocolGraph(TITLE,authority=SOURCE)
 dna=graph.input("community_dna","extracted total community DNA (10 ng in 30 µl)",MaterialKind.BIOLOGICAL,state=molecule("community DNA",(role("representative extracted community DNA molecule"),)),primary=True)
 adapter=graph.input("blt_adapter","bead-linked Nextera adapter transposomes (BLT)",MaterialKind.ADAPTER,state=blt_adapter_inventory(),source=SOURCE+" Public endpoint architecture is used; BLT reagent composition is proprietary.")
 tn5=graph.input("tn5","bead-linked Tn5 transposase",MaterialKind.CATALYST)
 tag=graph.input("tag_program","TB1; TAG program 55°C 15 min",MaterialKind.CONDITION)
 tagged=graph.product("tagged","bead-bound tagged community-DNA fragments",MaterialKind.MOLECULAR,state=_tagged_state())
 waste=graph.product("tagmentation_waste","non-library or same-end tagmentation products",MaterialKind.MOLECULAR)
 graph.react("tagment","Tagment community DNA with bead-linked transposomes",StepKind.TAGMENTATION,((dna,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST),(tag,InputRole.CONDITION)),(tagged,waste))
 graph.mark_terminal(Terminal(waste.id,Disposition.UNRESOLVED,"Failure to amplify does not establish physical disposal; exact partition through later bead washes is unassigned."))
 tsb=graph.input("tsb","Tagment Stop Buffer (TSB); PTC program 37°C 15 min",MaterialKind.CONDITION)
 stopped=graph.product("stopped","stopped bead-bound tagged fragments",MaterialKind.MOLECULAR,state=tagged.state)
 graph.react("stop_tagmentation","Stop tagmentation with TSB",StepKind.QUENCHING,((tagged,InputRole.TEMPLATE),(tsb,InputRole.CONDITION)),(stopped,),changes_molecule=False)
 wash=graph.input("wash_beads","magnet: discard supernatant; TWB 100 µl × 3, last TWB kept until PCR",MaterialKind.CONDITION)
 retained=graph.product("retained","washed bead-bound tagged fragments",MaterialKind.MOLECULAR,state=tagged.state)
 washed=graph.product("wash_waste","supernatant and TWB: unbound DNA and reaction components",MaterialKind.WASTE)
 graph.react("wash","Wash bead-bound tagged fragments",StepKind.CLEANUP,((stopped,InputRole.TEMPLATE),(wash,InputRole.CONDITION)),(retained,washed),changes_molecule=False)
 graph.discard(washed,"is removed with the post-tagmentation supernatant and TWB washes")
 i5=graph.input("i5_primer","i5 index adapter (IDT for Illumina DNA/RNA UD Indexes Set A well)",MaterialKind.OLIGO,state=index_primer("i5"))
 i7=graph.input("i7_primer","i7 index adapter (same pre-paired Set A well)",MaterialKind.OLIGO,state=index_primer("i7"))
 polymerase=graph.input("polymerase","Enhanced PCR Mix (EPM)",MaterialKind.CATALYST)
 program=graph.input("pcr_program",f"BLT PCR: 68°C 3 min, 98°C 3 min, {PCR_CYCLES_10_TO_24_NG} × (98°C 45 s, 62°C 30 s, 68°C 2 min), 68°C 1 min (cycle number derived from 10 ng input)",MaterialKind.CONDITION)
 amplified=graph.product("amplified","amplified dual-index library (unpurified)",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
 graph.react("index_pcr","Complete adapters and index the sample",StepKind.PCR,((retained,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),(i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST),(program,InputRole.CONDITION)),(amplified,),geometry=PCRGeometry.unavailable(retained.state,reason="Actual preceding bead-bound tagged fragments; native strand gaps and occupancy are unassigned, even though the public adapter and PCR-primer families are known."))
 ipb=graph.input("ipb","IPB: 45µl PCR supernatant +40µl water +45µl beads; transfer125µl supernatant to15µl IPB; two80%EtOH washes;32µl RSB, recover30µl",MaterialKind.SELECTION)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 large=graph.product("large_fraction","first-bead fraction (larger fragments)",MaterialKind.WASTE)
 small=graph.product("small_fraction","second supernatant (short products and reaction components)",MaterialKind.WASTE)
 graph.react("ipb_cleanup","Double-sided IPB cleanup, ethanol washes, elute in RSB",StepKind.CLEANUP,((amplified,InputRole.TEMPLATE),(ipb,InputRole.SELECTION)),(library,large,small),changes_molecule=False)
 graph.discard(large,"is discarded with the first bead plate")
 graph.discard(small,"is discarded with the second-bind supernatant")
 graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),read_lengths=dict(READ_LENGTHS)))
 graph.annotate(Highlight("pcr_program","PCR cycle number is derived",
  "The ©2024 application note states10ng and a350bp insert summary. The selected April2026 guide assigns8 PCR cycles at10–24ng; this is a current reference choice, not proof of the historical experiment's executed cycle count or insert extent."))
 graph.annotate(Highlight("library","Read configuration",
  "Separate source scopes: selected April2026 guide p32 recommends151×10×10×151; ©2024 application note reports2×251 on HiSeq2500 Rapid Run/NovaSeq6000 SP, with index cycles unassigned. No complete revision equivalence or historical151/10/10 acquisition is asserted.",slot=AnnotationSlot.SEQUENCING))
 graph.declare_reviewed_scope(("community_dna","blt_adapter","tn5","tag_program","tagment","tsb","stop_tagmentation","wash_beads","wash","i5_primer","i7_primer","polymerase","pcr_program","index_pcr","ipb","ipb_cleanup","library"),SOURCE)
 return graph.require_publishable()
