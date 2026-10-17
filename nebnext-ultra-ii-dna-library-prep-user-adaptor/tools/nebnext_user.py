"""NEBNext Ultra II DNA library prep with the E7600 USER-cleavable adaptor."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp
from dumbbell import Dumbbell, OpenedDumbbell
from endprep import dA_tailed_scene, repair_and_dA_tail
from adapter_ligation import reported_adapter_end
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, StepKind)
from topology import HairpinTopology

TITLE = "NEBNext Ultra II DNA Library Prep — USER-cleavable adaptor"
SOURCE = "Commercial authority: NEBNext Ultra II DNA Library Prep E7645/E7103 manual v7.0_9/22 and E7600 Dual Index Primers Set 1 manual v6.0_6/24."


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


INSERT = Construct([seg("DNA insert", "X" * 42, placeholder=True)], name="fragmented DNA")
END_PREP = repair_and_dA_tail(INSERT)
HAIRPIN = il.NEBNEXT_HAIRPIN_MODEL
HAIRPIN_TOPOLOGY = HairpinTopology(
    "NEBNext E7600 USER adaptor",
    HAIRPIN.sequence,
    first_stem_start=0,
    second_stem_start=(len(HAIRPIN.sequence)
                       - HAIRPIN.overhang_3_len - HAIRPIN.stem_len),
    stem_length=HAIRPIN.stem_len,
)
OPENED = il.NEBNEXT_USER_OPENED
PRE_USER = Dumbbell(INSERT, HAIRPIN.sequence, HAIRPIN.sequence,
                     name="hairpin-ligated NEBNext intermediate",
                     left_stem_length=HAIRPIN.stem_len,
                     right_stem_length=HAIRPIN.stem_len)
OPENED_LIBRARY = OpenedDumbbell(PRE_USER, OPENED, OPENED,
                                name="USER-opened adaptor library molecules")
I5_NAME, I7_NAME = "i501", "i701"
I5 = il.NEBNEXT_I5_SET1[I5_NAME]
I7_IN_OLIGO = il.NEBNEXT_I7_SET1[I7_NAME]
I7_READ = revcomp(I7_IN_OLIGO)
END_PREP_TIMES = ((20, 30), (65, 30))
LIGATION = (20, 15)
USER = (37, 15)
MINIMUM_PCR_CYCLES = 3
PCR_CYCLES = (3, 15)   # E7103/E7645 v7.0_9/22 Step 4.1.3 and Tables 4.1/4.2


def _cleanup(graph, reviewed, ident, label, material, beads_label, waste_text,
             kept_label, *, kind=MaterialKind.MOLECULAR):
    """A manual-specified SPRIselect / NEBNext Sample Purification Bead cleanup."""
    beads = graph.input(ident + "_beads", beads_label, MaterialKind.SELECTION)
    kept = graph.product(ident + "_kept", kept_label, kind, state=material.state)
    waste = graph.product(ident + "_waste", waste_text, MaterialKind.WASTE)
    graph.react(ident, label, StepKind.CLEANUP,
                ((material, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (kept, waste), changes_molecule=False)
    graph.discard(waste, "bead supernatant removed on the magnet")
    reviewed.extend((beads.id, ident, kept.id, waste.id))
    return kept


def hairpin_oligo() -> Construct:
    p = il.NEBNEXT_HAIRPIN_DU_POS
    n = HAIRPIN.stem_len
    end = len(HAIRPIN.sequence) - HAIRPIN.overhang_3_len
    return Construct([
        seg("5-prime stem", HAIRPIN.sequence[:n], "r2"),
        seg("Read-2-side loop", HAIRPIN.sequence[n:p], "r2"),
        seg("dU", HAIRPIN.sequence[p:p + 1], "w1"),
        seg("Read-1-side loop", HAIRPIN.sequence[p + 1:end - n], "r1"),
        seg("3-prime stem", HAIRPIN.sequence[end - n:end], "r1"),
        seg("3-prime dT", HAIRPIN.overhang_3, None),
    ], name="NEBNext Adaptor for Illumina")


def i5_primer() -> Construct:
    return Construct([seg("P5", il.P5, "p5"), seg("i501", I5, "cbc",
        feature=feature("sample_i5", "sample_index", "whitelist", whitelist="NEBNext index set 1")),
        seg("Read 1 arm", il.TRUSEQ_READ1, "r1")], name="NEBNext i501; terminal *T")


def i7_primer() -> Construct:
    return Construct([seg("P7", il.P7, "p7"), seg("i701", I7_IN_OLIGO, "cbc",
        feature=feature("sample_i7", "sample_index", "whitelist", whitelist="NEBNext index set 1")),
        seg("Read 2 arm", il.TRUSEQ_READ2, "r2")], name="NEBNext i701; terminal *T")


def final_library() -> Construct:
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5", I5, "cbc", feature=feature("sample_i5", "sample_index", "whitelist", whitelist="NEBNext index set 1")),
        seg("Read 1 arm", il.TRUSEQ_READ1, "r1"),
        seg("insert", "X" * 42, placeholder=True),
        seg("dA junction", "A"),
        seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        seg("i7 reverse complement", I7_READ, "cbc", feature=feature("sample_i7", "sample_index", "whitelist", whitelist="NEBNext index set 1")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="PCR-completed NEBNext library")
    problems = sp.verify(lib, sequencing_primers())
    if problems:
        raise ValueError("invalid NEBNext final library: " + "; ".join(problems))
    return lib


def sequencing_primers() -> tuple[sp.SeqPrimer, ...]:
    return (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])


def end_prep_scene():
    return dA_tailed_scene(END_PREP)


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("fragmented_dna","fragmented double-stranded DNA",MaterialKind.MOLECULAR,
                    state=INSERT,primary=True)
    reviewed=[]
    end_buffer=graph.input("endprep_buffer","NEBNext Ultra II End Prep Reaction Buffer",MaterialKind.CONDITION)
    repair_mix=graph.input("repair_mix","NEBNext Ultra II End Prep Enzyme Mix (composition not printed)",MaterialKind.CATALYST)
    repaired=graph.product("repaired_dna","blunt phosphorylated DNA",MaterialKind.MOLECULAR)
    graph.react("end_repair","End Prep: repair fragment ends; 5′ phosphorylation inferred (20 °C 30 min)",StepKind.END_REPAIR,
                ((dna,InputRole.TEMPLATE),(end_buffer,InputRole.CONDITION),(repair_mix,InputRole.CATALYST)),(repaired,))
    tailed=graph.product("tailed_dna","dA-tailed DNA",MaterialKind.MOLECULAR,
                         state=tuple(end_prep_scene().rows()))
    graph.react("da_tail","End Prep: add one 3′ dA to both insert ends (then 65 °C 30 min)",StepKind.TAILING,
                ((repaired,InputRole.TEMPLATE),(repair_mix,InputRole.CATALYST)),(tailed,))
    adaptor=graph.input("hairpin_adaptor","E7600 USER-cleavable hairpin adaptor",
                        MaterialKind.ADAPTER,state=hairpin_oligo())
    enhancer=graph.input("ligation_enhancer","NEBNext Ligation Enhancer",MaterialKind.CONDITION)
    ligase=graph.input("dna_ligase","NEBNext Ultra II Ligation Master Mix",MaterialKind.CATALYST)
    closed=graph.product("closed_dumbbell","hairpin-ligated closed intermediate",
                         MaterialKind.MOLECULAR,state=PRE_USER)
    graph.react("hairpin_ligation","Ligate one hairpin adaptor to each dA-tailed end (20 °C 15 min)",
                StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adaptor,InputRole.STRUCTURAL),
                 (enhancer,InputRole.CONDITION),(ligase,InputRole.CATALYST)),(closed,),
                geometry=reported_adapter_end(
                    "E7600 hairpin adaptor", kind="cohesive", polarity="3-prime",
                    terminal_sequence=HAIRPIN.overhang_3, body_sequence=HAIRPIN.sequence[:-1],
                    substrate_scenes=(end_prep_scene(),),
                    statement="E7600 v6.0_6/24 p.3 prints the full hairpin oligo with a terminal phosphorothioate-protected 3-prime T; the source-known fold and dA-tailed substrate determine the cohesive junction."))
    user=graph.input("user_enzyme","USER Enzyme (UDG plus Endo VIII), 3 µl added to the ligation mixture",
                     MaterialKind.CATALYST)
    opened=graph.product("opened_adaptor","USER-opened adaptor library molecules",
                         MaterialKind.MOLECULAR,state=OPENED_LIBRARY)
    cleavage_waste=graph.product("user_waste","excised uracil-derived adaptor residues",
                                 MaterialKind.WASTE)
    graph.react("user_open","Remove dU and cleave both hairpin loops (37 °C 15 min)",StepKind.DIGESTION,
                ((closed,InputRole.TEMPLATE),(user,InputRole.CATALYST)),
                (opened,cleavage_waste))
    graph.discard(cleavage_waste,"uracil-derived cleavage products do not enter the library")
    opened=_cleanup(graph,reviewed,"ligation_cleanup","Clean up the ligation without size selection (0.9X beads; input ≤ 50 ng, Section 3B)",
                    opened,"0.9X SPRIselect or NEBNext Sample Purification Beads",
                    "ligation supernatant (excess adaptor, reagents, small fragments)",
                    "purified USER-opened adaptor-ligated DNA")
    i5=graph.input("i5_primer","NEBNext i501 PCR primer",MaterialKind.OLIGO,state=i5_primer())
    i7=graph.input("i7_primer","NEBNext i701 PCR primer",MaterialKind.OLIGO,state=i7_primer())
    polymerase=graph.input("pcr_polymerase","NEBNext Ultra II Q5 Master Mix",MaterialKind.CATALYST)
    final=final_library()
    amplified=graph.product("amplified","PCR-amplified dual-index library",MaterialKind.MOLECULAR,state=final)
    graph.react("index_pcr","Complete the opened adaptor and install dual indexes (98 °C 30 s; 3–15 cycles 98 °C 10 s / 65 °C 75 s; 65 °C 5 min)",
                StepKind.PCR,
                ((opened,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(amplified,),
                geometry=PCRGeometry.unavailable(
                    final,reason="The opened hairpin is represented topologically rather than as a complete PCR duplex."))
    library=_cleanup(graph,reviewed,"pcr_cleanup","Clean up the PCR reaction (0.9X beads, elute 0.1X TE)",
                     amplified,"0.9X SPRIselect or NEBNext Sample Purification Beads",
                     "PCR supernatant (primers and small products)",final.name,
                     kind=MaterialKind.LIBRARY)
    graph.sequence(library,SequencingPlan(library=final,primers=sequencing_primers()))
    graph.annotate(Highlight("ligation_cleanup","Input-dependent cleanup branch",
        "E7103/E7645 v7.0_9/22 Section 3: inputs > 50 ng use double-sided size selection "
        "(about 0.4X beads, keep the supernatant, then 0.2X more beads, keep the beads) instead of "
        "this 0.9X cleanup. Either branch acts on the USER-treated ligation (96.5 µl); only the "
        "≤ 50 ng 0.9X branch is drawn."))
    graph.annotate(Highlight("index_pcr","Why at least three PCR cycles",
        "The manual states that NEBNext adaptors have a truncated design and require a minimum of "
        "3 amplification cycles to add the complete adaptor sequences (Tables 4.1/4.2 footnotes). "
        "Recommended cycles range from 3 (≥ 100 ng) to 15 (0.5 ng, target-enrichment yield)."))
    graph.declare_reviewed_scope(("fragmented_dna","repair_mix","end_repair",
      "da_tail","hairpin_adaptor","dna_ligase","hairpin_ligation",
      "user_enzyme","user_open","i5_primer","i7_primer","pcr_polymerase",
      "index_pcr","endprep_buffer","ligation_enhancer",*reviewed),SOURCE)
    return graph.require_publishable()
