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
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from topology import HairpinTopology

TITLE = "NEBNext Ultra II DNA Library Prep — USER-cleavable adaptor"
SOURCE = "Commercial authority: NEBNext Ultra II DNA Library Prep E7645/E7103 and E7600 adaptor documentation."


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
    repair_mix=graph.input("repair_mix","NEBNext end-repair enzymes",MaterialKind.CATALYST)
    repaired=graph.product("repaired_dna","blunt phosphorylated DNA",MaterialKind.MOLECULAR)
    graph.react("end_repair","Repair and phosphorylate fragment ends",StepKind.END_REPAIR,
                ((dna,InputRole.TEMPLATE),(repair_mix,InputRole.CATALYST)),(repaired,))
    tailing=graph.input("tailing_polymerase","NEBNext dA-tailing polymerase",MaterialKind.CATALYST)
    tailed=graph.product("tailed_dna","dA-tailed DNA",MaterialKind.MOLECULAR,
                         state=tuple(end_prep_scene().rows()))
    graph.react("da_tail","Add one 3′ dA to both insert ends",StepKind.TAILING,
                ((repaired,InputRole.TEMPLATE),(tailing,InputRole.CATALYST)),(tailed,))
    adaptor=graph.input("hairpin_adaptor","E7600 USER-cleavable hairpin adaptor",
                        MaterialKind.ADAPTER,state=hairpin_oligo())
    ligase=graph.input("dna_ligase","NEBNext DNA ligase",MaterialKind.CATALYST)
    closed=graph.product("closed_dumbbell","hairpin-ligated closed intermediate",
                         MaterialKind.MOLECULAR,state=PRE_USER)
    graph.react("hairpin_ligation","Ligate one hairpin adaptor to each dA-tailed end",
                StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adaptor,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(closed,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The dumbbell model guarantees two dA/dT junctions; the compact reaction geometry does not flatten both hairpins into one Scene."))
    user=graph.input("user_enzyme","USER enzyme (UDG plus Endonuclease VIII)",
                     MaterialKind.CATALYST)
    opened=graph.product("opened_adaptor","USER-opened adaptor library molecules",
                         MaterialKind.MOLECULAR,state=OPENED_LIBRARY)
    cleavage_waste=graph.product("user_waste","excised uracil-derived adaptor residues",
                                 MaterialKind.WASTE)
    graph.react("user_open","Remove dU and cleave both hairpin loops",StepKind.DIGESTION,
                ((closed,InputRole.TEMPLATE),(user,InputRole.CATALYST)),
                (opened,cleavage_waste))
    graph.discard(cleavage_waste,"uracil-derived cleavage products do not enter the library")
    i5=graph.input("i5_primer","NEBNext i501 PCR primer",MaterialKind.OLIGO)
    i7=graph.input("i7_primer","NEBNext i701 PCR primer",MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","NEBNext library PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    graph.react("index_pcr","Complete the opened adaptor and install dual indexes",
                StepKind.PCR,
                ((opened,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(
                    final,reason="The opened hairpin is represented topologically rather than as a complete PCR duplex."))
    graph.sequence(library,SequencingPlan(library=final,primers=sequencing_primers()))
    graph.declare_reviewed_scope(("fragmented_dna","repair_mix","end_repair",
      "tailing_polymerase","da_tail","hairpin_adaptor","dna_ligase","hairpin_ligation",
      "user_enzyme","user_open","i5_primer","i7_primer","pcr_polymerase",
      "index_pcr","library"),SOURCE)
    return graph.require_publishable()
