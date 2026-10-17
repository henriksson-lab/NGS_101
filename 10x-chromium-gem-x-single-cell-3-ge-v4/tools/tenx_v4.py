"""Construct model for 10x GEM-X Single Cell 3' Gene Expression v4 (CG000731)."""
from __future__ import annotations
import illumina as il
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,feature,revcomp
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

PARTIAL_R1=il.TRUSEQ_READ1[11:]
CELL_BARCODE=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x-gem-x-3prime-v4")
UMI_FEATURE=feature("umi","umi","random")
I5_FEATURE=feature("sample_index_i5","sample_index","unknown")
I7_FEATURE=feature("sample_index_i7","sample_index","unknown")
FEATURES={"cell barcode":CELL_BARCODE,"UMI":UMI_FEATURE,"i5":I5_FEATURE,"i7 reverse complement":I7_FEATURE}
def seg(name,top,tag=None,**kw):
    kw.setdefault("feature",FEATURES.get(name)); return Segment(name=name,top=top,tag=tag,**kw)
def inferred(name,top,tag=None,**kw): return Segment(name=name,top=top,tag=tag,inferred=True,**kw)
def bead_oligo(): return [seg("Partial Read 1",PARTIAL_R1,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("poly(dT)30","T"*30),seg("V","V",placeholder=True),seg("N","N",placeholder=True)]
def inferred_tso(): return [inferred("TSO handle",CDNA_REVERSE+"TACAT","tso"),inferred("rGrGrG","GGG","tso")]
def capture_scene():
    mrna=[seg("mRNA","X"*32,placeholder=True),seg("poly(A)","A"*30)]; sc=Scene(); sc.strand("mRNA",mrna,label="mRNA")
    sc.anneal("bead primer",bead_oligo(),to="mRNA",pair=("poly(dT)30","poly(A)"),label="released GEM-X bead primer",unpaired=("V","N")); sc.arrow("bead primer","reverse transcriptase"); return sc
def first_strand(): return Construct([*bead_oligo(),seg("antisense cDNA","X"*34,placeholder=True),inferred("CCC","CCC","tso")],name="v4 first strand")
def switch_scene():
    first=first_strand(); sc=Scene(); sc.strand("first",list(first),label="first strand"); sc.anneal("TSO",inferred_tso(),to="first",pair=("rGrGrG","CCC"),label="TSO",above=True); sc.arrow("first","RT switches template"); return sc
def amplified_cdna(): return Construct([*bead_oligo(),seg("antisense cDNA","X"*34,placeholder=True),inferred("copied predecessor-style TSO",revcomp("".join(part.top for part in inferred_tso())),"tso")],name="v4 amplified cDNA")
SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["I2"],sp.TRUSEQ["R2"])
def switched_strand():
    """Single-stranded template-switched first strand (before cDNA PCR makes a duplex)."""
    sc=Scene(); sc.strand("switched",list(amplified_cdna()),label="template-switched first-strand cDNA"); return sc
def final_library(insert_nt=36):
    lib=Construct([seg("P5",il.P5,"p5"),seg("i5","J"*10,"cbc",placeholder=True,feature=I5_FEATURE),seg("Read 1 arm",il.TRUSEQ_READ1,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL_BARCODE),seg("UMI","U"*12,"umi",placeholder=True,feature=UMI_FEATURE),seg("poly(dT)30","T"*30),seg("V","V",placeholder=True),seg("N","N",placeholder=True),seg("cDNA insert","X"*insert_nt,placeholder=True),seg("dA junction","A"),seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),seg("i7 reverse complement","I"*10,"cbc",placeholder=True,feature=I7_FEATURE),seg("P7 reverse complement",il.P7_RC,"p7")],name="10x 3' v4 library")
    problems=sp.verify(lib,SEQ_PRIMERS)
    if problems: raise ValueError("invalid v4 library: "+"; ".join(problems))
    return lib
def read1_layout(): return [(x.cycles,x.feature.label) for x in sp.feature_spans(final_library(),SEQ_PRIMERS,{"Read 1":28})]

# Same catalog reagent in this guide and the predecessor oligo appendix.
# Assignment is catalog-linked, not a direct reprint of the current supplied oligos.
CDNA_FORWARD = PARTIAL_R1
CDNA_REVERSE = "AAGCAGTGGTATCAACGCAGAG"

def cdna_primers():
    return (Construct([seg("cDNA forward", CDNA_FORWARD, "r1")], name="PN-2000089 forward"),
            Construct([seg("cDNA reverse", CDNA_REVERSE, "tso")], name="PN-2000089 reverse"))

def protocol_graph():
    authority = "10x Genomics CG000731 Rev A; manual singleplex GEX workflow"
    graph = ProtocolGraph("10x Chromium GEM-X Single Cell 3-prime Gene Expression v4", authority=authority)
    cells = graph.input("cells", "cells containing poly(A) RNA", MaterialKind.BIOLOGICAL,
                        state=molecule("cells", (role("poly(A) RNA"),)), primary=True)
    gems = graph.input("gems", "GEM-X partitioning condition", MaterialKind.CONDITION)
    occupied = graph.product("occupied", "occupied GEMs", MaterialKind.MOLECULAR, state=cells.state)
    failed = graph.product("failed", "empty, multiplet or unresolved GEM products", MaterialKind.WASTE)
    graph.react("partition", "Partition cells with gel beads", StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (gems, InputRole.CONDITION)), (occupied, failed), changes_molecule=False)
    graph.discard(failed, "does not retain a resolvable single-cell barcode")
    rt_primer = graph.input("rt_primer", "cell-barcode/UMI oligo(dT) gel-bead primer", MaterialKind.OLIGO,
                           state=Construct(bead_oligo(), name="RT primer"),
                           inferred=False, source=authority)
    enzyme = graph.input("rt", "RT Enzyme E", MaterialKind.CATALYST)
    first = graph.product("first_strand", "first-strand cDNA", MaterialKind.MOLECULAR,
                          state=first_strand(), inferred=True, source=authority + "; conditional intermediate reconstruction")
    graph.react("reverse_transcribe", "INFERRED — capture and reverse-transcribe poly(A) RNA", StepKind.REVERSE_TRANSCRIPTION,
                ((occupied, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (enzyme, InputRole.CATALYST)), (first,))
    tso = graph.input("tso", "Template Switch Oligo B PN-2001027", MaterialKind.OLIGO,
                      state=Construct(inferred_tso(), name="template-switch oligo"),
                      inferred=True, source=authority)
    switched = graph.product("switched", "template-switched cDNA", MaterialKind.MOLECULAR,
                             state=switched_strand(), inferred=True, source=authority + "; conditional end reconstruction")
    graph.react("template_switch", "INFERRED — copy the template-switch oligo", StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL), (enzyme, InputRole.CATALYST)), (switched,))
    recovery = graph.input("gem_recovery", "Recovery Agent (GEM break)", MaterialKind.CONDITION)
    broken = graph.product("broken_gems", "aqueous post-GEM-RT phase", MaterialKind.MOLECULAR, state=switched.state, inferred=True)
    oil = graph.product("recovery_oil", "Recovery Agent/Partitioning Oil phase", MaterialKind.WASTE)
    graph.react("recover", "Break GEMs with Recovery Agent; remove the oil phase (§2.1a–b)", StepKind.EXTRACTION,
                ((switched, InputRole.TEMPLATE), (recovery, InputRole.CONDITION)), (broken, oil), changes_molecule=False)
    graph.discard(oil, "pink Recovery Agent/Partitioning Oil phase aspirated from the tube bottom")
    recovered = cleanup(graph, broken, "dynabeads_cleanup", "Dynabeads MyOne SILANE cleanup (§2.1)",
                        "recovered", "bead-purified post-GEM-RT cDNA",
                        "leftover GEM-RT reagents and primers in the bead supernatant")
    forward, reverse = cdna_primers()
    pf = graph.input("preamp_forward", "cDNA Primers PN-2000089 forward (catalog-linked)", MaterialKind.OLIGO,
                     state=forward, source=authority + "; predecessor appendix with same catalog number")
    pr = graph.input("preamp_reverse", "cDNA Primers PN-2000089 reverse (catalog-linked)", MaterialKind.OLIGO,
                     state=reverse, source=authority + "; predecessor appendix with same catalog number")
    pol = graph.input("preamp_polymerase", "Amp Mix PN-2000047/2000103 (not Library Amp Mix)", MaterialKind.CATALYST)
    amplified = graph.product("amplified", "amplified cDNA", MaterialKind.MOLECULAR, state=amplified_cdna(), inferred=True)
    graph.react("preamp", "INFERRED — amplify recovered cDNA (98/63/72 °C, 11–13 total cycles by targeted cell recovery; §2.2)", StepKind.PCR,
                ((recovered, InputRole.TEMPLATE), (pf, InputRole.STRUCTURAL), (pr, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(amplified.state,
                    reason="PN-2000089 primer sequences are catalog-linked; complete current-kit cDNA ends are not established in this model."))
    cleaned = cleanup(graph, amplified, "cdna_cleanup", "SPRIselect 0.6X cDNA cleanup (§2.3)",
                      "cleaned_cdna", "SPRIselect-purified amplified cDNA",
                      "unbound short products and PCR reagents in the 0.6X supernatant")
    fragmentation = graph.input("fragmentase", "Fragmentation Mix: Fragmentation Enzyme PN-2000090 + Buffer PN-2000091", MaterialKind.CATALYST)
    fragments = graph.product("fragments", "fragmented cDNA", MaterialKind.MOLECULAR)
    graph.react("fragment", "Fragment 10 µl (25%) of purified cDNA at 32 °C, 5 min (§3.1)", StepKind.FRAGMENTATION,
                ((cleaned, InputRole.TEMPLATE), (fragmentation, InputRole.CATALYST)), (fragments,))
    end_mix = graph.input("end_mix", "same Fragmentation Mix, no reagent added", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "end-repaired A-tailed cDNA", MaterialKind.MOLECULAR)
    graph.react("end_prep", "End-repair and A-tail in the same tube at 65 °C, 30 min (§3.1)", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (end_mix, InputRole.CATALYST)), (tailed,))
    sized = double_sided(graph, tailed, "fragment_size_selection",
                         "Double-sided SPRIselect 0.6X/0.8X (§3.2)", "sized_fragments",
                         "size-selected A-tailed fragments")
    adapter = graph.input("read2_adapter", "Ligation Mix PN-2001109", MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "DNA Ligase PN-220110/220131", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "adapter-ligated cDNA", MaterialKind.MOLECULAR)
    graph.react("ligate", "Ligate adapters at 20 °C, 15 min (§3.3)", StepKind.LIGATION,
                ((sized, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The current Ligation Mix PN-2001109 oligos and physical ligating ends are not established in this declaration; the final library alone does not establish reagent geometry."))
    ligated = cleanup(graph, adapted, "ligation_cleanup", "SPRIselect 0.8X post-ligation cleanup (§3.4)",
                      "ligated", "purified adapter-ligated cDNA",
                      "unligated adapter and ligation reagents in the 0.8X supernatant")
    index_f = graph.input("index_forward", "Dual Index Kit TT Set A forward oligo", MaterialKind.OLIGO)
    index_r = graph.input("index_reverse", "Dual Index Kit TT Set A reverse oligo", MaterialKind.OLIGO)
    index_pol = graph.input("index_polymerase", "Library Amp Mix PN-2000531 or Amp Mix PN-2000047/2000103", MaterialKind.CATALYST)
    final = final_library()
    indexed = graph.product("indexed", "sample-index PCR product", MaterialKind.LIBRARY, state=final)
    graph.react("index_pcr", "Select the barcode-bearing end and add sample indexes (98/54/72 °C, 5–16 total cycles by cDNA input; §3.5)", StepKind.PCR,
                ((ligated, InputRole.TEMPLATE), (index_f, InputRole.STRUCTURAL), (index_r, InputRole.STRUCTURAL),
                 (index_pol, InputRole.CATALYST)), (indexed,), geometry=PCRGeometry.unavailable(final,
                    reason="The TT table supplies index pairs, not full physical PCR oligos; exact primer geometry is not declared."))
    library = double_sided(graph, indexed, "library_size_selection",
                           "Post-sample-index double-sided SPRIselect 0.6X/0.8X (§3.6)", "library",
                           "dual-index GEX library", kind=MaterialKind.LIBRARY, state=final)
    graph.sequence(library, SequencingPlan(library=final, primers=SEQ_PRIMERS,
                   acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
                   read_lengths={"Read 1":28,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":90}))
    graph.declare_reviewed_scope(["cells", "gems", "partition", "rt_primer", "rt", "reverse_transcribe",
        "tso", "template_switch", "gem_recovery", "recover", "dynabeads_cleanup", "dynabeads_cleanup_beads", "preamp_forward", "preamp_reverse",
        "preamp_polymerase", "preamp", "cdna_cleanup", "cdna_cleanup_beads", "fragmentase", "fragment", "end_mix", "end_prep",
        "fragment_size_selection", "fragment_size_selection_beads", "read2_adapter", "ligase", "ligate", "ligation_cleanup", "ligation_cleanup_beads",
        "index_forward", "index_reverse", "index_polymerase", "index_pcr", "library_size_selection", "library_size_selection_beads", "library"], authority)
    return graph.require_publishable()


def cleanup(graph, substrate, step_id, label, kept_id, kept_label, waste_label):
    """One bead cleanup: retained bead eluate plus the discarded supernatant."""
    beads = graph.input(step_id + "_beads", label, MaterialKind.SELECTION)
    kept = graph.product(kept_id, kept_label, MaterialKind.MOLECULAR, state=substrate.state,
                         inferred=substrate.inferred, source=substrate.source)
    waste = graph.product(step_id + "_waste", waste_label, MaterialKind.WASTE)
    graph.react(step_id, label, StepKind.CLEANUP,
                ((substrate, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (kept, waste))
    graph.discard(waste, "removed with the supernatant and ethanol washes")
    return kept


def double_sided(graph, substrate, step_id, label, kept_id, kept_label, *,
                 kind=MaterialKind.MOLECULAR, state=None):
    """Source double-sided SPRIselect: 0.6X binds large species, 0.8X (cumulative) keeps the middle."""
    beads = graph.input(step_id + "_beads", label, MaterialKind.SELECTION)
    kept = graph.product(kept_id, kept_label, kind, state=substrate.state if state is None else state)
    large = graph.product(step_id + "_large", "large fragments bound at 0.6X", MaterialKind.WASTE)
    small = graph.product(step_id + "_small", "small fragments left in the 0.8X supernatant", MaterialKind.WASTE)
    graph.react(step_id, label, StepKind.SELECTION,
                ((substrate, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (kept, large, small))
    graph.discard(large, "first-bead pellet is discarded; the supernatant is carried forward")
    graph.discard(small, "second-bead supernatant is removed; the bead-bound middle fraction is eluted")
    return kept
