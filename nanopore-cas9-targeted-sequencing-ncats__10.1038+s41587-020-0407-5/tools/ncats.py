"""Molecular model for nanopore Cas9-targeted sequencing (nCATS)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      annotation_rows, complement_segments)
from endprep import LigationJunction, dA_tailed_scene, repair_and_dA_tail
from selective_ligation import dephosphorylate_then_cut
from adapter_ligation import reported_adapter_end
from protocol_graph import (AnnotationSlot, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)

def seg(name, top, tag=None, **kw): return Segment(name, top, tag, **kw)

GENOME = Construct([
    seg("left genomic flank", "X"*22, placeholder=True),
    seg("left Cas9 target", "N"*20, "me", placeholder=True),
    seg("region of interest", "X"*44, placeholder=True, length_bp=20000),
    seg("right Cas9 target", "N"*20, "me", placeholder=True),
    seg("right genomic flank", "X"*22, placeholder=True),
], name="high-molecular-weight genomic DNA")

DEPHOSPHORYLATED = GENOME

CUT_TARGET = Construct([
    seg("left Cas9-cleavage flank", "N"*10, "r1", placeholder=True),
    seg("region of interest", "X"*44, placeholder=True, length_bp=20000),
    seg("right Cas9-cleavage flank", "N"*10, "r2", placeholder=True),
], name="Cas9-flanked target")

SELECTIVE_ENDS = dephosphorylate_then_cut(GENOME, CUT_TARGET)

# Taq adds dA only after Cas9 has exposed ligatable, phosphorylated ends.  The shared
# end-prep object derives opposed 3' A overhangs from the target rather than drawing them.
DA_TARGET = repair_and_dA_tail(CUT_TARGET)

TA_JUNCTION = LigationJunction("T", "A")
FINAL_LIBRARY = Construct([
    seg("left motor-loaded ONT adapter", "[LSK109 motor adapter]", "r1", placeholder=True),
    seg("left T:A junction", TA_JUNCTION.adapter_overhang,
        bottom=TA_JUNCTION.insert_overhang),
    *CUT_TARGET.segments,
    seg("right A:T junction", TA_JUNCTION.insert_overhang,
        bottom=TA_JUNCTION.adapter_overhang),
    seg("right ONT adapter", "[LSK109 adapter]", "r2", placeholder=True),
], name="nCATS target library")

def rows(con, label=""):
    unknown=tuple(s.name+"'" for s in con if s.is_role_token())
    return [*Scene.duplex(list(con),label=label,unpaired=unknown).rows(),*annotation_rows(con)]

def end_state_rows(con, label, mod5):
    top=list(con); bottom=complement_segments(top)
    sc=Scene(); sc.strand("top",top,label=label,mod5=mod5)
    sc.anneal("bottom",bottom,to="top",pair=(top[0].name+"'",top[0].name),
              label=label,mod5=mod5)
    return [*sc.rows(),*annotation_rows(con)]

def final_rows():
    unknown=tuple(s.name+"'" for s in FINAL_LIBRARY if s.is_role_token())
    sc=Scene.duplex(list(FINAL_LIBRARY),label="nCATS target",unpaired=unknown)
    sc.junction("top","left T:A junction","left Cas9-cleavage flank","adapter ligation")
    sc.junction("top","right Cas9-cleavage flank","right A:T junction","adapter ligation")
    return [*sc.rows(),*annotation_rows(FINAL_LIBRARY)]

INITIAL_NAME = "High-molecular-weight genomic DNA"
INITIAL_ROWS = tuple(end_state_rows(GENOME,"gDNA","p"))

def workflow():
    wf=Workflow(MolecularState(INITIAL_NAME,INITIAL_ROWS))
    wf.react("Quick CIP dephosphorylation",end_state_rows(DEPHOSPHORYLATED,"old ends","OH"),
             note="Pre-existing DNA ends lose the 5-prime phosphate required for ligation.")
    wf.react("Cas9–guide RNP cleavage around the target",end_state_rows(SELECTIVE_ENDS.selected_fragment,"fresh cuts","p"),
             note="Fresh cuts retain ligatable 5-prime phosphates; the old ends remain blocked.")
    wf.react("Taq plus dATP at 72 °C",dA_tailed_scene(DA_TARGET).rows(),
             note="A-tailing prepares the newly exposed target ends for ONT adapter ligation.")
    wf.react("Quick Ligase plus LSK109 adapters",final_rows(),
             note="Adapters preferentially ligate at Cas9-created, phosphorylated ends.")
    return wf

def protocol_graph():
    source="Defining source: https://doi.org/10.1038/s41587-020-0407-5"
    graph=ProtocolGraph("nCATS",authority=source)
    genome=graph.input("genome","high-molecular-weight genomic DNA",MaterialKind.BIOLOGICAL,
                       state=GENOME,primary=True)
    cip=graph.input("cip","Quick CIP phosphatase (NEB M0508) in CutSmart",MaterialKind.CATALYST)
    blocked=graph.product("dephosphorylated_genome","genomic DNA with blocked pre-existing ends",
                          MaterialKind.MOLECULAR,state=DEPHOSPHORYLATED)
    graph.react("dephosphorylate","Dephosphorylate all pre-existing DNA ends (37 °C 10 min; 80 °C 2 min inactivation)",
                StepKind.DEPHOSPHORYLATION,
                ((genome,InputRole.TEMPLATE),(cip,InputRole.CATALYST)),(blocked,))
    guides=graph.input("guides","pooled IDT crRNA:tracrRNA duplexes on complementary strands flanking the target (Supplementary Table 1)",
                       MaterialKind.OLIGO)
    cas9=graph.input("cas9","HiFi Cas9 Nuclease V3 (IDT) preassembled as RNP",MaterialKind.CATALYST)
    target=graph.product("cut_target","target fragment with fresh phosphorylated Cas9 ends",
                         MaterialKind.MOLECULAR,state=SELECTIVE_ENDS.selected_fragment)
    flanks=graph.product("blocked_flanks","off-target and old-ended genomic material",
                         MaterialKind.WASTE)
    graph.react("cas9_cut","Cut around the target with Cas9-guide RNPs",StepKind.DIGESTION,
                ((blocked,InputRole.TEMPLATE),(guides,InputRole.STRUCTURAL),
                 (cas9,InputRole.CATALYST)),(target,flanks))
    graph.discard(flanks,"old genomic ends remain dephosphorylated; at each cut Cas9 stays bound on the 5′ side of the gRNA, so ligation is preferential on the 3′ side (Results)")
    taq=graph.input("taq","Taq DNA polymerase (NEB M0267) with dATP, added with the RNP",MaterialKind.CATALYST)
    tailed=graph.product("a_tailed_target","dA-tailed Cas9 target fragment",MaterialKind.MOLECULAR,
                         state=dA_tailed_scene(DA_TARGET))
    graph.react("da_tail","Add one 3-prime dA to newly cut target ends (72 °C 5 min, same tube)",StepKind.TAILING,
                ((target,InputRole.TEMPLATE),(taq,InputRole.CATALYST)),(tailed,))
    adapter=graph.input("ont_adapter","LSK109 motor-loaded ONT adapter",MaterialKind.ADAPTER,
                        inferred=True,source=source+" (adapter body sequence proprietary)")
    ligase=graph.input("ligase","Quick Ligase (NEB M2200) in LSK109 ligation buffer",MaterialKind.CATALYST)
    ligated=graph.product("ligated",FINAL_LIBRARY.name+" in ligation reaction",MaterialKind.MOLECULAR,inferred=True,
                          source=source+" (terminal ONT adapter roles are proprietary)")
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,inferred=True,
                          source=source+" (terminal ONT adapter roles are proprietary)")
    geometry=reported_adapter_end(
        "LSK109 adapter",kind="cohesive",polarity="3-prime",terminal_sequence="T",
        substrate_scenes=(dA_tailed_scene(DA_TARGET),),
        statement="The nCATS workflow dA-tails fresh Cas9 ends before LSK109 adapter ligation, establishing a 3'-dT/3'-dA junction.")
    graph.react("adapter_ligation","Ligate LSK109 adapters preferentially to fresh Cas9 ends",
                StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(ligated,),geometry=geometry)
    beads=graph.input("ampure","0.3× AMPure XP with two LSK109 Long Fragment Buffer washes; Elution Buffer",
                      MaterialKind.SELECTION)
    bead_waste=graph.product("bead_supernatant","unbound adapters, short fragments and ligation reagents",
                             MaterialKind.WASTE)
    graph.react("bead_cleanup","Recover adapted DNA on 0.3× AMPure XP and elute",StepKind.CLEANUP,
                ((ligated,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(library,bead_waste))
    graph.discard(bead_waste,"is removed with the bead supernatant and LFB washes")
    graph.sequence(library,SequencingPlan(
        limitation="Nanopore sequencing is primer-free and the LSK109 motor-adapter bases are proprietary.",
        acquired_roles=("Sequencing output",),
        read_uses={"Sequencing output":"include"}))
    graph.declare_reviewed_scope((
        "genome","cip","dephosphorylate","guides","cas9","cas9_cut","taq",
        "da_tail","ont_adapter","ligase","adapter_ligation","ampure","bead_cleanup",
        "library"),source)
    graph.annotate(SourceDiscrepancy(
        target="cas9_cut",title="Cas9 cleavage time differs between article and detailed protocol",
        claims=(
            SourceClaim("Published article Online Methods, Cas9 Cleavage and Library Prep (PMC7145730)",
                        "Incubated at 37 °C for 20 min for Cas9 cleavage, then 72 °C 5 min."),
            SourceClaim("Cited protocols.io dx.doi.org/10.17504/protocols.io.68ihhue V.2 (10 Sep 2019), step 9",
                        "+37 °C for 15 min (Cas9 cleavage), then +72 °C 5 min.")),
        model_choice="The graph draws the cleavage and A-tailing order only and does not fix a cleavage time.",
        evidence="Both are the authors' own documents; the protocol revision predates the published article.",
        uncertainty="Neither source says which time was used for each reported experiment.",
        slot=AnnotationSlot.BEFORE_WORKFLOW))
    return graph.require_publishable()
