"""Construct model for 10x GEM-X Single Cell 5' Gene Expression v3 (CG000733)."""
from __future__ import annotations
from dataclasses import replace
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

PARTIAL_R1=il.TRUSEQ_READ1[11:]; SWITCH_SPACER="TTTCTTATAT"
CELL_BARCODE=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x-gem-x-5prime-v3")
UMI_FEATURE=feature("umi","umi","random")
I5_FEATURE=feature("sample_index_i5","sample_index","unknown")
I7_FEATURE=feature("sample_index_i7","sample_index","unknown")
FEATURES={"cell barcode":CELL_BARCODE,"UMI":UMI_FEATURE,"i5":I5_FEATURE,"i7 reverse complement":I7_FEATURE}
def seg(name,top,tag=None,**kw):
    kw.setdefault("feature",FEATURES.get(name)); return Segment(name=name,top=top,tag=tag,**kw)
def inferred(name,top,tag=None,**kw): return Segment(name=name,top=top,tag=tag,inferred=True,**kw)
def bead_tso(): return [seg("Partial Read 1",PARTIAL_R1,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("switch spacer",SWITCH_SPACER,"tso"),seg("rGrGrG","GGG","tso")]
def inferred_rt_primer(): return [inferred("RT handle",CDNA_REVERSE+"TAC","tso"),inferred("poly(dT)30","T"*30),inferred("V","V",placeholder=True),inferred("N","N",placeholder=True)]
def capture_scene():
    mrna=[seg("mRNA","X"*34,placeholder=True),seg("poly(A)","A"*30)]; sc=Scene(); sc.strand("mRNA",mrna,label="mRNA")
    sc.anneal("RT primer",inferred_rt_primer(),to="mRNA",pair=("poly(dT)30","poly(A)"),label="Poly-dT RT Primer B",unpaired=("V","N")); sc.arrow("RT primer","reverse transcriptase"); return sc
def first_strand(): return Construct([*inferred_rt_primer(),seg("antisense cDNA","X"*36,placeholder=True),inferred("CCC","CCC","tso")],name="v3 first strand")
def switch_scene():
    first=first_strand(); sc=Scene(); sc.strand("first",list(first),label="first strand"); sc.anneal("bead TSO",bead_tso(),to="first",pair=("rGrGrG","CCC"),label="barcoded bead TSO",above=True,mod5="bead"); sc.arrow("first","RT switches template"); return sc
def amplified_cdna(): return Construct([*bead_tso(),seg("5-prime cDNA","X"*40,placeholder=True),*[replace(part, name=part.name + " reverse complement", top=part.top if part.placeholder else revcomp(part.top)) for part in reversed(inferred_rt_primer())]],name="v3 amplified cDNA")
SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["I2"],sp.TRUSEQ["R2"])
def final_library(insert_nt=38):
    lib=Construct([seg("P5",il.P5,"p5"),seg("i5","J"*10,"cbc",placeholder=True,feature=I5_FEATURE),seg("Read 1 arm",il.TRUSEQ_READ1,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL_BARCODE),seg("UMI","U"*12,"umi",placeholder=True,feature=UMI_FEATURE),seg("switch spacer",SWITCH_SPACER,"tso"),seg("GGG","GGG","tso"),seg("5-prime cDNA insert","X"*insert_nt,placeholder=True),seg("dA junction","A"),seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),seg("i7 reverse complement","I"*10,"cbc",placeholder=True,feature=I7_FEATURE),seg("P7 reverse complement",il.P7_RC,"p7")],name="10x 5' GEM-X v3 library")
    problems=sp.verify(lib,SEQ_PRIMERS)
    if problems: raise ValueError("invalid 5' v3 library: "+"; ".join(problems))
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
    authority = "10x Genomics CG000733 Rev A; manual singleplex GEX workflow"
    graph = ProtocolGraph("10x Chromium GEM-X Single Cell 5-prime Gene Expression v3", authority=authority)
    cells = graph.input("cells", "cells containing poly(A) RNA", MaterialKind.BIOLOGICAL,
                        state=molecule("cells", (role("poly(A) RNA"),)), primary=True)
    gems = graph.input("gems", "GEM-X partitioning condition", MaterialKind.CONDITION)
    occupied = graph.product("occupied", "occupied GEMs", MaterialKind.MOLECULAR, state=cells.state)
    failed = graph.product("failed", "empty, multiplet or unresolved GEM products", MaterialKind.WASTE)
    graph.react("partition", "Partition cells with gel beads", StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (gems, InputRole.CONDITION)), (occupied, failed), changes_molecule=False)
    graph.discard(failed, "does not retain a resolvable single-cell barcode")
    rt_primer = graph.input("rt_primer", "Poly-dT RT Primer B PN-2001110", MaterialKind.OLIGO,
                           state=Construct(inferred_rt_primer(), name="RT primer"),
                           inferred=True, source=authority)
    enzyme = graph.input("rt", "RT Enzyme E", MaterialKind.CATALYST)
    first = graph.product("first_strand", "first-strand cDNA", MaterialKind.MOLECULAR,
                          state=first_strand(), inferred=True, source=authority + "; conditional intermediate reconstruction")
    graph.react("reverse_transcribe", "INFERRED — capture and reverse-transcribe poly(A) RNA", StepKind.REVERSE_TRANSCRIPTION,
                ((occupied, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (enzyme, InputRole.CATALYST)), (first,))
    tso = graph.input("tso", "cell-barcode/UMI gel-bead template-switch oligo", MaterialKind.OLIGO,
                      state=Construct(bead_tso(), name="template-switch oligo"),
                      inferred=False, source=authority)
    switched = graph.product("switched", "template-switched cDNA", MaterialKind.MOLECULAR,
                             state=amplified_cdna(), inferred=True, source=authority + "; conditional end reconstruction")
    graph.react("template_switch", "INFERRED — copy the template-switch oligo", StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL), (enzyme, InputRole.CATALYST)), (switched,))
    recovery = graph.input("gem_recovery", "GEM break and Dynabeads cleanup", MaterialKind.CONDITION)
    recovered = graph.product("recovered", "recovered cDNA", MaterialKind.MOLECULAR, state=switched.state, inferred=True)
    graph.react("recover", "Break GEMs and recover cDNA", StepKind.EXTRACTION,
                ((switched, InputRole.TEMPLATE), (recovery, InputRole.CONDITION)), (recovered,), changes_molecule=False)
    forward, reverse = cdna_primers()
    pf = graph.input("preamp_forward", "cDNA Primers PN-2000089 forward (catalog-linked)", MaterialKind.OLIGO,
                     state=forward, source=authority + "; predecessor appendix with same catalog number")
    pr = graph.input("preamp_reverse", "cDNA Primers PN-2000089 reverse (catalog-linked)", MaterialKind.OLIGO,
                     state=reverse, source=authority + "; predecessor appendix with same catalog number")
    pol = graph.input("preamp_polymerase", "Amp Mix", MaterialKind.CATALYST)
    amplified = graph.product("amplified", "amplified cDNA", MaterialKind.MOLECULAR, state=amplified_cdna(), inferred=True)
    graph.react("preamp", "INFERRED — amplify recovered cDNA", StepKind.PCR,
                ((recovered, InputRole.TEMPLATE), (pf, InputRole.STRUCTURAL), (pr, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(amplified.state,
                    reason="PN-2000089 primer sequences are catalog-linked; complete current-kit cDNA ends are not established in this model."))
    fragmentation = graph.input("fragmentase", "Fragmentation Enzyme", MaterialKind.CATALYST)
    fragments = graph.product("fragments", "fragmented cDNA", MaterialKind.MOLECULAR)
    graph.react("fragment", "Fragment amplified cDNA", StepKind.FRAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (fragmentation, InputRole.CATALYST)), (fragments,))
    end_mix = graph.input("end_mix", "Fragmentation Enzyme mix during end repair/A-tailing", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "end-repaired A-tailed cDNA", MaterialKind.MOLECULAR)
    graph.react("end_prep", "End-repair and A-tail cDNA fragments", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (end_mix, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("read2_adapter", "Ligation Mix PN-2001109", MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "DNA Ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "adapter-ligated cDNA", MaterialKind.MOLECULAR)
    graph.react("ligate", "Ligate adapters", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The current Ligation Mix PN-2001109 oligos and physical ligating ends are not established in this declaration; the final library alone does not establish reagent geometry."))
    index_f = graph.input("index_forward", "Dual Index Kit TT Set A forward oligo", MaterialKind.OLIGO)
    index_r = graph.input("index_reverse", "Dual Index Kit TT Set A reverse oligo", MaterialKind.OLIGO)
    index_pol = graph.input("index_polymerase", "Library Amp Mix", MaterialKind.CATALYST)
    final = final_library()
    library = graph.product("library", "dual-index GEX library", MaterialKind.LIBRARY, state=final)
    graph.react("index_pcr", "Select the barcode-bearing end and add sample indexes", StepKind.PCR,
                ((adapted, InputRole.TEMPLATE), (index_f, InputRole.STRUCTURAL), (index_r, InputRole.STRUCTURAL),
                 (index_pol, InputRole.CATALYST)), (library,), geometry=PCRGeometry.unavailable(final,
                    reason="The TT table supplies index pairs, not full physical PCR oligos; exact primer geometry is not declared."))
    graph.sequence(library, SequencingPlan(library=final, primers=SEQ_PRIMERS,
                   acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
                   read_lengths={"Read 1":28,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":90}))
    graph.declare_reviewed_scope(["cells", "gems", "partition", "rt_primer", "rt", "reverse_transcribe",
        "tso", "template_switch", "gem_recovery", "recover", "preamp_forward", "preamp_reverse",
        "preamp_polymerase", "preamp", "fragmentase", "fragment", "end_mix", "end_prep",
        "read2_adapter", "ligase", "ligate", "index_forward", "index_reverse", "index_polymerase", "index_pcr", "library"], authority)
    return graph.require_publishable()
