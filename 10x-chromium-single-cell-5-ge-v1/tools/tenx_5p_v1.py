"""Construct model for 10x Chromium Single Cell 5' Gene Expression v1 (CG000109)."""
from __future__ import annotations
import illumina as il
import rt
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,complement_segments,feature,revcomp
from protocol_graph import (AnnotationSlot, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import molecule, role
from pcr import PCRGeometry, PCRPrimer
from adapter_ligation import AdapterLigationGeometry

PARTIAL_R1=il.TRUSEQ_READ1[11:]
SWITCH_SPACER="TTTCTTATAT"
CELL_BARCODE=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x-chromium-5prime-v1")
UMI_FEATURE=feature("umi","umi","random")
SAMPLE_INDEX=feature("sample_index_i7","sample_index","unknown")
FEATURES={"cell barcode":CELL_BARCODE,"UMI":UMI_FEATURE,"i7 reverse complement":SAMPLE_INDEX}
def seg(name,top,tag=None,**kw):
    kw.setdefault("feature",FEATURES.get(name)); return Segment(name=name,top=top,tag=tag,**kw)
def bead_tso(): return [seg("Partial Read 1",PARTIAL_R1,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*10,"umi",placeholder=True),seg("switch spacer",SWITCH_SPACER,"tso"),seg("rGrGrG","GGG","tso")]
def rt_primer(): return [seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("AC","AC"),seg("poly(dT)30","T"*30),seg("V","V",placeholder=True),seg("N","N",placeholder=True)]
def capture_scene():
    mrna=[seg("mRNA","X"*34,placeholder=True),seg("poly(A)","A"*30)]; sc=Scene(); sc.strand("mRNA",mrna,label="mRNA")
    sc.anneal("RT primer",rt_primer(),to="mRNA",pair=("poly(dT)30","poly(A)"),label="universal RT primer",unpaired=("V","N")); sc.arrow("RT primer","reverse transcriptase"); return sc
def first_strand(): return Construct([*rt_primer(),seg("antisense cDNA","X"*36,placeholder=True),seg("CCC",rt.UNTEMPLATED_TAIL,"tso")],name="v1 first strand")
def switch_scene():
    first=first_strand(); sc=Scene(); sc.strand("first",list(first),label="first strand"); sc.anneal("bead TSO",bead_tso(),to="first",pair=("rGrGrG","CCC"),label="barcoded bead TSO",above=True,mod5="bead"); sc.arrow("first","RT switches template"); return sc
def amplified_cdna(): return Construct([*bead_tso(),seg("5-prime cDNA","X"*40,placeholder=True),*complement_segments(rt_primer())],name="v1 amplified cDNA")
def adapter_top(): return [seg("Read 2 adapter",il.INDEX1_PRIMER,"r2")]
def adapter_bottom(): return [seg("stem complement",il.STEM_COMPLEMENT,"r2"),seg("3' T","T","r2")]
def adapter_scene():
    sc=Scene(); sc.strand("top",adapter_top(),label="Read 2 adapter"); sc.anneal("bottom",adapter_bottom(),to="top",pair=("stem complement","Read 2 adapter"),label="short strand"); sc.mark("bottom","3' T","ligation overhang"); return sc
SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
def final_library(insert_nt=38):
    lib=Construct([seg("P5",il.P5,"p5"),seg("Read 1 arm",il.TRUSEQ_READ1[4:],"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*10,"umi",placeholder=True),seg("switch spacer",SWITCH_SPACER,"tso"),seg("GGG","GGG","tso"),seg("5-prime cDNA insert","X"*insert_nt,placeholder=True),seg("dA junction","A"),seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),seg("i7 reverse complement","I"*8,"cbc",placeholder=True),seg("P7 reverse complement",il.P7_RC,"p7")],name="10x 5' v1 library")
    problems=sp.verify(lib,SEQ_PRIMERS,required_roles=("Read 1","Index 1 (i7)","Read 2"))
    if problems: raise ValueError("invalid 5' v1 library: "+"; ".join(problems))
    return lib
def read1_layout():
    identifiers=[(x.cycles,x.feature.label) for x in sp.feature_spans(final_library(),SEQ_PRIMERS,{"Read 1":26})]
    return [*identifiers,("27 onward","switch spacer, then transcript 5' end")]

AUTHORITY=('Primary chemistry source: <a href="https://www.10xgenomics.com/support">'
           '10x Genomics Technical Note CG000109 Rev D</a>.')
TABLE_2B="CG000109 Rev D Table 2B (Protocol Option 2: 5' Gene Expression), p. 9"

# CG000109 Rev D Table 2B oligos.  Bases are slices of canonical handles; the
# selftest checks them letter-for-letter against the printed Table 2B strings.
CDNA_FORWARD=PARTIAL_R1                                  # Protocol 4.1 "Partial Read 1"
CDNA_REVERSE=rt.SMART_HANDLE[:20]                        # Protocol 4.1 "Non-Poly(dT) Primer", Table 2B form
CDNA_REVERSE_TABLE_2A=rt.SMART_HANDLE[:22]               # same PN-220106 as printed in Table 2A
INDEX_FORWARD=il.P5+il.TRUSEQ_READ1[4:24]                # Protocol 6.5 "P5 - Partial Read 1" PN-220111
INDEX_REVERSE=revcomp(il.P7_RC)+"N"*8+revcomp(il.INDEX1_PRIMER)[:21]  # "P7 - Sample Index - Partial Read 2" PN-220103

def _oligo(name,seq,tag): return Construct([Segment(name=name,top=seq,tag=tag)],name=name)
def rt_product():
    """CG000109 Table 2B Protocol Step 1.5 'Reverse Transcript Product' (one strand, 5'->3')."""
    return Construct([*rt_primer(),seg("antisense cDNA","X"*36,placeholder=True),seg("CCC",rt.UNTEMPLATED_TAIL,"tso"),
        seg("switch spacer complement",revcomp(SWITCH_SPACER),"tso"),seg("UMI complement","U"*10,"umi",placeholder=True,feature=UMI_FEATURE),
        seg("cell barcode complement","B"*16,"cbc",placeholder=True,feature=CELL_BARCODE),seg("Partial Read 1 complement",revcomp(PARTIAL_R1),"r1")],name="v1 reverse-transcript product")
def ligation_product(insert_nt=38):
    """CG000109 Table 2B Protocol 6.3 'Ligation Product' (bead end + insert + A + Read 2 adaptor)."""
    return Construct([seg("Partial Read 1",PARTIAL_R1,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*10,"umi",placeholder=True),seg("switch spacer",SWITCH_SPACER,"tso"),seg("GGG","GGG","tso"),seg("5-prime cDNA insert","X"*insert_nt,placeholder=True),seg("dA junction","A"),seg("Read 2 adapter",il.INDEX1_PRIMER,"r2")],name="v1 adaptor-ligation product")
def ligation_ready_adapter():
    sc=adapter_scene(); sc.ligation_end("bottom","3' T",kind="cohesive",polarity="3-prime",owner="adapter"); return sc

def protocol_graph():
    graph = ProtocolGraph("10x Chromium Single Cell 5-prime Gene Expression v1", authority=AUTHORITY)
    src = AUTHORITY
    reviewed = []
    cells = graph.input("cells", "cells containing poly(A) RNA", MaterialKind.BIOLOGICAL,
                        state=molecule("cells", (role("poly(A) RNA"),)), primary=True)
    reviewed.append("cells")
    gems = graph.input("gems", "Chromium GEM partitioning with Single Cell 5' gel beads", MaterialKind.CONDITION)
    occupied = graph.product("occupied", "occupied GEMs containing cells", MaterialKind.MOLECULAR, state=cells.state)
    failed = graph.product("failed", "empty, multiplet or unresolved GEM products", MaterialKind.WASTE)
    graph.react("partition", "Partition cells with barcoded gel beads", StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (gems, InputRole.CONDITION)), (occupied, failed),
                changes_molecule=False)
    graph.discard(failed, "does not retain a resolvable single-cell barcode")
    reviewed += ["gems", "partition"]
    primer = graph.input("rt_primer", "Poly-dT RT Primer PN-2000007", MaterialKind.OLIGO,
                         state=Construct(rt_primer(), name="Poly-dT RT Primer"),
                         source=TABLE_2B + ", Protocol Step 1.5")
    enzyme = graph.input("rt", "reverse transcriptase in the GEM Master Mix (not named in CG000109)",
                         MaterialKind.CATALYST)
    first = graph.product("first_strand", "first-strand cDNA with untemplated 3' C tail",
                          MaterialKind.MOLECULAR, state=first_strand(), inferred=True,
                          source=TABLE_2B + " (pre-switch intermediate and CCC tail inferred; the note prints only the switched product)")
    graph.react("reverse_transcribe", "INFERRED — capture and reverse-transcribe poly(A) RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((occupied, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                 (enzyme, InputRole.CATALYST)), (first,))
    tso = graph.input("tso", "Gel Bead Oligo Primer (TSO) PN-220112", MaterialKind.OLIGO,
                      state=Construct(bead_tso(), name="Gel Bead Oligo Primer"),
                      source=TABLE_2B + ", Protocol Step 1.5")
    switched = graph.product("switched", "reverse-transcript product carrying cell barcode and UMI",
                             MaterialKind.MOLECULAR, state=rt_product(),
                             source=TABLE_2B + ", Protocol Step 1.5 'Reverse Transcript Product'")
    graph.react("template_switch", "Template-switch onto the barcoded gel-bead oligo",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (enzyme, InputRole.CATALYST)), (switched,))
    recovery = graph.input("gem_recovery", "GEM break and cDNA cleanup (conditions in CG000086, not revision-matched)",
                           MaterialKind.CONDITION)
    recovered = graph.product("recovered", "recovered barcoded first-strand cDNA",
                              MaterialKind.MOLECULAR, state=switched.state)
    graph.react("recover", "Break GEMs and recover cDNA", StepKind.EXTRACTION,
                ((switched, InputRole.TEMPLATE), (recovery, InputRole.CONDITION)),
                (recovered,), changes_molecule=False)
    pf = graph.input("preamp_forward", "cDNA Primer Mix PN-220106 forward (Partial Read 1)",
                     MaterialKind.OLIGO, state=_oligo("Partial Read 1 primer", CDNA_FORWARD, "r1"),
                     source=TABLE_2B + ", Protocol 4.1")
    pr = graph.input("preamp_reverse", "cDNA Primer Mix PN-220106 reverse (Non-Poly(dT), 20-nt Table 2B form)",
                     MaterialKind.OLIGO, state=_oligo("Non-Poly(dT) primer", CDNA_REVERSE, "tso"),
                     source=TABLE_2B + ", Protocol 4.1")
    pol = graph.input("preamp_polymerase", "cDNA amplification polymerase (not named in CG000109)",
                      MaterialKind.CATALYST)
    amplified = graph.product("amplified", "amplified full-length barcoded cDNA",
                              MaterialKind.MOLECULAR, state=amplified_cdna(),
                              source=TABLE_2B + ", Protocol 4.1 'cDNA Amplified Product'")
    graph.react("preamp", "Amplify cDNA between Partial Read 1 and the Non-Poly(dT) handle", StepKind.PCR,
                ((recovered, InputRole.TEMPLATE), (pf, InputRole.STRUCTURAL),
                 (pr, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.exact(amplified_cdna(),
                    PCRPrimer("Partial Read 1", CDNA_FORWARD), PCRPrimer("Non-Poly(dT)", CDNA_REVERSE),
                    statement=TABLE_2B + ", Protocol 4.1 prints both cDNA Primer Mix oligos; they flank the bead end and the RT-primer handle."))
    fragmentase = graph.input("fragmentase", "enzymatic fragmentation (enzyme not named in CG000109)",
                              MaterialKind.CATALYST)
    fragments = graph.product("fragments", "fragmented cDNA (bead-end, internal and handle-end pieces)",
                              MaterialKind.MOLECULAR)
    graph.react("fragment", "Enzymatically fragment amplified cDNA", StepKind.FRAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (fragmentase, InputRole.CATALYST)), (fragments,))
    end_mix = graph.input("end_mix", "end-repair and A-tailing enzymes", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "end-repaired dA-tailed cDNA fragments", MaterialKind.MOLECULAR)
    graph.react("end_prep", "End-repair and dA-tail cDNA fragments", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (end_mix, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("read2_adapter", "Adaptor (Read 2) PN-220026, 3'-T partially double-stranded",
                          MaterialKind.ADAPTER, state=adapter_scene(), source=TABLE_2B + ", Protocol 6.3")
    ligase = graph.input("ligase", "DNA ligase (not named in CG000109)", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "Read-2-adaptor-ligated fragments", MaterialKind.MOLECULAR,
                            state=ligation_product(), source=TABLE_2B + ", Protocol 6.3 'Ligation Product'")
    graph.react("ligate", "Ligate the Read 2 adaptor to dA-tailed fragments", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.known("cohesive", ligation_ready_adapter(),
                    statement=TABLE_2B + ", Protocol 6.3 prints the 33-nt top strand paired over 12 nt by 3'-TCTAGCCTTCTCG-5', leaving a single 3'-T; adaptor 5' phosphorylation is not stated."))
    index_f = graph.input("index_forward", "Sample Index PCR Primer PN-220111 (P5 - Partial Read 1)",
                          MaterialKind.OLIGO, state=_oligo("P5 - Partial Read 1", INDEX_FORWARD, "p5"),
                          source=TABLE_2B + ", Protocol 6.5")
    index_r = graph.input("index_reverse", "Sample Index PCR Primer PN-220103 (P7 - Sample Index - Partial Read 2)",
                          MaterialKind.OLIGO, state=_oligo("P7 - Sample Index - Partial Read 2", INDEX_REVERSE, "p7"),
                          source=TABLE_2B + ", Protocol 6.5 (8-nt index well not specified)")
    index_pol = graph.input("index_polymerase", "sample-index PCR polymerase (not named in CG000109)",
                            MaterialKind.CATALYST)
    final = final_library()
    library = graph.product("library", "single-index 5' Gene Expression library", MaterialKind.LIBRARY,
                            state=final, source=TABLE_2B + ", Protocol 6.5 'Sample Index PCR Product'")
    graph.react("index_pcr", "Select the bead-derived end and add P5, i7 and P7", StepKind.PCR,
                ((adapted, InputRole.TEMPLATE), (index_f, InputRole.STRUCTURAL),
                 (index_r, InputRole.STRUCTURAL), (index_pol, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(ligation_product(),
                    PCRPrimer("P5 - Partial Read 1", INDEX_FORWARD, anneal_nt=13),
                    PCRPrimer("P7 - Sample Index - Partial Read 2", INDEX_REVERSE, anneal_nt=21),
                    statement=TABLE_2B + ", Protocol 6.5: the P5 primer's 3' 13 nt anneal to Partial Read 1 and the P7 primer's 3' 21 nt to the ligated adaptor (computed from the printed oligos)."))
    graph.sequence(library, SequencingPlan(library=final, primers=SEQ_PRIMERS,
                   acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
                   read_lengths={"Read 1":26,"Index 1 (i7)":8,"Read 2":98}))
    reviewed += ["rt_primer", "rt", "reverse_transcribe", "tso", "template_switch",
                 "gem_recovery", "recover", "preamp_forward", "preamp_reverse",
                 "preamp_polymerase", "preamp", "fragmentase", "fragment", "end_mix",
                 "end_prep", "read2_adapter", "ligase", "ligate", "index_forward",
                 "index_reverse", "index_polymerase", "index_pcr", "library"]
    graph.declare_reviewed_scope(reviewed, src)
    graph.annotate(SourceDiscrepancy(
        "preamp", "Two printed lengths for the cDNA Primer Mix reverse primer",
        (SourceClaim("CG000109 Rev D Table 2A, Protocol 4.1 (V(D)J from amplified cDNA), p. 8",
                     "PN-220106 Non-Poly(dT) primer 5'-AAGCAGTGGTATCAACGCAGAG-3' (22 nt)."),
         SourceClaim("CG000109 Rev D Table 2B, Protocol 4.1 (5' Gene Expression), p. 9",
                     "PN-220106 Non-Poly(dT) primer 5'-AAGCAGTGGTATCAACGCAG-3' (20 nt).")),
        "The model uses the 20-nt Table 2B form, because Table 2B is the 5' Gene Expression table; both forms anneal to the same RT-primer handle.",
        "Both tables name the same part number in the same Technical Note revision; the extra AG is the next two handle bases.",
        "Which oligo PN-220106 actually contains is not resolved by CG000109. The later CG000086 Rev K appendix (Protocol Step 3.1) prints the 22-nt form, but a later revision cannot establish the 2017 v1 kit, and the revision-matched CG000086 was not located.",
        slot=AnnotationSlot.INTRO))
    return graph.require_publishable()
