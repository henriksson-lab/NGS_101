"""Construct model for 10x Chromium Single Cell 3' Gene Expression v3 (CG000183)."""
from __future__ import annotations
import illumina as il
import nextera
import rt
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,feature,revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import Disposition, InputRole, MaterialKind, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role
from protocol_graph import ProtocolGraph, SourceClaim, SourceDiscrepancy

PARTIAL_R1=il.TRUSEQ_READ1[11:]
TSO_DNA=rt.SMART_HANDLE+"ACAT"
CELL_BARCODE=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x-chromium-3prime-v3")
UMI_FEATURE=feature("umi","umi","random")
SAMPLE_INDEX=feature("sample_index_i7","sample_index","unknown")
FEATURES={"cell barcode":CELL_BARCODE,"UMI":UMI_FEATURE,"i7 reverse complement":SAMPLE_INDEX}
def seg(name,top,tag=None,**kw):
    kw.setdefault("feature",FEATURES.get(name)); return Segment(name=name,top=top,tag=tag,**kw)
def bead_oligo(): return [seg("Partial Read 1",PARTIAL_R1,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("poly(dT)30","T"*30)]
def tso(): return [seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("ACAT","ACAT","tso"),seg("rGrGrG","GGG","tso")]
def capture_seq1_bead_oligo():
    """Second CG000183 gel-bead primer, printed as "not used for 3' Gene Expression"."""
    return [seg("Nextera Read 1 (Read 1N)",nextera.ADAPTOR_S5[11:],"r1"),
            seg("cell barcode","B"*16,"cbc",placeholder=True),
            seg("UMI","U"*12,"umi",placeholder=True),
            seg("Capture Sequence 1","TTGCTAGGACCGGCCTTAAAGC")]
def capture_scene():
    mrna=[seg("mRNA","X"*32,placeholder=True),seg("poly(A)","A"*30)]; sc=Scene(); sc.strand("mRNA",mrna,label="mRNA")
    sc.anneal("bead primer",bead_oligo(),to="mRNA",pair=("poly(dT)30","poly(A)"),label="released gel-bead primer"); sc.arrow("bead primer","reverse transcriptase"); return sc
def first_strand(): return Construct([*bead_oligo(),seg("antisense cDNA","X"*34,placeholder=True),seg("CCC",rt.UNTEMPLATED_TAIL,"tso")],name="v3 first strand")
def switch_scene():
    first=first_strand(); sc=Scene(); sc.strand("first",list(first),label="first strand"); sc.anneal("TSO",tso(),to="first",pair=("rGrGrG","CCC"),label="TSO",above=True); sc.arrow("first","RT switches template"); return sc
def amplified_cdna(): return Construct([*bead_oligo(),seg("antisense cDNA","X"*34,placeholder=True),seg("CCC",rt.UNTEMPLATED_TAIL,"tso"),seg("copied TSO",revcomp(TSO_DNA),"tso")],name="v3 amplified cDNA")
def adapter_top(): return [seg("Partial Read 2",il.INDEX1_PRIMER[:-1],"r2")]
def adapter_bottom(): return [seg("stem complement",il.STEM_COMPLEMENT,"r2"),seg("3' T","T","r2")]
def adapter_scene():
    """Adaptor Oligos PN-2000094 as printed in CG000183 Rev A Appendix p. 54, with the 3'-T end."""
    sc=Scene(); sc.strand("top",adapter_top(),label="Read 2 adapter"); sc.anneal("bottom",adapter_bottom(),to="top",pair=("stem complement","Partial Read 2"),label="short strand")
    sc.ligation_end("bottom","3' T",kind="cohesive",polarity="3-prime",owner="adapter",text="3'-T ligation overhang"); return sc
def read2_adapter_end():
    return AdapterLigationGeometry.known("cohesive",adapter_scene(),
        statement=("CG000183 Rev A Appendix p. 54 (Protocol Step 3.3) prints Adaptor Oligos "
                   "PN-2000094 as the 32-nt Partial Read 2 strand paired over 12 nt by "
                   "3'-TCTAGCCTTCTCG-5', leaving a single 3'-T for the dA-tailed fragment; "
                   "5' phosphorylation is not stated."))
SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
READ_LENGTHS={"Read 1":28,"Index 1 (i7)":8,"Read 2":91}
def final_library(insert_nt=36):
    lib=Construct([seg("P5",il.P5,"p5"),seg("Read 1 arm",il.TRUSEQ_READ1[4:],"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("poly(dT)30","T"*30),seg("cDNA insert","X"*insert_nt,placeholder=True),seg("dA junction","A"),seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),seg("i7 reverse complement","I"*8,"cbc",placeholder=True),seg("P7 reverse complement",il.P7_RC,"p7")],name="10x 3' v3 library")
    problems=sp.verify(lib,SEQ_PRIMERS,required_roles=("Read 1","Index 1 (i7)","Read 2"))
    if problems: raise ValueError("invalid v3 library: "+"; ".join(problems))
    return lib
def read1_layout(): return [(x.cycles,x.feature.label) for x in sp.feature_spans(final_library(),SEQ_PRIMERS,{"Read 1":28})]


def switched_cdna():
    """Single-stranded template-switched first strand; cDNA PCR makes the duplex later."""
    sc = Scene()
    sc.strand("switched", list(amplified_cdna()), label="v3 template-switched first-strand cDNA")
    return sc


def cdna_primers():
    return (PCRPrimer("cDNA forward primer", PARTIAL_R1),
            PCRPrimer("cDNA reverse primer", TSO_DNA[:22]))


def index_primers():
    # Source prints an N8 sample-index window; it represents one oligo in a well mix.
    return (PCRPrimer("sample-index forward primer", il.NEBNEXT_UNIVERSAL_PRIMER[:49],
                      anneal_nt=13),
            PCRPrimer("sample-index reverse primer",
                      il.P7+"N"*8+il.TRUSEQ_READ2[:21]))


def pre_index_library(insert_nt=36):
    return Construct([*bead_oligo(),
                      seg("cDNA insert", "X"*insert_nt, placeholder=True),
                      seg("dA junction", "A"), *adapter_top()],
                     name="v3 adapter-ligated cDNA fragment")


def _cleanup(g, reviewed, ident, label, material, reagent, wastes, *,
             kind=MaterialKind.MOLECULAR, kept_id=None, kept_label=None):
    """A CG000183 Rev A bead cleanup: the kept fraction and each removed fraction."""
    beads=g.input(ident+"_beads", reagent, MaterialKind.SELECTION)
    kept=g.product(kept_id or ident+"_kept", kept_label or label+" (retained DNA)", kind,
                   state=material.state)
    removed=[g.product(ident+"_waste%d" % i, text, MaterialKind.WASTE)
             for i, (text, _) in enumerate(wastes, 1)]
    g.react(ident, label, StepKind.CLEANUP,
            ((material,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
            (kept,*removed), changes_molecule=False)
    for waste, (_, why) in zip(removed, wastes):
        g.discard(waste, why)
    reviewed.extend((beads.id, ident, kept.id, *(w.id for w in removed)))
    return kept


SUPERNATANT=("bead supernatant (unbound primers, reagents and fragments below the cut-off)",
             "supernatant removed from the magnet-held SPRIselect beads")
LARGE=("fragments bound at 0.6X SPRIselect (above the upper size cut-off)",
       "first-SPRI beads discarded after the supernatant is transferred")
SMALL=("supernatant after 0.8X SPRIselect (below the lower size cut-off)",
       "supernatant removed from the magnet-held SPRIselect beads")


def _source_pcr(g, reviewed, prefix, label, template, product_state, primers, *, library=False,
                site_template=None):
    """The source prints both primers; locate their sites rather than claiming absence.

    ``site_template`` is the duplex on which both primer sites are drawn when the input
    material is single-stranded (the switched first strand gains its complement in cycle 1).
    """
    forward, reverse=primers
    f=g.input(prefix+"_forward", forward.name, MaterialKind.OLIGO,
              state=Construct([Segment(forward.name, forward.sequence)]))
    r=g.input(prefix+"_reverse", reverse.name, MaterialKind.OLIGO,
              state=Construct([Segment(reverse.name, reverse.sequence)]))
    polymerase=g.input(prefix+"_polymerase", label+" polymerase", MaterialKind.CATALYST)
    product=g.product("index_pcr_product" if library else "amplified", label+" product",
                      MaterialKind.MOLECULAR, state=product_state)
    g.react(prefix, label, StepKind.PCR,
            ((template,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),
             (r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)), (product,),
            geometry=PCRGeometry.exact(site_template or template.state,forward,reverse,
                statement="Both primer architectures are printed in CG000183 Rev A Appendix pp. 54–55."))
    reviewed.extend((f.id,r.id,polymerase.id,prefix,product.id))
    return product


def _source_ligation(g, reviewed, tailed):
    adapter=g.input("read2_adapter", "3-prime-dT Read-2 adapter", MaterialKind.ADAPTER,
                    state=adapter_scene())
    ligase=g.input("ligase", "adapter ligase", MaterialKind.CATALYST)
    adapted=g.product("adapted", "Read-2-adapted barcode-bearing fragment",
                      MaterialKind.MOLECULAR,state=pre_index_library())
    g.react("ligate", "Ligate the Read-2 adapter", StepKind.LIGATION,
            ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
             (ligase,InputRole.CATALYST)),(adapted,),geometry=read2_adapter_end())
    reviewed.extend((adapter.id,ligase.id,"ligate"))
    return adapted

def protocol_graph():
    authority=('Primary chemistry source: <a href="https://www.10xgenomics.com/support">'
               '10x Genomics User Guide CG000183 Rev A</a>.')
    graph = ProtocolGraph("10x Chromium Single Cell 3-prime Gene Expression v3", authority=authority)
    # Reaction choices and source coverage belong to this guide-specific declaration.
    reviewed = ["cells", "gems", "partition", "rt_primer", "rt", "reverse_transcribe",
                "tso", "template_switch", "gem_recovery", "recover", "fragmentase",
                "fragment", "end_mix", "end_prep"]
    final = final_library()
    cells = graph.input("cells", "cells containing poly(A) RNA", MaterialKind.BIOLOGICAL,
                        state=molecule("cells", (role("poly(A) RNA"),)), primary=True)
    gems = graph.input("gems", "Chromium GEM partitioning condition", MaterialKind.CONDITION)
    occupied = graph.product("occupied", "occupied GEMs containing cells", MaterialKind.MOLECULAR,
                             state=cells.state)
    failed = graph.product("failed", "empty, multiplet or unresolved GEM products", MaterialKind.MOLECULAR)
    graph.react("partition", "Partition cells with barcoded gel beads", StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (gems, InputRole.CONDITION)),
                (occupied, failed), changes_molecule=False)
    graph.mark_terminal(Terminal(failed.id, Disposition.UNRESOLVED,
        "GEM contents are pooled after GEM break. Barcode resolution is a downstream analysis "
        "criterion; the guide does not establish physical removal of these GEM classes here."))
    primer = graph.input("rt_primer", "cell-barcode/UMI oligo(dT) gel-bead primer", MaterialKind.OLIGO,
                         state=Construct(bead_oligo(), name="gel-bead primer"))
    enzyme = graph.input("rt", "template-switching reverse transcriptase", MaterialKind.CATALYST)
    first_state = Scene()
    first_state.strand("first", list(first_strand()), label="first-strand transcript cDNA")
    first = graph.product("first_strand", "first-strand transcript cDNA", MaterialKind.MOLECULAR,
                          state=first_state)
    graph.react("reverse_transcribe", "Capture and reverse-transcribe poly(A) RNA", StepKind.REVERSE_TRANSCRIPTION,
                ((occupied, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                 (enzyme, InputRole.CATALYST)), (first,))
    switch_oligo = graph.input("tso", "template-switch oligo", MaterialKind.OLIGO,
                              state=Construct(tso(), name="template-switch oligo"))
    switched = graph.product("switched", "template-switched cell-barcoded cDNA", MaterialKind.MOLECULAR,
                             state=switched_cdna())
    graph.react("template_switch", "Install the opposite cDNA amplification handle", StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (switch_oligo, InputRole.STRUCTURAL),
                 (enzyme, InputRole.CATALYST)), (switched,))
    recovery = graph.input("gem_recovery", "GEM break and cDNA recovery condition", MaterialKind.CONDITION)
    recovered = graph.product("recovered", "recovered cell-barcoded cDNA", MaterialKind.MOLECULAR,
                              state=switched.state)
    graph.react("recover", "Break GEMs and recover pooled cDNA", StepKind.EXTRACTION,
                ((switched, InputRole.TEMPLATE), (recovery, InputRole.CONDITION)),
                (recovered,), changes_molecule=False)
    recovered = _cleanup(graph, reviewed, "silane_cleanup", "Post GEM-RT cleanup with Dynabeads MyOne Silane",
        recovered, "Dynabeads MyOne SILANE beads",
        (("Silane-bead supernatant (leftover GEM-RT reagents and primers)",
          "removed from the magnet-held Silane beads"),))
    amplified = _source_pcr(graph, reviewed, "preamp", "Amplify recovered full-length cDNA",
                            recovered, amplified_cdna(), cdna_primers(),
                            site_template=amplified_cdna())
    amplified = _cleanup(graph, reviewed, "cdna_spri", "cDNA SPRIselect cleanup (0.6X)",
                         amplified, "SPRIselect reagent, 0.6X", (SUPERNATANT,))
    fragmentase = graph.input("fragmentase", "CG000183 Rev A Step 3.1 Fragmentation Enzyme", MaterialKind.CATALYST)
    fragments = graph.product("fragments", "fragmented barcode-bearing cDNA", MaterialKind.MOLECULAR)
    graph.react("fragment", "Fragment amplified cDNA", StepKind.FRAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (fragmentase, InputRole.CATALYST)), (fragments,))
    end_mix = graph.input("end_mix", "CG000183 Rev A Step 3.1 end-repair and dA-tailing mix", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "end-repaired dA-tailed cDNA fragments", MaterialKind.MOLECULAR)
    graph.react("end_prep", "End-repair and dA-tail cDNA fragments", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (end_mix, InputRole.CATALYST)), (tailed,))
    tailed = _cleanup(graph, reviewed, "fragment_size_selection",
                      "Double-sided SPRIselect size selection (0.6X, then 0.8X)",
                      tailed, "SPRIselect reagent, 0.6X then 0.8X", (LARGE, SMALL))
    adapted = _source_ligation(graph, reviewed, tailed)
    adapted = _cleanup(graph, reviewed, "ligation_spri", "Post-ligation SPRIselect cleanup (0.8X)",
                       adapted, "SPRIselect reagent, 0.8X", (SUPERNATANT,))
    indexed = _source_pcr(graph, reviewed, "index_pcr", "Select barcode-bearing end and add sample index",
                          adapted, final, index_primers(), library=True)
    library = _cleanup(graph, reviewed, "index_size_selection",
                       "Double-sided SPRIselect size selection (0.6X, then 0.8X)",
                       indexed, "SPRIselect reagent, 0.6X then 0.8X", (LARGE, SMALL),
                       kind=MaterialKind.LIBRARY, kept_id="library",
                       kept_label="size-selected indexed v3 library")
    graph.sequence(library, SequencingPlan(library=final, primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS), read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope(reviewed, authority)
    graph.require_publishable()
    graph.annotate(SourceDiscrepancy(
        "rt_primer", "CG000183 poly(dT) anchor disagreement",
        (SourceClaim("CG000183 Rev A, overview figure (Poly(dT)VN)",
                     "The bead primer is labelled Poly(dT)VN."),
         SourceClaim("CG000183 Rev A, Appendix, Oligonucleotide Sequences",
                     "The printed gel-bead primer ends with T30 and no VN bases.")),
        "The construct uses the appendix's printed T30 sequence.",
        "The appendix provides a base-resolved oligo, whereas the overview gives a role label.",
        "The physical bead-primer anchor may include VN bases omitted by the appendix."))
    graph.annotate(SourceDiscrepancy(
        "read2_adapter", "CG000183 adapter terminal-base disagreement",
        (SourceClaim("CG000183 Rev A, Appendix p. 54, Adaptor Oligos PN-2000094",
                     "The top-strand oligo is 32 bases and ends CAGTCA."),
         SourceClaim("CG000183 Rev A, Appendix p. 54, Ligation Product",
                     "The adapter region in the printed ligation product ends CAGTCAC.")),
        "The adapter input and pre-index state use the printed 32-base oligo; the final library contains the terminal C.",
        "The sample-index reverse primer on p. 55 restores this C when its full sequence is copied.",
        "The source does not establish whether the extra C in its ligation product is a drawing error or the oligo list omits a base."))
    return graph
