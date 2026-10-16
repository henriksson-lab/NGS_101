"""Molecular construct model for Seq-Well S3 (Hughes et al., 2020)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (AnnotationSlot, Highlight, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "random")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "whitelist",
                     whitelist="Nextera N700 index set; selected member not specified")
FEATURES = {"cell barcode": CELL_BARCODE, "UMI": UMI_FEATURE,
            "i7 index read": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


BEAD_PREFIX = "T" * 7
BEAD_SUFFIX = "AC"
P5_SPACER = "GCCTGTCCGCGG"
S3_RANDOM_TAIL = "NNNGGNNNB"
I7_INDEX_READ = "N" * 8
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"
P5_HYBRID = il.P5 + P5_SPACER + rt.SMART_HANDLE + BEAD_SUFFIX
CUSTOM_R1 = P5_SPACER + rt.SMART_HANDLE + BEAD_SUFFIX


def bead_oligo() -> Construct:
    return Construct([
        seg("bead prefix", BEAD_PREFIX), seg("SMART handle", rt.SMART_HANDLE, "tso"),
        seg("bead constant", BEAD_SUFFIX),
        seg("cell barcode", "B" * 12, "cbc", placeholder=True),
        seg("UMI", "N" * 8, "umi", placeholder=True), seg("poly(T)", "T" * 30),
    ], name="Seq-Well S3 bead oligo")


def first_strand() -> Construct:
    return Construct([
        *bead_oligo(), seg("copied transcript", "X" * 28, placeholder=True),
    ], name="Seq-Well S3 bead-extended first-strand cDNA")


def p5_primer() -> Construct:
    return Construct([
        seg("P5", il.P5, "p5"), seg("spacer", P5_SPACER),
        seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("AC", BEAD_SUFFIX),
    ], name="P5-SMART hybrid; terminal *A*C phosphorothioate bonds")


def n7xx_primer() -> Construct:
    return Construct([
        seg("P7", il.P7, "p7", inferred=True),
        seg("i7 index read", I7_INDEX_READ, "cbc", placeholder=True),
        seg("s7", nx.S7, "s7", inferred=True),
    ], name="N7xx library-PCR primer; unassigned eight-base index")


def tso_oligo() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
            seg("rGrGrG", "rGrGrG", placeholder=True)]


def randomer() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GA", "GA"),
            seg("random priming tail", S3_RANDOM_TAIL, placeholder=True)]


def rt_scene(insert_nt: int = 28) -> Scene:
    mrna = [seg("transcript", "X" * insert_nt, placeholder=True),
            seg("poly(A)", "A" * 30)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA", mod3="poly(A)")
    sc.anneal("bead", list(bead_oligo()), to="mRNA", pair=("poly(T)", "poly(A)"),
              label="bead primer")
    sc.arrow("bead", "reverse transcription; template switching may occur")
    return sc


def random_second_strand_scene() -> Scene:
    # Placeholder text is deliberately symmetric under the Scene's reverse placement:
    # the exact random landing sequence is molecule-specific and is not asserted here.
    template = [seg("first-strand cDNA", "X" * 20, placeholder=True),
                seg("randomer landing", "BNNNCCNNN", placeholder=True)]
    sc = Scene()
    sc.strand("first", template, label="bead-bound first strand")
    sc.anneal("randomer", randomer(), to="first",
              pair=("random priming tail", "randomer landing"), label="S3 randomer")
    sc.arrow("randomer", "Klenow exo-minus makes the second strand toward the bead")
    return sc


def wta_product(insert_nt: int = 28) -> Construct:
    return Construct([
        seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("bead constant",BEAD_SUFFIX),
        seg("cell barcode","B"*12,"cbc",placeholder=True),
        seg("UMI","N"*8,"umi",placeholder=True),seg("poly(T)","T"*30),
        seg("cDNA","X"*insert_nt,placeholder=True),
        seg("copied GA",revcomp("GA")),
        seg("copied SMART handle",revcomp(rt.SMART_HANDLE),"tso")],
        name="Seq-Well S3 handle-flanked WTA product")


def selected_tagmented_strand(insert_nt: int = 28) -> Construct:
    return Construct([
        seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("bead constant",BEAD_SUFFIX),
        seg("cell barcode","B"*12,"cbc",placeholder=True),
        seg("UMI","N"*8,"umi",placeholder=True),seg("poly(T)","T"*30),
        seg("cDNA","X"*insert_nt,placeholder=True),
        seg("mosaic end reverse complement",nx.ME_RC,"me",inferred=True),
        seg("s7 reverse complement",nx.S7_RC,"s7",inferred=True)],
        name="Seq-Well S3 selected bead-end fragment")


SEQ_PRIMERS = (
    sp.custom("Read 1", "Seq-Well S3 custom Read 1", CUSTOM_R1,
              "Hughes et al. 2020 Key Resources Table"),
    sp.NEXTERA["I1"], sp.NEXTERA["R2"],
)
RUN_ROLES = ("Read 1", "Index 1 (i7)", "Read 2")


def final_library(insert_nt: int = 28) -> Construct:
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("custom spacer", P5_SPACER),
        seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("bead constant", BEAD_SUFFIX),
        seg("cell barcode", "B" * 12, "cbc", placeholder=True),
        seg("UMI", "N" * 8, "umi", placeholder=True), seg("poly(T)", "T" * 30),
        seg("cDNA", "X" * insert_nt, placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me", inferred=True),
        seg("s7 reverse complement", nx.S7_RC, "s7", inferred=True),
        seg("i7 index read", I7_INDEX_READ, "cbc", placeholder=True),
        seg("P7 reverse complement", il.P7_RC, "p7", inferred=True),
    ], name="Seq-Well S3 sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=RUN_ROLES)
    if problems:
        raise ValueError("invalid final Seq-Well S3 construct: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def protocol_graph() -> ProtocolGraph:
    source="Defining source: https://doi.org/10.1016/j.immuni.2020.09.015"
    author_protocol="Hughes et al. Seq-Well S3 Master Protocol (June 2019), pp. 18–24, 39: https://shaleklab.com/wp-content/uploads/2024/02/Hughes.et_.al_.Master.Protocol.pdf"
    graph=ProtocolGraph("Seq-Well S3",authority=source)
    cells=graph.input("cells","single-cell suspension",MaterialKind.BIOLOGICAL,primary=True)
    beads=graph.input("barcoded_beads","barcoded capture beads",MaterialKind.ADAPTER,
                      state=bead_oligo(),source=author_protocol)
    membrane=graph.input("membrane","semipermeable nanowell membrane",MaterialKind.CONDITION)
    captured=graph.product("captured_rna","nanowell-confined bead-captured poly(A) RNA",
                           MaterialKind.MOLECULAR,state=rt_scene())
    empty=graph.product("empty_wells","unproductive cell/bead loading outcomes",
                        MaterialKind.WASTE)
    graph.react("partition_capture","Load cells and beads into membrane-sealed nanowells",
                StepKind.PARTITION,
                ((cells,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (membrane,InputRole.CONDITION)),(captured,empty))
    graph.discard(empty,"wells without a productive cell-bead pairing do not yield library")
    rt_enzyme=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    first=graph.product("first_strand","bead-barcoded first-strand cDNA",
                        MaterialKind.MOLECULAR,state=first_strand())
    graph.react("reverse_transcribe","Synthesize bead-barcoded first-strand cDNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (rt_enzyme,InputRole.CATALYST)),(first,))
    tso=graph.input("tso","template-switching oligo",MaterialKind.OLIGO,
                    state=tuple(tso_oligo()))
    switched=graph.product("switched_cdna","template-switched cDNA molecules",
                           MaterialKind.MOLECULAR)
    unswitched=graph.product("unswitched_cdna","cDNA molecules without a template switch",
                             MaterialKind.MOLECULAR)
    graph.react("template_switch","Resolve template-switched and unswitched cDNA fractions",
                StepKind.TEMPLATE_SWITCHING,
                ((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
                 (rt_enzyme,InputRole.CATALYST)),(switched,unswitched))
    exonuclease=graph.input("exonuclease_i","Exonuclease I",MaterialKind.CATALYST)
    cleaned=graph.product("exo_cleaned_cdna","bead-bound first strands after Exonuclease I",
                          MaterialKind.MOLECULAR)
    unused_primers=graph.product("unused_capture_primers","digested unextended capture primers",
                                 MaterialKind.WASTE)
    graph.react("exonuclease_cleanup","Digest unextended oligos with Exonuclease I",
                StepKind.DIGESTION,
                ((switched,InputRole.TEMPLATE),(unswitched,InputRole.TEMPLATE),
                 (exonuclease,InputRole.CATALYST)),(cleaned,unused_primers),
                source=author_protocol)
    graph.discard(unused_primers,"unextended capture primers cannot seed the S3 randomer")
    alkali=graph.input("alkali","RNA-stripping alkali",MaterialKind.CONDITION)
    te_wash=graph.input("te_wash","1xTE wash after NaOH",MaterialKind.CONDITION,
                        source=author_protocol)
    stripped=graph.product("stripped_first_strand","RNA-stripped bead-bound first strand",
                           MaterialKind.MOLECULAR)
    removed_rna=graph.product("removed_rna","hydrolysed RNA template",MaterialKind.WASTE)
    graph.react("strip_rna","Strip RNA with NaOH; wash bead cDNA with 1xTE",
                StepKind.STRAND_SEPARATION,
                ((cleaned,InputRole.TEMPLATE),(alkali,InputRole.CONDITION),
                 (te_wash,InputRole.CONDITION)),
                (stripped,removed_rna),source=author_protocol)
    graph.discard(removed_rna,"the RNA template is removed before random second-strand synthesis")
    randomer_input=graph.input("s3_randomer","S3 random second-strand primer",
                               MaterialKind.OLIGO,state=tuple(randomer()))
    klenow=graph.input("klenow","Klenow exo-minus polymerase",MaterialKind.CATALYST)
    rescued=graph.product("random_primed_cdna","random-primed second-strand cDNA molecules",
                          MaterialKind.MOLECULAR,state=wta_product())
    graph.react("random_second_strand","Make a handle-tagged second strand on bead-bound cDNA",
                StepKind.EXTENSION,
                ((stripped,InputRole.TEMPLATE),(randomer_input,InputRole.STRUCTURAL),
                 (klenow,InputRole.CATALYST)),(rescued,),source=author_protocol)
    smart=graph.input("smart_primer","SMART PCR primer",MaterialKind.OLIGO,
                      state=Construct([seg("SMART handle",rt.SMART_HANDLE,"tso")]),
                      source="Hughes et al. 2020 Key Resources, SMART PCR Primer")
    preamp_pol=graph.input("preamp_polymerase","WTA PCR polymerase",MaterialKind.CATALYST)
    amplified=graph.product("amplified_cdna","amplified S3 WTA cDNA",MaterialKind.MOLECULAR,
                            state=wta_product())
    graph.react("preamp","Amplify handle-flanked cDNA",StepKind.PCR,
                ((rescued,InputRole.TEMPLATE),(smart,InputRole.STRUCTURAL),
                 (preamp_pol,InputRole.CATALYST)),(amplified,),
                geometry=PCRGeometry.symmetric(
                    wta_product(),PCRPrimer("SMART PCR primer",rt.SMART_HANDLE),
                    statement="The Key Resources Table prints the SMART primer used at both handle boundaries."))
    transposome=graph.input("transposome","Nextera transposome",MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    selected=selected_tagmented_strand()
    tagged=graph.product("selected_tagmented",selected.name,MaterialKind.MOLECULAR,state=selected)
    waste=graph.product("tagment_waste","tagmentation products without the bead end",
                        MaterialKind.WASTE)
    graph.react("tagment","INFERRED — standard Nextera arms on tagmented cDNA",StepKind.TAGMENTATION,
                ((amplified,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,waste))
    graph.discard(waste,"fragments lacking the bead-derived handle are not PCR-selected")
    p5=graph.input("p5_primer","P5-SMART selective primer",MaterialKind.OLIGO,
                   state=p5_primer(),source="Hughes et al. 2020 Key Resources, P5-SMART Hybrid Oligo (*A*C)")
    p7=graph.input("p7_primer","N7xx indexed primer",MaterialKind.OLIGO,
                   state=n7xx_primer(),source="Hughes et al. 2020 Sequencing Library Preparation: N700 family; official Illumina Nextera PCR arms")
    library_pol=graph.input("library_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    graph.react("library_pcr","INFERRED — select the bead-end/standard Nextera fragment and add P5/P7",StepKind.PCR,
                ((tagged,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
                 (p7,InputRole.STRUCTURAL),(library_pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.exact(
                    selected,PCRPrimer("P5-SMART",P5_HYBRID,
                                       anneal_nt=len(rt.SMART_HANDLE)+len(BEAD_SUFFIX)),
                    PCRPrimer("N7xx",il.P7+I7_INDEX_READ+nx.S7,anneal_nt=len(nx.S7)),
                    statement="Hughes Key Resources prints P5-SMART; Methods names N700. The N7xx boundary uses the official Illumina Nextera PCR arm, with its eight-base member left unassigned."))
    lengths={"Read 1":20,"Index 1 (i7)":8,"Read 2":50}
    graph.sequence(library,SequencingPlan(
        library=final,primers=SEQ_PRIMERS,acquired_roles=RUN_ROLES,read_lengths=lengths,
        read_uses={"Read 1":"extract","Index 1 (i7)":"demultiplex","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "cells","barcoded_beads","membrane","partition_capture","reverse_transcriptase",
        "reverse_transcribe","tso","template_switch","exonuclease_i","exonuclease_cleanup",
        "alkali","te_wash","strip_rna","s3_randomer","klenow","random_second_strand",
        "smart_primer","preamp_polymerase",
        "preamp","transposome","tn5","tagment","p5_primer","p7_primer",
        "library_polymerase","library_pcr","library"),source)
    graph.annotate(SourceDiscrepancy(
        target="s3_randomer", title="Randomer length in the supplement",
        claims=(
            SourceClaim("Hughes et al. Supplementary oligo table, S3 Randomer",
                        "The tail is NNNGGNNNB: nine positions, seven degenerate and two fixed G."),
            SourceClaim("Supplementary Figure S1A legend", "Calls this a random octamer.")),
        model_choice="The molecular model uses the printed nine-position NNNGGNNNB tail.",
        evidence="The oligo table gives the actual synthesis sequence; the pairing scene uses that tail.",
        uncertainty="The legend's intended meaning of octamer is not explained.",
        slot=AnnotationSlot.INTRO))
    graph.annotate(SourceDiscrepancy(
        target="strip_rna", title="Ambiguous TE wash concentration",
        claims=(
            SourceClaim("Hughes et al. 2020 STAR Methods, second-strand synthesis; PMC JATS line 491",
                        "The text says 1 M TE."),
            SourceClaim("Hughes et al. 2020 Document S1, S3 Master Protocol p. 19; mmc1 extracted line 3644",
                        "The detailed wash instruction says 1xTE."),
            SourceClaim("Hughes et al. 2019 author S3 Master Protocol p. 19, step 4",
                        "The earlier author protocol also says 1xTE.")),
        model_choice="The graph follows the detailed 1xTE wash after 0.1 M NaOH stripping.",
        evidence="Both stepwise author protocols state 1xTE; the 2020 prose's 1 M TE is inconsistent with them.",
        uncertainty="The published prose has not been formally corrected; the reason for its 1 M wording is unknown.",
        slot=AnnotationSlot.INTRO))
    graph.annotate(Highlight(
        target="exonuclease_cleanup",title="Exonuclease I precedes S3 second-strand synthesis",
        explanation="The author S3 Master Protocol pp. 18–19 treats the beads with Exonuclease I, then NaOH strips RNA before dN-SMRT/Klenow second-strand synthesis.",
        slot=AnnotationSlot.AFTER_WORKFLOW))
    graph.annotate(SourceDiscrepancy(
        target="preamp",title="KAPA master-mix strength",
        claims=(
            SourceClaim("Hughes et al. 2020 STAR Methods; PMC JATS line 495",
                        "The WTA recipe says 25 µL of KAPA 5X Mastermix in a 40 µL mix."),
            SourceClaim("Hughes et al. 2020 Key Resources; PMC JATS lines 393 and 409",
                        "The named KK2602 reagent is KAPA 2x HiFi HotStart PCR mix."),
            SourceClaim("Hughes et al. 2020 Document S1, S3 Master Protocol p. 19; mmc1 line 3663",
                        "The stepwise recipe says 25 µL of 2X KAPA HiFi Hotstart Readymix.")),
        model_choice="The graph names KAPA HiFi PCR without assigning a stock concentration; the detailed 2X recipe is used in the note.",
        evidence="The named catalogue reagent and both stepwise author protocols agree on 2X.",
        uncertainty="The published STAR Methods 5X wording has not been formally corrected.",
        slot=AnnotationSlot.INTRO))
    graph.annotate(SourceDiscrepancy(
        target="tagment",title="Nextera TD/ATM volume split",
        claims=(
            SourceClaim("Hughes et al. 2020 STAR Methods; PMC JATS line 509",
                        "One ng WTA is combined with 10 µL Buffer TD and 5 µL Buffer ATM."),
            SourceClaim("Hughes et al. 2020 Document S1, S3 Master Protocol p. 22; mmc1 line 3810",
                        "The detailed recipe uses 11 µL Nextera TD and 4 µL ATM."),
            SourceClaim("Hughes et al. 2019 author S3 Master Protocol p. 22, step 3",
                        "The earlier author protocol also uses 11 µL TD and 4 µL ATM.")),
        model_choice="The graph records Nextera XT tagmentation without a numeric buffer/transposase split.",
        evidence="Both source recipes total 20 µL with 5 µL cDNA, but assign different TD and ATM volumes.",
        uncertainty="The authors do not explain whether the split changed between runs or is a transcription error.",
        slot=AnnotationSlot.INTRO))
    graph.annotate(Highlight(
        target="library", title="Displayed read configuration and Nextera context",
        explanation="The displayed 20/8/50 Read 1/Index 1/Read 2 configuration is the June 2019 author guide p. 24 MiSeq/NextSeq option and the 2023 chapter section 3.14 option. The 2019 guide separately lists NovaSeq Read 2 50–80 cycles, optional Index 2, and bead-lot-dependent 21-cycle Read 1. Those variants are not assigned to this single-index construct. The inferred Nextera right arm uses the official vendor family context; the historical kit revision and N700 member remain unassigned.",
        slot=AnnotationSlot.AFTER_WORKFLOW))
    return graph.require_publishable()
