"""Molecular model for Drop-seq.

Provenance.  Every library oligo is transcribed from the defining paper's own
Supplemental Information (Macosko et al., Cell 2015, Elsevier
1-s2.0-S0092867415005498-mmc1.pdf, Table S6 p. 29 and Supplemental Experimental
Procedures pp. 3-5).  Table S6 prints two bead batches:

    Barcoded Bead SeqA  5'-Bead-Linker-TTTTTTT + SMART handle + ACGT + J12 + N8 + T30
    Barcoded Bead SeqB  5'-Bead-Linker-TTTTTTT + SMART handle + AC   + J12 + N8 + T30
    Template_Switch_Oligo  SMART handle + GAAT + rGrGrG
    TSO_PCR                SMART handle
    P5-TSO_Hybrid          P5 + GCCTGTCCGCGG + SMART handle + *A*C
    Nextera_N701..N703     P7 + i7 + s7
    Read1CustomSeqA / B    GCCTGTCCGCGG + SMART handle + ACGT / AC

The Supplemental note "Bead Surface Primers and Custom Sequencing Primers" assigns SeqA
to the human-mouse experiments and retina replicates 1-3, SeqB to retina replicates 4-7,
and pairs each bead with its Read 1 primer.  The drawn library uses SeqB (the bead the
note says ChemGenes would manufacture); `final_library("SeqA")` builds the other batch.
The `SMART handle` is the canonical `rt.SMART_HANDLE`; P5, P7, s7 and the mosaic end are
canonical `illumina` / `nextera` constants.  Only the untemplated `CCC` added by the
reverse transcriptase is an inference (standard template-switching mechanism; its length
is not printed).

The scg_lib_structs drawn page is **not** a source here.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "random")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"cell barcode": CELL_BARCODE, "UMI": UMI_FEATURE,
            "i7 reverse complement": I7_FEATURE}
READ_LENGTHS = {"Read 1": 20, "Index 1 (i7)": 8, "Read 2": 50}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def inferred(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    """A base supplied by mechanism rather than printed by the source."""
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, inferred=True, **kw)


LINKER="T"*7
SPACER="GCCTGTCCGCGG"
# Cell 2015 Table S6: Barcoded Bead SeqA / SeqB and Read1CustomSeqA / SeqB.
BEAD_CONSTANT={"SeqA":"ACGT","SeqB":"AC"}
DRAWN_BEAD="SeqB"
N701_I7="TCGCCTTA"       # Table S6 Nextera_N701, between P7 and s7


def bead_primer(variant: str=DRAWN_BEAD) -> list[Segment]:
    return [seg("bead linker",LINKER),seg("SMART handle",rt.SMART_HANDLE,"tso"),
            seg("bead constant",BEAD_CONSTANT[variant]),
            seg("cell barcode","B"*12,"cbc",placeholder=True),
            seg("UMI","U"*8,"umi",placeholder=True),seg("dT30","T"*30)]


def tso() -> list[Segment]:
    return [seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("GAAT","GAAT"),
            seg("rGrGrG","GGG","tso")]


def tso_pcr() -> list[Segment]:
    return [seg("SMART handle",rt.SMART_HANDLE,"tso")]


def p5_hybrid() -> list[Segment]:
    """P5-TSO_Hybrid; its 3' AC (phosphorothioate-linked) matches both bead batches."""
    return [seg("P5",il.P5,"p5"),seg("spacer",SPACER,"r1"),
            seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("bead constant",BEAD_CONSTANT["SeqB"])]


def n701() -> list[Segment]:
    return [seg("P7",il.P7,"p7"),seg("N701 i7",N701_I7,"cbc"),seg("s7",nx.S7,"s7")]


def custom_r1(variant: str=DRAWN_BEAD) -> str:
    return SPACER+rt.SMART_HANDLE+BEAD_CONSTANT[variant]


def capture_scene() -> Scene:
    mrna=[seg("RNA body","X"*28,placeholder=True),seg("poly(A)","A"*30)]
    sc=Scene(); sc.strand("mRNA",mrna,label="captured mRNA")
    sc.anneal("bead primer",bead_primer(),to="mRNA",pair=("dT30","poly(A)"),
              label="barcoded bead",mod5="bead")
    sc.arrow("bead primer","reverse transcriptase"); return sc


def first_strand() -> Construct:
    return Construct([*bead_primer(),seg("antisense cDNA","X"*28,placeholder=True),
                      inferred("CCC","CCC","tso")],name="Drop-seq first strand")


def template_switch_scene() -> Scene:
    first=first_strand(); sc=Scene(); sc.strand("first strand",list(first),label="first strand",mod5="bead")
    sc.anneal("TSO",tso(),to="first strand",pair=("rGrGrG","CCC"),label="TSO",above=True)
    sc.arrow("first strand","RT switches template"); return sc


def wta_product() -> Construct:
    return Construct([*first_strand(),seg("copied GAAT",revcomp("GAAT")),
                      seg("copied SMART handle",revcomp(rt.SMART_HANDLE),"tso")],
                     name="Drop-seq template-switched cDNA")


def selected_tagmented_strand(insert_nt: int=34,variant: str=DRAWN_BEAD) -> Construct:
    return Construct([seg("SMART handle",rt.SMART_HANDLE,"tso"),
        seg("bead constant",BEAD_CONSTANT[variant]),seg("cell barcode","B"*12,"cbc",placeholder=True),
        seg("UMI","U"*8,"umi",placeholder=True),seg("poly(T)","T"*30),
        seg("antisense cDNA","X"*insert_nt,placeholder=True),seg("mosaic end reverse complement",nx.ME_RC,"me"),
        seg("s7 reverse complement",nx.S7_RC,"s7")],name="Drop-seq bead-end fragment")


def seq_primers(variant: str=DRAWN_BEAD):
    return (sp.custom("Read 1",f"Read1CustomSeq{variant[-1]}",custom_r1(variant),
                      f"Cell 2015 Table S6 Read1CustomSeq{variant[-1]}, paired with "
                      f"Barcoded Bead {variant} in the Supplemental note on bead primers"),
            sp.NEXTERA["I1"],sp.NEXTERA["R2"])


def final_library(insert_nt: int=34,variant: str=DRAWN_BEAD) -> Construct:
    selected=selected_tagmented_strand(insert_nt,variant)
    lib=Construct([seg("P5",il.P5,"p5"),seg("spacer",SPACER,"r1"),*list(selected),
                   seg("i7 reverse complement","I"*8,"cbc",placeholder=True),
                   seg("P7 reverse complement",il.P7_RC,"p7")],name=f"Drop-seq library ({variant} bead)")
    for p in seq_primers(variant):
        if sp.locate(lib,p) is None: raise ValueError(f"{p.name} does not land on {lib.name}")
    return lib


def primer_landings():
    lib=final_library(); out=[]
    for p in seq_primers():
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} does not land on {lib.name}")
        out.append((p,h))
    return out


def read1_layout() -> list[tuple[str,str]]:
    return [(x.cycles,x.feature.label) for x in
            sp.feature_spans(final_library(),seq_primers(),{"Read 1":20})]


SI=("Macosko et al. 2015 Cell Supplemental Information (Elsevier "
    "1-s2.0-S0092867415005498-mmc1.pdf)")
TS_INFERENCE=("the untemplated CCC added by the reverse transcriptase and paired by the "
              "TSO rGrGrG is the standard template-switching mechanism; its length is not "
              "printed in "+SI)


def protocol_graph() -> ProtocolGraph:
    source=("Defining source: https://doi.org/10.1016/j.cell.2015.05.002 and its "+SI+
            ": Table S6 oligos and Supplemental Experimental Procedures (Drop-Seq)")
    graph=ProtocolGraph("Drop-seq",authority=source)
    cells=graph.input("cells","single-cell suspension",MaterialKind.BIOLOGICAL,primary=True)
    beads=graph.input("barcoded_beads","Barcoded Bead SeqB in Drop-Seq Lysis Buffer + 50 mM DTT",
                      MaterialKind.ADAPTER,state=tuple(bead_primer()),source=source)
    droplets=graph.input("droplet_device","125 um co-flow device, ~1 nL droplets",MaterialKind.CONDITION)
    captured=graph.product("captured_rna","droplet-confined bead-captured poly(A) RNA",
                           MaterialKind.MOLECULAR,state=capture_scene(),source=source)
    empty=graph.product("empty_droplets","unproductive cell/bead loading outcomes",
                        MaterialKind.WASTE)
    graph.react("partition_capture","Co-encapsulate cells and barcoded beads; lyse and hybridize poly(A) RNA in droplets",StepKind.PARTITION,
                ((cells,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (droplets,InputRole.CONDITION)),(captured,empty),
                source=SI+" Droplet Generation and Droplet Breakage (perfluoro-1-octanol in 30 mL 6X SSC)")
    graph.discard(empty,"droplets without a productive cell-bead pairing do not yield library")
    reverse_transcriptase=graph.input("reverse_transcriptase","Maxima H- reverse transcriptase",
                                     MaterialKind.CATALYST)
    first=graph.product("first_strand",first_strand().name,MaterialKind.MOLECULAR,
                        state=first_strand(),inferred=True,source=TS_INFERENCE)
    graph.react("reverse_transcribe","INFERRED — Synthesize bead-barcoded first-strand cDNA (untemplated CCC by mechanism)",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (reverse_transcriptase,InputRole.CATALYST)),(first,),
                source=SI+" Reverse Transcription: RT mix with 4% Ficoll PM-400, 1 mM dNTPs, "
                       "10 U/uL Maxima H- RT; room temperature 30 min then 42 C 90 min")
    tso_input=graph.input("tso","Template_Switch_Oligo (2.5 uM)",MaterialKind.OLIGO,
                          state=tuple(tso()),source=source)
    switched=graph.product("template_switched",wta_product().name,MaterialKind.MOLECULAR,
                           state=wta_product(),inferred=True,source=TS_INFERENCE)
    graph.react("template_switch","INFERRED — Install the distal SMART handle by template switching",StepKind.TEMPLATE_SWITCHING,
                ((first,InputRole.TEMPLATE),(tso_input,InputRole.STRUCTURAL),
                 (reverse_transcriptase,InputRole.CATALYST)),(switched,))
    exo=graph.input("exonuclease_i","exonuclease I (1 U/uL)",MaterialKind.CATALYST)
    exo_beads=graph.product("exo_treated","bead-bound template-switched cDNA after exonuclease I",
                            MaterialKind.MOLECULAR,state=wta_product(),inferred=True,
                            source=TS_INFERENCE)
    exo_waste=graph.product("unextended_primer_digest",
                            "nucleotides from unextended single-stranded bead primers",
                            MaterialKind.WASTE)
    graph.react("exo_digest","Treat washed beads with exonuclease I to remove unextended primers",
                StepKind.DIGESTION,
                ((switched,InputRole.TEMPLATE),(exo,InputRole.CATALYST)),(exo_beads,exo_waste),
                source=(SI+" Exonuclease I Treatment: 1x Exonuclease I buffer, 1 U/uL "
                        "Exonuclease I, 37 C 45 min, then TE/SDS, TE/TW and water washes"))
    graph.discard(exo_waste,"the beads are washed after exonuclease I, before counting and PCR")
    smart=graph.input("smart_primer","TSO_PCR primer (0.8 uM)",MaterialKind.OLIGO,
                      state=tuple(tso_pcr()),source=source)
    preamp_pol=graph.input("preamp_polymerase","Kapa HiFi HotStart ReadyMix",
                           MaterialKind.CATALYST)
    amplified=graph.product("amplified_cdna","amplified full-length barcoded cDNA",
                            MaterialKind.MOLECULAR,state=wta_product(),inferred=True,
                            source=TS_INFERENCE)
    graph.react("preamp","INFERRED — Amplify template-switched cDNA from 1,000-bead aliquots (TSO junction CCC by mechanism)",StepKind.PCR,
                ((exo_beads,InputRole.TEMPLATE),(smart,InputRole.STRUCTURAL),
                 (preamp_pol,InputRole.CATALYST)),(amplified,),
                source=SI+" PCR: 4 cycles at 65 C annealing then X cycles at 67 C "
                       "(X = 8 human-mouse, 9 retina); 0.6x AMPure",
                geometry=PCRGeometry.symmetric(
                    wta_product(),PCRPrimer("TSO_PCR",rt.SMART_HANDLE),
                    statement="The molecule carries the same SMART handle at both amplification boundaries, so one primer amplifies it."))
    transposome=graph.input("transposome","Nextera XT transposome (600 pg cDNA per reaction)",MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    selected=selected_tagmented_strand()
    tagged=graph.product("selected_tagmented",selected.name,MaterialKind.MOLECULAR,
                         state=selected,source=source)
    waste=graph.product("tagment_waste","tagmentation products without the bead end",
                        MaterialKind.WASTE)
    graph.react("tagment","Tagment cDNA; the bead end nearest an s7 insertion is the selectable fragment",
                StepKind.TAGMENTATION,
                ((amplified,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,waste))
    graph.discard(waste,"fragments lacking the bead-derived handle are not PCR-selected")
    p5=graph.input("p5_primer","P5-TSO_Hybrid (200 nM)",MaterialKind.OLIGO,
                   state=tuple(p5_hybrid()),source=source)
    p7=graph.input("p7_primer","Nextera_N701 (200 nM; N702/N703 when multiplexed)",MaterialKind.OLIGO,
                   state=tuple(n701()),source=source)
    polymerase=graph.input("library_polymerase","Nextera XT PCR mix",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,
                          source=source)
    p5_sequence=il.P5+SPACER+rt.SMART_HANDLE+BEAD_CONSTANT["SeqB"]
    graph.react("library_pcr","Select the bead-end fragment and add P5/P7 (11 cycles)",StepKind.PCR,
                ((tagged,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
                 (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
                source=SI+" Preparation of Drop-Seq cDNA Library: custom P5_TSO_Hybrid and "
                       "Nextera_N701 replace the kit oligos; 95 C 30 s, 11 x (95 C 10 s, 55 C "
                       "30 s, 72 C 30 s), 72 C 5 min",
                geometry=PCRGeometry.exact(
                    selected,PCRPrimer("P5-TSO_Hybrid",p5_sequence,
                                       anneal_nt=len(rt.SMART_HANDLE)+len(BEAD_CONSTANT["SeqB"])),
                    PCRPrimer("Nextera N70X",il.P7+"N"*8+nx.S7,anneal_nt=len(nx.S7)),
                    statement="Cell 2015 Table S6 prints P5-TSO_Hybrid and Nextera_N701-N703; the hybrid's 3' AC matches the bead end, so only bead-end fragments with an s7 end amplify exponentially."))
    graph.sequence(library,SequencingPlan(
        library=final,primers=seq_primers(),acquired_roles=tuple(READ_LENGTHS),
        read_lengths=READ_LENGTHS,
        read_uses={"Read 1":"extract","Index 1 (i7)":"demultiplex","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "cells","barcoded_beads","droplet_device","partition_capture","reverse_transcriptase",
        "reverse_transcribe","tso","template_switch","exonuclease_i","exo_digest","smart_primer","preamp_polymerase",
        "preamp","transposome","tn5","tagment","p5_primer","p7_primer",
        "library_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
