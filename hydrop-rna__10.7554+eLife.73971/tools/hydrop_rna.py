"""Molecular construct model for HyDrop-RNA (De Rop et al., 2022)."""
from __future__ import annotations

import illumina as il
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from adapter_ligation import reported_adapter_end
from endprep import dA_tailed_scene,repair_and_dA_tail
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind)

UMI_FEATURE = feature("umi", "umi", "random")
I5_FEATURE = feature("sample_index_i5", "sample_index", "fixed")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    if name == "UMI": kw.setdefault("feature", UMI_FEATURE)
    elif name in ("i5", "i5 sample index"): kw.setdefault("feature", I5_FEATURE)
    elif name in ("i7", "i7 reverse complement"): kw.setdefault("feature", I7_FEATURE)
    return Segment(name=name, top=top, tag=tag, **kw)


T7_PROMOTER = "TAATACGACTCACTATAGGG"
ACRYDITE_PRIMER = "T" * 7 + T7_PROMOTER + rt.SMART_HANDLE + "AC"
LINK1 = "CAGCTACTGC"
LINK2 = "CGAGTACCCT"
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"
TSO_P = rt.SMART_HANDLE
I7_OLIGO_INDEX = "CGCTCAGTTC"
I5_OLIGO_INDEX = "TCGTGGAGCG"
HYI7_SPACER = "CTGTCCGCGG"
HYI7_SITE = HYI7_SPACER + rt.SMART_HANDLE + "AC"
CUSTOM_INDEX1 = revcomp(HYI7_SITE)
CUSTOM_READ2 = HYI7_SITE


def barcode(name: str) -> Segment:
    part = name.replace("barcode ", "round ")
    return seg(name, "B" * 10, "cbc", placeholder=True,
               feature=feature(name.replace(" ", "_"), "cell_barcode", "combinatorial",
                               group="cell_barcode", part=part))


def bead_oligo() -> Construct:
    """Final RNA bead oligo; linkage and barcode order are fixed by construction."""
    return Construct([
        seg("poly(T) before promoter", "T" * 7),
        seg("T7 promoter", T7_PROMOTER),
        seg("SMART handle", rt.SMART_HANDLE, "tso"),
        seg("AC", "AC"),
        barcode("barcode 1"), seg("linker 1", LINK1),
        barcode("barcode 2"), seg("linker 2", LINK2),
        barcode("barcode 3"),
        seg("UMI", "N" * 8, "umi", placeholder=True),
        seg("poly(T)", "T" * 30),
    ], name="HyDrop-RNA bead oligo")


def tso_oligo() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
            seg("rGrGrG", "rGrGrG", placeholder=True)]


def hyi7_primer() -> list[Segment]:
    return [seg("P7", il.P7, "p7"), seg("i7", I7_OLIGO_INDEX, "cbc"),
            seg("spacer", HYI7_SPACER), seg("SMART handle + AC", rt.SMART_HANDLE + "AC", "tso")]


def hyi5_primer() -> list[Segment]:
    return [seg("P5", il.P5, "p5"), seg("i5", I5_OLIGO_INDEX, "cbc"),
            seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1")]


def rt_scene(insert_nt: int = 28) -> Scene:
    """Released bead primer annealed to a poly(A)-tailed transcript."""
    primer = list(bead_oligo())
    mrna = [seg("transcript", "X" * insert_nt, placeholder=True),
            seg("poly(A)", "A" * 30)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA", mod3="poly(A)")
    sc.anneal("bead primer", primer, to="mRNA", pair=("poly(T)", "poly(A)"),
              label="released bead primer")
    sc.arrow("bead primer", "reverse transcription and template switching")
    return sc


def amplified_cdna(insert_nt: int = 28) -> Construct:
    """TSO-P amplicon; upstream acrydite/T7 prefix is excluded by the priming site."""
    bead = bead_oligo()
    keep = list(bead)[2:]
    return Construct([
        *keep,
        seg("cDNA", "X" * insert_nt, placeholder=True),
        seg("template-switch junction", "CCC"),
        seg("copied TSO suffix", revcomp("GAAT")),
        seg("opposite SMART handle", revcomp(rt.SMART_HANDLE), "tso"),
    ], name="HyDrop-RNA amplified cDNA")


def final_library(insert_nt: int = 28) -> Construct:
    """P5-to-P7' bead-end library selected by HYi5/HYi7 PCR."""
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("i5 sample index", I5_OLIGO_INDEX, "cbc"),
        seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
        seg("cDNA", "X" * insert_nt, placeholder=True), seg("poly(A)", "A" * 30),
        seg("UMI", "N" * 8, "umi", placeholder=True), barcode("barcode 3"),
        seg("linker 2, plate orientation", revcomp(LINK2)), barcode("barcode 2"),
        seg("linker 1, plate orientation", revcomp(LINK1)), barcode("barcode 1"),
        seg("HYi7 site reverse complement", revcomp(HYI7_SITE), "tso"),
        seg("i7 reverse complement", revcomp(I7_OLIGO_INDEX), "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="HyDrop-RNA sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final HyDrop-RNA construct: " + "; ".join(problems))
    return lib


SEQ_PRIMERS = (
    sp.TRUSEQ["R1"],
    sp.custom("Index 1 (i7)", "HyDrop custom Index 1", CUSTOM_INDEX1,
              "De Rop et al. 2022, Supplementary file 3"),
    sp.TRUSEQ["I2"],
    sp.custom("Read 2", "HyDrop custom Read 2", CUSTOM_READ2,
              "De Rop et al. 2022, Supplementary file 3"),
)


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def read2_layout() -> list[tuple[str, str]]:
    spans=sp.feature_spans(final_library(),SEQ_PRIMERS,{"Read 2":58})
    by_start={x.cycle_start:(x.cycles,x.feature.label+(
              f" ({x.feature.part})" if x.feature.part else "")) for x in spans}
    return [by_start[1],("11–20","linker 1"),by_start[21],
            ("31–40","linker 2"),by_start[41],by_start[51]]


def first_strand(insert_nt: int = 28) -> Construct:
    return Construct([*bead_oligo(),seg("cDNA","X"*insert_nt,placeholder=True),
                      seg("untemplated C","CCC")],name="HyDrop-RNA first-strand cDNA")


def protocol_graph() -> ProtocolGraph:
    source="Defining source: https://doi.org/10.7554/eLife.73971"
    graph=ProtocolGraph("HyDrop-RNA",authority=source)
    precursor=graph.input("bead_precursor","acrydite hydrogel bead primer",
                          MaterialKind.ADAPTER)
    barcode_plates=graph.input("barcode_plates","three 96-way barcode extension plates",
                               MaterialKind.OLIGO)
    bead_pol=graph.input("bead_polymerase","bead-extension polymerase",MaterialKind.CATALYST)
    bead=graph.product("barcoded_bead",bead_oligo().name,MaterialKind.ADAPTER,state=bead_oligo())
    graph.react("bead_synthesis","Build the three-round split-pool-barcoded RNA bead oligo",
                StepKind.EXTENSION,
                ((precursor,InputRole.TEMPLATE),(barcode_plates,InputRole.STRUCTURAL),
                 (bead_pol,InputRole.CATALYST)),(bead,))
    cells=graph.input("cells","single-cell suspension",MaterialKind.BIOLOGICAL,primary=True)
    droplets=graph.input("droplet_device","HyDrop droplet generator",MaterialKind.CONDITION)
    productive=graph.product("productive_droplets","droplets containing cell and barcoded bead",
                              MaterialKind.MOLECULAR)
    empty=graph.product("unproductive_droplets","droplets without a usable cell/bead pair",
                        MaterialKind.WASTE)
    graph.react("partition","Co-encapsulate cells with barcoded hydrogel beads",
                StepKind.PARTITION,
                ((cells,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
                 (droplets,InputRole.CONDITION)),(productive,empty))
    graph.discard(empty,"unproductive droplets do not contribute indexed library molecules")
    dtt=graph.input("dtt","DTT hydrogel-solubilization condition",MaterialKind.CONDITION)
    released=graph.product("released_capture","released bead primers with cellular RNA",
                           MaterialKind.MOLECULAR,state=rt_scene())
    matrix=graph.product("bead_matrix","insoluble hydrogel bead remnants",MaterialKind.WASTE)
    graph.react("release_bead","Dissolve the hydrogel bead and release capture primers",
                StepKind.SOLUBILIZATION,
                ((productive,InputRole.TEMPLATE),(dtt,InputRole.CONDITION)),(released,matrix))
    graph.discard(matrix,"the hydrogel matrix is not part of the library")
    reverse_transcriptase=graph.input("reverse_transcriptase","reverse transcriptase",
                                     MaterialKind.CATALYST)
    first=graph.product("first_strand",first_strand().name,MaterialKind.MOLECULAR,
                        state=first_strand())
    graph.react("reverse_transcribe","Reverse-transcribe poly(A) RNA from released bead primers",
                StepKind.REVERSE_TRANSCRIPTION,
                ((released,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
                 (reverse_transcriptase,InputRole.CATALYST)),(first,))
    tso=graph.input("tso","template-switching oligo",MaterialKind.OLIGO,
                    state=tuple(tso_oligo()))
    switched=graph.product("switched_cdna","template-switched HyDrop cDNA",
                           MaterialKind.MOLECULAR,state=amplified_cdna())
    graph.react("template_switch","Install the distal SMART handle",
                StepKind.TEMPLATE_SWITCHING,
                ((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
                 (reverse_transcriptase,InputRole.CATALYST)),(switched,))
    exonuclease=graph.input("exonuclease","Exonuclease I",MaterialKind.CATALYST)
    cleaned=graph.product("cleaned_cdna","template-switched cDNA after primer cleanup",
                          MaterialKind.MOLECULAR,state=amplified_cdna())
    oligo_waste=graph.product("oligo_waste","digested unused bead primers",MaterialKind.WASTE)
    graph.react("cleanup","Remove unused bead primers",StepKind.CLEANUP,
                ((switched,InputRole.TEMPLATE),(exonuclease,InputRole.CATALYST)),
                (cleaned,oligo_waste))
    graph.discard(oligo_waste,"unused capture primers are removed before cDNA amplification")
    tso_p=graph.input("tso_p","TSO-P cDNA-PCR primer",MaterialKind.OLIGO)
    preamp_pol=graph.input("preamp_polymerase","cDNA PCR polymerase",MaterialKind.CATALYST)
    amplified=graph.product("amplified_cdna",amplified_cdna().name,MaterialKind.MOLECULAR,
                            state=amplified_cdna())
    graph.react("preamp","Amplify SMART-handle-flanked cDNA with TSO-P",StepKind.PCR,
                ((cleaned,InputRole.TEMPLATE),(tso_p,InputRole.STRUCTURAL),
                 (preamp_pol,InputRole.CATALYST)),(amplified,),
                geometry=PCRGeometry.symmetric(
                    amplified_cdna(),PCRPrimer("TSO-P",TSO_P),
                    statement="Supplementary file 3 prints TSO-P and the copied SMART handles bound at both cDNA ends."))
    fragmentase=graph.input("fragmentase","library fragmentation enzyme mix",
                            MaterialKind.CATALYST)
    fragments=graph.product("fragments","fragmented HyDrop cDNA",MaterialKind.MOLECULAR)
    graph.react("fragment","Fragment amplified cDNA",StepKind.FRAGMENTATION,
                ((amplified,InputRole.TEMPLATE),(fragmentase,InputRole.CATALYST)),(fragments,))
    endprep=graph.input("endprep_mix","end-repair and dA-tailing mix",MaterialKind.CATALYST)
    prepared_model=repair_and_dA_tail(Construct([
        seg("HyDrop cDNA fragment","X"*36,placeholder=True)],name="HyDrop cDNA fragment"))
    prepared=graph.product("prepared_fragments","end-repaired, dA-tailed cDNA fragments",
                           MaterialKind.MOLECULAR,
                           state=dA_tailed_scene(prepared_model,label="ligation-ready HyDrop cDNA"))
    graph.react("end_prepare","End-repair and dA-tail fragmented cDNA",StepKind.END_REPAIR,
                ((fragments,InputRole.TEMPLATE),(endprep,InputRole.CATALYST)),(prepared,))
    adapter=graph.input("nebnext_hairpin","NEBNext hairpin adapter",MaterialKind.ADAPTER)
    ligase=graph.input("ligase","DNA ligase",MaterialKind.CATALYST)
    ligated=graph.product("hairpin_ligated","NEBNext-hairpin-ligated HyDrop fragments",
                          MaterialKind.MOLECULAR)
    graph.react("adapter_ligation","Ligate the NEBNext hairpin adapter",StepKind.LIGATION,
                ((prepared,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(ligated,),
                geometry=reported_adapter_end(
                    "NEBNext hairpin",kind="cohesive",polarity="3-prime",
                    terminal_sequence="T",
                    statement="The NEBNext hairpin presents a terminal 3′ T to dA-tailed fragments."))
    user=graph.input("user_enzyme","USER enzyme",MaterialKind.CATALYST)
    opened=graph.product("opened_library","USER-opened adapter-ligated fragments",
                         MaterialKind.MOLECULAR,state=final_library())
    loop=graph.product("released_hairpin","released hairpin loop",MaterialKind.WASTE)
    graph.react("user_digest","Open the dU-containing hairpin",StepKind.DIGESTION,
                ((ligated,InputRole.TEMPLATE),(user,InputRole.CATALYST)),(opened,loop))
    graph.discard(loop,"the cleaved hairpin loop is not part of the library")
    hyi5=graph.input("hyi5","HYi5 TruSeq library-PCR primer",MaterialKind.OLIGO,
                     state=tuple(hyi5_primer()))
    hyi7=graph.input("hyi7","HYi7 bead-end library-PCR primer",MaterialKind.OLIGO,
                     state=tuple(hyi7_primer()))
    library_pol=graph.input("library_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    graph.react("library_pcr","Select bead-end fragments and add indexed P5/P7 arms",
                StepKind.PCR,
                ((opened,InputRole.TEMPLATE),(hyi5,InputRole.STRUCTURAL),
                 (hyi7,InputRole.STRUCTURAL),(library_pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.exact(
                    final,PCRPrimer("HYi5","".join(s.top for s in hyi5_primer())),
                    PCRPrimer("HYi7","".join(s.top for s in hyi7_primer()),
                              anneal_nt=len(HYI7_SITE)),
                    statement="Supplementary file 3 prints HYi5/HYi7 and their TruSeq and bead-backbone landing regions."))
    lengths={"Read 1":50,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":58}
    graph.sequence(library,SequencingPlan(
        library=final,primers=SEQ_PRIMERS,acquired_roles=tuple(lengths),read_lengths=lengths,
        read_uses={"Read 1":"include","Index 1 (i7)":"demultiplex",
                   "Index 2 (i5)":"demultiplex","Read 2":"extract"}))
    graph.declare_reviewed_scope((
        "bead_precursor","barcode_plates","bead_polymerase","bead_synthesis","cells",
        "droplet_device","partition","dtt","release_bead","reverse_transcriptase",
        "reverse_transcribe","tso","template_switch","exonuclease","cleanup","tso_p",
        "preamp_polymerase","preamp","fragmentase","fragment","endprep_mix","end_prepare",
        "nebnext_hairpin","ligase","adapter_ligation","user_enzyme","user_digest","hyi5",
        "hyi7","library_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
