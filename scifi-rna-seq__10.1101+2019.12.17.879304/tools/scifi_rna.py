"""Molecular construct model for scifi-RNA-seq (Datlinger et al.)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind)

UMI_FEATURE = feature("umi", "umi", "random")
ROUND1_FEATURE = feature("cell_barcode_round1", "cell_barcode", "combinatorial",
                         group="cell_barcode", part="RT plate")
ROUND2_FEATURE = feature("cell_barcode_round2", "cell_barcode", "combinatorial",
                         group="cell_barcode", part="bead")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")
FEATURES = {"UMI": UMI_FEATURE, "round-1 barcode": ROUND1_FEATURE,
            "round-2 bead barcode": ROUND2_FEATURE, "i7 index read": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


RT_BARCODE_EXAMPLE = "AGTGATTAGCA"
P7_INDEX_READ = "TAAGGCGA"
BRIDGE = "CGTCGTGTAGGGAAAGAGTGTGACGCTGCCGACGA"
PARTIAL_P5 = il.P5[:22]
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"


def rt_primer() -> list[Segment]:
    return [seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
            seg("UMI", "N" * 8, "umi", placeholder=True), seg("fixed N", "A"),
            seg("round-1 barcode", "B" * 11, "cbc", placeholder=True),
            seg("fixed V", "A"), seg("poly(T)", "T" * 30),
            seg("V", "V", placeholder=True), seg("N", "N", placeholder=True)]


def bead_oligo() -> list[Segment]:
    return [seg("P5", il.P5, "p5"),
            seg("round-2 bead barcode", "C" * 16, "cbc", placeholder=True),
            seg("s5", nx.S5, "s5")]


def bridge_oligo() -> list[Segment]:
    return [seg("RT-primer bridge", BRIDGE[:21]), seg("bead bridge", BRIDGE[21:])]


def ligation_scene() -> Scene:
    """Bridge-enforced nick ligation between bead and phosphorylated RT oligos."""
    primer = rt_primer()
    joined = [*bead_oligo(),
              seg("Read-1 bridge site", il.TRUSEQ_READ1[:21], "r1"),
              seg("Read-1 remainder", il.TRUSEQ_READ1[21:], "r1"), *primer[1:]]
    sc = Scene()
    sc.strand("product", joined, label="bead + RT primer")
    sc.anneal("bridge", bridge_oligo(), to="product",
              pair=("RT-primer bridge", "Read-1 bridge site"), label="3'-ddC bridge",
              above=True)
    # The named pair fixes the 21-nt register; the adjacent 14 nt necessarily cover s5
    # because BRIDGE is the reverse complement of s5 + the first 21 nt of Read 1.
    if revcomp(BRIDGE) != nx.S5 + il.TRUSEQ_READ1[:21]:
        raise ValueError("bridge no longer joins s5 directly to TruSeq Read 1")
    sc.junction("product", "s5", "Read-1 bridge site")
    return sc


def ligation_ready_scene() -> Scene:
    sc=Scene();sc.strand("bead oligo",bead_oligo(),label="gel-bead oligo")
    sc.strand("RT primer",rt_primer(),label="phosphorylated RT primer",mod5="p")
    sc.strand("bridge",bridge_oligo(),label="3-prime-ddC bridge",mod3="ddC")
    sc.ligation_end("bead oligo","s5",kind="splinted",polarity="3-prime",
                    owner="adapter",text="bead s5 3′-OH at the nick")
    sc.ligation_end("RT primer","TruSeq Read 1",kind="splinted",polarity="5-prime",
                    owner="substrate",text="RT-primer 5′ phosphate at the nick")
    return sc


def ligated_first_strand(insert_nt: int = 28) -> Construct:
    return Construct([
        *bead_oligo(), *rt_primer(),
        seg("cDNA", "X" * insert_nt, placeholder=True),
        seg("untemplated C", "CCC"),
    ], name="scifi-RNA ligated first strand")


def rt_first_strand(insert_nt: int = 28) -> Construct:
    return Construct([*rt_primer(),seg("cDNA","X"*insert_nt,placeholder=True),
                      seg("untemplated C","CCC")],name="round-1-barcoded first strand")


def tso_scene(insert_nt: int = 28) -> Scene:
    first = ligated_first_strand(insert_nt)
    tso = [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
           seg("rGrGrG", "rGrGrG", placeholder=True)]
    sc = Scene()
    sc.strand("first", list(first), label="first-strand cDNA")
    sc.anneal("TSO", tso, to="first", pair=("rGrGrG", "untemplated C"), label="TSO")
    sc.arrow("TSO", "template switching; then TSO-primer enrichment")
    return sc


INDEX2_FROM_P5 = sp.custom(
    "Index 2 (i5)", "flow-cell P5 oligo (forward-strand workflow)", il.P5,
    "Illumina forward-strand indexed sequencing",
    "The bead barcode occupies the i5 position immediately after P5.")
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"], INDEX2_FROM_P5, sp.NEXTERA["R2"])


def final_library(insert_nt: int = 28) -> Construct:
    """Selected bead-end/i7-Tn5 product after partial-P5/P7 PCR."""
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("round-2 bead barcode", "C" * 16, "cbc", placeholder=True),
        seg("s5", nx.S5, "s5"), seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
        seg("UMI", "N" * 8, "umi", placeholder=True), seg("fixed N", "A"),
        seg("round-1 barcode", "B" * 11, "cbc", placeholder=True),
        seg("fixed V", "A"), seg("poly(T)", "T" * 30),
        seg("V", "V", placeholder=True), seg("N", "N", placeholder=True),
        seg("cDNA", "X" * insert_nt, placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 index read", P7_INDEX_READ, "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scifi-RNA-seq sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final scifi-RNA-seq construct: " + "; ".join(problems))
    return lib


def selected_tagmented_strand(insert_nt: int = 28) -> Construct:
    return Construct([
        *bead_oligo(),*rt_primer(),seg("cDNA","X"*insert_nt,placeholder=True),
        seg("mosaic end reverse complement",nx.ME_RC,"me"),
        seg("s7 reverse complement",nx.S7_RC,"s7")],
        name="scifi-RNA bead-end/i7-Tn5 fragment")


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def read1_layout() -> list[tuple[str, str]]:
    spans=sp.feature_spans(final_library(),SEQ_PRIMERS,{"Read 1":21})
    by_id={x.feature.id:(x.cycles,x.feature.label) for x in spans}
    return [by_id["umi"],("9","fixed N"),by_id["cell_barcode_round1"],
            ("21","fixed V")]


def protocol_graph() -> ProtocolGraph:
    source=("Defining source: https://doi.org/10.1101/2019.12.17.879304; "
            "expanded protocol: https://doi.org/10.1038/s41592-021-01153-z")
    graph=ProtocolGraph("scifi-RNA-seq",authority=source)
    cells=graph.input("fixed_cells","fixed permeabilized cells",MaterialKind.BIOLOGICAL,
                      primary=True)
    primer=graph.input("rt_primer","round-1 barcoded RT primer",MaterialKind.OLIGO,
                       state=tuple(rt_primer()))
    reverse_transcriptase=graph.input("reverse_transcriptase","reverse transcriptase",
                                     MaterialKind.CATALYST)
    first=graph.product("first_strand",rt_first_strand().name,MaterialKind.MOLECULAR,
                        state=rt_first_strand())
    graph.react("reverse_transcribe","Install the plate barcode and UMI by reverse transcription",
                StepKind.REVERSE_TRANSCRIPTION,
                ((cells,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),
                 (reverse_transcriptase,InputRole.CATALYST)),(first,))
    peer=graph.input("peer_wells","first-strand cDNA from the other RT wells",
                     MaterialKind.MOLECULAR)
    pooled=graph.product("pooled_cdna","pooled round-1-barcoded cells and cDNA",
                         MaterialKind.MOLECULAR,state=rt_first_strand())
    graph.react("pool","Pool the 384 RT wells",StepKind.POOL,
                ((first,InputRole.TEMPLATE),(peer,InputRole.TEMPLATE)),(pooled,),
                changes_molecule=False)
    droplets=graph.input("droplet_device","overloaded 10x scATAC droplet channel",
                         MaterialKind.CONDITION)
    partitioned=graph.product("productive_droplets","droplets containing round-1 cDNA and bead oligo",
                              MaterialKind.MOLECULAR,state=rt_first_strand())
    empty=graph.product("unproductive_droplets","droplets without usable cell/bead combinations",
                        MaterialKind.WASTE)
    graph.react("partition","Overload pooled cells into the bead-barcoding droplets",
                StepKind.PARTITION,
                ((pooled,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(partitioned,empty))
    graph.discard(empty,"unproductive droplets do not contribute indexed library molecules")
    bead=graph.input("bead_oligo","10x scATAC gel-bead oligo",MaterialKind.ADAPTER,
                     state=tuple(bead_oligo()))
    bridge=graph.input("bridge","3-prime-ddC bridge oligo",MaterialKind.OLIGO,
                       state=tuple(bridge_oligo()))
    ligase=graph.input("ampligase","thermostable DNA ligase",MaterialKind.CATALYST)
    ligated=graph.product("ligated_first_strand",ligated_first_strand().name,
                          MaterialKind.MOLECULAR,state=ligated_first_strand())
    graph.react("barcode_ligation","Ligate the bead barcode to the RT-primer end",
                StepKind.LIGATION,
                ((partitioned,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
                 (bridge,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(ligated,),
                geometry=AdapterLigationGeometry.known(
                    "splinted",ligation_ready_scene(),ligation_scene(),
                    statement="The expanded protocol prints the bridge, bead and phosphorylated RT-primer sequences that align the nick."))
    tso=graph.input("tso","template-switching oligo",MaterialKind.OLIGO)
    switched=graph.product("switched_cdna","template-switched two-level-barcoded cDNA",
                           MaterialKind.MOLECULAR,state=ligated_first_strand())
    graph.react("template_switch","Install the distal SMART amplification handle",
                StepKind.TEMPLATE_SWITCHING,
                ((ligated,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
                 (reverse_transcriptase,InputRole.CATALYST)),(switched,))
    preamp_forward=graph.input("preamp_forward","bead-end cDNA enrichment primer",
                               MaterialKind.OLIGO)
    preamp_reverse=graph.input("preamp_reverse","TSO-handle cDNA enrichment primer",
                               MaterialKind.OLIGO)
    preamp_pol=graph.input("preamp_polymerase","cDNA enrichment polymerase",
                           MaterialKind.CATALYST)
    amplified=graph.product("amplified_cdna","bulk-enriched full cDNA",MaterialKind.MOLECULAR,
                            state=ligated_first_strand())
    graph.react("preamp","Enrich full template-switched cDNA in bulk",StepKind.PCR,
                ((switched,InputRole.TEMPLATE),(preamp_forward,InputRole.STRUCTURAL),
                 (preamp_reverse,InputRole.STRUCTURAL),(preamp_pol,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    ligated_first_strand(),reason="The page model does not establish both bulk-enrichment primer landing scenes."))
    transposome=graph.input("i7_transposome","i7-only Tn5 transposome",MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    selected=selected_tagmented_strand()
    tagged=graph.product("selected_tagmented",selected.name,MaterialKind.MOLECULAR,state=selected)
    tag_waste=graph.product("tagment_waste","Tn5/Tn5 and TSO/Tn5 products lacking both PCR sites",
                            MaterialKind.WASTE)
    graph.react("tagment","Tagment enriched cDNA with i7-only Tn5",StepKind.TAGMENTATION,
                ((amplified,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,tag_waste))
    graph.discard(tag_waste,"fragments lacking the bead end and one s7 end are not PCR-selected")
    partial_p5=graph.input("partial_p5","partial-P5 selective primer",MaterialKind.OLIGO)
    indexed_p7=graph.input("indexed_p7","indexed P7/s7 primer",MaterialKind.OLIGO)
    library_pol=graph.input("library_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    graph.react("library_pcr","Select bead-end/s7 fragments and add indexed P7",StepKind.PCR,
                ((tagged,InputRole.TEMPLATE),(partial_p5,InputRole.STRUCTURAL),
                 (indexed_p7,InputRole.STRUCTURAL),(library_pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.exact(
                    selected,PCRPrimer("partial-P5",PARTIAL_P5),
                    PCRPrimer("indexed P7",il.P7+revcomp(P7_INDEX_READ)+nx.S7,
                              anneal_nt=len(nx.S7)),
                    statement="The expanded protocol defines partial-P5 selection of the bead end and the indexed P7/s7 primer."))
    lengths={"Read 1":21,"Index 1 (i7)":8,"Index 2 (i5)":16,"Read 2":78}
    graph.sequence(library,SequencingPlan(
        library=final,primers=SEQ_PRIMERS,acquired_roles=tuple(lengths),read_lengths=lengths,
        read_uses={"Read 1":"extract","Index 1 (i7)":"demultiplex",
                   "Index 2 (i5)":"extract","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "fixed_cells","rt_primer","reverse_transcriptase","reverse_transcribe","peer_wells",
        "pool","droplet_device","partition","bead_oligo","bridge","ampligase",
        "barcode_ligation","tso","template_switch","preamp_forward","preamp_reverse",
        "preamp_polymerase","preamp","i7_transposome","tn5","tagment","partial_p5",
        "indexed_p7","library_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
