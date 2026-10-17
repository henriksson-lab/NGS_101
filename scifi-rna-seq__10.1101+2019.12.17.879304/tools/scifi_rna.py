"""Molecular construct model for scifi-RNA-seq (Datlinger et al., Nature Methods 2021)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (AnnotationSlot, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)

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


def rt_primer() -> list[Segment]:
    return [seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
            seg("UMI", "N" * 8, "umi", placeholder=True), seg("fixed N", "A"),
            seg("round-1 barcode", "B" * 11, "cbc", placeholder=True),
            seg("fixed V", "A"), seg("poly(T)", "T" * 30),
            seg("V", "V", placeholder=True), seg("N", "N")]


# The gel-bead oligo appears only in the published article's Extended Data Fig. 3a
# (article PDF page 16, an image-only page): 5'-AATGATACGGCGACCACCGAGATCTACAC
# [10x_round2]TCGTCGGCAGCGTC, read from the rendered page. The page's OCR twin garbles
# the s5 part (e.g. "ICGTCGGCAGCOTC"), so the text guard cannot match it.
SOURCE_OLIGO_EXCEPTIONS = {
    "bead_oligo": ("Printed only as an image: published Extended Data Fig. 3a "
                   "(s41592-021-01153-z.pdf page 16, 'Gel Bead Oligo Primer with 737k "
                   "possible microfluidic barcodes (16b)') shows 5'-AATGATACGGCGACCACC"
                   "GAGATCTACAC[10x_round2]TCGTCGGCAGCGTC; verified on the rendered page "
                   "because the OCR text twin misreads the s5 bases."),
}


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


def tso_oligo() -> list[Segment]:
    """Supplementary Table 1 TSO; its printed rGrGrG kept as three RNA G residues."""
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
            seg("rGrGrG", "/rG//rG//rG/", placeholder=True, length_bp=3)]


def indexed_p7_primer() -> Construct:
    """Supplementary Table 1 P7-001...P7-016: P7 + 8-nt i7 window + s7."""
    return Construct([seg("P7", il.P7, "p7"),
                      seg("i7 index", "N" * len(P7_INDEX_READ), "cbc", placeholder=True,
                          feature=I7_FEATURE),
                      seg("s7", nx.S7, "s7")], name="scifi-RNA indexed P7 primer")


def i7_transposome_adapter() -> Scene:
    """Supplementary Table 1 Tn5-bottom_Read2N annealed to [5Phos] Tn5-top_ME."""
    sc = Scene()
    sc.strand("Read2N", [seg("s7", nx.S7, "s7"), seg("mosaic end", nx.ME, "me")],
              label="Tn5-bottom_Read2N")
    sc.anneal("ME", [seg("mosaic end reverse complement", nx.ME_RC, "me")], to="Read2N",
              pair=("mosaic end reverse complement", "mosaic end"), label="Tn5-top_ME",
              mod5="p")
    return sc


def tso_scene(insert_nt: int = 28) -> Scene:
    first = ligated_first_strand(insert_nt)
    tso = tso_oligo()
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
        seg("V", "V", placeholder=True), seg("N", "N"),
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
    source=("Defining source: Datlinger et al., Nature Methods 18:635 (2021), "
            "https://doi.org/10.1038/s41592-021-01153-z, Supplementary Table 1 and "
            "Supplementary Protocol v2021-01-25; gel-bead oligo drawn only in the "
            "article's Extended Data Fig. 3a (image)")
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
                    statement="Published Supplementary Table 1 prints the 3'-ddC bridge and the 5'-phosphorylated RT primers; the gel-bead oligo (P5 + 16-nt round-2 barcode + TCGTCGGCAGCGTC) is drawn only in the published article's Extended Data Fig. 3a, read from the rendered figure."))
    oligos="Datlinger et al. 2021 Supplementary Table 1 (scifi-RNA-seq oligos sheet)"
    tso=graph.input("tso","template-switching oligo",MaterialKind.OLIGO,source=oligos,
                    state=Construct(tso_oligo(),name="scifi-RNA TSO"))
    switched=graph.product("switched_cdna","template-switched two-level-barcoded cDNA",
                           MaterialKind.MOLECULAR,state=ligated_first_strand())
    graph.react("template_switch","Install the distal SMART amplification handle",
                StepKind.TEMPLATE_SWITCHING,
                ((ligated,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
                 (reverse_transcriptase,InputRole.CATALYST)),(switched,))
    preamp_forward=graph.input("preamp_forward","bead-end cDNA enrichment primer (Partial-P5)",
                               MaterialKind.OLIGO,source=oligos,
                               state=Construct([seg("Partial-P5",PARTIAL_P5,"p5")],
                                               name="Partial-P5"))
    preamp_reverse=graph.input("preamp_reverse","TSO-handle cDNA enrichment primer",
                               MaterialKind.OLIGO,source=oligos,
                               state=Construct([seg("TSO_enrichment_primer",rt.SMART_HANDLE,
                                                    "tso")],name="TSO_enrichment_primer"))
    preamp_pol=graph.input("preamp_polymerase","cDNA enrichment polymerase",
                           MaterialKind.CATALYST)
    amplified=graph.product("amplified_cdna","bulk-enriched full cDNA",MaterialKind.MOLECULAR,
                            state=ligated_first_strand())
    graph.react("preamp","Enrich full template-switched cDNA in bulk",StepKind.PCR,
                ((switched,InputRole.TEMPLATE),(preamp_forward,InputRole.STRUCTURAL),
                 (preamp_reverse,InputRole.STRUCTURAL),(preamp_pol,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    ligated_first_strand(),reason="The page model does not establish both bulk-enrichment primer landing scenes."))
    transposome=graph.input("i7_transposome","i7-only Tn5 transposome",MaterialKind.ADAPTER,
                            source=oligos,state=i7_transposome_adapter())
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    selected=selected_tagmented_strand()
    tagged=graph.product("selected_tagmented",selected.name,MaterialKind.MOLECULAR,state=selected)
    tag_waste=graph.product("tagment_waste","Tn5/Tn5 and TSO/Tn5 products lacking both PCR sites",
                            MaterialKind.WASTE)
    graph.react("tagment","Tagment enriched cDNA with i7-only Tn5",StepKind.TAGMENTATION,
                ((amplified,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,tag_waste))
    graph.discard(tag_waste,"fragments lacking the bead end and one s7 end are not PCR-selected")
    partial_p5=graph.input("partial_p5","partial-P5 selective primer",MaterialKind.OLIGO,
                           source=oligos,state=Construct([seg("Partial-P5",PARTIAL_P5,"p5")],
                                                         name="Partial-P5"))
    indexed_p7=graph.input("indexed_p7","indexed P7/s7 primer",MaterialKind.OLIGO,
                           source=oligos,state=indexed_p7_primer())
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
                    statement="Published Supplementary Table 1 and the Supplementary Protocol library-enrichment step define Partial-P5 and the indexed P7/s7 primers."))
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
    graph.annotate(SourceDiscrepancy(
        target="library_pcr", title="P7 index orientation differs within the published oligo table",
        claims=(
            SourceClaim("Datlinger et al. 2021 Supplementary Table 1 (oligos sheet), P7-001–P7-008",
                        "The primer contains the reverse complement of the listed P7_index_seq; P7-001 has TCGCCTTA for listed TAAGGCGA."),
            SourceClaim("The same table, P7-009–P7-016",
                        "The primer contains the listed bases as written; P7-009 has AAGAGGCA for listed N711 AAGAGGCA.")),
        model_choice="The representative library uses P7-001 and its printed primer sequence; it does not silently normalize P7-009–P7-016.",
        evidence="Direct comparison of each printed full primer and its index column shows the orientation change.",
        uncertainty="It is unknown whether the last eight primers were synthesized as printed or transcribed incorrectly.",
        slot=AnnotationSlot.INTRO))
    graph.annotate(SourceDiscrepancy(
        target="barcode_ligation", title="Thermoligation denaturation temperature",
        claims=(
            SourceClaim("Datlinger et al. 2021 Supplementary Protocol v2021-01-25, thermoligation (scATAC v1.0 and v1.1 chips)",
                        "Twelve cycles of 95 °C for 30 s and 59 °C for 2 min."),
            SourceClaim("Datlinger et al. 2021 Methods, 'Microfluidic thermoligation barcoding'",
                        "Twelve cycles of 98 °C for 30 s and 59 °C for 2 min.")),
        model_choice="The graph records thermoligation without fixing a temperature; the note presents both recipes.",
        evidence="Both parts of the published article describe the same round-2 bead barcode ligation with different denaturation settings.",
        uncertainty="Whether this was a deliberate protocol revision is not established.",
        slot=AnnotationSlot.INTRO))
    return graph.require_publishable()
