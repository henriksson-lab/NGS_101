"""Source-specific sci-ATAC-seq3 reaction graph (Supplementary Table S7)."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(ROOT / "sci-atac-seq__10.1126+science.aab1601" / "tools"))
import sciatac as S
import illumina as il
import nextera as nx
from sciatac import pool_indexed_nuclei
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def ligation_scene(which: str) -> Scene:
    if which == "N5":
        joined = [seg("head", S.N5_HEAD, "r2"), seg("barcode", "A" * 10, "cbc", placeholder=True), seg("tail", S.N5_TAIL, "r2"), seg("s5 first 8", nx.S5[:8], "s5"), seg("s5 rest", nx.S5[8:], "s5"), seg("ME", nx.ME, "me")]
        bridge = [seg("splint", revcomp(S.N5_TAIL + nx.S5[:8]), "r2")]
        pair = ("splint", "tail")
    else:
        joined = [seg("head", S.N7_HEAD, "r3"), seg("barcode", "B" * 10, "cbc", placeholder=True), seg("tail", S.N7_TAIL, "r3"), seg("s7 first 8", nx.S7[:8], "s7"), seg("s7 rest", nx.S7[8:], "s7"), seg("ME", nx.ME, "me")]
        bridge = [seg("splint", revcomp(S.N7_TAIL + nx.S7[:8]), "r3")]
        pair = ("splint", "tail")
    sc = Scene(); sc.strand("joined strand", joined); sc.anneal("splint", bridge, to="joined strand", pair=pair, unpaired=())
    sc.mark("joined strand", "tail", "ligation nick", through=("s5 first 8" if which == "N5" else "s7 first 8"))
    return sc


def ligation_reagent_scene(which: str) -> Scene:
    head, tail = ((S.N5_HEAD, S.N5_TAIL) if which == "N5"
                  else (S.N7_HEAD, S.N7_TAIL))
    adapter = [seg(f"{which} head", head),
               seg(f"{which} barcode", "N" * 10, "cbc", placeholder=True),
               seg(f"{which} tail", tail)]
    sc = Scene(); sc.strand("adapter", adapter, label=f"{which} ligation oligo")
    sc.ligation_end("adapter", f"{which} tail", kind="splinted",
                    polarity="3-prime", owner="adapter")
    sc.note("adapter", "end-specific splint aligns this 3′ end to the phosphorylated Tn5 end")
    return sc


def protocol_graph() -> ProtocolGraph:
    source = "https://doi.org/10.1126/science.aba7612"
    graph = ProtocolGraph(S.SCI3, authority=source)
    nuclei = graph.input(
        "nuclei", "accessible chromatin in nuclei", MaterialKind.BIOLOGICAL,
        state=Construct([seg("genomic DNA", "X" * 40, placeholder=True)]),
        primary=True)
    reviewed = ["nuclei"]
    plain_adapter = Construct([
        seg("s5", nx.S5, "s5"), seg("ME", nx.ME, "me"),
        seg("s7", nx.S7, "s7"), seg("ME-2", nx.ME, "me")])
    transposome = graph.input("plain_transposome", "plain Nextera transposome",
                              MaterialKind.ADAPTER, state=plain_adapter)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged_state = Construct([
        seg("s5", nx.S5, "s5"), seg("ME", nx.ME, "me"),
        seg("genomic insert", "X" * 30, placeholder=True),
        seg("ME'", nx.ME_RC, "me"), seg("s7'", nx.S7_RC, "s7")])
    tagged = graph.product("tagmented_fragments", "plain tagmented chromatin",
                           MaterialKind.MOLECULAR, state=tagged_state)
    tagment_waste = graph.product("tagment_waste", "non-library tagmentation products",
                                  MaterialKind.WASTE)
    graph.react("tagment", "Tagment nuclei with plain Nextera Tn5",
                StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE),
                 (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tagment_waste))
    graph.discard(tagment_waste, "does not retain the heterologous s5/s7 topology")
    pnk = graph.input("pnk", "T4 polynucleotide kinase", MaterialKind.CATALYST)
    phosphorylated = graph.product("phosphorylated_fragments",
                                   "5′-phosphorylated tagmented fragments",
                                   MaterialKind.MOLECULAR, state=tagged_state)
    graph.react("phosphorylate", "Phosphorylate transferred 5′ ends",
                StepKind.PHOSPHORYLATION,
                ((tagged, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)),
                (phosphorylated,))
    ligase = graph.input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
    n5_state = Construct([
        seg("N5 head", S.N5_HEAD), seg("N5 barcode", "N" * 10,
                                       placeholder=True), seg("N5 tail", S.N5_TAIL)])
    n5 = graph.input("n5_adapter", "N5 split-pool ligation oligo",
                     MaterialKind.ADAPTER, state=n5_state)
    first_ligated = graph.product("n5_ligated", "N5-barcoded s5 fragments",
                                  MaterialKind.MOLECULAR,
                                  state=ligation_scene("N5"))
    graph.react("n5_ligation", "Ligate the first-round N5 barcode through its splint",
                StepKind.LIGATION,
                ((phosphorylated, InputRole.TEMPLATE),
                 (n5, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (first_ligated,),
                geometry=AdapterLigationGeometry.known(
                    "splinted", ligation_reagent_scene("N5"),
                    statement="Table S7 prints the N5 oligo and its end-specific splint."))
    pooled1, pool1_reviewed = pool_indexed_nuclei(
        graph, first_ligated, prefix="n5_round", state=ligation_scene("N5"),
        label="N5-barcoded fragments from other first-round wells")
    n7_state = Construct([
        seg("N7 head", S.N7_HEAD), seg("N7 barcode", "N" * 10,
                                       placeholder=True), seg("N7 tail", S.N7_TAIL)])
    n7 = graph.input("n7_adapter", "N7 split-pool ligation oligo",
                     MaterialKind.ADAPTER, state=n7_state)
    final = S.sci3_library()
    second_ligated = graph.product("n7_ligated", "N5/N7-barcoded fragments",
                                   MaterialKind.MOLECULAR, state=final)
    graph.react("n7_ligation", "Ligate the second-round N7 barcode through its splint",
                StepKind.LIGATION,
                ((pooled1, InputRole.TEMPLATE),
                 (n7, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (second_ligated,),
                geometry=AdapterLigationGeometry.known(
                    "splinted", ligation_reagent_scene("N7"),
                    statement="Table S7 prints the N7 oligo and its end-specific splint."))
    pooled2, pool2_reviewed = pool_indexed_nuclei(
        graph, second_ligated, prefix="n7_round", state=final,
        label="N5/N7-barcoded fragments from other second-round wells")
    gap_polymerase = graph.input("gap_polymerase", "gap-fill DNA polymerase",
                                 MaterialKind.CATALYST)
    filled = graph.product("gap_filled_fragments", "gap-filled N5/N7 fragments",
                           MaterialKind.MOLECULAR, state=final)
    graph.react("gap_fill", "Fill tagmentation gaps", StepKind.EXTENSION,
                ((pooled2, InputRole.TEMPLATE),
                 (gap_polymerase, InputRole.CATALYST)), (filled,))
    p5_sequence = il.P5 + "N" * 10 + S.N5_HEAD
    p7_sequence = il.P7 + "N" * 10 + S.N7_HEAD
    p5_primer = PCRPrimer("P5/N5 index primer", p5_sequence,
                          anneal_nt=len(S.N5_HEAD))
    p7_primer = PCRPrimer("P7/N7 index primer", p7_sequence,
                          anneal_nt=len(S.N7_HEAD))
    p5 = graph.input("p5_index_primer", p5_primer.name, MaterialKind.OLIGO,
                     state=Construct([seg(p5_primer.name, p5_sequence,
                                          placeholder=True)]))
    p7 = graph.input("p7_index_primer", p7_primer.name, MaterialKind.OLIGO,
                     state=Construct([seg(p7_primer.name, p7_sequence,
                                          placeholder=True)]))
    polymerase = graph.input("pcr_polymerase", "PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("library", "sci-ATAC-seq3 sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("index_pcr", "Add the third-round P5/P7 indices by PCR",
                StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(
                    final, p5_primer, p7_primer,
                    statement="Table S7 prints both third-round PCR primers."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=S.SCI3_PRIMERS,
        acquired_roles=tuple(p.role for p in S.SCI3_PRIMERS)))
    reviewed.extend([
        "plain_transposome", "tn5", "tagment", "pnk", "phosphorylate",
        "dna_ligase", "n5_adapter", "n5_ligation", *pool1_reviewed,
        "n7_adapter", "n7_ligation", *pool2_reviewed, "gap_polymerase",
        "gap_fill", "p5_index_primer", "p7_index_primer", "pcr_polymerase",
        "index_pcr", "library",
    ])
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
