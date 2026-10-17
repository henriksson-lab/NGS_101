"""Molecular construct model for sci-RNA-seq (Cao et al., 2017)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Highlight, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind, Terminal, Disposition)

UMI_FEATURE = feature("umi", "umi", "random")
RT_BARCODE = feature("cell_barcode_rt", "cell_barcode", "combinatorial",
                     group="cell_barcode", part="RT round")
I5_FEATURE = feature("cell_barcode_i5", "cell_barcode", "combinatorial",
                     group="cell_barcode", part="PCR well i5")
I7_FEATURE = feature("cell_barcode_i7", "cell_barcode", "combinatorial",
                     group="cell_barcode", part="PCR well i7")
FEATURES = {"UMI": UMI_FEATURE, "RT barcode": RT_BARCODE, "i5": I5_FEATURE,
            "i7": I7_FEATURE, "i7 reverse complement": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


HANDLE = il.TRUSEQ_READ1[-18:]
RT_TEMPLATE = HANDLE + "N" * 8 + "N" * 10 + "T" * 30 + "VN"
P5_TEMPLATE = il.P5 + "N" * 10 + il.TRUSEQ_READ1
P7_TEMPLATE = il.P7 + "N" * 10 + nx.S7


def rt_primer() -> list[Segment]:
    return [seg("Read 1 handle", HANDLE, "r1"),
            seg("UMI", "U" * 8, "umi", placeholder=True),
            seg("RT barcode", "B" * 10, "cbc", placeholder=True),
            seg("dT30", "T" * 30), seg("V", "V", placeholder=True),
            seg("N", "N", placeholder=True)]


def p5_primer() -> list[Segment]:
    return [seg("P5", il.P5, "p5"), seg("i5", "I" * 10, "cbc", placeholder=True),
            seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1")]


def p7_primer() -> list[Segment]:
    return [seg("P7", il.P7, "p7"), seg("i7", "J" * 10, "cbc", placeholder=True),
            seg("s7", nx.S7, "s7")]


def rt_scene() -> Scene:
    mrna = [seg("RNA body", "X" * 28, placeholder=True), seg("poly(A)", "A" * 30)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("RT primer", rt_primer(), to="mRNA", pair=("dT30", "poly(A)"),
              label="barcoded RT primer", unpaired=("V", "N"))
    sc.arrow("RT primer", "reverse transcriptase")
    return sc


def first_strand() -> Construct:
    return Construct([*rt_primer(), seg("antisense cDNA", "X" * 34, placeholder=True)],
                     name="sci-RNA-seq first strand")


def second_strand_scene() -> Scene:
    first = first_strand()
    sc = Scene.duplex(list(first), label="cDNA")
    sc.labels("top")
    return sc


def selected_tagmented_strand(insert_nt: int = 34) -> Construct:
    """RT-primer end through the nearest s7 insertion after gap fill."""
    return Construct([
        seg("Read 1 handle", HANDLE, "r1"),
        seg("UMI", "U" * 8, "umi", placeholder=True),
        seg("RT barcode", "B" * 10, "cbc", placeholder=True),
        seg("poly(T)", "T" * 30), seg("V", "V", placeholder=True),
        seg("N", "N", placeholder=True),
        seg("antisense cDNA", "X" * insert_nt, placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me", inferred=True),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
    ], name="selected sci-RNA-seq tagmented strand")


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"],
               sp.TRUSEQ["I2"], sp.NEXTERA["R2"])


def final_library(insert_nt: int = 34) -> Construct:
    selected = selected_tagmented_strand(insert_nt)
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5", "I" * 10, "cbc", placeholder=True),
        seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
        *list(selected)[1:],
        seg("i7 reverse complement", "J" * 10, "cbc", placeholder=True),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="sci-RNA-seq sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final sci-RNA-seq construct: " + "; ".join(problems))
    return lib


def primer_landings(lib: Construct | None = None):
    lib = lib or final_library()
    out = []
    for primer in SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        if hit is None:
            raise ValueError(f"{primer.name} does not land on {lib.name}")
        out.append((primer, hit))
    return out


def read1_layout() -> list[tuple[str, str]]:
    return [(x.cycles, x.feature.label) for x in
            sp.feature_spans(final_library(), SEQ_PRIMERS, {"Read 1": 18})]


def protocol_graph() -> ProtocolGraph:
    source="Cao2017 final main: https://doi.org/10.1126/science.aam8940; selected implementation: original preprint Methods https://doi.org/10.1101/104844"
    graph=ProtocolGraph("sci-RNA-seq",authority=source)
    cells=graph.input("fixed_cells","fixed permeabilized cells",MaterialKind.BIOLOGICAL,
                      primary=True)
    primer=graph.input("rt_primer","barcoded anchored oligo-dT RT primer",
                       MaterialKind.OLIGO,state=tuple(rt_primer()))
    rt=graph.input("reverse_transcriptase","SuperScript IV reverse transcriptase",MaterialKind.CATALYST)
    rt_stop=graph.input("rt_stop_solution","2X stop solution (40 mM EDTA, 1 mM spermidine)",
                        MaterialKind.CONDITION)
    first=graph.product("first_strand",first_strand().name,MaterialKind.MOLECULAR,
                        state=first_strand())
    graph.react("reverse_transcribe","Install the RT-well barcode and UMI in fixed cells (55 °C 10 min, then stop solution)",
                StepKind.REVERSE_TRANSCRIPTION,
                ((cells,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST),(rt_stop,InputRole.CONDITION)),(first,))
    sorter=graph.input("facs","FACS well assignment",MaterialKind.SELECTION)
    dapi=graph.input("dapi","DAPI stain (3 µM) for singlet gating",MaterialKind.CONDITION)
    sorted_cdna=graph.product("sorted_cdna","first-strand cDNA retained in sorted cells",
                              MaterialKind.MOLECULAR,state=first_strand())
    rejected=graph.product("sorting_waste","cells not retained for library construction",
                           MaterialKind.WASTE)
    graph.react("sort","Pool cells, stain with DAPI and sort them into PCR wells",StepKind.SAMPLING,
                ((first,InputRole.TEMPLATE),(dapi,InputRole.CONDITION),(sorter,InputRole.SELECTION)),
                (sorted_cdna,rejected),changes_molecule=False)
    graph.discard(rejected,"cells outside the selected PCR wells do not enter the library")
    second_pol=graph.input("second_strand_polymerase","NEB mRNA Second Strand Synthesis enzyme mix",
                           MaterialKind.CATALYST)
    duplex=graph.product("double_stranded_cdna","double-stranded barcoded cDNA",
                         MaterialKind.MOLECULAR,state=second_strand_scene())
    graph.react("second_strand","Synthesize the second cDNA strand (16 °C 150 min; 75 °C 20 min stop)",StepKind.EXTENSION,
                ((sorted_cdna,InputRole.TEMPLATE),(second_pol,InputRole.CATALYST)),(duplex,))
    carrier=graph.input("carrier_dna","human genomic DNA carrier",MaterialKind.MOLECULAR)
    transposome=graph.input("nextera_transposome","Nextera DNA Sample Preparation TDE1 transposome",MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("selected_tagmented","RT-handle-bearing tagmentation products",
                         MaterialKind.MOLECULAR)
    other=graph.product("other_tagmentation_products",
                        "carrier and other tagmentation fragments",
                        MaterialKind.MOLECULAR)
    graph.react("tagment","Tagment double-stranded cDNA with TDE1 and human genomic carrier",
                StepKind.TAGMENTATION,
                ((duplex,InputRole.TEMPLATE),(carrier,InputRole.SPIKE_IN),
                 (transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),
                (tagged,other))
    graph.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,
        "The preprint states that carrier DNA is not appreciably amplified or sequenced; "
        "this does not establish physical removal or exact structures of all other fragments."))
    binding=graph.input("dna_binding_buffer","Zymo DNA binding buffer",MaterialKind.CONDITION)
    beads=graph.input("cleanup_beads","AMPure XP beads",MaterialKind.SELECTION)
    cleaned=graph.product("cleaned_tagmented","purified tagmented DNA",MaterialKind.MOLECULAR)
    waste=graph.product("cleanup_waste","material not recovered by bead purification",MaterialKind.WASTE)
    graph.react("tagmentation_cleanup","Stop tagmentation with DNA binding buffer and purify on AMPure XP",
                StepKind.CLEANUP,
                ((tagged,InputRole.TEMPLATE),(binding,InputRole.CONDITION),(beads,InputRole.SELECTION)),
                (cleaned,waste),changes_molecule=False)
    graph.discard(waste,"Material not recovered by the published bead cleanup; its molecular composition is unassigned.")
    p5=graph.input("p5_primer","indexed P5 PCR primer",MaterialKind.OLIGO,
                   state=tuple(p5_primer()))
    p7=graph.input("p7_primer","indexed P7 PCR primer",MaterialKind.OLIGO,
                   state=tuple(p7_primer()))
    polymerase=graph.input("pcr_polymerase","NEBNext High-Fidelity 2X PCR Master Mix",MaterialKind.CATALYST)
    final=final_library()
    amplified=graph.product("pcr_products","PCR-well products",MaterialKind.MOLECULAR,state=final)
    geometry=PCRGeometry.exact(
        selected_tagmented_strand(),
        PCRPrimer("indexed P5 PCR primer",P5_TEMPLATE,anneal_nt=len(HANDLE)),
        PCRPrimer("indexed P7 PCR primer",P7_TEMPLATE,anneal_nt=len(nx.S7)),
        statement="The original 2017 preprint prints both PCR-primer templates. The selected gap-filled strand is reconstructed from the RT handle and s7 pairing sites; the canonical Nextera mosaic end is inferred.")
    graph.react("library_pcr","INFERRED — initial extension then PCR selects the RT-handle/s7 fragment and adds PCR-well indexes",
                StepKind.PCR,
                ((cleaned,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
                 (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (amplified,),geometry=geometry)
    pcr_beads=graph.input("pcr_cleanup_beads","AMPure XP beads (0.8 volumes)",MaterialKind.SELECTION)
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    pcr_waste=graph.product("pcr_cleanup_waste","material not recovered by post-PCR bead purification",
                            MaterialKind.WASTE)
    graph.react("pcr_cleanup","Pool PCR wells and purify with 0.8 volumes AMPure XP (Qubit, 6% TBE-PAGE check)",
                StepKind.CLEANUP,((amplified,InputRole.TEMPLATE),(pcr_beads,InputRole.SELECTION)),
                (library,pcr_waste),changes_molecule=False)
    graph.discard(pcr_waste,"Material not recovered by the published post-PCR bead cleanup; its composition is unassigned.")
    lengths={"Read 1":18,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":52}
    graph.sequence(library,SequencingPlan(
        library=final,primers=SEQ_PRIMERS,acquired_roles=tuple(lengths),
        read_lengths=lengths,
        read_uses={"Read 1":"extract","Index 1 (i7)":"demultiplex",
                   "Index 2 (i5)":"demultiplex","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "fixed_cells","rt_primer","reverse_transcriptase","rt_stop_solution","reverse_transcribe",
        "facs","dapi","sort","second_strand_polymerase","second_strand",
        "carrier_dna","nextera_transposome","tn5","tagment","p5_primer","p7_primer",
        "dna_binding_buffer","cleanup_beads","tagmentation_cleanup",
        "pcr_polymerase","library_pcr","pcr_products","pcr_cleanup_beads","pcr_cleanup",
        "library"),source)
    graph.annotate(Highlight("tagmentation_cleanup","Final-article overview lists a lysis step",
        "The final Science overview (step iv, p. 1) lists second-strand synthesis, Tn5 transposition, "
        "lysis and PCR. The selected preprint Methods instead stop tagmentation with Zymo DNA binding "
        "buffer and purify on AMPure XP before PCR, with no separate lysis step. The final supplementary "
        "Methods were not available, so the published procedure may differ from the drawn preprint route."))
    return graph.require_publishable()
