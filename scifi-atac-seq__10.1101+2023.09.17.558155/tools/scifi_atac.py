"""Molecular construct model for scifi-ATAC-seq (Zhang et al., Genome Biology 2024)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (AnnotationSlot, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


SPACER_A = "GATATGTGATAATGAGGAC"
SPACER_B = "TGAATGTGTAGAAGACAGA"
BC_A_EXAMPLE = "CGTAT"
BC_B_EXAMPLE = "TCGGA"
# The paper names no i7 set; CG000209 (Rev G inspected) says each Single Index Kit N
# Set A well is a mix of four 8-nt index oligos, so the read is an unassigned placeholder.
I7_INDEX_NT = 8
R1_SEQ = nx.S5 + SPACER_A
R2_SEQ = nx.S7 + SPACER_B
I1_SEQ = revcomp(R2_SEQ)
I2_SEQ = revcomp(R1_SEQ)
TN5_A = feature("cell_tn5_a", "cell_barcode", "combinatorial", group="cell_id", part="Tn5 A")
TN5_B = feature("cell_tn5_b", "cell_barcode", "combinatorial", group="cell_id", part="Tn5 B")
SAMPLE_I7 = feature("sample_i7", "sample_index", "whitelist",
                    whitelist="10x Single Index Kit N Set A (well not stated; four-oligo mix per well)")
# CG000209 Rev G, Appendix oligo sequences, Step 4.1: SI-PCR Primer B is the 5'
# 22 nt of P5; the Single Index N primer is P7 + 8-nt i7 + s7 (as printed).
SI_PCR_PRIMER_B = il.P5[:22]


def adapter_a() -> list[Segment]:
    return [seg("s5", nx.S5, "s5"), seg("spacer A", SPACER_A),
            seg("Tn5 barcode A", "A" * 5, "cbc", placeholder=True, feature=TN5_A),
            seg("mosaic end", nx.ME, "me")]


def adapter_b() -> list[Segment]:
    return [seg("s7", nx.S7, "s7"), seg("spacer B", SPACER_B),
            seg("Tn5 barcode B", "B" * 5, "cbc", placeholder=True, feature=TN5_B),
            seg("mosaic end", nx.ME, "me")]


def bead_oligo() -> list[Segment]:
    """Unmodified 10x scATAC v1.1 gel-bead oligo from its vendor guide."""
    return [seg("P5", il.P5, "p5"),
            seg("GEM barcode", "C" * 16, "cbc", placeholder=True,
                feature=feature("cell_gem", "cell_barcode", "whitelist", whitelist="10x scATAC v1.1 whitelist", group="cell_id", part="GEM")),
            seg("s5", nx.S5, "s5")]


def indexed_fragment(insert_nt: int = 28) -> Construct:
    return Construct([*adapter_a(), seg("accessible DNA", "X" * insert_nt, placeholder=True),
                      seg("opposite ME reverse complement", nx.ME_RC, "me"),
                      seg("Tn5 barcode B reverse complement", "B" * 5, "cbc", placeholder=True,
                          feature=TN5_B),
                      seg("spacer B reverse complement", revcomp(SPACER_B)),
                      seg("s7 reverse complement", nx.S7_RC, "s7")],
                     name="A/B-indexed accessible fragment")


SEQ_PRIMERS = (
    sp.custom("Read 1", "scifi 1_Read1", R1_SEQ, "Zhang et al. 2024 Additional file 2 Table S1"),
    sp.custom("Index 1 (i7)", "scifi 2_Index1(i7)", I1_SEQ, "Zhang et al. 2024 Additional file 2 Table S1"),
    sp.custom("Index 2 (i5)", "scifi 3_Index2(i5)", I2_SEQ, "Zhang et al. 2024 Additional file 2 Table S1"),
    sp.custom("Read 2", "scifi 4_Read2", R2_SEQ, "Zhang et al. 2024 Additional file 2 Table S1"),
)


def final_library(insert_nt: int = 28) -> Construct:
    lib = Construct([
        *bead_oligo(), seg("spacer A", SPACER_A),
        seg("Tn5 barcode A", "A" * 5, "cbc", placeholder=True,
            feature=TN5_A),
        seg("mosaic end", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("opposite ME reverse complement", nx.ME_RC, "me"),
        seg("Tn5 barcode B reverse complement", "B" * 5, "cbc", placeholder=True,
            feature=TN5_B),
        seg("spacer B reverse complement", revcomp(SPACER_B)),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 sample index read", "N" * I7_INDEX_NT, "cbc", placeholder=True,
            feature=SAMPLE_I7),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scifi-ATAC-seq sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final scifi-ATAC-seq construct: " + "; ".join(problems))
    return lib


def droplet_indexed_fragment(insert_nt: int = 28) -> Construct:
    return Construct([
        *bead_oligo(),seg("spacer A",SPACER_A),
        seg("Tn5 barcode A","A"*5,"cbc",placeholder=True,
            feature=TN5_A),
        seg("mosaic end",nx.ME,"me"),seg("accessible genomic DNA","X"*insert_nt,placeholder=True),
        seg("opposite ME reverse complement",nx.ME_RC,"me"),
        seg("Tn5 barcode B reverse complement","B"*5,"cbc",placeholder=True,
            feature=TN5_B),
        seg("spacer B reverse complement",revcomp(SPACER_B)),
        seg("s7 reverse complement",nx.S7_RC,"s7")],
        name="GEM-indexed scifi-ATAC fragment")


def si_pcr_forward() -> Construct:
    return Construct([seg("SI-PCR Primer B (P5 5' 22 nt)", SI_PCR_PRIMER_B, "p5")],
                     name="10x SI-PCR Primer B")


def si_pcr_reverse() -> Construct:
    return Construct([seg("P7", il.P7, "p7"),
                      seg("i7 sample index", "N" * I7_INDEX_NT, "cbc", placeholder=True,
                          feature=SAMPLE_I7),
                      seg("s7", nx.S7, "s7")], name="10x Single Index N i7 primer")


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def protocol_graph() -> ProtocolGraph:
    source=("Defining source: Zhang et al., Genome Biology 25:90 (2024), "
            "https://doi.org/10.1186/s13059-024-03235-5 (Additional files 1-2); "
            "droplet chemistry: 10x CG000209 v1.1 guide (Rev E cited; Rev G inspected)")
    graph=ProtocolGraph("scifi-ATAC-seq",authority=source)
    nuclei=graph.input("nuclei","permeabilized nuclei",MaterialKind.BIOLOGICAL,primary=True)
    adapter_a_input=graph.input("adapter_a","Tn5 A barcode adapter",MaterialKind.ADAPTER,
                                state=tuple(adapter_a()))
    adapter_b_input=graph.input("adapter_b","Tn5 B barcode adapter",MaterialKind.ADAPTER,
                                state=tuple(adapter_b()))
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    indexed=graph.product("indexed_fragments",indexed_fragment().name,MaterialKind.MOLECULAR,
                          state=indexed_fragment())
    tag_waste=graph.product("tagment_waste","A/A, B/B and otherwise non-selected fragments",
                            MaterialKind.WASTE)
    graph.react("tagment","Tagment nuclei in plate wells with one A/B barcode pair",
                StepKind.TAGMENTATION,
                ((nuclei,InputRole.TEMPLATE),(adapter_a_input,InputRole.STRUCTURAL),
                 (adapter_b_input,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),
                (indexed,tag_waste))
    graph.discard(tag_waste,"only fragments retaining the A and B selection boundaries yield the modelled library")
    peers=graph.input("peer_wells","A/B-indexed nuclei from the other plate wells",
                      MaterialKind.MOLECULAR)
    pooled=graph.product("pooled_nuclei","pooled plate-indexed nuclei",MaterialKind.MOLECULAR,
                         state=indexed_fragment())
    graph.react("pool","Pool the 96 tagmentation wells",StepKind.POOL,
                ((indexed,InputRole.TEMPLATE),(peers,InputRole.TEMPLATE)),(pooled,),
                changes_molecule=False)
    droplets=graph.input("droplet_device","overloaded 10x scATAC v1.1 channel",
                         MaterialKind.CONDITION)
    productive=graph.product("productive_droplets","droplets with plate-indexed nuclei and gel bead",
                              MaterialKind.MOLECULAR,state=indexed_fragment())
    empty=graph.product("unproductive_droplets","droplets without usable nucleus/bead combinations",
                        MaterialKind.WASTE)
    graph.react("partition","Overload pooled nuclei into the 10x scATAC channel",
                StepKind.PARTITION,
                ((pooled,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(productive,empty))
    graph.discard(empty,"unproductive droplets do not contribute library molecules")
    bead=graph.input("bead_oligo","10x scATAC v1.1 gel-bead oligo",MaterialKind.ADAPTER,
                     state=tuple(bead_oligo()))
    polymerase=graph.input("gap_fill_polymerase","gap-fill and linear-extension polymerase",
                           MaterialKind.CATALYST)
    extended_state=droplet_indexed_fragment()
    extended=graph.product("droplet_indexed",extended_state.name,MaterialKind.MOLECULAR,
                           state=extended_state)
    graph.react("droplet_extension","Gap-fill and extend from bead s5 to add P5 and GEM barcode",
                StepKind.EXTENSION,
                ((productive,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
                 (polymerase,InputRole.CATALYST)),(extended,))
    vendor_oligos="10x CG000209 Rev G, Appendix oligo sequences, Protocol Step 4.1"
    forward=graph.input("pcr_forward","P5-side sample-index PCR primer (SI-PCR Primer B)",
                        MaterialKind.OLIGO,source=vendor_oligos,state=si_pcr_forward())
    reverse=graph.input("pcr_reverse","indexed P7-side sample-index PCR primer (Single Index N)",
                        MaterialKind.OLIGO,source=vendor_oligos,state=si_pcr_reverse())
    pcr_pol=graph.input("pcr_polymerase","sample-index PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    graph.react("library_pcr","Select the B end and add the i7/P7 sample index",
                StepKind.PCR,
                ((extended,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.exact(
                    extended_state,
                    PCRPrimer("SI-PCR Primer B",SI_PCR_PRIMER_B),
                    PCRPrimer("Single Index N i7 primer",
                              "".join(x.top for x in si_pcr_reverse()),
                              anneal_nt=len(nx.S7)),
                    statement=("CG000209 Rev G prints SI-PCR Primer B and the Single Index N "
                               "primer; Primer B primes from the bead-added P5 and the i7 "
                               "primer's s7 3' end anneals to the B-end s7, so only "
                               "GEM-indexed A/B fragments are amplified.")))
    graph.sequence(library,SequencingPlan(
        library=final,primers=SEQ_PRIMERS,
        acquired_roles=("Read 1","Index 1 (i7)","Index 2 (i5)","Read 2"),
        read_uses={"Read 1":"extract_then_include","Index 1 (i7)":"demultiplex",
                   "Index 2 (i5)":"extract","Read 2":"extract_then_include"}))
    graph.declare_reviewed_scope((
        "nuclei","adapter_a","adapter_b","tn5","tagment","peer_wells","pool",
        "droplet_device","partition","bead_oligo","gap_fill_polymerase",
        "droplet_extension","pcr_forward","pcr_reverse","pcr_polymerase",
        "library_pcr","library"),source)
    graph.annotate(SourceDiscrepancy(
        target="tagment", title="A/B plate axes in the published Methods",
        claims=(
            SourceClaim("Zhang et al. Methods 'Indexed Tn5 tagmentation and pooling' and Supplementary Methods III.5-6",
                        "Place 12 A adapters by rows and 8 B adapters by columns."),
            SourceClaim("Zhang et al. Results, two-sided barcode rationale",
                        "The plate has 8 rows by 12 columns, with 15 µL in 8 rows and 10 µL in 12 columns.")),
        model_choice="Twelve A adapters occupy columns and eight B adapters occupy rows.",
        evidence="That orientation fits the 8-by-12 plate and the published adapter counts.",
        uncertainty="The Methods axis wording is unresolved; individual well assignments need the experimental plate map.",
        slot=AnnotationSlot.INTRO))
    graph.annotate(SourceDiscrepancy(
        target="tagment", title="B-adapter count in the supplement",
        claims=(
            SourceClaim("Supplementary Methods step II.1", "States 18 Tn5-ME-Bx adapters."),
            SourceClaim("Supplementary Table S1 and the remaining Methods", "List eight B adapters, B1–B8.")),
        model_choice="The graph uses the eight B members printed in Table S1.",
        evidence="Only B1–B8 are listed, and eight B by twelve A yields the stated 96 combinations.",
        uncertainty="The source does not publish a formal correction of the number 18.",
        slot=AnnotationSlot.INTRO))
    graph.annotate(SourceDiscrepancy(
        target="tagment", title="Separate versus mixed Tn5 loading",
        claims=(
            SourceClaim("Main Methods, Tn5 loading", "Annealed A and B adapters are mixed 1:1 before addition to Tn5."),
            SourceClaim("Supplementary Methods, Tn5 assembly", "Load A and B adapters in separate Tn5 stocks (2 µL + 15 µL and 3 µL + 25 µL diluted Tn5).")),
        model_choice="The graph follows separate A- and B-loaded transposomes for combinatorial plating.",
        evidence="The supplement's separate stocks can preserve the intended A-by-B well identities.",
        uncertainty="The main Methods may describe a different preparation; its relation to the plated experiment is unclear.",
        slot=AnnotationSlot.INTRO))
    return graph.require_publishable()
