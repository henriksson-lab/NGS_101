"""Source-owned microSPLiT reaction order on shared split-pool molecular states."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "split-seq__10.1126+science.aam8999" / "tools"))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import splitseq as S
from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Scene, Segment, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind


# Mirrored 2021 Science Table S3, sheet General Oligonucleotides, rows 15–16.
# Table S3 prints C*T and describes * as a "phosphorothioated base".
BC_0243_CORE = "ACACTCTTTCCCTACACGACGCTCTTCCGATC"
BC_0243_TOP = BC_0243_CORE + "T"
BC_0244_BOTTOM = "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


def published_adapter_scene() -> Scene:
    if revcomp(BC_0244_BOTTOM) != BC_0243_CORE:
        raise ValueError("published BC_0243/BC_0244 paired cores disagree")
    sc = Scene()
    sc.strand("BC_0243", [Segment("paired core", BC_0243_CORE),
                           Segment("3′ T overhang", "T")],
              label="BC_0243", mod3="terminal C*T; *phosphorothioated base (2021 Table S3 mirror)")
    sc.anneal("BC_0244", [Segment("paired core, opposite strand", BC_0244_BOTTOM)],
              to="BC_0243", pair=("paired core, opposite strand", "paired core"),
              label="BC_0244")
    sc.ligation_end("BC_0243", "3′ T overhang", kind="cohesive",
                    polarity="3-prime", owner="adapter",
                    text="published 3′ T adapter overhang")
    return sc


# Mirrored 2021 Science Table S3, sheet General Oligonucleotides, rows 5 and 17–19
# (same bases in 2024 Nature Protocols Table 3 General Oligonucleotide).
BC_0335_ROUND2_LINKER = "CCACAGTCTCAAGCACG"
BC_0284_ROUND3_LINKER = "TACGCCGATGCGAAACATCG"
BC_0340_ROUND2_BLOCKER = "CGTGCTTGAGACTGTGG"
BC_0066_ROUND3_BLOCKER = "GTGGCCGATGTTTCGCATCGGCGTACGACT"
SPLINTS = {2: ("BC_0335", BC_0335_ROUND2_LINKER), 3: ("BC_0284", BC_0284_ROUND3_LINKER)}
BLOCKERS = {2: ("BC_0340", BC_0340_ROUND2_BLOCKER), 3: ("BC_0066", BC_0066_ROUND3_BLOCKER)}
# The splint scenes are computed from the barcode-oligo arms (splitseq.LINK2/LINK3);
# the transcribed linkers must be those splints, and each blocker must cover its linker.
for _n, (_name, _linker) in SPLINTS.items():
    if _linker != (S.LINK2 if _n == 2 else S.LINK3):
        raise ValueError(f"{_name} no longer equals the computed round-{_n} splint")
    if revcomp(_linker) not in BLOCKERS[_n][1]:
        raise ValueError(f"{BLOCKERS[_n][0]} no longer covers the {_name} linker")


def splint_state(n: int) -> Construct:
    name, linker = SPLINTS[n]
    return Construct([Segment(name, linker)])


# Mirrored 2021 Science Table S3, sheet Round 1 RT Barcodes, row E1 (Round1_49):
# \5Phos\ACTGTGGCTGCTTTGTNNNNNN. Published Methods (PMC8269303 JATS line 165)
# put random-hexamer and anchored dT15VN primers in every RT well at 1:2.
HEXAMER_R1_BARCODE = "CTGCTTTG"


def hexamer_round1():
    return (S._s("round-1 handle", S.R1_HANDLE),
            S._s("barcode 1", HEXAMER_R1_BARCODE, "cbc"),
            S._s("spacer T", "T"),
            S._s("random hexamer", "NNNNNN", placeholder=True))


def start_microsplit(graph: ProtocolGraph, *, source: str):
    """Published order: fix/permeabilize, in-cell PAP tailing, then RT barcoding."""
    cells = graph.input("fixed_cells", "formaldehyde-fixed, Tween-20/lysozyme-permeabilized bacteria",
                        MaterialKind.BIOLOGICAL, primary=True,
                        source=source + " Methods (JATS line 149)")
    pap = graph.input("poly_a_polymerase", "E. coli Poly(A) Polymerase I (NEB) + ATP",
                      MaterialKind.CATALYST, source=source + " Methods (JATS line 153)")
    tailed = graph.product("polyadenylated_cells", "cells with in-situ polyadenylated RNA",
                           MaterialKind.BIOLOGICAL, source=source + " Methods (JATS line 153)")
    graph.react("polyadenylate", "In-cell polyadenylation (PAP, ATP, 37 °C 30 min)",
                StepKind.TAILING, ((cells, InputRole.TEMPLATE), (pap, InputRole.CATALYST)),
                (tailed,), source=source + " Methods (JATS line 153)")
    dt_oligo = graph.input("round1_primer", "round-1 barcoded anchored dT15VN primer (5′ phosphate)",
                           MaterialKind.OLIGO, state=tuple(S.round1("microSPLiT")),
                           source=source + " Table S3 mirror, Round 1 RT Barcodes A1")
    hex_oligo = graph.input("round1_hexamer", "round-1 barcoded random-hexamer primer (5′ phosphate)",
                            MaterialKind.OLIGO, state=hexamer_round1(),
                            source=source + " Table S3 mirror, Round 1 RT Barcodes E1; "
                                            "1:2 hexamer:dT per well (Methods JATS line 165)")
    rt = graph.input("reverse_transcriptase", "Maxima H Minus reverse transcriptase",
                     MaterialKind.CATALYST, source=source + " Methods (JATS line 165)")
    first = graph.product("round1_cdna",
                          "round-1 barcoded first-strand cDNA (dT-primed form drawn; "
                          "hexamer-primed cDNA carries the same handle)",
                          MaterialKind.MOLECULAR, state=S.barcoded_first_strand("microSPLiT", 1),
                          source=source)
    graph.react("reverse_transcribe", "Install barcode 1 during in-cell RT (23 °C 10 min, 50 °C 50 min)",
                StepKind.REVERSE_TRANSCRIPTION,
                ((tailed, InputRole.TEMPLATE), (dt_oligo, InputRole.STRUCTURAL),
                 (hex_oligo, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (first,),
                source=source + " Methods (JATS line 165)")
    ligase = graph.input("dna_ligase", "T4 DNA ligase", MaterialKind.CATALYST,
                         source=source + " Methods (JATS line 169)")
    return first, rt, ligase


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1126/science.aba5257"
    graph = ProtocolGraph("microSPLiT split-pool RNA library", authority=source)
    first, rt, ligase = start_microsplit(graph, source=source)
    second = S.ligate_split_pool_round(graph, first, ligase, n=2,
                                       protocol="microSPLiT", splint_state=splint_state(2),
                                       source=source)
    third = S.ligate_split_pool_round(graph, second, ligase, n=3,
                                      protocol="microSPLiT", splint_state=splint_state(3),
                                      source=source)
    lysis = graph.input("lysis_proteinase_k", "SDS lysis buffer + proteinase K, 55 °C 2 h",
                        MaterialKind.CATALYST, source=source + " Methods (JATS line 175)")
    lysate = graph.product("sublibrary_lysate",
                           "sublibrary lysate with crosslink-reversed barcoded cDNA",
                           MaterialKind.MOLECULAR,
                           state=S.barcoded_first_strand("microSPLiT", 3),
                           source=source + " Methods (JATS line 175)")
    graph.react("lyse", "Split into sublibraries, freeze-thaw and lyse with proteinase K",
                StepKind.DEPROTEINIZATION,
                ((third, InputRole.TEMPLATE), (lysis, InputRole.CATALYST)), (lysate,),
                changes_molecule=False, relations=(S.MolecularRelation.STRAND_ORDER,),
                source=source + " Methods (JATS line 175)")
    captured = S.capture_split_pool_cdna(graph, lysate, protocol="microSPLiT")
    amplified = S.switch_and_amplify_split_pool(graph, captured, rt,
                                                protocol="microSPLiT", source=source)

    fragmentation = graph.input("fragmentation_mix", "WGS enzymatic fragmentation mix",
                                MaterialKind.CATALYST)
    fragmented = graph.product("fragmented_cdna", "enzymatically fragmented barcoded cDNA",
                               MaterialKind.MOLECULAR)
    graph.react("fragment", "Fragment amplified cDNA", StepKind.FRAGMENTATION,
                ((amplified, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CATALYST)), (fragmented,))
    prepared = graph.product("prepared_fragments", "end-repaired, dA-tailed cDNA fragments",
                             MaterialKind.MOLECULAR, inferred=True,
                             source=source + " (dA tail inferred from the disclosed 3′-T adapter)")
    graph.react("end_prepare", "Prepare fragment ends for the 3-prime-T adapter",
                StepKind.END_REPAIR,
                ((fragmented, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CATALYST)), (prepared,),
                source=source + " (end repair/dA tail inferred from adapter compatibility)")

    adapter = graph.input("library_adapter", "pre-annealed BC_0243/BC_0244 adapter",
                          MaterialKind.ADAPTER, state=published_adapter_scene(),
                          source=source + " Table S3 mirror (terminal BC_0243 C*T mark)")
    library_ligase = graph.input("library_ligase", "WGS DNA ligase", MaterialKind.CATALYST)
    template = graph.product("adapter_ligated", "adapter-ligated barcoded cDNA fragments",
                             MaterialKind.MOLECULAR)
    graph.react("adapter_ligation", "Ligate the 3-prime-T library adapter", StepKind.LIGATION,
                ((prepared, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (library_ligase, InputRole.CATALYST)), (template,),
                geometry=reported_adapter_end(
                    "BC_0243/BC_0244 adapter", kind="cohesive", polarity="3-prime",
                    terminal_sequence="T",
                    statement="The mirrored 2021 Table S3 discloses a duplex adapter with a terminal 3′ T and marks BC_0243 C*T as phosphorothioated; the original publisher attachment is gated."))

    forward = graph.input("pcr_forward", "library PCR forward primer", MaterialKind.OLIGO,
                          state=Construct([Segment("BC_0027", S.MICRO_LIBRARY_FWD)]),
                          source=source + " Supplementary Table S3")
    reverse = graph.input("pcr_reverse", "indexed library PCR reverse primer",
                          MaterialKind.OLIGO,
                          state=Construct([Segment("BC_0076 / TSBC07", S.MICRO_LIBRARY_REV)]),
                          source=source + " Supplementary Table S3")
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    final = S.final_boundary("microSPLiT")
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final,
                            source=source + " Supplementary Table S3")
    graph.react("library_pcr", "Build the indexed sequencing library", StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.exact(
                    final, PCRPrimer("BC_0027", S.MICRO_LIBRARY_FWD),
                    PCRPrimer("BC_0076 / TSBC07", S.MICRO_LIBRARY_REV),
                    statement="Supplementary Table S3 prints the complete terminal PCR pair."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=S.MICRO_SEQ_PRIMERS,
        acquired_roles=("Read 1", "Read 2", "Index 1 (i7)"),
        read_lengths={"Read 1": 74, "Read 2": 86, "Index 1 (i7)": 6},
        read_uses={"Read 1": "include", "Read 2": "extract", "Index 1 (i7)": "demultiplex"}))
    graph.declare_reviewed_scope((
        "fixed_cells", "poly_a_polymerase", "polyadenylate", "polyadenylated_cells",
        "round1_primer", "round1_hexamer", "reverse_transcriptase", "reverse_transcribe",
        "dna_ligase", "round2_adapter", "round2_splint", "ligate_round2",
        "round3_adapter", "round3_splint", "ligate_round3", "lysis_proteinase_k", "lyse",
        "sublibrary_lysate", "streptavidin", "capture",
        "tso", "template_switch", "preamp_forward", "preamp_reverse",
        "preamp_polymerase", "cdna_pcr", "fragmentation_mix", "fragment", "end_prepare",
        "library_adapter", "library_ligase", "adapter_ligation", "pcr_forward",
        "pcr_reverse", "pcr_polymerase", "library_pcr", "library"), source)
    return graph.require_publishable()
