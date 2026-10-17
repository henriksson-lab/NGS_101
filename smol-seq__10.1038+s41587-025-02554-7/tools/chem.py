"""smol-seq: bead-immobilized barcode structure-switching aptamers read by sequencing.

Source: Tan and Fraser, Nat Biotechnol 2026 (PMC12807852), Methods 'Ligand-binding
assay' and 'Multiplexing and sequencing of bSSAs'; oligos from Supplementary Data 1
(MOESM2 sheet 'SuppData1 Oligonucleotides').  The drawn sensor is the mefloquine (MQ)
bSSA with barcode 1, one of nine barcode SROs (MQ, ATP, PQ x three barcodes) pooled in
the sequenced experiment; the other sensors follow the same published SRO design.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind)
from protocol_model import molecule, role

TITLE = "smol-seq — small-molecule sequencing"
NOTES = "01_smol_seq.html"
SOURCE = ('Tan and Fraser, <i>Nature Biotechnology</i> (2026), '
          '<a href="https://doi.org/10.1038/s41587-025-02554-7">doi:10.1038/s41587-025-02554-7</a>; '
          'oligos from Supplementary Data 1.')
SUMMARY = ("A bead-immobilized ligand-binding oligo (LBO) holds a barcoded short release oligo (SRO) "
           "by a short stem; cognate ligand switches the LBO into a stem-loop and releases the SRO. "
           "Released and boiled-off SROs are amplified in two PCR rounds and counted by sequencing.")
CAVEAT = ("The drawn sensor is the mefloquine bSSA with barcode 1; nine barcode SROs (ATP, mefloquine, "
          "piperaquine; three barcodes each) were pooled in the sequenced experiment. The SRO's random "
          "N5–7 diversity region is drawn as 6 N. INFERRED — i5/i7 index lengths and the sequencing "
          "primers are not stated (the paper reports ~4 million 50-bp single-end MiSeq reads); the "
          "standard TruSeq Read 1 primer is inferred from the source-built Read 1 site and index "
          "windows are drawn as 8 nt.")

ANALYTE_BARCODE = feature("analyte_identity", "feature_barcode", "fixed")
I5_FEATURE = feature("sample_index_i5", "sample_index", "unknown")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")

# Supplementary Data 1, 'Barcode seq multiplexing' and 'Round 1/2 PCR' rows.
MQ_LBO = "ATACCGCTCTCGGGACGACGGCAGTCTATACCCCGTATCGCCGAAGGTTGTCGTCCC"   # /5Biosg/ 5'
MQ_SRO_STEM = "CGTCGTCCCGAGAGTAGC"
MQ_BARCODE_1 = "AGTCTGAC"
SRO_P5_HANDLE = "AGATCGGAAGAGCGTCGTGTAG"
ROUND1_REVERSE_P5 = "ACACTCTTTCCCTACACGACGCTCTTCCGATC"
ROUND1_FORWARD_HANDLE = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
ROUND1_FORWARD_MQ = ROUND1_FORWARD_HANDLE + MQ_SRO_STEM
ROUND2_I5_3PRIME = "ACACTCTTTCCCTACACGACG"
ROUND2_I7_3PRIME = "GTGACTGGAGTTCAGACGTG"
N_DRAWN = 6   # source: N5-7

# LBO bases 7-20 pair with the SRO 5' 14 nt (computed below by Scene.anneal).
LBO_PAIRED = MQ_LBO[6:20]


def lbo_segments():
    return [Segment("LBO 5' end", MQ_LBO[:6]),
            Segment("LBO stem", LBO_PAIRED),
            Segment("LBO ligand-binding/switch region", MQ_LBO[20:])]


def sro_segments():
    return [Segment("SRO stem", MQ_SRO_STEM[:14]),
            Segment("SRO stem tail", MQ_SRO_STEM[14:]),
            Segment("MQ barcode 1", MQ_BARCODE_1, "cbc", feature=ANALYTE_BARCODE),
            Segment("N5-7 diversity bases", "N" * N_DRAWN, placeholder=True),
            Segment("P5-side PCR handle", SRO_P5_HANDLE, "r1")]


def sensor_scene():
    sc = Scene()
    sc.strand("LBO", lbo_segments(), label="5'-biotin LBO", mod5="biotin")
    sc.anneal("SRO", sro_segments(), to="LBO", pair=("SRO stem", "LBO stem"), label="barcode SRO",
              unpaired=("SRO stem tail", "MQ barcode 1", "N5-7 diversity bases", "P5-side PCR handle"))
    return sc


SRO = Construct(sro_segments(), name="MQ barcode 1 SRO")

FINAL_LIBRARY = Construct([
    Segment("P5", il.P5, "p5"),
    Segment("i5 index", "I" * 8, "cbc", placeholder=True, inferred=True, feature=I5_FEATURE),
    Segment("round-1 reverse P5 primer / Read 1 site", ROUND1_REVERSE_P5 + "T", "r1"),
    Segment("N5-7 complement", "N" * N_DRAWN, placeholder=True),
    Segment("MQ barcode 1 (reverse complement)", revcomp(MQ_BARCODE_1), "cbc", feature=ANALYTE_BARCODE),
    Segment("MQ SRO stem (reverse complement)", revcomp(MQ_SRO_STEM)),
    Segment("round-1 forward handle (reverse complement)", revcomp(ROUND1_FORWARD_HANDLE), "r2"),
    Segment("i7 index (reverse complement)", "I" * 8, "cbc", placeholder=True, inferred=True, feature=I7_FEATURE),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
], name="smol-seq MQ barcode-1 library")

SEQ_PRIMERS = sp.standard_primers(FINAL_LIBRARY, ("Read 1",))
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=("Read 1",)):
    raise ValueError("smol-seq final library: " + "; ".join(problems))
# The round-2 primer 3' ends must land on the round-1 product they extend.
for _p in (ROUND2_I5_3PRIME, ROUND2_I7_3PRIME):
    _top = "".join(s.top for s in FINAL_LIBRARY)
    if _p not in _top and revcomp(_p) not in _top:
        raise ValueError(f"round-2 primer 3' end {_p} has no site on the round-1 product")

FINAL_CAPTION = ("INFERRED — 8-nt i5/i7 windows drawn; index lengths are not stated. MQ barcode-1 amplicon: "
                 "Read 1 crosses the N5–7 diversity bases, then the 8-nt analyte barcode and the MQ stem.")
SEQUENCING_INTRO = ("MiSeq, ~4 million 50-bp single-end reads (Methods). Read 1 primer identity and index-read "
                    "cycles are not stated; the TruSeq Read 1 primer is inferred from the source-built site.")


def oligos():
    return (
        f"<p>Supplementary Data 1, mefloquine LBO (5′→3′): /5Biosg/{MQ_LBO}</p>",
        f"<p>Supplementary Data 1, mefloquine barcode 1 SRO (5′→3′): {MQ_SRO_STEM} <cbc>{MQ_BARCODE_1}</cbc> N5-7 {SRO_P5_HANDLE}</p>",
        f"<p>Round 1: Reverse P5 primer {ROUND1_REVERSE_P5}; Forward MQ primer {ROUND1_FORWARD_HANDLE} {MQ_SRO_STEM}</p>",
        f"<p>Round 2: Forward i5 primer {il.P5} [i5 index] {ROUND2_I5_3PRIME}; Reverse i7 primer {il.P7} [i7 index] {ROUND2_I7_3PRIME}</p>",
    )


def _library_prep(graph, sample, prefix, label, pcr1_cycles, primers):
    rev1, fwd1, rev2i5, rev2i7, pol = primers
    kit = graph.input(prefix + "_oligo_cleanup", "Monarch PCR & DNA cleanup kit (oligonucleotide protocol)",
                      MaterialKind.SELECTION)
    clean = graph.product(prefix + "_clean", f"purified {label} SRO pool with spike-ins", MaterialKind.MOLECULAR,
                          chemical_state=SRO)
    w0 = graph.product(prefix + "_oligo_cleanup_waste", "column flow-through and washes", MaterialKind.WASTE)
    graph.react(prefix + "_purify", f"Pool replicates, add spike-in controls and purify the {label} single-stranded SROs",
                StepKind.CLEANUP, ((sample, InputRole.TEMPLATE), (kit, InputRole.SELECTION)), (clean, w0))
    graph.discard(w0, "column flow-through and washes")
    r1 = graph.product(prefix + "_round1", f"{label} round-1 product with partial P5/P7 handles", MaterialKind.MOLECULAR,
                       state=FINAL_LIBRARY, inferred=True)
    graph.react(prefix + "_pcr1", f"Round-1 PCR adds partial P5/P7 handles ({pcr1_cycles} cycles)", StepKind.PCR,
                ((clean, InputRole.TEMPLATE), (rev1, InputRole.STRUCTURAL), (fwd1, InputRole.STRUCTURAL),
                 (pol, InputRole.CATALYST)), (r1,),
                geometry=PCRGeometry.unavailable(FINAL_LIBRARY, reason="Round-1 product and full library are drawn as one final construct; the intermediate is not separately represented."))
    kit2 = graph.input(prefix + "_pcr_cleanup", "Monarch PCR & DNA cleanup kit (PCR protocol); 1:5 dilution",
                       MaterialKind.SELECTION)
    r1c = graph.product(prefix + "_round1_clean", f"purified {label} round-1 product", MaterialKind.MOLECULAR,
                        state=FINAL_LIBRARY, inferred=True)
    w1 = graph.product(prefix + "_pcr_cleanup_waste", "primers and flow-through", MaterialKind.WASTE)
    graph.react(prefix + "_pcr1_cleanup", "Clean up round-1 PCR and dilute 1:5", StepKind.CLEANUP,
                ((r1, InputRole.TEMPLATE), (kit2, InputRole.SELECTION)), (r1c, w1))
    graph.discard(w1, "column flow-through and washes")
    r2 = graph.product(prefix + "_round2", f"{label} dual-indexed library", MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react(prefix + "_pcr2", "Round-2 PCR adds i5/i7 indices and complete P5/P7 (100 µL, 8 cycles)", StepKind.PCR,
                ((r1c, InputRole.TEMPLATE), (rev2i5, InputRole.STRUCTURAL), (rev2i7, InputRole.STRUCTURAL),
                 (pol, InputRole.CATALYST)), (r2,),
                geometry=PCRGeometry.unavailable(FINAL_LIBRARY, reason="Index values are selected per library and not printed."))
    kit3 = graph.input(prefix + "_size_select", "Monarch PCR cleanup then NucleoMag NGS beads at 1.2×",
                       MaterialKind.SELECTION)
    lib = graph.product(prefix + "_library", f"{label} sequencing library", MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    w2 = graph.product(prefix + "_size_select_waste", "small fragments in the bead supernatant", MaterialKind.WASTE)
    graph.react(prefix + "_cleanup", "Clean up and size-select (NucleoMag 1.2× bead:sample)", StepKind.SELECTION,
                ((r2, InputRole.TEMPLATE), (kit3, InputRole.SELECTION)), (lib, w2))
    graph.discard(w2, "bead supernatant and washes")
    graph.sequence(lib, SequencingPlan(library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
                                       acquired_roles=("Read 1",), read_lengths={"Read 1": 50}))
    return [prefix + x for x in ("_oligo_cleanup", "_clean", "_oligo_cleanup_waste", "_purify", "_round1", "_pcr1",
                                 "_pcr_cleanup", "_round1_clean", "_pcr_cleanup_waste", "_pcr1_cleanup", "_round2",
                                 "_pcr2", "_size_select", "_library", "_size_select_waste", "_cleanup")]


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    lbo = graph.input("lbo", "5′-biotin mefloquine LBO (25 pmol)", MaterialKind.OLIGO, state=MQ_LBO,
                      source="Supplementary Data 1", primary=True)
    sro = graph.input("sro", "mefloquine barcode-1 SRO (125 pmol)", MaterialKind.OLIGO, chemical_state=SRO,
                      source="Supplementary Data 1")
    rows = tuple(sensor_scene().rows())
    duplex = graph.product("ssa", "LBO:SRO barcode SSA", MaterialKind.MOLECULAR, state=rows)
    graph.react("anneal", "Anneal LBO and SRO (95°C 5 min, slow cool to 25°C)", StepKind.HYBRIDIZATION,
                ((lbo, InputRole.TEMPLATE), (sro, InputRole.STRUCTURAL)), (duplex,))
    beads = graph.input("streptavidin_beads", "Dynabeads MyOne Streptavidin C1 (5 µL per sample)", MaterialKind.SELECTION)
    bound = graph.product("bead_ssa", "bead-immobilized barcode SSA (sensor pool: MQ, ATP, PQ × 3 barcodes)",
                          MaterialKind.MOLECULAR, state=rows)
    excess = graph.product("excess_sro", "unbound excess SRO", MaterialKind.WASTE)
    graph.react("immobilize", "Bind to streptavidin beads 30 min, wash twice, 45 min, wash once", StepKind.AFFINITY_CAPTURE,
                ((duplex, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (bound, excess))
    graph.discard(excess, "removed with the bead washes")
    ligand = graph.input("ligand", "sample: 125 µM ATP, 30 µM MQ, 25 µM PQ or PBS buffer", MaterialKind.BIOLOGICAL,
                         state=molecule("small-molecule sample", (role("cognate ligand"),)))
    released = graph.product("released", "'released' sample: SROs displaced into the supernatant", MaterialKind.MOLECULAR,
                             chemical_state=SRO)
    retained = graph.product("retained", "beads with SROs still paired to unswitched LBOs", MaterialKind.MOLECULAR, state=rows)
    graph.react("release", "Ligand binding (45 min, rotator) folds the LBO stem-loop and displaces the SRO; collect supernatant",
                StepKind.BINDING, ((bound, InputRole.TEMPLATE), (ligand, InputRole.ANALYTE)), (released, retained))
    heat = graph.input("strand_separation", "Mg-free strand separation buffer, 95°C", MaterialKind.CONDITION)
    boiled = graph.product("boiled", "'boiled' sample: SROs heat-released from the beads", MaterialKind.MOLECULAR,
                           chemical_state=SRO)
    spent = graph.product("spent_beads", "beads carrying biotin-LBO", MaterialKind.WASTE)
    graph.react("boil", "Boil off the remaining SROs to measure total sensor (percentage release)",
                StepKind.THERMAL_TREATMENT, ((retained, InputRole.TEMPLATE), (heat, InputRole.CONDITION)), (boiled, spent))
    graph.discard(spent, "magnetically retained; the supernatant is collected")
    rev1 = graph.input("round1_reverse_p5", "Round-1 Reverse P5 primer", MaterialKind.OLIGO,
                       state=ROUND1_REVERSE_P5, source="Supplementary Data 1")
    fwd1 = graph.input("round1_forward_mq", "Round-1 Forward MQ primer", MaterialKind.OLIGO,
                       state=ROUND1_FORWARD_MQ, source="Supplementary Data 1")
    i5 = graph.input("round2_i5", "Round-2 Forward i5 primer (index value selected per library)", MaterialKind.OLIGO,
                     source="Supplementary Data 1")
    i7 = graph.input("round2_i7", "Round-2 Reverse i7 primer (index value selected per library)", MaterialKind.OLIGO,
                     source="Supplementary Data 1")
    pol = graph.input("pcr_polymerase", "PCR polymerase (not named in Methods)", MaterialKind.CATALYST)
    primers = (rev1, fwd1, i5, i7, pol)
    scope = ["lbo", "sro", "anneal", "ssa", "streptavidin_beads", "immobilize", "bead_ssa", "excess_sro",
             "ligand", "release", "released", "retained", "strand_separation", "boil", "boiled", "spent_beads",
             "round1_reverse_p5", "round1_forward_mq", "round2_i5", "round2_i7", "pcr_polymerase"]
    scope += _library_prep(graph, released, "rel", "released", 14, primers)
    scope += _library_prep(graph, boiled, "boil", "boiled", 11, primers)
    graph.declare_reviewed_scope(scope, SOURCE)
    return graph.require_publishable()
