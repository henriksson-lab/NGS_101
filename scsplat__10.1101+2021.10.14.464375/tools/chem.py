"""Source-owned scSPLAT protocol definition and reaction graph.

Authority: Raine et al., Scientific Reports 12:5772 (2022), doi:10.1038/s41598-022-09798-2
(PMC8986790), Methods "scWGBS library preparation and sequencing" and Supplementary Table S1.
The bioRxiv preprint 10.1101/2021.10.14.464375 names the directory but is superseded.
"""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, oligo
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind, Terminal,
                            Disposition)
from protocol_specs import Definition, get

PUBLISHED = "Raine et al., Scientific Reports (2022), doi:10.1038/s41598-022-09798-2"
NEBNEXT = (PUBLISHED + " names NEBNext dual index primers but does not print them; "
           "the P5/i5/i7/P7 bases are the standard TruSeq-compatible arms (inferred)")

# Supplementary Table S1 (MOESM1): SSO_H_1..SSO_H_32 share this 5′ tail, then a
# six-base cell barcode and eight H (A/C/T) random bases. SSO_H_1 is drawn as the
# representative family member; every member carries 5AmMC6.
SSO_TAIL = "TTCCCTACACGACGCTCTTCCGATCT"
SSO_H_1_BARCODE = "TAGCTT"
H8 = "HHHHHHHH"
# Published Methods / Table S1 SPLAT-adapter oligos.
SS1A = "GACGTGTGCTCTTCCGATCTNNNNNN"     # 5′AmMC6 … 3′AmMo
SS1B = "AGATCGGAAGAGCACACGTC"           # 5′Phos … 3′AmMo

if il.TRUSEQ_READ1[-len(SSO_TAIL):] != SSO_TAIL:
    raise AssertionError("SSO_H tail is not the 3′ 26 nt of TruSeq Read 1")
_READ2_TOP = "A" + il.INDEX1_PRIMER        # top-strand Read 2 / Index 1 arm
if not _READ2_TOP.startswith(SS1B):
    raise AssertionError("ss1b is not the 5′ start of the TruSeq Read 2 arm")

CELL_BARCODE = feature("scsplat_cell_barcode", "cell_barcode", "whitelist",
                       whitelist="32 six-base SSO_H_X barcodes, Supplementary Table S1",
                       note="Methods: first six bases of Read 1; demultiplexed with je")

DEFINITION = Definition(
    "scSPLAT — scalable single-cell WGBS", "10.1101/2021.10.14.464375",
    "Raine et al., Scientific Reports (2022); bioRxiv 2021 preprint",
    "Bisulfite-converted single-cell DNA receives a P5-tail/cell-barcode random-primer copy "
    "(SSO_H_X, Klenow exo−) before pooling and bulk splinted ligation (SPLAT) of the opposite adapter.",
    "truseq", ("cell barcode", "bisulfite-converted insert"),
    "The first six bases of Read 1 are the inline cell barcode; reads were PE150 "
    "(HiSeqX or NovaSeq 6000). Index-read cycles are not stated.",
    "DNA",
    caveat=("INFERRED — the NEBNext dual-index PCR primers are not printed in the paper; the "
            "P5, i5, first seven Read 1 bases, Read 2 remainder, i7 and P7 bases shown are the "
            "standard TruSeq-compatible arms, not source transcriptions. The SSO_H tail, cell "
            "barcode/H8 primer and ss1a/ss1b adapter are printed in Methods and Table S1."))


def _final_library() -> tuple[Construct, tuple]:
    lib = Construct([
        Segment("P5", il.P5, "p5", inferred=True),
        Segment("i5", "J" * 8, "cbc", placeholder=True, inferred=True,
                feature=feature("sample_index_i5", "sample_index", "unknown")),
        Segment("Read 1 arm from i5 primer", il.TRUSEQ_READ1[:-len(SSO_TAIL)], "r1",
                inferred=True),
        Segment("SSO_H P5 tail", SSO_TAIL, "r1"),
        Segment("cell barcode", "X" * 6, "cbc", placeholder=True, feature=CELL_BARCODE),
        Segment("H8 random-primer copy", H8, placeholder=True),
        Segment("bisulfite-converted insert", "X" * 28, placeholder=True),
        Segment("ss1b SPLAT adapter", SS1B, "r2"),
        Segment("Read 2 arm from i7 primer", _READ2_TOP[len(SS1B):], "r2", inferred=True),
        Segment("i7 reverse complement", "I" * 8, "cbc", placeholder=True, inferred=True,
                feature=feature("sample_index_i7", "sample_index", "unknown")),
        Segment("P7 reverse complement", il.P7_RC, "p7", inferred=True),
    ], name="scSPLAT final library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers)
    if problems:
        raise ValueError("invalid scSPLAT endpoint: " + "; ".join(problems))
    return lib, primers


FINAL_LIBRARY, SEQ_PRIMERS = _final_library()


def sso_h_1() -> Scene:
    sc = Scene()
    sc.strand("SSO_H_1", [Segment("SSO_H P5 tail", SSO_TAIL, "r1"),
                          Segment("SSO_H_1 cell barcode", SSO_H_1_BARCODE, "cbc",
                                  feature=CELL_BARCODE),
                          Segment("H8 random primer", H8, placeholder=True)],
              label="SSO_H_1 (representative of SSO_H_1–32)", mod5="AmMC6")
    return sc


def splat_adapter() -> Scene:
    sc = Scene()
    sc.strand("ss1b", [Segment("ss1b ligating strand", SS1B)], label="SPLAT adapter",
              mod5="p", mod3="AmMo")
    sc.anneal("ss1a", [Segment("ss1a paired arm", SS1A[:-6]),
                       Segment("ss1a random splint", "N" * 6)],
              to="ss1b", pair=("ss1a paired arm", "ss1b ligating strand"),
              label="SPLAT adapter", mod5="AmMC6", mod3="",
              unpaired=("ss1a random splint",))
    sc.overhang("ss1a", "ss1a random splint", polarity="3-prime", owner="adapter",
                text="3′ random hexamer splint")
    sc.ligation_end("ss1b", "ss1b ligating strand", kind="splinted", polarity="5-prime",
                    owner="adapter", text="5′ phosphate joins the single-strand 3′-OH")
    return sc


def adapter_geometry() -> AdapterLigationGeometry:
    return AdapterLigationGeometry.known(
        "splinted", splat_adapter(),
        statement="Published Methods/Table S1 print ss1a (5′AmMC6/GACGTGTGCTCTTCCGATCTNNNNNN/3′AmMo) "
                  "and ss1b (5′Phos/AGATCGGAAGAGCACACGTC/3′AmMo); the protruding 3′ N6 of the bottom "
                  "strand hybridizes to the single-stranded 3′ end so T4 DNA ligase joins it to the "
                  "ss1b 5′ phosphate. The Scene verifies that the printed sequences pair.")


def oligos():
    return [
        oligo("SSO_H_1 (Table S1; SSO_H_2–32 differ only in the six-base barcode)",
              (Segment("tail", SSO_TAIL, "r1"), Segment("barcode", SSO_H_1_BARCODE, "cbc"),
               Segment("H8", H8, placeholder=True)), mods="/5AmMC6/"),
        oligo("ss1a (SPLAT adapter bottom strand)", (Segment("ss1a", SS1A),),
              mods="/5AmMC6/", three="/3AmMo/-3′"),
        oligo("ss1b (SPLAT adapter top strand; Methods also call the partner ss2b)",
              (Segment("ss1b", SS1B),), mods="/5Phos/", three="/3AmMo/-3′"),
    ]


def _barcoded_segments() -> list[Segment]:
    return [Segment("SSO_H P5 tail", SSO_TAIL, "r1"),
            Segment("cell barcode", "X" * 6, "cbc", placeholder=True, feature=CELL_BARCODE),
            Segment("H8 random-primer copy", H8, placeholder=True),
            Segment("Klenow copy of converted DNA", "X" * 28, placeholder=True)]


def _single(label: str, segs: list[Segment]) -> Scene:
    sc = Scene()
    sc.strand("strand", segs, label=label)
    return sc


def _primed_duplex() -> Scene:
    """Barcoded copy on its converted template; the SSO tail and barcode stay unpaired."""
    sc = _single("barcoded second strand", _barcoded_segments())
    sc.anneal("template", [Segment("converted template", "x" * 28, placeholder=True),
                           Segment("H8 binding site", "D" * 8, placeholder=True)],
              to="strand", pair=("converted template", "Klenow copy of converted DNA"),
              label="original converted strand")
    return sc


def _ligated_strand() -> Scene:
    sc = _single("SPLAT-ligated single strand",
                 [*_barcoded_segments(), Segment("ss1b SPLAT adapter", SS1B, "r2")])
    sc.junction("strand", "Klenow copy of converted DNA", "ss1b SPLAT adapter",
                "splinted ligation")
    return sc


def _ss(label: str) -> Scene:
    return _single(label, [Segment(label, "X" * 28, placeholder=True)])


SPEC = get("scsplat", definition=DEFINITION,
           final_override=(FINAL_LIBRARY, SEQ_PRIMERS), oligos_fn=oligos)


def protocol_graph():
    d = DEFINITION
    src = PUBLISHED + ", Methods scWGBS library preparation and sequencing"
    graph = ProtocolGraph(d.title, authority=src)
    cell = graph.input("cell", "FACS-sorted single cell or nucleus in 2 µl Zymo M-digestion buffer with proteinase K",
                       MaterialKind.BIOLOGICAL, state=Construct([Segment("single-cell genomic DNA", "X" * 28, placeholder=True)], name="genomic DNA"), primary=True)
    bisulfite = graph.input("bisulfite", "15 µl EZ Methylation Direct conversion reagent (manufacturer's program)",
                            MaterialKind.CONDITION)
    converted_raw = graph.product("converted_raw", "bisulfite-treated DNA in conversion reagent",
                                  MaterialKind.MOLECULAR, state=_ss("bisulfite-converted DNA"))
    graph.react("convert", "Lyse and bisulfite-convert DNA in each 384-plate well",
                StepKind.CHEMICAL_CONVERSION,
                ((cell, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)), (converted_raw,))
    column = graph.input("spin_plate", "Zymo 384-well spin plate: binding, M-wash, M-desulphonation 15 min RT, two washes",
                         MaterialKind.SELECTION)
    sso = graph.input("sso_h", "SSO_H_X barcoded strand-synthesis oligo, 500 nM in 7 µl elution buffer",
                      MaterialKind.OLIGO, state=sso_h_1(),
                      source=PUBLISHED + ", Supplementary Table S1 (SSO_H_1 drawn)")
    converted = graph.product("converted", "desulphonated single-stranded converted DNA eluted with SSO_H_X",
                              MaterialKind.MOLECULAR, state=_ss("single-stranded converted DNA"))
    column_waste = graph.product("column_waste", "spin-plate flow-through, desulphonation buffer and washes",
                                 MaterialKind.WASTE)
    graph.react("desulphonate", "Desulphonate and wash on the 384 spin plate; elute with SSO_H_X",
                StepKind.CLEANUP,
                ((converted_raw, InputRole.TEMPLATE), (column, InputRole.SELECTION)),
                (converted, column_waste))
    graph.discard(column_waste, "removed by centrifugation through the spin plate")
    klenow = graph.input("klenow", "Klenow exo− (2.5 U), 2× Blue buffer, 1 mM dNTP; 95 °C 2 min then ice, "
                         "4 °C 5 min, 25 °C 5 min, 37 °C 60 min", MaterialKind.CATALYST)
    barcoded = graph.product("barcoded", "SSO_H-primed barcoded second strand on the converted template",
                             MaterialKind.MOLECULAR, state=_primed_duplex())
    graph.react("strand_synthesis", "Random-prime second-strand synthesis with SSO_H_X and Klenow exo−",
                StepKind.EXTENSION,
                ((converted, InputRole.TEMPLATE), (sso, InputRole.STRUCTURAL),
                 (klenow, InputRole.CATALYST)), (barcoded,))
    exo = graph.input("exonuclease_i", "Exonuclease I (0.5 µl), 37 °C 30 min", MaterialKind.CATALYST)
    digested = graph.product("exo_treated", "barcoded duplex after excess-oligo digestion",
                             MaterialKind.MOLECULAR, state=barcoded.state)
    spent = graph.product("digested_sso", "Exo I-degraded excess SSO_H_X", MaterialKind.WASTE)
    graph.react("exo_digest", "Digest excess strand-synthesis oligos with Exonuclease I",
                StepKind.DIGESTION,
                ((barcoded, InputRole.TEMPLATE), (exo, InputRole.CATALYST)), (digested, spent))
    graph.discard(spent, "digestion products removed by the pooled 1.2× SPRI purification")
    spri1 = graph.input("spri_pool", "AMPure XP 1.2×; elute 10 µl 10 mM Tris–HCl pH 8.5",
                        MaterialKind.SELECTION)
    pooled = graph.product("pooled", "pooled barcoded reactions from up to 32 wells",
                           MaterialKind.MOLECULAR, state=barcoded.state)
    pool_waste = graph.product("pool_spri_waste", "1.2× SPRI supernatant and washes", MaterialKind.WASTE)
    graph.react("pool_purify", "Pool up to 32 single-cell reactions and SPRI-purify",
                StepKind.CLEANUP,
                ((digested, InputRole.TEMPLATE), (spri1, InputRole.SELECTION)), (pooled, pool_waste))
    graph.discard(pool_waste, "removed with the bead supernatant")
    ssb = graph.input("et_ssb", "95 °C 2 min with 50 ng ET-SSB, then ice", MaterialKind.CONDITION)
    single = graph.product("single_strand", "released barcoded single strand", MaterialKind.MOLECULAR,
                           state=_single("barcoded second strand", _barcoded_segments()))
    original = graph.product("original_strand", "released original converted strand (no SSO_H tail)",
                             MaterialKind.MOLECULAR, state=_ss("original converted strand"))
    graph.react("denature", "Heat-denature to separate the second strand from the original strand",
                StepKind.STRAND_SEPARATION,
                ((pooled, InputRole.TEMPLATE), (ssb, InputRole.CONDITION)), (single, original))
    graph.mark_terminal(Terminal("original_strand", Disposition.UNRESOLVED,
                                 "Carried into the bulk SPLAT reaction (not separated); without the "
                                 "SSO_H P5 tail it cannot form the two-primer PCR product. Its "
                                 "ligation and fate are not modelled."))
    adapter = graph.input("splat_adapter", "SPLAT 3′ adapter ss1a/ss1b, 5 µM final",
                          MaterialKind.ADAPTER, state=splat_adapter(),
                          source=PUBLISHED + ", Methods and Supplementary Table S1")
    ligase = graph.input("t4_ligase", "T4 DNA ligase HC (15 U), T4 PNK (10 U), 1× ligase buffer, PEG4000; 20 °C 1 h",
                         MaterialKind.CATALYST)
    ligated = graph.product("ligated", "SPLAT-ligated barcoded single strand", MaterialKind.MOLECULAR,
                            state=_ligated_strand())
    graph.react("splat", "SPLAT: splinted ligation of ss1b to the 3′ end in bulk",
                StepKind.LIGATION,
                ((single, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), geometry=adapter_geometry())
    spri2 = graph.input("spri_ligation", "SPRI 1:1; elute 10 µl", MaterialKind.SELECTION)
    clean = graph.product("ligated_clean", "purified SPLAT-ligated strands", MaterialKind.MOLECULAR,
                          state=ligated.state)
    lig_waste = graph.product("ligation_spri_waste", "1:1 SPRI supernatant, free adapter fraction and washes",
                              MaterialKind.WASTE)
    graph.react("ligation_cleanup", "SPRI-purify the ligation (1:1)", StepKind.CLEANUP,
                ((ligated, InputRole.TEMPLATE), (spri2, InputRole.SELECTION)), (clean, lig_waste))
    graph.discard(lig_waste, "removed with the bead supernatant")
    forward = graph.input("pcr_i5", "NEBNext i5 dual-index primer (sequence not printed)",
                          MaterialKind.OLIGO, inferred=True, source=NEBNEXT)
    reverse = graph.input("pcr_i7", "NEBNext i7 dual-index primer (sequence not printed)",
                          MaterialKind.OLIGO, inferred=True, source=NEBNEXT)
    kapa = graph.input("kapa", "KAPA HiFi 2× PCR mastermix, 12–16 cycles", MaterialKind.CATALYST)
    amplified = graph.product("amplified", "dual-indexed scSPLAT pool library", MaterialKind.LIBRARY,
                              state=FINAL_LIBRARY, inferred=True, source=NEBNEXT)
    graph.react("index_pcr", "INFERRED — Amplify with NEBNext dual-index primers (arms not printed)",
                StepKind.PCR,
                ((clean, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="The paper names NEBNext dual index primers without printing them."))
    spri3 = graph.input("spri_pcr", "two 0.8× SPRI purifications; elute 10 µl 10 mM Tris pH 8.0",
                        MaterialKind.SELECTION)
    library = graph.product("library", "purified scSPLAT sequencing library", MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY, inferred=True, source=NEBNEXT)
    pcr_waste = graph.product("pcr_spri_waste", "0.8× SPRI supernatants (primers, short products)",
                              MaterialKind.WASTE)
    graph.react("pcr_cleanup", "Purify the PCR twice with 0.8× SPRI", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (spri3, InputRole.SELECTION)), (library, pcr_waste))
    graph.discard(pcr_waste, "removed with the bead supernatants")
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=SEQ_PRIMERS, acquired_roles=("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2"),
        # PE150 is stated; index-read cycles are not. Pools carrying unique dual
        # indices were combined before sequencing, so the index reads demultiplex pools.
        read_lengths={"Read 1": 150, "Read 2": 150},
        read_uses={"Read 1": "extract_then_include", "Index 1 (i7)": "demultiplex",
                   "Index 2 (i5)": "demultiplex", "Read 2": "include"}))
    graph.declare_reviewed_scope(tuple(graph.materials) + tuple(graph.steps), src)
    graph.annotate(SourceDiscrepancy(
        "splat_adapter", "Published SPLAT-adapter partner name differs",
        (SourceClaim("Scientific Reports Methods, adapter preparation (PMC8986790 JATS line 131)",
                     "The SPLAT adapter was prepared from oligo ss1a and ss2b (Supplementary Table S1)."),
         SourceClaim("Same paragraph, printed oligo list; Supplementary Table S1",
                     "The second oligo is printed as ss1b: 5′Phos/AGATCGGAAGAGCACACGTC/3′AmMo; no ss2b entry exists.")),
        "The model draws the printed ss1b sequence as the partner of ss1a; the printed sequences pair over 20 nt.",
        "Both statements are in the published article; only the name differs.",
        "Whether ss2b is a typo for ss1b or a distinct reagent is not stated; the original SPLAT paper (Raine 2017, NAR) has not been checked."))
    return graph.require_publishable()


SPEC.protocol_graph = protocol_graph
