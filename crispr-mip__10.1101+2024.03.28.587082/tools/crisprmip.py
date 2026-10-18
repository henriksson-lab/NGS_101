"""
CRISPR-MIP: sgRNA quantification by molecular inversion probe.

Selinger M, Yakovenko I, Nazir I, Henriksson J. "CRISPR-MIP replaces PCR and reveals GC
and oversampling bias in pooled CRISPR screens." bioRxiv 2024.03.28.587082.
doi:10.1101/2024.03.28.587082

The idea: replace the readout PCR of a pooled CRISPR screen with a padlock capture. The
probe carries a UMI, so each captured genomic molecule is labelled before any
amplification -- which turns the readout from "count reads" into "count molecules", and
removes the PCR-derived GC and sequencing-depth biases the paper reports.
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from crispr import SCAFFOLD_V1, U6_PLUS1
from padlock import Capture, Padlock, capture
from plasmid import Plasmid, amplify
from topology import ClosedCircleTopology, PadlockCaptureTopology

# --------------------------------------------------------------- probe design
# All sequences verbatim from Table S2 of the preprint (supplementary media-2.xlsx).
EXT_ARM = "GTGGAAAGGACGAAACACC"        # 19 nt, published Tm 53.8 C -- ends at the U6 +1
LIG_ARM = "AGCTAGGTCTTGAAAGGAGTGGG"    # 23 nt, published Tm 58.3 C -- past the scaffold

UMI = "N" * 13                                             # 13 nt, published
# The backbone's two sequencing-primer sites are the Illumina TruSeq ones, referenced from
# lib/illumina.py rather than retyped. Table S2 (ref/TableS2_primers_and_probes.xlsx; the
# selftest checks every probe against it) shows the per-probe variable stretch is 8 nt, so
# the constant Read 1 site is TruSeq Read 1 less its 5' ACAC.
READ2_SITE = revcomp(il.TRUSEQ_READ2)                      # 34 nt
READ1_SITE = il.TRUSEQ_READ1[4:]                           # 29 nt
INDEX_LEN = 8

# The nine probes differ ONLY in this 8-nt i7 index (Table S2, MIP_probe_1..9).
PROBE_INDICES = ("ATCACGAC", "ACAGTGGT", "CAGATCCA", "ACAAACGG", "ACCCAGCA",
                 "AACCCCTC", "CCCAACCT", "CACCACAC", "GAAACCCA")

# Illustrative non-G spacer for drawings and tests; not an actual library member.
EXAMPLE_SPACER = "ATCGATCGATCGATCGATCG"

CAPTURED_SPAN = 112      # published, and computed from the real vector
# The 15 bases between the end of the v1 scaffold and the ligation arm.  This is the
# only local vector context needed to construct the captured molecule; it was extracted
# from Addgene #52961 and is cross-checked against that map whenever the optional source
# file is present.  Keeping this small, reviewed landmark in the model lets a clean clone
# build the page without bundling the third-party GenBank file.
CAPTURE_3PRIME_FLANK = "TTTTTTGAATTCGCT"
UMI_LEN = 13
READ1_CYCLES = 60        # published: reads the captured sgRNA
READ2_CYCLES = 15        # published: enough for the 13-nt UMI

# The amplification primers. These are SEPARATE oligos, and -- the key design point --
# they anneal inside the CAPTURED sequence (the tracrRNA/scaffold), not in the probe
# backbone. A probe that annealed but was never extended therefore has no primer sites.
# Their 5' ends are the Illumina P5/P7 flow-cell sequences.
P5_TRACR_TAIL = "GTGCTTTTTTGAATTCGC"          # Table S2: vector just past the scaffold
P7_TRACR_TAIL = revcomp(SCAFFOLD_V1[32:50])   # Table S2: GTTGATAACGGACTAGCC, on the scaffold
P5_TRACR_FWD = il.P5 + P5_TRACR_TAIL
P7_TRACR_REV = il.P7 + P7_TRACR_TAIL


def backbone(index: str = PROBE_INDICES[0]) -> str:
    """The 84-nt non-complementary middle of the probe."""
    return UMI + READ2_SITE + index + READ1_SITE


PROBE_LEN = len(LIG_ARM) + len(backbone()) + len(EXT_ARM)   # 126
BACKBONE_LEN = len(backbone())                              # 84

# The paper's Methods text says the probe is "134 bp". Table S2's actual sequences are
# 126 nt, and decompose exactly as above. We follow Table S2 and flag the discrepancy.
PROBE_LEN_IN_TEXT = 134

UMI_FEATURE = feature("probe_umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "whitelist",
                     whitelist="crispr-mip-nine-probe-index-set")

def probe_construct(index: str = PROBE_INDICES[0]) -> Construct:
    """The probe as ordered, 5'->3': lig arm - UMI - Read 2 site - i7 - Read 1 site - ext arm.

    One definition, used by the page diagram, the selftest and the notes (via mdfacts), so
    a length or an index width cannot be restated anywhere and go stale.
    """
    return Construct([
        Segment("ligation arm", LIG_ARM, "r3"),
        Segment("UMI", UMI, "umi", placeholder=True, feature=UMI_FEATURE),
        Segment("Read 2 site", READ2_SITE, "t7"),
        Segment("i7", index, None, feature=I7_FEATURE),
        Segment("Read 1 site", READ1_SITE, "s5"),
        Segment("extension arm", EXT_ARM, "r1"),
    ], name="CRISPR-MIP probe")


def probe(index: str = PROBE_INDICES[0]) -> Padlock:
    """One of the nine CRISPR-MIP probes, selected by its i7 index."""
    return Padlock(f"CRISPR-MIP-{index}", ext_arm=EXT_ARM, lig_arm=LIG_ARM,
                   backbone=backbone(index), phosphorylated=True)


def capture_fill(spacer: str = EXAMPLE_SPACER) -> str:
    """Sequence copied between the two probe arms for a cloned guide.

    This illustrative cassette adds a G before a non-G spacer; for a G-starting
    spacer the alternative model omits that separately added base. This convention
    is not an observed first-base policy for every guide in the screened library.
    The scaffold and 15-nt vector flank match the deposited lentiCRISPRv2 source.
    """
    spacer = spacer.upper()
    if len(spacer) != 20 or any(b not in "ACGT" for b in spacer):
        raise ValueError("CRISPR-MIP example spacer must be 20 DNA bases")
    return ("" if spacer.startswith("G") else U6_PLUS1) + spacer + SCAFFOLD_V1 + CAPTURE_3PRIME_FLANK


def modeled_capture(spacer: str = EXAMPLE_SPACER,
                    index: str = PROBE_INDICES[0]) -> Capture:
    """Build the capture from the protocol landmarks, independent of a plasmid file."""
    target = Plasmid("lentiCRISPRv2 guide cassette",
                     EXT_ARM + capture_fill(spacer) + LIG_ARM,
                     circular=False, features=[])
    hits = capture(target, probe(index))
    if len(hits) != 1:
        raise ValueError(f"modeled guide cassette gives {len(hits)} CRISPR-MIP captures")
    return hits[0]


def filled_probe_construct() -> Construct:
    """The selected example's extended probe, before its last nick is sealed.

    The capture constructor computes the copied gap from the supplied target.
    This represents one successful candidate, not a capture-efficiency claim.
    """
    captured = modeled_capture()
    product = Construct([
        *probe_construct(),
        Segment("copied guide-cassette gap", captured.fill,
                note="Illustrative guide cassette; not every screened clone"),
    ], name="extended CRISPR-MIP probe before nick sealing")
    origin = PROBE_LEN - len(EXT_ARM)
    if product.top()[origin:] + product.top()[:origin] != captured.circle.seq:
        raise ValueError("extended probe differs from the derived capture circle")
    return product


def library_construct() -> Construct:
    """The inverse-PCR product, segmented from the same captured molecule."""
    captured = modeled_capture()
    fill = captured.fill
    p5_start = fill.index(P5_TRACR_TAIL)
    p7_end = fill.index(revcomp(P7_TRACR_TAIL)) + len(P7_TRACR_TAIL)
    leading = len(fill) - len(EXAMPLE_SPACER) - len(SCAFFOLD_V1) - len(CAPTURE_3PRIME_FLANK)
    product = Construct([
        Segment("Illumina P5", il.P5, "p5"),
        Segment("captured vector after P5 site", fill[p5_start:], "r2"),
        *probe_construct(),
        Segment("illustrative U6 initial base", fill[:leading]),
        Segment("captured sgRNA example", fill[leading:leading + len(EXAMPLE_SPACER)]),
        Segment("captured scaffold through P7 site", fill[leading + len(EXAMPLE_SPACER):p7_end], "r2"),
        Segment("Illumina P7 reverse complement", il.P7_RC, "p7"),
    ], name="CRISPR-MIP library")
    expected = amplify(captured.circle, P5_TRACR_FWD, P7_TRACR_REV,
                       min_anneal=min(len(P5_TRACR_TAIL), len(P7_TRACR_TAIL)))
    if len(expected) != 1 or product.top() != expected[0].seq:
        raise ValueError("drawn CRISPR-MIP library differs from inverse-PCR derivation")
    return product


def _single_strand(parts, label, *, mod5=""):
    scene = Scene()
    scene.strand("strand", list(parts), label=label, mod5=mod5)
    return scene


def protocol_graph():
    """Source-owned capture, enzymatic treatment and inverse-PCR provenance."""
    from circular import CircularizationGeometry, circularize_ssdna
    from pcr import PCRGeometry, PCRPrimer
    from protocol_graph import (Disposition, InputRole, MaterialKind,
                                MolecularRelation, ProtocolGraph, SequencingPlan,
                                StepKind, Terminal)

    source = ("Selinger et al., bioRxiv 2024.03.28.587082 v1, Methods: "
              "Design of CRISPR-MIP probes; CRISPR-MIP sgRNA abundance quantification, "
              "Steps 1–4; Supplementary Table S2 (media-2.xlsx)")
    graph = ProtocolGraph("CRISPR-MIP — source-owned capture and inverse PCR", authority=source)
    target_sequence = EXT_ARM + capture_fill() + LIG_ARM
    target = Construct([Segment("illustrative guide-cassette target", target_sequence)], name="guide-cassette example")
    genomic = graph.input("genomic_dna", "genomic DNA containing the guide cassette (illustrative guide drawn)",
                          MaterialKind.BIOLOGICAL, state=target, primary=True)
    heat = graph.input("denaturation_program", PROTOCOL[0][1], MaterialKind.CONDITION)
    target_strand = graph.product("target_strand", "denatured strand complementary to the probe arms",
                                  MaterialKind.MOLECULAR,
                                  state=_single_strand([Segment("target genomic strand", revcomp(target_sequence))], "illustrative genomic target"))
    other_strand = graph.product("other_strand", "complementary genomic strand retained in the same reaction",
                                 MaterialKind.MOLECULAR,
                                 state=_single_strand(list(target), "complementary genomic strand"))
    graph.react("denature", "Denature genomic DNA; retain both strand outcomes", StepKind.STRAND_SEPARATION,
                ((genomic, InputRole.TEMPLATE), (heat, InputRole.CONDITION)),
                (target_strand, other_strand))
    probe_input = graph.input("mip_probe", "Table S2 MIP_probe_1, 5′-phosphorylated",
                              MaterialKind.OLIGO,
                              state=_single_strand(probe_construct(), "single-stranded MIP probe", mod5="P"))
    hybrid_view = Scene()
    hybrid_view.strand("target", [Segment("target genomic strand", revcomp(target_sequence))],
                       label="target genomic strand (5′→3′)")
    hybrid_view.strand("probe", list(probe_construct()), label="one intact padlock probe (5′→3′)", mod5="P")
    hybrid_view.note("probe", "Separate molecules: the two-arm loop pairing is shown in the checked capture topology below.")
    hybrid = graph.product("hybrid", "probe arms annealed around the guide-cassette gap",
                            MaterialKind.MOLECULAR,
                            state=hybrid_view)
    graph.react("hybridize", "Hybridize overnight at 60 °C after the source ramp", StepKind.HYBRIDIZATION,
                ((target_strand, InputRole.TEMPLATE), (probe_input, InputRole.STRUCTURAL)), (hybrid,))
    polymerase = graph.input("phusion", "4 U Phusion HF DNA polymerase (#F530S)", MaterialKind.CATALYST)
    nucleotides = graph.input("dntps", "0.2 mM dNTPs in the added extension/ligation mix", MaterialKind.MOLECULAR)
    reaction = graph.input("capture_program", PROTOCOL[1][1], MaterialKind.CONDITION)
    filled = filled_probe_construct()
    extended = graph.product("extended_probe", "gap-extended candidate probe, nick not yet sealed",
                              MaterialKind.MOLECULAR,
                              state=_single_strand(filled, "extended candidate probe", mod5="P"),
                              chemical_state=filled)
    graph.react("extend_gap", "Copy the gap from the selected example target (part of the same 60 °C reaction)",
                StepKind.EXTENSION, ((hybrid, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST),
                                    (nucleotides, InputRole.STRUCTURAL), (reaction, InputRole.CONDITION)),
                (extended,))
    ligase = graph.input("ampligase", "5 U Ampligase", MaterialKind.CATALYST)
    closed = circularize_ssdna(filled, five_prime_phosphate=True, three_prime_oh=True)
    circle = graph.product("capture_circle", "closed capture-circle candidate",
                            MaterialKind.MOLECULAR, state=closed, chemical_state=closed)
    linear = graph.product("linear_material", "target DNA and candidate unclosed/unreacted probes",
                            MaterialKind.MOLECULAR)
    graph.react("seal_circle", "Seal the gap-filled probe's remaining nick", StepKind.CIRCULARIZATION,
                ((extended, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (circle, linear),
                geometry=CircularizationGeometry.known("gap-sealing probe closure", statement=
                    "Table S2 supplies the probe 5′ phosphate; Methods Step 2 combines Phusion gap extension and Ampligase closure. The drawn split is conceptual within one reaction."),
                relations=(MolecularRelation.CIRCLE_CLOSURE,))
    exoi = graph.input("exoi", "10 U Exonuclease I (#EN0581)", MaterialKind.CATALYST)
    exoiii = graph.input("exoiii", "50 U Exonuclease III (#M0206L)", MaterialKind.CATALYST)
    exo_program = graph.input("exonuclease_program", "37 °C 45 min; 80 °C 20 min inactivation", MaterialKind.CONDITION)
    retained = graph.product("treated_circles", "closed capture circles after exonuclease treatment",
                              MaterialKind.MOLECULAR, state=closed, chemical_state=closed)
    hydrolysates = graph.product("hydrolysis_products", "linear-DNA hydrolysis products; no physical separation",
                                 MaterialKind.MOLECULAR)
    graph.react("digest_linear", "Treat linear DNA with Exonuclease I and III; heat-inactivate",
                StepKind.DIGESTION, ((circle, InputRole.TEMPLATE), (linear, InputRole.TEMPLATE),
                                    (other_strand, InputRole.TEMPLATE), (exoi, InputRole.CATALYST),
                                    (exoiii, InputRole.CATALYST), (exo_program, InputRole.CONDITION)),
                (retained, hydrolysates))
    graph.mark_terminal(Terminal(hydrolysates.id, Disposition.UNRESOLVED,
        "The whole 31-µl reaction enters PCR. Digestion is not a physical discard; hydrolysates remain in the reaction and their later AMPure partition is not reconstructed."))
    forward = graph.input("p5_tracr", "Table S2 P5_tracrRNA_fwd", MaterialKind.OLIGO, state=P5_TRACR_FWD)
    reverse = graph.input("p7_tracr", "Table S2 P7_tracrRNA_rev", MaterialKind.OLIGO, state=P7_TRACR_REV)
    kapa = graph.input("kapa", "KAPA HiFi HotStart ReadyMix", MaterialKind.CATALYST)
    pcr_program = graph.input("pcr_program", PROTOCOL[3][1], MaterialKind.CONDITION)
    library = library_construct()
    # The captured circle is ssDNA: only P7 can bind initially. P5 binds the
    # new complementary copy. These are conceptual phases of Methods Step 4,
    # not separate incubations or purification steps.
    fill = capture_fill()
    p5_start = fill.index(P5_TRACR_TAIL)
    p5_end = p5_start + len(P5_TRACR_TAIL)
    p7_start = fill.index(revcomp(P7_TRACR_TAIL))
    p7_end = p7_start + len(P7_TRACR_TAIL)
    opened = Construct([
        Segment("between sites", fill[p7_end:p5_start]),
        Segment("P5 site", fill[p5_start:p5_end]),
        Segment("fill end", fill[p5_end:]),
        *probe_construct(),
        Segment("captured gap before P7 site", fill[:p7_start]),
        Segment("P7 site", fill[p7_start:p7_end]),
    ], name="same ssDNA circle opened after the P7 site for coordinates")
    # Rotation is representational; every circle base appears exactly once,
    # and the UMI/index features survive it and the complementary copy.
    offset = PROBE_LEN + p7_end
    if opened.top() != closed.sequence[offset:] + closed.sequence[:offset]:
        raise ValueError("P7-first circle view differs from the captured circle")
    first_copy = Construct([Segment("P7 adapter", il.P7),
                            *complement_segments(list(opened))],
                           name="P7-primed one-traversal copy of the circle")
    # Copying the first P7 product to its adapter end must reproduce the
    # already derived inverse-PCR library, including both flow-cell tails.
    p5_on_copy = first_copy.top().upper().index(revcomp(P5_TRACR_TAIL))
    second_copy = il.P5 + revcomp(first_copy.top()[:p5_on_copy + len(P5_TRACR_TAIL)])
    if second_copy.upper() != library.top().upper():
        raise ValueError("P5 copying of the P7-first product differs from the inverse-PCR library")
    first_view = Scene()
    first_view.strand("circle", list(opened), label="ssDNA circle (opened for coordinates)")
    first_view.anneal("P7 primer", [Segment("P7 adapter", il.P7),
                                  Segment("P7 tail", P7_TRACR_TAIL)],
                     to="circle", pair=("P7 tail", "P7 site"),
                     unpaired=("P7 adapter",))
    first_view.arrow("P7 primer", "first copy; no P5 binding strand exists yet")
    first_product = graph.product("first_p7_copy", "P7-primed complementary copy (one-traversal candidate)",
                                   MaterialKind.MOLECULAR, state=_single_strand(first_copy, "first P7 copy"),
                                   chemical_state=first_copy)
    remaining_circle = graph.product("circle_after_first_copy", "circle retained after template-directed copying",
                                     MaterialKind.MOLECULAR, state=closed, chemical_state=closed)
    graph.react("first_p7_extension", "PCR initial phase: P7 primes the ssDNA circle", StepKind.EXTENSION,
                ((retained, InputRole.TEMPLATE), (reverse, InputRole.STRUCTURAL),
                 (kapa, InputRole.CATALYST), (pcr_program, InputRole.CONDITION)),
                (first_product, remaining_circle))
    # The first-copy primer view is the input geometry, while the material state
    # records the copied strand. The detailed cycle-1 panel shows the same sites.
    from protocol_graph import DetailView, Highlight
    graph.annotate(Highlight("first_p7_extension", "Initial ssDNA-circle priming",
        "Only P7 binds the initial circle. P5 requires the newly copied complementary strand. "
        "The one-traversal copy is the candidate used for the existing cycle-1/cycle-2 diagrams; "
        "longer synthesis and yields are not inferred. These phases share the same PCR incubation."))
    graph.annotate(DetailView("first_p7_extension", "P7 primes the initial single-stranded circle",
        "Primer-binding view of the first copying phase, not a pre-existing duplex.", state=first_view))
    graph.mark_terminal(Terminal(remaining_circle.id, Disposition.UNRESOLVED,
        "The circle is a retained copying template, not consumed by extension. It remains in the PCR mixture; "
        "its later AMPure partition is not reconstructed."))
    second_view = Scene()
    second_view.strand("first copy", list(first_copy), label="P7-primed complementary copy")
    second_view.anneal("P5 primer", [Segment("P5 adapter", il.P5), Segment("P5 tail", P5_TRACR_TAIL)],
                      to="first copy", pair=("P5 tail", "P5 site'"), unpaired=("P5 adapter",))
    second_view.arrow("P5 primer", "copy to the P7 adapter end")
    p5_copy = graph.product("first_p5_copy", "first bounded P5-bearing strand; its new P5 tail is not yet copied",
                            MaterialKind.MOLECULAR, state=_single_strand(library, "first P5-primed strand"),
                            chemical_state=library)
    retained_first_copy = graph.product("retained_first_p7_copy", "long first P7 copying template retained in the PCR mixture",
                                        MaterialKind.MOLECULAR, state=_single_strand(first_copy, "retained initial P7 copy"))
    graph.react("first_p5_extension", "PCR initial phase: P5 copies the P7-primed strand", StepKind.EXTENSION,
                ((first_product, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (kapa, InputRole.CATALYST), (pcr_program, InputRole.CONDITION)), (p5_copy, retained_first_copy))
    graph.mark_terminal(Terminal(retained_first_copy.id, Disposition.UNRESOLVED,
        "Copying does not consume the first P7 template. Its later amplification and cleanup partition are not enumerated."))
    first_duplex = graph.product("first_duplex", "bounded P5/P7 duplex after the new P5 tail has also been copied",
                                 MaterialKind.MOLECULAR, state=library)
    graph.react("recopy_p5_strand", "PCR initial phase: P7 copies the P5-bearing strand, including its newly added tail",
                StepKind.EXTENSION, ((p5_copy, InputRole.TEMPLATE), (reverse, InputRole.STRUCTURAL),
                                     (kapa, InputRole.CATALYST), (pcr_program, InputRole.CONDITION)), (first_duplex,))
    graph.annotate(Highlight("recopy_p5_strand", "Both adapter complements now exist",
        "P5 extension adds a new P5 tail to one strand. P7 must copy that strand before its full P5-primer complement "
        "exists on the opposite strand. Exact duplex PCR geometry below applies after this initial copying phase, "
        "within the same published PCR incubation."))
    graph.annotate(DetailView("first_p5_extension", "P5 primes the newly copied strand",
        "The complementary primer site exists after the first P7 copy. This priming view is not the first-copy product state.", state=second_view))
    pcr = graph.product("pcr_library", "inverse-PCR library before AMPure cleanup",
                         MaterialKind.LIBRARY, state=library)
    graph.react("inverse_pcr", "Later cycles amplify the P5/P7 duplex; total cycles set by qPCR",
                StepKind.PCR, ((first_duplex, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                               (reverse, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST),
                               (pcr_program, InputRole.CONDITION)), (pcr,),
                geometry=PCRGeometry.exact(first_duplex.state,
                    PCRPrimer("P5_tracrRNA_fwd", P5_TRACR_FWD, len(P5_TRACR_TAIL)),
                    PCRPrimer("P7_tracrRNA_rev", P7_TRACR_REV, len(P7_TRACR_TAIL)), statement=
                    "Exact inward geometry applies after P7-first copying, P5 extension and recopying of the newly added P5 tail. "
                    "It does not place a P5 primer on the initial ssDNA circle or validate copying efficiency."))
    beads = graph.input("ampure", "AMPure XP post-PCR cleanup", MaterialKind.SELECTION)
    final = graph.product("library", "AMPure-purified CRISPR-MIP library", MaterialKind.LIBRARY, state=library)
    washes = graph.product("ampure_waste", "AMPure supernatant and washes", MaterialKind.WASTE)
    graph.react("cleanup", "AMPure cleanup after PCR; Qubit and Bioanalyzer 2100 QC", StepKind.CLEANUP,
                ((pcr, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (final, washes),
                relations=(MolecularRelation.STRAND_ORDER,))
    graph.discard(washes, "Supernatant and washes are removed during post-PCR bead cleanup; no ratio-specific retention is claimed.")
    graph.sequence(final, SequencingPlan(library=library,
        primers=tuple(p for p in SEQ_PRIMERS if p.role in ("Read 1", "Index 1 (i7)", "Read 2")),
        acquired_roles=("Read 1", "Index 1 (i7)", "Read 2"),
        read_lengths={"Read 1": READ1_CYCLES, "Read 2": READ2_CYCLES},
        limitation="The model represents the Methods 60/15 read profile and the Table S2 8-base i7 sequence. Index-cycling configuration is outside this model; no i5 is assigned.",
        read_uses={"Read 1": "include", "Index 1 (i7)": "demultiplex", "Read 2": "extract"}))
    graph.declare_reviewed_scope((
        "genomic_dna", "denaturation_program", "denature", "mip_probe", "hybridize",
        "phusion", "dntps", "capture_program", "extend_gap", "ampligase", "seal_circle",
        "exoi", "exoiii", "exonuclease_program", "digest_linear", "p5_tracr", "p7_tracr",
        "kapa", "pcr_program", "first_p7_extension", "first_p5_extension", "recopy_p5_strand", "inverse_pcr", "ampure", "cleanup", "library",
    ), source)
    return graph.require_publishable()


# The target strand is the actual antiparallel genomic strand, written 5′→3′. The
# topology constructor locates both sites from the arm sequences and rejects overlap,
# incompatible pairing or a crossing layout.
CAPTURE_TOPOLOGY = PadlockCaptureTopology(
    "CRISPR-MIP probe hybridized around the guide-cassette gap",
    revcomp(EXT_ARM + capture_fill() + LIG_ARM),
    LIG_ARM,
    backbone(),
    EXT_ARM,
)
CLOSED_CAPTURE_TOPOLOGY = ClosedCircleTopology(
    "gap-filled CRISPR-MIP capture circle",
    modeled_capture().circle.seq,
    "var(--c-cbc)",
)


# ------------------------------------------------------------------ protocol
# Conditions from the preprint's four steps (bioRxiv v1 Methods, "CRISPR-MIP sgRNA
# abundance quantification", extracted-text lines 471-493).
PROTOCOL = (
    ("1. denature & hybridise",
     "20 ul reaction: 1-10 ug gDNA, 0.2-0.002 uM probe, 1.5x Ampligase buffer. 94 C 5 min, "
     "then ramp to 60 C at -0.1 C/s, hybridise overnight at 60 C."),
    ("2. extend & ligate",
     "Add 10 ul pre-heated (60 C) mix of 5 U Ampligase + 4 U Phusion HF + 0.2 mM dNTPs "
     "to the 20 ul hybridisation. 60 C, 1 h."),
    ("3. exonuclease",
     "37 C, 45 min with 10 U Exonuclease I (removes un-circularised probe) and 50 U "
     "Exonuclease III (removes genomic DNA). Inactivate 80 C, 20 min."),
    ("4. PCR",
     "Whole 31 ul reaction into one 100 ul KAPA HiFi HotStart ReadyMix PCR with "
     "P5_tracrRNA_fwd + P7_tracrRNA_rev: 95 C 3 min; n cycles of 98 C 20 s, 60 C 15 s, "
     "72 C 30 s (n set by qPCR per screen); 72 C 30 s. The P5/P7 primer sites lie WITHIN "
     "the captured sequence, so only probes that were both extended and ligated amplify."),
    ("5. cleanup & QC",
     "AMPure XP bead cleanup after the PCR; Qubit and Bioanalyzer 2100. The bead ratio is "
     "not stated; this cleanup follows the selection, it is not one of its steps."),
)

# Why the method is specific, in one place: three independent filters, none of which is a
# size selection or a cleanup.
SELECTION_STEPS = (
    "both arms must anneal to the same molecule, a fixed distance apart",
    "the polymerase must cross the gap and the ligase must close the circle",
    "exonuclease I/III destroy everything still linear -- probe, gDNA, failed captures",
    "the PCR primer sites sit inside the captured region, not the backbone",
)


# ------------------------------------------------------------ sequencing primers
# Sequences come from lib/ via seqprimers; where each lands, and what it reads first, is
# computed from the final library on the page. The library has no i5 index: P5 runs
# straight into captured vector sequence.
SEQ_PRIMERS = [
    replace(sp.TRUSEQ["R1"], note="Its 5' ACAC sits over the last bases of the 8-nt i7, "
            "so how many of those 4 pair varies by probe; the 3' 29 nt pair on all nine."),
    sp.TRUSEQ["I1"],
    sp.mismatching(sp.TRUSEQ["I2"],
                   "its 3' ...GTGT needs the ACAC that starts a full TruSeq Read 1 site; "
                   "the probe has only the 3' 29 nt of that site, preceded by the i7's last "
                   "bases, so the primer's 3' end is mismatched -- and there is no i5 to read"),
    sp.custom("Index 2 (i5)", "P5 flow-cell oligo (forward-strand i5 workflow)", il.P5,
              'Illumina "Indexed Sequencing Overview Guide" #15057455',
              "No i5 in this library: these cycles read constant captured-vector sequence."),
    sp.TRUSEQ["R2"],
]
