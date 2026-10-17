"""ASAP-seq molecular model.

Authority: Mimitou et al., Nat. Biotechnol. 39:1246 (2021), PMC8763625 Methods
"Cell fixation and permeabilization for ASAP-seq" and "Transposition and barcoding
for ASAP-seq"; publisher supplementary protocol MOESM1 (ASAP-seq protocol, Sept 2020).
The 10x scATAC v1.1 constant arms come from CG000209 Rev A Appendix p48 (the article
cites Rev D, which was not inspected).  The drawn antibody-tag library is the
TotalSeq-A/BOA branch; TotalSeq-B/BOB and hashtag (D7xx) libraries are not drawn.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

ARTICLE = '<a href="https://doi.org/10.1038/s41587-021-00927-2">Mimitou et al., <i>Nature Biotechnology</i> (2021)</a>'
SOURCE = "Defining source: " + ARTICLE + "."
GUIDE = ('<a href="https://assets.ctfassets.net/an68im79xiti/3SEbiWVegS8it73fsQNXr0/458c7ee408f249e1762233ef2f437c9e/'
         'CG000209_Chromium_NextGEM_SingleCell_ATAC_ReagentKits_v1.1_UserGuide_RevA.pdf">10x CG000209 Rev A, Appendix p. 48</a>')

CELL = feature("asap_cell", "cell_barcode", "whitelist", whitelist="10x Chromium scATAC v1.1 gel-bead barcode")
UBI = feature("asap_ubi", "umi", "random")
ANTIBODY = feature("asap_antibody", "feature_barcode", "whitelist",
                   whitelist="TotalSeq-A antibody barcodes (MOESM3 Supplementary Table 1)")
ATAC_I7 = feature("asap_atac_i7", "sample_index", "unknown")
ADT_I7 = feature("asap_adt_i7", "sample_index", "unknown")

# Article Methods / MOESM1 "Oligos": BOA, 3'-inverted-dT blocked.
BOA_HANDLE = nx.READ1_PRIMER                      # TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
BOA_UBI = "NNNNNNNNNV"
BOA_POLY_T = "T" * 30
BOA = BOA_HANDLE + BOA_UBI + BOA_POLY_T
BOB = BOA_HANDLE + "TTGCTAGGACCGGCCTTAAAGC"
# RPI-x primer 3' arm; Methods: its 3' part is the TruSeq Small RNA handle present in TSA tags.
RPI_ARM = "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
TSA_HANDLE = RPI_ARM[12:]                         # CCTTGGCACCCGAGAATTCCA
D7_ARM = "GTGACTGGAGTTCAGACGTGTGC"
# CG000209 Rev A Appendix p48.
BEAD_ARM = nx.READ1_PRIMER[:14]                   # TCGTCGGCAGCGTC
SI_PCR_B_FORWARD = il.P5[:22]                     # AATGATACGGCGACCACCGAGA
SI_REVERSE_ARM = revcomp(nx.INDEX1_PRIMER)[:15]   # GTCTCGTGGGCTCGG
assert BEAD_ARM == "TCGTCGGCAGCGTC" and SI_PCR_B_FORWARD == "AATGATACGGCGACCACCGAGA"
assert SI_REVERSE_ARM == "GTCTCGTGGGCTCGG" and TSA_HANDLE == "CCTTGGCACCCGAGAATTCCA"


def _cell_barcode():
    return Segment("cell barcode", "N" * 16, "cbc", placeholder=True, feature=CELL)


def bead_primer():
    return Construct([Segment("P5", il.P5, "p5"), _cell_barcode(),
                      Segment("partial Read 1N", BEAD_ARM, "r1")], name="scATAC v1.1 gel-bead primer")


def boa_state():
    return Construct([Segment("Read 1N handle", BOA_HANDLE, "r1"),
                      Segment("UBI (N9V)", BOA_UBI, "umi", placeholder=True, feature=UBI),
                      Segment("poly(dT) tag-annealing arm", BOA_POLY_T)], name="BOA /3InvdT/")


def atac_precursor():
    """CG000209 Rev A 'Linear Amplification DNA Product'."""
    return Construct([Segment("P5", il.P5, "p5"), _cell_barcode(),
                      Segment("Read 1N", nx.READ1_PRIMER, "r1"),
                      Segment("accessible genomic DNA", "X" * 42, placeholder=True),
                      Segment("Read 2N reverse complement", nx.INDEX1_PRIMER, "r2")],
                     name="ASAP-seq cell-barcoded ATAC fragment")


def atac_library():
    """CG000209 Rev A 'Sample Index PCR Product'."""
    return Construct([*atac_precursor().segments,
                      Segment("i7 reverse complement", "N" * 8, "cbc", placeholder=True, feature=ATAC_I7),
                      Segment("P7 reverse complement", il.P7_RC, "p7")], name="ASAP-seq ATAC library")


def _adt_core():
    """Bead-barcoded TSA tag: bead arm, BOA-derived handle/UBI/dT, tag barcode, handle."""
    return [Segment("P5", il.P5, "p5"), _cell_barcode(),
            Segment("Read 1N (from bead and BOA)", nx.READ1_PRIMER, "r1"),
            Segment("UBI (N9V)", BOA_UBI, "umi", placeholder=True, feature=UBI),
            Segment("BOA poly(dT)", BOA_POLY_T),
            Segment("antibody barcode reverse complement", "N" * 15, "cbc", placeholder=True,
                    feature=ANTIBODY, inferred=True),
            Segment("TSA small-RNA handle reverse complement", revcomp(TSA_HANDLE), "r2")]


def adt_precursor():
    return Construct(_adt_core(), name="ASAP-seq bead-barcoded TotalSeq-A tag")


def adt_library():
    core = _adt_core()[:-1]
    return Construct([*core, Segment("RPI-x arm reverse complement", revcomp(RPI_ARM), "r2"),
                      Segment("i7 reverse complement", "N" * 8, "cbc", placeholder=True, feature=ADT_I7),
                      Segment("P7 reverse complement", il.P7_RC, "p7")],
                     name="ASAP-seq TotalSeq-A antibody-tag library")


ATAC = atac_library()
ADT = adt_library()
ATAC_P = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["I2"], sp.NEXTERA["R2"])
ADT_P = (sp.NEXTERA["R1"], sp.NEXTERA["I2"])
for _library, _primers in ((ATAC, ATAC_P), (ADT, ADT_P)):
    _problems = sp.verify(_library, _primers, required_roles=tuple(p.role for p in _primers))
    if _problems:
        raise ValueError("; ".join(_problems))

TITLE = "ASAP-seq — ATAC and surface-protein profiling"
NOTES = "01_asap-seq.html"
SUMMARY = ("Antibody-stained cells are fixed, permeabilized and tagmented; in the 10x scATAC GEM a "
           "3′-blocked bridge oligo templates extension of each antibody tag so the bead primer barcodes "
           "tags together with chromatin fragments, which are then recovered as separate libraries.")
CAVEAT = ("The article and its supplementary protocol print the BOA/BOB bridges (3′ inverted dT), RPI-x and "
          "D7xx index-primer families, and the reaction modifications. The scATAC constant arms, linear "
          "amplification product, SI-PCR Primer B and final ATAC library are transcribed from the historical "
          + GUIDE + "; the article cites CG000209 Rev D (and CG000168 Rev D for v1), not inspected, so revision "
          "equivalence is unverified. <strong>INFERRED:</strong> the TotalSeq-A conjugate sequence is not printed "
          "in the inspected sources; the antibody-tag library places the 15-nt MOESM3 barcode directly after the "
          "BOA poly(dT) and before the small-RNA handle named in Methods, omitting the conjugate's unresolved "
          "poly(A) length and any linker. The P5 primer used for tag PCR is named but not printed, and the "
          "acquired read cycles and tag-library Read 2/Index 1 primers are not stated in the inspected sources; "
          "no cycle counts are assigned. TotalSeq-B (BOB) and hashtag (D7xx) libraries are not drawn.")
ADT_LIMITATION = ("Read 2 and Index 1 primers for the small-RNA arm of the TotalSeq-A tag library are not "
                  "printed in the inspected sources, so their binding sites and antibody-barcode cycles are not asserted.")
FINAL_LIBRARIES = (
    ("ATAC library", ATAC, ATAC_P,
     "CG000209 Rev A Sample Index PCR product: the 10x barcode sits between P5 and Read 1N.",
     "Read 1N and Read 2N read the insert; the i5 read acquires the cell barcode."),
    ("Antibody-tag library", ADT, ADT_P,
     "INFERRED — TotalSeq-A/BOA tag library; antibody-barcode adjacency to the BOA poly(dT) is inferred.",
     "Read 1N reads the UBI then poly(T); the i5 read acquires the cell barcode."),
)


def oligos():
    return [
        "<b>BOA bridge, TotalSeq-A (Methods; MOESM1)</b> — 5′-" + BOA_HANDLE + BOA_UBI + BOA_POLY_T + "/3InvdT/-3′",
        "<b>BOB bridge, TotalSeq-B (Methods; MOESM1)</b> — 5′-" + BOB + "/3InvdT/-3′",
        "<b>RPI-x, TotalSeq-A tag index primer (Methods)</b> — 5′-" + il.P7 + "xxxxxxxx" + RPI_ARM + "-3′ (x: sample index)",
        "<b>D7xx_s, TotalSeq-B tag / hashtag index primer (Methods)</b> — 5′-" + il.P7 + "xxxxxxxx" + D7_ARM + "-3′ (x: sample index)",
        "<b>scATAC v1.1 gel-bead primer PN-2000210 (CG000209 Rev A)</b> — 5′-" + il.P5 + "-NNNNNNNNNNNNNNNN-" + BEAD_ARM + "-3′ (N16: cell barcode)",
        "<b>SI-PCR Primer B forward (CG000209 Rev A)</b> — 5′-" + SI_PCR_B_FORWARD + "-3′",
        "<b>Sample Index N reverse (CG000209 Rev A)</b> — 5′-" + il.P7 + "-NNNNNNNN-" + SI_REVERSE_ARM + "-3′",
    ]


def _endpoint(g, prefix, template, final, primers, reviewed, *, forward_state, reverse_state,
              geometry, limitation=""):
    forward = g.input(prefix + "_forward", prefix + " indexing forward primer", MaterialKind.OLIGO,
                      state=forward_state, inferred=forward_state is None)
    reverse = g.input(prefix + "_reverse", prefix + " indexing reverse primer", MaterialKind.OLIGO,
                      state=reverse_state, inferred=reverse_state is None)
    polymerase = g.input(prefix + "_polymerase", prefix + " indexing polymerase", MaterialKind.CATALYST)
    library = g.product(prefix + "_library", final.name, MaterialKind.LIBRARY, state=final)
    g.react(prefix + "_pcr", "Index and amplify the " + prefix + " library", StepKind.PCR,
            ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,), geometry=geometry)
    g.sequence(library, SequencingPlan(library=final, primers=tuple(primers),
                                       acquired_roles=tuple(p.role for p in primers), limitation=limitation))
    reviewed.extend((forward.id, reverse.id, polymerase.id, prefix + "_pcr", library.id))


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    cells = g.input("cells", "single-cell suspension with surface epitopes and chromatin", MaterialKind.BIOLOGICAL,
                    state=molecule("ASAP-seq cells", (role("surface epitopes"), role("chromatin"))), primary=True)
    antibodies = g.input("antibodies", "TotalSeq-A/-B oligo-conjugated antibodies after Fc block", MaterialKind.SELECTION,
                         inferred=True, source=SOURCE + " (conjugate sequences not printed; MOESM3 lists barcodes)")
    stained = g.product("stained", "stained cells retaining antibody-derived tags", MaterialKind.MOLECULAR,
                        state=molecule("stained cells", (role("antibody tags"), role("chromatin"))))
    unbound = g.product("unbound", "unbound antibody conjugates", MaterialKind.WASTE)
    g.react("stain", "Stain cells, then wash three times", StepKind.BINDING,
            ((cells, InputRole.TEMPLATE), (antibodies, InputRole.SELECTION)), (stained, unbound))
    g.discard(unbound, "removed in the three staining-buffer washes")
    fixative = g.input("fixative", "1% formaldehyde 10 min RT, glycine quench", MaterialKind.CONDITION)
    fixed = g.product("fixed", "fixed stained cells", MaterialKind.MOLECULAR, state=stained.state)
    g.react("fix", "Fix stained cells", StepKind.CROSS_LINKING, ((stained, InputRole.TEMPLATE), (fixative, InputRole.CONDITION)),
            (fixed,), changes_molecule=False)
    lysis = g.input("lysis", "LLL (mtDNA retained) or OMNI lysis buffer on ice", MaterialKind.CONDITION)
    permeable = g.product("permeabilized", "permeabilized fixed cells in Diluted Nuclei Buffer", MaterialKind.MOLECULAR,
                          state=stained.state)
    g.react("permeabilize", "Permeabilize fixed cells (optional intracellular staining follows)", StepKind.PERMEABILIZATION,
            ((fixed, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)), (permeable,), changes_molecule=False)
    transposome = g.input("transposome", "10x scATAC Transposition Mix", MaterialKind.ADAPTER, inferred=True,
                          source=GUIDE + " (Read 1N/Read 2N ends printed; loaded-complex state undeclared here)")
    tn5 = g.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = g.product("tagged", "cells with tagmented chromatin and antibody tags", MaterialKind.MOLECULAR,
                       state=molecule("tagmented stained cells", (role("antibody tags"), role("Read 1N/Read 2N tagmented chromatin"))))
    tn5_waste = g.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    g.react("tagment", "Tagment accessible chromatin (10x step 1)", StepKind.TAGMENTATION,
            ((permeable, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)),
            (tagged, tn5_waste))
    g.discard(tn5_waste, "does not retain an amplifiable accessibility fragment")
    droplets = g.input("droplets", "10x scATAC GEM partitioning", MaterialKind.CONDITION)
    occupied = g.product("occupied", "GEMs with cell, bead, barcoding mix and bridge", MaterialKind.MOLECULAR, state=tagged.state)
    failed = g.product("failed", "empty or multiply occupied GEM products", MaterialKind.WASTE)
    g.react("partition", "Partition cells with gel beads and barcoding mix", StepKind.PARTITION,
            ((tagged, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)), (occupied, failed), changes_molecule=False)
    g.discard(failed, "filtered computationally; not a separate library")
    bridge = g.input("bridge", "0.5 µl 1 µM BOA (and/or BOB) 3′-inverted-dT bridge in barcoding mix", MaterialKind.OLIGO,
                     state=boa_state(), source=SOURCE + " (Methods; MOESM1 Oligos)")
    bead = g.input("bead_primer", "scATAC v1.1 gel-bead primer PN-2000210", MaterialKind.OLIGO, state=bead_primer(), source=GUIDE)
    polymerase = g.input("gem_polymerase", "10x barcoding enzyme", MaterialKind.CATALYST)
    barcoded = g.product("barcoded", "cell-barcoded ATAC fragments and bridge-extended antibody tags", MaterialKind.MOLECULAR,
                         state=molecule("barcoded pool", (role("cell-barcoded ATAC fragments"), role("cell-barcoded TSA tags"))))
    g.react("gem_barcoding", "40°C 5 min bridge annealing, 72°C extension, 12-cycle linear barcoding", StepKind.EXTENSION,
            ((occupied, InputRole.TEMPLATE), (bridge, InputRole.STRUCTURAL), (bead, InputRole.STRUCTURAL),
             (polymerase, InputRole.CATALYST)), (barcoded,))
    silane = g.input("silane", "silane bead cleanup, 43.5 µl elution; 3 µl set aside", MaterialKind.SELECTION)
    eluate = g.product("eluate", "silane eluate", MaterialKind.MOLECULAR, state=barcoded.state)
    silane_waste = g.product("silane_waste", "GEM reagents not bound to silane beads", MaterialKind.WASTE)
    g.react("silane_cleanup", "Break GEMs and purify with silane beads", StepKind.CLEANUP,
            ((barcoded, InputRole.TEMPLATE), (silane, InputRole.SELECTION)), (eluate, silane_waste), changes_molecule=False)
    g.discard(silane_waste, "removed with the silane supernatant")
    spri = g.input("spri", "1.2× SPRI; supernatant raised to 2.0× (32 µl beads)", MaterialKind.SELECTION)
    atac = g.product("atac_branch", "1.2× SPRI-bound ATAC fragments", MaterialKind.MOLECULAR, state=atac_precursor())
    adt = g.product("adt_branch", "2.0× supernatant fraction (+3 µl silane eluate) with short tag products",
                    MaterialKind.MOLECULAR, state=adt_precursor())
    g.react("spri_split", "Separate ATAC fragments from short tag products by SPRI", StepKind.SELECTION,
            ((eluate, InputRole.TEMPLATE), (spri, InputRole.SELECTION)), (atac, adt), changes_molecule=False)
    reviewed = ["cells", "antibodies", "stain", "fixative", "fix", "lysis", "permeabilize", "transposome", "tn5",
                "tagment", "droplets", "partition", "bridge", "bead_primer", "gem_polymerase", "gem_barcoding",
                "silane", "silane_cleanup", "spri", "spri_split"]
    si_forward = PCRPrimer("SI-PCR Primer B forward", SI_PCR_B_FORWARD, anneal_nt=len(SI_PCR_B_FORWARD))
    si_reverse = PCRPrimer("Sample Index N reverse", il.P7 + "N" * 8 + SI_REVERSE_ARM, anneal_nt=len(SI_REVERSE_ARM))
    _endpoint(g, "atac", atac, ATAC, ATAC_P, reviewed,
              forward_state=Construct([Segment("SI-PCR Primer B forward", SI_PCR_B_FORWARD, "p5")]),
              reverse_state=Construct([Segment("P7", il.P7, "p7"),
                                       Segment("i7 index", "N" * 8, "cbc", placeholder=True, feature=ATAC_I7),
                                       Segment("partial Read 2N", SI_REVERSE_ARM, "r2")]),
              geometry=PCRGeometry.exact(atac.state, si_forward, si_reverse,
                                         statement="CG000209 Rev A Appendix p48 prints SI-PCR Primer B and Sample Index N primers."))
    _endpoint(g, "adt", adt, ADT, ADT_P, reviewed, forward_state=None,
              reverse_state=Construct([Segment("P7", il.P7, "p7"),
                                       Segment("i7 index", "N" * 8, "cbc", placeholder=True, feature=ADT_I7),
                                       Segment("RPI-x arm", RPI_ARM, "r2")]),
              geometry=PCRGeometry.unavailable(adt.state, reason=(
                  "Methods/MOESM1 print RPI-x but name the forward primer only as 'P5'; KAPA, 14–16 cycles.")),
              limitation=ADT_LIMITATION)
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
