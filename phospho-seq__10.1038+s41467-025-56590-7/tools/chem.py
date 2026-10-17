"""Phospho-seq (10x scATAC route) declaration, workflow and reaction graph.

Source authority: Blair et al., Nat. Commun. 16, 1346 (2025),
doi:10.1038/s41467-025-56590-7, open article PMC11794950 (JATS Methods "Antibody
conjugation and pooling", "Phospho-seq", "Sequencing and alignment").  The bioRxiv
preprint 10.1101/2023.03.27.534442 (which named the directory until 2026-10-17) is
superseded for active claims.  BOA/BOB sequences are not reprinted by Phospho-seq; they are taken from the
cited ASAP-seq method (Mimitou et al. 2021, PMC8763625 JATS line 234), which the
Phospho-seq Methods name as the protocol it follows for the bridge step.
"""
from chemdraw import Construct, Segment, feature
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role
from multimodal_protocols import ANTIBODY, CELL, UMI, Definition, ident

DEFINITION = Definition(
    "Phospho-seq — intracellular protein dynamics with chromatin accessibility",
    "10.1038/s41467-025-56590-7",
    "Blair et al., Nature Communications (2025); preprint bioRxiv 10.1101/2023.03.27.534442",
    "DNA-tagged intracellular antibodies are bridged into a 10x scATAC capture reaction; tag "
    "and ATAC branches are separated after GEM barcoding while retaining cell identity.",
    "The ATAC library and antibody/hash tag libraries share a droplet cell barcode; antibody "
    "barcodes identify total or phosphorylated targets.", "01_phospho_seq.html")

# Published PMC11794950 Methods (Antibody conjugation; Phospho-seq), 5'->3'.
PHOSPHO_HASH_A_TAIL = "A" * 28            # then *A*A: 30 A in total as printed
PHOSPHO_BLOCKER_N = "N" * 29              # then /3ddC/ ("30-nt 3' blocked" oligo)
PHOSPHO_P5 = "AATGATACGGCGACCACCGA"
PHOSPHO_P7_HEAD = "CAAGCAGAAGACGGCATACGAGAT"
PHOSPHO_P7_TSA_TAIL = "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
PHOSPHO_P7_TSB_TAIL = "GTGACTGGAGTTCAGACGTGTGC"
PHOSPHO_P7_TSA = PHOSPHO_P7_HEAD + "XXXXXXXX" + PHOSPHO_P7_TSA_TAIL
PHOSPHO_P7_TSB = PHOSPHO_P7_HEAD + "XXXXXXXX" + PHOSPHO_P7_TSB_TAIL
# Phospho-seq-multi additives / P5 (not drawn in the graph; listed only).
PHOSPHO_HTO_ADDITIVE = "CCTTGGCACCCGAGAATTCC"
PHOSPHO_ADT_ADDITIVE = "GTGACTGGAGTTCAGACGTGTGCTC"
PHOSPHO_MULTIOME_P5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC"
# ASAP-seq PMC8763625 JATS line 234 (inherited, see module docstring); both /3InvdT/.
PHOSPHO_BOA = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNVTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
PHOSPHO_BOB = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTTGCTAGGACCGGCCTTAAAGC"

# BOB is printed by the cited ASAP-seq method, not by the Phospho-seq bundle; verified
# exactly at sources/asap-seq__10.1038+s41587-021-00927-2 PMC8763625 JATS line 234 and
# MOESM1 line 256 (2026-10-17 pass).
SOURCE_OLIGO_EXCEPTIONS = {
    "bridge_bob": "Inherited from ASAP-seq (PMC8763625 JATS line 234); Phospho-seq names "
                  "BOB without reprinting it.",
}

HASH = ident("hash tag", "sample_hash", "sample_index")
I7 = role("8-nt i7 index", tag="cbc", length_bp=8,
          molecular_feature=feature("sample_index", "sample_index", "unknown"))

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = (f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a> '
          '(PMC11794950 Methods). BOA/BOB from the cited ASAP-seq method (PMC8763625).')
SUMMARY = DEFINITION.summary
CAVEAT = ("The published Methods print the home-conjugated TotalSeq-A-style hash and "
          "TotalSeq-B-style antibody tags, the 3′-blocked blocking oligo, the P5 primer and two "
          "indexed P7 families, and the 34/8/16/34 run. BOA and BOB are named but not reprinted; "
          "their sequences are taken from the cited ASAP-seq method. The 10x scATAC v1.1 gel-bead "
          "oligo and ATAC library steps follow the vendor kit (not printed), so the cell-barcode "
          "end and the ATAC branch are drawn at role level. The free TotalSeq-A spike-in, SSB "
          "pre-binding and the Phospho-seq-multi (10x Multiome) variant are described in the notes, "
          "not drawn. The author step-by-step PDFs were unreachable in this pass.")

TAG_LIBRARY = Construct([
    role("P5 end (P5 primer)"), CELL,
    role("Read 1 / bridge-oligo capture region", inferred=True),
    UMI, ANTIBODY, role("TotalSeq-A or -B tag handle"), I7, role("P7 end"),
], name="Phospho-seq antibody/hash tag library roles")
ATAC_LIBRARY = Construct([role("P5 end"), CELL, role("Tn5 accessible-chromatin insert"),
                          role("i7 sample index (10x kit)"), role("P7 end")],
                         name="Phospho-seq ATAC library roles")
FINAL_LIBRARY = TAG_LIBRARY
SEQ_PRIMERS = ()
FINAL_CAPTION = ("INFERRED — tag-library role order: the bridge-capture junction and 10x "
                 "bead-oligo end are not printed. TotalSeq-B tags carry UMIs; TotalSeq-A hash "
                 "tags do not.")
READ_LENGTHS = {"Read 1": 34, "Index 1 (i7)": 8, "Index 2 (i5)": 16, "Read 2": 34}
ACQUIRED_ROLES = ("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2")
SEQUENCING_INTRO = ("Libraries pooled at 25% ADT, 10% HTO, 65% ATAC; NextSeq 550, 75-cycle High "
                    "Output kit: 34 bp Read 1, 8 bp i7, 16 bp i5, 34 bp Read 2.")
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = ("Run lengths are published; the 10x scATAC bead oligo that supplies "
                          "the cell-barcode end is not printed, so primer landing is not asserted.")


def oligos():
    return (
        "<b>TotalSeq-A-style hashing tag</b> — 5′-/5AmMC12/CCTTGGCACCCGAGAATTCCA"
        "[15-nt antibody barcode]B" + PHOSPHO_HASH_A_TAIL + "*A*A-3′",
        "<b>TotalSeq-B-style antibody tag</b> — "
        "5′-/5AmMC12/GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT[10-nt UMI]"
        "[15-nt antibody barcode][9-nt UMI]GCTTTAAGGCCGGTCCTAGC*A*A-3′",
        "<b>30-nt blocking oligo</b> — 5′-" + PHOSPHO_BLOCKER_N + "/3ddC/-3′",
        "<b>Bridge oligo A (BOA, ASAP-seq)</b> — 5′-" + PHOSPHO_BOA + "/3InvdT/-3′",
        "<b>Bridge oligo B (BOB, ASAP-seq)</b> — 5′-" + PHOSPHO_BOB + "/3InvdT/-3′",
        "<b>P5 tag-library primer</b> — 5′-" + PHOSPHO_P5 + "-3′",
        "<b>TotalSeq-A P7 index-primer family</b> — 5′-" + PHOSPHO_P7_TSA + "-3′",
        "<b>TotalSeq-B P7 index-primer family</b> — 5′-" + PHOSPHO_P7_TSB + "-3′",
        "<b>Multi: HTO additive</b> — 5′-CCTTGGCACCCGAGAATT*C*C-3′",
        "<b>Multi: ADT additive</b> — 5′-GTGACTGGAGTTCAGACGTGTGC*T*C-3′",
        "<b>Multi: SI-PCR P5</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGC*T*C-3′",
    )


def workflow():
    return chain(molecule("fixed permeabilized cell", (
        role("accessible chromatin"), role("intracellular phosphoproteins"))), (
        ("Block (with hashing antibodies and 3′-blocked N29 oligo), then stain with DNA-tagged antibodies",
         molecule("tagged cell", (ANTIBODY, role("phosphoprotein epitope"), HASH)),
         "Click-conjugated (mTz/TCO) antibodies carry TotalSeq-A or TotalSeq-B style tags."),
        ("Tagment; add BOA and BOB to the barcoding mix; 5 min 40 °C before GEM incubation",
         molecule("capturable branches", (role("tagmented accessible DNA"),
                                           role("bridge-bound antibody tag"))),
         "Bridge oligos make the antibody tags compatible with the gel-bead capture sequence."),
        ("Save the first-SPRI supernatant for tags; ATAC DNA stays on the beads",
         molecule("matched Phospho-seq libraries", (CELL, role("ATAC insert"), UMI,
                                                     ANTIBODY, role("phosphoprotein tag"))),
         "Tags are amplified from the cleaned supernatant plus 3 µl set-aside eluate."),
    ))


def _oligo(name, seq):
    return Construct([Segment(name, seq)], name=name)


def _boa():
    head, tail = PHOSPHO_BOA.split("NNNNNNNNNV")
    return Construct([Segment("BOA Read 1/ME", head),
                      Segment("N9V", "NNNNNNNNNV", placeholder=True, length_bp=10),
                      Segment("T30", tail)], name="BOA")


def _p7(name, tail):
    return Construct([Segment("P7 head", PHOSPHO_P7_HEAD),
                      Segment("i7", "X" * 8, "cbc", placeholder=True, length_bp=8,
                              feature=feature("sample_index", "sample_index", "unknown")),
                      Segment(name + " tail", tail)], name=name)


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    source = f"{d.citation}; doi:{d.doi}; PMC11794950 Methods"
    vendor = ("10x Chromium Single Cell ATAC v1.1 user guide (kit 100175), named in the Methods; "
              "not read in this pass")
    asap = "ASAP-seq, Mimitou et al. 2021, PMC8763625 JATS line 234 (cited bridge method)"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cells", "fixed, permeabilized cells", MaterialKind.BIOLOGICAL,
                        state=molecule("fixed permeabilized cell", (
                            role("accessible chromatin"), role("intracellular phosphoproteins"))),
                        primary=True)
    hashing = graph.input("hash_antibodies", "TotalSeq-A-style CD298/B2M hashing antibodies (1 µg)",
                          MaterialKind.SELECTION)
    blocker = graph.input("blocking_oligo", "100 µg 3′-blocked N29/3ddC/ single-stranded oligo",
                          MaterialKind.OLIGO,
                          state=Construct([Segment("N29", PHOSPHO_BLOCKER_N, placeholder=True,
                                                   length_bp=29)], name="blocking oligo"))
    hashed = graph.product("hashed_cells", "blocked, hashed cells", MaterialKind.MOLECULAR,
                           state=molecule("hashed cell", (HASH, role("accessible chromatin"))))
    unbound_hash = graph.product("unbound_hash", "unbound hashing antibody", MaterialKind.WASTE)
    graph.react("block_hash", "Block 30 min RT with hashing antibodies and blocking oligo; wash",
                StepKind.BINDING,
                ((cells, InputRole.TEMPLATE), (hashing, InputRole.SELECTION),
                 (blocker, InputRole.STRUCTURAL)), (hashed, unbound_hash))
    graph.discard(unbound_hash, "removed by centrifugation and wash buffer #2")
    antibodies = graph.input("antibodies", "SSB-preincubated DNA-tagged antibody pool (0.25–0.5 µg each)",
                             MaterialKind.SELECTION)
    stained = graph.product("stained_cells", "antibody-tagged hashed cells",
                            MaterialKind.MOLECULAR,
                            state=molecule("tagged cell", (ANTIBODY, HASH,
                                                           role("accessible chromatin"))))
    unbound = graph.product("unbound_antibodies", "unbound antibody tags", MaterialKind.WASTE)
    graph.react("stain", "Pool hashed cells; stain 1 h RT; wash twice; 1x nuclei buffer; 40 µm filter",
                StepKind.BINDING,
                ((hashed, InputRole.TEMPLATE), (antibodies, InputRole.SELECTION)),
                (stained, unbound))
    graph.discard(unbound, "removed by two washes in wash buffer #2")
    transposome = graph.input("transposome", "10x ATAC enzyme (3 µl) in ATAC buffer B (7 µl)",
                              MaterialKind.ADAPTER, source=vendor)
    tn5 = graph.input("tn5", "10x ATAC Tn5 transposase", MaterialKind.CATALYST, source=vendor)
    tagmented = graph.product("tagmented_cells", "tagged cells with tagmented accessible chromatin",
                              MaterialKind.MOLECULAR, state=stained.state)
    tn5_waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment 30,000 cells (10x scATAC v1.1)", StepKind.TAGMENTATION,
                ((stained, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagmented, tn5_waste), source=vendor)
    graph.discard(tn5_waste, "does not retain an amplifiable accessibility fragment")
    boa = graph.input("bridge_boa", "Bridge oligo A /3InvdT/ (0.5 µl 1 µM)", MaterialKind.OLIGO,
                      state=_boa(), source=asap)
    bob = graph.input("bridge_bob", "Bridge oligo B /3InvdT/ (0.5 µl 1 µM)", MaterialKind.OLIGO,
                      state=_oligo("BOB", PHOSPHO_BOB), source=asap)
    droplets = graph.input("droplets", "10x scATAC GEM generation with barcoding mix",
                           MaterialKind.CONDITION, source=vendor)
    gems = graph.product("gems", "GEMs with tagmented cells, bead oligos and bridge oligos",
                         MaterialKind.MOLECULAR,
                         state=molecule("capturable branches", (role("tagmented accessible DNA"),
                                                                role("bridge-bound antibody tag"))))
    failed = graph.product("failed_partitions", "empty or ambiguous GEMs", MaterialKind.WASTE)
    graph.react("partition", "Add BOA/BOB to the barcoding reaction and generate GEMs",
                StepKind.PARTITION,
                ((tagmented, InputRole.TEMPLATE), (boa, InputRole.STRUCTURAL),
                 (bob, InputRole.STRUCTURAL), (droplets, InputRole.CONDITION)),
                (gems, failed), changes_molecule=False)
    graph.discard(failed, "does not retain an interpretable cell barcode")
    bead = graph.input("bead_oligo", "10x scATAC gel-bead barcode oligo (not printed)",
                       MaterialKind.OLIGO, source=vendor)
    gem_pol = graph.input("gem_polymerase", "10x GEM barcoding enzyme", MaterialKind.CATALYST,
                          source=vendor)
    barcoded = graph.product("barcoded_products", "cell-barcoded ATAC and antibody-tag products",
                             MaterialKind.MOLECULAR,
                             state=molecule("barcoded Phospho-seq products", (
                                 CELL, role("ATAC insert"), UMI, ANTIBODY,
                                 role("phosphoprotein or hash tag"))))
    graph.react("gem_barcoding", "5 min 40 °C hybridization, then GEM incubation (bridge-mediated tag capture)",
                StepKind.EXTENSION,
                ((gems, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
                 (gem_pol, InputRole.CATALYST)), (barcoded,),
                reason="Bridge-oligo capture geometry with the unprinted bead oligo is not "
                       "derived here.")
    silane = graph.input("silane_beads", "silane beads; 43.5 µl Elution Solution I",
                         MaterialKind.SELECTION, source=vendor)
    eluate = graph.product("silane_eluate", "40 µl eluate for SPRI cleanup", MaterialKind.MOLECULAR,
                           state=barcoded.state)
    aside = graph.product("set_aside", "3 µl eluate set aside for tag PCR", MaterialKind.MOLECULAR,
                          state=barcoded.state)
    silane_waste = graph.product("silane_waste", "broken-emulsion supernatant", MaterialKind.WASTE)
    graph.react("silane_cleanup", "Silane cleanup (step 3.1o); elute 43.5 µl; keep 40 µl + 3 µl aside",
                StepKind.CLEANUP, ((barcoded, InputRole.TEMPLATE), (silane, InputRole.SELECTION)),
                (eluate, aside, silane_waste))
    graph.discard(silane_waste, "removed from bead-bound barcoded products")
    spri = graph.input("first_spri", "first SPRI cleanup (step 3.2d)", MaterialKind.SELECTION,
                       source=vendor)
    atac = graph.product("atac_branch", "bead-bound ATAC DNA", MaterialKind.MOLECULAR,
                         state=barcoded.state)
    sup = graph.product("tag_supernatant", "saved supernatant with hash and antibody tags",
                        MaterialKind.MOLECULAR, state=barcoded.state)
    graph.react("first_spri_split", "First SPRI: ATAC DNA on beads; supernatant saved (not discarded)",
                StepKind.CLEANUP, ((eluate, InputRole.TEMPLATE), (spri, InputRole.SELECTION)),
                (atac, sup))

    # ATAC branch: continues per the 10x kit.
    atac_fwd = graph.input("atac_forward", "10x SI-PCR primer", MaterialKind.OLIGO, source=vendor)
    atac_rev = graph.input("atac_reverse", "10x i7 sample-index primer", MaterialKind.OLIGO, source=vendor)
    atac_pol = graph.input("atac_polymerase", "10x amplification mix", MaterialKind.CATALYST, source=vendor)
    atac_lib = graph.product("atac_library", "Phospho-seq ATAC library", MaterialKind.LIBRARY,
                             state=ATAC_LIBRARY)
    graph.react("atac_pcr", "ATAC sample-index PCR as written in the 10x kit", StepKind.PCR,
                ((atac, InputRole.TEMPLATE), (atac_fwd, InputRole.STRUCTURAL),
                 (atac_rev, InputRole.STRUCTURAL), (atac_pol, InputRole.CATALYST)), (atac_lib,),
                source=vendor,
                geometry=PCRGeometry.unavailable(
                    ATAC_LIBRARY, reason="ATAC library proceeds 'normally as written' in the "
                                         "10x scATAC v1.1 kit; primers are not printed."))
    graph.sequence(atac_lib, SequencingPlan(limitation=SEQUENCING_UNAVAILABLE,
                                            acquired_roles=ACQUIRED_ROLES,
                                            read_lengths=READ_LENGTHS))

    # Tag branch: published Methods.
    spri32 = graph.input("tag_spri_beads", "32 µl SPRI beads; 80% EtOH; elute 42 µl EB",
                         MaterialKind.SELECTION)
    tags = graph.product("tag_clean", "cleaned tag fraction", MaterialKind.MOLECULAR,
                         state=barcoded.state)
    tag_waste = graph.product("tag_spri_waste", "unbound supernatant", MaterialKind.WASTE)
    graph.react("tag_spri", "Bind supernatant to 32 µl SPRI 5 min; wash; elute 42 µl EB",
                StepKind.CLEANUP, ((sup, InputRole.TEMPLATE), (spri32, InputRole.SELECTION)),
                (tags, tag_waste))
    graph.discard(tag_waste, "removed from the SPRI-bound tag fraction")
    pooled = graph.product("tag_input", "42 µl tag eluate + 3 µl set-aside eluate",
                           MaterialKind.MOLECULAR, state=barcoded.state)
    graph.react("combine_tag_input", "Combine with the 3 µl set-aside eluate", StepKind.POOL,
                ((tags, InputRole.TEMPLATE), (aside, InputRole.TEMPLATE)), (pooled,),
                changes_molecule=False)
    tsa_in = graph.product("tsa_input", "22.5 µl input for TotalSeq-A PCR", MaterialKind.MOLECULAR,
                           state=barcoded.state)
    tsb_in = graph.product("tsb_input", "22.5 µl input for TotalSeq-B PCR", MaterialKind.MOLECULAR,
                           state=barcoded.state)
    graph.react("split_tag_pcr", "Split into TotalSeq-A and TotalSeq-B PCRs", StepKind.SPLIT,
                ((pooled, InputRole.TEMPLATE),), (tsa_in, tsb_in), changes_molecule=False)
    for key, label, tail, inp in (("tsa", "TotalSeq-A (hash)", PHOSPHO_P7_TSA_TAIL, tsa_in),
                                  ("tsb", "TotalSeq-B (antibody)", PHOSPHO_P7_TSB_TAIL, tsb_in)):
        p5 = graph.input(key + "_p5", "P5 primer (2.5 µl 10 µM)", MaterialKind.OLIGO,
                         state=_oligo("P5 primer", PHOSPHO_P5))
        p7 = graph.input(key + "_p7", label + " indexed P7 primer", MaterialKind.OLIGO,
                         state=_p7(label + " P7", tail))
        kapa = graph.input(key + "_kapa", "2× KAPA HotStart ReadyMix (50 µl)", MaterialKind.CATALYST)
        amp = graph.product(key + "_amplicon", label + " tag amplicon", MaterialKind.MOLECULAR,
                            state=TAG_LIBRARY)
        graph.react(key + "_pcr", label + " PCR: 95 °C 3 min; 14–18 × (95 °C 20 s, 60 °C 30 s, 72 °C 20 s); 72 °C 5 min",
                    StepKind.PCR,
                    ((inp, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                     (p7, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)), (amp,),
                    geometry=PCRGeometry.unavailable(
                        TAG_LIBRARY,
                        reason="P5 and P7 are published, but the P5-side template end comes "
                               "from the unprinted 10x bead oligo."))
        beads = graph.input(key + "_spri", "1.6× SPRI beads; elute 20 µl EB", MaterialKind.SELECTION)
        lib = graph.product(key + "_library", "Phospho-seq " + label + " library",
                            MaterialKind.LIBRARY, state=TAG_LIBRARY)
        waste = graph.product(key + "_waste", "unbound primers and reaction components",
                              MaterialKind.WASTE)
        graph.react(key + "_cleanup", "1.6× SPRI; elute 20 µl EB", StepKind.CLEANUP,
                    ((amp, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (lib, waste))
        graph.discard(waste, "left in the SPRI supernatant")
        graph.sequence(lib, SequencingPlan(limitation=SEQUENCING_UNAVAILABLE,
                                           acquired_roles=ACQUIRED_ROLES,
                                           read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope(list(graph.materials) + list(graph.steps), source)
    return graph.require_publishable()
