"""Molecular model for Oxford Nanopore SQK-NBD114.24 native barcoding."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from endprep import EndPreparedDuplex, dA_tailed_scene
from protocol_model import role
from protocol_graph import SequencingPlan, Terminal, Disposition

# ONT CHTD_500 v1 revAV section15, source-qualified locally.
LA_TOP = "CCTGTACTTCGTTCAGTTACGTATTGCT"
NATIVE_ADAPTER_BOTTOM = "ACGTAACTGAACGAAGTACAGG"
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role as molecular_role

TITLE = "Oxford Nanopore Native Barcoding V14 — SQK-NBD114.24"
NOTES = "01_native-barcoding-v14.html"
SOURCE = ('Commercial defining source: Oxford Nanopore '
          '<a href="https://nanoporetech.com/document/ligation-sequencing-gdna-native-barcoding-v14-sqk-nbd114-24">'
          'NBE_9169_v114 rev X (15 Jul 2026)</a>; barcode flanks and Native Adapter strands from the '
          '<a href="https://nanoporetech.com/document/chemistry-technical-document">Chemistry Technical '
          'Document</a> CHTD_500 v1 rev AV, sections 14–15.')
SUMMARY = ("PCR-free sample multiplexing: repair and dA-tail each DNA sample, clean up and "
           "normalise to equimolar mass, ligate a dT-tailed native barcode, stop with EDTA, "
           "pool and clean up, then ligate the motor-bearing Native Adapter through the "
           "barcode adapter's cohesive end.")
CAVEAT = ("The selected genomic-DNA route is NBE_9169_v114_revX_15Jul2026. "
          "CHTD_500_v1_revAV_24July2026 publishes barcode analysis contexts and Native "
          "Adapter DNA strands. Analysis contexts do not establish a complete physical "
          "barcode adapter. The native graph retains the known dT-tail ligation and "
          "cohesive-hook operation but leaves the complete assembly geometry and "
          "motor association unassigned. Separate rows denote distinct populations or "
          "unplaced components, never a covalent DNA–protein polymer.")

NATIVE_BARCODES = {
    "NB01": "CACAAAGACACCGACAACTTTCTT", "NB02": "ACAGACGACTACAAACGGAATCGA",
    "NB03": "CCTGGTAACTGGGACACAAGACTC", "NB04": "TAGGGAAACACGATAGAATCCGAA",
    "NB05": "AAGGTTACACAAACCCTGGACAAG", "NB06": "GACTACTTTCTGCCTTTGCGAGAA",
    "NB07": "AAGGATTCATTCCCACGGTAACAC", "NB08": "ACGTAACTTGGTTTGTTCCCTGAA",
    "NB09": "AACCAAGACTCGCTGTGCCTAGTT", "NB10": "GAGAGGACAAAGGTTTCAACGCTT",
    "NB11": "TCCATTCCCTCCGATAGATGAAAC", "NB12": "TCCGATTCTGCTTCTTTCTACCTG",
    "NB13": "AGAACGACTTCCATACTCGTGTGA", "NB14": "AACGAGTCTCTTGGGACCCATAGA",
    "NB15": "AGGTCTACCTCGCTAACACCACTG", "NB16": "CGTCAACTGACAGTGGTTCGTACT",
    "NB17": "ACCCTCCAGGAAAGTACCTCTGAT", "NB18": "CCAAACCCAACAACCTAGATAGGC",
    "NB19": "GTTCCTCGTGCAGTGTCAAGAGAT", "NB20": "TTGCGTCCTGTTACGAGAACTCAT",
    "NB21": "GAGCCTCTCATTGTCCGTTCTCTA", "NB22": "ACCACTGCCATGTATCAAAGTACG",
    "NB23": "CTTACTACCCAGTGAACCTCCTCG", "NB24": "GCATAGTTCTGCATGATGGGTTAG",
}
NB01 = NATIVE_BARCODES["NB01"]
# CHTD_500 v1 rev AV section 14 ("Barcode flanking sequences"), written 5'->3' around
# the forward barcode and its reverse complement respectively.
NB_FLANK_FORWARD = ("AAGGTTAA", "CAGCACCT")
NB_FLANK_REVERSE = ("GGTGCTG", "TTAACCTTAGCAAT")
INSERT = Construct([Segment("genomic DNA", "X" *42, placeholder=True)], name="genomic DNA")
# The product guide specifies dA-tailing;5′ phosphate identity is not established
# by that phrase alone. Native barcode joining requires prepared ligatable ends.
END_PREP = EndPreparedDuplex(INSERT, five_prime_phosphate=True, three_prime_overhang="A")

def prepared_population():
    sc=dA_tailed_scene(END_PREP)
    sc.strand("control population", [role("independent end-prepared DCS DNA; same sample barcode", inferred=True)], label="DCS control population")
    sc.note("top", "Ligatable5′phosphates are inferred from nativeTA ligation; exact native enzyme components unassigned")
    return sc
BARCODE_FEATURE = feature("native_barcode", "sample_index", "whitelist",
                          whitelist="Oxford Nanopore NB01–NB24 published barcode set")
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = ("A motor-associated Native Adapter captures DNA at a nanopore; "
                     "there is no synthesis sequencing primer. The selected sample barcode "
                     "is used for demultiplexing; read length and barcode-cycle placement "
                     "are not assigned by a short-read run recipe.")
MOTOR_LIMITATION = ("Native Adapter DNA strands and barcode analysis contexts are published; "
                    "the full motor-bearing assembly, native hook bases in the complete "
                    "barcode adapter, and strand-entry polarity remain unassigned in this scope.")
END_PREP_PROGRAM = ((20,5),(65,5))
BARCODE_LIGATION_MIN =20
NATIVE_ADAPTER_LIGATION_MIN =20
BARCODE_END = reported_adapter_end(
    "Native Barcode adapter", kind="cohesive", polarity="3-prime", terminal_sequence="T",
    substrate_scenes=(dA_tailed_scene(END_PREP),),
    statement="The kit24 product guide specifies a dT-tailed barcode adapter for dA-tailed DNA; full barcode-adapter strands are not assigned here.")
NATIVE_HOOK = AdapterLigationGeometry.unavailable(
    "The product guide specifies the barcode-adapter cohesive hook and the protocol "
    "names Native Adapter plus QuickT4 ligase. The complete barcode-adapter strands "
    "are not established by analysis-context flanks; their native hook alignment remains unassigned.")


def native_adapter_scene():
    core = revcomp(NATIVE_ADAPTER_BOTTOM)
    if not LA_TOP.startswith(core):
        raise ValueError("source Native Adapter strands do not share an antiparallel core")
    sc = Scene()
    sc.strand("NA top", [Segment("NA paired core", core),
                        Segment("NA3′ unpaired extension", LA_TOP[len(core):])], label="NA top")
    sc.anneal("NA bottom", [Segment("NA bottom", NATIVE_ADAPTER_BOTTOM)], to="NA top",
              pair=("NA bottom", "NA paired core"), label="NA bottom")
    sc.note("NA top", "Disclosed DNA component; full native-barcode hook alignment unassigned")
    return sc


def barcoded_scene(*, pooled=False, adapted=False):
    sc=Scene()
    sc.strand("DNA", [role("barcoded genomic DNA; full end architecture unassigned", inferred=True)],
              label="barcoded DNA population")
    sc.strand("barcode component", [Segment("barcode identity; placement unassigned", "N"*24,
              "cbc", placeholder=True, inferred=True, feature=BARCODE_FEATURE)],
              label="NB01–NB24 identity, not a separate free molecule")
    sc.note("barcode component", "Unplaced feature of the barcode adapter; one independently assigned barcode per sample")
    sc.strand("DCS control", [role("independently barcoded DCS control DNA; same sample barcode", inferred=True)], label="independent control molecules")
    if pooled:
        sc.strand("other samples", [role("other independently barcoded genomic DNA populations", inferred=True)],
                  label="other pool members; membership unassigned")
    if adapted:
        sc.strand("motor association", [role("Native Adapter plus associated motor; connection unassigned", inferred=True)],
                  label="associated component; no covalent placement")
    return sc


def oligos():
    yield oligo("Native Adapter top", [Segment("Native Adapter top", LA_TOP, "r2")])
    yield oligo("Native Adapter bottom", [Segment("Native Adapter bottom", NATIVE_ADAPTER_BOTTOM, "r2")])
    for name, forward in NATIVE_BARCODES.items():
        yield oligo(f"{name} forward", [Segment("native barcode", forward, "cbc",
                                                feature=BARCODE_FEATURE)])
        yield oligo(f"{name} reverse", [Segment("native barcode reverse", revcomp(forward),
                                                "cbc", feature=BARCODE_FEATURE)])


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    dna = g.input("dna", "genomic DNA sample (400 ng for >4 barcodes, 1000 ng for ≤4)",
                  MaterialKind.MOLECULAR,
                  state=INSERT,
                  primary=True)
    control = g.input("dna_control", "DNA Control Sample (DCS),1 µl per end-prep recipe",
                      MaterialKind.MOLECULAR,
                      state=molecule("control DNA", (molecular_role("independent control DNA"),)))
    repair = g.input("repair_mix", "NEBNext FFPE DNA Repair Mix and Ultra II End-prep enzymes",
                     MaterialKind.CATALYST)
    tailed = g.product("tailed", "repaired dA-tailed genomic DNA",
                       MaterialKind.MOLECULAR, state=prepared_population())
    g.react("end_prep", "Repair DNA ends and add one 3-prime dA (20 °C 5 min, 65 °C 5 min)",
            StepKind.END_REPAIR,
            ((dna, InputRole.TEMPLATE), (control, InputRole.SPIKE_IN), (repair, InputRole.CATALYST)), (tailed,))
    beads_1 = g.input("end_prep_beads", "15 µl AMPure XP Beads, 80% ethanol washes, water elution",
                      MaterialKind.SELECTION)
    tailed_clean = g.product("tailed_clean", "bead-purified dA-tailed DNA",
                             MaterialKind.MOLECULAR, state=prepared_population())
    end_prep_waste = g.product("end_prep_waste", "end-prep supernatant and ethanol washes",
                               MaterialKind.WASTE)
    g.react("end_prep_cleanup", "Clean up each end-prep reaction on AMPure XP beads",
            StepKind.CLEANUP,
            ((tailed, InputRole.TEMPLATE), (beads_1, InputRole.SELECTION)),
            (tailed_clean, end_prep_waste))
    g.discard(end_prep_waste, "unbound end-prep reaction components removed with the supernatant")
    normalise = g.input("equimolar_rule", "Qubit quantification; equimolar mass per sample",
                        MaterialKind.CONDITION)
    normalised = g.product("normalised", "equimolar end-prepped DNA for barcode ligation",
                           MaterialKind.MOLECULAR, state=prepared_population())
    leftover = g.product("leftover", "end-prepped DNA/control not taken into barcode ligation",
                         MaterialKind.MOLECULAR, state=prepared_population())
    g.react("normalise", "Take an equimolar mass of each sample into barcode ligation",
            StepKind.SAMPLING,
            ((tailed_clean, InputRole.TEMPLATE), (normalise, InputRole.CONDITION)),
            (normalised, leftover), changes_molecule=False)
    g.mark_terminal(Terminal(leftover.id, Disposition.UNRESOLVED,
                             "Not sampled into barcode ligation; source does not prescribe disposal."))
    barcode = g.input("barcode_adapter", "dT-tailed Native Barcode adapter (one of NB01–NB24)",
                      MaterialKind.ADAPTER, inferred=True,
                      source=SOURCE + " (barcode identities and analysis contexts published)")
    ligase_1 = g.input("barcode_ligase", "NEB Blunt/TA Ligase Master Mix", MaterialKind.CATALYST)
    barcoded = g.product("barcoded", "native-barcoded DNA sample",
                         MaterialKind.MOLECULAR, state=barcoded_scene(), inferred=True,
                         source=SOURCE + " (complete physical barcode adapter unassigned)")
    barcode_residual = g.product("barcode_residual", "unligated/partially barcoded DNA/control and off-products",
                                 MaterialKind.MOLECULAR, state=normalised.state, inferred=True)
    g.react("barcode", "Ligate one Native Barcode adapter per sample (20 min, room temperature)",
            StepKind.LIGATION,
            ((normalised, InputRole.TEMPLATE), (barcode, InputRole.STRUCTURAL),
             (ligase_1, InputRole.CATALYST)), (barcoded, barcode_residual), geometry=BARCODE_END)
    g.mark_terminal(Terminal(barcode_residual.id, Disposition.UNRESOLVED,
        "Barcode ligation is not quantitative adaptation; later pool/bead recovery of this class is unassigned."))
    edta = g.input("edta", "EDTA (4 µl per well)", MaterialKind.CONDITION)
    stopped = g.product("stopped", "EDTA-stopped barcode ligation",
                        MaterialKind.MOLECULAR, state=barcoded.state, inferred=True,
                        source=SOURCE + " (complete physical barcode adapter unassigned)")
    g.react("stop", "Stop each barcode ligation with EDTA", StepKind.QUENCHING,
            ((barcoded, InputRole.TEMPLATE), (edta, InputRole.CONDITION)), (stopped,),
            changes_molecule=False)
    others = g.input("other_barcoded", "other EDTA-stopped barcode ligations "
                     "(different NB barcodes, processed in parallel)", MaterialKind.MOLECULAR)
    pooled = g.product("pooled", "pool of all barcoded samples",
                       MaterialKind.MOLECULAR, state=barcoded_scene(pooled=True), inferred=True,
                       source=SOURCE + " (distinct sample populations)")
    g.react("pool", "Pool the entire volume of every barcoded sample", StepKind.POOL,
            ((stopped, InputRole.TEMPLATE), (others, InputRole.TEMPLATE)), (pooled,),
            changes_molecule=False)
    beads_2 = g.input("pool_beads", "0.4X AMPure XP Beads, 700 µl SFB washes, "
                      "water elution at 37 °C", MaterialKind.SELECTION)
    pooled_clean = g.product("pooled_clean", "bead-purified barcoded pool",
                             MaterialKind.MOLECULAR, state=barcoded_scene(pooled=True), inferred=True,
                             source=SOURCE + " (distinct sample populations)")
    pool_waste = g.product("pool_waste", "pool supernatant (unligated barcode adapters, "
                           "ligase, EDTA) and SFB washes", MaterialKind.WASTE)
    g.react("pool_cleanup", "Clean up the pooled barcode ligations", StepKind.CLEANUP,
            ((pooled, InputRole.TEMPLATE), (beads_2, InputRole.SELECTION)),
            (pooled_clean, pool_waste))
    g.discard(pool_waste, "unbound material removed with the bead supernatant and washes")
    native = g.input("native_adapter", "motor-loaded Native Adapter",
                     MaterialKind.ADAPTER, state=native_adapter_scene(), inferred=True,
                     source=SOURCE + " (NA strands published; motor association unassigned)")
    ligase_2 = g.input("native_ligase", "NEBNext Quick Ligation Reaction Buffer and Quick T4 DNA Ligase",
                       MaterialKind.CATALYST)
    library = g.product("library", "motor-adapted native-barcode library",
                        MaterialKind.LIBRARY, state=barcoded_scene(pooled=True, adapted=True), inferred=True,
                        source=SOURCE + " (complete motor association unassigned)")
    hook_residual = g.product("hook_residual", "unligated/partially Native-Adapter-adapted pool and off-products",
                              MaterialKind.MOLECULAR, state=pooled_clean.state, inferred=True)
    g.react("native_adapter_ligation", "Attach Native Adapter to the pooled library "
            "(20 min, room temperature)", StepKind.LIGATION,
            ((pooled_clean, InputRole.TEMPLATE), (native, InputRole.STRUCTURAL),
             (ligase_2, InputRole.CATALYST)), (library,hook_residual), geometry=NATIVE_HOOK)
    g.mark_terminal(Terminal(hook_residual.id, Disposition.UNRESOLVED,
        "Cohesive-hook ligation is not physical purification; later bead partition of unadapted classes is unassigned."))
    beads_3 = g.input("library_beads", "20 µl AMPure XP Beads, 125 µl LFB or SFB washes "
                      "(not ethanol), Elution Buffer at 37 °C", MaterialKind.SELECTION)
    library_clean = g.product("library_clean", "bead-purified native-barcode library",
                              MaterialKind.LIBRARY, state=barcoded_scene(pooled=True, adapted=True), inferred=True,
                              source=SOURCE + " (complete motor association unassigned)")
    library_waste = g.product("library_waste", "ligation supernatant (unligated Native "
                              "Adapter) and LFB/SFB washes", MaterialKind.WASTE)
    g.react("library_cleanup", "Clean up the adapted library; LFB enriches ≥3 kb, SFB keeps all sizes",
            StepKind.CLEANUP,
            ((library, InputRole.TEMPLATE), (beads_3, InputRole.SELECTION)),
            (library_clean, library_waste))
    g.discard(library_waste, "unbound material removed with the bead supernatant and washes")
    g.sequence(library_clean, SequencingPlan(
        acquired_roles=("Nanopore strand",), read_uses={"Nanopore strand":"include"},
        limitation=MOTOR_LIMITATION))
    g.declare_reviewed_scope(("dna", "dna_control", "repair_mix", "end_prep", "end_prep_beads",
                              "end_prep_cleanup", "equimolar_rule", "normalise",
                              "barcode_adapter", "barcode_ligase", "barcode", "edta", "stop",
                              "other_barcoded", "pool", "pool_beads", "pool_cleanup",
                              "native_adapter", "native_ligase", "native_adapter_ligation",
                              "library_beads", "library_cleanup", "library_clean"), SOURCE)
    return g.require_publishable()
