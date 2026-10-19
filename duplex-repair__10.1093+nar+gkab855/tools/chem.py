"""Protocol-owned source, oligos, workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Workflow, feature
from endprep import ligation_product_scene
from protocol_model import molecule, role
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from accuracy_protocols import _workflow

# Supplementary Figures S9-S10 legends: in those capillary-electrophoresis assays the NGS
# adapters add 39-40 or 37-38 bp to the exposed 3' and 5' ends because "unique molecular
# indices can be either 3 or 4 base pairs".  The legends do not say these are the custom
# IDT dual-index adapters used for the sequenced libraries (inferred, not stated).
# Neither the UMI composition (random or a fixed set) nor the position of the UMI or the
# dual sample indexes within the custom adapters is published, so each feature marks the
# whole adapter role, unsized, with unknown encoding.
_UMI_NOTE = ("UMI within this custom IDT dual-index adapter; position, bases and composition "
             "unpublished. Supplementary Figures S9-S10 legends give 3- or 4-bp UMIs for the "
             "adapters in those CE assays; that these are the library adapters is inferred")


def _umi(side):
    return feature(f"duplex_umi_{side}", "umi", "unknown", note=_UMI_NOTE)


def workflow() -> Workflow:
    damaged = molecule("damaged or nicked DNA duplex", (role("damaged duplex with ragged ends"),))
    excised = molecule("lesion-excised duplex", (role("short gaps"), role("trimmed overhangs")))
    blunt = molecule("blunt repaired duplex", (role("restricted fill-in"), role("remaining nicks")))
    sealed = molecule("nick-sealed duplex", (role("original duplex interior"),))
    ata = molecule("dA-tailed duplex", (role("original duplex interior"), role("3′ dA", length_bp=1)))
    adapted_construct = Construct((
        role("custom duplex UMI adapter A", tag="umi",
             molecular_feature=_umi("a")),
        role("repaired insert"),
        role("custom duplex UMI adapter B", tag="umi",
             molecular_feature=_umi("b"))),
        name="duplex-UMI library")
    adapted = ligation_product_scene(
        adapted_construct, "repaired insert", adapter_group="duplex UMI adapter A",
        right_adapter_group="duplex UMI adapter B")
    return _workflow(damaged, (
        ("Excise damaged bases and remove single-stranded overhangs", excised,
         "EndoIV, Fpg, UDG, T4 PDG and EndoVIII act with ExoVII."),
        ("Blunt ends and restrict fill-in", blunt,
         "T4 PNK and T4 DNA polymerase fill short gaps and residual short 5′ overhangs."),
        ("Seal nicks with HiFi Taq DNA ligase", sealed,
         "Nick sealing precedes dA-tailing so polymerase does not start from open nicks."),
        ("Restricted dA-tail with dATP as the only nucleotide", ata,
         "Klenow fragment exo− and Taq add the ligation-compatible 3′ dA."),
        ("Ligate custom dual-index duplex-UMI adapters and PCR-amplify", adapted,
         "Complementary UMI families recover both original strands."),
    ))


FINAL_LIBRARY = Construct([
    role("custom left duplex-UMI adapter", tag="umi", molecular_feature=_umi("a")),
    role("repaired insert"),
    role("custom right duplex-UMI adapter", tag="umi", molecular_feature=_umi("b")),
], name="Duplex-Repair sequencing library")
SEQ_PRIMERS = ()
TITLE = 'Duplex-Repair — damage-aware end repair for duplex sequencing'
NOTES = '01_duplex_repair.html'
SOURCE = '<a href="https://doi.org/10.1093/nar/gkab855">Xiong et al., <i>Nucleic Acids Research</i> (2022)</a>.'
SUMMARY = 'Lesion excision, controlled blunting, nick sealing and dATP-only tailing are separated to minimize copying of one-strand damage into an apparent duplex mutation.'
CAVEAT = 'The paper identifies custom dual-index duplex-UMI adapters but does not publish their sequences; adapter geometry and primer landing positions therefore remain unspecified.'
FINAL_CAPTION = 'Published roles of the custom duplex-UMI library; exact adapter bases are unavailable.'
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'After two rounds of hybrid selection and amplification, the authors sequenced custom duplex-UMI libraries on a HiSeq 2500 rapid run (100 bp paired-end) or HiSeq X (151 bp paired-end).'
SEQUENCING_UNAVAILABLE = 'The custom dual-index duplex-UMI adapter and primer oligos are not given in the accessible paper or supplement.'
def oligos():
    return ()

def _cleanup(graph, ident, label, material, kept_label, reviewed, *, state=None):
    """A source-named AMPure cleanup; the article does not itemize removed species."""
    beads = graph.input(ident + "_beads", "Ampure beads", MaterialKind.SELECTION)
    kept = graph.product(ident + "_kept", kept_label, MaterialKind.MOLECULAR,
                         state=material.state if state is None else state)
    waste = graph.product(ident + "_waste", "Bead supernatant and washes", MaterialKind.WASTE)
    graph.react(ident, label, StepKind.CLEANUP,
                ((material, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (kept, waste), changes_molecule=False)
    graph.discard(waste, "Removed from the bead-bound DNA; the source does not list the removed species.")
    reviewed.extend((beads.id, ident, kept.id, waste.id))
    return kept


def protocol_graph():
    final = FINAL_LIBRARY
    source = "Xiong et al., Nucleic Acids Research (2022); doi:10.1093/nar/gkab855"
    graph = ProtocolGraph("Duplex-Repair", authority=source)
    reviewed = []
    damaged = graph.input(
        "damaged_dna", "10–20 ng cfDNA or sheared FFPE/buffy-coat gDNA (damaged or nicked duplexes)",
        MaterialKind.BIOLOGICAL,
        state=molecule("damaged or nicked DNA duplex", (
            role("damaged duplex with ragged ends"),)), primary=True)
    lesion_enzymes = graph.input(
        "lesion_enzymes", "EndoIV/Fpg/UDG/T4-PDG/EndoVIII/ExoVII cocktail in NEBuffer 2 + BSA, 37°C 30 min",
        MaterialKind.CATALYST)
    excised = graph.product(
        "excised_dna", "lesion-excised duplex with trimmed overhangs",
        MaterialKind.MOLECULAR,
        state=molecule("lesion-excised duplex", (role("short gaps"), role("trimmed overhangs"))))
    graph.react(
        "lesion_excision", "Step 1: excise damaged bases and remove single-stranded overhangs",
        StepKind.CHEMICAL_CONVERSION,
        ((damaged, InputRole.TEMPLATE), (lesion_enzymes, InputRole.CATALYST)),
        (excised,))
    blunting = graph.input("blunting_enzymes", "T4 PNK and T4 DNA polymerase (added to step 1 mix)",
                           MaterialKind.CATALYST)
    fill_mix = graph.input("atp_dntps", "0.8 mM ATP and 0.5 mM each dNTP, 37°C 30 min",
                           MaterialKind.CONDITION)
    blunt = graph.product(
        "blunt_dna", "blunt repaired duplex with restricted fill-in",
        MaterialKind.MOLECULAR,
        state=molecule("blunt repaired duplex", (
            role("restricted fill-in"), role("remaining nicks"))))
    graph.react("blunt_ends", "Step 2: blunt ends and fill remaining short gaps/overhangs", StepKind.END_REPAIR,
                ((excised, InputRole.TEMPLATE), (blunting, InputRole.CATALYST),
                 (fill_mix, InputRole.CONDITION)),
                (blunt,))
    nick_ligase = graph.input("nick_ligase", "HiFi Taq DNA ligase + buffer (added to step 2 mix), 35→65°C over 45 min",
                              MaterialKind.CATALYST)
    sealed = graph.product(
        "sealed_dna", "nick-sealed repaired duplex", MaterialKind.MOLECULAR,
        state=molecule("nick-sealed duplex", (role("original duplex interior"),)))
    graph.react("seal_nicks", "Step 3: seal nicks before tailing", StepKind.NICK_SEALING,
                ((blunt, InputRole.TEMPLATE), (nick_ligase, InputRole.CATALYST)),
                (sealed,))
    sealed_clean = _cleanup(graph, "sealed_cleanup", "3X Ampure cleanup; elute in 17 µL 10 mM Tris",
                            sealed, "purified nick-sealed duplex", reviewed)
    tailing = graph.input(
        "restricted_tailing_enzyme", "Klenow fragment exo-minus and Taq polymerase in NEBuffer 2",
        MaterialKind.CATALYST)
    datp = graph.input("datp_only", "0.2 mM dATP as the only nucleotide; RT 30 min then 65°C 30 min",
                       MaterialKind.CONDITION)
    tailed = graph.product(
        "tailed_dna", "restricted dA-tailed repaired duplex",
        MaterialKind.MOLECULAR,
        state=molecule("dA-tailed duplex", (
            role("original duplex interior"), role("3-prime dA", length_bp=1))))
    graph.react("restricted_a_tail", "Step 4: restricted dA-tail with dATP only",
                StepKind.TAILING,
                ((sealed_clean, InputRole.TEMPLATE), (datp, InputRole.CONDITION),
                 (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input("duplex_umi_adapters", "custom dual-index duplex-UMI adapters (IDT)",
                           MaterialKind.ADAPTER)
    adapter_ligase = graph.input("adapter_ligase", "T4 DNA ligase (1000 U) and 5′-deadenylase",
                                 MaterialKind.CATALYST)
    peg = graph.input("ligation_peg", "10% (w/v) PEG 8000 in the step 4 mix (55 µL), RT 1 h",
                      MaterialKind.CONDITION)
    adapted = graph.product("adapted_dna", "duplex-UMI adapter-ligated DNA",
                            MaterialKind.MOLECULAR, state=final)
    graph.react("adapter_ligation", "Ligate custom duplex-UMI adapters in the step 4 mix",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (adapter_ligase, InputRole.CATALYST), (peg, InputRole.CONDITION)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The paper identifies custom duplex-UMI adapters but does not publish their ligating-end sequences."))
    adapted_clean = _cleanup(graph, "ligation_cleanup", "1.2X Ampure cleanup after ligation",
                             adapted, "purified adapter-ligated DNA", reviewed)
    forward = graph.input("pcr_forward", "custom duplex-library PCR forward primer (sequence not published)",
                          MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "custom duplex-library PCR reverse primer (sequence not published)",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase (not specified in this article)",
                             MaterialKind.CATALYST)
    amplified = graph.product("amplified_library", "amplified pre-capture duplex-UMI library",
                              MaterialKind.MOLECULAR, state=final)
    graph.react("library_pcr", "PCR-amplify the duplex-UMI library", StepKind.PCR,
                ((adapted_clean, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The accessible source does not disclose the custom duplex-library primer sequences."))
    panel = graph.input("capture_panel", "IDT pan-cancer hybrid-selection panel; xGen hybridization and wash kit with xGen Universal blockers",
                        MaterialKind.SELECTION)
    captured = graph.product("captured_library", "library after two rounds of hybrid selection",
                             MaterialKind.MOLECULAR, state=final)
    off_target = graph.product("off_target", "uncaptured library molecules", MaterialKind.WASTE)
    graph.react("hybrid_selection", "Two rounds of hybrid selection (between-round handling not stated in this article)",
                StepKind.AFFINITY_CAPTURE,
                ((amplified, InputRole.TEMPLATE), (panel, InputRole.SELECTION)), (captured, off_target))
    graph.discard(off_target, "Molecules not retained by the hybrid-selection capture are washed away.")
    post_forward = graph.input("post_capture_forward", "post-capture PCR forward primer (not specified in this article)",
                               MaterialKind.OLIGO)
    post_reverse = graph.input("post_capture_reverse", "post-capture PCR reverse primer (not specified in this article)",
                               MaterialKind.OLIGO)
    post_polymerase = graph.input("post_capture_pcr", "post-capture PCR polymerase (not specified in this article)",
                                  MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "Duplex-Repair sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("post_capture_pcr_step", "Amplify after the second hybrid selection", StepKind.PCR,
                ((captured, InputRole.TEMPLATE), (post_forward, InputRole.STRUCTURAL),
                 (post_reverse, InputRole.STRUCTURAL), (post_polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The article states that libraries were amplified after the second hybrid selection but gives no primers."))
    graph.sequence(library, SequencingPlan(
        limitation="HiSeq 2500 rapid run (100 bp paired-end) or HiSeq X (151 bp paired-end) is reported, but custom adapter and sequencing-primer bases are unavailable, so read starts and UMI cycles are not derived."))
    graph.declare_reviewed_scope((
        "damaged_dna", "lesion_enzymes", "lesion_excision", "blunting_enzymes", "atp_dntps",
        "blunt_ends", "nick_ligase", "seal_nicks", "restricted_tailing_enzyme",
        "datp_only", "restricted_a_tail", "duplex_umi_adapters", "adapter_ligase", "ligation_peg",
        "adapter_ligation", "pcr_forward", "pcr_reverse", "pcr_polymerase", "amplified_library",
        "library_pcr", "capture_panel", "captured_library", "off_target", "hybrid_selection",
        "post_capture_forward", "post_capture_reverse", "post_capture_pcr", "post_capture_pcr_step", "sequencing_library", *reviewed,
    ), source)
    return graph.require_publishable()
