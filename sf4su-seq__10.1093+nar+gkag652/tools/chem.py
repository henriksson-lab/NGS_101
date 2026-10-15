"""sf4sU-Seq source definition and size-selection route."""
from nascent_protocols import Definition, FourSUReadout, make_spec, role_library
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

# Supplementary Methods, "Kits and Reagents", identifies E7300S. The official
# NEB E7300 manual v9.0_07/25 prints these molecules on page 16. The study does
# not identify which multiplex-index tube was used, so no sample index is fixed.
KIT_SOURCE = "NEBNext Small RNA Library Prep Kit E7300S; NEB manual v9.0_07/25, p. 16"
ADAPTER_3P = "AGATCGGAAGAGCACACGTCT"      # 5′-rApp…-NH2-3′
ADAPTER_5P = "GUUCAGAGUUCUACAGUCCGACGAUC"  # RNA, 5′→3′
RT_PRIMER = "AGACGTGTGCTCTTCCGATCT"
PCR_FORWARD = "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"
PCR_REVERSE_5P = "CAAGCAGAAGACGGCATACGAGAT"
PCR_REVERSE_3P = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"


def kit_oligos():
    return (
        "NEBNext E7300 3′ SR adapter — 5′-rApp" + ADAPTER_3P + "-NH2-3′",
        "NEBNext E7300 5′ SR adapter — 5′-" + "".join("r" + base for base in ADAPTER_5P) + "-3′",
        "NEBNext E7300 SR RT primer — 5′-" + RT_PRIMER + "-3′",
        "NEBNext E7300 SR PCR primer — 5′-" + PCR_FORWARD[:-1] + "*" + PCR_FORWARD[-1] + "-3′",
        "NEBNext E7300 index PCR primer family — 5′-" + PCR_REVERSE_5P +
        "[six-base index]" + PCR_REVERSE_3P[:-1] + "*" + PCR_REVERSE_3P[-1] +
        "-3′; the study does not identify its index tube",
    )

DEFINITION = Definition(
    "sf4sU-Seq — size-fractionated 4sU sequencing", "10.1093/nar/gkag652",
    "Gerlevik et al., Nucleic Acids Research (2026)",
    "Biotin-enriched 4sU RNA is gel-selected around 100 nt and sequenced directionally from both ends to map initiation and promoter-proximal pausing.",
    "The two read ends report the native 5′ and 3′ coordinates of the selected short labelled RNA.",
    "01_sf4su_seq.html",
)
READOUT = FourSUReadout(
    "thiol biotinylation, streptavidin purification and ≈100-nt denaturing-gel selection",
    "native RNA ends", True,
)
FINAL_LIBRARY = role_library(DEFINITION.title + " final library",
    (role("assay-derived cDNA"),), indexed=True, platform="Illumina-compatible")
SPEC = make_spec(DEFINITION, FINAL_LIBRARY,
    unavailable="The authoritative paper establishes the library roles and read direction, but does not expose one complete immutable sequencing-primer geometry; adapters are therefore shown by role.")
SPEC.oligos = kit_oligos
SPEC.CAVEAT = ("Supplementary Methods identify NEBNext Small RNA kit E7300S. "
               "The official E7300 manual supplies adapter and primer families; "
               "the study does not identify its index tube or verify the kit lot. "
               "The final amplicon remains role-level pending base-by-base reconstruction.")


def _rna_ligation_geometry(adapter: str, substrate: Construct, *,
                           adapter_end: str, substrate_end: str,
                           statement: str) -> AdapterLigationGeometry:
    scene = Scene()
    scene.strand("substrate", list(substrate), label="selected 4sU RNA")
    scene.strand("adapter", [Segment("kit adapter", adapter)],
                 label="NEBNext E7300 adapter",
                 mod5="rApp" if adapter_end == "5-prime" else "",
                 mod3="NH2" if adapter_end == "5-prime" else "")
    scene.ligation_end("substrate", substrate_end, kind="single-strand",
                       polarity="3-prime" if adapter_end == "5-prime" else "5-prime",
                       owner="substrate")
    scene.ligation_end("adapter", "kit adapter", kind="single-strand",
                       polarity=adapter_end, owner="adapter")
    return AdapterLigationGeometry.known("single-strand", scene, statement=statement)


def protocol_graph():
    source = f"{DEFINITION.citation}; doi:{DEFINITION.doi}"
    graph = ProtocolGraph(DEFINITION.title, authority=source)
    cells = graph.input(
        "cells", "living cells during the 4sU pulse", MaterialKind.BIOLOGICAL,
        state=molecule("living cells", (role("RNA synthesis during pulse"),)),
        primary=True)
    four_su = graph.input("four_su", "4-thiouridine pulse-label nucleotide",
                          MaterialKind.MOLECULAR)
    polymerase = graph.input("cellular_rna_polymerase", "cellular RNA polymerases",
                             MaterialKind.CATALYST)
    labelled = graph.product("labelled_rna", "4sU pulse-labelled RNA",
                             MaterialKind.MOLECULAR, state=READOUT.labelled())
    graph.react("pulse_label", "Pulse-label newly synthesized RNA with 4sU",
                StepKind.TRANSCRIPTION,
                ((cells, InputRole.TEMPLATE), (four_su, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (labelled,))
    reviewed = ["cells", "four_su", "cellular_rna_polymerase", "pulse_label"]

    biotin = graph.input("thiol_biotin", "4sU thiol-biotinylation reagent",
                         MaterialKind.MOLECULAR)
    condition = graph.input("biotinylation_condition", "thiol-biotinylation condition",
                            MaterialKind.CONDITION)
    biotinylated = graph.product("biotinylated_rna", "biotinylated 4sU RNA",
                                MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("biotinylate", "Thiol-biotinylate incorporated 4sU",
                StepKind.CHEMICAL_CONVERSION,
                ((labelled, InputRole.TEMPLATE), (biotin, InputRole.STRUCTURAL),
                 (condition, InputRole.CONDITION)), (biotinylated,))
    streptavidin = graph.input("streptavidin", "streptavidin RNA-capture reagent",
                              MaterialKind.SELECTION)
    enriched = graph.product("enriched_rna", "streptavidin-enriched 4sU RNA",
                             MaterialKind.MOLECULAR, state=READOUT.product())
    unlabelled = graph.product("unlabelled_rna", "unlabelled RNA background",
                               MaterialKind.WASTE)
    graph.react("enrich", "Affinity-purify biotinylated 4sU RNA",
                StepKind.AFFINITY_CAPTURE,
                ((biotinylated, InputRole.TEMPLATE),
                 (streptavidin, InputRole.SELECTION)), (enriched, unlabelled))
    graph.discard(unlabelled, "does not carry the pulse-installed 4sU affinity label")
    reviewed.extend(("thiol_biotin", "biotinylation_condition", "biotinylate",
                     "streptavidin", "enrich"))

    size_rule = graph.input("size_selection", "approximately 100-nt gel-selection rule",
                            MaterialKind.SELECTION)
    selected_rna = Construct((Segment("selected short RNA", "[size-selected 4sU RNA]",
                                      placeholder=True),), name="size-selected 4sU RNA")
    selected = graph.product("selected_rna", "approximately 100-nt enriched 4sU RNA",
                             MaterialKind.MOLECULAR, state=enriched.state,
                             chemical_state=selected_rna)
    rejected = graph.product("rejected_sizes", "enriched RNA outside the selected size",
                             MaterialKind.WASTE)
    graph.react("select_size", "Gel-select the approximately 100-nt RNA fraction",
                StepKind.SELECTION,
                ((enriched, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "falls outside the size-resolved endpoint population")
    reviewed.extend(("size_selection", "select_size"))

    rpph = graph.input("rpph", "RNA 5-prime pyrophosphohydrolase", MaterialKind.CATALYST)
    decapped = graph.product("decapped_rna", "5-prime-monophosphate short RNA",
                            MaterialKind.MOLECULAR, chemical_state=selected_rna)
    graph.react("decap", "Remove the RNA 5-prime cap with RppH", StepKind.DEPHOSPHORYLATION,
                ((selected, InputRole.TEMPLATE), (rpph, InputRole.CATALYST)), (decapped,))
    reviewed.extend(("rpph", "decap"))

    adapter3 = graph.input("adapter_3p", "NEBNext E7300 preadenylated, 3-prime-blocked SR adapter",
                           MaterialKind.ADAPTER,
                           state=Construct((Segment("3-prime SR adapter", ADAPTER_3P),)),
                           source=KIT_SOURCE)
    ligase3 = graph.input("ligase_3p", "NEBNext 3-prime RNA-adapter ligase",
                          MaterialKind.CATALYST)
    three_prime = Construct((*selected_rna, Segment("3-prime SR adapter", ADAPTER_3P)),
                            name="3-prime-adapter-ligated short RNA")
    three_ligated = graph.product("three_prime_ligated", "3-prime-adapter-ligated RNA",
                                  MaterialKind.MOLECULAR, chemical_state=three_prime)
    geometry3 = _rna_ligation_geometry(
        ADAPTER_3P, selected_rna, adapter_end="5-prime", substrate_end="selected short RNA",
        statement=("NEB E7300 manual p. 16 prints the 5′-rApp and 3′-NH2 adapter; "
                   "its activated 5′ end joins the selected RNA 3′ end."))
    graph.react("ligate_3p", "Ligate the preadenylated 3-prime SR adapter",
                StepKind.LIGATION,
                ((decapped, InputRole.TEMPLATE), (adapter3, InputRole.STRUCTURAL),
                 (ligase3, InputRole.CATALYST)), (three_ligated,), geometry=geometry3)
    reviewed.extend(("adapter_3p", "ligase_3p", "ligate_3p"))

    adapter5 = graph.input("adapter_5p", "NEBNext E7300 5-prime RNA SR adapter",
                           MaterialKind.ADAPTER,
                           state=Construct((Segment("5-prime SR adapter", ADAPTER_5P),)),
                           source=KIT_SOURCE)
    ligase5 = graph.input("ligase_5p", "NEBNext 5-prime RNA-adapter ligase",
                          MaterialKind.CATALYST)
    doubly_adapted = Construct((Segment("5-prime SR adapter", ADAPTER_5P),
                               *selected_rna, Segment("3-prime SR adapter", ADAPTER_3P)),
                              name="doubly-adapter-ligated short RNA")
    ligated = graph.product("library_adapter_product", "doubly-adapter-ligated RNA",
                            MaterialKind.MOLECULAR, chemical_state=doubly_adapted)
    geometry5 = _rna_ligation_geometry(
        ADAPTER_5P, three_prime, adapter_end="3-prime", substrate_end="selected short RNA",
        statement=("NEB E7300 manual p. 16 prints the 5′ RNA adapter; "
                   "its 3′ end joins the monophosphorylated RNA 5′ end."))
    graph.react("ligate_5p", "Ligate the 5-prime SR adapter to decapped RNA",
                StepKind.LIGATION,
                ((three_ligated, InputRole.TEMPLATE), (adapter5, InputRole.STRUCTURAL),
                 (ligase5, InputRole.CATALYST)), (ligated,), geometry=geometry5)
    reviewed.extend(("adapter_5p", "ligase_5p", "ligate_5p"))

    rt_primer = graph.input("library_rt_primer", "NEBNext E7300 SR RT primer",
                            MaterialKind.OLIGO,
                            state=Construct((Segment("SR RT primer", RT_PRIMER),)),
                            source=KIT_SOURCE)
    rt = graph.input("library_reverse_transcriptase", "reverse transcriptase",
                     MaterialKind.CATALYST)
    cdna = graph.product("library_cdna", "adapter-primed first-strand cDNA",
                         MaterialKind.MOLECULAR)
    graph.react("library_rt", "Reverse-transcribe the doubly-adapted RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    reviewed.extend(("library_rt_primer", "library_reverse_transcriptase", "library_rt"))

    forward = graph.input("library_pcr_forward", "NEBNext E7300 SR PCR primer",
                          MaterialKind.OLIGO,
                          state=Construct((Segment("SR PCR primer", PCR_FORWARD),)),
                          source=KIT_SOURCE)
    reverse = graph.input("library_pcr_reverse", "NEBNext E7300 indexed PCR primer family",
                          MaterialKind.OLIGO,
                          state=Construct((Segment("P7 and index landing", PCR_REVERSE_5P),
                                           Segment("unidentified six-base index", "N" * 6,
                                                   placeholder=True),
                                           Segment("SR reverse annealing arm", PCR_REVERSE_3P))),
                          source=KIT_SOURCE, inferred=True)
    pcr_polymerase = graph.input("library_pcr_polymerase", "library PCR polymerase",
                                 MaterialKind.CATALYST)
    library = graph.product("library_library", "indexed sequencing library",
                            MaterialKind.LIBRARY, state=SPEC.FINAL_LIBRARY)
    graph.react("library_pcr", "Index and amplify the cDNA library", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (pcr_polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    SPEC.FINAL_LIBRARY, reason="The E7300 primer family is published, but the study's exact index tube and base-level pre-PCR cDNA are not identified."))
    graph.sequence(library, SequencingPlan(limitation=DEFINITION.endpoint))
    reviewed.extend(("library_pcr_forward", "library_pcr_reverse", "library_pcr_polymerase",
                     "library_pcr", "library_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


SPEC.protocol_graph = protocol_graph
