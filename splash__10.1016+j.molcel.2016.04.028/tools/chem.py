"""SPLASH source-specific page specification and psoralen choices."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Scene, Segment
from circular import CircularizationGeometry
from interaction_capture import ReversibleDuplexCrosslink
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_specs import Definition, get
from rna_library_operations import molecule
from splash_oligos import (INDEX1_PCR_PRIMER, RNA_LINKER,
                           UNIVERSAL_PCR_PRIMER, oligos)

DEFINITION = Definition(
    "SPLASH — sequencing psoralen-crosslinked RNA hybrids",
    "10.1016/j.molcel.2016.04.028", "Aw et al., Molecular Cell (2016)",
    "Biotinylated psoralen photocrosslinks RNA duplexes for streptavidin selection, proximity ligation and sequencing.",
    "custom", ("RNA A-derived cDNA", "RNA B-derived cDNA"),
    "The authors' 2017 JoVE protocol specifies paired-end 150-base sequencing. Its official Materials table prints the activated 3-prime adapter and 96-nt RT/circularization primer, and names NEB E7335L for library PCR. Exact circle-to-amplicon and read-primer geometry remain unresolved.",
    "RNA",
    "The official 2017 JoVE Materials table prints the custom adapter and RT primer, including 5rApp, 3ddC and two internal iSp18 modifications, and names the NEB E7335L PCR primer kit. The available NEB manual is revision 8.0 (June 2024), so its printed primer family is a later kit reference, not proof of the 2017 lot's selected index. Full circle-derived read geometry is unresolved. The recovered 2016 Molecular Cell supplement (Figure S2F and Supplemental Experimental Procedures) independently prints the same adapter and two-spacer RT primer, labels the cDNA 5′-phosphate, and shows circularization/PCR. Its 3′ adapter string does not annotate the activation/blocking groups supplied by the later author protocol. The original main article remains inaccessible; base-resolved circle-to-amplicon/read geometry is still a model gap.",
)

# Workflow reagent, workflow selection, graph reagent and graph selection.
PSORALEN = ("biotinylated psoralen",
            "streptavidin selects biotin-psoralen RNA duplexes",
            "biotinylated psoralen", "streptavidin capture")

# The author materials list names NEB E7335L; the available NEB Rev 8.0
# manual prints these Set 1 primer bases with a phosphorothioate bond to
# terminal T. Historical 2017 lot identity remains unverified.
# The reaction's exact circle-derived PCR template remains unavailable.
NEB_SET1_SOURCE = "https://www.neb.com/en/-/media/nebus/files/manuals/manuale7335_e7500_-e7710_e7730.pdf"

SPEC = get("splash", definition=DEFINITION, psoralen=PSORALEN)
SPEC.oligos = oligos
SPEC.SEQUENCING_UNAVAILABLE = (
    "The activated adapter and 96-nt RT/circularization-primer bases are "
    "published. The detailed protocol names NEB E7335L Universal PCR Primer and "
    "Index (X) Primer; the exact circle junction and completed amplicon have "
    "not yet been derived, so read-primer sites remain unresolved.")
SOURCE_OLIGO_EXCEPTIONS = {
    "rna_linker": ("The retrieved JoVE author PDF contains the reaction protocol but "
                   "omits the web article's materials-list sequences; the author's "
                   "HTML materials list prints 5rApp-CTGTAGGCACCATCAAT-3ddC."),
}


def rna_linker_geometry() -> AdapterLigationGeometry:
    """The author protocol's activated linker meets a repaired RNA 3′-OH."""
    sc = Scene()
    sc.strand("RNA chimera", [Segment("RNA 3-prime end", "X" * 18,
                                      placeholder=True)], mod3="OH")
    sc.ligation_end("RNA chimera", "RNA 3-prime end", kind="single-strand",
                    polarity="3-prime", owner="substrate", text="RNA 3′-OH")
    sc.strand("RNA linker", [Segment("activated RNA linker", RNA_LINKER)],
              mod5="rApp", mod3="ddC")
    sc.ligation_end("RNA linker", "activated RNA linker", kind="single-strand",
                    polarity="5-prime", owner="adapter", text="linker 5′-rApp")
    return AdapterLigationGeometry.known(
        "single-strand", sc,
        statement="Aw et al. JoVE 2017 section 5 repairs RNA 3′ ends to OH; section 7 "
                  "ligates the printed 5′-rApp, 3′-ddC RNA linker with T4 Rnl2(tr).")


def protocol_graph():
    source = f"{DEFINITION.source_label}; doi:{DEFINITION.doi}"
    graph = ProtocolGraph(DEFINITION.title, authority=source)
    reagent, selection_label = PSORALEN[2:]
    topology = ReversibleDuplexCrosslink(reagent, selection_label)
    rna = graph.input("rna", "cellular RNA containing paired regions",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("cellular RNA", ("RNA molecules",)), primary=True)
    crosslinker = graph.input("crosslinker", reagent + " and long-wave UV condition",
                              MaterialKind.CONDITION)
    crosslinked = graph.product("crosslinked", reagent + "-crosslinked RNA duplexes",
                                MaterialKind.MOLECULAR, state=topology.duplex())
    graph.react("crosslink", "Photocrosslink paired RNA with " + reagent,
                StepKind.CROSS_LINKING,
                ((rna, InputRole.TEMPLATE), (crosslinker, InputRole.CONDITION)),
                (crosslinked,))
    fragmentation = graph.input("fragmentation_condition", "RNA fragmentation condition",
                                MaterialKind.CONDITION)
    fragments = graph.product("fragments", "fragmented crosslinked RNA duplexes",
                              MaterialKind.MOLECULAR, state=crosslinked.state)
    graph.react("fragment", "Fragment the crosslinked RNA", StepKind.FRAGMENTATION,
                ((crosslinked, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CONDITION)), (fragments,))
    selector = graph.input("selector", selection_label, MaterialKind.SELECTION)
    selected = graph.product("selected", "enriched crosslinked RNA duplexes",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("rejected", "nonselected RNA fragments", MaterialKind.WASTE)
    graph.react("select_duplexes", "Enrich crosslinked RNA duplexes",
                StepKind.SELECTION,
                ((fragments, InputRole.TEMPLATE), (selector, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "does not retain the selected crosslinked duplex population")
    phosphatase = graph.input("end_phosphatase", "T4 PNK 3′-cyclic-phosphate repair",
                              MaterialKind.CATALYST)
    three_oh = graph.product("three_oh", "crosslinked RNA duplexes with 3′-OH ends",
                             MaterialKind.MOLECULAR, state=selected.state)
    graph.react("repair_three_prime", "Convert RNA 3′ cyclic phosphates to 3′-OH",
                StepKind.DEPHOSPHORYLATION,
                ((selected, InputRole.TEMPLATE), (phosphatase, InputRole.CATALYST)),
                (three_oh,))
    kinase = graph.input("end_kinase", "T4 PNK plus ATP for RNA 5′ phosphorylation",
                         MaterialKind.CATALYST)
    five_phosphate = graph.product("five_phosphate", "ligatable RNA duplex ends",
                                   MaterialKind.MOLECULAR, state=three_oh.state)
    graph.react("phosphorylate_five_prime", "Convert RNA 5′-OH ends to 5′ phosphate",
                StepKind.PHOSPHORYLATION,
                ((three_oh, InputRole.TEMPLATE), (kinase, InputRole.CATALYST)),
                (five_phosphate,))
    ligase = graph.input("proximity_ligase", "T4 RNA ligase 1", MaterialKind.CATALYST)
    chimera = graph.product("chimera", "proximity-ligated RNA chimera",
                            MaterialKind.MOLECULAR, state=topology.chimera())
    graph.react("proximity_ligation", "Proximity-ligate the paired RNA fragments",
                StepKind.PROXIMITY_LIGATION,
                ((five_phosphate, InputRole.TEMPLATE),
                 (ligase, InputRole.CATALYST)), (chimera,))
    reversal = graph.input("photoreversal", reagent + " photoreversal condition",
                           MaterialKind.CONDITION)
    reversed_chimera = graph.product("reversed_chimera", "photoreversed RNA chimera",
                                     MaterialKind.MOLECULAR, state=chimera.state)
    graph.react("photoreverse", "Photoreverse the psoralen crosslink",
                StepKind.CHEMICAL_CONVERSION,
                ((chimera, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)),
                (reversed_chimera,))
    linker = graph.input("rna_linker", "published activated 3′ RNA linker",
                         MaterialKind.ADAPTER, state=RNA_LINKER,
                         source="Aw et al. JoVE 2017, materials list")
    linker_ligase = graph.input("linker_ligase", "T4 RNA ligase 2 truncated",
                                MaterialKind.CATALYST)
    linked_state = Scene()
    linked_state.strand("linked RNA", [
        Segment("RNA chimera", "X" * 36, placeholder=True),
        Segment("3-prime RNA linker", RNA_LINKER)], label="3′-linker-ligated RNA")
    linked_state.junction("linked RNA", "RNA chimera", "3-prime RNA linker",
                          "RNA 3′-OH to linker 5′-rApp")
    linked_rna = graph.product("linked_rna", "3′-linker-ligated RNA chimera",
                               MaterialKind.MOLECULAR, state=linked_state)
    graph.react("ligate_rna_linker", "Ligate the activated 3′ RNA linker",
                StepKind.LIGATION,
                ((reversed_chimera, InputRole.TEMPLATE),
                 (linker, InputRole.STRUCTURAL),
                 (linker_ligase, InputRole.CATALYST)),
                (linked_rna,), geometry=rna_linker_geometry())
    rt_primer = graph.input("library_rt_primer", "published 96-nt RT primer",
                            MaterialKind.OLIGO,
                            source="Aw et al. JoVE 2017, materials list")
    rt_enzyme = graph.input("library_rt_enzyme", "reverse transcriptase",
                            MaterialKind.CATALYST)
    cdna = graph.product("library_cdna", "RT-primer-extended SPLASH cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule("SPLASH cDNA", ("cDNA insert", "96-nt RT primer")))
    graph.react("library_rt", "Reverse-transcribe the linker-ligated RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((linked_rna, InputRole.TEMPLATE),
                 (rt_primer, InputRole.STRUCTURAL),
                 (rt_enzyme, InputRole.CATALYST)), (cdna,))
    circle_ligase = graph.input("circle_ligase", "CircLigase II ssDNA ligase",
                                MaterialKind.CATALYST)
    circle = graph.product("circular_cdna", "circularized SPLASH cDNA",
                           MaterialKind.MOLECULAR,
                           source="Aw et al. Molecular Cell 2016 Figure S2F; first-strand cDNA is labelled 5′P")
    graph.react("circularize_cdna", "Circularize single-stranded cDNA",
                StepKind.CIRCULARIZATION,
                ((cdna, InputRole.TEMPLATE),
                 (circle_ligase, InputRole.CATALYST)), (circle,),
                geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="Aw et al. 2016 original Figure S2F labels first-strand cDNA 5′P and shows its circularization; "
                              "the 2017 detailed protocol section 7/Table 6 names CircLigase II. "
                              "A base-resolved circular junction is not assembled in this model."))
    reviewed = ["rna", "crosslinker", "crosslink", "fragmentation_condition", "fragment",
                "selector", "select_duplexes", "end_phosphatase", "repair_three_prime",
                "end_kinase", "phosphorylate_five_prime", "proximity_ligase",
                "proximity_ligation", "photoreversal", "photoreverse", "rna_linker",
                "linker_ligase", "ligate_rna_linker", "library_rt_primer",
                "library_rt_enzyme", "library_rt", "circle_ligase", "circularize_cdna"]
    forward = graph.input("library_forward", "NEB E7335L universal PCR primer",
                          MaterialKind.OLIGO, state=UNIVERSAL_PCR_PRIMER,
                          source=NEB_SET1_SOURCE + " (3′ terminal T on phosphorothioate)")
    reverse = graph.input("library_reverse", "NEB E7335L Index 1 PCR primer",
                          MaterialKind.OLIGO, state=INDEX1_PCR_PRIMER,
                          source=NEB_SET1_SOURCE + " (representative index 1; 3′ terminal T on phosphorothioate)")
    polymerase = graph.input("library_polymerase", "SPLASH library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("library_library", "SPLASH sequencing library",
                            MaterialKind.LIBRARY, state=SPEC.FINAL_LIBRARY)
    graph.react("library_pcr", "Amplify and index the SPLASH cDNA circle",
                StepKind.PCR,
                ((circle, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    SPEC.FINAL_LIBRARY,
                    reason="NEB E7335L primer sequences are known, but the circular "
                           "SPLASH cDNA junction and both inward sites have not been "
                           "derived base-for-base."))
    graph.sequence(library, SequencingPlan(limitation=DEFINITION.sequencing))
    reviewed.extend(("library_forward", "library_reverse", "library_polymerase",
                     "library_pcr", "library_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


SPEC.protocol_graph = protocol_graph
