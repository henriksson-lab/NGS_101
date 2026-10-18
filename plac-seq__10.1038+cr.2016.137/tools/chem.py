"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from types import SimpleNamespace
from chemdraw import Construct, Scene
from pcr import PCRGeometry
from protocol_model import molecule, role
from protocol_graph import SequencingPlan
from protocol_graph import Disposition
from protocol_graph import InputRole
from protocol_graph import Terminal
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind

DEFINITION = SimpleNamespace(
    title='PLAC-seq — proximity ligation-assisted ChIP-seq',
    doi='10.1038/cr.2016.137',
    source_label='Fang et al., Cell Research (2016)',
    summary='In-nucleus restriction digestion and proximity ligation precede shearing and chromatin immunoprecipitation, enriching contacts anchored at a chosen protein or mark.',
    sequencing='Paired-end Illumina reads identify both genomic sides of a protein-anchored proximity-ligation product. Data S1 maps each end separately with BWA-MEM, keeps uniquely mapped cis pairs within 500 bp of an MboI site and removes PCR duplicates; it does not state read lengths.',
    caveat='The endpoint records library-adapter roles without assigning a TruSeq shell or index layout. Data S1 performs the biotin pull-down and library construction according to the in situ Hi-C protocol (Rao et al. 2014; not inspected here) with 20 µL T1 beads and a qPCR-chosen minimal PCR cycle number, and prints no adapter or index oligo. Adapter identities, index layout and read lengths remain unresolved within this source scope; this is not an accepted publication-absence claim.',
)


# Data S1 delegates post-ChIP library chemistry to Rao 2014. No generic platform
# selector supplies adapter/index sequences or end geometry at that boundary.
FINAL_LIBRARY = Construct([
    role("left library adapter; exact sequence and index assignment unresolved"),
    role("locus A"), role("biotinylated contact junction"), role("locus B"),
    role("right library adapter; exact sequence and index assignment unresolved"),
], name="PLAC-seq contact library (role-level endpoint)")
SPEC = SimpleNamespace(
    TITLE=DEFINITION.title, NOTES="01_plac_seq.html",
    SOURCE='<a href="https://doi.org/10.1038/cr.2016.137">Fang et al. 2016</a>; published Data S1.',
    SUMMARY=DEFINITION.summary, CAVEAT=DEFINITION.caveat,
    FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=(),
    FINAL_CAPTION=DEFINITION.sequencing, SEQUENCING_INTRO=DEFINITION.sequencing,
    SEQUENCING_ENDING=DEFINITION.sequencing,
    SEQUENCING_UNAVAILABLE="Paired-end sequencing is established by Data S1; primer identities, index assignment and acquired cycle counts remain unassigned pending the delegated library method and run metadata.",
    READ_LENGTHS=None, oligos=lambda: (),
)


def filled_ends_scene():
    """Separate restriction-fragment ends before proximity ligation."""
    scene = Scene()
    for name in ("A", "B"):
        strand = "locus_" + name
        scene.strand(strand, [role("filled locus " + name + " end")],
                     label="separate filled locus " + name + " end")
        scene.labels(strand)
    return scene


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    data_s1 = f"{source}; Supplementary information Data S1 (Methods & Materials)"
    graph = ProtocolGraph(d.title, authority=data_s1)
    nuclei = graph.input("nuclei", "0.5–5 M 1% formaldehyde-fixed cells, lysed and SDS-permeabilized (62 °C), Triton-quenched",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("crosslinked nuclei", (role("chromatin contacts"),)), primary=True)
    restriction = graph.input("restriction_enzyme", "MboI, 100 U in NEBuffer 2 (37 °C, 2 h; inactivated 62 °C, 20 min)",
                              MaterialKind.CATALYST)
    digested = graph.product("digested", "restriction-digested chromatin contacts",
                             MaterialKind.MOLECULAR, state=nuclei.state)
    graph.react("digest", "MboI-digest chromatin in nuclei, then heat-inactivate MboI", StepKind.DIGESTION,
                ((nuclei, InputRole.TEMPLATE), (restriction, InputRole.CATALYST)), (digested,))
    fill_primer = graph.input("biotin_nucleotides", "15 nmol each dCTP, dGTP, dTTP and biotin-14-dATP",
                              MaterialKind.OLIGO)
    polymerase = graph.input("fill_polymerase", "Klenow fragment, 40 U (37 °C, 1.5 h)",
                             MaterialKind.CATALYST)
    filled = graph.product("filled", "biotin-filled restriction fragment ends",
                           MaterialKind.MOLECULAR,
                           state=filled_ends_scene())
    graph.react("fill", "Fill restriction ends with biotinylated nucleotides",
                StepKind.EXTENSION,
                ((digested, InputRole.TEMPLATE), (fill_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (filled,))
    ligase = graph.input("proximity_ligase", "T4 DNA ligase, 4000 U in 1.2 mL with 1% Triton X-100 and BSA (RT)", MaterialKind.CATALYST)
    contacts = graph.product("contacts", "biotinylated proximity-ligation products",
                             MaterialKind.MOLECULAR,
                             state=molecule("contact product", (role("locus A"), role("locus B")),
                                            junctions=(("locus A", "locus B",
                                                        "biotinylated proximity ligation"),)))
    graph.react("proximity_ligation", "Proximity-ligate inside intact nuclei",
                StepKind.PROXIMITY_LIGATION,
                ((filled, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (contacts,))
    shear = graph.input("shear_condition", "RIPA lysis on ice, Covaris M220 sonication (75 W, 10% duty, 200 cycles/burst, 10 min, 7 °C)",
                        MaterialKind.CONDITION)
    sheared = graph.product("sheared", "soluble sheared contact chromatin", MaterialKind.MOLECULAR,
                            state=contacts.state)
    insoluble = graph.product("insoluble_chromatin", "pellet after 14,000 rpm clearing", MaterialKind.WASTE)
    graph.react("shear", "Pellet nuclei, lyse in RIPA, sonicate and clear the lysate", StepKind.FRAGMENTATION,
                ((contacts, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (sheared, insoluble))
    graph.discard(insoluble, "is pelleted; the cleared supernatant is collected")
    antibody = graph.input("antibody", "H3K27ac (ab4729, 2.5 µg), H3K4me3 (04-745, 2.5 µg) or Pol II (ab817, 5 µg) antibody with BSA-blocked Protein G Sepharose",
                           MaterialKind.SELECTION)
    chip = graph.product("chip", "bead-bound protein-anchored contact chromatin",
                         MaterialKind.MOLECULAR, state=sheared.state)
    chip_waste = graph.product("chip_waste", "pre-clearing bead material, unbound lysate and washes", MaterialKind.WASTE)
    input_control = graph.product("input_control", "~5% pre-cleared lysate saved as input control",
                                  MaterialKind.MOLECULAR, state=sheared.state)
    graph.react("chip_capture", "Pre-clear with Protein G (3 h), save ~5% input, immunoprecipitate overnight and wash (RIPA, high-salt RIPA, LiCl, TE)",
                StepKind.AFFINITY_CAPTURE,
                ((sheared, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (chip, chip_waste, input_control))
    graph.discard(chip_waste, "does not contain the selected protein-anchored contact")
    graph.mark_terminal(Terminal(input_control.id, Disposition.STORED,
                                 "Data S1 saves it as an input control without describing its library processing"))
    reversal = graph.input("reverse_crosslinking", "RNase A (10 µg, 37 °C 1 h) then proteinase K (20 µg, 65 °C overnight or ≥2 h)",
                           MaterialKind.CATALYST)
    purified = graph.product("purified_chip_dna", "purified ChIP contact DNA",
                             MaterialKind.MOLECULAR, state=sheared.state)
    protein_waste = graph.product("protein_waste", "beads, digested protein, RNA and organic phase", MaterialKind.WASTE)
    graph.react("reverse_crosslink", "Reverse crosslinks, phenol/chloroform/isoamyl alcohol extract and ethanol-precipitate",
                StepKind.DEPROTEINIZATION,
                ((chip, InputRole.TEMPLATE), (reversal, InputRole.CATALYST)),
                (purified, protein_waste))
    graph.discard(protein_waste, "is removed by extraction and precipitation")
    streptavidin = graph.input("streptavidin", "Dynabeads MyOne Streptavidin T1, 20 µL per sample",
                               MaterialKind.SELECTION)
    enriched = graph.product("enriched", "biotin-junction-enriched PLAC-seq contacts",
                             MaterialKind.MOLECULAR, state=chip.state)
    nonjunction = graph.product("nonjunction", "non-biotinylated ChIP DNA", MaterialKind.WASTE)
    graph.react("junction_capture", "Enrich biotinylated contact junctions (in situ Hi-C pull-down, 20 µL T1 beads)",
                StepKind.AFFINITY_CAPTURE,
                ((purified, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (enriched, nonjunction))
    graph.discard(nonjunction, "does not retain a biotinylated proximity-ligation junction")
    adapted = graph.product("adapter_ready", "bead-bound PLAC-seq contacts with Illumina adapter ends",
                            MaterialKind.MOLECULAR, state=Construct([
                                role("left adapter-bearing end; unresolved"), role("locus A"),
                                role("biotinylated contact junction"), role("locus B"),
                                role("right adapter-bearing end; unresolved"),
                            ], name="PLAC-seq pre-PCR substrate (unresolved ends)"))
    graph.react("hic_library_construction", "Construct the library as in the in situ Hi-C protocol",
                StepKind.MAGIC, ((enriched, InputRole.TEMPLATE),), (adapted,),
                reason="Data S1 delegates library construction after biotin pull-down to the cited in situ Hi-C protocol (Rao et al. 2014), which was not inspected in this pass; adapter, end-processing and enzyme details are therefore not represented.")
    reviewed = ["nuclei", "restriction_enzyme", "digest", "biotin_nucleotides",
                "fill_polymerase", "fill", "proximity_ligase", "proximity_ligation",
                "shear_condition", "shear", "antibody", "chip_capture",
                "reverse_crosslinking", "reverse_crosslink", "streptavidin",
                "junction_capture", "hic_library_construction"]
    # Data S1 lines 57–61 establishes amplification and qPCR cycle selection,
    # without assigning the delegated primer/enzyme choices in this review.
    forward = graph.input("library_forward", "delegated in situ Hi-C forward PCR primer (identity unresolved)", MaterialKind.OLIGO)
    reverse = graph.input("library_reverse", "delegated in situ Hi-C reverse PCR primer (identity unresolved)", MaterialKind.OLIGO)
    polymerase = graph.input("library_polymerase", "library PCR polymerase (identity unresolved in inspected Data S1)", MaterialKind.CATALYST)
    library = graph.product("library_library", "PLAC-seq contact sequencing library", MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Amplify PLAC-seq contacts with a qPCR-chosen minimal cycle number", StepKind.PCR,
                ((adapted, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    adapted.state, reason="Data S1 lines 57–61 delegates library construction to Rao 2014 and specifies qPCR cycle selection; exact primer sites on this pre-PCR substrate remain unresolved pending that method."))
    graph.sequence(library, SequencingPlan(library=final, acquired_roles=("Read 1", "Read 2"),
        limitation=SPEC.SEQUENCING_UNAVAILABLE))
    reviewed.extend(("library_forward", "library_reverse", "library_polymerase", "library_pcr", "library_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
