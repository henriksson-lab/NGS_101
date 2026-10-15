"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_pcr
from rna_library_operations import molecule

DEFINITION = Definition("PLAC-seq — proximity ligation-assisted ChIP-seq", "10.1038/cr.2016.137", "Fang et al., Cell Research (2016)", "In-nucleus restriction digestion and proximity ligation precede shearing and chromatin immunoprecipitation, enriching contacts anchored at a chosen protein or mark.", "truseq", ("locus A", "biotinylated contact junction", "locus B"), "Paired-end Illumina reads identify both genomic sides of a protein-anchored proximity-ligation product.", "DNA")


from protocol_specs import workflow, two_branch


def _source_workflow():
    return workflow(molecule("cross-linked nuclei", ("chromatin contacts",)), (
        ("Restriction-digest chromatin and fill ends with biotinylated nucleotides",
         two_branch("one chromatin contact", "filled locus A end", "filled locus B end"),
         "Biotin marks the ends that can become contact junctions."),
        ("Proximity-ligate inside intact nuclei",
         molecule("contact product", ("locus A", "locus B"), joined=True,
                  junctions=(("locus A", "locus B", "biotinylated proximity ligation"),)),
         "The junction preserves the same in-nucleus contact geometry as in situ Hi-C."),
        ("Shear, immunoprecipitate the chosen chromatin target, and enrich biotin",
         molecule("protein-anchored contact library", ("left Illumina adapter", "locus A", "locus B", "right Illumina adapter"), joined=True,
                  junctions=(("left Illumina adapter", "locus A", "adapter ligation"),
                             ("locus A", "locus B", "contact ligation"),
                             ("locus B", "right Illumina adapter", "adapter ligation"))),
         "ChIP occurs after proximity ligation, defining the distinction from HiChIP."),
    ))

SPEC = get('plac-seq', definition=DEFINITION,
           workflow_fn=_source_workflow)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input("nuclei", "crosslinked nuclei containing chromatin contacts",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("crosslinked nuclei", ("chromatin contacts",)), primary=True)
    restriction = graph.input("restriction_enzyme", "Hi-C restriction endonuclease",
                              MaterialKind.CATALYST)
    digested = graph.product("digested", "restriction-digested chromatin contacts",
                             MaterialKind.MOLECULAR, state=nuclei.state)
    graph.react("digest", "Restriction-digest chromatin in nuclei", StepKind.DIGESTION,
                ((nuclei, InputRole.TEMPLATE), (restriction, InputRole.CATALYST)), (digested,))
    fill_primer = graph.input("biotin_nucleotides", "biotinylated fill-in nucleotides",
                              MaterialKind.OLIGO)
    polymerase = graph.input("fill_polymerase", "restriction-end fill-in polymerase",
                             MaterialKind.CATALYST)
    filled = graph.product("filled", "biotin-filled restriction fragment ends",
                           MaterialKind.MOLECULAR,
                           state=molecule("filled contact ends", ("filled locus A end", "filled locus B end")))
    graph.react("fill", "Fill restriction ends with biotinylated nucleotides",
                StepKind.EXTENSION,
                ((digested, InputRole.TEMPLATE), (fill_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (filled,))
    ligase = graph.input("proximity_ligase", "in-nucleus DNA ligase", MaterialKind.CATALYST)
    contacts = graph.product("contacts", "biotinylated proximity-ligation products",
                             MaterialKind.MOLECULAR,
                             state=molecule("contact product", ("locus A", "locus B"),
                                            junctions=(("locus A", "locus B",
                                                        "biotinylated proximity ligation"),)))
    graph.react("proximity_ligation", "Proximity-ligate inside intact nuclei",
                StepKind.PROXIMITY_LIGATION,
                ((filled, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (contacts,))
    shear = graph.input("shear_condition", "chromatin shearing condition", MaterialKind.CONDITION)
    sheared = graph.product("sheared", "sheared contact products", MaterialKind.MOLECULAR,
                            state=contacts.state)
    graph.react("shear", "Shear proximity-ligated chromatin", StepKind.FRAGMENTATION,
                ((contacts, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (sheared,))
    antibody = graph.input("antibody", "chosen chromatin-target ChIP antibody",
                           MaterialKind.SELECTION)
    chip = graph.product("chip", "protein-anchored contact products",
                         MaterialKind.MOLECULAR, state=sheared.state)
    chip_waste = graph.product("chip_waste", "non-target chromatin contacts", MaterialKind.WASTE)
    graph.react("chip_capture", "Immunoprecipitate the selected chromatin target",
                StepKind.AFFINITY_CAPTURE,
                ((sheared, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (chip, chip_waste))
    graph.discard(chip_waste, "does not contain the selected protein-anchored contact")
    streptavidin = graph.input("streptavidin", "streptavidin junction-capture beads",
                               MaterialKind.SELECTION)
    enriched = graph.product("enriched", "biotin-junction-enriched PLAC-seq contacts",
                             MaterialKind.MOLECULAR, state=chip.state)
    nonjunction = graph.product("nonjunction", "non-biotinylated ChIP material", MaterialKind.WASTE)
    graph.react("junction_capture", "Enrich biotinylated contact junctions",
                StepKind.AFFINITY_CAPTURE,
                ((chip, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (enriched, nonjunction))
    graph.discard(nonjunction, "does not retain a biotinylated proximity-ligation junction")
    reviewed = ["nuclei", "restriction_enzyme", "digest", "biotin_nucleotides",
                "fill_polymerase", "fill", "proximity_ligase", "proximity_ligation",
                "shear_condition", "shear", "antibody", "chip_capture", "streptavidin",
                "junction_capture"]
    finish_pcr(graph, "library", enriched, final, primers, d.sequencing, reviewed,
                label="PLAC-seq contact")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
