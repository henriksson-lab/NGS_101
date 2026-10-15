"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_rt_pcr
from rna_library_operations import molecule

DEFINITION = Definition("MARIO — mapping RNA interactome in vivo", "10.1038/ncomms12023", "Nguyen et al., Nature Communications (2016)", "UV-crosslinked RNPs receive a selectable biotinylated RNA linker that joins two RNAs co-bound by one protein.", "truseq", ("UMI", "RNA A cDNA", "MARIO linker", "RNA B cDNA"), "Standard paired-end Illumina reads traverse the two RNA-derived arms and linker. The paper publishes the full PE1/PE2 primers, three inline sample-barcode/UMI RT primers, and every custom linker/adapter oligo.", "RNA", "The original Methods section has a complete source-exact oligo list. Three RT-primer variants differ in their four-base sample barcode and random-barcode length; all three are shown.")


from protocol_specs import workflow
from interaction_capture import ProteinMediatedRNAContact
from mario_chem import oligos as mario_oligos


def _source_workflow():
    x = ProteinMediatedRNAContact("biotinylated RNA linker", True)
    return workflow(molecule("cellular RNPs", ("protein-bound RNAs",)), (
        ("UV-cross-link RNA to its bound proteins", x.complex(),
         "Protein, rather than RNA base pairing, preserves the contact."),
        ("Attach the biotinylated linker and proximity-ligate co-bound RNAs", x.chimera(),
         "The linker creates an explicit RNA1–linker–RNA2 boundary."),
        ("Select linker-containing chimeras and convert them to cDNA", x.chimera(),
         "Streptavidin enrichment precedes paired-end library preparation."),
    ))

SPEC = get('mario', definition=DEFINITION,
           workflow_fn=_source_workflow,
           oligos_fn=mario_oligos)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rnps = graph.input("rnps", "cellular protein-bound RNA assemblies",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("cellular RNPs", ("protein-bound RNAs",)), primary=True)
    uv = graph.input("uv", "UV crosslinking condition", MaterialKind.CONDITION)
    fixed = graph.product("fixed", "UV-crosslinked RNP assemblies", MaterialKind.MOLECULAR,
                          state=rnps.state)
    graph.react("crosslink", "UV-crosslink RNA to bound proteins", StepKind.CROSS_LINKING,
                ((rnps, InputRole.TEMPLATE), (uv, InputRole.CONDITION)), (fixed,))
    linker = graph.input("linker", "biotinylated RNA linker", MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "RNA proximity ligase", MaterialKind.CATALYST)
    chimera = graph.product("chimera", "RNA1-linker-RNA2 chimera", MaterialKind.MOLECULAR,
                            state=molecule("MARIO chimera", ("RNA A", "biotinylated linker", "RNA B"),
                                           junctions=(("RNA A", "biotinylated linker", "linker ligation"),
                                                      ("biotinylated linker", "RNA B", "linker ligation"))))
    graph.react("ligate_linker", "Attach the linker and ligate co-bound RNAs",
                StepKind.LIGATION,
                ((fixed, InputRole.TEMPLATE), (linker, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (chimera,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The public schematic establishes the selectable linker boundary but not one universal end geometry."))
    beads = graph.input("streptavidin", "streptavidin capture beads", MaterialKind.SELECTION)
    captured = graph.product("captured", "linker-containing RNA chimeras", MaterialKind.MOLECULAR,
                             state=chimera.state)
    rejected = graph.product("rejected", "linker-negative RNA", MaterialKind.WASTE)
    graph.react("capture", "Select linker-containing chimeras", StepKind.AFFINITY_CAPTURE,
                ((chimera, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not contain the biotinylated interaction linker")
    reviewed = ["rnps", "uv", "crosslink", "linker", "ligase", "ligate_linker",
                "streptavidin", "capture"]
    finish_rt_pcr(graph, "library", captured, final, primers, d.sequencing, reviewed,
                   label="MARIO chimera")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
