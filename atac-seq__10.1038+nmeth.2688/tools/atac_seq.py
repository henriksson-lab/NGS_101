"""Buenrostro et al. 2013 bulk ATAC-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import nextera_library, seg
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from protocol_families.atac import ATAC_AD1_NOMX, ATAC_AD2_1
from chemdraw import Construct, Row, Scene
from chromatin_epigenetics import tagmented_library, tagmentation_scene
import nextera as nx
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, Input, InputRole, Material, MaterialKind,
                            ProtocolGraph, SequencingPlan, SourceCoverage, SourceItem,
                            Step, StepKind, Terminal)

TITLE = "ATAC-seq — transposition into accessible chromatin"
NOTES = "01_atac-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.2688">Buenrostro et al., <i>Nature Methods</i> (2013)</a>.'
SUMMARY = "Tn5 simultaneously cleaves accessible native chromatin and transfers sequencing-adapter ends; extension and indexed PCR complete the paired-end library."
CAVEAT = "The paper uses Nextera transposase and prints its PCR primers. The final construct is assembled from the repository’s canonical Nextera components."
FINAL_LIBRARY, SEQ_PRIMERS = nextera_library(
    [seg("accessible chromatin fragment", "X" * 36, placeholder=True)],
    "accessible chromatin fragment library", dual_index=False)
FINAL_CAPTION = "Productive heterologous Tn5 product after 72 °C extension and limited indexed PCR."
SEQUENCING_INTRO = "Paired reads enter the accessible chromatin fragment from its two transposition sites; the index reads identify the PCR-added sample indexes."
SHOW_READ_CONTENT = True

def sections():
    return [
        ("Expose native chromatin", [Row(chunks=[("nucleosome -- linker -- nucleosome     [open regulatory DNA]", None, False)])],
         "Mild lysis releases nuclei while retaining chromatin architecture."),
        ("Tagment accessible DNA", tagmentation_scene("accessible chromatin fragment").rows(),
         "Tn5 cleavage and mosaic-end transfer are coupled; ** marks both transferred adapter–insert boundaries."),
        ("Extend and amplify", [Row(chunks=[("72 °C fill/extension → indexed PCR → size-distributed library", None, False)])],
         "The initial extension makes both transposase-tagged ends amplifiable; PCR cycle number is limited by a qPCR side reaction."),
    ]


def _tagmented_template() -> Construct:
    """Productive heterologous S5/S7 molecule before PCR adds flow-cell arms."""
    return Construct([
        seg("S5", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"),
        seg("accessible chromatin fragment", "X" * 36, placeholder=True),
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("S7 reverse complement", nx.S7_RC, "s7"),
    ], name="productive ATAC-seq tagmentation product")


def protocol_graph() -> ProtocolGraph:
    source = "Buenrostro et al. 2013, doi:10.1038/nmeth.2688"
    graph = ProtocolGraph(TITLE, authority=source)
    graph.add_input(Material("native_chromatin", "native chromatin",
                             MaterialKind.BIOLOGICAL, primary=True, source=source))
    graph.add_input(Material("loaded_mosaic_end", "loaded Nextera mosaic-end adapters",
                             MaterialKind.ADAPTER, source=source))
    graph.add_input(Material("tn5", "Tn5 transposase", MaterialKind.CATALYST,
                             source=source))
    tagged = Material("tagged_fragment", "productive heterologous S5/S7 fragment",
                      MaterialKind.MOLECULAR, state=_tagmented_template(), source=source)
    off_product = Material(
        "non_library_tagmentation_products",
        "S5/S5, S7/S7, untagged and otherwise non-library fragments",
        MaterialKind.WASTE, source=source)
    graph.add_step(Step(
        "tagment", "Tagment accessible chromatin", StepKind.TAGMENTATION,
        (Input("native_chromatin", InputRole.TEMPLATE),
         Input("loaded_mosaic_end", InputRole.STRUCTURAL),
         Input("tn5", InputRole.CATALYST)),
        (tagged.id, off_product.id), source=source), (tagged, off_product))
    graph.mark_terminal(Terminal(off_product.id, Disposition.DISCARD,
                                 "not a complete heterologous sequencing molecule"))

    graph.add_input(Material("extension_polymerase", "72 °C extension polymerase",
                             MaterialKind.CATALYST, source=source))
    extended = Material("extended_fragment", "gap-filled tagmented fragment",
                        MaterialKind.MOLECULAR, state=_tagmented_template(), source=source)
    graph.add_step(Step(
        "extend", "Extend the transposase-tagged ends", StepKind.EXTENSION,
        (Input(tagged.id, InputRole.TEMPLATE),
         Input("extension_polymerase", InputRole.CATALYST)),
        (extended.id,), source=source), (extended,))

    graph.add_input(Material("ad1_nomx", "Ad1_noMX", MaterialKind.OLIGO, source=source))
    graph.add_input(Material("ad2_indexed", "Ad2.1_TAAGGCGA", MaterialKind.OLIGO,
                             source=source))
    geometry = PCRGeometry.exact(
        FINAL_LIBRARY, PCRPrimer("Ad1_noMX", ATAC_AD1_NOMX),
        PCRPrimer("Ad2.1_TAAGGCGA", ATAC_AD2_1),
        statement="Supplementary Table 1 prints both amplification primers; their 3-prime mosaic-end regions face inward across the accessible-DNA insert.")
    library = Material("atac_library", "indexed ATAC-seq library", MaterialKind.LIBRARY,
                       state=FINAL_LIBRARY, source=source)
    graph.add_step(Step(
        "indexed_pcr", "Limited indexed PCR", StepKind.PCR,
        (Input(extended.id, InputRole.TEMPLATE),
         Input("ad1_nomx", InputRole.STRUCTURAL),
         Input("ad2_indexed", InputRole.STRUCTURAL)),
        (library.id,), geometry=geometry, source=source), (library,))
    graph.mark_terminal(Terminal(library.id, Disposition.SEQUENCING), SequencingPlan(
        FINAL_LIBRARY, tuple(SEQ_PRIMERS), tuple(p.role for p in SEQ_PRIMERS)))

    graph.declare_reviewed_scope(dict((
            ("native_nuclei", "native_chromatin"),
            ("loaded_adapter_reagent", "loaded_mosaic_end"),
            ("tn5_reagent", "tn5"),
            ("tagmentation", "tagment"),
            ("productive_and_nonproductive_ends", "non_library_tagmentation_products"),
            ("extension_polymerase_reagent", "extension_polymerase"),
            ("extension", "extend"),
            ("ad1_primer_reagent", "ad1_nomx"),
            ("ad2_primer_reagent", "ad2_indexed"),
            ("indexed_pcr", "indexed_pcr"),
            ("paired_end_sequencing", "atac_library"))), source)
    return graph.require_publishable()
