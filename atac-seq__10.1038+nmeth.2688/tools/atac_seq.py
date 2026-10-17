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
SEQUENCING_INTRO = ("Paired reads enter the accessible chromatin fragment from its two transposition sites; "
                    "one 8-nt index read identifies the Ad2.N sample index (Ad1_noMX carries no index). "
                    "Online Methods (Primary data processing) report 34 × 8 × 34 reads on MiSeq or "
                    "50 × 8 × 50 on HiSeq; the HiSeq configuration is drawn.")
READ_LENGTHS = {"Read 1": 50, "Index 1 (i7)": 8, "Read 2": 50}
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

    graph.add_input(Material("minelute", "Qiagen MinElute column", MaterialKind.SELECTION,
                             source=source))
    purified = Material("purified_fragment", "MinElute-purified tagmented DNA",
                        MaterialKind.MOLECULAR, state=_tagmented_template(), source=source)
    column_waste = Material("minelute_waste", "column flow-through: transposase and reaction mix",
                            MaterialKind.WASTE, source=source)
    graph.add_step(Step(
        "purify_tagmented", "Purify transposed DNA (MinElute; skipped for 500/5,000 cells)",
        StepKind.CLEANUP,
        (Input(tagged.id, InputRole.TEMPLATE), Input("minelute", InputRole.SELECTION)),
        (purified.id, column_waste.id), source=source, changes_molecule=False),
        (purified, column_waste))
    graph.mark_terminal(Terminal(column_waste.id, Disposition.DISCARD,
                                 "column flow-through and washes are discarded"))

    graph.add_input(Material("extension_polymerase", "NEBNext PCR master mix polymerase",
                             MaterialKind.CATALYST, source=source))
    extended = Material("extended_fragment", "gap-filled tagmented fragment",
                        MaterialKind.MOLECULAR, state=_tagmented_template(), source=source)
    graph.add_step(Step(
        "extend", "Extend the transposase-tagged ends, 72 °C 5 min", StepKind.EXTENSION,
        (Input(purified.id, InputRole.TEMPLATE),
         Input("extension_polymerase", InputRole.CATALYST)),
        (extended.id,), source=source), (extended,))

    # Supplementary Table 1 prints both primers unmodified; Ad2.1 stands for Ad2.1-Ad2.24.
    graph.add_input(Material("ad1_nomx", "Ad1_noMX", MaterialKind.OLIGO, source=source,
                             state=(seg("Ad1_noMX", ATAC_AD1_NOMX),)))
    graph.add_input(Material("ad2_indexed", "Ad2.1_TAAGGCGA", MaterialKind.OLIGO,
                             source=source, state=(seg("Ad2.1_TAAGGCGA", ATAC_AD2_1),)))
    geometry = PCRGeometry.exact(
        FINAL_LIBRARY, PCRPrimer("Ad1_noMX", ATAC_AD1_NOMX),
        PCRPrimer("Ad2.1_TAAGGCGA", ATAC_AD2_1),
        statement="Supplementary Table 1 prints both amplification primers; their 3-prime mosaic-end regions face inward across the accessible-DNA insert.")
    amplified = Material("amplified_library", "amplified indexed ATAC-seq library",
                         MaterialKind.MOLECULAR, state=FINAL_LIBRARY, source=source)
    graph.add_step(Step(
        "indexed_pcr", "Indexed PCR: 5 cycles, qPCR-chosen extra cycles (10–12 total)",
        StepKind.PCR,
        (Input(extended.id, InputRole.TEMPLATE),
         Input("ad1_nomx", InputRole.STRUCTURAL),
         Input("ad2_indexed", InputRole.STRUCTURAL),
         Input("extension_polymerase", InputRole.CATALYST)),
        (amplified.id,), geometry=geometry, source=source), (amplified,))
    graph.add_input(Material("pcr_cleanup", "Qiagen PCR cleanup column", MaterialKind.SELECTION,
                             source=source))
    library = Material("atac_library", "indexed ATAC-seq library", MaterialKind.LIBRARY,
                       state=FINAL_LIBRARY, source=source)
    pcr_waste = Material("pcr_cleanup_waste", "primers and PCR mix removed by the column",
                         MaterialKind.WASTE, source=source)
    graph.add_step(Step(
        "purify_library", "Purify the library (Qiagen PCR cleanup; no size selection)",
        StepKind.CLEANUP,
        (Input(amplified.id, InputRole.TEMPLATE), Input("pcr_cleanup", InputRole.SELECTION)),
        (library.id, pcr_waste.id), source=source, changes_molecule=False),
        (library, pcr_waste))
    graph.mark_terminal(Terminal(pcr_waste.id, Disposition.DISCARD,
                                 "column flow-through and washes are discarded"))
    graph.mark_terminal(Terminal(library.id, Disposition.SEQUENCING), SequencingPlan(
        FINAL_LIBRARY, tuple(SEQ_PRIMERS), tuple(p.role for p in SEQ_PRIMERS),
        read_lengths=READ_LENGTHS))

    graph.declare_reviewed_scope(dict((
            ("native_nuclei", "native_chromatin"),
            ("loaded_adapter_reagent", "loaded_mosaic_end"),
            ("tn5_reagent", "tn5"),
            ("tagmentation", "tagment"),
            ("minelute_reagent", "minelute"),
            ("minelute_cleanup", "purify_tagmented"),
            ("productive_and_nonproductive_ends", "non_library_tagmentation_products"),
            ("extension_polymerase_reagent", "extension_polymerase"),
            ("extension", "extend"),
            ("ad1_primer_reagent", "ad1_nomx"),
            ("ad2_primer_reagent", "ad2_indexed"),
            ("indexed_pcr", "indexed_pcr"),
            ("pcr_cleanup_reagent", "pcr_cleanup"),
            ("library_cleanup", "purify_library"),
            ("paired_end_sequencing", "atac_library"))), source)
    return graph.require_publishable()
