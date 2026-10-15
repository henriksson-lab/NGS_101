"""Bulk Methyl-HiC molecular model (Li et al. 2019)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from restriction import DPNII
from single_cell_hic import ContactWorkflow, SEQ_PRIMERS, inferred_truseq_library
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

WORKFLOW = ContactWorkflow("Methyl-HiC", DPNII, "bulk nuclei", "biotin-fill", "bisulfite-PCR", capture=True)
TITLE = "Methyl-HiC — chromatin contacts plus DNA methylation"
NOTES = "01_methyl-hic.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/s41592-019-0502-z">Li et al., <i>Nature Methods</i> (2019)</a>.'
SUMMARY = "In situ DpnII Hi-C preserves proximity junctions; on-bead adapter ligation followed by bisulfite conversion turns each paired-end contact molecule into a joint contact-and-methylation readout."
STEPS = (
    ("Digest nuclei and fill ends", "DpnII digestion is followed by fill-in with biotinylated nucleotides."),
    ("Proximity ligation", "T4 DNA ligase joins spatially adjacent filled restriction fragments in intact nuclei."),
    ("Shear and enrich", "Reverse crosslinks, sonicate to about 400 bp and capture biotin-marked junctions on streptavidin beads."),
    ("Ligate adapters, then convert", "Ligate sequencing adapters on beads before sodium-bisulfite conversion."),
    ("Uracil-tolerant PCR", "KAPA HiFi HotStart Uracil+ mix amplifies the converted library."),
)
READOUT = "Paired-end reads span the two DpnII-derived contact arms while C-to-T conversion reports cytosine methylation."
JUNCTION_CAPTION = "DpnII GATC ends are filled with biotinylated nucleotide and proximity-ligated; ** marks the contact junction."
LIBRARY_CAPTION = "Bisulfite-converted contact insert in the TruSeq-compatible shell named by the supplementary protocol. The particular indexed adapter member is not reported."
LIBRARY_CAVEAT = ("INFERRED ADAPTER BASES — the supplement names TruSeq PCR primers and "
                   "an indexed Illumina adapter, but does not print the ligated adapter "
                   "oligos or selected index. Canonical TruSeq-compatible arms therefore "
                   "show primer placement without claiming a source transcription.")
LIBRARY_ADAPTER_GEOMETRY = AdapterLigationGeometry.unavailable(
    "The supplementary protocol names indexed Illumina adapters and TruSeq PCR "
    "primers, but does not print the ligated adapter strands or their ligation-ready "
    "end geometry.")
LIBRARY_ADAPTER_CAPTION = "Source limitation before on-bead adapter ligation."
REQUIRED_ROLES = ("Read 1", "Index 1 (i7)", "Read 2")
ADAPTER_BEFORE_CONVERSION = True

def contact_product(): return WORKFLOW.junction()
def final_library():
    return inferred_truseq_library(contact_product(), "bisulfite-converted Methyl-HiC library")


def protocol_graph() -> ProtocolGraph:
    """The bulk DpnII/contact/capture/convert order in the defining protocol."""
    source = "https://doi.org/10.1038/s41592-019-0502-z"
    graph = ProtocolGraph(TITLE, authority=source)
    sample = graph.input(
        "sample", "input nuclei or cells", MaterialKind.BIOLOGICAL,
        state=molecule("chromatin sample", (role("chromosomal DNA in nuclei"),)),
        primary=True)
    nuclease = graph.input("restriction_enzyme", "DpnII restriction enzyme",
                           MaterialKind.CATALYST)
    digested = graph.product("digested", "DpnII-digested chromatin",
                             MaterialKind.MOLECULAR, state=sample.state)
    graph.react("digest", "Digest chromatin with DpnII", StepKind.DIGESTION,
                ((sample, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)),
                (digested,))

    nucleotide = graph.input("biotin_nucleotide", "biotinylated fill-in nucleotide",
                             MaterialKind.MOLECULAR)
    fill_enzyme = graph.input("fill_enzyme", "restriction-end fill-in polymerase",
                              MaterialKind.CATALYST)
    marked = graph.product("marked", "biotin-marked restriction ends",
                           MaterialKind.MOLECULAR, state=digested.state)
    graph.react("mark_ends", "Fill and biotin-mark restriction ends", StepKind.EXTENSION,
                ((digested, InputRole.TEMPLATE),
                 (nucleotide, InputRole.STRUCTURAL),
                 (fill_enzyme, InputRole.CATALYST)), (marked,))
    ligase = graph.input("contact_ligase", "proximity-ligation DNA ligase",
                         MaterialKind.CATALYST)
    contacts = graph.product("contacts", "proximity-ligated chromatin contacts",
                             MaterialKind.MOLECULAR, state=contact_product())
    graph.react("proximity_ligation", "Ligate spatially adjacent restriction fragments",
                StepKind.PROXIMITY_LIGATION,
                ((marked, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
                (contacts,))

    shear = graph.input("fragmentation_condition", "post-ligation fragmentation condition",
                        MaterialKind.CONDITION)
    fragments = graph.product("fragments", "fragmented proximity-ligation products",
                              MaterialKind.MOLECULAR, state=contacts.state)
    graph.react("fragment", "Fragment contact products for capture",
                StepKind.FRAGMENTATION,
                ((contacts, InputRole.TEMPLATE), (shear, InputRole.CONDITION)),
                (fragments,))
    beads = graph.input("streptavidin", "streptavidin capture beads",
                        MaterialKind.SELECTION)
    captured = graph.product("captured", "biotin-enriched contact fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("uncaptured", "contact-poor unbound fragments",
                             MaterialKind.WASTE)
    graph.react("capture", "Capture biotin-bearing contact junctions",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not retain a biotin-marked contact junction")

    adapter = graph.input("library_adapter", "indexed Illumina library adapter",
                          MaterialKind.ADAPTER, inferred=True,
                          source=source + " (ligated adapter strands not disclosed)")
    adapter_ligase = graph.input("adapter_ligase", "library-adapter ligase",
                                 MaterialKind.CATALYST)
    ligated = graph.product("adapter_ligated", "adapter-ligated contact DNA",
                            MaterialKind.MOLECULAR, state=captured.state)
    graph.react("ligate_library_adapter", "Install sequencing-library adapters",
                StepKind.LIGATION,
                ((captured, InputRole.TEMPLATE),
                 (adapter, InputRole.STRUCTURAL),
                 (adapter_ligase, InputRole.CATALYST)), (ligated,),
                geometry=LIBRARY_ADAPTER_GEOMETRY)
    bisulfite = graph.input("bisulfite", "sodium-bisulfite conversion condition",
                            MaterialKind.CONDITION)
    converted = graph.product("converted", "bisulfite-converted contact DNA",
                              MaterialKind.MOLECULAR, state=ligated.state)
    graph.react("convert", "Bisulfite-convert contact DNA",
                StepKind.CHEMICAL_CONVERSION,
                ((ligated, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
                (converted,))

    forward = graph.input("pcr_forward", "TruSeq library PCR forward primer",
                          MaterialKind.OLIGO, inferred=True,
                          source=source + " (primer sequence not printed)")
    reverse = graph.input("pcr_reverse", "TruSeq library PCR reverse primer",
                          MaterialKind.OLIGO, inferred=True,
                          source=source + " (primer sequence not printed)")
    polymerase = graph.input("pcr_polymerase", "uracil-tolerant KAPA HiFi polymerase",
                             MaterialKind.CATALYST)
    final = final_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY,
                            state=final, inferred=True,
                            source=source + " (indexed adapter member not printed)")
    graph.react("library_pcr", "Amplify the final contact library", StepKind.PCR,
                ((converted, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The source does not print the PCR primers or indexed adapter strands."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=SEQ_PRIMERS, acquired_roles=REQUIRED_ROLES,
        limitation="The indexed adapter member and ligated end geometry are not published."))
    graph.declare_reviewed_scope((
        "sample", "restriction_enzyme", "digest", "biotin_nucleotide", "fill_enzyme",
        "mark_ends", "contact_ligase", "proximity_ligation", "fragmentation_condition",
        "fragment", "streptavidin", "capture", "library_adapter", "adapter_ligase",
        "ligate_library_adapter", "bisulfite", "convert", "pcr_forward", "pcr_reverse",
        "pcr_polymerase", "library_pcr", "library"), source)
    return graph.require_publishable()
