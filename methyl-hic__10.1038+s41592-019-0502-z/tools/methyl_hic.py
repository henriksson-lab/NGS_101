"""Bulk Methyl-HiC molecular model (Li et al. 2019)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Segment
from restriction import DPNII
from single_cell_hic import ContactWorkflow, SEQ_PRIMERS, inferred_truseq_library
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind, Terminal)
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
LIBRARY_CAPTION = "Bisulfite-converted contact insert in an inferred TruSeq-compatible shell. The particular indexed adapter member is unresolved from the checked sources."
LIBRARY_CAVEAT = ("INFERRED ADAPTER BASES — the supplement names TruSeq PCR primers and "
                   "an indexed Illumina adapter, but does not print the ligated adapter "
                   "oligos or selected index. Canonical TruSeq-compatible arms therefore "
                   "show primer placement without claiming a source transcription. The "
                   "article's explicit TruSeq adapter is for a separate WGBS control. "
                   "The bulk source and linked-material review remains in progress.")
LIBRARY_ADAPTER_GEOMETRY = AdapterLigationGeometry.unavailable(
    "The supplementary protocol names indexed Illumina adapters and TruSeq PCR "
    "primers, but does not print the ligated adapter strands or their ligation-ready "
    "end geometry.")
LIBRARY_ADAPTER_CAPTION = "Source limitation before on-bead adapter ligation."
REQUIRED_ROLES = ("Read 1", "Index 1 (i7)", "Read 2")
ADAPTER_BEFORE_CONVERSION = True
LIBRARY_CAVEAT += (" Published enzyme-label difference: the main Methyl-HiC Methods and "
    "supplementary reagent list specify DpnII; the shared supplementary protocol's "
    "Biotin Fill-in step1 says Inactivate MboI. The model retains DpnII from the main "
    "Methods and leaves the intended heat-inactivation instruction unresolved. "
    "Bulk post-conversion handling is also ambiguous: the protocol transfers supernatant "
    "to a purification column but later says PCR directly on beads. This graph follows "
    "the article's purified converted-DNA route, without asserting a resolved bead/column mechanism.")

def contact_product(): return WORKFLOW.junction()
def final_library():
    converted = Construct([
        Segment(s.name, s.top if s.placeholder else s.top.replace("C", "Y"), s.tag,
                placeholder=s.placeholder or "C" in s.top, inferred=not s.placeholder,
                note="C/T class after conversion and PCR; methylation assignment unresolved")
        for s in contact_product()
    ], name="converted contact sequence classes")
    return inferred_truseq_library(converted, "bisulfite-converted Methyl-HiC library")


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

    nucleotide = graph.input("biotin_nucleotide", "biotin-14-dATP fill-in nucleotide",
                             MaterialKind.MOLECULAR)
    fill_enzyme = graph.input("fill_enzyme", "Klenow DNA polymerase I large fragment",
                              MaterialKind.CATALYST)
    marked = graph.product("marked", "biotin-marked restriction ends",
                           MaterialKind.MOLECULAR, state=digested.state)
    graph.react("mark_ends", "Fill and biotin-mark restriction ends", StepKind.EXTENSION,
                ((digested, InputRole.TEMPLATE),
                 (nucleotide, InputRole.STRUCTURAL),
                 (fill_enzyme, InputRole.CATALYST)), (marked,))
    ligase = graph.input("contact_ligase", "T4 DNA ligase",
                         MaterialKind.CATALYST)
    contacts = graph.product("contacts", "proximity-ligated chromatin contacts",
                             MaterialKind.MOLECULAR, state=contact_product())
    graph.react("proximity_ligation", "Ligate spatially adjacent restriction fragments",
                StepKind.PROXIMITY_LIGATION,
                ((marked, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
                (contacts,))

    protease = graph.input("proteinase", "Proteinase K/SDS and reverse-crosslinking conditions", MaterialKind.CATALYST)
    released = graph.product("released", "reverse-crosslinked contact DNA", MaterialKind.MOLECULAR)
    graph.react("reverse_crosslinks", "Proteinase treatment and reverse crosslinks before DNA shearing", StepKind.DEPROTEINIZATION,
                ((contacts, InputRole.TEMPLATE), (protease, InputRole.CATALYST)), (released,))
    ethanol = graph.input("ethanol", "ethanol/sodium acetate DNA precipitation", MaterialKind.SELECTION)
    recovered = graph.product("recovered", "precipitated contact DNA", MaterialKind.MOLECULAR)
    precip_waste = graph.product("precip_waste", "unrecovered precipitation material", MaterialKind.WASTE)
    graph.react("precipitate", "Recover contact DNA by ethanol precipitation", StepKind.CLEANUP,
                ((released, InputRole.TEMPLATE), (ethanol, InputRole.SELECTION)), (recovered, precip_waste))
    graph.discard(precip_waste, "not recovered with the DNA pellet")
    shear = graph.input("fragmentation_condition", "Covaris sonication to approximately 400 bp",
                        MaterialKind.CONDITION)
    fragments = graph.product("fragments", "fragmented proximity-ligation products",
                              MaterialKind.MOLECULAR, state=contacts.state)
    graph.react("fragment", "Fragment contact products for capture",
                StepKind.FRAGMENTATION,
                ((recovered, InputRole.TEMPLATE), (shear, InputRole.CONDITION)),
                (fragments,))
    spri = graph.input("spri", "SPRI beads, 1x post-sonication cleanup", MaterialKind.SELECTION)
    cleaned = graph.product("cleaned", "SPRI-cleaned sonicated contact DNA", MaterialKind.MOLECULAR)
    spri_waste = graph.product("spri_waste", "unrecovered SPRI cleanup material", MaterialKind.WASTE)
    graph.react("spri_cleanup", "Purify sonicated products with 1x SPRI beads", StepKind.CLEANUP,
                ((fragments, InputRole.TEMPLATE), (spri, InputRole.SELECTION)), (cleaned, spri_waste))
    graph.discard(spri_waste, "not recovered in SPRI eluate")
    beads = graph.input("streptavidin", "streptavidin capture beads",
                        MaterialKind.SELECTION)
    captured = graph.product("captured", "biotin-enriched contact fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("uncaptured", "capture supernatant retained as optional methylation background control",
                             MaterialKind.MOLECULAR)
    graph.react("capture", "Capture biotin-bearing contact junctions",
                StepKind.AFFINITY_CAPTURE,
                ((cleaned, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, rejected))
    graph.mark_terminal(Terminal(rejected.id, Disposition.STORED, "Optional background-control route selected: protocol retains and ethanol-precipitates the capture supernatant; its subsequent library preparation is outside this graph."))

    repair = graph.input("repair", "T4 PNK, T4 DNA polymerase and Klenow end-repair mix", MaterialKind.CATALYST)
    repaired = graph.product("repaired", "end-repaired bead-bound DNA", MaterialKind.MOLECULAR)
    graph.react("end_repair", "Repair sheared ends and remove biotin from unligated ends", StepKind.END_REPAIR,
                ((captured, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("tailing", "Klenow exo-minus with dATP", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "dA-tailed bead-bound DNA", MaterialKind.MOLECULAR)
    graph.react("da_tail", "Add dA tails before adapter ligation", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("library_adapter", "indexed Illumina library adapter",
                          MaterialKind.ADAPTER, inferred=True,
                          source=source + " (ligated adapter strands not disclosed)")
    adapter_ligase = graph.input("adapter_ligase", "DNA Quick Ligase",
                                 MaterialKind.CATALYST)
    ligated = graph.product("adapter_ligated", "adapter-ligated contact DNA",
                            MaterialKind.MOLECULAR)
    graph.react("ligate_library_adapter", "Install sequencing-library adapters",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE),
                 (adapter, InputRole.STRUCTURAL),
                 (adapter_ligase, InputRole.CATALYST)), (ligated,),
                geometry=LIBRARY_ADAPTER_GEOMETRY)
    bisulfite = graph.input("bisulfite", "sodium-bisulfite conversion condition",
                            MaterialKind.CONDITION)
    lambda_control = graph.input("lambda_control", "0.5% same-adapter-ligated unmethylated lambda DNA", MaterialKind.MOLECULAR)
    converted = graph.product("converted", "bisulfite-converted contact DNA",
                              MaterialKind.MOLECULAR)
    graph.react("convert", "Bisulfite-convert contact DNA",
                StepKind.CHEMICAL_CONVERSION,
                ((ligated, InputRole.TEMPLATE), (lambda_control, InputRole.SPIKE_IN), (bisulfite, InputRole.CONDITION)),
                (converted,))

    purified = graph.product("purified_converted", "purified bisulfite-converted library", MaterialKind.MOLECULAR)
    graph.react("conversion_recovery", "Recover converted DNA for uracil-tolerant PCR", StepKind.MAGIC,
                ((converted, InputRole.TEMPLATE),), (purified,),
                reason="Main Methods use purified converted DNA; supplement directs supernatant to a column but subsequently mentions bead PCR. Exact transfer/retention geometry is unresolved.")
    forward = graph.input("pcr_forward", "TruSeq library PCR forward primer",
                          MaterialKind.OLIGO, inferred=True,
                          source=source + " (primer sequence not printed)")
    reverse = graph.input("pcr_reverse", "TruSeq library PCR reverse primer",
                          MaterialKind.OLIGO, inferred=True,
                          source=source + " (primer sequence not printed)")
    polymerase = graph.input("pcr_polymerase", "uracil-tolerant KAPA HiFi polymerase",
                             MaterialKind.CATALYST)
    final = final_library()
    library = graph.product("amplified", final.name, MaterialKind.MOLECULAR,
                            state=final, inferred=True,
                            source=source + " (indexed adapter member not printed)")
    graph.react("library_pcr", "Amplify the final contact library", StepKind.PCR,
                ((purified, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The source does not print the PCR primers or indexed adapter strands."))
    size_rule = graph.input("size_rule", "double AMPure selection, mainly 400–700 bp", MaterialKind.SELECTION)
    selected = graph.product("library", final.name, MaterialKind.LIBRARY, state=final, inferred=True)
    offsize = graph.product("offsize", "products outside selected library-size interval", MaterialKind.WASTE)
    graph.react("size_selection", "Double-size-select amplified libraries with AMPure beads", StepKind.SELECTION,
                ((library, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)), (selected, offsize))
    graph.discard(offsize, "not recovered in selected library interval")
    graph.sequence(selected, SequencingPlan(
        library=final, primers=SEQ_PRIMERS, acquired_roles=REQUIRED_ROLES,
        limitation="The indexed adapter member and ligated end geometry are unresolved from the checked Li and Rao sources; independent source review continues."))
    graph.declare_reviewed_scope((
        "sample", "restriction_enzyme", "digest", "biotin_nucleotide", "fill_enzyme",
        "mark_ends", "contact_ligase", "proximity_ligation", "fragmentation_condition",
        "fragment", "streptavidin", "capture", "library_adapter", "adapter_ligase",
        "ligate_library_adapter", "bisulfite", "convert", "pcr_forward", "pcr_reverse",
        "pcr_polymerase", "library_pcr", "library", "proteinase", "reverse_crosslinks",
        "ethanol", "precipitate", "spri", "spri_cleanup", "repair", "end_repair", "tailing", "da_tail",
        "lambda_control", "conversion_recovery", "size_rule", "size_selection"), source)
    return graph.require_publishable()
