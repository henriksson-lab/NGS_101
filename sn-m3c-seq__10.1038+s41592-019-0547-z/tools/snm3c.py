"""Lee et al. 2019 sn-m3C-seq chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Segment, feature
from restriction import DPNII
from single_cell_hic import ContactWorkflow, unresolved_illumina_library
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
WORKFLOW = ContactWorkflow("sn-m3C-seq", DPNII, "after ligation", "none", "bisulfite-PCR")
TITLE = "sn-m3C-seq — contacts and DNA methylation from one nucleus"
NOTES = "01_sn-m3c-seq.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/s41592-019-0547-z">Lee et al., <i>Nature Methods</i> (2019)</a>.'
SUMMARY = "Single-nucleus chromosome-conformation capture is followed by bisulfite conversion and snmC-seq2-style indexed random-priming library construction."
STEPS = (("Digest and ligate nuclei", "Generate proximity-ligation products in fixed nuclei without a biotin pull-down."),
         ("Isolate one nucleus", "Sort nuclei so the contact and methylation observations retain a common cellular identity."),
         ("Bisulfite conversion", "Convert unmethylated cytosine to uracil while methylated cytosine remains protected."),
         ("Indexed random priming and PCR", "Copy converted fragments with an indexed P5-bearing random primer, add the opposite adapter and amplify."))
READOUT = "Paired bisulfite reads jointly encode the two contact fragments and cytosine-conversion state."
JUNCTION_CAPTION = "The unconverted contact architecture is shown here; after ligation, bisulfite treatment changes sequence content but not the contact boundary. ** marks ligation."
LIBRARY_CAPTION = "Bisulfite-compatible snmC-seq2 product. The 25-base random-primer-derived prefix and 3-base Adaptase tail are experimentally trimmed; the sn-m3C paper does not print the complete adapter oligos."
LIBRARY_CAVEAT = "INFERRED — the trim-defined protocol-specific regions are established, but the complete outer adapter bases are not reported in the defining sn-m3C source."
SEQUENCING_UNAVAILABLE = "The source specifies paired-end reads and trimming 25 bases from the start and 3 bases from the end of each read, but does not print a complete sequencing-primer set on which base-level placement can be verified."
def contact_product(): return WORKFLOW.junction()
RANDOM_PRIMER_TRIM_NT = 25
ADAPTASE_TAIL_TRIM_NT = 3
def final_library():
    converted = Construct([
        Segment("indexed random-primer-derived sequence", "N"*RANDOM_PRIMER_TRIM_NT,
                "cbc", placeholder=True,
                feature=feature("cell_index", "cell_barcode", "unknown",
                                note="exact barcode sub-boundary is unresolved")),
        *list(contact_product()),
        Segment("Adaptase low-complexity tail", "N"*ADAPTASE_TAIL_TRIM_NT,
                "me", placeholder=True),
    ], name="converted and tagged contact")
    return unresolved_illumina_library(converted, "sn-m3C-seq bisulfite library", indexed=True)


def protocol_graph() -> ProtocolGraph:
    """Contact ligation, nucleus isolation, then bisulfite-compatible readout."""
    g = ProtocolGraph(TITLE, authority=CITATION)
    sample = g.input("sample", "fixed chromatin nuclei", MaterialKind.BIOLOGICAL,
                     state=molecule("fixed chromatin", (role("chromosomal DNA"),)),
                     primary=True)
    dpnii = g.input("dpnii", "DpnII", MaterialKind.CATALYST)
    cut = g.product("cut", "DpnII-cut chromatin", MaterialKind.MOLECULAR,
                    state=sample.state, chemical_state=WORKFLOW.digest())
    g.react("digest", "Digest fixed chromatin with DpnII", StepKind.DIGESTION,
            ((sample, InputRole.TEMPLATE), (dpnii, InputRole.CATALYST)), (cut,))
    ligase = g.input("contact_ligase", "contact-ligation DNA ligase",
                     MaterialKind.CATALYST)
    contacts = g.product("contacts", "unselected proximity-ligation products",
                         MaterialKind.MOLECULAR, state=contact_product())
    g.react("contact_ligation", "Ligate DpnII contact ends without biotin capture",
            StepKind.PROXIMITY_LIGATION,
            ((cut, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (contacts,))
    isolation = g.input("isolation", "single-nucleus sorting", MaterialKind.CONDITION)
    nucleus = g.product("nucleus", "one ligated nucleus", MaterialKind.MOLECULAR,
                        state=contacts.state)
    excluded = g.product("excluded", "other ligated nuclei", MaterialKind.WASTE)
    g.react("isolate", "Isolate one contact-bearing nucleus", StepKind.SAMPLING,
            ((contacts, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)),
            (nucleus, excluded), changes_molecule=False)
    g.discard(excluded, "outside the retained single-nucleus library")
    bisulfite = g.input("bisulfite", "bisulfite conversion reagent",
                        MaterialKind.CONDITION)
    converted = g.product("converted", "bisulfite-converted contact DNA",
                          MaterialKind.MOLECULAR, state=nucleus.state)
    g.react("convert", "Bisulfite-convert contact DNA", StepKind.CHEMICAL_CONVERSION,
            ((nucleus, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
            (converted,))
    primer = g.input("random_primer", "indexed P5-bearing random primer",
                     MaterialKind.OLIGO, inferred=True, source=LIBRARY_CAVEAT)
    polymerase = g.input("extension_polymerase", "random-priming polymerase",
                         MaterialKind.CATALYST)
    copied = g.product("copied", "indexed random-primed contact copies",
                       MaterialKind.MOLECULAR, state=converted.state)
    g.react("random_prime", "Copy converted fragments from indexed random primers",
            StepKind.EXTENSION,
            ((converted, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
             (polymerase, InputRole.CATALYST)), (copied,))
    adaptase = g.input("adaptase", "Adaptase adapter-tagging reagent",
                       MaterialKind.CATALYST)
    opposite_adapter = g.input("opposite_adapter", "unreported opposite adapter",
                               MaterialKind.ADAPTER, inferred=True,
                               source=LIBRARY_CAVEAT)
    tagged = g.product("tagged", "Adaptase-tagged copies", MaterialKind.MOLECULAR,
                       state=copied.state)
    g.react("adaptase_tag", "Add the Adaptase low-complexity tail and opposite adapter",
            StepKind.ADAPTER_ATTACHMENT,
            ((copied, InputRole.TEMPLATE), (opposite_adapter, InputRole.STRUCTURAL),
             (adaptase, InputRole.CATALYST)),
            (tagged,))
    forward = g.input("pcr_forward", "library PCR forward primer", MaterialKind.OLIGO,
                      inferred=True, source=LIBRARY_CAVEAT)
    reverse = g.input("pcr_reverse", "library PCR reverse primer", MaterialKind.OLIGO,
                      inferred=True, source=LIBRARY_CAVEAT)
    pcr = g.input("pcr_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    final = final_library()
    library = g.product("library", final.name, MaterialKind.LIBRARY,
                        state=final, inferred=True, source=LIBRARY_CAVEAT)
    g.react("library_pcr", "Amplify the snmC-seq2-style indexed library", StepKind.PCR,
            ((tagged, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)),
            (library,), geometry=PCRGeometry.unavailable(
                final, reason=SEQUENCING_UNAVAILABLE))
    g.sequence(library, SequencingPlan(acquired_roles=("Read 1", "Read 2"),
                                       limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(tuple(g.steps) + tuple(g.materials), CITATION)
    return g.require_publishable()
