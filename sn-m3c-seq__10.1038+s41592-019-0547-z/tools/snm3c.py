"""Lee et al. 2019 sn-m3C-seq chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Segment
from restriction import DPNII
from single_cell_hic import ContactWorkflow
from pcr import PCRGeometry
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role
import snmc2_core
WORKFLOW = ContactWorkflow("sn-m3C-seq", DPNII, "after ligation", "none", "bisulfite-PCR")
TITLE = "sn-m3C-seq — contacts and DNA methylation from one nucleus"
NOTES = "01_sn-m3c-seq.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/s41592-019-0547-z">Lee et al., <i>Nature Methods</i> (2019)</a>.'
SUMMARY = "Single-nucleus chromosome-conformation capture is followed by bisulfite conversion and snmC-seq2-style indexed random-priming library construction."
STEPS = (("Digest and ligate nuclei", "Generate proximity-ligation products in fixed nuclei without a biotin pull-down."),
         ("Isolate one nucleus", "Sort nuclei so the contact and methylation observations retain a common cellular identity."),
         ("Bisulfite conversion", "Convert unmethylated cytosine to uracil while methylated cytosine remains protected."),
         ("Indexed random priming and PCR", "Copy converted fragments with an inline-indexed, partial Read 1 random primer inherited from snmC-seq2, add the opposite adapter and amplify."))
READOUT = "Paired bisulfite reads jointly encode the two contact fragments and cytosine-conversion state."
JUNCTION_CAPTION = "The unconverted contact architecture is shown here; after ligation, bisulfite treatment changes sequence content but not the contact boundary. ** marks ligation."
LIBRARY_CAPTION = "INFERRED — a representative snmC-seq2-family endpoint uses the source-printed P5L_AD002_H inline barcode and dual-index primer pair 1. Sample-specific primer/index assignments and Adaptase-derived bases remain unresolved in the inspected source scope."
LIBRARY_CAVEAT = "The sn-m3C Methods explicitly inherit the snmC-seq2 library procedure. Luo et al. Supplementary Methods print its random primers, and Supplementary Table 2 prints the indexed PCR primers. This panel uses one published representative from each family; This capped pass does not establish sample-specific assignments or vendor-derived bases; those remain unresolved. The inherited endpoint is provisional and its schematic index widths are not acquired cycle declarations."
SEQUENCING_UNAVAILABLE = "PE150 acquisition is source-stated. Both reads are computationally trimmed by 25 bases at the start and 3 bases at the end. Complete primer/Adaptase junctions have not been established in this capped pass, so no base-level read map is claimed."
READ_LENGTHS = {"Read 1": 150, "Read 2": 150}
REQUIRED_ROLES = ("Read 1", "Read 2")
def contact_product(): return WORKFLOW.junction()
RANDOM_PRIMER_TRIM_NT = 25
ADAPTASE_TAIL_TRIM_NT = 3
def final_library():
    inherited, _ = snmc2_core.library(
        "bisulfite-converted contact fragments", insert_nt=40,
        library_name="representative sn-m3C-seq bisulfite library",
        feature_id="snm3c_cell_inline",
    )
    # The inherited source leaves Adaptase's span unspecified. Give its
    # diagram a schematic width so the duplex renderer can pair placeholders.
    return Construct([
        Segment(segment.name, "X" * 8, segment.tag, placeholder=True,
                inferred=True, note="schematic width; actual Adaptase-added length unknown")
        if segment.name == "Adaptase-added bases" else segment
        for segment in inherited
    ], name=inherited.name)


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
    inactivation = g.input("dpnii_inactivation", "62 C for 10 min DpnII heat inactivation", MaterialKind.CONDITION)
    inactivated = g.product("inactivated", "DpnII-digested nuclei after enzyme heat inactivation", MaterialKind.MOLECULAR)
    g.react("inactivate_dpnii", "Heat-inactivate DpnII before contact ligation", StepKind.THERMAL_TREATMENT,
            ((cut, InputRole.TEMPLATE), (inactivation, InputRole.CONDITION)), (inactivated,), changes_molecule=False)
    ligase = g.input("contact_ligase", "T4 DNA ligase",
                     MaterialKind.CATALYST)
    contacts = g.product("contacts", "unselected proximity-ligation products",
                         MaterialKind.MOLECULAR, state=contact_product())
    g.react("contact_ligation", "Ligate DpnII contact ends without biotin capture",
            StepKind.PROXIMITY_LIGATION,
            ((inactivated, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (contacts,))
    isolation = g.input("isolation", "FANS singlet/Hoechst DNA-content gates, 384-well deposition", MaterialKind.CONDITION)
    nucleus = g.product("nucleus", "one ligated nucleus", MaterialKind.MOLECULAR,
                        state=contacts.state)
    excluded = g.product("excluded", "other ligated nuclei", MaterialKind.MOLECULAR)
    g.react("isolate", "Isolate one contact-bearing nucleus", StepKind.SAMPLING,
            ((contacts, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)),
            (nucleus, excluded), changes_molecule=False)
    g.mark_terminal(Terminal(excluded.id, Disposition.UNRESOLVED, "Other nuclei lie outside this sorted-nucleus lineage; their disposition is unrepresented."))
    protease = g.input("protease", "proteinase digestion reagent; conditions unresolved here", MaterialKind.CATALYST)
    released = g.product("released_dna", "contact DNA after sorted-nucleus proteinase digestion", MaterialKind.MOLECULAR)
    g.react("proteinase_digest", "Digest proteins before bisulfite conversion", StepKind.DEPROTEINIZATION,
            ((nucleus, InputRole.TEMPLATE), (protease, InputRole.CATALYST)), (released,))
    bisulfite = g.input("bisulfite", "bisulfite conversion reagent",
                        MaterialKind.CONDITION)
    converted = g.product("converted", "bisulfite-converted contact DNA",
                          MaterialKind.MOLECULAR)
    g.react("convert", "Bisulfite-convert contact DNA", StepKind.CHEMICAL_CONVERSION,
            ((released, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
            (converted,))
    primer = g.input("random_primer", "inline-indexed partial Read 1 random primer",
                     MaterialKind.OLIGO, inferred=True, source=LIBRARY_CAVEAT)
    polymerase = g.input("extension_polymerase", "random-priming polymerase",
                         MaterialKind.CATALYST)
    copied = g.product("copied", "indexed random-primed contact copies",
                       MaterialKind.MOLECULAR)
    g.react("random_prime", "Copy converted fragments from indexed random primers",
            StepKind.EXTENSION,
            ((converted, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
             (polymerase, InputRole.CATALYST)), (copied,))
    adaptase = g.input("adaptase", "Adaptase adapter-tagging reagent",
                       MaterialKind.CATALYST)
    opposite_adapter = g.input("opposite_adapter", "opposite adapter, sequence assignment unresolved",
                               MaterialKind.ADAPTER, inferred=True,
                               source=LIBRARY_CAVEAT)
    tagged = g.product("tagged", "Adaptase-tagged copies", MaterialKind.MOLECULAR,
                       state=None)
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
    g.sequence(library, SequencingPlan(acquired_roles=REQUIRED_ROLES, read_lengths=READ_LENGTHS,
                                       limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(tuple(g.steps) + tuple(g.materials), CITATION)
    return g.require_publishable()
