"""Reusable, deliberately small model for single-cell chromosome-conformation pages.

The workflow object is the policy boundary: it rejects combinations that are not real
protocol classes (for example, streptavidin capture without biotin marking), while each
paper-specific module remains responsible for its enzyme, ordering and amplification.
"""
from __future__ import annotations

from dataclasses import dataclass

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, MolecularRelation, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role
from restriction import Digest, RestrictionEnzyme


@dataclass(frozen=True)
class ContactWorkflow:
    name: str
    enzyme: RestrictionEnzyme
    isolate: str
    marking: str                 # biotin-fill, biotin-blunt, none, or bridge-adaptor
    amplification: str           # PCR, MDA, META, or bisulfite-PCR
    barcoding_rounds: int = 0
    capture: bool = False

    def __post_init__(self) -> None:
        if self.marking not in {"biotin-fill", "biotin-blunt", "none", "bridge-adaptor"}:
            raise ValueError(f"unknown end-marking mode: {self.marking}")
        if self.amplification not in {"PCR", "MDA", "META", "bisulfite-PCR"}:
            raise ValueError(f"unknown amplification mode: {self.amplification}")
        if self.capture and not self.marking.startswith("biotin"):
            raise ValueError("streptavidin capture requires a biotin-marked workflow")
        if self.marking == "bridge-adaptor" and self.barcoding_rounds != 2:
            raise ValueError("indexed bridge-adaptor Hi-C requires two barcode rounds")

    def digest(self) -> Digest:
        return Digest(self.enzyme, "X" * 16, "X" * 16)

    def junction(self) -> Construct:
        d = self.digest()
        if self.marking == "biotin-fill":
            middle = d.fill_in(biotin_base="A").junction
        elif self.enzyme.end == "blunt":
            middle = self.enzyme.site
        else:
            # Cohesive ends re-form one recognition site when compatible ends ligate.
            middle = self.enzyme.site
        return Construct([
            Segment("locus A", "X" * 22, placeholder=True),
            Segment("contact junction", middle, tag="me"),
            Segment("locus B", "X" * 22, placeholder=True),
        ], name=f"{self.name} contact product")


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


def inferred_truseq_library(insert: Construct, name: str) -> Construct:
    """Canonical paired-end shell when a paper names, but does not print, the kit arms."""
    lib = Construct([
        Segment("P5", il.P5, "p5", inferred=True),
        Segment("Read 1 arm", il.TRUSEQ_READ1[4:], "r1", inferred=True),
        *list(insert), Segment("dA junction", "A", inferred=True),
        Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=True),
        Segment("i7 reverse complement", "I" * 8, "cbc", placeholder=True, inferred=True,
                feature=feature("sample_i7", "sample_index", "unknown")),
        Segment("P7 reverse complement", il.P7_RC, "p7", inferred=True),
    ], name=name)
    errors = sp.verify(lib, SEQ_PRIMERS,
                       required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if errors:
        raise ValueError("invalid Hi-C library: " + "; ".join(errors))
    return lib


def unresolved_illumina_library(insert: Construct, name: str, *, indexed: bool | None) -> Construct:
    """Library whose source identifies Illumina prep but not the adapter oligos.

    Role-labelled placeholders deliberately prevent a named kit generation, index count,
    or modern primer set from being smuggled into a historical protocol.
    """
    index = ("indexed adapter (sequence not reported)" if indexed is True else
             "adapter/index region (not reported)" if indexed is None else
             "adapter arm (sequence not reported)")
    semantic = (feature("sample_index_region", "sample_index", "unknown",
                        note="index sub-boundary is not reported")
                if indexed is True else None)
    return Construct([
        Segment("Illumina left arm (not reported)", "X" * 20,
                placeholder=True, inferred=True),
        *list(insert),
        Segment(index, "X" * 20, placeholder=True, inferred=True, feature=semantic),
    ], name=name)


# Nagano 2013 and Stevens 2017 print this historical paired-end primer in their
# supplementary methods.  It is not the modern TruSeq Read 2 primer.
HISTORICAL_PE_READ2 = "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"


def historical_inline_pe_library(insert: Construct, name: str, barcode: str = "CAA") -> Construct:
    """Exact custom 3-bp inline-tag architecture shared by Nagano and Stevens."""
    if len(barcode) != 3 or any(b not in "ACGT" for b in barcode):
        raise ValueError("historical Hi-C identification tag must be three DNA bases")
    return Construct([
        Segment("P5", il.P5, "p5"),
        Segment("historical Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        Segment("3-bp identification tag", barcode, "cbc",
                feature=feature("cell_inline_a", "cell_barcode", "whitelist",
                                whitelist="published identification-tag set", group="cell_id", part="read 1 end")),
        *list(insert),
        Segment("opposite 3-bp identification tag", revcomp(barcode), "cbc",
                feature=feature("cell_inline_b", "cell_barcode", "whitelist",
                                whitelist="published identification-tag set", group="cell_id", part="read 2 end")),
        Segment("historical Read 2 site", revcomp(HISTORICAL_PE_READ2), "r2"),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name=name)


def protocol_graph(module) -> ProtocolGraph:
    """Build the typed process implied by a paper module's ``ContactWorkflow``.

    The paper modules remain responsible for source-specific restriction enzymes,
    isolation timing, marking, capture, amplification, final construct and caveats.
    This constructor owns the shared causal order and terminal accounting.
    """
    design = module.WORKFLOW
    source = module.CITATION
    graph = ProtocolGraph(module.TITLE, authority=source)
    sample = graph.input("sample", "input nuclei or cells", MaterialKind.BIOLOGICAL,
                         state=molecule("chromatin sample", (role("chromosomal DNA in nuclei"),)),
                         primary=True)
    reviewed = ["sample"]

    def isolate(template, prefix):
        condition = graph.input(prefix + "_condition", design.isolate,
                                MaterialKind.CONDITION)
        retained = graph.product(prefix, "isolated single nucleus or nuclear fraction",
                                 MaterialKind.MOLECULAR, state=template.state)
        excluded = graph.product(prefix + "_excluded", "excluded nuclei or bulk material",
                                 MaterialKind.WASTE)
        graph.react(prefix + "_sampling", "Isolate the protocol-defined nuclear unit",
                    StepKind.SAMPLING,
                    ((template, InputRole.TEMPLATE), (condition, InputRole.CONDITION)),
                    (retained, excluded), changes_molecule=False)
        graph.discard(excluded, "is outside the retained single-nucleus or nuclear fraction")
        reviewed.extend((condition.id, prefix + "_sampling"))
        return retained

    current = sample
    before = "before" in design.isolate.lower()
    bulk = "bulk" in design.isolate.lower()
    if before and not bulk:
        current = isolate(current, "isolated_nucleus")

    nuclease = graph.input("restriction_enzyme", design.enzyme.name + " restriction enzyme",
                           MaterialKind.CATALYST)
    digested = graph.product("digested", design.enzyme.name + "-digested chromatin",
                             MaterialKind.MOLECULAR, state=current.state,
                             chemical_state=getattr(module, "DEMO", None))
    graph.react("digest", "Digest chromatin with " + design.enzyme.name,
                StepKind.DIGESTION,
                ((current, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)),
                (digested,))
    reviewed.extend(("restriction_enzyme", "digest"))
    current = digested

    if design.marking.startswith("biotin"):
        nucleotide = graph.input("biotin_nucleotide", "biotinylated fill-in nucleotide",
                                 MaterialKind.MOLECULAR)
        fill_enzyme = graph.input("fill_enzyme", "restriction-end fill-in polymerase",
                                  MaterialKind.CATALYST)
        marked = graph.product("marked", "biotin-marked restriction ends",
                               MaterialKind.MOLECULAR, state=current.state,
                               chemical_state=getattr(module, "FILLED", None))
        graph.react("mark_ends", "Fill and biotin-mark restriction ends",
                    StepKind.EXTENSION,
                    ((current, InputRole.TEMPLATE), (nucleotide, InputRole.STRUCTURAL),
                     (fill_enzyme, InputRole.CATALYST)), (marked,),
                    relations=(MolecularRelation.RESTRICTION_FILL,)
                    if marked.chemical_state is not None else ())
        reviewed.extend(("biotin_nucleotide", "fill_enzyme", "mark_ends"))
        current = marked
    elif design.marking == "bridge-adaptor":
        bridge = graph.input("bridge_adapter", "indexed contact bridge adapter",
                             MaterialKind.ADAPTER)
        bridge_ligase = graph.input("bridge_ligase", "bridge-adapter ligase",
                                    MaterialKind.CATALYST)
        bridged = graph.product("bridged", "bridge-adapter-bearing restriction fragments",
                                MaterialKind.MOLECULAR, state=current.state)
        graph.react("install_bridge", "Ligate indexed bridge adapters",
                    StepKind.LIGATION,
                    ((current, InputRole.TEMPLATE), (bridge, InputRole.STRUCTURAL),
                     (bridge_ligase, InputRole.CATALYST)), (bridged,),
                    geometry=AdapterLigationGeometry.unavailable(
                        "The paper-specific bridge geometry is represented by the final contact construct, not a universal adapter end."))
        reviewed.extend(("bridge_adapter", "bridge_ligase", "install_bridge"))
        current = bridged

    contact_ligase = graph.input("contact_ligase", "proximity-ligation DNA ligase",
                                 MaterialKind.CATALYST)
    contacts = graph.product("contacts", "proximity-ligated chromatin contacts",
                             MaterialKind.MOLECULAR, state=module.contact_product())
    graph.react("proximity_ligation", "Ligate spatially adjacent restriction fragments",
                StepKind.PROXIMITY_LIGATION,
                ((current, InputRole.TEMPLATE), (contact_ligase, InputRole.CATALYST)),
                (contacts,), relations=(MolecularRelation.CONTACT_JUNCTION,)
                if getattr(module, "FILLED", None) is not None else ())
    reviewed.extend(("contact_ligase", "proximity_ligation"))
    current = contacts

    if not before and not bulk:
        current = isolate(current, "isolated_nucleus")

    if design.capture:
        shear = graph.input("fragmentation_condition", "post-ligation fragmentation condition",
                            MaterialKind.CONDITION)
        fragments = graph.product("fragments", "fragmented proximity-ligation products",
                                  MaterialKind.MOLECULAR, state=current.state)
        graph.react("fragment", "Fragment contact products for capture",
                    StepKind.FRAGMENTATION,
                    ((current, InputRole.TEMPLATE), (shear, InputRole.CONDITION)),
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
        reviewed.extend(("fragmentation_condition", "fragment", "streptavidin", "capture"))
        current = captured

    if design.amplification == "MDA":
        phi29 = graph.input("phi29", "phi29 strand-displacing polymerase",
                            MaterialKind.CATALYST)
        amplified = graph.product("amplified", "MDA-amplified contact DNA",
                                  MaterialKind.MOLECULAR, state=current.state)
        graph.react("amplify", "Multiple-displacement amplify the contact mixture",
                    StepKind.EXTENSION,
                    ((current, InputRole.TEMPLATE), (phi29, InputRole.CATALYST)),
                    (amplified,))
        reviewed.extend(("phi29", "amplify")); current = amplified
    elif design.amplification == "META":
        adapter = graph.input("meta_adapter", "META transposon adapter",
                              MaterialKind.ADAPTER, inferred=True,
                              source=source + " (transposon bases unavailable)")
        transposase = graph.input("meta_transposase", "META transposase",
                                  MaterialKind.CATALYST)
        tagged = graph.product("meta_tagged", "META-tagged contact DNA",
                               MaterialKind.MOLECULAR, state=current.state)
        waste = graph.product("meta_waste", "unusable META tagging products",
                              MaterialKind.WASTE)
        graph.react("meta_tag", "Tag contact molecules for META amplification",
                    StepKind.TAGMENTATION,
                    ((current, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                     (transposase, InputRole.CATALYST)), (tagged, waste))
        graph.discard(waste, "does not yield an amplifiable META-tagged molecule")
        reviewed.extend(("meta_adapter", "meta_transposase", "meta_tag")); current = tagged

    def install_library_adapter(template):
        library_adapter = graph.input(
            "library_adapter", "paper-defined Illumina library adapter",
            MaterialKind.ADAPTER,
            inferred=not bool(getattr(module, "SEQ_PRIMERS", ())), source=source)
        adapter_ligase = graph.input("adapter_ligase", "library-adapter ligase",
                                     MaterialKind.CATALYST)
        ligated_product = graph.product("adapter_ligated", "adapter-ligated contact DNA",
                                        MaterialKind.MOLECULAR, state=template.state)
        graph.react("ligate_library_adapter", "Install sequencing-library adapters",
                    StepKind.LIGATION,
                    ((template, InputRole.TEMPLATE),
                     (library_adapter, InputRole.STRUCTURAL),
                     (adapter_ligase, InputRole.CATALYST)), (ligated_product,),
                    geometry=AdapterLigationGeometry.unavailable(
                        getattr(module, "LIBRARY_CAVEAT",
                                "The source does not define one universal ligation-end geometry.")
                        or "The paper-specific adapter is represented in the final construct."))
        reviewed.extend(("library_adapter", "adapter_ligase", "ligate_library_adapter"))
        return ligated_product

    adapter_before_conversion = bool(getattr(module, "ADAPTER_BEFORE_CONVERSION", False))
    if adapter_before_conversion:
        current = install_library_adapter(current)

    if design.amplification == "bisulfite-PCR":
        bisulfite = graph.input("bisulfite", "sodium-bisulfite conversion condition",
                                MaterialKind.CONDITION)
        converted = graph.product("converted", "bisulfite-converted contact DNA",
                                  MaterialKind.MOLECULAR, state=current.state)
        graph.react("convert", "Bisulfite-convert contact DNA",
                    StepKind.CHEMICAL_CONVERSION,
                    ((current, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
                    (converted,))
        reviewed.extend(("bisulfite", "convert")); current = converted

    if not adapter_before_conversion:
        current = install_library_adapter(current)

    forward = graph.input("pcr_forward", "library PCR forward primer", MaterialKind.OLIGO,
                          inferred=not bool(getattr(module, "SEQ_PRIMERS", ())), source=source)
    reverse = graph.input("pcr_reverse", "library PCR reverse primer", MaterialKind.OLIGO,
                          inferred=not bool(getattr(module, "SEQ_PRIMERS", ())), source=source)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    final_construct = module.final_library()
    library = graph.product("library", final_construct.name, MaterialKind.LIBRARY,
                            state=final_construct)
    graph.react("library_pcr", "Amplify the final contact library", StepKind.PCR,
                ((current, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final_construct, reason="The pre-PCR contact population is not represented as one exact Construct."))
    primers = tuple(getattr(module, "SEQ_PRIMERS", ()))
    if primers:
        graph.sequence(library, SequencingPlan(
            library=final_construct, primers=primers,
            acquired_roles=tuple(getattr(module, "REQUIRED_ROLES", ())) or tuple(p.role for p in primers),
            read_lengths=getattr(module, "READ_LENGTHS", {}) or {}))
    else:
        graph.sequence(library, SequencingPlan(
            acquired_roles=tuple(getattr(module, "REQUIRED_ROLES", ("Read 1", "Read 2"))),
            read_lengths=getattr(module, "READ_LENGTHS", {}) or {},
            limitation=module.SEQUENCING_UNAVAILABLE))
    reviewed.extend(("pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr", "library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
