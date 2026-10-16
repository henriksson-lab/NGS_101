"""Molecular construct model for Drop-ChIP (Rotem et al., 2015)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

# Representative row 1 from Supplementary Table 2.  Every other barcode adaptor has
# the same core and replaces the eight-base barcode and its reversed copy.
REPRESENTATIVE_ADAPTOR = (
    "TTAAGGGCTTTCGTATCCGGGGGACCTTAATTAAGGTGGGGGGGATACCTTTCGGGTTAA"
)
REPRESENTATIVE_BARCODE = "GGGCTTTC"
ADAPTOR_CORE = "GTATCCGGGGGACCTTAATTAAGGTGGGGGGGATAC"
END = "TTAA"

SC_PCR1 = "TAAGGTGGGGGGGATAC"
SC_PCR2 = "TAAGGTCCCCCGGATAC"
PACI = "TTAATTAA"
BCIVI = "GTATCC"


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def barcode_adaptor(barcode: str = REPRESENTATIVE_BARCODE) -> Construct:
    """One published 60-nt adaptor strand, with symmetry enforced by construction."""
    if len(barcode) != 8 or any(b not in "ACGT" for b in barcode):
        raise ValueError("Drop-ChIP barcode must be eight DNA bases")
    return Construct([
        _seg("left blunt end", END),
        _seg("barcode", barcode, "cbc",
             feature=feature("cell_barcode", "cell_barcode", "whitelist", whitelist="published Drop-ChIP adapter set")),
        _seg("symmetric adaptor core", ADAPTOR_CORE),
        _seg("reversed barcode", barcode[::-1], "cbc",
             feature=feature("cell_barcode", "cell_barcode", "whitelist", whitelist="published Drop-ChIP adapter set")),
        _seg("right blunt end", END),
    ], name="Drop-ChIP barcode adaptor")


def adaptor_scene(barcode: str = REPRESENTATIVE_BARCODE) -> Scene:
    return Scene.duplex(list(barcode_adaptor(barcode)))


def blunt_ligation_geometry() -> AdapterLigationGeometry:
    sc = Scene()
    sc.strand("insert", [_seg("repaired blunt nucleosomal end", "X" * 18,
                               placeholder=True)], label="end-repaired nucleosomal DNA")
    sc.ligation_end("insert", "repaired blunt nucleosomal end", kind="blunt",
                    owner="substrate", text="blunt repaired DNA end")
    sc.strand("adapter", list(barcode_adaptor()), label="published barcode adapter")
    sc.ligation_end("adapter", "left blunt end", kind="blunt", owner="adapter",
                    text="published TTAA blunt adapter end")
    return AdapterLigationGeometry.known(
        "blunt", sc,
        statement="The defining source prints the double-stranded barcode adapter and specifies blunt ligation to end-repaired nucleosomal DNA.")


def blunt_ligated_fragment(barcode: str = REPRESENTATIVE_BARCODE) -> Construct:
    """Representative one of four ligation orientations: same top-strand orientation.

    Both adaptor/insert junctions are TTAA-to-genome blunt junctions.  The duplex itself
    is derived by Scene at render time rather than independently entered.
    """
    left = barcode_adaptor(barcode)
    right = barcode_adaptor(barcode)
    return Construct([
        *[_seg("left " + s.name, s.top, s.tag, placeholder=s.placeholder, feature=s.feature) for s in left],
        _seg("nucleosomal DNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        *[_seg("right " + s.name, s.top, s.tag, placeholder=s.placeholder, feature=s.feature) for s in right],
    ], name="blunt-ligated Drop-ChIP fragment")


def paci_product(barcode: str = REPRESENTATIVE_BARCODE) -> Construct:
    """Adaptor-labelled fragment after PacI removes the outer adaptor halves."""
    full = barcode_adaptor(barcode).top()
    site = full.index(PACI)
    # PacI cuts TTAAT^TAA.  The left molecule retains sequence after that top-strand cut;
    # at the right junction it retains the reciprocal prefix.
    left = full[site + 5:]
    right = full[:site + 3]
    if not left.startswith(SC_PCR1) or revcomp(SC_PCR2) not in right:
        raise ValueError("PacI product no longer exposes both published SC-PCR sites")
    if left != SC_PCR1 + barcode[::-1] + END \
            or right != END + barcode + revcomp(SC_PCR2):
        raise ValueError("PacI product no longer preserves the symmetric barcode junctions")
    return Construct([
        _seg("SC-PCR1 site", SC_PCR1),
        _seg("left barcode", barcode[::-1], "cbc",
             feature=feature("cell_barcode", "cell_barcode", "whitelist", whitelist="published Drop-ChIP adapter set")),
        _seg("left PacI half", END),
        _seg("nucleosomal DNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        _seg("right PacI half", END),
        _seg("right barcode", barcode, "cbc",
             feature=feature("cell_barcode", "cell_barcode", "whitelist", whitelist="published Drop-ChIP adapter set")),
        _seg("SC-PCR2 site'", revcomp(SC_PCR2)),
    ], name="PacI-trimmed Drop-ChIP fragment")


def sc_pcr_scene(barcode: str = REPRESENTATIVE_BARCODE) -> Scene:
    """PacI product with both SC-PCR primers placed on their derived sites."""
    top = list(paci_product(barcode))
    bottom = complement_segments(top)
    sc = Scene()
    sc.strand("top", top, label="template")
    sc.anneal("bottom", bottom, to="top", pair=(top[0].name + "'", top[0].name),
              label="template")
    sc.anneal("SC-PCR1", [_seg("SC-PCR1 primer", SC_PCR1)], to="bottom",
              pair=("SC-PCR1 primer", "SC-PCR1 site'"), label="SC-PCR1")
    sc.anneal("SC-PCR2", [_seg("SC-PCR2 primer", SC_PCR2)], to="top",
              pair=("SC-PCR2 primer", "SC-PCR2 site'"), label="SC-PCR2", above=True)
    sc.arrow("SC-PCR1", "")
    sc.arrow("SC-PCR2", "")
    sc.stack("SC-PCR2", "top", "bottom", "SC-PCR1")
    return sc


def bcivi_product(barcode: str = REPRESENTATIVE_BARCODE) -> tuple[Scene, Construct]:
    """BciVI-cut insert, with a single 3'-A overhang at each end.

    The exact retained sequence is sliced from the symmetric adaptor constructor.  The
    right terminal A is explicit; the reciprocal bottom-strand A is derived as its own
    unpaired segment in the Scene.
    """
    full = barcode_adaptor(barcode).top()
    left_site = full.index(BCIVI)
    right_site = full.index(revcomp(BCIVI))
    left = full[right_site - 5:]       # GGGGGGGATAC · reversed barcode · TTAA
    right = full[:left_site + 12]      # TTAA · barcode · GTATCCGGGGGA
    if not right.endswith("A"):
        raise ValueError("BciVI right product must end in the published single-A overhang")
    paired = [
        _seg("left 11-nt constant", left[:11]),
        _seg("left barcode", left[11:19], "cbc",
             feature=feature("cell_barcode", "cell_barcode", "whitelist", whitelist="published Drop-ChIP adapter set")),
        _seg("left PacI half", left[19:]),
        _seg("nucleosomal DNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        _seg("right PacI half", right[:4]),
        _seg("right barcode", right[4:12], "cbc",
             feature=feature("cell_barcode", "cell_barcode", "whitelist", whitelist="published Drop-ChIP adapter set")),
        _seg("right 11-nt constant", right[12:-1]),
    ]
    product = Construct([*paired, _seg("right 3' A", right[-1])],
                        name="BciVI-cut Drop-ChIP insert")
    sc = Scene()
    sc.strand("top", list(product), label="BciVI product")
    sc.anneal("bottom", [*complement_segments(paired), _seg("left 3' A", "A")],
              to="top", pair=(paired[0].name + "'", paired[0].name),
              label="BciVI product", unpaired=("left 3' A",))
    sc.mark("top", "right 3' A", "3'-A")
    sc.mark("bottom", "left 3' A", "3'-A")
    return sc, product


def inferred_final_library(barcode: str = REPRESENTATIVE_BARCODE) -> Construct:
    """Only the supported inner architecture; outer Illumina identity is unresolved."""
    _, inner = bcivi_product(barcode)
    return Construct([
        _seg("left unresolved Illumina arm", "[ILLUMINA ARM]", inferred=True,
             placeholder=True),
        *list(inner),
        _seg("right unresolved Illumina arm", "[ILLUMINA ARM]", inferred=True,
             placeholder=True),
    ], name="Drop-ChIP final library, authoritative boundary")


def _validate_source_transcription() -> None:
    if barcode_adaptor().top() != REPRESENTATIVE_ADAPTOR:
        raise ValueError("representative adaptor no longer matches Supplementary Table 2")
    if barcode_adaptor().top().count(PACI) != 1:
        raise ValueError("barcode adaptor must contain one central PacI site")
    if barcode_adaptor().top().count(BCIVI) != 1 \
            or barcode_adaptor().top().count(revcomp(BCIVI)) != 1:
        raise ValueError("barcode adaptor must expose reciprocal BciVI sites")
    adaptor_scene().rows()
    sc_pcr_scene().rows()
    bcivi_product()[0].rows()


_validate_source_transcription()


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1038/nbt.3383"
    inferred_source = source + "; downstream Illumina adapter and library primers unresolved from checked article and supplements"
    graph = ProtocolGraph("Drop-ChIP", authority=source)
    dna = graph.input("nucleosomal_dna", "nucleosomal DNA fragments",
                      MaterialKind.BIOLOGICAL, primary=True)
    repair = graph.input("repair_enzymes", "DNA end-repair enzymes", MaterialKind.CATALYST)
    repaired = graph.product("repaired_dna", "blunt end-repaired nucleosomal DNA",
                             MaterialKind.MOLECULAR)
    graph.react("end_repair", "End-repair nucleosomal DNA", StepKind.END_REPAIR,
                ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    adapter = graph.input("barcode_adapter", "published symmetric barcode adapter",
                          MaterialKind.ADAPTER, state=adaptor_scene())
    ligase = graph.input("ligase", "DNA ligase", MaterialKind.CATALYST)
    ligated_construct = blunt_ligated_fragment()
    ligated = graph.product("barcoded_fragments", ligated_construct.name,
                            MaterialKind.MOLECULAR, state=ligated_construct)
    graph.react("barcode_ligation", "Blunt-ligate barcode adapters to both fragment ends",
                StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), geometry=blunt_ligation_geometry())
    paci = graph.input("paci", "PacI restriction enzyme", MaterialKind.CATALYST)
    trimmed_construct = paci_product()
    trimmed = graph.product("paci_trimmed", trimmed_construct.name, MaterialKind.MOLECULAR,
                            state=trimmed_construct)
    graph.react("paci_digest", "Remove concatemers and outer adapter halves with PacI",
                StepKind.DIGESTION,
                ((ligated, InputRole.TEMPLATE), (paci, InputRole.CATALYST)), (trimmed,))
    sc1 = graph.input("sc_pcr1", "published SC-PCR1 primer", MaterialKind.OLIGO)
    sc2 = graph.input("sc_pcr2", "published SC-PCR2 primer", MaterialKind.OLIGO)
    polymerase = graph.input("polymerase", "SC-PCR and library PCR polymerase",
                             MaterialKind.CATALYST)
    amplified = graph.product("sc_amplicon", "SC-PCR-amplified barcode-labelled fragment",
                              MaterialKind.MOLECULAR, state=sc_pcr_scene())
    graph.react("sc_pcr", "Amplify the PacI-trimmed barcode-labelled fragment", StepKind.PCR,
                ((trimmed, InputRole.TEMPLATE), (sc1, InputRole.STRUCTURAL),
                 (sc2, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.exact(
                    trimmed_construct, PCRPrimer("SC-PCR1", SC_PCR1, 17),
                    PCRPrimer("SC-PCR2", SC_PCR2, 17),
                    statement="The defining source prints both 17-nt SC-PCR primers and the adapter sequence from which their inward landing sites are derived."))
    bcivi = graph.input("bcivi", "BciVI type-IIS restriction enzyme", MaterialKind.CATALYST)
    bcivi_scene, bcivi_construct = bcivi_product()
    da_insert = graph.product("da_insert", bcivi_construct.name, MaterialKind.MOLECULAR,
                              state=bcivi_scene)
    graph.react("bcivi_digest", "Remove SC-PCR ends and expose 3′-A overhangs with BciVI",
                StepKind.DIGESTION,
                ((amplified, InputRole.TEMPLATE), (bcivi, InputRole.CATALYST)), (da_insert,))
    illumina_adapter = graph.input("illumina_adapter", "unresolved Illumina adapter",
                                   MaterialKind.ADAPTER, inferred=True,
                                   source=inferred_source)
    precursor = graph.product("adapter_ligated", "Illumina-adapter-ligated Drop-ChIP insert",
                              MaterialKind.MOLECULAR, inferred=True,
                              source=inferred_source)
    graph.react("illumina_ligation", "Ligate Illumina adapters to the 3′-A-ended insert",
                StepKind.LIGATION,
                ((da_insert, InputRole.TEMPLATE), (illumina_adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (precursor,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The source establishes adapter ligation; selected Illumina adapter bases and end geometry remain unresolved from the checked article and supplements."),
                source=inferred_source)
    library_primers = graph.input("library_primers", "unresolved library PCR-primer pair",
                                  MaterialKind.OLIGO, inferred=True, source=inferred_source)
    final = inferred_final_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final,
                            inferred=True, source=inferred_source)
    graph.react("library_pcr", "PCR-complete the Illumina library", StepKind.PCR,
                ((precursor, InputRole.TEMPLATE), (library_primers, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The selected library PCR-primer bases remain unresolved from the checked article and supplements.",
                    symmetric_species=True), source=inferred_source)
    graph.sequence(library, SequencingPlan(
        limitation="The paper publishes barcode coordinates but not the Illumina adapter or sequencing-primer sequences.",
        acquired_roles=("Read 1", "Read 2", "Index 1 (i7)"),
        read_uses={"Read 1": "extract_then_include", "Read 2": "extract_then_include",
                   "Index 1 (i7)": "demultiplex"}))
    graph.declare_reviewed_scope((
        "nucleosomal_dna", "repair_enzymes", "end_repair", "barcode_adapter",
        "ligase", "barcode_ligation", "paci", "paci_digest", "sc_pcr1",
        "sc_pcr2", "polymerase", "sc_pcr", "bcivi", "bcivi_digest",
        "illumina_adapter", "illumina_ligation", "library_primers", "library_pcr",
        "library"), source)
    return graph.require_publishable()
