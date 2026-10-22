"""Published Yeast Toolkit: direct primary-type cassette assembly."""
from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from cohesive_assembly import CircularCohesiveAssembly, FivePrimeFragment
from protocol_graph import (Disposition, InputRole as R, MaterialKind as K,
                            MolecularRelation, ProtocolGraph, StepKind as S)
from protocol_model import closed_circle, role
from protocol_views import WorkflowSection, assemble_page
from page import head, info
import seqprimers as sp

PAGE_FILE = "yeast_toolkit.html"
TITLE = "Yeast Toolkit — direct cassette assembly"
SOURCE = "Lee et al. (2015), doi:10.1021/sb500366v; published SI Figures S6–S7, S15–S17 and cassette assembly pp20–22; published YTK_Parts.xls"

# Primary types only. The selected identities are workbook records, not a claim
# that the source tested this particular combination as one experimental clone.
PARTS = (
    ("pYTK002", "ConLS", "CCCT", "AACG"),
    ("pYTK009", "pTDH3", "AACG", "TATG"),
    ("pYTK032", "mTurquoise2", "TATG", "ATCC"),
    ("pYTK051", "tENO1", "ATCC", "GCTG"),
    ("pYTK067", "ConR1", "GCTG", "TACA"),
    ("pYTK074", "URA3", "TACA", "GAGT"),
    ("pYTK081", "CEN6/ARS4", "GAGT", "CCGA"),
    ("pYTK083", "AmpR–ColE1", "CCGA", "CCCT"),
)
BSAI_SITE = "GGTCTC"
BSMBI_SITE = "CGTCTC"


def bp(name, bases):
    return Segment(name, bases)


def connector_barcode(identifier, label, length):
    return role(label, length_bp=length, molecular_feature=feature(
        identifier, "inline_barcode", "fixed",
        note="Source-defined assembly connector identifier; not a UMI or sample-sheet index. Bases are not expanded in this scope."))


def part_core(plasmid, label):
    # SI Figure S7B: the 21-bp barcode is inside the cassette. The later
    # multigene scar length is not assigned from this drawing alone.
    if plasmid == "pYTK002":
        return (bp("EcoRI and XbaI region", "GAATTCGCATCTAGA"),
                connector_barcode("conls.concatenated", "ConLS concatenated barcodes", 143),
                bp("BsmBI forward", BSMBI_SITE), bp("BsmBI spacer", "A"),
                bp("ConLS higher-level fusion", "CTGA"),
                connector_barcode("conls.scar", "ConLS cassette barcode", 21))
    # SI Figure S15B, with the specific ConR1 fusion independently read from
    # pYTK067. This site points upstream; it is not a second ConLS geometry.
    if plasmid == "pYTK067":
        return (bp("ConR1 higher-level fusion", "CCAA"), bp("BsmBI spacer", "T"),
                bp("BsmBI reverse", revcomp(BSMBI_SITE)),
                connector_barcode("conr1.concatenated", "ConR1 concatenated barcodes", 143),
                bp("SpeI and PstI region", "ACTAGTGCACTGCAG"))
    return (role(label),)


def digest_remainders():
    """Independent donor backbones and RFP dropout, never a covalent join."""
    sc = Scene()
    for identifier, label in (("donors", "chloramphenicol donor-backbone population"),
                              ("dropout", "Type 8 mRFP1 dropout")):
        top = (role(label),)
        bottom = complement_segments(top)
        sc.strand(identifier, top, label=label)
        sc.anneal(identifier + " bottom", bottom, to=identifier,
                  pair=(bottom[0].name, top[0].name), unpaired=(bottom[0].name,))
    return sc


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE, required_output=Disposition.CONSTRUCT)
    enzyme = g.input("bsai", "BsaI", K.CATALYST)
    ligase = g.input("ligase", "DNA ligase (identity unassigned in this scope)", K.CATALYST)
    selection = g.input("selection", "ampicillin and red/white colony screening", K.SELECTION)
    sources, fragments, geometry_parts = [], [], []
    for n, (plasmid, label, left, right) in enumerate(PARTS, 1):
        core = part_core(plasmid, label)
        fragment = FivePrimeFragment(f"Type {n}: {label}", left, core, right)
        remainder = "mRFP1 dropout" if n == 8 else "chloramphenicol donor backbone"
        donor = Construct((bp("BsaI forward", BSAI_SITE), bp("left spacer", "A"),
                           bp("left fusion", left), *core, bp("right fusion", right),
                           bp("right spacer", "T"), bp("BsaI reverse", revcomp(BSAI_SITE)),
                           role(remainder)), plasmid)
        source = g.input(f"type{n}_plasmid", f"{plasmid}: {label}", K.MOLECULAR,
                         state=closed_circle(plasmid, tuple(donor)), chemical_state=donor,
                         primary=True)
        cut = g.product(f"type{n}_fragment", fragment.name, K.MOLECULAR,
                        state=fragment, chemical_state=fragment)
        sources.append(source)
        fragments.append(cut)
        geometry_parts.append((cut.id, fragment))
    remainders = g.product("remainders", "donor backbones and Type 8 dropout", K.MOLECULAR,
                           state=digest_remainders())
    g.react("release", "BsaI cleavage within the direct cassette assembly", S.DIGESTION,
            [(s, R.TEMPLATE) for s in sources] + [(enzyme, R.CATALYST)],
            [*fragments, remainders])
    geometry = CircularCohesiveAssembly(tuple(geometry_parts),
        "Published SI Figure S6: one primary part of each Type; the Type 8 AmpR–ColE1 fragment is the receiving backbone.")
    product = geometry.product("eight-part yeast cassette")
    # Labels for the circle are derived from the exact checked product; no bases
    # or barcode identities are reconstructed for rendering.
    aliases = {}
    for identifier, fragment in geometry_parts:
        aliases[identifier + " fusion"] = fragment.name + " fusion"
        for seg in fragment.core:
            aliases[identifier + ": " + seg.name] = fragment.name + ": " + seg.name
    displayed = tuple(replace(seg, name=aliases[seg.name]) for seg in product)
    joined = g.product("joined", product.name, K.MOLECULAR,
                       state=closed_circle(product.name, displayed), chemical_state=product)
    g.react("seal", "Seal the declared eight-part circle", S.LIGATION,
            [(fragments[0], R.TEMPLATE)] + [(p, R.STRUCTURAL) for p in fragments[1:]]
            + [(ligase, R.CATALYST)], [joined], geometry=geometry,
            relations=(MolecularRelation.COHESIVE_ASSEMBLY,))
    selected = g.product("selected_construct", "selected cassette candidate", K.MOLECULAR,
                         state=joined.state, chemical_state=product)
    rejected = g.product("rejected", "material not carried into selected colonies", K.WASTE)
    g.react("select", "Transform and select candidate cassette clones", S.SELECTION,
            [(joined, R.TEMPLATE), (remainders, R.TEMPLATE), (selection, R.SELECTION)],
            [selected, rejected], changes_molecule=False,
            relations=(MolecularRelation.STRAND_ORDER,))
    g.discard(rejected, "Unselected material is not carried forward; colony colour does not prove the expected genotype or assembly efficiency.")
    g.constructed(selected, "Expected direct cassette construct among selected candidates; expression, clone genotype and pooled-library representation are not inferred.")
    g.declare_reviewed_scope((
        "bsai", "ligase", "selection", "type1_plasmid", "type2_plasmid", "type3_plasmid",
        "type4_plasmid", "type5_plasmid", "type6_plasmid", "type7_plasmid", "type8_plasmid",
        "release", "seal", "select", "selected_construct",
    ), SOURCE)
    return g.require_publishable()


def render_page():
    return assemble_page(protocol_graph(), (
        head(TITLE), '<div class="wrap">', f"<h1>{TITLE}</h1>",
        '<p class="research-notes"><a href="01_yeast_toolkit.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1021/sb500366v">Lee et al., ACS Synthetic Biology (2015)</a>; '
             '<a href="https://doi.org/10.1021/sb500366v.s001">published assembly supplement</a> and '
             '<a href="https://doi.org/10.1021/sb500366v.s002">sequence archive</a>.'),
        info("Direct assembly from one primary part of each Type. The selected workbook parts illustrate the source architecture; no particular experimental clone is asserted. Part bodies and unexpanded barcodes retain their declared roles and lengths. Checks establish the supplied cohesive ends and circular product, with symbolic cores remaining conditional; they do not establish complete plasmid digestion, clone genotype, expression or pool coverage. Part entry, the retained-BsaI GFP cassette-vector branch and multigene BsmBI assembly remain outside this schematic."),
        info("The author-hosted article copy specifies T7 DNA ligase, but its exact publication-version identity is unverified. Ligase identity remains unassigned in this published-SI-scoped model pending that check; this is an open source question, not a claim of non-disclosure."),
        info("Published discrepancy: SI Figure S16B shows Type 6 boundaries TGGC…GCTG (the right end is printed on the complementary strand as CGAC), whereas Figure S16A, its text and overview Figure S6 give TACA…GAGT. This model uses TACA…GAGT, independently supported by BsaI cuts in the published pYTK074 sequence. Panel B is not used to assign the boundary; the authors' intended correction to that panel is unconfirmed."),
        WorkflowSection("Direct cassette assembly", paths=(
            ("type1_fragment", ()), ("remainders", ("release",)),
            ("selected_construct", ("release",)),
        )),
        "<h2>Readout</h2>",
        info("The output is constructed DNA. Connector barcodes identify assembly components and are distinct from sample-sheet indexes. No NGS preparation or run is assigned to this scope."),
        sp.read_spec_table([sp.read_spec("Sequencing output", None,
            "Not applicable to this construct-assembly scope", "not_applicable")]),
        "</div>",
    ))


def selftest():
    from checks import Check, Source, have
    from chemistry_validation import Status, check_graph
    from plasmid import BsaI, BsmBI, Plasmid, digest
    from storage import sources_root
    check = Check()
    graph = protocol_graph()
    check("direct cassette graph has no contradictory molecular claims",
          [c for c in check_graph(graph) if c.status is Status.CONTRADICTORY], [])
    source = Source(sources_root() / Path(__file__).resolve().parents[1].name
                    / "3625419_YTK_Parts.xls.txt",
                    "Published supplement doi:10.1021/sb500366v.s002; YTK_Parts.xls, converted with tools/doctext.py")
    if have(check, source, label="published primary-part fusion transcription"):
        records = {r["Plasmid Name"]: r for r in csv.DictReader(
            source.path.read_text().split("\n")[1:], delimiter="\t") if r.get("Plasmid Name")}
        for name, _, left, right in PARTS:
            plasmid = Plasmid(name, records[name]["Plasmid Sequence"].upper())
            check(name + " published BsaI cut ends", [c.overhang for c in digest(plasmid, BsaI)], [left, right])
        for name, expected in (("pYTK002", "CTGA"), ("pYTK067", "CCAA")):
            plasmid = Plasmid(name, records[name]["Plasmid Sequence"].upper())
            check(name + " connector BsmBI end", [c.overhang for c in digest(plasmid, BsmBI)], [expected])
            cuts = digest(plasmid, BsaI)
            actual_core = plasmid.sub(cuts[0].bottom_cut, cuts[1].top_cut)
            pattern = "".join(f"[ACGT]{{{seg.length_bp}}}" if seg.is_role_token()
                              else re.escape(seg.top) for seg in part_core(name, ""))
            check(name + " retained connector architecture matches source sequence",
                  re.fullmatch(pattern, actual_core) is not None)
    check.report()
