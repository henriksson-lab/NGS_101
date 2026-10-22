"""Source-owned MoClo module-to-transcription-unit and multigene assembly."""
from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, complement_segments, revcomp
from cohesive_assembly import CircularCohesiveAssembly, FivePrimeFragment
from protocol_graph import (Disposition, InputRole as R, MaterialKind as K,
                            MolecularRelation, ProtocolGraph, StepKind as S)
from protocol_model import closed_circle, role
from protocol_views import WorkflowSection, assemble_page
from page import head, info
import seqprimers as sp

PAGE_FILE = "moclo.html"
TITLE = "MoClo — modular DNA assembly"
SOURCE = "Weber et al., PLoS ONE (2011), doi:10.1371/journal.pone.0016765; Figures 1, 3 and 5; Results and MoClo cloning protocol"

# Figure 1B–C. Names identify source-defined part types, not biological inserts.
PARTS = {
    "promoter": ("GGAG", "TACT"),
    "UTR": ("TACT", "AATG"),
    "signal peptide": ("AATG", "AGGT"),
    "secreted CDS": ("AGGT", "GCTT"),
    "cytosolic CDS": ("AATG", "GCTT"),
    "terminator": ("GCTT", "CGCT"),
}
# Original Figure 5, forward positions 1 and 2 only. Other positions are not inferred.
POSITIONS = {1: ("TGCC", "GCAA"), 2: ("GCAA", "ACTA")}
BSAI_SITE = "GGTCTC"  # Figure 3
BPII_SITE = "GAAGAC"  # Figure 3


def bp(name, bases):
    return Segment(name, bases)


def part_module(name: str, kind: str) -> Construct:
    left, right = PARTS[kind]
    return Construct((bp("BsaI forward", BSAI_SITE),
                      Segment("left source spacer", "N", placeholder=True),
                      bp("left fusion", left), role(name), bp("right fusion", right),
                      Segment("right source spacer", "N", placeholder=True),
                      bp("BsaI reverse", revcomp(BSAI_SITE)), role("level 0 donor backbone")), name)


def digest_remainder(backbones: str, dropout: str) -> Scene:
    """Two independent DNA classes in the mixture, with no covalent fusion."""
    sc = Scene()
    for identifier, label in (("backbones", backbones), ("dropout", dropout)):
        top = (role(label),)
        bottom = complement_segments(top)
        sc.strand(identifier, top, label=label)
        sc.anneal(identifier + " bottom", bottom, to=identifier,
                  pair=(bottom[0].name, top[0].name), unpaired=(bottom[0].name,))
    return sc


def product_view(product: Construct, geometry: CircularCohesiveAssembly):
    """Draw the sealed circle from the same product used by the checker."""
    # Opening the circle for a second view allows every internal nick boundary to
    # be marked on both strands; the last-to-first seam is the circle closure.
    aliases = {}
    for identifier, fragment in geometry.fragments:
        aliases[identifier + " fusion"] = fragment.name + " fusion"
        for seg in fragment.core:
            aliases[identifier + ": " + seg.name] = fragment.name + ": " + seg.name
    view = Construct(tuple(replace(seg, name=aliases[seg.name]) for seg in product), product.name)
    sc = Scene.duplex(tuple(view), label="circle opened at drawing origin",
                      unpaired=tuple(seg.name + "'" for seg in view if seg.is_role_token()))
    for (identifier, fragment), (next_id, _) in zip(
            geometry.fragments[:-1], geometry.fragments[1:]):
        last = identifier + ": " + fragment.core[-1].name
        fusion = next_id + " fusion"
        sc.junction("top", aliases[last], aliases[fusion])
        first = next_id + ": " + dict(geometry.fragments)[next_id].core[0].name
        sc.junction("bottom", aliases[fusion] + "'", aliases[first] + "'")
    sc.footer("Last-to-first seam is closed in the companion circle view.", "top")
    circle = closed_circle(product.name, tuple(view),
                           closure_label="sealed last-to-first cohesive junction")
    from protocol_model import Drawing
    return Drawing(tuple(circle.rows()) + tuple(sc.rows()))


def protocol_graph() -> ProtocolGraph:
    g = ProtocolGraph(TITLE, authority=SOURCE, required_output=Disposition.CONSTRUCT)
    bsai = g.input("bsai", "BsaI", K.CATALYST)
    bpii = g.input("bpii", "BpiI", K.CATALYST)
    ligase = g.input("ligase", "T4 DNA ligase", K.CATALYST)
    blue_white = g.input("blue_white", "ampicillin and blue/white clone screening", K.SELECTION)
    red_white = g.input("red_white", "kanamycin and red/white clone screening", K.SELECTION)
    units = []

    for position, kinds in ((1, ("promoter", "UTR", "cytosolic CDS", "terminator")),
                            (2, ("promoter", "UTR", "signal peptide", "secreted CDS", "terminator"))):
        prefix = f"tu{position}"
        parts = []
        for n, kind in enumerate(kinds):
            identifier = f"{prefix}_part{n}"
            label = f"TU{position} {kind}"
            source = g.input(identifier + "_plasmid", label + " level 0 plasmid", K.MOLECULAR,
                             state=part_module(label, kind), primary=True)
            left, right = PARTS[kind]
            fragment = FivePrimeFragment(label, left, (role(label),), right)
            parts.append((source, identifier, fragment))
        left_pos, right_pos = POSITIONS[position]
        # Figure 3 places BpiI fusion sites outside the BsaI transcription-unit
        # cassette. The variable backbone is a role, not an asserted full map.
        vector_core = (bp("right BpiI fusion", right_pos),
                       Segment("right BpiI spacer", "NN", placeholder=True),
                       bp("right BpiI site", revcomp(BPII_SITE)), role(f"pL1F-{position} backbone"),
                       bp("left BpiI site", BPII_SITE),
                       Segment("left BpiI spacer", "NN", placeholder=True),
                       bp("left BpiI fusion", left_pos))
        vector_fragment = FivePrimeFragment(f"pL1F-{position} retained vector", "CGCT", vector_core, "GGAG")
        vector = g.input(prefix + "_vector", f"pL1F-{position} destination plasmid", K.MOLECULAR,
                         state=Construct((*vector_core, bp("left TU boundary", "GGAG"),
                                          Segment("left BsaI spacer", "N", placeholder=True),
                                          bp("left outward BsaI", revcomp(BSAI_SITE)), role("lacZ dropout"),
                                          bp("right outward BsaI", BSAI_SITE),
                                          Segment("right BsaI spacer", "N", placeholder=True),
                                          bp("right TU boundary", "CGCT")), f"pL1F-{position}"))
        parts.append((vector, prefix + "_vector_fragment", vector_fragment))
        cut_products = [g.product(identifier, fragment.name, K.MOLECULAR,
                                  state=fragment, chemical_state=fragment)
                        for _, identifier, fragment in parts]
        remainder = g.product(prefix + "_remainder", "donor backbones and lacZ dropout", K.MOLECULAR,
                              state=digest_remainder("donor-backbone fragment population", "lacZ dropout"))
        g.react(prefix + "_cut", "BsaI cleavage within the one-pot assembly", S.DIGESTION,
                      [(source, R.TEMPLATE) for source, _, _ in parts] + [(bsai, R.CATALYST)],
                      [*cut_products, remainder])
        geometry = CircularCohesiveAssembly(tuple((identifier, fragment) for _, identifier, fragment in parts),
                    "Figures 1 and 3: source-defined BsaI fusion sites join selected level 0 modules to the retained level 1 vector.")
        product = geometry.product(f"TU{position} level 1 plasmid")
        assembled = g.product(prefix + "_assembled", product.name, K.MOLECULAR,
                              state=product_view(product, geometry), chemical_state=product)
        g.react(prefix + "_ligate", "Seal the declared level 1 assembly", S.LIGATION,
                       [(cut_products[0], R.TEMPLATE)] + [(p, R.STRUCTURAL) for p in cut_products[1:]]
                       + [(ligase, R.CATALYST)], [assembled], geometry=geometry,
                       relations=(MolecularRelation.COHESIVE_ASSEMBLY,))
        selected = g.product(prefix + "_selected", f"TU{position} selected candidate plasmid", K.MOLECULAR,
                             state=assembled.state, chemical_state=product)
        rejected = g.product(prefix + "_rejected", "material not carried into selected clones", K.WASTE)
        g.react(prefix + "_select", "Transform and select level 1 candidate clones", S.SELECTION,
                         [(assembled, R.TEMPLATE), (remainder, R.TEMPLATE), (blue_white, R.SELECTION)],
                         [selected, rejected], changes_molecule=False,
                         relations=(MolecularRelation.STRAND_ORDER,))
        g.discard(rejected, "Unselected material is not carried into the chosen clone preparation; no genotype or recovery proportion is inferred from colony colour.")
        # The released TU core is derived from the chosen part fragments, not a
        # second transcription of the assembled molecule. Vector sequence stays out.
        transferred = []
        for _, identifier, fragment in parts[:-1]:
            transferred.append(bp(identifier + " fusion", fragment.left))
            transferred.extend(fragment.core)
        transferred.append(bp(prefix + " terminal fusion", "CGCT"))
        unit_fragment = FivePrimeFragment(f"TU{position}", left_pos, tuple(transferred), right_pos)
        unit = g.product(prefix + "_released", f"BpiI-released TU{position}", K.MOLECULAR,
                         state=unit_fragment, chemical_state=unit_fragment)
        backbone = g.product(prefix + "_backbone", "level 1 donor backbone", K.MOLECULAR,
                             state=Construct((role(f"pL1F-{position} donor backbone"),)))
        g.react(prefix + "_release", "BpiI cleavage exposes the position-specific TU ends", S.DIGESTION,
                          [(selected, R.TEMPLATE), (bpii, R.CATALYST)], [unit, backbone])
        units.append((unit, unit_fragment, backbone))

    # Figure 5: pELE-2 joins the position-2 end ACTA to the pL2-1 end GGGA.
    end_linker = FivePrimeFragment("pELE-2 end-linker", "ACTA", (role("basic end-linker 2"),), "GGGA")
    vector_fragment = FivePrimeFragment("pL2-1 retained vector", "GGGA", (role("pL2-1 backbone"),), "TGCC")
    linker = g.input("end_linker", "pELE-2 plasmid", K.MOLECULAR,
                     state=Construct((bp("left BpiI", BPII_SITE),
                                      Segment("left BpiI spacer", "NN", placeholder=True),
                                      bp("left fusion", end_linker.left), *end_linker.core,
                                      bp("right fusion", end_linker.right),
                                      Segment("right BpiI spacer", "NN", placeholder=True),
                                      bp("right BpiI", revcomp(BPII_SITE)),
                                      role("end-linker donor backbone"))))
    vector = g.input("level2_vector", "pL2-1 destination plasmid", K.MOLECULAR,
                     state=Construct((bp("right vector fusion", vector_fragment.left),
                                      *vector_fragment.core,
                                      bp("left vector fusion", vector_fragment.right),
                                      Segment("left BpiI spacer", "NN", placeholder=True),
                                      bp("left outward BpiI", revcomp(BPII_SITE)), role("CRed dropout"),
                                      bp("right outward BpiI", BPII_SITE),
                                      Segment("right BpiI spacer", "NN", placeholder=True))))
    linker_cut = g.product("linker_cut", end_linker.name, K.MOLECULAR, state=end_linker, chemical_state=end_linker)
    vector_cut = g.product("vector_cut", vector_fragment.name, K.MOLECULAR,
                           state=vector_fragment, chemical_state=vector_fragment)
    remainder = g.product("level2_remainder", "end-linker backbone and CRed dropout", K.MOLECULAR,
                          state=digest_remainder("end-linker donor backbone", "CRed dropout"))
    g.react("level2_cut", "BpiI cleavage of end-linker and destination plasmids", S.DIGESTION,
            [(linker, R.TEMPLATE), (vector, R.TEMPLATE), (bpii, R.CATALYST)],
            [linker_cut, vector_cut, remainder])
    geometry = CircularCohesiveAssembly(tuple((unit.id, fragment) for unit, fragment, _ in units)
                                       + ((linker_cut.id, end_linker), (vector_cut.id, vector_fragment)),
                  "Figures 3 and 5: forward positions 1 and 2 join via GCAA; pELE-2 closes ACTA to GGGA in pL2-1.")
    product = geometry.product("two-transcription-unit level 2 plasmid")
    joined = g.product("level2_joined", product.name, K.MOLECULAR,
                       state=product_view(product, geometry), chemical_state=product)
    g.react("level2_ligate", "Seal the declared level 2 assembly", S.LIGATION,
            [(units[0][0], R.TEMPLATE), (units[1][0], R.STRUCTURAL),
             (linker_cut, R.STRUCTURAL), (vector_cut, R.STRUCTURAL), (ligase, R.CATALYST)],
            [joined], geometry=geometry, relations=(MolecularRelation.COHESIVE_ASSEMBLY,))
    final = g.product("selected_construct", "selected level 2 candidate construct", K.MOLECULAR,
                      state=joined.state, chemical_state=product)
    rejected = g.product("level2_rejected", "material not carried into selected level 2 clones", K.WASTE)
    g.react("level2_select", "Transform and select level 2 candidate clones", S.SELECTION,
            [(joined, R.TEMPLATE), (remainder, R.TEMPLATE)]
            + [(backbone, R.TEMPLATE) for _, _, backbone in units] + [(red_white, R.SELECTION)],
            [final, rejected], changes_molecule=False, relations=(MolecularRelation.STRAND_ORDER,))
    g.discard(rejected, "Material not carried into chosen clones; colony selection does not establish the complete expected genotype.")
    g.constructed(final, "Expected assembled DNA among selected candidates; no sequencing acquisition or experimental coverage is asserted.")
    # Source coverage is an independently written inventory, not generated from
    # graph construction. Figures 1/3/5 and Results supply this bounded scope.
    g.declare_reviewed_scope((
        "bsai", "bpii", "ligase", "blue_white", "red_white",
        "tu1_part0_plasmid", "tu1_part1_plasmid", "tu1_part2_plasmid", "tu1_part3_plasmid",
        "tu1_vector", "tu1_cut", "tu1_ligate", "tu1_select", "tu1_release",
        "tu2_part0_plasmid", "tu2_part1_plasmid", "tu2_part2_plasmid", "tu2_part3_plasmid", "tu2_part4_plasmid",
        "tu2_vector", "tu2_cut", "tu2_ligate", "tu2_select", "tu2_release",
        "end_linker", "level2_vector", "level2_cut", "level2_ligate", "level2_select", "selected_construct",
    ), SOURCE)
    return g.require_publishable()


def render_page() -> str:
    g = protocol_graph()
    return assemble_page(g, (
        head(TITLE), '<div class="wrap">',
        f"<h1>{TITLE}</h1>", '<p class="research-notes"><a href="01_moclo.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1371/journal.pone.0016765">Weber et al., PLoS ONE (2011)</a>'),
        info("From sequenced level 0 modules to a level 2 construct. Part bodies and vector backbones are role placeholders, without asserted sequences or lengths. The two TU choices illustrate the published cytosolic and secreted part architectures; they do not name experimental clones. Restriction cleavage and sealing are drawn separately within the one-pot reactions. Checks establish the supplied cohesive-end and product relations, not complete vector digestion, clone genotype or library coverage."),
        WorkflowSection("Molecular assembly", paths=(
            ("tu1_released", ()), ("tu2_released", ()),
            ("linker_cut", ()), ("selected_construct", ("tu1_release",)),
        )),
        "<h2>Readout</h2>", info("The output is assembled DNA. The paper uses clone examination and expression demonstrations; this construction scope assigns no NGS library or run."),
        sp.read_spec_table([sp.read_spec("Sequencing output", None,
            "Not applicable to this construct-assembly scope", "not_applicable")]),
        "</div>",
    ))


def selftest() -> None:
    from checks import Check
    from chemistry_validation import Status, check_graph
    check = Check()
    graph = protocol_graph()
    check("source-scoped assembly graph has no contradictory molecular claims",
          [c for c in check_graph(graph) if c.status is Status.CONTRADICTORY], [])
    check.report()
