"""Source-scoped DropSynth2018 construction graph; detailed chemistry remains open."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Row, Scene, Segment, complement_segments, feature, revcomp
from page import head, info
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import Disposition, InputRole as R, MaterialKind as K, ProtocolGraph, StepKind as S, Terminal
from protocol_model import Drawing, closed_circle, role
from protocol_views import OligoEntry, OligoSection, WorkflowSection, assemble_page

PAGE_FILE = "dropsynth.html"
TITLE = "DropSynth — capture groups and pooled construction"
SOURCE = "Plesa et al.2018 doi:10.1126/science.aao5167; printed final Science article science.aao5167.pdf, assembly description and original Figure1A–D; the PMC XML is an associated author manuscript, not asserted final-equivalent. Published SI generic DropSynth protocol pp3–4 now checked; no DropSynth2.0 or superseded preprint recipe is substituted. Target-specific design and functional-assay branches are outside this graph."


# Published SI TableS6, rendered original p52 checked including wrapped anchor.
BEAD_LIGATION = "TCCGCGAGTAAACCTAACAA"
BEAD_ANCHOR = "TTGTTAGGTTTACTCGCGGAACACGTGCTATTAGATGCCT"
BEAD_COMMON = "AGGCATCTAATAGCACGTGT"
SOURCE_OLIGO_EXCEPTIONS = {
    "bead_anchor": "Published SI TableS6 rendered original p52 verified: the40nt anchor wraps after31 bases through the modifications/amount columns, which the generic text matcher cannot rejoin. selftest independently reconstructs the two source sequence tokens and checks the full anchor plus source end modifications.",
    "bead_barcodes": "Published SI TableS6 original p52 prints12-mer microbead barcode reverse complement + AGGCATCTAATAGCACGTGT across table cells, rather than a contiguous32-base oligo. selftest checks that exact shorthand and common-arm token. No actual barcode values or equal-N physical pairing are asserted.",
}


def bead_barcode(group):
    return (Segment("bead " + group.identifier + " complementary capture barcode", "N" * 12,
        placeholder=True, feature=group.barcode.feature,
        note="TableS6 specifies the reverse complement of the fragment capture barcode. Values remain unknown; this is a declared group association, not a sequence equality or binding proof."),
        Segment("bead " + group.identifier + " common arm", BEAD_COMMON))


def bead_geometry(groups, *, sealed=False, phosphorylated=False):
    if BEAD_ANCHOR != revcomp(BEAD_COMMON + BEAD_LIGATION):
        raise ValueError("source anchor must complement both common arms in the stated order")
    sc = Scene()
    for group in groups:
        name = group.identifier
        anchor_lig = Segment("anchor " + name + " ligation complement", revcomp(BEAD_LIGATION))
        anchor_bc = Segment("anchor " + name + " barcode-arm complement", revcomp(BEAD_COMMON))
        sc.strand(name + " anchor", (anchor_lig, anchor_bc), label="group " + name + " anchor", mod5="dual biotin")
        barcode = bead_barcode(group)
        lig = Segment("bead " + name + " ligation oligo", BEAD_LIGATION)
        if sealed:
            sc.anneal(name + " sealed", (*barcode, lig), to=name + " anchor",
                pair=(barcode[-1].name, anchor_bc.name), label="group " + name + " ligated capture strand",
                mod5="phosphate" if phosphorylated else "", mod3="biotin")
            sc.junction(name + " sealed", barcode[-1].name, lig.name, "ligation")
        else:
            sc.anneal(name + " barcode", barcode, to=name + " anchor",
                pair=(barcode[-1].name, anchor_bc.name), label="group " + name + " barcode oligo")
            sc.anneal(name + " ligation", (lig,), to=name + " anchor",
                pair=(lig.name, anchor_lig.name), label="group " + name + " ligation oligo",
                mod5="phosphate", mod3="biotin")
            sc.ligation_end(name + " ligation", lig.name, kind="splinted", polarity="5-prime", owner="adapter",
                text="source5-prime phosphate at splinted nick")
            sc.ligation_end(name + " barcode", barcode[-1].name, kind="splinted", polarity="3-prime", owner="substrate",
                text="source barcode3-prime end at splinted nick")
        sc.note(name + " anchor", "Common arms pair exactly; unknown12mer remains unpaired. Bead occupancy and physical capture of gene fragments are unassigned.")
    return sc


def capture_feature(group):
    return Segment("unplaced capture barcode " + group, "N" * 12, placeholder=True,
        feature=feature("dropsynth.capture." + group, "feature_barcode", "unknown", group=group,
            note="Logical assembly-group capture feature. Bases, strand and position on the oligo remain unassigned; equal unknown strings do not establish sequence identity or binding."))


@dataclass(frozen=True)
class Fragment:
    identifier: str
    capture: Segment

    def __post_init__(self):
        if not self.identifier.strip() or self.capture.feature is None or len(self.capture.top) != 12:
            raise ValueError("fragment class needs an identity and its structured twelve-base capture feature")


@dataclass(frozen=True)
class CaptureGroup:
    identifier: str
    barcode: Segment
    fragments: tuple[Fragment, ...]

    def __post_init__(self):
        if not self.identifier or not self.fragments:
            raise ValueError("capture group needs an identity and explicit fragment classes")
        identity = self.barcode.feature
        if identity is None or identity.group != self.identifier or len(self.barcode.top) != 12:
            raise ValueError("capture-group identity must match its structured feature")
        ids = [f.identifier for f in self.fragments]
        if len(set(ids)) != len(ids):
            raise ValueError("a fragment class cannot be counted twice in one capture group")
        if any(f.capture.feature != identity for f in self.fragments):
            raise ValueError("a declared capture group cannot silently accept another group's feature")


def schematic_groups():
    # A/B are schematic identities, not published barcode values or target genes.
    # Four-fragment model assemblies are shown in original Figure1D.
    return tuple(CaptureGroup(name, capture_feature(name), tuple(
        Fragment(name + ".fragment" + str(i), capture_feature(name)) for i in range(1, 5)))
        for name in ("a", "b"))


@dataclass(frozen=True)
class SeparateFragments:
    groups: tuple[CaptureGroup, ...]
    phase: str
    duplex: bool
    capture_present: bool = True
    exposed_capture: bool = False
    forward_biotin: bool = False

    def __post_init__(self):
        names = [g.identifier for g in self.groups]
        ids = [f.identifier for g in self.groups for f in g.fragments]
        if self.exposed_capture and (not self.duplex or not self.capture_present):
            raise ValueError("exposed capture state needs a duplex body and retained capture feature")
        if len(set(names)) != len(names) or len(set(ids)) != len(ids):
            raise ValueError("molecular classes must remain separate across capture groups")

    def scene(self):
        sc = Scene()
        for group in self.groups:
            for fragment in group.fragments:
                label = fragment.identifier + ": " + self.phase
                part = role(fragment.identifier + " molecular layout unassigned")
                parts = (part,)
                if self.exposed_capture:
                    overhang = Segment(fragment.identifier + " source5-prime capture overhang", group.barcode.top,
                        placeholder=True, feature=group.barcode.feature,
                        note="SI p3–4: twelve-base exposed top-strand capture region; unknown bases do not establish binding or identity across groups.")
                    parts = (overhang, part)
                sc.strand(fragment.identifier, parts, label=label, mod5="biotin" if self.forward_biotin else "")
                if self.duplex:
                    bottom = complement_segments((part,))
                    sc.anneal(fragment.identifier + " bottom", bottom, to=fragment.identifier,
                        pair=(bottom[0].name, part.name), unpaired=(bottom[0].name,))
        return sc

    def rows(self):
        rows = self.scene().rows()
        for group in self.groups:
            meaning = (("exposed5-prime capture feature " if self.exposed_capture else "contains unplaced capture feature ") + group.barcode.feature.id
                       if self.capture_present else "capture identity is provenance, not a retained DNA feature")
            rows.append(Row(chunks=[("Group " + group.identifier + ": " + meaning + "; binding geometry/occupancy unassigned", None, False)]))
        return rows


def bead_state(groups):
    # Only unplaced barcode-feature annotations are shown, not a fabricated bead
    # oligo sequence, attachment polarity or duplex with an unknown N-mer.
    return Drawing(tuple(Row(chunks=[(
        "Bead group " + g.identifier + ": intended complementary 12-nt capture region; sequence/attachment geometry unassigned",
        None, False)]) for g in groups))


def expected_assembly(group):
    # Unknown overlaps prevent concatenating fragment strings or counting their
    # lengths as a gene length. This is the intended assembly class only.
    return Construct((role("group " + group.identifier + " intended assembly; sequence and overlaps unassigned", inferred=True),),
                     "group " + group.identifier + " intended assembled DNA")


def clone_barcode():
    return Segment("later clone/assembly barcode", "N" * 20, placeholder=True, inferred=True,
        feature=feature("dropsynth.clone", "feature_barcode", "random", group="post-assembly clone tracking",
            note="Later random twenty-base assembly/clone barcode from the published main text. Independent of the twelve-base bead capture feature; no particular clone sequence or sample index is assigned."))


def expected_clone():
    # Order is an explicitly inferred drawing arrangement, not a vector map.
    return Construct((role("assembled insert; exact position unassigned", inferred=True), clone_barcode(),
        role("source expression-plasmid context; order/junctions unassigned", inferred=True)),
        "intended barcode-bearing expression-plasmid class")


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE, required_output=Disposition.CONSTRUCT,
        representation="illustrative_example", example_description="Two invented A/B group labels illustrate the published four-fragment assembly class. Barcode values, target genes, overlaps and recipient arrangement are unassigned; this is not a reconstruction of a particular experimental clone.")
    groups = schematic_groups()
    initial = SeparateFragments(groups, "synthesis pool", False)
    pool = g.input("pool", "source oligo pool; two schematic capture groups", K.OLIGO, state=initial, primary=True)
    forward = g.input("pool_f", "source pool-amplification forward reagent; identity unassigned", K.OLIGO)
    reverse = g.input("pool_r", "source pool-amplification reverse reagent; identity unassigned", K.OLIGO)
    polymerase = g.input("pool_polymerase", "KAPA HiFi for oligo-pool amplification", K.CATALYST)
    amplified_state = SeparateFragments(groups, "amplified oligo class", True)
    amplified = g.product("amplified", "amplified oligo classes remain separate molecules", K.MOLECULAR, state=amplified_state)
    g.react("amplify", "Amplify the array-derived oligo pool", S.PCR,
        [(pool, R.TEMPLATE), (forward, R.STRUCTURAL), (reverse, R.STRUCTURAL), (polymerase, R.CATALYST)], [amplified],
        geometry=PCRGeometry.unavailable(initial.scene(), reason="Published SI p4 and detailed protocol Step2 specify initial subpool amplification. Selected15mer primer values, initial strand conversion and full oligo states remain unassigned."))
    bulk_f = g.input("bulk_f", "5-prime biotinylated subpool forward primer; selected bases unassigned", K.OLIGO)
    bulk_r = g.input("bulk_r", "subpool reverse primer; selected bases unassigned", K.OLIGO)
    bulk_state = SeparateFragments(groups, "bulk-amplified oligo; biotin on source forward end", True, forward_biotin=True)
    bulk = g.product("bulk", "bulk-amplified biotinylated oligo classes", K.MOLECULAR, state=bulk_state)
    g.react("bulk_amplify", "Bulk amplify with the biotinylated forward primer", S.PCR,
        [(amplified, R.TEMPLATE), (bulk_f, R.STRUCTURAL), (bulk_r, R.STRUCTURAL), (polymerase, R.CATALYST)], [bulk],
        geometry=PCRGeometry.unavailable(amplified_state.scene(), reason="SI p4 and detailed emulsion protocol Step3 distinguish the second PCR and biotinylated forward end. Exact selected subpool primers and end coordinates remain unassigned; no qPCR-derived cycle count or yield is inferred."))
    nickase = g.input("nickase", "Nt.BspQI; exact selected-oligo cut coordinates unassigned", K.CATALYST)
    exposed_state = SeparateFragments(groups, "source5-prime capture region exposed; exact cut bases unassigned", True, exposed_capture=True)
    exposed = g.product("exposed", "oligo classes with exposed capture regions", K.MOLECULAR, state=exposed_state)
    short_state = Scene(); short_state.strand("short", (role("displaced short-strand class"),), mod5="biotin")
    short = g.product("short_strands", "short displaced strand classes; geometry unassigned", K.MOLECULAR,
                      state=short_state)
    g.react("expose", "Nick and expose the capture region", S.DIGESTION,
        [(bulk, R.TEMPLATE), (nickase, R.CATALYST)], [exposed, short])
    short_capture = g.input("short_capture", "streptavidin M270 Dynabeads for biotinylated short-strand removal", K.SELECTION)
    processed = g.product("processed", "processed oligo classes after short-strand removal", K.MOLECULAR, state=exposed_state)
    removed_short = g.product("removed_short", "captured short biotinylated strand class", K.WASTE)
    g.react("remove_short", "Remove the short biotinylated cleavage product", S.AFFINITY_CAPTURE,
        [(exposed, R.TEMPLATE), (short, R.TEMPLATE), (short_capture, R.SELECTION)], [processed, removed_short])
    g.discard(removed_short, "Published SI p4 specifies M270 capture followed by column cleanup before assembly-bead association. No selected target sequence is assigned.")
    anchor = g.input("bead_anchor", "TableS6 anchor;5-prime dual biotin", K.OLIGO,
        state=Construct((Segment("source anchor", BEAD_ANCHOR),), "bead anchor"))
    ligation_oligo = g.input("bead_ligation", "TableS6 ligation oligo;5-prime phosphate and3-prime biotin", K.OLIGO,
        state=Construct((Segment("source ligation oligo", BEAD_LIGATION),), "bead ligation oligo"))
    barcode_scene = Scene()
    for group in groups:
        barcode_scene.strand(group.identifier, bead_barcode(group), label="separate source barcode class " + group.identifier)
    barcodes = g.input("bead_barcodes", "TableS6 barcode oligos; two schematic group identities", K.OLIGO, state=barcode_scene)
    ready_scene = bead_geometry(groups)
    bead_ready = g.product("bead_ready", "independent source bead-oligo annealing classes", K.MOLECULAR, state=ready_scene)
    g.react("anneal_bead_oligos", "Anneal source anchor, barcode and ligation oligos", S.HYBRIDIZATION,
        [(anchor, R.TEMPLATE), (barcodes, R.STRUCTURAL), (ligation_oligo, R.STRUCTURAL)], [bead_ready])
    bead_ligase = g.input("bead_ligase", "T4 DNA ligase for bead preparation", K.CATALYST)
    ligated_beads = g.product("ligated_beads", "ligated capture strands on source anchors", K.MOLECULAR,
        state=bead_geometry(groups, sealed=True))
    g.react("ligate_bead_oligos", "Seal the common-arm nick", S.NICK_SEALING,
        [(bead_ready, R.TEMPLATE), (bead_ligase, R.CATALYST)], [ligated_beads],
        geometry=AdapterLigationGeometry.known("splinted", ready_scene,
            statement="SI TableS6 and bead protocol pp10–11: barcode common arm3-prime end meets5-prime-phosphorylated ligation oligo on the dual-biotin anchor. Both20nt common regions are verified by Scene. No unknown barcode is treated as a base-pair match."))
    pnk = g.input("bead_pnk", "T4 PNK", K.CATALYST)
    phosphorylated = g.product("phosphorylated_beads", "phosphorylated source capture-strand classes", K.MOLECULAR,
        state=bead_geometry(groups, sealed=True, phosphorylated=True))
    g.react("phosphorylate_bead_oligos", "Phosphorylate the barcode5-prime end after ligation", S.PHOSPHORYLATION,
        [(ligated_beads, R.TEMPLATE), (pnk, R.CATALYST)], [phosphorylated])
    bead_support = g.input("bead_support", "M270 streptavidin beads for capture-oligo immobilization", K.SELECTION)
    beads = g.product("beads", "source capture-bead classes; physical biotin attachments not positioned", K.MOLECULAR,
        state=phosphorylated.state)
    bead_wash = g.product("bead_wash", "unbound bead-preparation material", K.WASTE)
    g.react("immobilize_bead_oligos", "Immobilize the biotinylated duplex classes on streptavidin beads", S.AFFINITY_CAPTURE,
        [(phosphorylated, R.TEMPLATE), (bead_support, R.SELECTION)], [beads, bead_wash])
    g.discard(bead_wash, "SI p3 and detailed bead protocol specify washing after immobilization; no recovery fraction or exact composition of unbound material is assigned.")
    capture_ligase = g.input("capture_ligase", "source Taq ligase for capture", K.CATALYST)
    associated_state = SeparateFragments(groups, "bead-associated; attachment unexpanded", True, exposed_capture=True)
    associated = g.product("associated", "bead-associated oligo classes; no covalent gene assembly", K.MOLECULAR, state=associated_state)
    g.react("capture", "Associate oligos with their declared capture-group bead class", S.BINDING,
        [(processed, R.TEMPLATE), (beads, R.STRUCTURAL), (capture_ligase, R.CATALYST)], [associated])
    oil = g.input("emulsion", "Bio-Rad Droplet Generation Oil; compartment membership schematic", K.CONDITION)
    droplets = [g.product("droplet_" + group.identifier, "schematic droplet class " + group.identifier,
        K.MOLECULAR, state=SeparateFragments((group,), "assigned droplet class; actual occupancy unresolved", True, exposed_capture=True)) for group in groups]
    bsa = g.input("bsa", "BSA in the source emulsion assembly mixture; dose unassigned", K.MOLECULAR)
    mix = g.input("pca_polymerase", "Actual KAPA2G Robust Mastermix supply; class-level allocations unassigned", K.MOLECULAR)
    pca = g.product("pca_enzyme", "Polymerase component of that same KAPA2G supply", K.CATALYST)
    nts = g.product("pca_nt", "Nucleotide component of that same KAPA2G supply", K.MOLECULAR)
    buf = g.product("pca_buffer", "Buffer component of that same KAPA2G supply", K.CONDITION)
    g.react("pca_components", "Expose components of the one named assembly Mastermix", S.MAGIC,
        [(mix, R.TEMPLATE)], [pca, nts, buf], changes_molecule=False,
        reason="Inventory decomposition of the actual source Mastermix, not three fresh additions. Exact formulation, total amounts and allocations among schematic droplet classes are unassigned.")
    g.react("partition", "Place the source assembly mixture into schematic compartments", S.PARTITION,
        [(associated, R.TEMPLATE), (oil, R.CONDITION), (bsa, R.STRUCTURAL), (nts, R.STRUCTURAL), (buf, R.CONDITION)], droplets, changes_molecule=False,
        reason="SI p4 mixes captured beads with Mastermix, primers, BSA and BtsI before emulsification. Reagent ports on the ensuing release/PCA are co-present inventory views, not additions to sealed droplets; quantitative distribution and occupancy remain unassigned.")
    nuclease = g.input("release_enzyme", "Co-present BtsI in source emulsion mix; cut coordinates unassigned", K.CATALYST)
    assembly_f = g.input("assembly_f", "Co-present source20mer assembly forward primer; selected bases unassigned", K.OLIGO)
    assembly_r = g.input("assembly_r", "Co-present source20mer assembly reverse primer; selected bases unassigned", K.OLIGO)
    break_reagent = g.input("break_reagent", "chloroform for emulsion break", K.CONDITION)
    recovered = []
    for group, droplet in zip(groups, droplets):
        name = group.identifier
        released = g.product("released_" + name, "separate released fragments " + name, K.MOLECULAR,
            state=SeparateFragments((group,), "released fragment; overlap geometry unassigned", True, False))
        residual = g.product("residual_" + name, "bead/capture remnants " + name, K.MOLECULAR, state=bead_state((group,)))
        g.react("release_" + name, "Release DNA within droplet class " + name, S.DIGESTION,
            [(droplet, R.TEMPLATE), (nuclease, R.CATALYST)], [released, residual])
        product = expected_assembly(group)
        assembled = g.product("assembled_" + name, "INFERRED intended assembly " + name, K.MOLECULAR, state=product)
        g.react("pca_" + name, "INFERRED intended product — polymerase cycling assembly " + name, S.EXTENSION,
            [(released, R.TEMPLATE), (assembly_f, R.STRUCTURAL), (assembly_r, R.STRUCTURAL), (pca, R.CATALYST)], [assembled])
        recovered.append(g.product("recovered_" + name, "INFERRED recovered assembly class " + name, K.MOLECULAR, state=product))
        remainder = g.product("recovery_remainder_" + name, "recovery remainder; physical fate unresolved " + name, K.MOLECULAR)
        inputs = [(assembled, R.TEMPLATE), (residual, R.TEMPLATE), (break_reagent, R.CONDITION)]
        g.react("recover_" + name, "INFERRED intended product — emulsion break/recovery " + name, S.EXTRACTION,
            inputs, [recovered[-1], remainder], changes_molecule=False)
        g.mark_terminal(Terminal(remainder.id, Disposition.UNRESOLVED,
            "The generic SI specifies emulsion break, purification and re-amplification; bead/capture-remnant phase partition and quantitative distributions are not assigned in this graph."))
    pooled = g.product("pooled", "pool of distinct intended assembly classes", K.MOLECULAR,
        state=Drawing(tuple(row for group in groups for row in Scene.duplex(tuple(expected_assembly(group)),
            label=expected_assembly(group).name, unpaired=tuple(seg.name + "'" for seg in expected_assembly(group))).rows())))
    g.react("pool_assemblies", "INFERRED intended products — combine recovered classes without concatenation", S.POOL,
        [(node, R.TEMPLATE) for node in recovered], [pooled], changes_molecule=False)
    size = g.input("size_selection", "source tested post-assembly size-selection option; not universal", K.SELECTION)
    selected = g.product("selected", "INFERRED size-selected assembly classes", K.MOLECULAR, state=pooled.state)
    rejected = g.product("size_waste", "unselected size classes", K.WASTE)
    g.react("size_select", "INFERRED intended products — optional tested size-selection branch", S.SELECTION,
        [(pooled, R.TEMPLATE), (size, R.SELECTION)], [selected, rejected], changes_molecule=False)
    g.discard(rejected, "SI p5 names presence/absence of size selection as a tested factor. This is an optional class-level branch, not a universal source instruction or sequencing-based correctness filter; no yield or specificity follows.")
    reamplified = g.product("reamplified", "INFERRED re-amplified assembly classes", K.MOLECULAR, state=selected.state)
    g.react("reamplify", "INFERRED intended products — re-amplify recovered assemblies", S.PCR,
        [(selected, R.TEMPLATE), (assembly_f, R.STRUCTURAL), (assembly_r, R.STRUCTURAL), (polymerase, R.CATALYST)], [reamplified],
        geometry=PCRGeometry.unavailable(Scene(), reason="SI p4 and detailed emulsion protocol Step10 name assembly-primer re-amplification with KAPA HiFi. Selected primer sequences, assembled product spans and exact binding coordinates remain unassigned. The two schematic classes remain independent molecules."))
    recipient = g.input("recipient", "source expression-plasmid recipient with random twenty-base clone barcode", K.MOLECULAR,
        state=Construct((role("INFERRED recipient drawing flank; placement unassigned", inferred=True), clone_barcode(),
            role("INFERRED remaining recipient context; placement unassigned", inferred=True))))
    expected = expected_clone()
    cloned = g.product("constructed", "INFERRED intended barcode-bearing construct class", K.MOLECULAR,
        state=closed_circle(expected.name, tuple(expected)), chemical_state=expected)
    g.react("clone", "INFERRED construct arrangement — unexpanded source cloning", S.MAGIC,
        [(reamplified, R.TEMPLATE), (recipient, R.STRUCTURAL)], [cloned],
        reason="The main article names cloning into a barcode-bearing expression plasmid. Its exact vector, reaction mechanism and molecular junctions are outside this selected construction scope; a ligase recipe is not inferred from the word cloning.")
    g.constructed(cloned, "Scoped source-expected construct class. Exact reaction products and the published sequencing/QC branch remain unresolved; this is not a completed DropSynth review.")
    g.declare_reviewed_scope(SOURCE_INVENTORY, SOURCE)
    return g.require_publishable()


SOURCE_INVENTORY = (
    "pool", "pool_f", "pool_r", "pool_polymerase", "amplify", "nickase", "expose", "beads", "capture_ligase", "capture",
    "bead_anchor", "bead_ligation", "bead_barcodes", "anneal_bead_oligos", "bead_ligase",
    "ligate_bead_oligos", "bead_pnk", "phosphorylate_bead_oligos", "bead_support", "immobilize_bead_oligos",
    "bulk_f", "bulk_r", "bulk_amplify", "assembly_f", "assembly_r", "reamplify",
    "short_capture", "remove_short", "emulsion", "bsa", "pca_components", "partition", "release_enzyme", "pca_polymerase", "break_reagent", "release_a", "release_b", "pca_a", "pca_b",
    "recover_a", "recover_b", "pool_assemblies", "size_selection", "size_select", "recipient", "clone", "constructed",
)


def render_page():
    graph = protocol_graph()
    return assemble_page(graph, (
        head(TITLE), '<div class="wrap" data-protocol-complete="true">', f"<h1>{TITLE}</h1>",
        '<p class="research-notes"><a href="01_dropsynth.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1126/science.aao5167">Plesa et al., Science (2018)</a>'),
        "<h2>From pooled oligos to separate gene constructs</h2>",
        "<p>DropSynth gathers the fragments for each intended gene on a matching barcoded bead, separates the captured material into water-in-oil compartments, then releases and assembles the fragments. The emulsion is broken and recovered assemblies are amplified and cloned into expression plasmids.</p>",
        info("ILLUSTRATIVE CONSTRUCTION EXAMPLE: A/B are invented group labels within the published four-fragment model class, not experimental target genes or barcode values. The endpoint is a source-expected DNA construct class, not a sequencing-ready library. The full workflow exposes unresolved chemistry and material fates; its typed construct endpoint does not establish sequence correctness, clone recovery or yield."),
        "<h2>What is known and what remains open</h2>",
        "<p>The published bead oligos, their end modifications and common-arm pairing are explicit below. Initial, bulk and recovery amplification stay separate. The assembly mix includes BSA and components of the same named KAPA2G Mastermix before partitioning; distribution among droplets is unassigned. Selected gene sequences, fragment overlaps, release coordinates, primer bodies and cloning junctions remain unassembled.</p>",
        "<p>The twelve-base capture barcode groups fragments for assembly. A distinct random twenty-base barcode later tracks clones. Neither is assigned as a sequencing sample index. The source sequencing and functional assays are outside this construction view; no read lengths or sequencing-primer geometry are invented.</p>",
        info("INFERRED: intended assembled products, clone arrangement and donor junctions are marked. Capture association, compartment membership and covalent assembly remain separate. No equal-N comparison proves binding; no graph step establishes droplet occupancy, specificity, coverage, abundance, protein function or yield. Full oligo layouts, nick/release sites and cloning geometry remain unassigned."),
        WorkflowSection("Capture, compartmentalization and construction", paths=(("pca_enzyme", ()), ("beads", ()), ("associated", ()),
            ("recovered_a", ("capture",)), ("recovered_b", ("capture",)),
            ("constructed", ("recover_a", "recover_b")))),
        OligoSection("Source bead oligos", (
            OligoEntry("anchor (5-prime dual biotin)", (Segment("anchor", BEAD_ANCHOR),)),
            OligoEntry("ligation oligo (5-prime phosphate;3-prime biotin)", (Segment("ligation", BEAD_LIGATION),)),
            OligoEntry("barcode common3-prime arm (preceded by source12mer reverse-complement barcode)", (Segment("common barcode arm", BEAD_COMMON),)),
        ), note_html="Published SI TableS6, original p52. Unknown barcode values are schematic; source anchor and ligation modifications are retained. T4 PNK is applied after common-arm ligation as specified by the detailed bead protocol."),
        OligoSection("Unplaced identifier features", tuple(OligoEntry("capture group " + group.identifier,
            (group.barcode,)) for group in schematic_groups()) + (OligoEntry("later clone barcode", (clone_barcode(),)),),
            note_html="These are feature annotations, not orderable full bead oligos or fragment sequences. The twelve-base capture feature groups input fragments; the later random twenty-base barcode tracks clones. Neither is assigned as a sequencing sample index."),
        "</div>",
    ))


def selftest():
    import re
    from checks import Check, Source, have
    from chemistry_validation import Status, check_graph
    from storage import sources_root
    check = Check(); graph = protocol_graph(); groups = schematic_groups()
    check("no contradictory molecular claims in the scoped construction graph", [c for c in check_graph(graph) if c.status is Status.CONTRADICTORY], [])
    # Identical unknown strings must not merge distinct logical capture groups.
    state = SeparateFragments(groups, "fixture", False)
    check("unknown barcode values preserve eight separate fragment classes", len(state.scene().strands), 8)
    wrong = Fragment("foreign", groups[1].barcode)
    for label, fragments in (("foreign capture feature", (wrong,)),
        ("duplicate fragment class", (groups[0].fragments[0], groups[0].fragments[0]))):
        try: CaptureGroup("a", groups[0].barcode, fragments)
        except ValueError: rejected = True
        else: rejected = False
        check(label + " is rejected during construction", rejected, True)
    final = expected_clone()
    check("capture barcode is not silently retained as the clone barcode",
        [(seg.feature.id, seg.feature.encoding) for seg in final if seg.feature], [("dropsynth.clone", "random")])
    source = Source(sources_root() / Path(__file__).resolve().parents[1].name /
        "science.aao5167.pdf.txt", "Printed final Plesa2018 Science article; external manifest")
    if have(check, source, label="published identifier meanings"):
        text = " ".join(source.path.read_text().split())
        capture = re.search(r"microbead barcodes are distinct (\d+)-nt sequences", text)
        clone = re.search(r"random (\d+)-bp barcode \(assembly barcode\)", text)
        check("capture-feature length matches its own published source clause",
            int(capture[1]) if capture else None, len(groups[0].barcode.top))
        check("later clone-feature length matches its separate published source clause",
            int(clone[1]) if clone else None, len(clone_barcode().top))
    supplement = Source(sources_root() / Path(__file__).resolve().parents[1].name /
        "aao5167-plesa-sm.pdf.txt", "User-supplied published SI, TableS6 rendered original p52")
    if have(check, supplement, label="published bead oligos and modifications"):
        table = supplement.path.read_text().rsplit("\nTable S6.", 1)[-1].split("Table S7.", 1)[0]
        tokens = re.findall(r"(?<![A-Za-z])[ACGT]{8,}(?![A-Za-z])", table)
        check("source ligation oligo transcription", BEAD_LIGATION, tokens[0] if tokens else None)
        check("source wrapped anchor transcription", BEAD_ANCHOR, "".join(tokens[1:3]))
        check("source barcode common arm transcription", BEAD_COMMON, tokens[3] if len(tokens)>3 else None)
        check("source barcode prefix explicitly names twelve-base reverse complement",
            "12-mer microbead barcode reverse complement +" in table, True)
        check("source table explicitly retains end modifications",
            all(value in table for value in ("3’ biotin", "5’ phosphorylation", "5' dual biotin")), True)
    check("only the intended construct endpoint is required", graph.required_output, Disposition.CONSTRUCT)
    check("invented group labels are declared illustrative", graph.representation, "illustrative_example")
    check("BSA enters the mixture before compartmentalization", any(port.material == "bsa" for port in graph.steps["partition"].inputs), True)
    check("named Mastermix is decomposed once, not three fresh supplies", [port.material for port in graph.steps["pca_components"].inputs], ["pca_polymerase"])
    check.report()
