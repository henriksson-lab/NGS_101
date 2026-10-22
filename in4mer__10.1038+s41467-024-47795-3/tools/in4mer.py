"""Published four-guide prototype and Inzolia construction/readout declarations."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Row, Scene, Segment, complement_segments, feature, revcomp
from cohesive_assembly import FivePrimeFragment
import illumina as il
from page import head, info
from pcr import PCRGeometry
from protocol_graph import (InputRole as R, MaterialKind as K, MolecularRelation,
                            ProtocolGraph, SequencingPlan, StepKind as S)
from protocol_model import closed_circle, role
from protocol_views import (OligoEntry, OligoSection, PrimerBindingSection,
                            ReadLayoutSection, SourceRead, WorkflowSection, assemble_page)
import seqprimers as sp

PAGE_FILE = "in4mer.html"
TITLE = "in4mer — four-guide array construction and readout"
SOURCE = "Esmaeili Anvar et al. 2024, doi:10.1038/s41467-024-47795-3; Methods prototype/Inzolia library production and genomic-DNA sequencing; published Supplementary Data10 (MOESM13)"
GUIDE_LENGTH = 20
REPEATS = ("TAATTTCTACTATTGTAGAT", "AAATTTCTACTCTAGTAGAT", "TAATTTCTACTGTCGTAGAT")
TERMINATOR = "TTTTTT"
LEFT_END, RIGHT_END = "AGAT", "GAAT"
TYPE_IIS = "CGTCTC"
PROTOTYPE_F = il.P5[:20]
PROTOTYPE_R = il.P7
INZOLIA_F = "AGGCACTTGCTCGTACGACG"
INZOLIA_R = "ATGTGGGCCCGGCACCTTAA"
GENOMIC_F_SITE = "CTTGTGGAAAGGACGAAACACCG"
INZOLIA_F_SITE = "TGTGGAAAGGACGAAACACCG"
PROTOTYPE_R_SITE = "ACCGACTCGGTGCCACTTTTTCAAGACCAG"
INZOLIA_550_R_SITE = "ACCGACTCGGTGCCACTTTTTCAA"
INZOLIA_052_R_SITE = "TCTACTATTCTTTCCCCTGCACTGT"
CUSTOM_READ1 = "CTTGTGGAAAGGACGAAACACCGGTAATTTCTACTCTTGTAGAT"
SOURCE_TOKEN = "[s]"
PROTOTYPE_CYCLES = {"Read 1": 151, "Index 1 (i7)": 8, "Index 2 (i5)": 8}


def bp(name, bases, **kw):
    return Segment(name, bases, **kw)


def array_core():
    parts = []
    for position in range(1, 5):
        parts.append(Segment(f"guide {position}", "N" * GUIDE_LENGTH, placeholder=True,
            feature=feature(f"in4mer.guide{position}", "guide_barcode", "unknown",
                note=f"Ordered guide-identity window {position}; no target sequence or editing efficiency is selected.")))
        if position <= len(REPEATS):
            parts.append(bp(f"direct repeat {position}", REPEATS[position - 1]))
    return (*parts, bp("source terminator", TERMINATOR))


def synthesis(branch):
    if branch not in {"prototype", "inzolia"}:
        raise ValueError("pool chemistry must be selected from the two source-defined branches")
    forward, reverse = ((PROTOTYPE_F, PROTOTYPE_R) if branch == "prototype" else (INZOLIA_F, INZOLIA_R))
    return Construct((bp("pool amplification forward", forward), bp("forward type-IIS site", TYPE_IIS),
        bp("left cut spacer", "G"), bp("left insert end", LEFT_END), *array_core(),
        bp("right insert end", RIGHT_END), bp("right cut spacer", "G"),
        bp("reverse type-IIS site", revcomp(TYPE_IIS)), bp("pool amplification reverse complement", revcomp(reverse))),
        branch + " synthesized pool architecture")


def physical_length(product):
    return sum(seg.length_bp or len(seg) for seg in product)


def selected_release(branch, template):
    left_cut = template.span("forward type-IIS site")[1] + 1
    right_cut = template.span("reverse type-IIS site")[0] - 5
    names = [s.name for s in template]
    core = tuple(template.segments[names.index("left insert end") + 1:names.index("right insert end")])
    sequence = template.top()
    left, right = sequence[left_cut:left_cut + 4], sequence[right_cut:right_cut + 4]
    if left != LEFT_END or right != RIGHT_END or "".join(s.top for s in core) != sequence[left_cut + 4:right_cut]:
        raise ValueError("the selected source type-IIS cuts must release the declared array and ends")
    return FivePrimeFragment(branch + " released four-guide insert", left, core, right)


def ss(product, label):
    sc = Scene(); sc.strand("strand", tuple(product), label=label)
    return sc


@dataclass(frozen=True)
class UnassignedRecipient:
    """A known released insert and an independently unassigned vector end state."""
    insert: FivePrimeFragment
    vector: str
    availability = "unavailable"
    ready_label = "released insert; recipient ends unassigned"
    action_label = "Present the known insert and unresolved recipient"

    @property
    def statement(self):
        return (f"Methods name assembly into {self.vector}. The insert's ends derive from the printed pool; "
                "the complete receiving-vector map and opened cassette are not assigned in this review. "
                "No recipient sequence or source-absence conclusion is inferred from the insert design.")

    def rows(self):
        return self.insert.rows() + [Row(chunks=[(
            self.vector + " receiving cassette: end geometry unassigned", None, False)])]


def expected_construct(vector, core):
    # These junctions are an explicitly inferred intended product, not source-exact
    # receiving-vector ends supplied to the neutral cohesive assembly checker.
    return Construct((bp("INFERRED left recipient junction", LEFT_END, inferred=True), *core,
                      bp("INFERRED right recipient junction", RIGHT_END, inferred=True),
                      role(vector + " remaining backbone")), "intended " + vector + " array construct")


def sample_index(branch, channel):
    return Segment(channel + " sample index", "N" * 8, placeholder=True,
        feature=feature(f"in4mer.{branch}.{channel}", "sample_index", "unknown",
                        note="Source names eight-base Illumina sample-index sets; no actual sample assignment is selected."))


def genomic_primers(branch):
    if branch not in {"prototype", "inzolia550", "inzolia052"}:
        raise ValueError("genomic primer use must name a declared source/candidate branch")
    is_proto = branch == "prototype"
    fsite = GENOMIC_F_SITE if is_proto else INZOLIA_F_SITE
    forward = [bp("P5", il.P5), sample_index(branch, "i5"), bp("TruSeq Read1 arm", il.TRUSEQ_READ1)]
    if not is_proto:
        forward.append(Segment("literal unassigned source token", SOURCE_TOKEN, placeholder=True))
    forward.append(bp("source genomic forward end", fsite))
    rsite = (PROTOTYPE_R_SITE if is_proto else
             INZOLIA_550_R_SITE if branch.endswith("550") else INZOLIA_052_R_SITE)
    reverse = (bp("P7", il.P7), sample_index(branch, "i7"), bp("TruSeq Read2 arm", il.TRUSEQ_READ2),
               bp("source genomic reverse end", rsite))
    return tuple(forward), reverse


def library(branch):
    forward, reverse = genomic_primers(branch)
    # Prototype's custom sequencing primer supplies evidence for an expected
    # cassette prefix. Its complete vector-site placement is not verified.
    if branch == "prototype":
        prefix = CUSTOM_READ1[len(GENOMIC_F_SITE):-len(LEFT_END)]
        context = (bp("INFERRED pre-guide vector context", prefix, inferred=True),)
    else:
        context = (role("unassigned vector-specific pre-guide interval"),)
    return Construct((*forward, *context,
        bp("INFERRED left cassette junction", LEFT_END, inferred=True), *array_core(),
        bp("INFERRED right cassette junction", RIGHT_END, inferred=True),
        role("unassigned downstream vector interval"), *complement_segments(reverse)), branch + " genomic array library")


def assert_array(product):
    guides = [seg for seg in product if seg.feature and seg.feature.role == "guide_barcode"]
    if [(s.feature.id, len(s.top)) for s in guides] != [
            (f"in4mer.guide{i}", GUIDE_LENGTH) for i in range(1, 5)]:
        raise ValueError("four distinct ordered guide windows must be preserved without deduplication")


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    bsmbi = g.input("bsmbi", "prototype BsmBI-v2 assembly kit", K.CATALYST)
    kit_ligase = g.input("kit_ligase", "prototype kit ligase; exact identity unassigned", K.CATALYST)
    esp3i = g.input("esp3i", "Inzolia Esp3I", K.CATALYST)
    t7 = g.input("t7", "Inzolia T7 DNA ligase", K.CATALYST)
    kapa = g.input("kapa", "prototype KAPA HiFi pool PCR", K.CATALYST)
    neb = g.input("neb", "Inzolia NEBNext high-fidelity pool PCR", K.CATALYST)
    taq = g.input("taq", "Titanium Taq genomic PCR", K.CATALYST)
    mono = g.input("monarch", "prototype Monarch cleanup", K.SELECTION)
    qia = g.input("qiaquick", "Inzolia Qiaquick PCR cleanup", K.SELECTION)
    egel = g.input("egel", "prototype E-Gel SizeSelect", K.SELECTION)
    spri = g.input("spri", "Inzolia AMPure XP cleanup", K.SELECTION)
    carb = g.input("carbenicillin", "source carbenicillin clone selection", K.SELECTION)
    branches = (("prototype", "pRDA-550", bsmbi, kit_ligase, kapa, mono),
                ("inzolia550", "pRDA-550", esp3i, t7, neb, qia),
                ("inzolia052", "pRDA-052", esp3i, t7, neb, qia))
    # Inzolia's shared molecular preparation is shown twice for its two receiving
    # vectors; these are alternative source branches, not two successive reactions.
    for branch, vector_name, enzyme, ligase, polymerase, clean_reagent in branches:
        pool_kind = "prototype" if branch == "prototype" else "inzolia"
        template = synthesis(pool_kind)
        source = g.input(branch + "_synthesis", branch + " synthesized four-guide pool", K.OLIGO,
                         state=ss(template, "synthesized array strand"), chemical_state=template, primary=True)
        primer_seqs = ((PROTOTYPE_F, PROTOTYPE_R) if pool_kind == "prototype" else (INZOLIA_F, INZOLIA_R))
        primers = [g.input(branch + "_" + key, branch + " pool " + key, K.OLIGO,
                    state=Construct((bp(key, seq),))) for key, seq in zip(("forward", "reverse"), primer_seqs)]
        amplicon = g.product(branch + "_amplified", branch + " pool amplicon", K.MOLECULAR, state=template)
        g.react(branch + "_pool_pcr", "Amplify the " + branch + " pool with its own handles", S.PCR,
                [(source, R.TEMPLATE), *[(p, R.STRUCTURAL) for p in primers], (polymerase, R.CATALYST)], [amplicon],
                geometry=PCRGeometry.unavailable(source.state,
                    reason="The intended pool primer sites are explicit; the initial complementary copy of the synthesized strand is not expanded into a pre-existing duplex geometry."))
        cleaned = g.product(branch + "_cleaned", branch + " purified pool amplicon", K.MOLECULAR, state=template)
        pcr_waste = g.product(branch + "_pool_waste", branch + " unretained pool-PCR material", K.WASTE)
        g.react(branch + "_pool_clean", "Purify the " + branch + " pool amplicon", S.CLEANUP,
                [(amplicon, R.TEMPLATE), (clean_reagent, R.SELECTION)], [cleaned, pcr_waste],
                changes_molecule=False, relations=(MolecularRelation.STRAND_ORDER,))
        g.discard(pcr_waste, "Material not retained in source cleanup; no recovery proportion is inferred.")
        recipient = g.input(branch + "_recipient", branch + " " + vector_name + " recipient", K.MOLECULAR,
                            state=Construct((role(vector_name + " source cassette and backbone unassigned"),)))
        insert = selected_release(branch, template)
        insert_node = g.product(branch + "_insert", insert.name, K.MOLECULAR, state=insert, chemical_state=insert)
        opened = g.product(branch + "_opened", vector_name + " opened recipient, ends unassigned", K.MOLECULAR,
                           state=Construct((role("source opened recipient; ends unassigned"),)))
        cut_waste = g.product(branch + "_cut_waste", branch + " unretained terminal/cassette fragments", K.WASTE)
        g.react(branch + "_cut", "Cleavage within the " + branch + " source assembly", S.DIGESTION,
                [(cleaned, R.TEMPLATE), (recipient, R.TEMPLATE), (enzyme, R.CATALYST)],
                [insert_node, opened, cut_waste])
        g.discard(cut_waste, "Unretained digest classes; complete vector cleavage and fragment removal are not established.")
        product = expected_construct(vector_name, array_core()); assert_array(product)
        assembled = g.product(branch + "_construct", "INFERRED recipient junctions — " + branch + " expected construct",
                              K.MOLECULAR, state=closed_circle(product.name, tuple(product)), chemical_state=product)
        g.react(branch + "_join", "INFERRED recipient junctions — " + branch + " assembly", S.LIGATION,
                [(insert_node, R.TEMPLATE), (opened, R.STRUCTURAL), (ligase, R.CATALYST)], [assembled],
                geometry=UnassignedRecipient(insert, vector_name))
        recovered = g.product(branch + "_plasmid", branch + " expected recovered plasmid class", K.MOLECULAR,
                              state=assembled.state, chemical_state=product)
        rejected = g.product(branch + "_clone_waste", branch + " unretained clone-selection material", K.WASTE)
        g.react(branch + "_recover", "INFERRED recipient junctions — recover " + branch + " selected plasmids", S.SELECTION,
                [(assembled, R.TEMPLATE), (carb, R.SELECTION)], [recovered, rejected], changes_molecule=False)
        g.discard(rejected, "Source carbenicillin selection/recovery is represented at class level; no copy fidelity, population dynamics, genotype, efficiency or representation is established.")
        samples = g.input(branch + "_samples", branch + " source assay samples; biological state unassigned",
                          K.BIOLOGICAL, primary=True)
        genomic = g.product(branch + "_genomic", branch + " post-assay genomic DNA", K.MOLECULAR,
                            state=Construct((role("unassigned recovered genomic DNA"),)))
        g.react(branch + "_assay_boundary", "Unexpanded source assay and genomic recovery boundary", S.MAGIC,
                [(recovered, R.TEMPLATE), (samples, R.TEMPLATE)], [genomic],
                reason="The reported downstream biological assay and genomic recovery are an explicit unmodelled boundary. No delivery, integration, guide processing, editing, target specificity, genomic identity or molecular copy relation is asserted.")
        fparts, rparts = genomic_primers(branch)
        candidate = branch != "prototype"
        prefix = "candidate experiment assignment — " if candidate else ""
        forward = g.input(branch + "_genomic_f", prefix + branch + " indexed genomic forward", K.OLIGO,
                          state=Construct(fparts), inferred=candidate)
        reverse = g.input(branch + "_genomic_r", prefix + branch + " indexed genomic reverse", K.OLIGO,
                          state=Construct(rparts), inferred=candidate)
        lib = library(branch); assert_array(lib)
        amplified = g.product(branch + "_genomic_amp", prefix + "INFERRED cassette context — " + lib.name,
                              K.LIBRARY, state=lib, inferred=candidate)
        g.react(branch + "_genomic_pcr", prefix + "INFERRED cassette context — make " + branch + " genomic library", S.PCR,
                [(genomic, R.TEMPLATE), (forward, R.STRUCTURAL), (reverse, R.STRUCTURAL), (taq, R.CATALYST)], [amplified],
                geometry=PCRGeometry.unavailable(genomic.state,
                    reason="Full vector/genomic binding coordinates remain unassigned. Inzolia forward primers also retain the literal unresolved [s] token; no meaning, sequence or length is supplied for it."))
        ready = g.product(branch + "_library", prefix + branch + " purified genomic array library", K.LIBRARY,
                          state=lib, inferred=candidate)
        waste = g.product(branch + "_genomic_waste", branch + " unretained genomic-PCR material", K.WASTE)
        g.react(branch + "_genomic_clean", "INFERRED cassette context — purify " + branch + " genomic products", S.CLEANUP,
                [(amplified, R.TEMPLATE), (egel if branch == "prototype" else spri, R.SELECTION)], [ready, waste],
                changes_molecule=False, relations=(MolecularRelation.STRAND_ORDER,))
        g.discard(waste, "Unretained cleanup material; yields and abundance distributions are not modelled.")
        if branch == "prototype":
            read1 = sp.custom("Read 1", "source prototype custom forward sequencing", CUSTOM_READ1, SOURCE)
            seqprimers = (read1, sp.TRUSEQ["I1"], sp.TRUSEQ["I2"])
            spans = [(s.feature.id, s.cycle_start, s.cycle_end) for s in sp.feature_spans(lib, seqprimers, PROTOTYPE_CYCLES)
                     if s.feature.role == "guide_barcode"]
            expected = [(f"in4mer.guide{i}", (i - 1) * GUIDE_LENGTH + sum(map(len, REPEATS[:i - 1])) + 1,
                         i * GUIDE_LENGTH + sum(map(len, REPEATS[:i - 1]))) for i in range(1, 5)]
            if spans != expected:
                raise ValueError("the reported prototype read must cover all four ordered guide windows")
            g.sequence(ready, SequencingPlan(library=lib, primers=seqprimers,
                acquired_roles=tuple(PROTOTYPE_CYCLES), read_lengths=PROTOTYPE_CYCLES,
                read_uses={"Read 1":"extract", "Index 1 (i7)":"demultiplex", "Index 2 (i5)":"demultiplex"},
                limitation="151–8–8 single-end prototype case. The custom Read1 template context is inferred from its reported primer; index sites derive from canonical PCR arms. Full genomic binding coordinates and sample assignments remain unassigned. The 7mer paired-end recipe is a different source branch."))
        else:
            seqprimers = sp.standard_primers(lib)
            g.sequence(ready, SequencingPlan(library=lib, primers=seqprimers, acquisition="unspecified",
                limitation="Inzolia source names HiSeq2500 Rapid Run, 200-cycle kit; this review does not assign reads or per-read cycles from that capacity. Displayed primer sites are potential canonical sites only. The repeated workbook heading, vector-specific assignment, literal [s] token and cassette/genomic context remain unresolved; no guide-cycle coordinates or sample-index values are inferred."))
    g.declare_reviewed_scope((
        "bsmbi", "kit_ligase", "esp3i", "t7", "kapa", "neb", "taq", "monarch", "qiaquick", "egel", "spri", "carbenicillin",
        "prototype_synthesis", "prototype_forward", "prototype_reverse", "prototype_pool_pcr", "prototype_pool_clean",
        "prototype_recipient", "prototype_cut", "prototype_join", "prototype_recover", "prototype_plasmid",
        "prototype_samples", "prototype_assay_boundary", "prototype_genomic_f", "prototype_genomic_r", "prototype_genomic_pcr", "prototype_genomic_clean", "prototype_library",
        "inzolia550_synthesis", "inzolia550_forward", "inzolia550_reverse", "inzolia550_pool_pcr", "inzolia550_pool_clean",
        "inzolia550_recipient", "inzolia550_cut", "inzolia550_join", "inzolia550_recover", "inzolia550_plasmid",
        "inzolia550_samples", "inzolia550_assay_boundary", "inzolia550_genomic_f", "inzolia550_genomic_r", "inzolia550_genomic_pcr", "inzolia550_genomic_clean", "inzolia550_library",
        "inzolia052_synthesis", "inzolia052_forward", "inzolia052_reverse", "inzolia052_pool_pcr", "inzolia052_pool_clean",
        "inzolia052_recipient", "inzolia052_cut", "inzolia052_join", "inzolia052_recover", "inzolia052_plasmid",
        "inzolia052_samples", "inzolia052_assay_boundary", "inzolia052_genomic_f", "inzolia052_genomic_r", "inzolia052_genomic_pcr", "inzolia052_genomic_clean", "inzolia052_library",
    ), SOURCE)
    return g.require_publishable()


def render_page():
    g = protocol_graph()
    content = [head(TITLE), '<div class="wrap" data-protocol-complete="true">', f"<h1>{TITLE}</h1>",
        '<p class="research-notes"><a href="01_in4mer.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1038/s41467-024-47795-3">Esmaeili Anvar et al., Nature Communications (2024)</a>'),
        info("Four-guide prototype and two Inzolia recipient alternatives. Guide targets are unassigned; ordered identity windows and distinct repeats are retained. INFERRED: intended recipient junctions and the prototype pre-guide sequencing context are marked. Full vector maps, opened recipient ends and genomic binding coordinates remain unassigned. The downstream assay/recovery boundary is explicitly unmodelled; no delivery, guide processing, editing, genotype fidelity, abundance or coverage is inferred. Inzolia's duplicated preparation represents alternative recipients, not sequential assembly."),
        info("Supplementary Data10 (MOESM13) repeats the prototype/7mer heading above its vector-specific primer block. Those pRDA-550/pRDA-052 rows are shown as candidates for the matching Inzolia recipients; experiment assignment remains unresolved. Their literal [s] token is preserved with unassigned meaning and length. The 7mer-only reverse sequencing primer and paired-end acquisition are not assigned to the four-guide prototype."),
        WorkflowSection("Construction and genomic recovery", paths=(
            ("prototype_plasmid", ()), ("prototype_library", ("prototype_recover",)),
            ("inzolia550_plasmid", ()), ("inzolia550_library", ("inzolia550_recover",)),
            ("inzolia052_plasmid", ()), ("inzolia052_library", ("inzolia052_recover",)),
        )),
        OligoSection("Pool amplification primers", tuple(OligoEntry(name, (bp(name, sequence),)) for name, sequence in (
            ("prototype forward", PROTOTYPE_F), ("prototype reverse", PROTOTYPE_R),
            ("Inzolia forward", INZOLIA_F), ("Inzolia reverse", INZOLIA_R)))),
        OligoSection("Indexed genomic primers", tuple(entry for branch in ("prototype", "inzolia550", "inzolia052")
            for entry in (OligoEntry(branch + " forward", genomic_primers(branch)[0]),
                          OligoEntry(branch + " reverse", genomic_primers(branch)[1])))),
        PrimerBindingSection("prototype_library", "Prototype sequencing", "INFERRED pre-guide context — prototype primer sites"),
        ReadLayoutSection("prototype_library", (
            SourceRead("Read 1", "151", "four ordered guide windows and intervening repeats", "Read 1", "extract"),
            SourceRead("Index 1 (i7)", "8", "sample index", "Index 1 (i7)", "demultiplex"),
            SourceRead("Index 2 (i5)", "8", "sample index", "Index 2 (i5)", "demultiplex")),
            "Reported NextSeq500 single-end prototype format; guide positions derive from the declared construct, not from an editing or representation measurement."),
    ]
    for branch in ("inzolia550", "inzolia052"):
        content.append(PrimerBindingSection(branch + "_library", branch + " potential primer sites",
                       "Potential canonical sites; acquisition unassigned", show_read_content=False))
    content.append("</div>")
    return assemble_page(g, content)


def selftest():
    import csv
    from checks import Check, Source, have
    from chemistry_validation import Status, check_graph
    from plasmid import Enzyme, Plasmid, digest
    from storage import sources_root
    check = Check(); graph = protocol_graph()
    check("source-scoped array declarations have no contradictory molecular claims",
          [c for c in check_graph(graph) if c.status is Status.CONTRADICTORY], [])
    original = list(array_core())
    original[0], original[2] = original[2], original[0]
    check.raises("guide positions cannot be swapped by treating a combination as a set",
                 lambda: assert_array(Construct(original)))
    for branch in ("prototype", "inzolia"):
        # N windows leave hidden internal sites unresolved. This checks the known
        # terminal recognition sites, not a complete digest of every library member.
        cuts = digest(Plasmid(branch, synthesis(branch).top(), circular=False), Enzyme("source type-IIS", TYPE_IIS, 1, 5))
        check(branch + " source-declared terminal cuts give the intended ends",
              [cut.overhang for cut in cuts], [LEFT_END, RIGHT_END])
    for branch in ("inzolia550", "inzolia052"):
        primer = graph.materials[branch + "_genomic_f"].state
        check(branch + " preserves the literal source token", SOURCE_TOKEN in primer.top(), True)
        check(branch + " does not assign a run recipe from kit capacity",
              (graph.sequencing[branch + "_library"].acquisition, graph.sequencing[branch + "_library"].read_lengths),
              ("unspecified", {}))
    source = Source(sources_root() / Path(__file__).resolve().parents[1].name /
        "s41467-024-47795-3_41467_2024_47795_MOESM13_ESM.xlsx.txt",
        "Published Supplementary Data10 (MOESM13), doi:10.1038/s41467-024-47795-3; tools/get_sources.py and tools/doctext.py")
    if have(check, source, label="literal candidate primer transcription"):
        rows = list(csv.reader(source.path.read_text().split("\n")[1:], delimiter="\t"))
        def normalise(value):
            return value.replace(" ", "").replace("n", "N").removeprefix("5’-").removeprefix("5'-").removesuffix("-3’").removesuffix("-3'")
        entries = {row[0].strip(): normalise(row[1]) for row in rows if len(row) > 1}
        for branch, vector in (("inzolia550", "pRDA-550"), ("inzolia052", "pRDA-052")):
            for direction, parts in zip(("forward", "reverse"), genomic_primers(branch)):
                key = vector + " compatible " + direction + " primer:"
                check(key + " literal row is preserved, without resolving experiment assignment",
                      Construct(parts).top(), entries[key])
    check.report()
