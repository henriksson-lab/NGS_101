"""CHyMErA pooled cloning and readout, explicitly the 2021 detailed method."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Row, Scene, Segment, complement_segments, feature, revcomp
from cohesive_assembly import FivePrimeFragment
import illumina as il
from page import head, info
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole as R, MaterialKind as K,
                            ProtocolGraph, SequencingPlan, StepKind as S)
from protocol_model import closed_circle, role
from protocol_views import (OligoEntry, OligoSection, PrimerBindingSection,
                            WorkflowSection, assemble_page)
import seqprimers as sp

TITLE = "CHyMErA — pooled cloning and sequencing"
PAGE_FILE = "chymera.html"
SOURCE = "Aregger et al. 2021, doi:10.1038/s41596-021-00595-1; pooled cloning Steps34–69, NGS Steps126–140, Table1 and Figures4–5 captions. Gonatopoulos-Pournatzis et al.2020 doi:10.1038/s41587-020-0437-z Supplementary Table14 independently supplies the named primer rows. The newly supplied 2020 main Methods confirm two assembly rounds but report a distinct31/30 imaged-cycle recipe, not the2021 27/27 recipe modelled here."
TKO_F = "TGTCAGTTGTCATTCGCGAAAAAGAGAACCTGCAGAGACC"
TKO_R = "GTCACTGACGCGGTTTTGTAGATCTCTACCTGCCTCTAAA"
OUTER_F = "TGAAGAATCGCAAAACCAGCA"
OUTER_R = "GAGGGCCTATTTCCCATGATTC"
F_SITE = "CCGAGGGGACCCAGAGAATTC"
R_SITE = "GTGGAAAGGACGAGGTACCG"
LEFT_SYNTH = "agagaACCTGCagagaccg"
STUFFER = "gtttaGAGACGgctaaatccgCGTCTCgagat"
RIGHT_SYNTH = "ttttagagGCAGGTagaga"
I5 = "TGCTAAGT"  # A502-Cas12a-F1, called A502-Cpf1-F1 in the 2020 table.
I7 = "ATCACGAC"  # A701-Cas9-R1.


def bp(name, bases, **kw):
    return Segment(name, bases.upper(), **kw)


def guide(which):
    if which not in ("Cas9", "Cas12a"):
        raise ValueError("guide identity must name one of the two source spacer roles")
    return Segment(which + " spacer", "N" * (20 if which == "Cas9" else 23), placeholder=True,
        feature=feature("chymera." + which.lower(), "guide_barcode", "unknown", group="encoded guide pair",
            part=which, note="Ordered source spacer-identity window; no target or editing outcome is selected."))


def synthesis():
    return Construct((bp("synthesis left", LEFT_SYNTH), guide("Cas9"), bp("removable stuffer", STUFFER),
        guide("Cas12a"), bp("synthesis right", RIGHT_SYNTH)), "2021 synthesized pooled-oligo architecture")


def synthesis_scene():
    sc = Scene()
    sc.strand("oligo", tuple(synthesis()), label="source single-stranded oligo pool")
    return sc


def pool_geometry():
    return PCRGeometry.exact(synthesis(), PCRPrimer("TKO_F1", TKO_F, 15),
        PCRPrimer("TKO_R1", TKO_R, 15), statement="2021 Step41 supplies TKO_F1/TKO_R1 for the Step34 oligo template. Exact terminal sites only; spacer values remain unassigned.")


def amplified_pool():
    geometry = pool_geometry(); template = synthesis()
    start = template.span("Cas9 spacer")[0]; end = template.span("Cas12a spacer")[1]
    left = template.top()[geometry.forward_hit.end:start]
    right = template.top()[end:geometry.reverse_hit.start]
    if geometry.forward_hit.end > start or geometry.reverse_hit.start < end:
        raise ValueError("source amplification must preserve both complete spacer windows")
    parts = [bp("TKO_F1 written end", TKO_F)]
    if left: parts.append(bp("retained left oligo flank", left))
    parts.extend((guide("Cas9"), bp("removable stuffer", STUFFER), guide("Cas12a")))
    if right: parts.append(bp("retained right oligo flank", right))
    parts.append(bp("TKO_R1 written end", revcomp(TKO_R)))
    return Construct(tuple(parts), "amplified paired-spacer pool")


def paired_insert():
    return (guide("Cas9"), bp("removable stuffer", STUFFER), guide("Cas12a"))


def assert_pair(product):
    observed = [(s.feature.id, len(s.top)) for s in product
                if s.feature and s.feature.role == "guide_barcode"]
    if observed != [("chymera.cas9", 20), ("chymera.cas12a", 23)]:
        raise ValueError("both distinct ordered spacer windows must survive construction and readout")


def intermediate():
    product = Construct((*paired_insert(), role("right recipient junction"),
        role("pLCHKO remaining backbone"), role("left recipient junction")), "intended first-round construct")
    assert_pair(product)
    return product


def opened_intermediate():
    # CGTCTC(1/5) is the shared BsmBI/Esp3I operation. Only the two
    # source-written stuffer sites are selected; a full vector digest is not claimed.
    sequence = STUFFER.upper()
    left_site = sequence.index(revcomp("CGTCTC"))
    right_site = sequence.index("CGTCTC") + len("CGTCTC")
    left_end = sequence[left_site - 5:left_site - 1]
    right_end = sequence[right_site + 1:right_site + 5]
    return FivePrimeFragment("selected BsmBI-opened intermediate", right_end,
        (guide("Cas12a"), role("right recipient junction"), role("pLCHKO remaining backbone"),
         role("left recipient junction"), guide("Cas9")), left_end)


def completed_cassette():
    opened = opened_intermediate()
    return (guide("Cas9"), bp("INFERRED donor-side left fusion", opened.right, inferred=True),
        role("SpCas9 tracrRNA and LbCas12a direct-repeat core; donor sequence unassigned"),
        bp("INFERRED donor-side right fusion", opened.left, inferred=True), guide("Cas12a"))


def completed_construct():
    product = Construct((*completed_cassette(), role("right recipient junction"),
        role("pLCHKO remaining backbone"), role("left recipient junction")), "intended completed hybrid-guide construct")
    assert_pair(product)
    return product


def index_segment(channel, value):
    return Segment(channel + " source example index", value,
        feature=feature("chymera.example." + channel, "sample_index", "fixed",
            note="A source-listed primer index. This displayed i5/i7 combination is illustrative, not an assigned experimental sample."))


def index_primers():
    forward = (bp("P5", il.P5), index_segment("i5", I5), bp("TruSeq Read1 arm", il.TRUSEQ_READ1),
               bp("U6-side annealing end", F_SITE))
    reverse = (bp("P7", il.P7), index_segment("i7", I7), bp("TruSeq Read2 arm", il.TRUSEQ_READ2),
               bp("termination-side annealing end", R_SITE))
    return forward, reverse


def library():
    forward, reverse = index_primers()
    product = Construct((*forward, role("unassigned pre-spacer vector interval"), *completed_cassette(),
        role("unassigned post-spacer vector interval"), *complement_segments(reverse)), "illustrative indexed paired-guide library")
    assert_pair(product)
    return product


@dataclass(frozen=True)
class UnassignedAssembly:
    statement: str
    known: FivePrimeFragment | None = None
    availability = "unavailable"
    ready_label = "source assembly classes; incomplete end assignment"
    action_label = "Present source assembly inputs"

    def rows(self):
        return (self.known.rows() if self.known else []) + [Row(chunks=[(
            "Receiving vector/donor maps and complementary end assignment remain unresolved", None, False)])]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    oligo = synthesis(); amp = amplified_pool()
    pool = g.input("pool", oligo.name, K.OLIGO, state=synthesis_scene(), primary=True)
    tf = g.input("tko_f", "TKO_F1", K.OLIGO, state=Construct((bp("TKO_F1", TKO_F),)))
    tr = g.input("tko_r", "TKO_R1", K.OLIGO, state=Construct((bp("TKO_R1", TKO_R),)))
    q5 = g.input("q5", "source Q5 high-fidelity polymerase", K.CATALYST)
    ngs_q5 = g.input("ngs_q5", "source NEBNext Ultra II Q5 Master Mix", K.CATALYST)
    cleanup = g.input("cleanup", "source PCR/ligation cleanup", K.SELECTION)
    gel = g.input("gel", "source gel-band isolation", K.SELECTION)
    antibiotic = g.input("antibiotic", "source ampicillin or carbenicillin clone selection", K.SELECTION)
    bfu = g.input("bfuai", "BfuAI/BveI", K.CATALYST)
    bsm = g.input("bsmbi", "BsmBI/Esp3I", K.CATALYST)
    ligase = g.input("ligase", "source T4 DNA ligase", K.CATALYST)
    sap = g.input("sap", "source rSAP phosphatase", K.CATALYST)
    recipient = g.input("recipient", "source pLCHKO; complete map/ends unassigned", K.MOLECULAR,
                        state=Construct((role("pLCHKO recipient"),)))
    opened = g.product("opened_recipient", "BfuAI-opened recipient; exact ends unassigned", K.MOLECULAR,
                       state=Construct((role("opened pLCHKO backbone"),)))
    stuffer = g.product("recipient_stuffer", "source recipient stuffer fallout", K.WASTE)
    g.react("open_recipient", "Open the source recipient independently of the oligo PCR", S.DIGESTION,
            [(recipient, R.TEMPLATE), (bfu, R.CATALYST)], [opened, stuffer])
    g.discard(stuffer, "Recipient stuffer removed; complete vector cleavage pattern remains unassigned.")
    dephos = g.product("recipient_dephos", "rSAP-treated receiving backbone", K.MOLECULAR, state=opened.state)
    g.react("recipient_dephosphorylate", "Treat the receiving backbone with source rSAP", S.DEPHOSPHORYLATION,
            [(opened, R.TEMPLATE), (sap, R.CATALYST)], [dephos])
    def select(identifier, source, reagent, state, *, label, kind=S.CLEANUP, material_kind=K.MOLECULAR):
        retained = g.product(identifier, label, material_kind, state=state)
        waste = g.product(identifier + "_waste", label + " unretained material", K.WASTE)
        g.react(identifier + "_step", label, kind, [(source, R.TEMPLATE), (reagent, R.SELECTION)],
                [retained, waste], changes_molecule=False)
        g.discard(waste, "Source physical selection/cleanup; no recovery fraction or genotype is established.")
        return retained
    ready = select("recipient_ready", dephos, gel, dephos.state, label="Gel-isolate the opened recipient", kind=S.SELECTION)
    amplified = g.product("amplified_pool", amp.name, K.MOLECULAR, state=amp)
    g.react("amplify_pool", "Amplify the paired-spacer pool", S.PCR,
            [(pool, R.TEMPLATE), (tf, R.STRUCTURAL), (tr, R.STRUCTURAL), (q5, R.CATALYST)],
            [amplified], geometry=PCRGeometry.unavailable(pool.state,
                reason="The source pool starts as single-stranded oligos. Complementary-strand synthesis within the reported PCR is unexpanded. Product-end calculation uses the known source reference axis and primer tails; a starting duplex or complete first-cycle binding geometry is not asserted."))
    clean = select("pool_clean", amplified, cleanup, amp, label="Purify the oligo amplicon")
    insert = g.product("first_insert", "BfuAI-released spacer insert; ends unassigned", K.MOLECULAR,
                       state=Construct(paired_insert(), "paired-spacer insert class"))
    removed = g.product("pool_ends", "unretained synthesis/PCR terminal classes", K.WASTE)
    g.react("cut_pool", "Release the pooled insert for the first assembly", S.DIGESTION,
            [(clean, R.TEMPLATE), (bfu, R.CATALYST)], [insert, removed])
    g.discard(removed, "Terminal synthesis/amplification material excluded from the intended insert; exact released end chemistry remains unassigned.")
    first = intermediate()
    joined = g.product("first_construct", first.name, K.MOLECULAR, state=closed_circle(first.name, tuple(first)))
    g.react("first_assembly", "First assembly: insert the paired spacers and internal stuffer", S.LIGATION,
            [(ready, R.TEMPLATE), (insert, R.STRUCTURAL), (ligase, R.CATALYST)], [joined],
            geometry=UnassignedAssembly("2021 pooled assembly1; recipient maps and complete release geometry are unassigned. No complete cohesive-product relation is claimed."))
    purified = select("first_cleanup", joined, cleanup, joined.state, label="Purify first-assembly products")
    first_pool = select("first_pool", purified, antibiotic, joined.state, label="Recover the first selected plasmid class", kind=S.SELECTION)
    fragment = opened_intermediate()
    cut = g.product("opened_intermediate", fragment.name, K.MOLECULAR, state=fragment, chemical_state=fragment)
    fallout = g.product("internal_stuffer", "internal BsmBI stuffer fallout", K.WASTE)
    g.react("second_cut", "Open the intermediate at its source-written BsmBI sites", S.DIGESTION,
            [(first_pool, R.TEMPLATE), (bsm, R.CATALYST)], [cut, fallout])
    g.discard(fallout, "Internal oligo stuffer is removed before the scaffold insert; no complete vector-digest or off-product distribution is asserted.")
    cut_dephos = g.product("intermediate_dephos", "rSAP-treated opened intermediate", K.MOLECULAR, state=fragment)
    g.react("intermediate_dephosphorylate", "Treat the opened intermediate with source rSAP", S.DEPHOSPHORYLATION,
            [(cut, R.TEMPLATE), (sap, R.CATALYST)], [cut_dephos])
    cut_ready = select("intermediate_ready", cut_dephos, cleanup, fragment, label="Purify the opened intermediate")
    donor = g.input("scaffold_donor", "source TOPO-tracrRNA-LbDR donor; map/ends unassigned", K.MOLECULAR,
                    state=Construct((role("TOPO donor carrying tracrRNA and LbCas12a direct repeat"),)))
    released = g.product("scaffold_insert", "released scaffold class; exact donor ends unassigned", K.MOLECULAR,
                         state=Construct((role("tracrRNA and LbCas12a direct-repeat donor interval"),)))
    donor_waste = g.product("donor_waste", "unretained TOPO donor classes", K.WASTE)
    g.react("release_scaffold", "Release the separate scaffold donor in assembly round two", S.DIGESTION,
            [(donor, R.TEMPLATE), (bsm, R.CATALYST)], [released, donor_waste])
    g.discard(donor_waste, "Source donor backbone is distinct from the retained scaffold interval; exact map and fragment fates unresolved.")
    final = completed_construct()
    assembled = g.product("completed", "INFERRED donor junctions — " + final.name, K.MOLECULAR,
                          state=closed_circle(final.name, tuple(final)))
    g.react("second_assembly", "INFERRED donor junctions — introduce the scaffold in the second assembly", S.LIGATION,
            [(cut_ready, R.TEMPLATE), (released, R.STRUCTURAL), (ligase, R.CATALYST)], [assembled],
            geometry=UnassignedAssembly("The source names the scaffold donor. Recipient-side stuffer ends are computable; matching donor-side junctions remain inferred pending its original map.", fragment))
    second_clean = select("second_cleanup", assembled, cleanup, assembled.state, label="INFERRED donor junctions — purify second-assembly products")
    plasmids = select("plasmids", second_clean, antibiotic, assembled.state, label="INFERRED donor junctions — recover the completed plasmid class", kind=S.SELECTION)
    sample = g.input("sample", "source downstream sample context; biological state unassigned", K.BIOLOGICAL, primary=True)
    genomic = g.product("genomic", "post-assay genomic DNA; integration/genotype unassigned", K.MOLECULAR,
                        state=Construct((role("source recovered genomic DNA"),)))
    g.react("assay_boundary", "Unexpanded downstream assay and genomic recovery boundary", S.MAGIC,
            [(plasmids, R.TEMPLATE), (sample, R.TEMPLATE)], [genomic],
            reason="Downstream delivery, integration, guide processing, cellular assay and genomic recovery are unmodelled. No cellular genotype, editing or quantitative coverage follows from this graph.")
    of = g.input("outer_f", "PCR1_Hybrid_Outer_F2", K.OLIGO, state=Construct((bp("outer forward", OUTER_F),)))
    orev = g.input("outer_r", "PCR1_Hybrid_Outer_R1", K.OLIGO, state=Construct((bp("outer reverse", OUTER_R),)))
    outer = Construct((bp("outer forward-written end", OUTER_F), role("unassigned outer-to-inner interval"),
        bp("U6-side inner site", F_SITE), role("unassigned pre-spacer vector interval"), *completed_cassette(),
        role("unassigned post-spacer vector interval"), bp("termination-side inner site", revcomp(R_SITE)),
        role("unassigned inner-to-outer interval"), bp("outer reverse-written end", revcomp(OUTER_R))), "intended outer-PCR cassette family")
    first_pcr = g.product("outer_amplicon", outer.name, K.MOLECULAR, state=outer)
    g.react("outer_pcr", "INFERRED genomic context — amplify the cassette with outer primers", S.PCR,
            [(genomic, R.TEMPLATE), (of, R.STRUCTURAL), (orev, R.STRUCTURAL), (ngs_q5, R.CATALYST)], [first_pcr],
            geometry=PCRGeometry.unavailable(genomic.state, reason="The source outer primers are exact, but full genomic sites and outer/inner intervals are unassigned. Intended cassette context is shown without a genomic-coordinate claim."))
    fparts, rparts = index_primers()
    f = g.input("index_f", "A502-Cas12a-F1 (2020 A502-Cpf1-F1)", K.OLIGO, state=Construct(fparts))
    r = g.input("index_r", "A701-Cas9-R1", K.OLIGO, state=Construct(rparts))
    lib = library()
    pcr2 = g.product("indexed", lib.name, K.LIBRARY, state=lib)
    g.react("index_pcr", "INFERRED cassette context — add the illustrative source i5/i7 pair", S.PCR,
            [(first_pcr, R.TEMPLATE), (f, R.STRUCTURAL), (r, R.STRUCTURAL), (ngs_q5, R.CATALYST)], [pcr2],
            geometry=PCRGeometry.exact(outer, PCRPrimer("A502-Cas12a-F1", f.state.top(), len(F_SITE)),
                PCRPrimer("A701-Cas9-R1", r.state.top(), len(R_SITE)), statement="The source annealing ends are represented on an intended outer-PCR family. This establishes primer/declared-site consistency, not genomic mapping or source read-side assignment."))
    selected = select("sequencing_library", pcr2, gel, lib, label="INFERRED cassette context — select the indexed sequencing-library band", kind=S.SELECTION, material_kind=K.LIBRARY)
    g.sequence(selected, SequencingPlan(library=lib, primers=tuple(sp.TRUSEQ.values()),
        acquired_roles=tuple(p.role for p in sp.TRUSEQ.values()),
        limitation="2021 Step139 gives dark29/Read1 27/Index1 8/Index2 8/dark20/Read2 27, and names Read1 Cas12a/Read2 Cas9. That assignment has not been reconciled with the Step34 oligo order and Table1 primer orientations. This model shows potential canonical primer sites without guide-cycle spans or automatic dark-cycle coordinates. Recipient/donor/genomic maps and interval lengths are unassigned. The displayed A502/A701 combination is an illustrative pair, not an experimental sample assignment."))
    g.declare_reviewed_scope(SOURCE_INVENTORY, SOURCE)
    return g.require_publishable()


SOURCE_INVENTORY = (
    "pool", "tko_f", "tko_r", "q5", "ngs_q5", "cleanup", "gel", "antibiotic", "bfuai", "bsmbi", "ligase", "sap", "recipient",
    "open_recipient", "recipient_dephosphorylate", "recipient_ready_step", "amplify_pool", "pool_clean_step", "cut_pool",
    "first_assembly", "first_cleanup_step", "first_pool_step", "second_cut", "intermediate_dephosphorylate",
    "intermediate_ready_step", "scaffold_donor", "release_scaffold", "second_assembly", "second_cleanup_step", "plasmids_step",
    "sample", "assay_boundary", "outer_f", "outer_r", "outer_pcr", "index_f", "index_r", "index_pcr",
    "sequencing_library_step", "sequencing_library",
)


def render_page():
    graph = protocol_graph(); f, r = index_primers()
    return assemble_page(graph, (
        head(TITLE), '<div class="wrap" data-protocol-complete="true">', f"<h1>{TITLE}</h1>",
        '<p class="research-notes"><a href="01_chymera.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1038/s41587-020-0437-z">Defining method: Gonatopoulos-Pournatzis et al. (2020)</a>. '
             '<a href="https://doi.org/10.1038/s41596-021-00595-1">Recipe modelled here: Aregger et al., detailed protocol (2021)</a>.'),
        info("This graph models the 2021 pooled two-round construction and nested-PCR readout. The 2020 Table14 independently matches the named primer rows; the now-reviewed 2020 main Methods independently describe two rounds but do not establish an identical recipe. In particular the2020 read lengths differ from2021; the versions remain separate. Focused-cloning alternatives and higher-order arrays are outside this graph. Recipient, scaffold-donor and genomic maps remain unassigned."),
        info("INFERRED: intended scaffold-side junctions and downstream cassette context are marked. Spacer windows retain their distinct Cas9 and Cas12a identities; no target is selected. Downstream delivery and the cellular assay are an explicit unmodelled boundary. No sequence check establishes genotype, editing, functional effect, balance or experimental coverage."),
        WorkflowSection("Two assemblies and downstream readout", paths=(("recipient_ready", ()),
            ("first_insert", ()), ("scaffold_insert", ()),
            ("plasmids", ("recipient_ready_step", "cut_pool", "release_scaffold")),
            ("sequencing_library", ("plasmids_step",)))),
        OligoSection("Source primers used in this example", (
            OligoEntry("TKO_F1", (bp("pool forward", TKO_F),)), OligoEntry("TKO_R1", (bp("pool reverse", TKO_R),)),
            OligoEntry("PCR1_Hybrid_Outer_F2", (bp("outer forward", OUTER_F),)),
            OligoEntry("PCR1_Hybrid_Outer_R1", (bp("outer reverse", OUTER_R),)),
            OligoEntry("A502-Cas12a-F1 / 2020 A502-Cpf1-F1", f), OligoEntry("A701-Cas9-R1", r)),
            note_html="The displayed source-listed i5/i7 combination is illustrative. It does not assign those indexes to an experimental sample; these F1/R1 examples carry no added stagger bases."),
        "<h2>Unresolved read-side assignment</h2>",
        info("2021 Step139 names Read1 as Cas12a and Read2 as Cas9. Step34 orders the encoded spacer oligo Cas9 then Cas12a, while the Table1 i5 primer binds the U6 side and i7 binds the termination side. These statements have not been reconciled against the complete original maps and Figure5 graphic. The model retains the source statements and source oligos, shows potential canonical primer sites, and assigns no guide-cycle spans or corrective recipe. Dark-cycle offsets are recorded below rather than treated as imaged sequence."),
        info("Version-specific acquisition: the2020 defining article Methods, Preparation of CRISPR sequencing libraries and Illumina sequencing, reports29 dark cycles then31 imaged bases (Cas12a), and20 dark cycles then30 imaged bases (Cas9), with two8-base indexes. The2021 detailed protocol Step139 reports27 imaged bases for each side. This page retains the2021 recipe; the earlier cycle counts are not silently substituted. Read-side orientation remains unresolved."),
        sp.read_spec_table((sp.read_spec("Read 1", 27, "source Step139: Cas12a; preceded by 29 dark cycles; side unresolved", "include"),
            sp.read_spec("Index 1 (i7)", 8, "i7 sample index", "demultiplex"),
            sp.read_spec("Index 2 (i5)", 8, "i5 sample index", "demultiplex"),
            sp.read_spec("Read 2", 27, "source Step139: Cas9; preceded by 20 dark cycles; side unresolved", "include"))),
        PrimerBindingSection("sequencing_library", "Potential sequencing-primer sites", "Declared source adapter arms",
            intro_html="Source adapter arms support these potential sites. No run-length or dark-cycle offset is applied to infer guide coverage.", show_read_content=False),
        "</div>",
    ))


def selftest():
    import csv
    from checks import Check, Source, have
    from chemistry_validation import Status, check_graph
    from plasmid import BsmBI, Plasmid, amplify, digest
    from storage import sources_root
    check = Check(); graph = protocol_graph(); amp = amplified_pool()
    check("source-scoped graph has no contradictory molecular claims", [c for c in check_graph(graph) if c.status is Status.CONTRADICTORY], [])
    independent = amplify(Plasmid("source synthetic pool", synthesis().top().upper(), circular=False), TKO_F, TKO_R, min_anneal=15)
    check("independent PCR reconstructs the declared pooled amplicon", [p.seq for p in independent], [amp.top().upper()])
    check("pooled amplicon agrees with source Step42 band length", len(amp.top()), 157)
    cuts = digest(Plasmid("known stuffer fixture", "N" * 30 + STUFFER.upper() + "N" * 30, circular=False), BsmBI)
    opened = opened_intermediate()
    check("source-written stuffer ends agree with independent BsmBI cuts", (opened.right, opened.left), tuple(c.overhang for c in cuts))
    source = Source(sources_root() / Path(__file__).resolve().parents[1].name /
        "s41587-020-0437-z_41587_2020_437_MOESM16_ESM.xlsx.txt", "Published 2020 Supplementary Table14; tools/get_sources.py and tools/doctext.py")
    if have(check, source, label="source primer identity across the published table"):
        rows = dict(r[:2] for r in csv.reader(source.path.read_text().splitlines(), delimiter="\t") if len(r)>1)
        f, r = index_primers()
        for name, bases in (("TKO_F1", TKO_F), ("TKO_R1", TKO_R), ("PCR1_Hybrid_Outer_F2", OUTER_F),
            ("PCR1_Hybrid_Outer_R1", OUTER_R), ("A502-Cpf1-F1", Construct(f).top()), ("A701-Cas9-R1", Construct(r).top())):
            check(name + " matches the independently published 2020 row", rows[name].strip(), bases)
    check.report()
