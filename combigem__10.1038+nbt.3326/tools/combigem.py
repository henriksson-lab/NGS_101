"""CombiGEM 2015: declared insert/barcode order and barcode-library readout."""
from __future__ import annotations

from dataclasses import dataclass, replace
from collections import Counter
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from cohesive_assembly import CircularCohesiveAssembly, FivePrimeFragment
import illumina as il
from page import head, info, table
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole as R, MaterialKind as K,
                            MolecularRelation, ProtocolGraph, SequencingPlan,
                            StepKind as S, Terminal)
from protocol_model import closed_circle, role
from protocol_views import (OligoEntry, OligoSection, PlanLimitationSection,
                            WorkflowSection, assemble_page)
import seqprimers as sp

PAGE_FILE = "combigem.html"
TITLE = "CombiGEM — ordered inserts and concatenated barcodes"
SOURCE = "Wong et al. 2015, doi:10.1038/nbt.3326; Online Methods, barcoded library creation/pooled assembly/barcode sequencing; published Supplementary Tables 1 and 13"
BGLII = "AGATCT"
BAMHI = "GGATCC"
ECORI = "GAATTC"
MFEI = "CAATTG"
PRECURSOR_JUNCTION = BAMHI[0] + BGLII[1:]
BARCODE_JUNCTION = MFEI[0] + ECORI[1:]
SPACER = "GCAACG"  # Table 1 reverse primers, reverse-complemented; also PCR forward.
FORWARD = il.P5 + BAMHI + SPACER + ECORI
REVERSE_SITE = "GGTTGCGTCAGCAAACACAG"
REVERSE = il.P7 + "N" * 8 + REVERSE_SITE
BARCODE_PRIMER = il.P5[13:] + BAMHI + SPACER + ECORI
INDEX_PRIMER = "GTGGCGTGGTGTGCACTGTGTTTGCTGACGCAACC"
QPCR_FORWARD = "AATGATACGGCGACCACCGA"
QPCR_REVERSE = "CAAGCAGAAGACGGCATACGA"


@dataclass(frozen=True)
class Member:
    key: str
    precursor: str
    reverse_primer_barcode: str

    @property
    def barcode(self):
        # Table 1 writes the barcode on the reverse PCR primer. Its molecule's
        # reference strand therefore carries the reverse complement.
        return revcomp(self.reverse_primer_barcode)

    def identifier_segment(self):
        return Segment(self.key + " barcode", self.barcode, feature=feature(
            "member." + self.key, "feature_barcode", "fixed",
            note=f"Table 1 source association: {self.precursor}; not a sample index."))


# First three original Table 1 rows, checked in the rendered original p3.
# One illustrative member trajectory, not an exhaustive enumeration of the pool.
MEMBERS = (Member("a", "let-7a", "TCGCGTTT"),
           Member("b", "let-7b/miR-4763 cluster", "ATGACCTT"),
           Member("c", "let-7c", "CGTTCTGT"))


def bp(name, bases):
    return Segment(name, bases)


def donor_core(member):
    return (bp("BglII downstream remnant", BGLII[-1]), role(member.precursor),
            bp("active BamHI", BAMHI), bp("inter-site spacer", SPACER),
            bp("active EcoRI", ECORI), member.identifier_segment(),
            bp("MfeI upstream remnant", MFEI[0]))


def donor_plasmid(member):
    return Construct((bp("BglII", BGLII), role(member.precursor),
                      bp("active BamHI", BAMHI), bp("inter-site spacer", SPACER),
                      bp("active EcoRI", ECORI), member.identifier_segment(),
                      bp("MfeI", MFEI), role("pBT264 storage backbone")),
                     "stored " + member.precursor)


def active_interval(product):
    """Selected known cassette, not a search for sites in unknown part bodies."""
    left = [i for i, seg in enumerate(product) if seg.top == BAMHI]
    right = [i for i, seg in enumerate(product) if seg.top == ECORI]
    if len(left) != 1 or len(right) != 1 or left[0] >= right[0]:
        raise ValueError("the declared active BamHI/EcoRI cassette must be unique and ordered")
    return left[0], right[0]


def retained_destination(product, number):
    left, right = active_interval(product)
    # EcoRI cuts G|AATTC, BamHI cuts G|GATCC. Open the retained circular
    # fragment at EcoRI; keep its final C and BamHI's first G in the core.
    core = (bp(f"round {number} EcoRI downstream remnant", ECORI[-1]),
            *product.segments[right + 1:], *product.segments[:left],
            bp(f"round {number} BamHI upstream remnant", BAMHI[0]))
    return FivePrimeFragment(f"round {number} retained destination", "AATT", core, "GATC")


def barcode_segments(product):
    _, right = active_interval(product)
    parts = []
    for seg in product.segments[right + 1:] + product.segments[:right + 1]:
        if seg.is_role_token():
            break  # unassigned backbone is outside the concatenated barcode tract
        parts.append(seg)
    return tuple(parts)


def declared_orders(product, members):
    """Check source-assigned identities in this declared product, not clone genotype.

    Barcode identity is structured data. Precursor cores remain source-labelled roles;
    matching these roles checks declaration bookkeeping, not their unknown DNA sequence.
    """
    labels = {"[" + m.precursor + "]": m.key for m in members}
    precursors = [labels[s.top] for s in product if s.top in labels]
    barcodes = [s.feature.id.removeprefix("member.") for s in barcode_segments(product)
                if s.feature and s.feature.id.startswith("member.")]
    if precursors != [m.key for m in members] or barcodes != list(reversed(precursors)):
        raise ValueError("declared precursor/barcode association order does not agree")
    multiplicity = Counter(m.key for m in members)
    for m in members:
        matching = [s for s in product if s.feature and s.feature.id == "member." + m.key]
        if len(matching) != multiplicity[m.key] or any(s.top != m.barcode for s in matching):
            raise ValueError("source-assigned barcode identity was changed or duplicated")
    return tuple(precursors), tuple(barcodes)


def sample_index(sample):
    return Segment(sample + " sample index", "N" * 8, placeholder=True,
                   feature=feature("sample." + sample, "sample_index", "unknown",
                                   note="Source assigns an eight-base sample index; no particular sample value is selected here."))


def reverse_oligo(sample):
    return (bp("P7", il.P7), sample_index(sample), bp("reverse barcode-PCR site", REVERSE_SITE))


def amplicon(sample, barcodes):
    # End sequence is supplied by the written PCR primers. The internal vector
    # interval is deliberately not filled from an unverified destination map.
    index = sample_index(sample)
    product = Construct((bp("P5", il.P5), bp("BamHI and spacer", BAMHI + SPACER),
                      bp("EcoRI", ECORI), *barcodes,
                      role("unassigned internal vector interval"),
                      bp("reverse-primer-derived site", revcomp(REVERSE_SITE)),
                      replace(index, top=revcomp(index.top)), bp("P7 reverse complement", il.P7_RC)),
                     sample + " barcode amplicon")
    problems = sp.verify(product, PRIMERS, required_roles=("Read 1", "Index 1 (i7)"))
    if problems:
        raise ValueError("barcode amplicon primer geometry: " + "; ".join(problems))
    return product


def genomic_barcode_tract():
    parts = []
    for n in range(3):
        parts.extend((Segment(f"genomic member {n + 1}", "N" * 8, placeholder=True,
                              feature=feature(f"genomic.member{n + 1}", "feature_barcode", "whitelist",
                                  whitelist="Wong2015 published Supplementary Table1",
                                  note="Unassigned member of the source's 39-barcode set; no recovered genomic identity is asserted.")),
                      bp(f"genomic junction {n + 1}", BARCODE_JUNCTION)))
    return tuple(parts)


PRIMERS = (sp.custom("Read 1", "CombiGEM barcode sequencing", BARCODE_PRIMER, SOURCE),
           sp.custom("Index 1 (i7)", "CombiGEM index sequencing", INDEX_PRIMER, SOURCE,
               "The model contains the primer's 20-base 3′ site supplied by reverse PCR. Its 15-base 5′ extension lies in an unassigned interval; a full 35-base template match is not established.",
               min_anneal=len(REVERSE_SITE)))


def pool_scene(libraries):
    sc = Scene()
    for sample, product in libraries.items():
        top = tuple(product)
        bottom = complement_segments(top)
        sc.strand(sample, top, label=sample)
        sc.anneal(sample + " bottom", bottom, to=sample,
                  pair=(bottom[-1].name, top[0].name),
                  unpaired=tuple(s.name for s in bottom if s.is_role_token()))
    return sc


def _graph_for_members(members):
    if len(members) != 3 or any(member not in MEMBERS for member in members):
        raise ValueError("this source scope covers three rounds using the inspected Table 1 members")
    g = ProtocolGraph(TITLE, authority=SOURCE)
    enzymes = {key: g.input(key, label, K.CATALYST) for key, label in (
        ("bglii", "BglII"), ("mfei", "MfeI"), ("bamhi", "BamHI"), ("ecori", "EcoRI"))}
    ligase = g.input("ligase", "DNA ligase (pooled-reaction identity unassigned)", K.CATALYST)
    polymerase = g.input("polymerase", "Kapa HiFi barcode PCR", K.CATALYST)
    gel = g.input("gel", "source barcode-amplicon gel isolation", K.SELECTION)
    forward = g.input("forward", "barcode PCR forward", K.OLIGO,
                      state=Construct((bp("P5", il.P5), bp("barcode-PCR binding end", BAMHI + SPACER + ECORI))))
    qforward = g.input("qforward", "quantification forward", K.OLIGO,
                       state=Construct((bp("quantification forward", QPCR_FORWARD),)))
    qreverse = g.input("qreverse", "quantification reverse", K.OLIGO,
                       state=Construct((bp("quantification reverse", QPCR_REVERSE),)))
    initial = Construct((bp("active BamHI", BAMHI), role("initial excised interval"),
                         bp("active EcoRI", ECORI), role("unassigned destination backbone")),
                        "source destination cassette")
    current = g.input("destination", "destination plasmid: declared cassette only", K.MOLECULAR,
                      state=closed_circle(initial.name, tuple(initial)), chemical_state=initial)
    product = initial
    for number, member in enumerate(members, 1):
        prefix = f"r{number}"
        donor = donor_plasmid(member)
        stored = g.input(prefix + "_stored", "stored source member: " + member.precursor,
                         K.MOLECULAR, state=closed_circle(donor.name, tuple(donor)),
                         chemical_state=donor, primary=True)
        insert_geometry = FivePrimeFragment(member.precursor + " and barcode", "GATC", donor_core(member), "AATT")
        insert = g.product(prefix + "_insert", insert_geometry.name, K.MOLECULAR,
                           state=insert_geometry, chemical_state=insert_geometry)
        donor_waste = g.product(prefix + "_donor_waste", "storage-vector backbone", K.WASTE,
                                state=Construct((role("unassigned storage-vector backbone"),)))
        g.react(prefix + "_release", "Release the stored insert/barcode member", S.DIGESTION,
                [(stored, R.TEMPLATE), (enzymes["bglii"], R.CATALYST), (enzymes["mfei"], R.CATALYST)],
                [insert, donor_waste])
        g.discard(donor_waste, "Storage-vector backbone is outside the intended assembly product; removal efficiency is not inferred.")
        vector_geometry = retained_destination(product, number)
        vector = g.product(prefix + "_vector", vector_geometry.name, K.MOLECULAR,
                           state=vector_geometry, chemical_state=vector_geometry)
        removed = g.product(prefix + "_excised", "excised BamHI–EcoRI cassette interval", K.WASTE,
                            state=Construct((role("excised selected-site interval"),)))
        g.react(prefix + "_open", "Open the destination at its selected BamHI/EcoRI cassette", S.DIGESTION,
                [(current, R.TEMPLATE), (enzymes["bamhi"], R.CATALYST), (enzymes["ecori"], R.CATALYST)],
                [vector, removed])
        g.discard(removed, "Selected-site excision product is not part of the intended retained circle.")
        geometry = CircularCohesiveAssembly(((vector.id, vector_geometry), (insert.id, insert_geometry)),
            "Online Methods: BamHI/BglII and EcoRI/MfeI compatible ends. Only the selected cassette sites and supplied fragment states are checked.")
        product = geometry.product(f"round {number} declared construct")
        declared_orders(product, members[:number])
        current = g.product(prefix + "_joined", product.name, K.MOLECULAR,
                            state=closed_circle(product.name, tuple(product)), chemical_state=product)
        g.react(prefix + "_join", f"Round {number}: seal the declared insert and retained vector", S.LIGATION,
                [(vector, R.TEMPLATE), (insert, R.STRUCTURAL), (ligase, R.CATALYST)], [current],
                geometry=geometry, relations=(MolecularRelation.COHESIVE_ASSEMBLY,))

    genomic = g.input("genomic_dna", "genomic DNA recovered after the reported assay", K.BIOLOGICAL,
                       state=Construct((role("unassigned genomic DNA"),)), primary=True)
    templates = {"plasmid": current, "genomic": genomic}
    libraries = {"plasmid": amplicon("plasmid", barcode_segments(product)),
                 "genomic": amplicon("genomic", genomic_barcode_tract())}
    cleaned = []
    for sample, template in templates.items():
        reverse = g.input(sample + "_reverse", sample + " barcode PCR reverse", K.OLIGO,
                           state=Construct(reverse_oligo(sample)))
        amplified = g.product(sample + "_amplicon", libraries[sample].name, K.LIBRARY, state=libraries[sample])
        geometry = PCRGeometry.unavailable(template.chemical_state or template.state,
            reason="Written primers and intended barcode target are known. Full source-vector/genomic primer-site coordinates are unassigned in this source scope; no complete template PCR relation is asserted.")
        g.react(sample + "_pcr", "Amplify " + sample + " barcode targets with sample-indexed anchors", S.PCR,
                [(template, R.TEMPLATE), (forward, R.STRUCTURAL), (reverse, R.STRUCTURAL),
                 (polymerase, R.CATALYST)], [amplified], geometry=geometry)
        clean = g.product(sample + "_clean", "gel-isolated " + sample + " barcode amplicon", K.LIBRARY,
                          state=libraries[sample])
        waste = g.product(sample + "_gel_waste", "unretained " + sample + " PCR material", K.WASTE)
        g.react(sample + "_isolate", "Isolate the source barcode-amplicon band", S.SELECTION,
                [(amplified, R.TEMPLATE), (gel, R.SELECTION)], [clean, waste],
                changes_molecule=False, relations=(MolecularRelation.STRAND_ORDER,))
        g.discard(waste, "Material outside the isolated target band; neither its identity nor its recovery fraction is inferred.")
        cleaned.append(clean)
    quantified = [g.product(sample + "_quantified", "quantified " + sample + " amplicon", K.LIBRARY,
                            state=libraries[sample]) for sample in templates]
    measurement = g.product("concentrations", "qPCR concentration estimates", K.QC)
    g.react("quantify", "Measure library concentrations for pooling", S.GENERIC,
            [(s, R.TEMPLATE) for s in cleaned] + [(qforward, R.STRUCTURAL), (qreverse, R.STRUCTURAL)],
            [*quantified, measurement], changes_molecule=False,
            reason="Measurement of aliquots informs pooling. The retained library sequences are unchanged; qPCR products and calibration chemistry are not modelled here.")
    g.mark_terminal(Terminal(measurement.id, Disposition.QC, "Concentration estimates do not establish barcode identity, balance or coverage."))
    pooled = g.product("indexed_pool", "pooled sample-indexed barcode libraries", K.LIBRARY,
                       state=pool_scene(libraries))
    g.react("pool", "Pool quantified samples for custom-primer HiSeq readout", S.POOL,
            [(s, R.TEMPLATE) for s in quantified], [pooled], changes_molecule=False)
    readout_primers = [g.input(identifier, label, K.OLIGO,
                              state=Construct((bp(label, sequence),)))
                      for identifier, label, sequence in (
                          ("barcode_read_primer", "custom barcode sequencing primer", BARCODE_PRIMER),
                          ("index_read_primer", "custom index sequencing primer", INDEX_PRIMER))]
    ready = g.product("sequencing_library", "barcode libraries for custom-primer HiSeq acquisition", K.LIBRARY,
                      state=pooled.state)
    g.react("readout_setup", "Use the source custom barcode and index sequencing primers", S.GENERIC,
            [(pooled, R.TEMPLATE)] + [(p, R.STRUCTURAL) for p in readout_primers], [ready],
            changes_molecule=False,
            reason="Acquisition setup records the source-declared readout oligos; no conversion of the retained template sequences or instrument-specific loading recipe is modelled.")
    g.sequence(ready, SequencingPlan(acquired_roles=("Read 1", "Index 1 (i7)"),
        limitation="The final pool retains distinct sample molecular classes. Primer geometry is computed separately for each class below. Source sample-index assignments, pooling proportions, complete internal vector intervals and acquisition cycle counts remain unassigned in this review. The reported 359-bp band is not used to invent read cycles or a common size for every assembly order."))
    g.declare_reviewed_scope((
        "bglii", "mfei", "bamhi", "ecori", "ligase", "polymerase", "gel", "forward", "qforward", "qreverse",
        "destination", "r1_stored", "r1_release", "r1_open", "r1_join",
        "r2_stored", "r2_release", "r2_open", "r2_join", "r3_stored", "r3_release", "r3_open", "r3_join",
        "genomic_dna", "plasmid_reverse", "genomic_reverse", "plasmid_pcr", "genomic_pcr",
        "plasmid_isolate", "genomic_isolate", "quantify", "pool", "barcode_read_primer",
        "index_read_primer", "readout_setup", "sequencing_library",
    ), SOURCE)
    return g.require_publishable()


def protocol_graph():
    return _graph_for_members(MEMBERS)


def render_page():
    graph = protocol_graph()
    order_rows = [(m.precursor, m.barcode, m.reverse_primer_barcode) for m in MEMBERS]
    content = [head(TITLE), '<div class="wrap" data-protocol-complete="true">',
        f"<h1>{TITLE}</h1>", '<p class="research-notes"><a href="01_combigem.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1038/nbt.3326">Wong et al., Nature Biotechnology (2015)</a>'),
        info("One illustrative member path through the source's pooled assembly: three source-assigned precursor/barcode pairs, followed by plasmid barcode QC. A separate genomic-DNA input represents recovery after the reported assay; delivery, integration and cellular phenotype are not reconstructed. Precursor bodies, full vector maps and some internal amplicon sequence remain unassigned. Checks establish declared cohesive-junction and barcode-order relations, not complete restriction digestion, clone genotype, association fidelity, library coverage or balance."),
        info("The Methods' shorthand BamHI–EcoRI interval includes GCAACG, verified in the published Table 1 reverse primers and barcode-PCR forward primer. Adding a precursor at the active cassette appends it to the precursor tract and prepends its barcode to the barcode tract. Thus the illustrative A→B→C precursor order has C→B→A barcode order, separated by derived CAATTC hybrid junctions. Genomic barcode identities remain unassigned; the diagram does not assert that these particular three precursors were recovered from cells."),
        table(("Illustrative member", "Construct-strand barcode", "Table 1 reverse-primer barcode"), order_rows),
        WorkflowSection("Declared assembly and barcode recovery", paths=(
            ("r1_insert", ()), ("r1_joined", ("r1_release",)),
            ("r2_insert", ()), ("r2_joined", ("r1_join", "r2_release")),
            ("r3_insert", ()), ("r3_joined", ("r2_join", "r3_release")),
            ("plasmid_clean", ("r3_join",)), ("genomic_clean", ()),
            ("sequencing_library", ("plasmid_isolate", "genomic_isolate")),
        )),
        OligoSection("Barcode readout oligos", (
            OligoEntry("barcode PCR forward", (bp("P5", il.P5), bp("barcode-PCR end", BAMHI + SPACER + ECORI))),
            OligoEntry("barcode PCR reverse", reverse_oligo("plasmid")),
            OligoEntry("custom barcode sequencing", (bp("barcode sequencing", BARCODE_PRIMER),)),
            OligoEntry("custom index sequencing", (bp("index sequencing", INDEX_PRIMER),)),
            OligoEntry("quantification forward", (bp("quantification forward", QPCR_FORWARD),)),
            OligoEntry("quantification reverse", (bp("quantification reverse", QPCR_REVERSE),)),
        )),
        "<h2>Sequencing readout</h2>",
        info("The two diagrams describe distinct sample classes before pooling. The custom index primer's 20-base 3′ match comes from the reverse PCR primer; its remaining 15 bases are not assigned a full template match. The index is separate from the precursor-identifying barcodes. Read cycles and actual sample indexes remain unassigned; the source's PCR-product size does not determine acquisition lengths."),
    ]
    for sample in ("plasmid", "genomic"):
        library = graph.materials[sample + "_clean"].state
        content.append(sp.section(library, PRIMERS, heading=sample + " sample primer sites",
            required_roles=("Read 1", "Index 1 (i7)"), read_lengths={},
            read_uses={"Read 1": "extract", "Index 1 (i7)": "demultiplex"}))
    content.extend((sp.read_spec_table([
        sp.read_spec("Read 1", None, "concatenated precursor-identifying barcodes", "extract"),
        sp.read_spec("Index 1 (i7)", None, "sample multiplexing index", "demultiplex")]),
        PlanLimitationSection("sequencing_library"), "</div>"))
    return assemble_page(graph, content)


def selftest():
    import re
    from checks import Check, Source, have
    from chemistry_validation import Status, check_graph
    from storage import sources_root
    check = Check()
    graph = protocol_graph()
    check("declared assembly and readout have no contradictory molecular claims",
          [c for c in check_graph(graph) if c.status is Status.CONTRADICTORY], [])
    final = graph.materials["r3_joined"].chemical_state
    # The source's iterative placement is the independent expectation here, not
    # a second transcription of the assembly-construction function.
    check("three rounds append precursors and prepend their associated barcodes",
          declared_orders(final, MEMBERS), (("a", "b", "c"), ("c", "b", "a")))
    check("ligation destroys the intended outer sites while retaining the next active cassette",
          (final.top().count(BGLII), final.top().count(MFEI), final.top().count(BAMHI), final.top().count(ECORI)),
          (0, 0, 1, 1))
    changed = Construct([replace(s, top="AAAAAAAA") if s.feature and s.feature.id == "member.b" else s for s in final])
    check.raises("changed barcode-to-member identity is rejected", lambda: declared_orders(changed, MEMBERS))
    repeated = (MEMBERS[0],) * 3
    repeat_product = _graph_for_members(repeated).materials["r3_joined"].chemical_state
    check("source Cartesian combinations allow repeated member identities",
          declared_orders(repeat_product, repeated), (("a", "a", "a"), ("a", "a", "a")))
    for sample in ("plasmid", "genomic"):
        product = graph.materials[sample + "_clean"].state
        hit = sp.locate(product, PRIMERS[1])
        check(sample + " index primer exposes exactly its unassigned 5-prime extension", hit.free5, 15)
    source = Source(sources_root() / Path(__file__).resolve().parents[1].name
                    / "nbt.3326_41587_2015_BFnbt3326_MOESM33_ESM.pdf.txt",
                    "Published Supplementary Table 1 for doi:10.1038/nbt.3326, original PDF p3; tools/get_sources.py and tools/doctext.py")
    if have(check, source, label="Table 1 barcode/spacer transcription"):
        text = source.path.read_text().split("Supplementary Table 1 List", 1)[1].split("Supplementary Table 2", 1)[0]
        joined = re.sub(r"\s+", "", text)
        for member in MEMBERS:
            prefix = "ATCG" + MFEI + member.reverse_primer_barcode + ECORI + revcomp(SPACER) + BAMHI
            check(member.precursor + " assigned reverse-primer prefix remains source-exact", prefix in joined)
    check.report()
