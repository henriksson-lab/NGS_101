"""Published Kir2.1 DIMPLE sublibraries and source-qualified sequencing scope."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Row, Segment, feature, revcomp
from cohesive_assembly import FivePrimeFragment
from page import head, info
from pcr import PCRGeometry
from protocol_graph import (InputRole as R, MaterialKind as K, MolecularRelation,
                            ProtocolGraph, SequencingPlan, StepKind as S)
from protocol_model import closed_circle, role
from protocol_views import OligoEntry, OligoSection, WorkflowSection, assemble_page
import seqprimers as sp

TITLE = "DIMPLE — Kir2.1 sublibrary assembly and sequencing"
PAGE_FILE = "dimple.html"
SOURCE = "Macdonald et al. 2023, doi:10.1186/s13059-023-02880-6; Methods library generation/cloning, sequencing, landing-pad construction; published Additional file3 FigureS2 and Additional file6 TableS3"
BSAI = "GGTCTC"
TARGET_F = "CCCTATCAGTGATAGAGAACGTATGACC"
TARGET_R = "TGGTCCTGGATTTTCCTCCACATC"

# Published Table S3: nine Kir2.1 sublibrary primer sets, source case retained.
# (sublibrary number, oligo F, oligo R, gene F, gene R)
SPECS = (
    (1, 'AAGACTCAACCAATGACCCTTC', 'CAATACTGAGCCGTGGGTT', 'ATAGGTCTCAAGATGGCCACTGTAATGTGCA', 'ATAGGTCTCGTGAGACGGCcggatctg'),
    (2, 'ACAATAACTATGGGTGCGGGG', 'CCCGATGGATCTTTATACTCGG', 'ATAGGTCTCGCTGCGTTTTTTGGTTGATCGC', 'ATAGGTCTCAAGCGTGACCGGCACTG'),
    (3, 'TGGTTGTATAGCTGAGCGATG', 'GAGAGTTGCGCCGAGTTTG', 'ATAGGTCTCCCACTATCGGATACGGGTTC', 'ATAGGTCTCAGCCAAGACAGTACGAAGGC'),
    (4, 'AGCTAAGACAATACCCTCGGAGG', 'CCGTAAGCTCTAGAGGTTGTC', 'ATAGGTCTCATGAGACATTGGTTTTCAGTCACAACG', 'ATAGGTCTCTCGATAGAGAAAAGGAAAGCTGCG'),
    (5, 'TGGTGATAGGTAAGGATGGCA', 'CATAATATCTGACGTTAGGCGTTGG', 'ATAGGTCTCCAGAAGGAGAGTACATACCACTCG', 'ATAGGTCTCTTGGGCTTAGCCATTTTTGCC'),
    (6, 'GTCGATGGCTCCCCTTATTAC', 'GACAGTCTCTCGTTTGTTCCTTAGG', 'ATAGGTCTCACATCGACAACGCTGATTTCGAA', 'ATAGGTCTCCGTGATTTCAGCAGTTGTGCTCG'),
    (7, 'AGTTAGTGGTTGCCTCAGAATC', 'ATCGCATCCCATTATCTAGCG', 'ATAGGTCTCAGCACTATTATAAGGTCGACTACTCTAG', 'ATAGGTCTCGAGAGGTCATAAAGAGGAGAATCCTC'),
    (8, 'GCCTTTATGTTCCATGGTCCT', 'CGTCTAGCCCAGAGAGAGTC', 'ATAGGTCTCTCACATCAAAGGAAGAAGAAGAAGATT', 'ATAGGTCTCAACAAGACAGGTTCGTAGCGATG'),
    (9, 'GAGCGGATCTGGATACGTAATTGG', 'TTAGACTTGCAGTCTCATAATGGG', 'ATAGGTCTCCGAGACGtGCAACTAATTTTAGTCTAC', 'ATAGGTCTCTCGTTCTCGTAACAGAAAGAGTTTGC'),
)


def bp(name, sequence, **kw):
    return Segment(name, sequence, **kw)


def handle(number, side, sequence):
    return Segment(f"subpool{number} {side} amplification handle", sequence,
        feature=feature(f"dimple.subpool{number}.{side}", "feature_barcode", "fixed",
            group=f"subpool{number}", part=side,
            note="Sublibrary amplification handle outside the transferred variant-bearing region; not a sequencing sample index or UMI."))


def oligo_architecture(number, forward, reverse):
    return Construct((handle(number, "forward", forward), role("left variable length buffer"),
        bp("left BsaI recognition", BSAI), role("left source cut spacer and fusion unassigned"),
        role("left unmutated overlap", length_bp=4), role(f"subpool{number} variant-bearing gene region"),
        role("right unmutated overlap", length_bp=4), role("right source fusion and cut spacer unassigned"),
        bp("right BsaI recognition", revcomp(BSAI)), role("right variable length buffer"),
        handle(number, "reverse", revcomp(reverse))), f"subpool{number} oligo architecture")


def backbone_fragment(number, forward, reverse):
    """Selected BsaI ends written into the inverse-PCR primers, not a full map."""
    ends, tails = [], []
    for primer in (forward.upper(), reverse.upper()):
        start = primer.index(BSAI) + len(BSAI) + 1
        ends.append(primer[start:start + 4]); tails.append(primer[start + 4:])
    return FivePrimeFragment(f"subpool{number} backbone fragment", ends[0],
        (bp("forward primer-written gene flank", tails[0]), role("remaining Kir2.1 and source backbone"),
         bp("reverse primer-written gene flank", revcomp(tails[1]))), revcomp(ends[1]))


@dataclass(frozen=True)
class SourceAssemblyGeometry:
    label: str
    statement: str
    backbone: FivePrimeFragment | None = None
    availability = "unavailable"
    ready_label = "source-scoped assembly inputs"
    action_label = "Present source inputs without inventing unassigned ends"

    def rows(self):
        rows = self.backbone.rows() if self.backbone is not None else []
        return rows + [Row(chunks=[(self.label, None, False)])]


def expected_product(number, fragment):
    product = Construct((bp("INFERRED variant-side left junction", fragment.right, inferred=True),
        role("left retained unmutated overlap", length_bp=4), role(f"subpool{number} retained variant region"),
        role("right retained unmutated overlap", length_bp=4),
        bp("INFERRED variant-side right junction", fragment.left, inferred=True), *fragment.core),
        f"subpool{number} intended Kir2.1 construct")
    if any(seg.feature and seg.feature.id.startswith("dimple.subpool") for seg in product):
        raise ValueError("sublibrary handles and length buffers must not become retained identifier regions")
    return product


def pool_view(products):
    # Each independently drawn circle remains one molecular class in the pool.
    # No concatemer or combined gene sequence is constructed by pooling.
    from protocol_model import Drawing
    rows = []
    for product in products:
        rows.extend(closed_circle(product.name, tuple(product)).rows())
        rows.append(Row())
    return Drawing(tuple(rows))


def indexed_family():
    return Construct((role("kit-derived left library adapter; sequence unassigned"),
        role("sample i5", molecular_feature=feature("dimple.i5", "sample_index", "unknown",
             note="IDT for Illumina UD index; exact product/version, length and sample value unassigned.")),
        role("tagmented Kir2.1 variant-bearing fragment; endpoints unassigned"),
        role("sample i7", molecular_feature=feature("dimple.i7", "sample_index", "unknown",
             note="Partner UD index; no index length is inferred from the kit family name.")),
        role("kit-derived right library adapter; sequence unassigned")), "indexed fragment family")


def protocol_graph():
    if [row[0] for row in SPECS] != list(range(1, 10)):
        raise ValueError("published Kir2.1 inventory must contain its nine distinct sublibraries")
    g = ProtocolGraph(TITLE, authority=SOURCE)
    ols = g.input("ols", "source Kir2.1 oligo pool; concrete member sequences unassigned", K.OLIGO, primary=True)
    reference = g.input("reference", "source Kir2.1 template vector; full sequence unassigned", K.MOLECULAR,
                        state=Construct((role("Kir2.1 source vector"),)))
    polymerase = g.input("primestar", "PrimeStar GXL", K.CATALYST)
    bsai = g.input("bsai", "BsaI-HF v2 Golden Gate mixture", K.CATALYST)
    ligase = g.input("ligation", "source assembly ligation component; identity unassigned", K.CATALYST)
    cleanup = g.input("cleanup", "source Zymo DNA cleanup", K.SELECTION)
    gel = g.input("gel", "source gel isolation", K.SELECTION)
    kan = g.input("kanamycin", "source kanamycin clone selection", K.SELECTION)
    selected, products = [], []
    for number, of, orev, gf, gr in SPECS:
        prefix = f"s{number}"
        oligo = oligo_architecture(number, of, orev)
        of_node = g.input(prefix + "_of", prefix + " oligo forward", K.OLIGO, state=Construct((handle(number, "forward", of),)))
        or_node = g.input(prefix + "_or", prefix + " oligo reverse", K.OLIGO, state=Construct((handle(number, "reverse", orev),)))
        gf_node = g.input(prefix + "_gf", prefix + " inverse-PCR forward", K.OLIGO, state=Construct((bp("gene forward", gf),)))
        gr_node = g.input(prefix + "_gr", prefix + " inverse-PCR reverse", K.OLIGO, state=Construct((bp("gene reverse", gr),)))
        amplified = g.product(prefix + "_oligo_amp", oligo.name, K.MOLECULAR, state=oligo)
        g.react(prefix + "_oligo_pcr", "Amplify " + prefix + " with its sublibrary handles", S.PCR,
                [(ols, R.TEMPLATE), (of_node, R.STRUCTURAL), (or_node, R.STRUCTURAL), (polymerase, R.CATALYST)],
                [amplified], geometry=PCRGeometry.unavailable(oligo,
                    reason="Table S3 assigns the handles; concrete Kir2.1 oligo members, buffering lengths and target-site coordinates are unassigned. Figure S2 supplies the generic architecture only."))
        clean = g.product(prefix + "_oligo_clean", prefix + " purified oligo amplicon", K.MOLECULAR, state=oligo)
        waste = g.product(prefix + "_oligo_waste", prefix + " unretained oligo-PCR material", K.WASTE)
        g.react(prefix + "_oligo_cleanup", "Purify " + prefix + " oligo PCR", S.CLEANUP,
                [(amplified, R.TEMPLATE), (cleanup, R.SELECTION)], [clean, waste], changes_molecule=False,
                relations=(MolecularRelation.STRAND_ORDER,))
        g.discard(waste, "Source cleanup rejects material; no recovery fraction is inferred.")
        backbone_amp = Construct((bp("inverse-PCR forward", gf), role("retained reference gene/vector interval"),
                                  bp("inverse-PCR reverse complement", revcomp(gr))), prefix + " inverse-PCR product")
        backbone = g.product(prefix + "_gene_amp", backbone_amp.name, K.MOLECULAR, state=backbone_amp)
        g.react(prefix + "_gene_pcr", "Make the matching " + prefix + " backbone by inverse PCR", S.PCR,
                [(reference, R.TEMPLATE), (gf_node, R.STRUCTURAL), (gr_node, R.STRUCTURAL), (polymerase, R.CATALYST)],
                [backbone], geometry=PCRGeometry.unavailable(backbone_amp,
                    reason="Primer-written terminal sequences are known; full reference-vector identity and primer landing coordinates remain unassigned."))
        gene_clean = g.product(prefix + "_gene_clean", prefix + " gel-purified backbone PCR", K.MOLECULAR, state=backbone_amp)
        gene_waste = g.product(prefix + "_gene_waste", prefix + " unselected backbone-PCR material", K.WASTE)
        g.react(prefix + "_gene_cleanup", "Gel-isolate " + prefix + " backbone PCR", S.SELECTION,
                [(backbone, R.TEMPLATE), (gel, R.SELECTION)], [gene_clean, gene_waste], changes_molecule=False,
                relations=(MolecularRelation.STRAND_ORDER,))
        g.discard(gene_waste, "Material outside the selected backbone band; no efficiency is inferred.")
        fragment = backbone_fragment(number, gf, gr)
        released_backbone = g.product(prefix + "_backbone", fragment.name, K.MOLECULAR, state=fragment, chemical_state=fragment)
        insert = g.product(prefix + "_insert", prefix + " selected variant-bearing insert, ends unassigned", K.MOLECULAR,
                           state=Construct((role("selected variant-bearing region; terminal sequences unassigned"),)))
        terminals = g.product(prefix + "_terminal_waste", prefix + " discarded handles/buffers/terminal fragments", K.WASTE)
        g.react(prefix + "_cut", "Release the paired " + prefix + " BsaI assembly substrates", S.DIGESTION,
                [(clean, R.TEMPLATE), (gene_clean, R.TEMPLATE), (bsai, R.CATALYST)], [released_backbone, insert, terminals])
        g.discard(terminals, "Amplification handles, variable external buffers and recognition-site terminal material stay outside the intended transferred variant region; exact fragment sequences are unassigned.")
        product = expected_product(number, fragment)
        joined = g.product(prefix + "_joined", "INFERRED insert junctions — " + product.name, K.MOLECULAR,
                           state=closed_circle(product.name, tuple(product)), chemical_state=product)
        g.react(prefix + "_join", "INFERRED insert junctions — assemble " + prefix, S.LIGATION,
                [(released_backbone, R.TEMPLATE), (insert, R.STRUCTURAL), (ligase, R.CATALYST)], [joined],
                geometry=SourceAssemblyGeometry("Variant-bearing oligo insert: actual terminal sequence unassigned",
                    "Table S3 defines backbone primer-written ends. Compatible insert junctions are intended/inferred until actual Kir2.1 oligos and the full template map are assigned; no complete assembly-product relation is claimed.", fragment))
        retained = g.product(prefix + "_selected", "INFERRED insert junctions — selected " + prefix + " plasmid class", K.MOLECULAR,
                             state=joined.state, chemical_state=product)
        rejected = g.product(prefix + "_clone_waste", prefix + " unselected clone material", K.WASTE)
        g.react(prefix + "_recover", "INFERRED insert junctions — recover " + prefix + " selected plasmids", S.SELECTION,
                [(joined, R.TEMPLATE), (kan, R.SELECTION)], [retained, rejected], changes_molecule=False)
        g.discard(rejected, "Represents source bacterial selection/recovery without proving genotype, population growth, coverage or abundance balance.")
        selected.append(retained); products.append(product)
    pooled = g.product("sublibrary_pool", "pool of nine distinct Kir2.1 sublibrary classes", K.MOLECULAR, state=pool_view(products))
    g.react("pool_sublibraries", "INFERRED insert junctions — pool separately constructed sublibraries", S.POOL,
            [(p, R.TEMPLATE) for p in selected], [pooled], changes_molecule=False)
    bsm = g.input("bsmbi", "source BsmBI transfer", K.CATALYST)
    recipient = g.input("landing_recipient", "source landing-pad-compatible recipient; sequence/ends unassigned", K.MOLECULAR,
                         state=Construct((role("source landing-pad recipient"),)))
    transfer = g.product("transfer_donor", "BsmBI-released donor classes, exact ends unassigned", K.MOLECULAR,
                         state=Construct((role("source variant-bearing transfer donor"),)))
    opened_recipient = g.product("transfer_recipient", "BsmBI-opened recipient, exact ends unassigned", K.MOLECULAR,
                         state=Construct((role("source opened transfer recipient"),)))
    transfer_waste = g.product("transfer_waste", "unretained transfer-backbone/cassette material", K.WASTE)
    g.react("transfer_cut", "Prepare the separate BsmBI landing-pad transfer", S.DIGESTION,
            [(pooled, R.TEMPLATE), (recipient, R.TEMPLATE), (bsm, R.CATALYST)], [transfer, opened_recipient, transfer_waste])
    g.discard(transfer_waste, "Unretained transfer material; target-specific Kir2.1 end sequences are not borrowed from OPRM1/TRPV1.")
    transferred = g.product("transfer_construct", "intended landing-pad-compatible variant pool", K.MOLECULAR,
                            state=Construct((role("intended Kir2.1 variant-bearing constructs in source transfer backbone"),)))
    g.react("transfer_join", "Assemble the source landing-pad-compatible pool", S.LIGATION,
            [(transfer, R.TEMPLATE), (opened_recipient, R.STRUCTURAL), (ligase, R.CATALYST)], [transferred],
            geometry=SourceAssemblyGeometry("Kir2.1 transfer ends and complete donor/recipient maps unassigned",
                "The published source names a separate BsmBI transfer. CATG/GGGT values explicitly assigned to OPRM1/TRPV1 are not assigned to Kir2.1 by similarity."))
    samples = g.input("samples", "source post-assay sample context; biological state unassigned", K.BIOLOGICAL, primary=True)
    genomic = g.product("genomic", "recovered genomic DNA; genotype unassigned", K.MOLECULAR,
                        state=Construct((role("source genomic DNA"),)))
    g.react("assay_boundary", "Unexpanded cell assay and genomic recovery boundary", S.MAGIC,
            [(transferred, R.TEMPLATE), (samples, R.TEMPLATE)], [genomic],
            reason="The source downstream assay/sorting and DNA extraction are an explicit unmodelled boundary. No integration, cellular selection mechanism, genotype, functional effect or quantitative recovery is established by this graph.")
    tf = g.input("target_f", "cell_line_for_3", K.OLIGO, state=Construct((bp("target forward", TARGET_F),)))
    tr = g.input("target_r", "P2A_cell_line_rev", K.OLIGO, state=Construct((bp("target reverse", TARGET_R),)))
    target = Construct((bp("forward primer-written end", TARGET_F), role("enriched Kir2.1 variant-bearing region"),
                        bp("reverse primer-written end", revcomp(TARGET_R))), "target-enriched amplicon family")
    enriched = g.product("target_amplicon", target.name, K.MOLECULAR, state=target)
    g.react("enrich_target", "Amplify the source Kir2.1 genomic target", S.PCR,
            [(genomic, R.TEMPLATE), (tf, R.STRUCTURAL), (tr, R.STRUCTURAL), (polymerase, R.CATALYST)], [enriched],
            geometry=PCRGeometry.unavailable(genomic.state,
                reason="Written target primers are source-exact; their complete genomic coordinates and recovered target genotypes remain unassigned."))
    clean = g.product("target_clean", "gel-isolated target amplicons", K.MOLECULAR, state=target)
    target_waste = g.product("target_waste", "unselected target-PCR material", K.WASTE)
    g.react("target_select", "Gel-isolate the enriched target", S.SELECTION,
            [(enriched, R.TEMPLATE), (gel, R.SELECTION)], [clean, target_waste], changes_molecule=False,
            relations=(MolecularRelation.STRAND_ORDER,))
    g.discard(target_waste, "Physical gel rejection is distinct from later computational rejection of unexpected sequence variants.")
    tn5 = g.input("nextera_xt", "source Nextera XT library-preparation component", K.CATALYST)
    adapters = g.input("kit_adapters", "source kit adapters; exact sequences/version unassigned", K.ADAPTER)
    tagged = g.product("tagged", "tagmented fragment population; endpoints/adapters unassigned", K.MOLECULAR,
                       state=Construct((role("tagged variant-bearing fragment classes"),)))
    g.react("tagment", "Prepare Nextera XT fragments from enriched amplicons", S.TAGMENTATION,
            [(clean, R.TEMPLATE), (adapters, R.STRUCTURAL), (tn5, R.CATALYST)], [tagged])
    i5 = g.input("index_i5", "IDT for Illumina UD i5 reagent; exact product unassigned", K.OLIGO,
                 state=Construct((role("i5 index reagent sequence/length unassigned"),)))
    i7 = g.input("index_i7", "IDT for Illumina UD i7 reagent; exact product unassigned", K.OLIGO,
                 state=Construct((role("i7 index reagent sequence/length unassigned"),)))
    kit_pcr = g.input("kit_pcr", "source sequencing-library PCR component", K.CATALYST)
    indexed = g.product("indexed", "source-indexed fragment family", K.LIBRARY, state=indexed_family())
    g.react("index_library", "Generate the source UD-indexed sequencing library", S.PCR,
            [(tagged, R.TEMPLATE), (i5, R.STRUCTURAL), (i7, R.STRUCTURAL), (kit_pcr, R.CATALYST)], [indexed],
            geometry=PCRGeometry.unavailable(tagged.state,
                reason="The source names Nextera XT and IDT for Illumina UD indexes, without an assigned exact kit/index product in this model. Adapter sequences, index lengths, PCR sites and fragment endpoints are not filled from another kit recipe."))
    spri = g.input("spriselect", "source SPRISelect library selection", K.SELECTION)
    ready = g.product("sequencing_library", "selected indexed Kir2.1 fragment family", K.LIBRARY, state=indexed.state)
    removed = g.product("library_waste", "unretained sequencing-library material", K.WASTE)
    g.react("select_library", "Select the indexed sequencing-library products", S.CLEANUP,
            [(indexed, R.TEMPLATE), (spri, R.SELECTION)], [ready, removed], changes_molecule=False,
            relations=(MolecularRelation.STRAND_ORDER,))
    g.discard(removed, "Unretained cleanup classes; fragment-size distributions, end-combination fates and coverage are not assigned from kit family names.")
    g.sequence(ready, SequencingPlan(acquired_roles=("Read 1", "Read 2"),
        limitation="Source NovaSeq6000 SP300 paired-end acquisition is recorded. SP300 capacity is not converted into per-read cycles. Exact Nextera XT/IDT UD product versions, molecular sequences, index lengths, sample values, fragment endpoints and sequencing-primer geometry remain unassigned; no per-variant coverage or fitness follows from these declarations."))
    # Independent, curator-written source inventory, populated as literal IDs below.
    g.declare_reviewed_scope(SOURCE_INVENTORY, SOURCE)
    return g.require_publishable()


SOURCE_INVENTORY = (
    's1_of',
    's1_or',
    's1_gf',
    's1_gr',
    's1_oligo_pcr',
    's1_oligo_cleanup',
    's1_gene_pcr',
    's1_gene_cleanup',
    's1_cut',
    's1_join',
    's1_recover',
    's2_of',
    's2_or',
    's2_gf',
    's2_gr',
    's2_oligo_pcr',
    's2_oligo_cleanup',
    's2_gene_pcr',
    's2_gene_cleanup',
    's2_cut',
    's2_join',
    's2_recover',
    's3_of',
    's3_or',
    's3_gf',
    's3_gr',
    's3_oligo_pcr',
    's3_oligo_cleanup',
    's3_gene_pcr',
    's3_gene_cleanup',
    's3_cut',
    's3_join',
    's3_recover',
    's4_of',
    's4_or',
    's4_gf',
    's4_gr',
    's4_oligo_pcr',
    's4_oligo_cleanup',
    's4_gene_pcr',
    's4_gene_cleanup',
    's4_cut',
    's4_join',
    's4_recover',
    's5_of',
    's5_or',
    's5_gf',
    's5_gr',
    's5_oligo_pcr',
    's5_oligo_cleanup',
    's5_gene_pcr',
    's5_gene_cleanup',
    's5_cut',
    's5_join',
    's5_recover',
    's6_of',
    's6_or',
    's6_gf',
    's6_gr',
    's6_oligo_pcr',
    's6_oligo_cleanup',
    's6_gene_pcr',
    's6_gene_cleanup',
    's6_cut',
    's6_join',
    's6_recover',
    's7_of',
    's7_or',
    's7_gf',
    's7_gr',
    's7_oligo_pcr',
    's7_oligo_cleanup',
    's7_gene_pcr',
    's7_gene_cleanup',
    's7_cut',
    's7_join',
    's7_recover',
    's8_of',
    's8_or',
    's8_gf',
    's8_gr',
    's8_oligo_pcr',
    's8_oligo_cleanup',
    's8_gene_pcr',
    's8_gene_cleanup',
    's8_cut',
    's8_join',
    's8_recover',
    's9_of',
    's9_or',
    's9_gf',
    's9_gr',
    's9_oligo_pcr',
    's9_oligo_cleanup',
    's9_gene_pcr',
    's9_gene_cleanup',
    's9_cut',
    's9_join',
    's9_recover',
    'ols',
    'reference',
    'primestar',
    'bsai',
    'ligation',
    'cleanup',
    'gel',
    'kanamycin',
    'pool_sublibraries',
    'bsmbi',
    'landing_recipient',
    'transfer_cut',
    'transfer_join',
    'samples',
    'assay_boundary',
    'target_f',
    'target_r',
    'enrich_target',
    'target_select',
    'nextera_xt',
    'kit_adapters',
    'tagment',
    'index_i5',
    'index_i7',
    'kit_pcr',
    'index_library',
    'spriselect',
    'select_library',
    'sequencing_library',
)


def render_page():
    graph = protocol_graph()
    paths = []
    for number in range(1, 10):
        paths.extend(((f"s{number}_gene_clean", ()), (f"s{number}_selected", (f"s{number}_gene_cleanup",))))
    paths.append(("sequencing_library", tuple(f"s{i}_recover" for i in range(1, 10))))
    return assemble_page(graph, (
        head(TITLE), '<div class="wrap" data-protocol-complete="true">', f"<h1>{TITLE}</h1>",
        '<p class="research-notes"><a href="01_dimple.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1186/s13059-023-02880-6">Macdonald et al., Genome Biology (2023)</a>'),
        info("Kir2.1 source branch, with all nine Table S3 primer sets. Concrete oligo variants, reference/transfer vector maps and some kit details remain unassigned. INFERRED: intended variant-side junctions are marked; actual insert ends are not inferred as source facts from a matching backbone. The cell-assay boundary is explicitly unmodelled. No sequence check establishes genotype, library balance, editing, protein function, coverage or fitness."),
        info("Amplification handles and variable length buffers stay outside the transferred variant region. BsaI sublibrary assembly, pooling and BsmBI landing-pad transfer are separate operations. Target-specific CATG/GGGT values reported for other genes are not assigned to Kir2.1. Nextera XT/IDT UD names alone do not assign exact index oligos, lengths or sequencing cycles."),
        WorkflowSection("Nine sublibraries and downstream readout", paths=tuple(paths)),
        OligoSection("Published Kir2.1 primers", tuple(OligoEntry(f"subpool{n} {label}", (bp(label, seq),))
            for n, of, orev, gf, gr in SPECS for label, seq in
            (("oligo forward", of), ("oligo reverse", orev), ("gene forward", gf), ("gene reverse", gr))) + (
                OligoEntry("cell_line_for_3", (bp("target forward", TARGET_F),)),
                OligoEntry("P2A_cell_line_rev", (bp("target reverse", TARGET_R),)))),
        "<h2>Readout</h2>", info("Paired-end sequence alignment identifies variants in the source analysis. Sublibrary amplification handles are distinct from sequencing sample indexes and are not retained as variant tags. Physical selection/recovery is separate from the source's computational filtering of undesigned variants. Exact kit/index primer details and programmed cycle counts remain unassigned in this review."),
        sp.read_spec_table([sp.read_spec("Read 1", None, "variant-bearing fragment sequence", "include"),
                            sp.read_spec("Read 2", None, "paired fragment sequence", "include")]),
        "</div>",
    ))


def selftest():
    import csv
    from checks import Check, Source, have
    from chemistry_validation import Status, check_graph
    from plasmid import BsaI, Plasmid, digest
    from storage import sources_root
    check = Check(); graph = protocol_graph()
    check("source-scoped Kir2.1 graph has no contradictory molecular claims",
          [c for c in check_graph(graph) if c.status is Status.CONTRADICTORY], [])
    for n, of, orev, gf, gr in SPECS:
        known = gf.upper() + "N" * 30 + revcomp(gr.upper())
        cuts = digest(Plasmid("synthetic primer-end fixture", known, circular=False), BsaI)
        fragment = backbone_fragment(n, gf, gr)
        check(f"subpool{n} known primer-written BsaI ends agree with independent cuts",
              (fragment.left, fragment.right), tuple(c.overhang for c in cuts))
    source = Source(sources_root() / Path(__file__).resolve().parents[1].name /
        "s13059-023-02880-6_supp_13059_2023_2880_MOESM6_ESM.xlsx.txt",
        "Published Additional file6 TableS3, doi:10.1186/s13059-023-02880-6; tools/get_sources.py and tools/doctext.py")
    if have(check, source, label="Kir2.1 source primer assignment"):
        rows = list(csv.reader(source.path.read_text().split("\n")[1:], delimiter="\t"))
        for n, of, orev, gf, gr in SPECS:
            for prefix, sequence in ((f"Kir21_oligoP_DMS-{n}_F",of),(f"Kir21_oligoP_DMS-{n}_R",orev),
                                     (f"Kir21_geneP_Mut-{n}_F",gf),(f"Kir21_geneP_Mut-{n}_R",gr)):
                original = [r[1] for r in rows if len(r)>1 and r[0].split(" ")[0] == prefix]
                check(prefix + " source assignment/case preserved", original, [sequence])
    check.report()
