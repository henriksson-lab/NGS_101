"""Original STARR-seq published-method construction/readout graph; notes scope."""
from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Row, Scene, Segment, revcomp
from page import head, info
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import AnnotationSlot, InputRole as R, MaterialKind as K, ProtocolGraph, SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind as S
from protocol_model import Drawing, closed_circle, role
from protocol_views import OligoEntry, OligoSection, WorkflowSection, assemble_page
import seqprimers as sp
import illumina as il

PAGE_FILE = "starr-seq.html"
TITLE = "STARR-seq — OSC reporter library and paired readout"
SOURCE = "Arnold et al.2013 doi:10.1126/science.1232542; published article Fig1A and corrected Supporting Online Material (6 February2013), Methods pp5–8 and Mapping p11; original OSC deposits GSE40739/GSM1000403/GSM1000405 remain identity checks. No human or later STARR-seq recipe is substituted."

# Verbatim Drosophila primer bases, corrected published SOM pp6,8.
CLONING_F = "TAGAGCATGCACCGG" + il.TRUSEQ_READ1
CLONING_R = "GGCCGAATTCGTCGA" + revcomp(il.INDEX1_PRIMER)
RT_PRIMER = "CAAACTCATCAATGTATCTTATCATG"
REPORTER_F = "AAGCCACCATGGAAAAGGCCAT"
REPORTER_R = "TATCATGTCTGCTCGAAGCGG"
INPUT_F = "CTAGGCACACCGAAACGACTAAC"


def primer(name, sequence):
    return Construct((Segment(name, sequence),), name)


def source_oligos():
    return OligoSection("Published Drosophila primers", (
        OligoEntry("cloning forward", tuple(primer("cloning forward", CLONING_F))),
        OligoEntry("cloning reverse", tuple(primer("cloning reverse", CLONING_R))),
        OligoEntry("reporter RT", tuple(primer("reporter RT", RT_PRIMER))),
        OligoEntry("reporter first-PCR forward (3-prime terminal5nt phosphorothioate protected)", tuple(primer("reporter first-PCR forward", REPORTER_F))),
        OligoEntry("reporter first-PCR reverse (3-prime terminal3nt phosphorothioate protected); shared input reverse bases", tuple(primer("shared first-PCR reverse", REPORTER_R))),
        OligoEntry("input first-PCR forward", tuple(primer("input first-PCR forward", INPUT_F))),
    ), note_html="Corrected SOM pp6,8. Reporter first-PCR forward and reverse have phosphorothioate protection at the terminal five and three nucleotides respectively, as stated in the source; these modifications are not absent merely because the base strings above omit bond symbols. The input forward is a different primer; its protection is not assigned by analogy. E7335L second-PCR oligos and experimental index values remain unassigned.")



def unknown(name):
    return Construct((role(name),), name)


def recipient_readout_tail():
    """Source-qualified short pGL3 tail, not the complete author GenBank map.

    Corrected SOM p5 supplies the Drosophila cassette through its FseI boundary.
    The author-linked DSCP original matches that entire2502-base block byte-for-byte
    and its downstream tail contains both exact SOM p8 primer sites. Only that tail
    is used here. The inconsistent LOCUS length is never accepted as map identity.
    """
    reverse_site, rt_site = revcomp(REPORTER_R), revcomp(RT_PRIMER)
    overlaps = [n for n in range(1, min(len(reverse_site), len(rt_site)))
                if reverse_site.endswith(rt_site[:n])]
    if not overlaps:
        raise ValueError("source RT and first-PCR reverse contexts have no overlap")
    overlap = max(overlaps)
    # GG is the two-base source interval after the15nt right homology arm.
    return (Segment("source right-arm to readout interval", "GG", inferred=True),
            Segment("source shared first-PCR reverse site", reverse_site, inferred=True),
            Segment("source RT-site continuation", rt_site[overlap:], inferred=True))


def expected_reporter():
    """Partial intended DNA context; unresolved intervals are never fabricated bases."""
    return Construct((role("upstream DSCP context; extent unassigned", inferred=True),
        Segment("source input-forward site in intron", INPUT_F, inferred=True),
        role("remaining intron and sgGFP interval; length unassigned", inferred=True),
        Segment("source15nt left homology arm", CLONING_F[:15], inferred=True),
        Segment("cloning-forward written adapter context", CLONING_F[15:], inferred=True),
        role("retained genomic and adapter interior; end details unassigned", inferred=True),
        Segment("cloning-reverse written adapter context", revcomp(CLONING_R[15:]), inferred=True),
        Segment("source15nt right homology arm", revcomp(CLONING_R[:15]), inferred=True),
        *recipient_readout_tail(),
        role("remaining SV40 polyA and pGL3 backbone; map unassigned", inferred=True)),
        "INFERRED partial intended reporter DNA context")


def input_pcr_geometry():
    return PCRGeometry.exact(expected_reporter(),
        PCRPrimer("published input forward", INPUT_F, len(INPUT_F)),
        PCRPrimer("published shared reverse", REPORTER_R, len(REPORTER_R)),
        statement="Exact inward sites within the declared partial, inferred reporter context. Input forward is verified in SOM p5; the shared reverse and RT tail are verified in the matching author-linked DSCP record. Unknown interior intervals are not molecular coordinates, a full plasmid map or an experimental product guarantee.")


def input_amplicon():
    template = expected_reporter()
    geometry = input_pcr_geometry()
    product = template.slice("source input-forward site in intron",
        "source shared first-PCR reverse site", name="INFERRED partial input first-PCR product")
    if (geometry.forward_hit.free5 or geometry.reverse_hit.free5
            or not product.top().startswith(INPUT_F)
            or not product.top().endswith(revcomp(REPORTER_R))):
        raise ValueError("input amplicon must follow both complete source primer ends")
    return product


def source_splice_candidate():
    # SOM p5 matches the author record's full cassette; the author record has two
    # competing intron annotations. This is ONLY the402..468 candidate, not a
    # selected experimental transcript or a repaired primer.
    return "AAGCCACCATGGAAAAG" + "GCCAA"


def reporter_rna():
    sc = Scene()
    sc.strand("rna", (role("RNA population containing reporter RNA; splice, ends and primer sites unassigned"),),
              label="RNA class, single strand")
    return sc


class UnassignedAdapterGeometry:
    availability = "unavailable"
    ready_label = "source genomic fragments and PE-400-1001 adapters"
    action_label = "Present source kit-defined inputs without inventing end chemistry"
    statement = "The author deposit names E6000L and PE-400-1001. Their exact applied kit versions, end geometry and subsequent adapter retention are not assigned from the product names alone."

    def rows(self):
        return [Row(chunks=[("Fragment and adapter ends unassigned; no complete ligation-product relation is claimed", None, False)])]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    genomic = g.input("genomic", "source Drosophila genomic DNA", K.MOLECULAR,
                      state=unknown("source genomic DNA population"), primary=True)
    shear = g.input("shear", "Covaris S220 sonication; exact fragment sequences unassigned", K.CONDITION)
    fragmented = g.product("fragmented", "sheared genomic fragment population", K.MOLECULAR,
                           state=unknown("source genomic fragments"))
    g.react("fragment", "Shear the source genomic DNA", S.FRAGMENTATION,
        [(genomic, R.TEMPLATE), (shear, R.CONDITION)], [fragmented])
    size = g.input("size_selection", "published 500–700-bp agarose size selection", K.SELECTION)
    selected = g.product("selected", "source selected fragment population; exact lengths unassigned", K.MOLECULAR, state=fragmented.state)
    other = g.product("size_waste", "unselected genomic size classes", K.WASTE)
    g.react("select_fragments", "Select the source approximate fragment-size class", S.SELECTION,
        [(fragmented, R.TEMPLATE), (size, R.SELECTION)], [selected, other], changes_molecule=False)
    g.discard(other, "Physical size selection; no exact distribution, recovery fraction or genome coverage is inferred.")
    kit = g.input("kit", "NEBNext DNA Library Prep Reagent Set E6000L", K.CATALYST)
    prepared = g.product("prepared", "kit-prepared fragments; end chemistry unassigned", K.MOLECULAR, state=selected.state)
    g.react("kit_prepare", "Unexpanded source E6000L preparation", S.MAGIC,
        [(selected, R.TEMPLATE), (kit, R.CATALYST)], [prepared],
        reason="The deposit refers to the kit instructions. Applied original kit version and molecular end-preparation steps await the corresponding authoritative original; no modern kit recipe is substituted.")
    adapters = g.input("adapters", "source Illumina Multiplexing Adapters PE-400-1001", K.ADAPTER)
    ligase = g.input("adapter_ligase", "source adapter-ligation component; identity unassigned", K.CATALYST)
    adapted = g.product("adapted", "source adapter-ligated fragment classes; complete structure unassigned", K.MOLECULAR,
                        state=unknown("source adapter-ligated genomic fragments"))
    g.react("attach_adapters", "Ligate the source-named adapters", S.LIGATION,
        [(prepared, R.TEMPLATE), (adapters, R.STRUCTURAL), (ligase, R.CATALYST)], [adapted], geometry=UnassignedAdapterGeometry())
    cloning_f = g.input("cloning_f", "published cloning forward", K.OLIGO, state=primer("cloning forward", CLONING_F))
    cloning_r = g.input("cloning_r", "published cloning reverse", K.OLIGO, state=primer("cloning reverse", CLONING_R))
    cloning_pcr = g.input("cloning_pcr", "KAPA HiFi HotStart ReadyMix (KK2602)", K.CATALYST)
    homology = g.product("homology", "source PCR product with In-Fusion arms; exact sites unassigned", K.MOLECULAR,
                         state=unknown("candidate fragments with source homology-arm context unassigned"))
    g.react("add_homology", "Add source In-Fusion homology arms by PCR", S.PCR,
        [(adapted, R.TEMPLATE), (cloning_f, R.STRUCTURAL), (cloning_r, R.STRUCTURAL), (cloning_pcr, R.CATALYST)], [homology],
        geometry=PCRGeometry.unavailable(adapted.state,
            reason="Corrected SOM p6 gives the two exact cloning primers and15nt homology extensions. The complete ligation-ready adapter/template ends and exact applied kit version remain unassigned; no precise landing coordinates are fabricated."))
    recipient = g.input("recipient", "source DSCP/mhc16-intron/sgGFP/ccdB pGL3-based recipient", K.MOLECULAR,
                        state=unknown("source recipient; complete map unassigned"))
    digest = g.input("recipient_digest", "AgeI-HF and SalI-HF", K.CATALYST)
    opened = g.product("opened_recipient", "digested reporter backbone; exact cut coordinates unassigned", K.MOLECULAR,
        state=unknown("AgeI/SalI prepared recipient; exact ends unassigned"))
    removed = g.product("removed_cassette", "excised recipient cassette", K.WASTE)
    g.react("open_recipient", "Digest recipient and isolate the backbone", S.DIGESTION,
        [(recipient, R.TEMPLATE), (digest, R.CATALYST)], [opened, removed])
    g.discard(removed, "SOM p6 specifies agarose separation and gel extraction of linearized recipient; exact excised sequence unassigned.")
    product = expected_reporter()
    assembled = g.product("assembled", "INFERRED intended source reporter construct", K.MOLECULAR,
                          state=closed_circle(product.name, tuple(product)), chemical_state=product)
    g.react("infusion", "INFERRED intended arrangement — unexpanded In-Fusion construction", S.MAGIC,
        [(homology, R.TEMPLATE), (opened, R.STRUCTURAL)], [assembled],
        reason="The source names In-Fusion recombination into the reporter. AgeI-HF/SalI-HF preparation and the two cloning oligos are source-defined. Full recipient/adapter end states and exact product junctions remain unassigned.")
    plasmids = g.product("plasmids", "INFERRED recovered reporter-plasmid class", K.MOLECULAR, state=assembled.state, chemical_state=product)
    g.react("recover_plasmids", "INFERRED intended arrangement — source clone propagation and isolation", S.MAGIC,
        [(assembled, R.TEMPLATE)], [plasmids],
        reason="The deposit describes transformation, growth and plasmid isolation. These population processes are unexpanded and do not prove genotype, diversity or quantitative representation.")
    cells = g.input("osc_cells", "source OSC sample context; biological state unassigned", K.BIOLOGICAL, primary=True)
    total_rna = g.product("total_rna", "extracted OSC total RNA", K.MOLECULAR, state=reporter_rna())
    g.react("reporter_assay", "Unexpanded OSC reporter assay and total-RNA recovery", S.MAGIC,
        [(plasmids, R.TEMPLATE), (cells, R.TEMPLATE)], [total_rna],
        reason="SOM p7 specifies electroporation,24h expression and RNeasy extraction. This graph records the assay boundary, not delivery efficiency, precise transcription/splice boundaries or reporter activity.")
    selection = g.input("poly_a_selection", "Dynabeads Oligo-dT25", K.SELECTION)
    poly_a = g.product("poly_a", "selected poly(A) RNA population", K.MOLECULAR, state=reporter_rna())
    non_poly_a = g.product("non_poly_a", "unselected RNA", K.WASTE)
    g.react("select_rna", "Select poly(A) RNA", S.AFFINITY_CAPTURE,
        [(total_rna, R.TEMPLATE), (selection, R.SELECTION)], [poly_a, non_poly_a])
    g.discard(non_poly_a, "Oligo-dT selection removes the unselected population; it is not the gene-specific RT priming step.")
    dnase = g.input("dnase", "Ambion Turbo DNase", K.CATALYST)
    rna = g.product("reporter_rna", "DNase-treated poly(A) RNA; reporter sequences present", K.MOLECULAR, state=reporter_rna())
    dna_waste = g.product("dna_waste", "digested DNA contamination", K.WASTE)
    g.react("digest_dna", "Treat selected RNA with Turbo DNase", S.DIGESTION,
        [(poly_a, R.TEMPLATE), (dnase, R.CATALYST)], [rna, dna_waste])
    g.discard(dna_waste, "SOM p7 specifies Turbo DNase followed by RNeasy MinElute cleanup; complete background elimination is not established by this graph.")
    # The first strand is explicit, without assigning unverified splice coordinates.
    rt_primer = g.input("rt_primer", "published reporter-specific RT primer", K.OLIGO,
        state=primer("reporter RT", RT_PRIMER))
    rt = g.input("rt", "SuperScript III", K.CATALYST)
    cdna_scene = Scene()
    cdna_scene.strand("cdna", (role("reporter first-strand cDNA; splice and RT extent unassigned"),), label="5-prime to3-prime cDNA")
    cdna = g.product("cdna", "reporter first-strand cDNA", K.MOLECULAR, state=Drawing(tuple([*cdna_scene.rows(), *reporter_rna().rows()])))
    g.react("make_cdna", "Copy reporter RNA with the source-specific primer", S.REVERSE_TRANSCRIPTION,
        [(rna, R.TEMPLATE), (rt_primer, R.STRUCTURAL), (rt, R.CATALYST)], [cdna])
    rnase = g.input("rnase_a", "RNase A", K.CATALYST)
    clean_cdna = g.product("clean_cdna", "RNase-treated reporter cDNA", K.MOLECULAR, state=cdna_scene)
    degraded_rna = g.product("degraded_rna", "digested RNA", K.WASTE)
    g.react("remove_rna", "Digest RNA after first-strand synthesis", S.DIGESTION,
        [(cdna, R.TEMPLATE), (rnase, R.CATALYST)], [clean_cdna, degraded_rna])
    g.discard(degraded_rna, "SOM p8 specifies RNase A followed by QIAquick purification; RNA is not the PCR template.")
    ngs_polymerase = g.input("ngs_polymerase", "KAPA HiFi HotStart ReadyMix", K.CATALYST)
    for branch, template, fseq in (("rna", clean_cdna, REPORTER_F), ("input", plasmids, INPUT_F)):
        forward = g.input(branch + "_f", branch + " published first-PCR forward", K.OLIGO,
            state=primer(branch + " first-PCR forward", fseq))
        reverse = g.input(branch + "_r", branch + " published first-PCR reverse", K.OLIGO,
            state=primer(branch + " first-PCR reverse", REPORTER_R))
        amplified = g.product(branch + "_amplicon", branch + " first-PCR product; exact span unassigned", K.MOLECULAR,
            state=(input_amplicon() if branch == "input" else unknown(branch + " first-PCR product; primer-to-template coordinates unresolved")))
        g.react(branch + "_pcr", branch + " reporter-specific first PCR", S.PCR,
            [(template, R.TEMPLATE), (forward, R.STRUCTURAL), (reverse, R.STRUCTURAL), (ngs_polymerase, R.CATALYST)], [amplified],
            geometry=(input_pcr_geometry() if branch == "input" else PCRGeometry.unavailable(cdna_scene,
                reason="Corrected SOM p8 gives the exact primers. The reporter forward spans the mhc16 splice junction; input uses a distinct forward primer. The published forward ends in T; the author-record402..468 splice candidate gives A at that base. No candidate is silently selected or primer repaired.")))
        f2 = g.input(branch + "_f2", branch + " NEBNext E7335L second-PCR forward; exact oligo/version unassigned", K.OLIGO)
        r2 = g.input(branch + "_r2", branch + " NEBNext E7335L second-PCR reverse; exact oligo/index unassigned", K.OLIGO)
        state = unknown(branch + " source sequenced-library family; final adapters and index architecture unassigned")
        library = g.product(branch + "_library", branch + " source paired-end readout family", K.LIBRARY, state=state)
        g.react(branch + "_index_pcr", branch + " NEBNext second PCR", S.PCR,
            [(amplified, R.TEMPLATE), (f2, R.STRUCTURAL), (r2, R.STRUCTURAL), (ngs_polymerase, R.CATALYST)], [library],
            geometry=PCRGeometry.unavailable(amplified.state,
                reason="SOM p8 names NEBNext Multiplex Oligos E7335L; historical component sequences, index assignments and exact template geometry remain unassigned."))
        g.sequence(library, SequencingPlan(acquired_roles=("Read 1", "Read 2"),
            read_lengths={"Read 1": 36, "Read 2": 36},
            limitation="Corrected SOM Mapping p11 reports36mer paired-end Drosophila STARR-seq, separate from50mer single-end human STARR-seq. Sequencing primers, final retained adapter ends, experimental sample indexes and source-specific read-to-insert coordinates remain unassigned."))
    g.annotate(SourceDiscrepancy("add_homology", "Published cloning-PCR reaction count",
        (SourceClaim("Corrected SOM p6, start of cloning-PCR paragraph", "Eight PCR reactions"),
         SourceClaim("Same paragraph, pooling sentence", "The five PCR reactions were pooled")),
        "The graph models the shared reaction class without assigning either batch count.",
        "1232542.arnold.sm.revision1.pdf, printed p6, Generation of input libraries.",
        "The source does not reconcile eight versus five; no quantitative pooling or representation claim is made."))
    g.annotate(SourceDiscrepancy("rna_pcr", "Reporter splice-primer context remains unresolved",
        (SourceClaim("Corrected SOM p8, reporter first-PCR forward", REPORTER_F),
         SourceClaim("SOM p5 cassette, also identical in author DSCP record; deleting its annotated402..468 intron", source_splice_candidate())),
        "Retain the exact published primer; show no exact reporter-RNA PCR landing or corrected sequence.",
        "SOM p5 and p8 rendered originals inspected. The cassette matches the author DSCP record, which separately annotates399..469 and402..468 intron intervals. Neither candidate gives a full exact match to this forward primer.",
        "The402..468 candidate differs at the primer terminal base (A versus T). Experimental splice choice, mismatch tolerance and the source discrepancy remain unresolved; this is not an absence conclusion or proof the experiment failed."))
    g.declare_reviewed_scope(SOURCE_INVENTORY, SOURCE)
    return g.require_publishable()


SOURCE_INVENTORY = (
    "genomic", "shear", "fragment", "size_selection", "select_fragments", "kit", "kit_prepare", "adapters", "adapter_ligase", "attach_adapters",
    "cloning_f", "cloning_r", "cloning_pcr", "add_homology", "recipient", "infusion", "recover_plasmids", "osc_cells",
    "reporter_assay", "make_cdna", "ngs_polymerase", "rna_f", "rna_r", "rna_pcr", "rna_library",
    "input_f", "input_r", "input_pcr", "input_library",
    "poly_a_selection", "select_rna", "dnase", "digest_dna",
    "recipient_digest", "open_recipient", "rt_primer", "rt", "rnase_a", "remove_rna",
    "rna_f2", "rna_r2", "rna_index_pcr", "input_f2", "input_r2", "input_index_pcr",
)


def render_page():
    return assemble_page(protocol_graph(), (
        head(TITLE), '<div class="wrap" data-protocol-complete="true">', f"<h1>{TITLE}</h1>",
        '<p class="research-notes"><a href="01_starr_seq.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1126/science.1232542">Arnold et al., Science (2013)</a>. '
             '<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1000403">OSC transcript sample GSM1000403</a>; '
             '<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1000405">input library GSM1000405</a>.'),
        info("WORK IN PROGRESS: the published article and corrected supplement now define this scoped graph. Complete vector/adapter, splice and primer-site assignments remain unresolved. This entry remains notes and is not a finished public schematic. The OSC paired-end branch is kept separate from the human single-end experiment and later STARR-seq variants."),
        info("INFERRED: the intended reporter arrangement is conceptual; it is not a complete plasmid map. The author DSCP GenBank has a LOCUS/ORIGIN length discrepancy and is not imported wholesale: only the short downstream readout context is used after confirming the printed2502-base cassette matches exactly. Kit names do not assign a version-specific adapter recipe. The published primer bases, RT/PCR reagents, splice-selective amplification and36-base paired acquisition are recorded; exact reporter splice/RT geometry and sample-index structure remain unassigned. Input-DNA primer landing is checked on an inferred partial context with unknown interior intervals; no full plasmid or amplicon length is assigned."),
        AnnotationSlot.INTRO,
        source_oligos(),
        WorkflowSection("Reporter construction and separate RNA/input readouts", paths=(("opened_recipient", ()), ("plasmids", ()),
            ("rna_library", ("recover_plasmids",)), ("input_library", ("recover_plasmids",)))),
        "<h2>Readout</h2>", info("Source Drosophila acquisition is36-base paired-end; sequencing-primer sites remain unassigned. The insert itself identifies the tested region; no independent MPRA reporter tag or UMI is added by analogy. Molecular consistency does not establish reporter activity, expression, library balance, coverage or an endogenous enhancer-target relationship."),
        sp.read_spec_table((sp.read_spec("Read 1", 36, "source paired fragment; exact layout unassigned", "include"),
                            sp.read_spec("Read 2", 36, "paired fragment; exact layout unassigned", "include"))),
        "</div>",
    ))


def selftest():
    from checks import Check, Source, have
    from chemistry_validation import Status, check_graph
    from storage import sources_root
    check = Check(); graph = protocol_graph()
    check("no contradictory claims in the source-scoped published-method graph", [c for c in check_graph(graph) if c.status is Status.CONTRADICTORY], [])
    root = sources_root() / Path(__file__).resolve().parents[1].name
    for accession, label in (("GSM1000403", "PCR amplified cDNA (STARR-seq transcript)"),
                              ("GSM1000405", "PCR amplified Plasmid library, used for OSC STARR-seq")):
        source = Source(root / ("geo-" + accession + "-self-full.txt"), "Original GSE40739 deposited sample; tools/get_sources.py via catalogue/extra_sources.tsv")
        if have(check, source, label=accession + " source identity"):
            lines = source.path.read_text().splitlines()
            check(accession + " is the original record, not a shared CGI cache alias", lines[0], "^SAMPLE = " + accession)
            check(accession + " assigns the correct readout population",
                [x.split(" = ", 1)[1] for x in lines if x.startswith("!Sample_description = ")], [label])
            check(accession + " belongs to the defining series",
                [x for x in lines if x.startswith("!Sample_series_id = ")], ["!Sample_series_id = GSE40739"])
    supplement = Source(root / "1232542.arnold.sm.revision1.pdf.txt", "User-supplied published corrected SOM for doi:10.1126/science.1232542")
    if have(check, supplement, label="corrected published primer transcription"):
        text = supplement.path.read_text()
        for name, sequence in (("cloning forward", CLONING_F), ("cloning reverse", CLONING_R),
            ("reporter RT", RT_PRIMER), ("reporter first forward", REPORTER_F),
            ("reporter first reverse", REPORTER_R), ("input first forward", INPUT_F)):
            check(name + " occurs verbatim in corrected published supplement", sequence in text, True)
        check("selected supplement identifies its correction", "corrected 6 February 2013" in text, True)
    vector = Source(root / "STARR-seq_DSCP.gbk", "Original author-linked DSCP record; see external manifest")
    if have(check, supplement, vector, label="partial source vector context"):
        import re
        source_text = supplement.path.read_text()
        vector_text = vector.path.read_text()
        vector_bases = "".join(re.findall("[ACGTacgt]+", vector_text.split("ORIGIN",1)[1].split("//",1)[0])).upper()
        # Extract the first explicitly printed DNA block, not a name-matched map.
        groups = []; current = []
        for line in source_text.splitlines():
            if re.fullmatch("[ACGT]{60,}", line.strip()): current.append(line.strip())
            elif current: groups.append("".join(current)); current = []
        if current: groups.append("".join(current))
        cassette = groups[0] if groups else ""
        check("author cassette exactly matches the whole printed Drosophila block", len(cassette) >= 2000 and cassette in vector_bases, True)
        check("input forward is in that printed cassette", INPUT_F in cassette, True)
        left_arm = CLONING_F[:15]; right_arm = revcomp(CLONING_R[:15])
        check("both15nt cloning homologies match printed cassette", left_arm in cassette and right_arm in cassette, True)
        tail = "".join(part.top for part in recipient_readout_tail())
        check("short declared right context is contiguous in original vector", right_arm + tail in vector_bases, True)
        candidate = vector_bases[:401] + vector_bases[468:]
        start = candidate.find("AAGCCACCATGGAAAAG")
        check("record-derived402..468 candidate gives the displayed primer context", candidate[start:start+len(REPORTER_F)] if start >= 0 else None, source_splice_candidate())
    check.report()
