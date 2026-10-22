"""Melnikov 2012 MPRA: source-scoped construction and DNA/RNA Tag-Seq."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from circular import CircularizationGeometry
from cohesive_assembly import CircularCohesiveAssembly, CohesiveFragment, FivePrimeFragment
import illumina as il
from page import head, info
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (AnnotationSlot, DetailView, InputRole as R, MaterialKind as K, MolecularRelation,
                            ProtocolGraph, SequencingPlan, StepKind as S)
from protocol_model import closed_circle, role
from protocol_views import (OligoEntry, OligoSection, PrimerBindingSection,
                            ReadLayoutSection, SourceRead, WorkflowSection, assemble_page)
import seqprimers as sp

TITLE = "MPRA — reporter construction and Tag-Seq"
PAGE_FILE = "mpra.html"
SOURCE = "Melnikov et al. 2012, doi:10.1038/nbt.2137; Methods oligonucleotide design/plasmid construction/Tag-Seq; published Supplementary Figure1"
UNIVERSAL_F = "ACTGGCCGCTTCACTG"
UNIVERSAL_R = "CGACGCTCTTCCGATCT"
SFI_F = "GCTAAGGGCCTAACTGGCCGCTTCACTG"
SFI_R = "GTTTAAGGCCTCCGAGGCCGACGCTCTTC"
SFI_LEFT = "GGCCTAACTGGCC"
SFI_RIGHT = "GGCCTCGGAGGCC"
SFI_TOP_CUT, SFI_BOTTOM_CUT = 8, 5
KPNI = "GGTACC"
XBAI = "TCTAGA"
HINDIII = "AAGCTT"
TAG_F = il.TRUSEQ_P5_FULL
TAG_R_SITE = "CGAGGTGCCTAAAGG"
SAMPLE_INDEX_LENGTH = 8
TAG_R = il.P7 + "N" * SAMPLE_INDEX_LENGTH + il.TRUSEQ_READ2 + TAG_R_SITE
VARIANT_LENGTH = 87
TAG_LENGTH = 10
READ1_CYCLES = 36
REPORTER_INTERVAL = 90
INFERRED_SITE = "INFERRED gene-side primer site — "


def bp(name, bases, **kw):
    return Segment(name, bases, **kw)


def reporter_tag(sequence="N" * TAG_LENGTH):
    if len(sequence) != TAG_LENGTH or set(sequence) - set("ACGTN"):
        raise ValueError("this MPRA design uses a ten-base tag (or an unassigned N window)")
    return Segment("design-assigned reporter tag", sequence, placeholder="N" in sequence,
                   feature=feature("mpra.reporter_tag", "feature_barcode", "unknown" if "N" in sequence else "fixed",
                       note="Identifier of the associated regulatory design, not a UMI or sample index. The generic schematic selects no particular published design/tag pair."))


def oligo_template(tag):
    return Construct((bp("universal upstream site", UNIVERSAL_F),
                      role("associated regulatory variant", length_bp=VARIANT_LENGTH),
                      bp("KpnI", KPNI), bp("XbaI", XBAI), tag,
                      bp("universal downstream site", revcomp(UNIVERSAL_R))), "synthesized design/tag oligo")


def synthesis_length():
    return sum(s.length_bp or len(s) for s in oligo_template(reporter_tag()))


def part_range(template, start, end):
    """Copy a source-declared character interval, preserving region metadata.

    Display-role tokens are indivisible; no nucleotide coordinate is invented inside
    them. Every PCR endpoint in this scope lies in a concrete primer-binding region.
    """
    parts, offset = [], 0
    for seg in template:
        lo, hi = max(start, offset) - offset, min(end, offset + len(seg)) - offset
        if lo < hi:
            if seg.is_role_token() and (lo != 0 or hi != len(seg)):
                raise ValueError("a copied interval cannot split an unexpanded molecular role")
            whole = lo == 0 and hi == len(seg)
            parts.append(replace(seg, top=seg.top[lo:hi],
                                 name=seg.name if whole else seg.name + f" retained interval {lo}:{hi}",
                                 bottom=seg.bottom[lo:hi] if seg.bottom is not None else None,
                                 length_bp=seg.length_bp if whole else None))
        offset += len(seg)
    return tuple(parts)


def amplified_product(template, forward_parts, reverse_parts, geometry, label):
    if geometry.availability != "exact":
        raise ValueError("product derivation needs both declared PCR sites")
    f, r = geometry.forward_hit, geometry.reverse_hit
    # A PCR product incorporates the complete primers, including their paired
    # ends. Known primer bases replace any inference on the template landing;
    # the geometry still retains that uncertainty about the source template.
    middle = part_range(template, f.end, r.start)
    return Construct((*forward_parts, *middle, *complement_segments(reverse_parts)), label)


def single_strand(product, label):
    sc = Scene(); sc.strand("strand", tuple(product), label=label)
    return sc


def blunt_view(product, *, left_overhang=""):
    core = tuple(product)
    top = ((bp("left 5-prime extension", left_overhang),) if left_overhang else ()) + core
    bottom = complement_segments(core)
    sc = Scene(); sc.strand("top", top)
    sc.anneal("bottom", bottom, to="top", pair=(bottom[-1].name, core[0].name),
              unpaired=tuple(s.name for s in bottom if s.is_role_token()))
    if left_overhang:
        sc.overhang("top", "left 5-prime extension", polarity="5-prime", owner="substrate")
    else:
        sc.ligation_end("top", core[0].name, kind="blunt", owner="substrate", text="left blunt end")
    sc.ligation_end("top", core[-1].name, kind="blunt", owner="substrate", text="right blunt end")
    return sc


def circular_geometry(statement):
    return CircularizationGeometry.known("double-strand intramolecular ligation", statement=statement)


def vector_cassette():
    return Construct((bp("source left SfiI site", SFI_LEFT), role("initial SfiI-excised interval"),
                      bp("source right SfiI site", SFI_RIGHT), role("remaining pGL4.10M backbone")),
                     "source-defined pGL4.10M cassette")


def prepare_vector(g, enzymes, ligase):
    original = Construct((bp("HindIII", HINDIII), role("luc2 region removed from pGL4.10"),
                          bp("XbaI", XBAI), role("retained pGL4.10 backbone")), "pGL4.10 starting vector")
    raw = g.input("pgl410", "pGL4.10[luc2]", K.MOLECULAR,
                  state=closed_circle(original.name, tuple(original)), chemical_state=original)
    backbone = FivePrimeFragment("HindIII/XbaI retained backbone", XBAI[1:5],
        (bp("XbaI last base", XBAI[-1]), role("retained pGL4.10 backbone"),
         bp("HindIII first base", HINDIII[0])), HINDIII[1:5])
    opened = g.product("vector_open", backbone.name, K.MOLECULAR, state=backbone, chemical_state=backbone)
    discarded = g.product("luc2_removed", "excised original luc2 region", K.WASTE)
    g.react("remove_original_luc2", "Remove the original luc2 region with HindIII/XbaI", S.DIGESTION,
            [(raw, R.TEMPLATE), (enzymes["hindiii"], R.CATALYST), (enzymes["xbai"], R.CATALYST)],
            [opened, discarded])
    g.discard(discarded, "The original luc2 region is removed before the source-defined reporter cassette is installed.")
    klenow = g.input("klenow", "Klenow fragment", K.CATALYST)
    mung = g.input("mung", "Mung Bean nuclease", K.CATALYST)
    filled_core = Construct((*backbone.core, bp("HindIII end fill", HINDIII[1:5])), "right-end-filled backbone")
    filled = g.product("vector_filled", filled_core.name, K.MOLECULAR,
                       state=blunt_view(filled_core, left_overhang=backbone.left))
    g.react("fill_hindiii", "Fill the HindIII 5-prime extension", S.END_REPAIR,
            [(opened, R.TEMPLATE), (klenow, R.CATALYST)], [filled])
    blunted = g.product("vector_blunted", "backbone with the XbaI extension removed", K.MOLECULAR,
                        state=blunt_view(filled_core), chemical_state=filled_core)
    g.react("remove_xbai_extension", "Remove the XbaI 5-prime extension", S.END_REPAIR,
            [(filled, R.TEMPLATE), (mung, R.CATALYST)], [blunted])
    vector = g.product("pgl410m", "pGL4.10M: repaired source backbone", K.MOLECULAR,
                       state=closed_circle("reclosed repaired pGL4.10M backbone", tuple(filled_core)),
                       chemical_state=filled_core)
    g.react("close_vector", "Self-ligate to form the source pGL4.10M backbone", S.CIRCULARIZATION,
            [(blunted, R.TEMPLATE), (ligase, R.CATALYST)], [vector],
            geometry=circular_geometry("Methods: filled/trimmed backbone is self-ligated. The SfiI region is separately transcribed from Supplementary Figure1, not proved from a complete pGL4.10 map."),
            relations=(MolecularRelation.STRAND_ORDER,))
    g.annotate(DetailView("pgl410m", "Source-defined SfiI region of pGL4.10M",
        "Supplementary Figure1 supplies these selected SfiI flanks. Their placement in a complete prepared-vector map has not been checked. This regional source declaration is not DNA created during self-ligation.",
        state=vector_cassette(), slot=AnnotationSlot.AFTER_WORKFLOW))
    return vector


def site_index(product, sequence):
    sites = [i for i, s in enumerate(product) if s.top == sequence]
    if len(sites) != 1:
        raise ValueError("the selected source site must be represented exactly once")
    return sites[0]


def retained_kpn_xba(product):
    k, x = site_index(product, KPNI), site_index(product, XBAI)
    if k >= x:
        raise ValueError("source KpnI/XbaI order is reversed")
    core = (bp("XbaI retained last base", XBAI[-1]), *product.segments[x + 1:],
            *product.segments[:k], bp("KpnI retained first base", KPNI[0]))
    return CohesiveFragment("design/tag plasmid retained around the reporter gap", XBAI[1:5], core,
                            KPNI[1:5], "5-prime", "3-prime")


def reporter_fragment():
    # The gene-side Tag-Seq landing is inferred from the written PCR primer.
    # Its exact source-vector coordinate has not been checked against a full map.
    core = (bp("KpnI reporter last base", KPNI[-1]), role("minimal promoter and upstream luc2 region"),
            bp("INFERRED reporter Tag-Seq site", TAG_R_SITE, inferred=True),
            role("remaining reporter interval", length_bp=REPORTER_INTERVAL), bp("XbaI reporter first base", XBAI[0]))
    return CohesiveFragment("minimal-promoter/luc2 KpnI–XbaI fragment", KPNI[1:5], core,
                            XBAI[1:5], "3-prime", "5-prime")


def target_template(product):
    """Open the source circle in the orientation used by Read 1 PCR."""
    reverse = Construct(complement_segments(tuple(product)))
    hit = sp.locate(reverse, sp.SeqPrimer("PCR", "Tag-Seq forward", TAG_F,
                                        min_anneal=len(UNIVERSAL_R)))
    if hit is None:
        raise ValueError("the source downstream universal region is not represented")
    return Construct((*part_range(reverse, hit.start, len(reverse)), *part_range(reverse, 0, hit.start)),
                     "reporter template opened in Tag-Seq orientation")


def reverse_indexed_primer(sample):
    return (bp("P7", il.P7),
            Segment(sample + " sample index", "N" * SAMPLE_INDEX_LENGTH, placeholder=True,
                    feature=feature("mpra.sample." + sample, "sample_index", "unknown",
                                    note="Study-assigned sample index; value unassigned in this generic schematic.")),
            bp("TruSeq Read2 arm", il.TRUSEQ_READ2), bp("reporter priming end", TAG_R_SITE))


def assert_tag(product, expected, *, variant_expected, reverse=False):
    hits = [s for s in product if s.feature and s.feature.id == expected.feature.id]
    expected_bases = revcomp(expected.top) if reverse else expected.top
    if len(hits) != 1 or hits[0].feature != expected.feature or hits[0].top.upper() != expected_bases.upper():
        raise ValueError("the declared regulatory-design/tag association was changed or lost")
    variants = sum(s.top == f"[{VARIANT_LENGTH}-nt associated regulatory variant]" for s in product)
    if variants != int(variant_expected):
        raise ValueError("regulatory variant and transcript/readout regions were confused")


def _graph_for_tag(tag_sequence):
    g = ProtocolGraph(TITLE, authority=SOURCE)
    tag = reporter_tag(tag_sequence)
    enzymes = {key: g.input(key, name, K.CATALYST) for key, name in (
        ("sfii", "SfiI"), ("kpni", "KpnI"), ("xbai", "XbaI"), ("hindiii", "HindIII"))}
    ligase = g.input("ligase", "DNA ligase (identity unassigned in this source scope)", K.CATALYST)
    phusion = g.input("phusion", "Phusion high-fidelity PCR", K.CATALYST)
    pfu = g.input("pfu", "Pfu Ultra HS Tag-Seq PCR", K.CATALYST)
    gel = g.input("gel", "source gel size selection", K.SELECTION)
    universal = [g.input(key, label, K.OLIGO, state=Construct((bp(label, seq),)))
                 for key, label, seq in (("universal_f", "universal forward", UNIVERSAL_F),
                                         ("universal_r", "universal reverse", UNIVERSAL_R))]
    sfi_primers = [g.input(key, label, K.OLIGO, state=Construct((bp(label, seq),)))
                   for key, label, seq in (("sfi_f", "SfiI-adding forward", SFI_F),
                                           ("sfi_r", "SfiI-adding reverse", SFI_R))]
    synthesis = oligo_template(tag)
    oligo = g.input("synthesis", "pooled synthesized regulatory-design/tag oligos", K.OLIGO,
                    state=single_strand(synthesis, "synthesized single-strand oligo"),
                    chemical_state=synthesis, primary=True)
    duplex = g.product("initial_duplex", "amplified design/tag cassette", K.MOLECULAR, state=synthesis)
    g.react("initial_pcr", "Amplify the synthesized oligo pool", S.PCR,
            [(oligo, R.TEMPLATE), *[(p, R.STRUCTURAL) for p in universal], (phusion, R.CATALYST)],
            [duplex], geometry=PCRGeometry.unavailable(oligo.state,
                reason="The source starts from single-stranded synthesis oligos. Both written sites are known, but the duplex PCR geometry does not model the first reverse-primer copy of that initial strand."))
    selected = g.product("initial_selected", "size-selected design/tag cassette", K.MOLECULAR, state=synthesis)
    initial_waste = g.product("initial_waste", "unselected initial PCR material", K.WASTE)
    g.react("initial_select", "Size-select the initial PCR product", S.SELECTION,
            [(duplex, R.TEMPLATE), (gel, R.SELECTION)], [selected, initial_waste],
            changes_molecule=False, relations=(MolecularRelation.STRAND_ORDER,))
    g.discard(initial_waste, "Material outside the source target band; no recovery fraction is inferred.")
    sfi_geometry = PCRGeometry.exact(synthesis, PCRPrimer("SfiI forward", SFI_F, len(UNIVERSAL_F)),
        PCRPrimer("SfiI reverse", SFI_R, 11), statement="Written source primers: the reverse primer retains an 11-base 3-prime site on the original downstream universal sequence.")
    amplified = amplified_product(synthesis, (bp("SfiI forward", SFI_F),), (bp("SfiI reverse", SFI_R),),
                                   sfi_geometry, "SfiI-flanked design/tag cassette")
    amplified_node = g.product("sfi_amplified", amplified.name, K.MOLECULAR, state=amplified)
    g.react("add_sfii", "Re-amplify with the source SfiI-adding primers", S.PCR,
            [(selected, R.TEMPLATE), *[(p, R.STRUCTURAL) for p in sfi_primers], (phusion, R.CATALYST)],
            [amplified_node], geometry=sfi_geometry)
    vector = prepare_vector(g, enzymes, ligase)
    left_end = SFI_LEFT[SFI_BOTTOM_CUT:SFI_TOP_CUT]
    right_end = SFI_RIGHT[SFI_BOTTOM_CUT:SFI_TOP_CUT]
    left_site = amplified.top().index(SFI_LEFT)
    right_site = amplified.top().index(SFI_RIGHT)
    insert = CohesiveFragment("SfiI-released regulatory design/tag", left_end,
        part_range(amplified, left_site + SFI_TOP_CUT, right_site + SFI_BOTTOM_CUT),
        right_end, "3-prime", "3-prime")
    cassette = vector_cassette()
    left_site = cassette.top().index(SFI_LEFT)
    right_site = cassette.top().index(SFI_RIGHT)
    receiving = CohesiveFragment("pGL4.10M retained SfiI backbone", right_end,
        (*part_range(cassette, right_site + SFI_TOP_CUT, len(cassette)),
         *part_range(cassette, 0, left_site + SFI_BOTTOM_CUT)),
        left_end, "3-prime", "3-prime")
    insert_node = g.product("sfi_insert", insert.name, K.MOLECULAR, state=insert, chemical_state=insert)
    vector_node = g.product("sfi_vector", receiving.name, K.MOLECULAR, state=receiving, chemical_state=receiving)
    sfi_waste = g.product("sfi_waste", "unretained SfiI terminal and vector-cassette fragments", K.WASTE)
    g.react("sfi_release", "Release insert and receiving backbone at the source SfiI sites", S.DIGESTION,
            [(amplified_node, R.TEMPLATE), (vector, R.TEMPLATE), (enzymes["sfii"], R.CATALYST)],
            [insert_node, vector_node, sfi_waste])
    g.discard(sfi_waste, "Unretained source digest fragments; this does not prove complete plasmid digestion.")
    sfi_join = CircularCohesiveAssembly(((insert_node.id, insert), (vector_node.id, receiving)),
        "Supplementary Figure1: directional SfiI assembly uses distinct AAC and CGG 3-prime ends.")
    intermediate = sfi_join.product("design/tag plasmid before reporter insertion")
    assert_tag(intermediate, tag, variant_expected=True)
    intermediate_node = g.product("design_plasmid", intermediate.name, K.MOLECULAR,
                                  state=closed_circle(intermediate.name, tuple(intermediate)), chemical_state=intermediate)
    g.react("sfi_join", "Join the source SfiI insert and retained vector", S.LIGATION,
            [(insert_node, R.TEMPLATE), (vector_node, R.STRUCTURAL), (ligase, R.CATALYST)],
            [intermediate_node], geometry=sfi_join, relations=(MolecularRelation.COHESIVE_ASSEMBLY,))
    reporter = reporter_fragment()
    reporter_source = g.input("reporter_source", "pGL4.23 minimal-promoter/luc2 donor", K.MOLECULAR,
        state=closed_circle("pGL4.23 source regions", (bp("KpnI", KPNI), *reporter.core[1:-1],
            bp("XbaI", XBAI), role("pGL4.23 donor backbone"))))
    opened_geometry = retained_kpn_xba(intermediate)
    opened = g.product("reporter_gap", opened_geometry.name, K.MOLECULAR,
                       state=opened_geometry, chemical_state=opened_geometry)
    reporter_node = g.product("reporter_insert", reporter.name, K.MOLECULAR,
                              state=reporter, chemical_state=reporter)
    reporter_waste = g.product("reporter_waste", "donor backbone and excised cassette interval", K.WASTE)
    g.react("reporter_release", "Open the design/tag gap and release the KpnI–XbaI reporter", S.DIGESTION,
            [(intermediate_node, R.TEMPLATE), (reporter_source, R.TEMPLATE),
             (enzymes["kpni"], R.CATALYST), (enzymes["xbai"], R.CATALYST)],
            [opened, reporter_node, reporter_waste])
    g.discard(reporter_waste, "Independent unretained digest classes; no covalent backbone/cassette fusion or recovery fraction is inferred.")
    reporter_join = CircularCohesiveAssembly(((reporter_node.id, reporter), (opened.id, opened_geometry)),
        "Methods and Figure1: KpnI supplies a GTAC 3-prime end; XbaI supplies a CTAG 5-prime end. The promoter/luc2 fragment lies between regulatory variant and tag.")
    final = reporter_join.product("MPRA reporter construct")
    assert_tag(final, tag, variant_expected=True)
    assembled = g.product("reporter_plasmid", "INFERRED gene-side primer site — assembled reporter", K.MOLECULAR,
                          state=closed_circle(final.name, tuple(final)), chemical_state=final)
    g.react("reporter_join", "INFERRED gene-side primer site — insert the source reporter fragment", S.LIGATION,
            [(opened, R.TEMPLATE), (reporter_node, R.STRUCTURAL), (ligase, R.CATALYST)], [assembled],
            geometry=reporter_join, relations=(MolecularRelation.COHESIVE_ASSEMBLY,))
    # Background removal is a physical source operation, not a read filter.
    # The assembly drawing opens at the KpnI fusion. Cutting the restored site
    # leaves GTAC on each physical 3-prime end, not a single-stranded plasmid.
    if final.segments[0].top != KPNI[1:5]:
        raise ValueError("reporter circle drawing origin must be its declared KpnI fusion")
    linear_geometry = CohesiveFragment("KpnI-linearized reporter", KPNI[1:5], tuple(final.segments[1:]),
                                       KPNI[1:5], "3-prime", "3-prime")
    cut = g.product("background_cut", "source KpnI-linearized reporter target", K.MOLECULAR,
                    state=linear_geometry, chemical_state=linear_geometry)
    g.react("background_digest", INFERRED_SITE + "digest with KpnI for background removal", S.DIGESTION,
            [(assembled, R.TEMPLATE), (enzymes["kpni"], R.CATALYST)], [cut])
    selected_target = g.product("background_selected", "gel-selected linear reporter target", K.MOLECULAR, state=cut.state)
    background = g.product("background_waste", "unselected vector background and other material", K.WASTE)
    g.react("background_select", INFERRED_SITE + "size-select the reporter target", S.SELECTION,
            [(cut, R.TEMPLATE), (gel, R.SELECTION)], [selected_target, background], changes_molecule=False)
    g.discard(background, "Physical material outside the selected reporter target; this is separate from computational Tag-Seq read filtering.")
    recovered = g.product("recovered_plasmid", "reclosed reporter pool for DNA Tag-Seq", K.MOLECULAR,
                          state=assembled.state, chemical_state=final)
    g.react("background_close", INFERRED_SITE + "self-ligate the selected reporter target", S.CIRCULARIZATION,
            [(selected_target, R.TEMPLATE), (ligase, R.CATALYST)], [recovered],
            geometry=circular_geometry("Source KpnI digest, size selection and self-ligation; the full cut map and off-product frequencies remain unassigned."))
    # RNA begins after the unmodelled cellular assay, not as a proved transformation
    # of this illustrative plasmid into an experimentally verified transcript.
    mrna = g.input("total_rna", "total RNA collected after the reported reporter assay", K.BIOLOGICAL,
                  state=single_strand(Construct((role("total RNA including reporter transcripts"),)), "total RNA"), primary=True)
    polya = g.input("polya", "MicroPoly(A)Purist poly(A) selection", K.SELECTION)
    preparation = Scene()
    preparation.strand("rna", (role("poly(A) RNA with design-assigned reporter tags"),), label="selected RNA")
    preparation.strand("dna", (role("DNA targeted by the DNase step, amount unassigned"),),
                       label="potential DNA carryover; no quantity or strandedness assigned")
    selected_rna = g.product("mrna_selected", "poly(A)-selected reporter-containing preparation", K.MOLECULAR,
                             state=preparation)
    other_rna = g.product("other_rna", "unselected RNA", K.WASTE)
    g.react("select_mrna", "Select poly(A) RNA", S.SELECTION,
            [(mrna, R.TEMPLATE), (polya, R.SELECTION)], [selected_rna, other_rna], changes_molecule=False)
    g.discard(other_rna, "Unselected RNA; no target recovery or transcript abundance is established.")
    dnase = g.input("dnase", "DNase I / Turbo DNA-free", K.CATALYST)
    treated = g.product("dnase_rna", "DNase-treated poly(A) RNA", K.MOLECULAR,
                        state=single_strand(Construct((role("DNase-treated poly(A) RNA"),)), "retained RNA"))
    dna_products = g.product("dnase_products", "DNA-degradation products targeted for removal", K.WASTE)
    g.react("remove_dna", "DNase-treat the selected preparation", S.DIGESTION,
            [(selected_rna, R.TEMPLATE), (dnase, R.CATALYST)], [treated, dna_products])
    g.discard(dna_products, "DNase-targeted DNA is outside the intended RNA/cDNA readout; its original amount, complete removal and exact fragment sizes are unassigned.")
    rt = g.input("rt", "High Capacity RNA-to-cDNA reverse transcription kit", K.CATALYST)
    rt_primer = g.input("rt_primer", "source kit RT priming reagent, exact species unassigned", K.OLIGO,
                         state=Construct((role("unassigned kit RT primer species"),)))
    template = target_template(final)
    # Only the source target region is carried into the cDNA view; the upstream
    # regulatory variant controls expression but is not the sequenced transcript.
    f_query = sp.SeqPrimer("PCR", "Tag-Seq forward", TAG_F, min_anneal=len(UNIVERSAL_R))
    r_query = sp.SeqPrimer("PCR", "Tag-Seq reverse", TAG_R, min_anneal=len(TAG_R_SITE))
    f_hit = sp.locate(template, f_query)
    r_hit = sp.locate(template, r_query)
    if f_hit is None or r_hit is None:
        raise ValueError("source Tag-Seq target sites do not land")
    cdna = Construct(part_range(template, f_hit.start, r_hit.end), "reporter target within first-strand cDNA")
    assert_tag(cdna, tag, variant_expected=False, reverse=True)
    cdna_node = g.product("cdna", "INFERRED gene-side primer site — reporter cDNA target", K.MOLECULAR,
                          state=single_strand(cdna, "reporter-target region in cDNA"), chemical_state=cdna)
    g.react("reverse_transcribe", "INFERRED gene-side primer site — synthesize first-strand cDNA", S.REVERSE_TRANSCRIPTION,
            [(treated, R.TEMPLATE), (rt_primer, R.STRUCTURAL), (rt, R.CATALYST)], [cdna_node])
    tag_forward = g.input("tag_f", "Tag-Seq forward / P5 primer", K.OLIGO,
                         state=Construct((bp("P5/Read1 primer", TAG_F),)))
    for sample, source_node, source_template in (("dna", recovered, template), ("rna", cdna_node, cdna)):
        reverse_parts = reverse_indexed_primer(sample)
        reverse_sequence = "".join(s.top for s in reverse_parts)
        reverse_node = g.input(sample + "_tag_r", sample + " indexed Tag-Seq reverse", K.OLIGO,
                               state=Construct(reverse_parts))
        geometry = PCRGeometry.exact(source_template, PCRPrimer("Tag-Seq forward", TAG_F, len(UNIVERSAL_R)),
            PCRPrimer("Tag-Seq reverse", reverse_sequence, len(TAG_R_SITE)),
            statement="Primer ends are from Methods. Gene-side placement is inferred from their reported use and Figure1's output layout; a complete reporter map has not independently established that site.")
        library = amplified_product(source_template, (bp("P5/Read1 primer", TAG_F),), reverse_parts,
                                    geometry, sample + " Tag-Seq library")
        assert_tag(library, tag, variant_expected=False, reverse=True)
        amplified_node = g.product(sample + "_amplified", "INFERRED gene-side site — " + library.name, K.LIBRARY, state=library)
        drawn_geometry = geometry if sample == "dna" else PCRGeometry.unavailable(source_node.state,
            reason="The source RNA path starts from first-strand cDNA. Product derivation uses the expected duplex-cycle target and written primer tails; the initial complementary copy is not expanded into a false pre-existing duplex geometry.")
        g.react(sample + "_tag_pcr", "INFERRED gene-side primer site — make " + sample + " Tag-Seq library", S.PCR,
                [(source_node, R.TEMPLATE), (tag_forward, R.STRUCTURAL), (reverse_node, R.STRUCTURAL), (pfu, R.CATALYST)],
                [amplified_node], geometry=drawn_geometry)
        ready = g.product(sample + "_library", "gel-selected " + sample + " Tag-Seq library", K.LIBRARY, state=library)
        waste = g.product(sample + "_tag_waste", "unselected " + sample + " Tag-Seq products", K.WASTE)
        g.react(sample + "_tag_select", INFERRED_SITE + "select the " + sample + " Tag-Seq band", S.SELECTION,
                [(amplified_node, R.TEMPLATE), (gel, R.SELECTION)], [ready, waste], changes_molecule=False,
                relations=(MolecularRelation.STRAND_ORDER,))
        g.discard(waste, "Unselected PCR material; no recovery proportion or coverage is inferred.")
        primers = sp.standard_primers(library, ("Read 1", "Index 1 (i7)"))
        tag_cycles = [(span.cycle_start, span.cycle_end) for span in sp.feature_spans(
            library, primers, {"Read 1": READ1_CYCLES}) if span.feature.id == tag.feature.id]
        if tag_cycles != [(1, TAG_LENGTH)]:
            raise ValueError("the source's first ten Read1 cycles must identify the reporter tag")
        g.sequence(ready, SequencingPlan(library=library, primers=primers,
            acquired_roles=("Read 1", "Index 1 (i7)"), read_lengths={"Read 1": READ1_CYCLES},
            read_uses={"Read 1": "extract", "Index 1 (i7)": "demultiplex"},
            limitation="Indexed-run case. Methods also allow individually sequenced libraries. Index acquisition cycles and actual sample-index assignments remain unassigned; the eight-base index window alone does not set a run recipe. The gene-side PCR site remains inferred as marked. The regulatory design/tag association is declared, not experimentally proved."))
    g.declare_reviewed_scope((
        "sfii", "kpni", "xbai", "hindiii", "ligase", "phusion", "pfu", "gel",
        "universal_f", "universal_r", "sfi_f", "sfi_r", "synthesis", "initial_pcr", "initial_select", "add_sfii",
        "pgl410", "remove_original_luc2", "klenow", "mung", "fill_hindiii", "remove_xbai_extension", "close_vector",
        "sfi_release", "sfi_join", "reporter_source", "reporter_release", "reporter_join",
        "background_digest", "background_select", "background_close", "total_rna", "polya", "select_mrna",
        "dnase", "remove_dna", "rt", "rt_primer", "reverse_transcribe", "tag_f", "dna_tag_r", "rna_tag_r",
        "dna_tag_pcr", "rna_tag_pcr", "dna_tag_select", "rna_tag_select", "dna_library", "rna_library",
    ), SOURCE)
    return g.require_publishable()


def protocol_graph():
    return _graph_for_tag("N" * TAG_LENGTH)


def render_page():
    g = protocol_graph()
    content = [head(TITLE), '<div class="wrap" data-protocol-complete="true">', f"<h1>{TITLE}</h1>",
        '<p class="research-notes"><a href="01_mpra.html">Research notes</a></p>',
        info('<a href="https://doi.org/10.1038/nbt.2137">Melnikov et al., Nature Biotechnology (2012)</a>'),
        info("Generic design/tag architecture from the published Methods and Supplementary Figure1. Published design tables exist; this schematic selects no particular regulatory variant or tag. Full vector maps and RT primer species remain unassigned. INFERRED: the reporter-side Tag-Seq landing is inferred from the written primer and final-library figure; its template bases are marked. PCR products incorporate the known written primer ends. Cohesive-end, primer-binding and identifier checks do not establish full restriction maps, clone genotype, reporter activity, association fidelity, library coverage or balance. Cellular stimulation and phenotype experiments are outside this construction/readout graph."),
        info("The design/tag insert uses distinct SfiI AAC/CGG 3′ ends. Reporter insertion uses KpnI GTAC 3′ and XbaI CTAG 5′ ends. Background removal is a later physical KpnI digest, gel selection and self-ligation. The plasmid and RNA inputs remain separate: the graph does not infer measured transcript abundance from a construct's sequence."),
        WorkflowSection("Construction and Tag-Seq", paths=(
            ("sfi_amplified", ()), ("pgl410m", ()), ("design_plasmid", ("add_sfii", "close_vector")),
            ("recovered_plasmid", ("sfi_join",)), ("cdna", ()),
            ("dna_library", ("background_close",)), ("rna_library", ("reverse_transcribe",)),
        )),
        AnnotationSlot.AFTER_WORKFLOW,
        OligoSection("Written source primers", (
            OligoEntry("universal forward", (bp("universal forward", UNIVERSAL_F),)),
            OligoEntry("universal reverse", (bp("universal reverse", UNIVERSAL_R),)),
            OligoEntry("SfiI-adding forward", (bp("SfiI forward", SFI_F),)),
            OligoEntry("SfiI-adding reverse", (bp("SfiI reverse", SFI_R),)),
            OligoEntry("Tag-Seq forward", (bp("Tag-Seq forward", TAG_F),)),
            OligoEntry("Tag-Seq reverse", reverse_indexed_primer("dna")),
        )),
        info("The source writes the sample-index window as eight X positions; N denotes the same unassigned window here. The ten-base reporter identifier is a distinct structured feature. Read1 reports the reverse-complement orientation of the synthesis-strand tag; its first ten cycles identify the design and the following cycles check upstream construct sequence."),
    ]
    for sample in ("dna", "rna"):
        content.extend((PrimerBindingSection(sample + "_library", sample.upper() + " sample sequencing primers",
            "Computed primer sites for " + sample.upper()),
            ReadLayoutSection(sample + "_library", (
                SourceRead("Read 1", str(READ1_CYCLES), "reporter tag + upstream construct sequence", "Read 1", "extract"),
                SourceRead("Index 1 (i7)", "unassigned", "sample index", "Index 1 (i7)", "demultiplex")),
                "The source counts only reads with a designed ten-base tag and the expected subsequent construct sequence; molecular consistency does not establish the observed RNA/DNA ratio.")))
    content.append("</div>")
    return assemble_page(g, content)


def selftest():
    from checks import Check
    from chemistry_validation import Status, check_graph
    check = Check()
    graph = protocol_graph()
    check("source-scoped MPRA graph has no contradictory molecular claims",
          [c for c in check_graph(graph) if c.status is Status.CONTRADICTORY], [])
    # Non-palindromic synthetic tag probes orientation; this is not a published
    # enhancer design or a recommendation for a source library.
    tag = "ACGTTGCAAG"
    fixture = _graph_for_tag(tag)
    for sample in ("dna", "rna"):
        plan = fixture.sequencing[sample + "_library"]
        first = sp.locate(plan.library, plan.primers[0])
        check(sample + " Read1 reports the reverse-complemented synthesis tag", first.reads[:TAG_LENGTH], revcomp(tag))
        written_end = plan.library.get("reporter priming end'")
        check(sample + " product end is supplied by the known primer rather than an inferred template base",
              (written_end.top, written_end.inferred), (revcomp(TAG_R_SITE), False))
    original = fixture.materials["design_plasmid"].chemical_state
    inverted_tag = Construct([replace(s, top=revcomp(s.top)) if s.feature else s for s in original])
    check.raises("inverting a synthesis-strand tag breaks its declared source association",
                 lambda: assert_tag(inverted_tag, reporter_tag(tag), variant_expected=True))
    check("background removal remains physical and precedes DNA Tag-Seq",
          [graph.steps[s].kind for s in ("background_digest", "background_select", "background_close")],
          [S.DIGESTION, S.SELECTION, S.CIRCULARIZATION])
    check.report()
