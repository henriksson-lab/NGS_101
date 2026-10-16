"""Reusable molecular states for non-standard RNA amplification workflows.

These helpers describe transformations of RNA itself.  They deliberately stop before
platform adapters are installed; protocol modules remain responsible for the exact final
library and sequencing-primer declarations.
"""
from __future__ import annotations

from chemdraw import Construct, Row, Scene, Segment, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, StepKind


def barcoded_ivt_front(graph, primary, *, rt_primer_label, rt_state,
                       first_state, duplex_state, arna_state, prefix="",
                       pool_before_second=False, exonuclease_when=None,
                       final_pool_after_exonuclease=False):
    """Build a source-selected well-barcoded RT and T7-IVT provenance path.

    One representative well is modelled from its primary RNA.  ``peer_wells`` is an
    explicit second pooled template, so pooling cannot be represented as a cosmetic
    label on one molecule.  The helper stops at amplified RNA because ligation versus
    random-primed second RT is protocol-specific.
    """
    def key(name: str) -> str:
        return f"{prefix}_{name}" if prefix else name

    rt_primer = graph.input(key("rt_primer"), rt_primer_label,
                            MaterialKind.OLIGO, state=rt_state)
    reverse_transcriptase = graph.input(
        key("reverse_transcriptase"), "reverse transcriptase",
        MaterialKind.CATALYST)
    first = graph.product(key("first_strand"), "well-barcoded first-strand cDNA",
                          MaterialKind.MOLECULAR, state=first_state)
    free_primers = (graph.product(key("free_rt_primers"),
                                  "unincorporated RT primers after reverse transcription",
                                  MaterialKind.MOLECULAR)
                    if exonuclease_when is not None else None)
    graph.react(
        key("reverse_transcription"),
        "Reverse-transcribe with the well-barcoded anchored oligo(dT) primer",
        StepKind.REVERSE_TRANSCRIPTION,
        ((primary, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
         (reverse_transcriptase, InputRole.CATALYST)),
        (first, free_primers) if free_primers is not None else (first,))

    if exonuclease_when not in (None, "before_pool", "after_pool"):
        raise ValueError("exonuclease_when must be before_pool, after_pool, or None")
    if exonuclease_when is not None and not pool_before_second:
        raise ValueError("Exo I position is only modeled for first-strand pooling")
    if final_pool_after_exonuclease and exonuclease_when != "after_pool":
        raise ValueError("final pooling requires Exo I after first pooling")

    reviewed = [key("rt_primer"), key("reverse_transcriptase"),
                key("reverse_transcription")]

    def remove_primers(template, leftover_primers):
        exonuclease = graph.input(key("exonuclease_i"), "Exonuclease I",
                                  MaterialKind.CATALYST)
        cleaned = graph.product(key("exo_cleaned"),
                                "barcoded cDNA after unincorporated-primer removal",
                                MaterialKind.MOLECULAR, state=first_state)
        waste = graph.product(key("exo_primer_fragments"),
                              "digested unincorporated RT-primer fragments",
                              MaterialKind.WASTE)
        graph.react(key("remove_free_primers"),
                    "Remove unincorporated RT primers with Exonuclease I",
                    StepKind.DIGESTION,
                    ((template, InputRole.TEMPLATE),
                     (leftover_primers, InputRole.TEMPLATE),
                     (exonuclease, InputRole.CATALYST)), (cleaned, waste))
        graph.discard(waste, "digested free RT primers are removed before second-strand synthesis")
        reviewed.extend((key("free_rt_primers"), key("exonuclease_i"),
                         key("remove_free_primers"), key("exo_primer_fragments")))
        return cleaned

    if pool_before_second:
        if exonuclease_when == "before_pool":
            first = remove_primers(first, free_primers)
        peers = graph.input(key("peer_wells"),
                            "corresponding barcoded first-strand cDNA from other wells",
                            MaterialKind.MOLECULAR, state=first_state)
        pooled_first = graph.product(key("pooled_first"),
                                     "pooled well-barcoded first-strand cDNA",
                                     MaterialKind.MOLECULAR, state=first_state)
        pooled_free = (graph.product(key("pooled_free_rt_primers"),
                                     "pooled unincorporated RT primers",
                                     MaterialKind.MOLECULAR)
                       if exonuclease_when == "after_pool" else None)
        graph.react(key("pool"), "Pool barcoded first-strand cDNA",
                    StepKind.POOL,
                    ((first, InputRole.TEMPLATE), (peers, InputRole.TEMPLATE)),
                    (pooled_first,), changes_molecule=False)
        reviewed.extend((key("peer_wells"), key("pool")))
        if pooled_free is not None:
            peer_free = graph.input(key("peer_free_rt_primers"),
                                    "free RT primers from corresponding wells",
                                    MaterialKind.MOLECULAR)
            graph.react(key("pool_free_rt_primers"),
                        "Free RT primers accompany the pooled first-strand cDNA",
                        StepKind.POOL,
                        ((free_primers, InputRole.TEMPLATE),
                         (peer_free, InputRole.TEMPLATE)),
                        (pooled_free,), changes_molecule=False)
            reviewed.extend((key("peer_free_rt_primers"),
                             key("pool_free_rt_primers"),
                             key("pooled_free_rt_primers")))
        first = (remove_primers(pooled_first, pooled_free)
                 if exonuclease_when == "after_pool" else pooled_first)
        if final_pool_after_exonuclease:
            peer_pools = graph.input(key("peer_pools"),
                                     "other Exo I-treated first-strand cDNA pools",
                                     MaterialKind.MOLECULAR, state=first_state)
            final_pool = graph.product(key("final_pool"),
                                       "final first-strand cDNA pool for amplification",
                                       MaterialKind.MOLECULAR, state=first_state)
            graph.react(key("pool_after_exonuclease"),
                        "Pool cleaned first-strand cDNA for second-strand synthesis",
                        StepKind.POOL,
                        ((first, InputRole.TEMPLATE),
                         (peer_pools, InputRole.TEMPLATE)),
                        (final_pool,), changes_molecule=False)
            reviewed.extend((key("peer_pools"), key("pool_after_exonuclease")))
            first = final_pool

    second_polymerase = graph.input(
        key("second_strand_polymerase"), "second-strand DNA polymerase",
        MaterialKind.CATALYST)
    duplex = graph.product(key("duplex"), "double-stranded barcoded T7 template",
                           MaterialKind.MOLECULAR, state=duplex_state)
    graph.react(key("second_strand"), "Synthesize the second cDNA strand",
                StepKind.EXTENSION,
                ((first, InputRole.TEMPLATE),
                 (second_polymerase, InputRole.CATALYST)), (duplex,))

    reviewed.extend((key("second_strand_polymerase"), key("second_strand")))
    if pool_before_second:
        pooled = duplex
    else:
        peers = graph.input(
            key("peer_wells"), "corresponding barcoded duplex cDNA from other wells",
            MaterialKind.MOLECULAR, state=duplex_state)
        pooled = graph.product(key("pooled"), "pooled well-barcoded duplex cDNA",
                               MaterialKind.MOLECULAR, state=duplex_state)
        graph.react(key("pool"), "Pool wells after molecular cell indexing",
                    StepKind.POOL,
                    ((duplex, InputRole.TEMPLATE), (peers, InputRole.TEMPLATE)),
                    (pooled,), changes_molecule=False)
        reviewed.extend((key("peer_wells"), key("pool")))

    t7 = graph.input(key("t7_polymerase"), "T7 RNA polymerase",
                     MaterialKind.CATALYST)
    arna = graph.product(key("arna"), "linearly amplified antisense RNA",
                         MaterialKind.MOLECULAR, state=arna_state)
    graph.react(key("ivt"), "Linearly amplify pooled barcoded cDNA by T7 IVT",
                StepKind.TRANSCRIPTION,
                ((pooled, InputRole.TEMPLATE), (t7, InputRole.CATALYST)),
                (arna,))
    reviewed.extend((key("t7_polymerase"), key("ivt")))
    return arna, reviewed


def poly_a_tagging_amplification(
        graph, template, *, tailed_state, tagging_primer_label,
        tagging_primer_state, duplex_state, suppression_primer_label,
        suppression_primer_sequence, amplified_state, prefix=""):
    """Model TdT tailing, second-strand tagging and symmetric suppression PCR.

    Quartz-family workflows share this molecular mechanism even though barcode and
    pooling occur at different points.  Exact opposite-facing PCR sites are checked on
    ``amplified_state``; a same-orientation pair cannot pass as symmetric amplification.
    """
    def key(name: str) -> str:
        return f"{prefix}_{name}" if prefix else name

    terminal_transferase = graph.input(
        key("terminal_transferase"), "terminal deoxynucleotidyl transferase",
        MaterialKind.CATALYST)
    tailed = graph.product(key("tailed_first_strand"),
                           "poly(A)-tailed first-strand cDNA",
                           MaterialKind.MOLECULAR, state=tailed_state)
    graph.react(key("poly_a_tail"), "Add a short poly(A) tail to first-strand cDNA",
                StepKind.TAILING,
                ((template, InputRole.TEMPLATE),
                 (terminal_transferase, InputRole.CATALYST)), (tailed,))

    tagging_primer = graph.input(key("tagging_primer"), tagging_primer_label,
                                 MaterialKind.OLIGO, state=tagging_primer_state)
    second_strand_polymerase = graph.input(
        key("second_strand_polymerase"), "second-strand DNA polymerase",
        MaterialKind.CATALYST)
    duplex = graph.product(key("tagged_duplex"),
                           "M-tagged double-stranded cDNA",
                           MaterialKind.MOLECULAR, state=duplex_state)
    graph.react(key("second_strand"),
                "Prime the poly(A) tail and synthesize the second strand",
                StepKind.EXTENSION,
                ((tailed, InputRole.TEMPLATE),
                 (tagging_primer, InputRole.STRUCTURAL),
                 (second_strand_polymerase, InputRole.CATALYST)), (duplex,))

    primer = PCRPrimer(suppression_primer_label, suppression_primer_sequence)
    suppression_primer = graph.input(
        key("suppression_primer"), suppression_primer_label,
        MaterialKind.OLIGO,
        state=Construct([Segment(suppression_primer_label,
                                 suppression_primer_sequence)]))
    polymerase = graph.input(key("suppression_polymerase"),
                             "suppression-PCR polymerase",
                             MaterialKind.CATALYST)
    amplified = graph.product(key("amplified_cdna"),
                              "suppression-amplified cDNA",
                              MaterialKind.MOLECULAR, state=amplified_state)
    graph.react(key("suppression_pcr"), "Suppress short byproducts and amplify cDNA",
                StepKind.PCR,
                ((duplex, InputRole.TEMPLATE),
                 (suppression_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.symmetric(
                    amplified_state, primer,
                    statement=("The single suppression primer has exact inward-facing "
                               "sites on the long cDNA product.")))
    reviewed = [
        key("terminal_transferase"), key("poly_a_tail"),
        key("tagging_primer"), key("second_strand_polymerase"),
        key("second_strand"), key("suppression_primer"),
        key("suppression_polymerase"), key("suppression_pcr"),
    ]
    return amplified, reviewed


def seg(name: str, bases: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, bases, tag, **kw)


def poly_a_capture(primer: list[Segment], *, annealing_segment: str,
                   transcript_width: int = 30, poly_a_length: int | None = None,
                   primer_label: str = "capture primer") -> Scene:
    """Anneal a compound 3-prime capture primer to an RNA poly(A) tail.

    The primer's annealing segment determines the tail length, preventing a page from
    drawing a partial overlap whose geometry was maintained separately by hand.
    """
    anneal = next((s for s in primer if s.name == annealing_segment), None)
    if anneal is None:
        raise ValueError(f"capture primer has no segment {annealing_segment!r}")
    n = poly_a_length or len(anneal)
    if n != len(anneal):
        raise ValueError("the drawn poly(A) tract must equal the annealing segment length")
    rna = [seg("transcript body", "X" * transcript_width, placeholder=True),
           seg("poly(A)", "A" * n)]
    sc = Scene(); sc.strand("mRNA", rna, label="mRNA")
    sc.anneal("primer", primer, to="mRNA", pair=(annealing_segment, "poly(A)"),
              label=primer_label,
              unpaired=tuple(s.name for s in primer if s.name != annealing_segment))
    sc.labels("primer")
    return sc


def terminal_tail(template: Construct, base: str, count: int, *, enzyme: str,
                  name: str | None = None) -> tuple[Construct, Scene]:
    """Append a homopolymer tail and draw it, deriving the product from the input."""
    if base.upper() not in "ACGT" or len(base) != 1 or count < 1:
        raise ValueError("terminal_tail needs one DNA base and a positive length")
    tail = seg(name or f"poly(d{base.upper()}) tail", base.upper() * count, "tso")
    out = Construct([*template.segments, tail], name=f"{template.name} + d{base.upper()} tail")
    sc = Scene(); sc.strand("product", list(out), label="tailed single strand")
    sc.mark("product", tail.name, f"{enzyme} adds d{base.upper()}")
    sc.labels("product")
    return out, sc


def direct_ssrna_ivt(template: list[Segment], promoter: Segment) -> Scene:
    """Draw a duplex T7 promoter attached to an otherwise single-stranded RNA template."""
    sc = Scene(); sc.strand("RNA template", [promoter, *template], label="original ssRNA")
    # ``Scene.anneal`` accepts every strand 5′→3′ and reverses its drawing, so the
    # partner must be the reverse complement, not the left-to-right plain complement.
    partner = Segment("T7 promoter partner", revcomp(promoter.top))
    sc.anneal("promoter partner", [partner], to="RNA template",
              pair=(partner.name, promoter.name), label="short promoter strand")
    sc.mark("RNA template", promoter.name, "duplex T7 promoter")
    sc.arrow("RNA template", "T7 RNA polymerase transcribes along the original ssRNA")
    return sc


def linear_amplification_rows(*, source: str, product: str) -> list[Row]:
    """Proportional-font explanation beneath a molecular IVT state."""
    return [Row(chunks=[(f"one {source} → many antisense {product} copies", None, False)])]
