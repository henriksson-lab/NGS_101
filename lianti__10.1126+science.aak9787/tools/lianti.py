"""Construct model for Chen et al. LIANTI (Science, 2017).

The defining paper establishes the architecture and reaction order (Fig. 1B/C legends:
a 19-bp double-stranded transposase binding site closed by a single-stranded T7-promoter
loop; gap extension to double-stranded promoters; overnight IVT; self-primed RT; RNase
digestion; UMI-bearing second strand), but prints no oligo.  Its Supplementary Materials
and Methods, which would, is a DOCX behind the PMC download gate.

Nothing here comes from the scg_lib_structs drawn page, which is inspiration only.  The
consequences, all visible in the model:

* the transposon's internal spacer between the mosaic end and the T7 promoter is an
  unverified placeholder -- neither its bases nor its length are established;
* the mosaic ends are the canonical Tn5 ME (``nextera``) and the promoter the canonical
  T7 class III promoter plus its GGG start, both inferred identifications of parts the
  paper names only by role;
* the UMI and index widths are drawing placeholders; the available main article does not give their lengths;
* no library-prep kit is named by the paper, so the adapter is an unspecified indexed
  Illumina adapter with unavailable ligation geometry, the library PCR geometry is
  unavailable, and every adapter-derived base in the drawn library is inferred.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp

UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"UMI": UMI_FEATURE, "i7": I7_FEATURE}
from endprep import dA_tailed_scene, repair_and_dA_tail
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (Disposition,Terminal,InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind)
from topology import HairpinTopology


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


# The paper gives the transposon's architecture, not its sequence.  The interval between
# the mosaic end and the T7 promoter is unverified in both bases and length: it is drawn
# as a placeholder of SPACER_NT positions purely so the loop has a drawable extent.
SPACER_NT = 10
SPACER = "N" * SPACER_NT
# Canonical T7 class III promoter and its first transcribed bases.  The paper names a
# "T7 promoter"; that this is the sequence used is an inference, so every segment built
# from it is created inferred.
T7_PROMOTER = "TAATACGACTCACTATA"
T7_START = "GGG"
TRANSPOSON = nx.ME_RC + SPACER + T7_PROMOTER + T7_START + nx.ME
TRANSPOSON_TOPOLOGY = HairpinTopology(
    "LIANTI hairpin transposon", TRANSPOSON,
    first_stem_start=0, second_stem_start=len(TRANSPOSON) - len(nx.ME),
    stem_length=len(nx.ME),
)
# Drawing widths for unverified identifier lengths.  The paper establishes that the
# second-strand primer carries a unique molecular barcode and that the library was
# indexed; neither width is published.
UMI_NT = 8
I7_NT = 6


def transposon_segments(*, inferred: bool = True) -> list[Segment]:
    """The single hairpin-forming transposon, written 5' to 3'.

    Its total length remains unverified because the spacer length is not given in the available main article.
    """
    return [
        seg("ME reverse complement", nx.ME_RC, "me", inferred=inferred),
        seg("unverified spacer", SPACER, placeholder=True, inferred=inferred,
            note="bases and length not verified from the available main article"),
        seg("T7 promoter", T7_PROMOTER, "t7", inferred=inferred),
        seg("T7 start", T7_START, "t7", inferred=inferred),
        seg("mosaic end", nx.ME, "me", inferred=inferred),
    ]


def second_strand_primer_segments() -> list[Segment]:
    return [
        seg("UMI", "N" * UMI_NT, "umi", placeholder=True, inferred=True),
        seg("T7 start", T7_START, "t7", inferred=True),
        seg("mosaic end", nx.ME, "me", inferred=True),
    ]


def gap_filled_fragment(insert_nt: int = 28) -> Construct:
    """A representative ds fragment after polymerase gap fill.

    The second end is the reverse complement of the first, so the two T7 promoters point
    inward.  All LIANTI-specific bases stay dotted because their exact sequence is not in
    the accessible defining paper.
    """
    insert = "X" * insert_nt if insert_nt <= 30 else "XXXXXXXXXXXX...XXXXXXXXXXXX"
    right = revcomp(TRANSPOSON)
    lib = Construct([
        *transposon_segments(),
        seg("genomic insert", insert, placeholder=True),
        seg("ME reverse complement, opposite", right[:len(nx.ME)], "me", inferred=True),
        seg("T7 start, opposite", right[len(nx.ME):len(nx.ME) + len(T7_START)],
            "t7", inferred=True),
        seg("T7 promoter, opposite",
            right[len(nx.ME) + len(T7_START):len(nx.ME) + len(T7_START) + len(T7_PROMOTER)],
            "t7", inferred=True),
        seg("unverified spacer, opposite",
            right[len(nx.ME) + len(T7_START) + len(T7_PROMOTER):
                  len(nx.ME) + len(T7_START) + len(T7_PROMOTER) + len(SPACER)],
            placeholder=True, inferred=True,
            note="bases and length not verified from the available main article"),
        seg("mosaic end, opposite", right[-len(nx.ME):], "me", inferred=True),
    ], name="gap-filled LIANTI fragment")
    return lib


def transcript(insert_nt: int = 28) -> Construct:
    """One IVT product in DNA letters (T represents U in the actual RNA)."""
    insert = "X" * insert_nt if insert_nt <= 30 else "XXXXXXXXXXXX...XXXXXXXXXXXX"
    return Construct([
        seg("5' GGG", T7_START, "t7", inferred=True),
        seg("5' mosaic end", nx.ME, "me", inferred=True),
        seg("genomic RNA", insert, placeholder=True),
        seg("self-primer site", nx.ME_RC, "me", inferred=True),
        seg("CCC", revcomp(T7_START), "t7", inferred=True),
        seg("promoter reverse complement", revcomp(T7_PROMOTER), "t7", inferred=True),
        seg("unverified spacer reverse complement", revcomp(SPACER),
            placeholder=True, inferred=True,
            note="bases and length not verified from the available main article"),
        seg("3' self-primer", nx.ME, "me", inferred=True),
    ], name="LIANTI RNA transcript")


def transcript_topology(insert_nt: int = 28) -> HairpinTopology:
    """Fold the terminal mosaic ends that self-prime reverse transcription."""
    tx = transcript(insert_nt)
    sequence = tx.top().replace("X", "N")
    first = len(T7_START) + len(nx.ME) + insert_nt
    return HairpinTopology(
        "self-priming LIANTI RNA", sequence,
        first_stem_start=first,
        second_stem_start=len(sequence) - len(nx.ME),
        stem_length=len(nx.ME),
        colour="var(--c-r2)",
    )


def umi_amplicon(insert_nt: int = 28) -> Construct:
    insert = "X" * insert_nt if insert_nt <= 30 else "XXXXXXXXXXXX...XXXXXXXXXXXX"
    return Construct([
        *second_strand_primer_segments(),
        seg("genomic copy", insert, placeholder=True),
    ], name="UMI-tagged LIANTI amplicon")


def end_prepared_amplicon(insert_nt: int = 28):
    """The standard NEBNext-polished, dA-tailed LIANTI amplicon."""
    return repair_and_dA_tail(umi_amplicon(insert_nt), inferred=True)


def end_prep_scene(insert_nt: int = 28):
    return dA_tailed_scene(end_prepared_amplicon(insert_nt), label="LIANTI amplicon")


SEQ_PRIMERS = [sp.TRUSEQ[k] for k in ("R1", "I1", "I2", "R2")]


def final_library(insert_nt: int = 28) -> Construct:
    """One of the two possible adapter-ligation orientations, P5 to P7'.

    The paper names no library-prep kit, so no adapter base here is established: the
    drawn arms are the canonical single-index Illumina/TruSeq parts, every one of them
    created inferred, and the index width is a placeholder.
    """
    insert = "X" * insert_nt if insert_nt <= 30 else "XXXXXXXXXXXX...XXXXXXXXXXXX"
    d_a = end_prepared_amplicon(insert_nt).three_prime_overhang
    lib = Construct([
        # P5 overlaps the first four bases of the Read 1 site (ACAC).
        seg("P5", il.P5[:-4], "p5", inferred=True),
        seg("Read 1 site", il.TRUSEQ_READ1, "r1", inferred=True),
        seg("UMI", "N" * UMI_NT, "umi", placeholder=True, inferred=True),
        seg("T7 start", T7_START, "t7", inferred=True),
        seg("mosaic end", nx.ME, "me", inferred=True),
        seg("genomic insert", insert, placeholder=True),
        seg("dA", d_a, inferred=True),
        seg("index-read site", il.INDEX1_PRIMER, "r2", inferred=True),
        seg("i7", "I" * I7_NT, "cbc", placeholder=True, inferred=True),
        seg("P7'", il.P7_RC, "p7", inferred=True),
    ], name="LIANTI sequencing library, UMI at P5 end")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid LIANTI final library: " + "; ".join(problems))
    return lib


def protocol_graph() -> ProtocolGraph:
    source=("Defining source: https://doi.org/10.1126/science.aak9787; "
            "the paper establishes the architecture but prints no oligo, and its "
            "Supplementary Materials and Methods is behind the PMC download gate")
    graph=ProtocolGraph("LIANTI",authority=source)
    genome=graph.input("genome","single-cell genomic DNA",MaterialKind.BIOLOGICAL,
                       primary=True)
    transposon=graph.input("hairpin_transposon","T7-promoter hairpin transposon",
                           MaterialKind.ADAPTER,state=tuple(transposon_segments()),
                           inferred=True,source=source)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagmented_genome","hairpin-transposon-tagmented genomic fragments",
                         MaterialKind.MOLECULAR,inferred=True,source=source)
    tag_waste=graph.product("tagment_waste","non-library tagmentation products",
                            MaterialKind.WASTE)
    graph.react("tagment","INFERRED — insert inward-facing T7-promoter transposons",
                StepKind.TAGMENTATION,
                ((genome,InputRole.TEMPLATE),(transposon,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,tag_waste))
    graph.discard(tag_waste,"fragments without usable paired promoter ends do not yield the modelled library")
    fill_pol=graph.input("gap_fill_polymerase","gap-fill DNA polymerase",MaterialKind.CATALYST)
    dntps=graph.input("gap_fill_dntps","gap-fill dNTPs",MaterialKind.OLIGO)
    filled=graph.product("gap_filled",gap_filled_fragment().name,MaterialKind.MOLECULAR,
                         state=gap_filled_fragment(),inferred=True,source=source)
    graph.react("gap_fill","INFERRED — fill the 9-nt Tn5 gaps and complete both T7 promoters",
                StepKind.EXTENSION,
                ((tagged,InputRole.TEMPLATE),(dntps,InputRole.STRUCTURAL),
                 (fill_pol,InputRole.CATALYST)),(filled,))
    t7=graph.input("t7_polymerase","T7 RNA polymerase",MaterialKind.CATALYST)
    rna=graph.product("rna_transcript",transcript().name,MaterialKind.MOLECULAR,
                      state=transcript(),inferred=True,source=source)
    graph.react("ivt","INFERRED — linearly amplify each fragment by T7 transcription",
                StepKind.TRANSCRIPTION,
                ((filled,InputRole.TEMPLATE),(t7,InputRole.CATALYST)),(rna,))
    reverse_transcriptase=graph.input("reverse_transcriptase","reverse transcriptase",
                                     MaterialKind.CATALYST)
    self_primed=graph.product("self_primed_cdna","self-primed LIANTI first-strand cDNA",
                              MaterialKind.MOLECULAR,inferred=True,source=source)
    graph.react("self_primed_rt","Self-prime reverse transcription from the RNA 3-prime end",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna,InputRole.TEMPLATE),(rna,InputRole.SELF_PRIMER),
                 (reverse_transcriptase,InputRole.CATALYST)),(self_primed,))
    rna_removal=graph.input("rna_removal","RNase digestion reagent; exact study choice unassigned",MaterialKind.CATALYST)
    cleaned=graph.product("cleaned_cdna","RNA-free self-primed cDNA",MaterialKind.MOLECULAR,
                          inferred=True,source=source)
    rna_waste=graph.product("rna_waste","LIANTI RNA digestion products; recovery unassigned",MaterialKind.MOLECULAR)
    graph.react("remove_rna","Digest the RNA template",StepKind.DIGESTION,
                ((self_primed,InputRole.TEMPLATE),(rna_removal,InputRole.CATALYST)),
                (cleaned,rna_waste))
    graph.mark_terminal(Terminal(rna_waste.id,Disposition.UNRESOLVED,"Original Figure1C describes RNase digestion before second-strand synthesis. Individual digestion products and their recovery/removal are not assigned by this model."))
    umi_primer=graph.input("umi_primer","UMI-bearing second-strand primer",MaterialKind.OLIGO,
                           state=tuple(second_strand_primer_segments()),inferred=True,source=source)
    second_pol=graph.input("second_strand_polymerase","second-strand DNA polymerase",
                           MaterialKind.CATALYST)
    amplicon=graph.product("umi_amplicon",umi_amplicon().name,MaterialKind.MOLECULAR,
                           state=umi_amplicon(),inferred=True,source=source)
    graph.react("second_strand","INFERRED — synthesize the UMI-bearing second strand",StepKind.EXTENSION,
                ((cleaned,InputRole.TEMPLATE),(umi_primer,InputRole.STRUCTURAL),
                 (second_pol,InputRole.CATALYST)),(amplicon,))
    repair=graph.input("end_repair_mix","end-repair enzyme mix",MaterialKind.CATALYST)
    repaired=graph.product("repaired_amplicon","end-repaired LIANTI amplicon",
                           MaterialKind.MOLECULAR,state=umi_amplicon(),inferred=True,source=source)
    graph.react("end_repair","INFERRED — end-repair LIANTI amplicons",StepKind.END_REPAIR,
                ((amplicon,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
    tailing=graph.input("tailing_polymerase","dA-tailing polymerase",MaterialKind.CATALYST)
    datp=graph.input("datp","dATP",MaterialKind.OLIGO)
    prepared_scene=end_prep_scene()
    prepared=graph.product("da_tailed","dA-tailed LIANTI amplicon",MaterialKind.MOLECULAR,
                           state=prepared_scene,inferred=True,source=source)
    graph.react("da_tail","INFERRED — add 3-prime dA overhangs",StepKind.TAILING,
                ((repaired,InputRole.TEMPLATE),(datp,InputRole.STRUCTURAL),
                 (tailing,InputRole.CATALYST)),(prepared,))
    adapter=graph.input("indexed_adapter","unspecified indexed Illumina adapter",
                        MaterialKind.ADAPTER,inferred=True,
                        source=source+" (study-specific kit and adapter remain unassigned; original SOM uninspected)")
    ligase=graph.input("ligase","DNA ligase",MaterialKind.CATALYST)
    opened=graph.product("adapter_ligated","indexed-adapter-ligated LIANTI amplicon",
                         MaterialKind.MOLECULAR,inferred=True,source=source)
    graph.react("adapter_ligation","INFERRED — ligate an indexed Illumina adapter to the dA-tailed amplicon",
                StepKind.LIGATION,
                ((prepared,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(opened,),
                geometry=AdapterLigationGeometry.unavailable(
                    "Study-specific kit and adapter are unassigned pending original SOM review; the end geometry, "
                    "polarity and bases have not been assembled in this model."))
    forward=graph.input("pcr_forward","library PCR primer completing the P5 arm",MaterialKind.OLIGO,
                        inferred=True,source=source+" (ordered study primer remains unassigned)")
    reverse=graph.input("pcr_reverse","indexed library PCR primer completing the P7 arm",
                        MaterialKind.OLIGO,inferred=True,
                        source=source+" (ordered study primer remains unassigned)")
    pcr_pol=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,
                          inferred=True,source=source)
    graph.react("library_pcr","INFERRED — amplify and complete the indexed library",StepKind.PCR,
                ((opened,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(
                    final,reason="The study-specific kit, ordered primers and their template sites remain unassigned pending original SOM review; this is not an absence conclusion."))
    graph.sequence(library,SequencingPlan(
        library=final,primers=tuple(SEQ_PRIMERS),
        limitation="The inspected main article and Figure1 establish amplified DNA ready for library preparation. Actual kit, indexing, physical acquisition and read lengths are unassigned pending the original SOM. Drawn TruSeq arms and identifier widths are reference assumptions, not source absence or an acquired run."))
    graph.declare_reviewed_scope((
        "genome","hairpin_transposon","tn5","tagment","gap_fill_polymerase",
        "gap_fill_dntps","gap_fill","t7_polymerase","ivt","reverse_transcriptase",
        "self_primed_rt","rna_removal","remove_rna","umi_primer","second_strand_polymerase",
        "second_strand","end_repair_mix","end_repair","tailing_polymerase","datp",
        "da_tail","indexed_adapter","ligase","adapter_ligation",
        "pcr_forward","pcr_reverse","pcr_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
