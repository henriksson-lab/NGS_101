"""Construct model for Chen et al. LIANTI (Science, 2017).

The defining paper establishes the architecture and reaction order, but the accessible
main text does not print oligo sequences.  LIANTI-specific bases below are transcribed
from the upstream scg_lib_structs page and are therefore kept visibly secondary-source
in the schematic.  Standard Tn5 and Illumina/NEBNext parts come from the shared library.
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
from adapter_ligation import reported_adapter_end
from pcr import PCRGeometry
from protocol_graph import (InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind)
from topology import HairpinTopology


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


# LIANTI-specific sequence: secondary source (upstream scg_lib_structs), not the paper.
SPACER = "GAACAGAATT"
T7_PROMOTER = "TAATACGACTCACTATA"
T7_START = "GGG"
TRANSPOSON = nx.ME_RC + SPACER + T7_PROMOTER + T7_START + nx.ME
TRANSPOSON_TOPOLOGY = HairpinTopology(
    "LIANTI hairpin transposon", TRANSPOSON,
    first_stem_start=0, second_stem_start=len(TRANSPOSON) - len(nx.ME),
    stem_length=len(nx.ME),
)
UMI_NT = 8

# Library-prep assignment is also only available in the upstream page.  These bases are
# canonical vendor sequences, imported rather than copied; the upstream specifies 6-nt i7.
I7_NT = 6
NEBNEXT_HAIRPIN = il.NEBNEXT_HAIRPIN_MODEL
NEBNEXT_OPENED = il.NEBNEXT_USER_OPENED


def transposon_segments(*, inferred: bool = True) -> list[Segment]:
    """The single 68-nt hairpin-forming transposon, written 5' to 3'."""
    return [
        seg("ME reverse complement", nx.ME_RC, "me", inferred=inferred),
        seg("spacer", SPACER, inferred=inferred),
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
        seg("spacer, opposite",
            right[len(nx.ME) + len(T7_START) + len(T7_PROMOTER):
                  len(nx.ME) + len(T7_START) + len(T7_PROMOTER) + len(SPACER)],
            inferred=True),
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
        seg("spacer reverse complement", revcomp(SPACER), inferred=True),
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
    """One of the two possible adapter-ligation orientations, P5 to P7'."""
    insert = "X" * insert_nt if insert_nt <= 30 else "XXXXXXXXXXXX...XXXXXXXXXXXX"
    d_a = end_prepared_amplicon(insert_nt).three_prime_overhang
    lib = Construct([
        # P5 overlaps the first four bases of the Read 1 site (ACAC).
        seg("P5", il.P5[:-4], "p5"),
        seg("Read 1 site", NEBNEXT_OPENED.right_arm, "r1"),
        seg("UMI", "N" * UMI_NT, "umi", placeholder=True, inferred=True),
        seg("T7 start", T7_START, "t7", inferred=True),
        seg("mosaic end", nx.ME, "me", inferred=True),
        seg("genomic insert", insert, placeholder=True),
        seg("dA", d_a, inferred=True),
        seg("index-read site", NEBNEXT_OPENED.left_arm
            + il.INDEX1_PRIMER[len(NEBNEXT_OPENED.left_arm):], "r2"),
        seg("i7", "I" * I7_NT, "cbc", placeholder=True, inferred=True),
        seg("P7'", il.P7_RC, "p7"),
    ], name="LIANTI sequencing library, UMI at P5 end")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid LIANTI final library: " + "; ".join(problems))
    return lib


def protocol_graph() -> ProtocolGraph:
    source=("Defining source: https://doi.org/10.1126/science.aak9787; "
            "LIANTI-specific oligo bases reconstructed from scg_lib_structs because the supplement is unavailable")
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
    graph.react("tagment","Insert inward-facing T7-promoter transposons",
                StepKind.TAGMENTATION,
                ((genome,InputRole.TEMPLATE),(transposon,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,tag_waste))
    graph.discard(tag_waste,"fragments without usable paired promoter ends do not yield the modelled library")
    fill_pol=graph.input("gap_fill_polymerase","gap-fill DNA polymerase",MaterialKind.CATALYST)
    dntps=graph.input("gap_fill_dntps","gap-fill dNTPs",MaterialKind.OLIGO)
    filled=graph.product("gap_filled",gap_filled_fragment().name,MaterialKind.MOLECULAR,
                         state=gap_filled_fragment(),inferred=True,source=source)
    graph.react("gap_fill","Fill the 9-nt Tn5 gaps and complete both T7 promoters",
                StepKind.EXTENSION,
                ((tagged,InputRole.TEMPLATE),(dntps,InputRole.STRUCTURAL),
                 (fill_pol,InputRole.CATALYST)),(filled,))
    t7=graph.input("t7_polymerase","T7 RNA polymerase",MaterialKind.CATALYST)
    rna=graph.product("rna_transcript",transcript().name,MaterialKind.MOLECULAR,
                      state=transcript(),inferred=True,source=source)
    graph.react("ivt","Linearly amplify each fragment by T7 transcription",
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
    rna_removal=graph.input("rna_removal","RNA-removal condition",MaterialKind.CONDITION)
    cleaned=graph.product("cleaned_cdna","RNA-free self-primed cDNA",MaterialKind.MOLECULAR,
                          inferred=True,source=source)
    rna_waste=graph.product("rna_waste","removed LIANTI RNA",MaterialKind.WASTE)
    graph.react("remove_rna","Remove the RNA template",StepKind.STRAND_SEPARATION,
                ((self_primed,InputRole.TEMPLATE),(rna_removal,InputRole.CONDITION)),
                (cleaned,rna_waste))
    graph.discard(rna_waste,"RNA is removed before UMI-bearing second-strand synthesis")
    umi_primer=graph.input("umi_primer","UMI-bearing second-strand primer",MaterialKind.OLIGO,
                           state=tuple(second_strand_primer_segments()),inferred=True,source=source)
    second_pol=graph.input("second_strand_polymerase","second-strand DNA polymerase",
                           MaterialKind.CATALYST)
    amplicon=graph.product("umi_amplicon",umi_amplicon().name,MaterialKind.MOLECULAR,
                           state=umi_amplicon(),inferred=True,source=source)
    graph.react("second_strand","Synthesize the UMI-bearing second strand",StepKind.EXTENSION,
                ((cleaned,InputRole.TEMPLATE),(umi_primer,InputRole.STRUCTURAL),
                 (second_pol,InputRole.CATALYST)),(amplicon,))
    repair=graph.input("end_repair_mix","end-repair enzyme mix",MaterialKind.CATALYST)
    repaired=graph.product("repaired_amplicon","end-repaired LIANTI amplicon",
                           MaterialKind.MOLECULAR,state=umi_amplicon(),inferred=True,source=source)
    graph.react("end_repair","End-repair LIANTI amplicons",StepKind.END_REPAIR,
                ((amplicon,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
    tailing=graph.input("tailing_polymerase","dA-tailing polymerase",MaterialKind.CATALYST)
    datp=graph.input("datp","dATP",MaterialKind.OLIGO)
    prepared_scene=end_prep_scene()
    prepared=graph.product("da_tailed","dA-tailed LIANTI amplicon",MaterialKind.MOLECULAR,
                           state=prepared_scene,inferred=True,source=source)
    graph.react("da_tail","Add 3-prime dA overhangs",StepKind.TAILING,
                ((repaired,InputRole.TEMPLATE),(datp,InputRole.STRUCTURAL),
                 (tailing,InputRole.CATALYST)),(prepared,))
    adapter=graph.input("nebnext_hairpin","NEBNext hairpin adapter",MaterialKind.ADAPTER)
    ligase=graph.input("ligase","DNA ligase",MaterialKind.CATALYST)
    ligated=graph.product("hairpin_ligated","NEBNext-hairpin-ligated LIANTI amplicon",
                          MaterialKind.MOLECULAR,inferred=True,source=source)
    geometry=reported_adapter_end(
        "NEBNext hairpin",kind="cohesive",polarity="3-prime",terminal_sequence="T",
        body_sequence=NEBNEXT_HAIRPIN.sequence[:-1],
        statement="The canonical NEBNext hairpin model has the disclosed terminal 3′ T used with dA-tailed inserts.")
    graph.react("adapter_ligation","Ligate the NEBNext hairpin adapter in either orientation",
                StepKind.LIGATION,
                ((prepared,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(ligated,),geometry=geometry)
    user=graph.input("user_enzyme","USER enzyme",MaterialKind.CATALYST)
    opened=graph.product("opened_library","USER-opened adapter-ligated LIANTI molecules",
                         MaterialKind.MOLECULAR,state=final_library(),inferred=True,source=source)
    released=graph.product("released_hairpin","released hairpin loop",MaterialKind.WASTE)
    graph.react("user_digest","Open the dU-containing hairpin",StepKind.DIGESTION,
                ((ligated,InputRole.TEMPLATE),(user,InputRole.CATALYST)),(opened,released))
    graph.discard(released,"the cleaved hairpin loop is not part of the final library")
    forward=graph.input("pcr_forward","NEBNext library PCR forward primer",MaterialKind.OLIGO,
                        inferred=True,source=source+" (exact LIANTI primer assignment secondary)")
    reverse=graph.input("pcr_reverse","indexed NEBNext library PCR reverse primer",
                        MaterialKind.OLIGO,inferred=True,
                        source=source+" (exact LIANTI primer assignment secondary)")
    pcr_pol=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,
                          inferred=True,source=source)
    graph.react("library_pcr","Amplify the USER-opened indexed library",StepKind.PCR,
                ((opened,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(
                    final,reason="The accessible source does not establish the exact LIANTI library-PCR oligo assignment."))
    graph.sequence(library,SequencingPlan(
        library=final,primers=tuple(SEQ_PRIMERS),
        acquired_roles=("Read 1","Index 1 (i7)","Index 2 (i5)","Read 2"),
        read_lengths={"Index 1 (i7)":I7_NT},
        read_uses={"Read 1":"include","Index 1 (i7)":"demultiplex",
                   "Index 2 (i5)":"demultiplex","Read 2":"include"},
        limitation="The 6-nt i7 placement and LIANTI-specific oligos are secondary-source reconstructions."))
    graph.declare_reviewed_scope((
        "genome","hairpin_transposon","tn5","tagment","gap_fill_polymerase",
        "gap_fill_dntps","gap_fill","t7_polymerase","ivt","reverse_transcriptase",
        "self_primed_rt","rna_removal","remove_rna","umi_primer","second_strand_polymerase",
        "second_strand","end_repair_mix","end_repair","tailing_polymerase","datp",
        "da_tail","nebnext_hairpin","ligase","adapter_ligation","user_enzyme",
        "user_digest","pcr_forward","pcr_reverse","pcr_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
