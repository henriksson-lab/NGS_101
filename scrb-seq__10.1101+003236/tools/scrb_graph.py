"""Source-owned reaction graphs for SCRB-seq and its mcSCRB successor."""
from __future__ import annotations

import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]

import scrb as S
import illumina as il
import nextera as nx
from chemdraw import Construct, Scene, Segment, complement_segments
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

PAPERS = {
    S.SCRB: ("10.1101/003236", "SCRB-seq"),
    S.MCSCRB: ("10.1038/s41467-018-05347-6", "mcSCRB-seq"),
}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def rt_scenes(protocol: str) -> tuple[Scene, Scene]:
    mrna = [seg("body", "XXXXXXXX...XXXXXXXX", placeholder=True), seg("poly(A)", "A" * S.DT_LEN)]
    primer = [seg("R1 handle", S.FULL_HANDLE, "r1"),
              seg("cell barcode", "B" * S.BARCODE_LEN, "cbc", placeholder=True),
              seg("UMI", "N" * S.UMI_LEN, "umi", placeholder=True),
              seg("dT", "T" * S.DT_LEN)]
    first = [*primer, *complement_segments([mrna[0]]), seg("CCC", "CCC", "tso")]
    tso = [seg("PCR handle", S.PCR_HANDLE, "r1"), seg("GGG", "GGG", "tso")]
    p = Scene(); p.strand("mRNA", mrna); p.anneal("barcoded oligo-dT", primer, to="mRNA", pair=("dT", "poly(A)"), mod5="biotin"); p.arrow("barcoded oligo-dT", "reverse transcriptase")
    sw = Scene(); sw.strand("mRNA", mrna); sw.anneal("first-strand cDNA", first, to="mRNA", pair=("dT", "poly(A)"), mod5="biotin"); sw.anneal("TSO", tso, to="first-strand cDNA", pair=("GGG", "CCC"), above=True); sw.arrow("first-strand cDNA", "template switch")
    return p, sw


def graph_for(protocol: str) -> ProtocolGraph:
    if protocol not in PAPERS:
        raise ValueError(f"unknown SCRB protocol: {protocol!r}")
    source = "https://doi.org/" + PAPERS[protocol][0]
    graph = ProtocolGraph(protocol, authority=source)
    prime, switched_state = rt_scenes(protocol)
    mrna = graph.input("mrna", "polyadenylated RNA from one well",
                       MaterialKind.BIOLOGICAL, state=prime, primary=True)
    oligo_dt = graph.input("barcoded_oligo_dt", "barcoded oligo-dT E3V6NEXT",
                           MaterialKind.OLIGO,
                           state=S.FULL_HANDLE + "N" * (S.BARCODE_LEN + S.UMI_LEN)
                                 + "T" * S.DT_LEN + "VN",
                           source=source + "; E3V6NEXT in Supplementary Methods")
    reverse_transcriptase = graph.input("reverse_transcriptase",
                                       "template-switching reverse transcriptase",
                                       MaterialKind.CATALYST)
    first = graph.product("first_strand", "barcoded first-strand cDNA",
                          MaterialKind.MOLECULAR, state=prime)
    graph.react("reverse_transcription",
                "Reverse-transcribe from the barcoded oligo-dT primer",
                StepKind.REVERSE_TRANSCRIPTION,
                ((mrna, InputRole.TEMPLATE),
                 (oligo_dt, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (first,))
    tso = graph.input("tso", "E5V6NEXT template-switching oligo",
                      MaterialKind.OLIGO, state=Construct(S.tso_segments(protocol)))
    switched = graph.product("template_switched_cdna",
                             "template-switched barcoded cDNA",
                             MaterialKind.MOLECULAR, state=switched_state)
    graph.react("template_switch", "Install the opposite preamplification handle",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (switched,))

    peers = graph.input("peer_wells", "barcoded cDNA from the other wells",
                        MaterialKind.MOLECULAR, state=switched_state)
    pooled = graph.product("pooled_cdna", "pooled barcoded cDNA",
                           MaterialKind.MOLECULAR, state=switched_state)
    graph.react("pool", "Pool wells after cell-barcode installation", StepKind.POOL,
                ((switched, InputRole.TEMPLATE), (peers, InputRole.TEMPLATE)),
                (pooled,), changes_molecule=False)
    exonuclease = graph.input("exonuclease_i", "exonuclease I",
                              MaterialKind.CATALYST)
    cleaned = graph.product("cleaned_cdna", "pooled cDNA after primer removal",
                            MaterialKind.MOLECULAR, state=switched_state)
    primer_waste = graph.product("unused_primer_waste",
                                 "digested unused primer products",
                                 MaterialKind.WASTE)
    graph.react("remove_unused_primers", "Remove unused primers with exonuclease I",
                StepKind.CLEANUP,
                ((pooled, InputRole.TEMPLATE),
                 (exonuclease, InputRole.CATALYST)), (cleaned, primer_waste))
    graph.discard(primer_waste, "unused primers are digested before preamplification")

    preamp_primer = PCRPrimer("SINGV6", S.PCR_HANDLE)
    preamp_oligo = graph.input("preamp_primer", "SINGV6 preamplification primer",
                               MaterialKind.OLIGO,
                               state=Construct([seg("SINGV6", S.PCR_HANDLE)]))
    preamp_polymerase = graph.input("preamp_polymerase", "preamplification polymerase",
                                    MaterialKind.CATALYST)
    preamp_state = S.preamp_product()
    amplified = graph.product("amplified_cdna", "amplified full-length cDNA",
                              MaterialKind.MOLECULAR, state=preamp_state)
    graph.react("preamp", "Preamplify full-length cDNA with SINGV6",
                StepKind.PCR,
                ((cleaned, InputRole.TEMPLATE),
                 (preamp_oligo, InputRole.STRUCTURAL),
                 (preamp_polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.symmetric(
                    preamp_state, preamp_primer,
                    statement="SINGV6 binds the shared handle at both cDNA ends."))

    transposome = graph.input("transposome", "s7-loaded Nextera transposome",
                              MaterialKind.ADAPTER,
                              state=Construct([seg("s7", nx.S7, "s7"),
                                               seg("ME", nx.ME, "me")]))
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    enrichment_template = graph.product(
        "enrichment_template", "oligo-dT/s7 tagmentation product",
        MaterialKind.MOLECULAR, state=S.enrichment_template())
    tagment_waste = graph.product(
        "tagment_waste", "other tagmentation topologies",
        MaterialKind.WASTE)
    graph.react("tagment", "Tagment amplified full-length cDNA",
                StepKind.TAGMENTATION,
                ((amplified, InputRole.TEMPLATE),
                 (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)),
                (enrichment_template, tagment_waste))
    graph.discard(tagment_waste,
                  "lacks the selected oligo-dT-handle/s7 primer combination")

    p5_primer = PCRPrimer("P5NEXTPT5", S.P5_ENRICH)
    n7_sequence = il.P7 + "N" * 8 + nx.S7
    n7_primer = PCRPrimer("N7xx", n7_sequence, anneal_nt=len(nx.S7))
    p5 = graph.input("p5_enrichment_primer", "P5NEXTPT5 enrichment primer",
                     MaterialKind.OLIGO,
                     state=Construct([seg("P5NEXTPT5", S.P5_ENRICH)]))
    n7 = graph.input("n7_index_primer", "N7xx indexed enrichment primer",
                     MaterialKind.OLIGO,
                     state=Construct([seg("N7xx", n7_sequence, placeholder=True)]))
    library_polymerase = graph.input("library_polymerase", "library-PCR polymerase",
                                     MaterialKind.CATALYST)
    final = S.final_library(protocol)
    library = graph.product("library", f"{protocol} sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("enrichment_pcr", "Select the 3′ cDNA end and add P5, P7 and i7",
                StepKind.PCR,
                ((enrichment_template, InputRole.TEMPLATE),
                 (p5, InputRole.STRUCTURAL), (n7, InputRole.STRUCTURAL),
                 (library_polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(
                    final, p5_primer, n7_primer,
                    statement="P5NEXTPT5 and N7xx select the oligo-dT/s7 product."))
    lengths = ({"Read 1": 17, "Index 1 (i7)": 8, "Read 2": 34}
               if protocol == S.SCRB else
               {"Read 1": 16, "Index 1 (i7)": 8, "Read 2": 50})
    graph.sequence(library, SequencingPlan(
        library=final, primers=S.READ_PRIMERS,
        acquired_roles=tuple(p.role for p in S.READ_PRIMERS),
        read_lengths=lengths))
    reviewed = [
        "mrna", "barcoded_oligo_dt", "reverse_transcriptase",
        "reverse_transcription", "tso", "template_switch", "peer_wells", "pool",
        "exonuclease_i", "remove_unused_primers", "preamp_primer",
        "preamp_polymerase", "preamp", "transposome", "tn5", "tagment",
        "p5_enrichment_primer", "n7_index_primer", "library_polymerase",
        "enrichment_pcr", "library",
    ]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph() -> ProtocolGraph:
    """The named original SCRB-seq method, for independent audits."""
    return graph_for(S.SCRB)
