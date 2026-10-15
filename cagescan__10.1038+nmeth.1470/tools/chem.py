from nascent_protocols import make_spec, role_library

from nascent_protocols import Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, molecule, role
CAGE_OLIGOS = (
    "5′-TAGTCGAACTGAAGGTCTCCAGCAXnrGrGrG-3′ — CAGEscan template-switch primer (Xn is a 3- or 6-base DNA barcode; terminal rGrGrG is RNA)",
    "5′-TAGTCGAACGAAGGTCTCCGAACCGCTCTTCCGATCTN6-3′ — CAGEscan RT primer",
    "5′-TAGTCGAACTGAAGGTCTCCAGC-3′ — CAGEscan semi-suppressive PCR forward primer",
    "5′-TGACGTCGTCTAGTCGAACTGAAGGTCTCCGAACC-3′ — CAGEscan semi-suppressive PCR reverse primer",
    "5′-AATGATACGGCGACCACCGAGATCTACACTAGTCGAACTGAAGG-3′ — CAGEscan forward PCR primer",
    "5′-CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT-3′ — CAGEscan reverse PCR primer",
    "5′-TAGTCGAACTGAAGGTCTCCAGCA-3′ — CAGEscan forward paired-end sequencing primer",
    "5′-CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT-3′ — CAGEscan reverse paired-end sequencing primer",
)

DEFINITION = Definition("CAGEscan — paired-end cap-to-transcript sequencing", "10.1038/nmeth.1470", "Plessy et al., Nature Methods (2010)", "Cap-dependent template switching marks the transcript 5′ end while a random RT primer supplies a downstream mate, linking promoters to transcript bodies without tag cleavage.", "One mate enters at the cap-derived 5′ end and the other enters from the random-primer-defined downstream end.", "01_cagescan.html")

FINAL_LIBRARY = role_library(DEFINITION.title + " final library", (role("cap-linked transcript cDNA"), role("downstream transcript cDNA")),
                             indexed=True, platform="Illumina-compatible")
globals().update(vars(make_spec(DEFINITION, FINAL_LIBRARY,
    unavailable='The original CAGEscan supplement publishes the RT and PCR oligos, but the historical instrument-primer geometry is not asserted as a modern invariant.', oligos=lambda: CAGE_OLIGOS)))

def _cagescan_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input(
        "total_rna", "total RNA containing capped transcripts",
        MaterialKind.BIOLOGICAL,
        state=molecule("total RNA", (role("capped transcript"),)), primary=True)
    n6 = graph.input("n6_rt_primer", "published CAGEscan N6 reverse-transcription primer",
                     MaterialKind.OLIGO)
    reverse_transcriptase = graph.input("reverse_transcriptase", "template-switching reverse transcriptase",
                                       MaterialKind.CATALYST)
    first = graph.product(
        "first_strand", "random-primed first-strand cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("first-strand cDNA", (role("RT handle"),
                                              role("random-primed transcript cDNA"))))
    graph.react("reverse_transcription", "Random-prime reverse transcription",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (n6, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (first,))
    tso = graph.input("cap_tso", "cap-dependent CAGEscan template-switch oligo",
                      MaterialKind.OLIGO)
    switched = graph.product(
        "cap_linked_cdna", "cap-linked CAGEscan cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("cap-linked cDNA", (role("template-switch 5-prime tag"),
                                            role("transcript cDNA"), role("random-primer end"))))
    graph.react("template_switch", "Template-switch at the capped RNA 5-prime end",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (switched,))
    forward = graph.input("pcr_forward", "published CAGEscan forward PCR primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "published CAGEscan reverse PCR primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "paired-end CAGEscan library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "PCR with paired-end Illumina primers", StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The historical PCR oligos are published, but this role-only construct does not assert modern instrument-primer landing geometry."))
    graph.sequence(library, SequencingPlan(limitation=d.endpoint))
    graph.declare_reviewed_scope((
        "total_rna", "n6_rt_primer", "reverse_transcriptase", "reverse_transcription",
        "cap_tso", "template_switch", "pcr_forward", "pcr_reverse", "pcr_polymerase",
        "library_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'cagescan'
    return _cagescan_graph(d, final)
