"""Protocol-owned source, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import (Definition, CELL, UMI, SAMPLE, TARGET, GUIDE, RT_WELL, LIG_WELL, PCR_WELL, _finish_pcr, _rt_then_pcr)

DEFINITION = Definition('EasySci-RNA', '10.1038/s41588-023-01572-y', 'Sziraki et al., Nature Genetics (2023)', 'Barcoded oligo(dT) and random-hexamer RT, barcode ligation and indexed PCR create a scalable three-round combinatorial single-nucleus RNA library.', 'RT-well, ligation-well and PCR-well indexes jointly identify each nucleus; a UMI identifies each cDNA molecule.', '01_easysci.html')
PAYLOAD = (RT_WELL, LIG_WELL, PCR_WELL, UMI, role("transcript insert"))

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = ""
FINAL_LIBRARY = role_library(TITLE + " final library roles", PAYLOAD, index=SAMPLE)
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = None
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = "The source defines the read roles, but the primer sequence is panel-, kit- or implementation-dependent; no universal binding oligo is asserted."

def workflow():
    return chain(molecule("fixed nuclei distributed into RT wells", (role("nuclear RNA"),)), (
        ("Reverse-transcribe with matched barcoded oligo(dT) and random-hexamer primers",
         molecule("first-index cDNA in nuclei", (RT_WELL, UMI, role("cDNA"))),
         "The paired priming modes recover polyadenylated and broader nuclear RNA from the same RT well."),
        ("Pool, redistribute and ligate the second well barcode",
         molecule("two-round indexed nuclei", (RT_WELL, LIG_WELL, UMI, role("cDNA")),
                  junctions=(("RT-well barcode", "ligation-well barcode", "ligation"),)),
         "The ordered RT-plus-ligation pair narrows cell identity before PCR."),
        ("Pool, make second strand, tagment and add the PCR-well index",
         molecule("EasySci-RNA library", (RT_WELL, LIG_WELL, PCR_WELL, UMI, role("transcript insert"), SAMPLE)),
         "Three combinatorial rounds identify a nucleus; the P5 barcode is read separately from the P7 sample index."),
    ))


def oligos():
    return ()


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input(
        "rt_well_nuclei", "fixed nuclei distributed into reverse-transcription wells",
        MaterialKind.BIOLOGICAL,
        state=molecule("fixed nuclei in RT wells", (role("nuclear RNA"),)), primary=True)
    dt = graph.input("barcoded_dt", "barcoded oligo(dT) RT primer", MaterialKind.OLIGO)
    random = graph.input("barcoded_random", "matched barcoded random-hexamer RT primer",
                         MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    first = graph.product(
        "first_index_cdna", "RT-well-indexed cDNA retained in nuclei",
        MaterialKind.MOLECULAR,
        state=molecule("first-index cDNA in nuclei", (RT_WELL, UMI, role("cDNA"))))
    graph.react("reverse_transcription", "Reverse-transcribe with barcoded dT and random primers",
                StepKind.REVERSE_TRANSCRIPTION,
                ((nuclei, InputRole.TEMPLATE), (dt, InputRole.STRUCTURAL),
                 (random, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (first,))
    repartition = graph.input("repartition_condition", "pool-and-redistribute ligation-well condition",
                              MaterialKind.CONDITION)
    redistributed = graph.product("redistributed_nuclei", "nuclei redistributed into ligation wells",
                                  MaterialKind.MOLECULAR, state=first.state)
    collisions = graph.product("repartition_losses", "empty or collision-prone redistributed wells",
                               MaterialKind.WASTE)
    graph.react("repartition", "Pool and redistribute nuclei into ligation wells",
                StepKind.PARTITION,
                ((first, InputRole.TEMPLATE), (repartition, InputRole.CONDITION)),
                (redistributed, collisions), changes_molecule=False)
    graph.discard(collisions, "does not provide an interpretable combinatorial cell identity")
    lig_barcode = graph.input("ligation_barcode", "second-round ligation-well barcode adapter",
                              MaterialKind.ADAPTER)
    ligase = graph.input("barcode_ligase", "barcode ligase", MaterialKind.CATALYST)
    indexed = graph.product(
        "two_round_nuclei", "two-round-indexed nuclei",
        MaterialKind.MOLECULAR,
        state=molecule("two-round indexed nuclei", (RT_WELL, LIG_WELL, UMI, role("cDNA")),
                       junctions=(("RT-well barcode", "ligation-well barcode", "ligation"),)))
    graph.react("barcode_ligation", "Ligate the second well barcode",
                StepKind.LIGATION,
                ((redistributed, InputRole.TEMPLATE), (lig_barcode, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (indexed,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The combinatorial barcode role is defined, but the complete plate-specific ligation-end family is not a universal sequence."))
    second_primer = graph.input("second_strand_primer", "second-strand synthesis primer",
                                MaterialKind.OLIGO)
    second_polymerase = graph.input("second_strand_polymerase", "second-strand DNA polymerase",
                                    MaterialKind.CATALYST)
    double = graph.product("double_stranded_cdna", "double-stranded two-index cDNA",
                           MaterialKind.MOLECULAR, state=indexed.state)
    graph.react("second_strand", "Synthesize the second cDNA strand",
                StepKind.EXTENSION,
                ((indexed, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL),
                 (second_polymerase, InputRole.CATALYST)), (double,))
    transposome = graph.input("transposome", "loaded cDNA tagmentation adapter",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagmented = graph.product("tagmented_cdna", "tagmented combinatorially indexed cDNA",
                              MaterialKind.MOLECULAR, state=double.state)
    losses = graph.product("tn5_losses", "non-library Tn5 products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment the double-stranded cDNA", StepKind.TAGMENTATION,
                ((double, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagmented, losses))
    graph.discard(losses, "does not retain an amplifiable two-ended tagmentation product")
    reviewed = [
        "rt_well_nuclei", "barcoded_dt", "barcoded_random", "reverse_transcriptase",
        "reverse_transcription", "repartition_condition", "repartition",
        "ligation_barcode", "barcode_ligase", "barcode_ligation",
        "second_strand_primer", "second_strand_polymerase", "second_strand",
        "transposome", "tn5", "tagment",
    ]
    _finish_pcr(graph, "library", tagmented, final, d.endpoint, reviewed,
                label="EasySci-RNA PCR-well-indexed")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
