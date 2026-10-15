"""Protocol-owned source, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import (Definition, CELL, UMI, SAMPLE, TARGET, GUIDE, RT_WELL, LIG_WELL, PCR_WELL, _finish_pcr, _rt_then_pcr)

DEFINITION = Definition('sc-end5-seq', '10.1101/2021.04.04.438388', 'Moody et al., bioRxiv (2021)', 'Oligo(dT) or random priming followed by 10x template switching retains nucleotide-resolution RNA 5′ ends with cell barcodes and UMIs.', "The insert's genomic 5′ boundary reports a candidate transcription start site; cell barcode and UMI support single-cell counting.", '01_sc_end5_seq.html')
PAYLOAD = (CELL, UMI, role("RNA 5-prime insert"))

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
    return chain(molecule("single-cell RNA in a 10x droplet", (role("poly(A) and non-poly(A) RNA"),)), (
        ("Prime first-strand synthesis with oligo(dT) or random hexamers",
         molecule("nascent first strand", (role("transcript cDNA"),)),
         "sc-end5-dT and sc-end5-rand are explicit priming variants of the named method."),
        ("Template-switch at the RNA 5′ end onto the gel-bead oligo",
         molecule("5-prime-end cDNA", (CELL, UMI, role("template-switch junction"), role("transcript cDNA"))),
         "Barcode and UMI become associated with a nucleotide-resolution transcript 5′ end."),
        ("Amplify and make the sequencing library",
         molecule("sc-end5-seq library", (CELL, UMI, role("RNA 5-prime insert"), SAMPLE)),
         "The genomic end of the insert reports a candidate transcription start site."),
    ))


def oligos():
    return ()


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input(
        "single_cell_rna", "single-cell poly(A) and non-poly(A) RNA",
        MaterialKind.BIOLOGICAL,
        state=molecule("single-cell RNA in a 10x droplet", (
            role("poly(A) and non-poly(A) RNA"),)), primary=True)
    dt_rna = graph.product("dt_input", "RNA assigned to sc-end5-dT",
                           MaterialKind.MOLECULAR, state=rna.state)
    random_rna = graph.product("random_input", "RNA assigned to sc-end5-rand",
                               MaterialKind.MOLECULAR, state=rna.state)
    graph.react("split_priming_variants", "Separate oligo(dT) and random-hexamer variants",
                StepKind.SPLIT, ((rna, InputRole.TEMPLATE),), (dt_rna, random_rna),
                changes_molecule=False)
    reviewed = ["single_cell_rna", "split_priming_variants"]

    def branch(prefix: str, template, primer_label: str, label: str):
        primer = graph.input(prefix + "_primer", primer_label, MaterialKind.OLIGO)
        rt = graph.input(prefix + "_rt_enzyme", "10x template-switching reverse transcriptase",
                         MaterialKind.CATALYST)
        first = graph.product(prefix + "_first_strand", label + " nascent first strand",
                              MaterialKind.MOLECULAR,
                              state=molecule("nascent first strand", (
                                  role("transcript cDNA"),)))
        graph.react(prefix + "_rt", "Prime first-strand synthesis with " + primer_label,
                    StepKind.REVERSE_TRANSCRIPTION,
                    ((template, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                     (rt, InputRole.CATALYST)), (first,))
        bead_tso = graph.input(prefix + "_bead_tso", "cell-barcode/UMI gel-bead template-switch oligo",
                               MaterialKind.OLIGO)
        switched = graph.product(
            prefix + "_switched", label + " 5-prime-end cDNA",
            MaterialKind.MOLECULAR,
            state=molecule("5-prime-end cDNA", (
                CELL, UMI, role("template-switch junction"), role("transcript cDNA"))))
        graph.react(prefix + "_template_switch", "Template-switch at the RNA 5-prime end",
                    StepKind.TEMPLATE_SWITCHING,
                    ((first, InputRole.TEMPLATE), (bead_tso, InputRole.STRUCTURAL),
                     (rt, InputRole.CATALYST)), (switched,))
        reviewed.extend((prefix + "_primer", prefix + "_rt_enzyme", prefix + "_rt",
                         prefix + "_bead_tso", prefix + "_template_switch"))
        _finish_pcr(graph, prefix, switched, final, label + ": " + d.endpoint,
                    reviewed, label=label)

    branch("dt", dt_rna, "oligo(dT) RT primer", "sc-end5-dT")
    branch("random", random_rna, "random-hexamer RT primer", "sc-end5-rand")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
