from types import SimpleNamespace
from chemdraw import feature
from protocol_model import chain, role_library

from readout_protocols import Construct, Definition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, _authority, _unavailable, molecule, role
DEFINITION = Definition("SDR-seq — single-cell DNA–RNA sequencing",
        "10.1038/s41592-025-02805-0", "Lindenhofer et al., Nature Methods (2025)",
        "In-situ barcoded poly(dT) RT is followed by Tapestri droplet multiplex PCR; distinct R2N and R2 overhangs split targeted gDNA and RNA into matched sequencing libraries.",
        "The gDNA branch spans targeted variants; the RNA branch reports transcript sequence plus cell barcode, RT sample barcode and UMI.",
        "01_sdr_seq.html")


def _sdr_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = _authority(d)
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input(
        "fixed_cells", "fixed permeabilized cells containing genomic DNA and poly(A) RNA",
        MaterialKind.BIOLOGICAL,
        state=molecule("fixed permeabilized cell", (role("genomic DNA"), role("poly(A) RNA"))),
        primary=True)
    rt_primer = graph.input("barcoded_rt_primer", "UMI/sample-barcoded poly(dT) RT primer",
                            MaterialKind.OLIGO)
    reverse_transcriptase = graph.input("reverse_transcriptase", "reverse transcriptase",
                                       MaterialKind.CATALYST)
    tagged = graph.product(
        "tagged_cells", "cells retaining gDNA and tagged cDNA", MaterialKind.MOLECULAR,
        state=molecule("cell retaining gDNA and tagged cDNA", (role("genomic DNA"), UMI,
            RT_SAMPLE_BARCODE, role("capture sequence"), role("target cDNA"))))
    graph.react(
        "in_situ_rt", "Reverse-transcribe RNA in situ with barcoded poly(dT)",
        StepKind.REVERSE_TRANSCRIPTION,
        ((cells, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
         (reverse_transcriptase, InputRole.CATALYST)), (tagged,))

    droplets = graph.input("droplet_condition", "Tapestri droplet partitioning condition",
                           MaterialKind.CONDITION)
    occupied = graph.product("occupied_droplets", "cell-containing Tapestri droplets",
                             MaterialKind.MOLECULAR, state=tagged.state)
    empty = graph.product("empty_droplets", "empty or unusable droplets", MaterialKind.WASTE)
    graph.react("partition", "Partition cells with barcoding beads", StepKind.PARTITION,
                ((tagged, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (occupied, empty))
    graph.discard(empty, "does not contain a usable cell template")

    target_forward = graph.input("target_forward_pool", "multiplex target-forward primer pool",
                                 MaterialKind.OLIGO)
    branch_reverse = graph.input(
        "branch_reverse_pool", "gDNA R2N and RNA R2 reverse-primer pool",
        MaterialKind.OLIGO)
    bead_oligos = graph.input("barcode_bead_oligos", "cell-barcode capture-sequence bead oligos",
                              MaterialKind.OLIGO)
    pcr_polymerase = graph.input("droplet_pcr_polymerase", "droplet PCR polymerase",
                                 MaterialKind.CATALYST)
    dna_amplicons = graph.product(
        "dna_amplicons", "cell-barcoded targeted gDNA amplicons", MaterialKind.MOLECULAR,
        state=molecule("gDNA target amplicons", (CELL_BARCODE, role("capture sequence"),
            role("gDNA target"), role("R2N overhang"))))
    rna_amplicons = graph.product(
        "rna_amplicons", "cell-barcoded targeted cDNA amplicons", MaterialKind.MOLECULAR,
        state=molecule("cDNA target amplicons", (CELL_BARCODE, role("capture sequence"),
            role("cDNA target"), UMI, RT_SAMPLE_BARCODE, role("R2 overhang"))))
    graph.react(
        "droplet_multiplex_pcr", "Multiplex-amplify targeted gDNA and cDNA in droplets",
        StepKind.PCR,
        ((occupied, InputRole.TEMPLATE), (target_forward, InputRole.STRUCTURAL),
         (branch_reverse, InputRole.STRUCTURAL), (bead_oligos, InputRole.STRUCTURAL),
         (pcr_polymerase, InputRole.CATALYST)), (dna_amplicons, rna_amplicons),
        geometry=_unavailable(
            final, "The panel-specific target primers and commercial bead oligo are not one fixed public amplicon."))

    index_polymerase = graph.input("index_polymerase", "library indexing PCR polymerase",
                                   MaterialKind.CATALYST)
    dna_forward = graph.input("dna_index_forward", "gDNA-library index forward primer",
                              MaterialKind.OLIGO)
    dna_reverse = graph.input("dna_index_reverse", "R2N-library index reverse primer",
                              MaterialKind.OLIGO)
    dna_library = graph.product("dna_library", "SDR-seq targeted DNA library",
                                MaterialKind.LIBRARY, state=final)
    graph.react(
        "dna_library_pcr", "Build the R2N-marked DNA library", StepKind.PCR,
        ((dna_amplicons, InputRole.TEMPLATE), (dna_forward, InputRole.STRUCTURAL),
         (dna_reverse, InputRole.STRUCTURAL), (index_polymerase, InputRole.CATALYST)),
        (dna_library,), geometry=_unavailable(
            final, "The public role-level description does not establish a universal R2N indexing pair."))
    rna_forward = graph.input("rna_index_forward", "RNA-library index forward primer",
                              MaterialKind.OLIGO)
    rna_reverse = graph.input("rna_index_reverse", "R2-library index reverse primer",
                              MaterialKind.OLIGO)
    rna_library = graph.product("rna_library", "SDR-seq targeted RNA library",
                                MaterialKind.LIBRARY, state=final)
    graph.react(
        "rna_library_pcr", "Build the R2-marked RNA library", StepKind.PCR,
        ((rna_amplicons, InputRole.TEMPLATE), (rna_forward, InputRole.STRUCTURAL),
         (rna_reverse, InputRole.STRUCTURAL), (index_polymerase, InputRole.CATALYST)),
        (rna_library,), geometry=_unavailable(
            final, "The public role-level description does not establish a universal R2 indexing pair."))
    graph.sequence(dna_library, SequencingPlan(
        limitation="DNA branch: the genomic-target read is acquired together with its cell barcode; exact panel-dependent cycles are not universal."))
    graph.sequence(rna_library, SequencingPlan(
        limitation="RNA branch: transcript, cell barcode, RT sample barcode and UMI are acquired; exact panel-dependent cycles are not universal."))
    graph.declare_reviewed_scope((
        "fixed_cells", "barcoded_rt_primer", "reverse_transcriptase", "in_situ_rt",
        "droplet_condition", "partition", "target_forward_pool", "branch_reverse_pool",
        "barcode_bead_oligos", "droplet_pcr_polymerase", "droplet_multiplex_pcr",
        "index_polymerase", "dna_index_forward", "dna_index_reverse", "dna_library_pcr",
        "rna_index_forward", "rna_index_reverse", "rna_library_pcr", "dna_library",
        "rna_library",
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'sdr-seq'
    return _sdr_graph(d, final)

CELL_BARCODE = role("cell barcode", tag="cbc",
                    molecular_feature=feature("cell_barcode", "cell_barcode", "unknown"))
UMI = role("UMI", tag="umi",
           molecular_feature=feature("umi", "umi", "unknown"))
RT_SAMPLE_BARCODE = role("RT sample barcode", tag="cbc",
                         molecular_feature=feature("rt_sample_barcode", "sample_index", "unknown"))

def _sdr():
    return chain(molecule("fixed permeabilized cell", (role("genomic DNA"),
                                                        role("poly(A) RNA"))), (
        ("Reverse-transcribe RNA in situ with barcoded poly(dT)",
         molecule("cell retaining gDNA and tagged cDNA", (role("genomic DNA"), UMI,
                                                           RT_SAMPLE_BARCODE,
                                                           role("capture sequence"),
                                                           role("target cDNA"))),
         "The RT primer installs a UMI and sample barcode before droplet processing."),
        ("Multiplex-amplify targeted gDNA and cDNA in Tapestri droplets",
         molecule("target amplicons", (CELL_BARCODE, role("capture sequence"),
                                        role("gDNA or cDNA target"),
                                        role("R2N or R2 overhang"))),
         "Capture-sequence complementarity copies the bead cell barcode onto both molecular branches."),
        ("Use distinct reverse overhangs to build separate DNA and RNA libraries",
         molecule("matched SDR-seq libraries", (role("gDNA: Nextera-R2N library"),
                                                 CELL_BARCODE,
                                                 role("RNA: TruSeq-R2 library"), UMI,
                                                 RT_SAMPLE_BARCODE)),
         "The DNA read spans the variant target; the RNA read retains cell, sample and UMI identity."),
    ))

d = DEFINITION
FINAL_LIBRARY = role_library("SDR-seq matched library roles", (CELL_BARCODE, UMI, RT_SAMPLE_BARCODE, role("targeted gDNA or cDNA insert")))
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT="", FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=(),
    FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
    SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE=("The defining source establishes adapter and identifier roles, but the page does not promote panel-specific or supplementary oligo sequences into universal protocol constants."),
    workflow=_sdr, oligos=lambda: (),
)
SPEC.protocol_graph = protocol_graph
globals().update(vars(SPEC))
