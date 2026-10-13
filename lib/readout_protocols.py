"""Shared role-level models for targeted immune, DNA–RNA and sensor readouts."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from chemdraw import Construct, Segment, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, molecule, role, role_library


SAMPLE_INDEX = role("sample index", tag="cbc", length_bp=8,
                    molecular_feature=feature("sample_index", "sample_index", "unknown"))
PLATE_BARCODE = role("plate barcode", tag="cbc",
                     molecular_feature=feature("plate_barcode", "sample_index", "unknown"))
WELL_UDI = role("well unique dual index", tag="cbc",
                molecular_feature=feature("well_udi", "sample_index", "unknown"))
CELL_BARCODE = role("cell barcode", tag="cbc",
                    molecular_feature=feature("cell_barcode", "cell_barcode", "unknown"))
RT_SAMPLE_BARCODE = role("RT sample barcode", tag="cbc",
                         molecular_feature=feature("rt_sample_barcode", "sample_index", "unknown"))
UMI = role("UMI", tag="umi",
           molecular_feature=feature("umi", "umi", "unknown"))
ANALYTE_BARCODE = role("aptamer identity barcode", tag="cbc",
                       molecular_feature=feature("analyte_identity", "feature_barcode", "unknown"))


def _tirtl():
    return chain(molecule("cells distributed across replicate wells", (
        role("T-cell RNA"),)), (
        ("Lyse and reverse-transcribe in each miniaturized well",
         molecule("TCR cDNA", (role("TCR-alpha cDNA"), role("TCR-beta cDNA"))),
         "Alpha/beta pairing is reconstructed from shared occurrence patterns across replicate wells."),
        ("Multiplex-PCR V-to-C segments with partial Nextera-tailed primers",
         molecule("PCR-I products", (role("partial Nextera end"), role("TCR V(D)J insert"),
                                      PLATE_BARCODE, role("partial Nextera end"))),
         "The C-segment reverse primer carries the plate barcode; alpha and beta primer pools are amplified together."),
        ("PCR-II adds full adapters and a well-specific UDI pair",
         molecule("TIRTL-seq library", (role("P5 / Read 1 end"), WELL_UDI,
                                         role("TCR amplicon"), PLATE_BARCODE,
                                         role("Read 2 / P7 end"))),
         "The second PCR makes 384 well-addressable libraries that can be pooled for Illumina sequencing."),
    ))


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


def _smol():
    return chain(molecule("structure-switching aptamer sensor", (
        role("target-binding aptamer"), ANALYTE_BARCODE, role("quenched capture strand"))), (
        ("Bind the cognate small molecule and switch aptamer structure",
         molecule("analyte-bound aptamer", (role("small-molecule analyte"),
                                             role("folded aptamer"), ANALYTE_BARCODE)),
         "Target binding destabilizes the duplex that retained the DNA identity barcode."),
        ("Release the analyte-specific DNA barcode",
         molecule("released reporter", (ANALYTE_BARCODE,)),
         "Reporter abundance is a proxy for the amount of the metabolite recognized by that aptamer."),
        ("PCR-index the released reporters and sequence them",
         molecule("smol-seq library", (role("P5 end"), SAMPLE_INDEX,
                                        ANALYTE_BARCODE, role("P7 end"))),
         "Multiplexed barcode counts convert chemically diverse analytes into one DNA-sequencing readout."),
    ))


@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    citation: str
    summary: str
    endpoint: str
    notes: str


DEFS = {
    "tirtl-seq": Definition("TIRTL-seq — throughput-intensive rapid TCR library sequencing",
        "10.1038/s41592-025-02907-9", "Pogorelyy et al., Nature Methods (2026)",
        "Miniaturized plate RT, multiplex V/C amplification, plate barcodes and well-specific UDIs generate deep alpha/beta TCR repertoires whose chain pairs are inferred across replicate wells.",
        "TCR insert reads are associated with the PCR-I plate barcode and PCR-II well UDI pair.",
        "01_tirtl_seq.html"),
    "sdr-seq": Definition("SDR-seq — single-cell DNA–RNA sequencing",
        "10.1038/s41592-025-02805-0", "Lindenhofer et al., Nature Methods (2025)",
        "In-situ barcoded poly(dT) RT is followed by Tapestri droplet multiplex PCR; distinct R2N and R2 overhangs split targeted gDNA and RNA into matched sequencing libraries.",
        "The gDNA branch spans targeted variants; the RNA branch reports transcript sequence plus cell barcode, RT sample barcode and UMI.",
        "01_sdr_seq.html"),
    "smol-seq": Definition("smol-seq — small-molecule sequencing",
        "10.1038/s41587-025-02554-7", "Tan and Fraser, Nature Biotechnology (2026)",
        "A structure-switching aptamer releases an analyte-specific DNA barcode when it binds its small-molecule target; barcode sequencing converts metabolite abundance into counts.",
        "The sequenced feature barcode identifies the aptamer—and therefore the analyte—not the small molecule itself.",
        "01_smol_seq.html"),
}


def _authority(d: Definition) -> str:
    return f"{d.citation}; doi:{d.doi}"


def _unavailable(final: Construct, reason: str) -> PCRGeometry:
    return PCRGeometry.unavailable(final, reason=reason)


def _tirtl_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = _authority(d)
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input(
        "replicate_well_cells", "T cells distributed across replicate wells",
        MaterialKind.BIOLOGICAL,
        state=molecule("cells distributed across replicate wells", (role("T-cell RNA"),)),
        primary=True)
    lysis = graph.input("lysis_condition", "miniaturized-well lysis condition",
                        MaterialKind.CONDITION)
    rt_primers = graph.input("rt_primer_pool", "TCR reverse-transcription primer pool",
                             MaterialKind.OLIGO)
    reverse_transcriptase = graph.input("reverse_transcriptase", "reverse transcriptase",
                                       MaterialKind.CATALYST)
    cdna = graph.product(
        "tcr_cdna", "TCR alpha/beta cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("TCR cDNA", (role("TCR-alpha cDNA"), role("TCR-beta cDNA"))))
    graph.react(
        "reverse_transcription", "Lyse and reverse-transcribe each replicate well",
        StepKind.REVERSE_TRANSCRIPTION,
        ((cells, InputRole.TEMPLATE), (lysis, InputRole.CONDITION),
         (rt_primers, InputRole.STRUCTURAL),
         (reverse_transcriptase, InputRole.CATALYST)), (cdna,))

    vc_forward = graph.input("vc_forward_pool", "Nextera-tailed V-segment primer pool",
                             MaterialKind.OLIGO)
    vc_reverse = graph.input("vc_reverse_pool", "plate-barcoded C-segment primer pool",
                             MaterialKind.OLIGO)
    pcr1_polymerase = graph.input("pcr1_polymerase", "first-PCR polymerase",
                                  MaterialKind.CATALYST)
    pcr1 = graph.product(
        "pcr1_product", "plate-barcoded TCR V-to-C amplicons",
        MaterialKind.MOLECULAR,
        state=molecule("PCR-I products", (role("partial Nextera end"),
            role("TCR V(D)J insert"), PLATE_BARCODE, role("partial Nextera end"))))
    graph.react(
        "pcr1", "Multiplex-PCR TCR V-to-C segments",
        StepKind.PCR,
        ((cdna, InputRole.TEMPLATE), (vc_forward, InputRole.STRUCTURAL),
         (vc_reverse, InputRole.STRUCTURAL), (pcr1_polymerase, InputRole.CATALYST)),
        (pcr1,), geometry=_unavailable(
            pcr1.state,
            "V/C primer identities depend on the receptor segment and the source does not define one universal amplicon."))

    pcr2_forward = graph.input("pcr2_forward", "full-adapter PCR-II forward primer",
                               MaterialKind.OLIGO)
    pcr2_reverse = graph.input("pcr2_reverse", "well-UDI PCR-II reverse primer",
                               MaterialKind.OLIGO)
    pcr2_polymerase = graph.input("pcr2_polymerase", "second-PCR polymerase",
                                  MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "well-indexed TIRTL-seq library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(
        "pcr2", "Add full adapters and a well-specific UDI pair",
        StepKind.PCR,
        ((pcr1, InputRole.TEMPLATE), (pcr2_forward, InputRole.STRUCTURAL),
         (pcr2_reverse, InputRole.STRUCTURAL), (pcr2_polymerase, InputRole.CATALYST)),
        (library,), geometry=_unavailable(
            final, "The page records the UDI and adapter roles but not a universal fixed primer pair."))
    graph.sequence(library, SequencingPlan(limitation=d.endpoint))
    graph.declare_reviewed_scope((
        "replicate_well_cells", "lysis_condition", "rt_primer_pool",
        "reverse_transcriptase", "reverse_transcription", "vc_forward_pool",
        "vc_reverse_pool", "pcr1_polymerase", "pcr1", "pcr2_forward",
        "pcr2_reverse", "pcr2_polymerase", "pcr2", "sequencing_library",
    ), source)
    return graph.require_publishable()


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


def _smol_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = _authority(d)
    graph = ProtocolGraph(d.title, authority=source)
    sensor = graph.input(
        "aptamer_sensor", "structure-switching aptamer sensor",
        MaterialKind.MOLECULAR,
        state=molecule("structure-switching aptamer sensor", (
            role("target-binding aptamer"), ANALYTE_BARCODE, role("quenched capture strand"))),
        primary=True)
    analyte = graph.input("small_molecule_analyte", "cognate small-molecule analyte",
                          MaterialKind.MOLECULAR)
    bound = graph.product(
        "bound_sensor", "analyte-bound switched aptamer", MaterialKind.MOLECULAR,
        state=molecule("analyte-bound aptamer", (role("small-molecule analyte"),
            role("folded aptamer"), ANALYTE_BARCODE)))
    unbound = graph.product("unbound_sensor", "unbound aptamer sensor", MaterialKind.WASTE)
    graph.react(
        "analyte_binding", "Bind cognate analyte and switch aptamer structure",
        StepKind.BINDING,
        ((sensor, InputRole.TEMPLATE), (analyte, InputRole.ANALYTE)), (bound, unbound))
    graph.discard(unbound, "does not release an analyte-responsive reporter")

    release = graph.input("release_condition", "reporter-release condition",
                          MaterialKind.CONDITION)
    reporter = graph.product(
        "released_reporter", "released analyte-identity reporter", MaterialKind.MOLECULAR,
        state=molecule("released reporter", (ANALYTE_BARCODE,)))
    spent = graph.product("spent_sensor", "spent aptamer/capture-strand complex",
                          MaterialKind.WASTE)
    graph.react(
        "reporter_release", "Release the analyte-specific DNA reporter",
        StepKind.STRAND_SEPARATION,
        ((bound, InputRole.TEMPLATE), (release, InputRole.CONDITION)), (reporter, spent))
    graph.discard(spent, "the released DNA identity reporter is the sequenced analyte proxy")

    forward = graph.input("index_forward", "reporter index-PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("index_reverse", "reporter index-PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "reporter PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "indexed smol-seq reporter library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(
        "index_pcr", "PCR-index released reporter molecules", StepKind.PCR,
        ((reporter, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
         (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
        (library,), geometry=_unavailable(
            final, "The aptamer panel determines the reporter family; no one reporter sequence is universal."))
    graph.sequence(library, SequencingPlan(limitation=d.endpoint))
    graph.declare_reviewed_scope((
        "aptamer_sensor", "small_molecule_analyte", "analyte_binding",
        "release_condition", "reporter_release", "index_forward", "index_reverse",
        "pcr_polymerase", "index_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()


def get(key: str):
    d = DEFS[key]
    if key == "tirtl-seq":
        final = role_library("TIRTL-seq final library", (WELL_UDI, role("TCR V(D)J amplicon"), PLATE_BARCODE))
        workflow = _tirtl
    elif key == "sdr-seq":
        final = role_library("SDR-seq matched library roles", (CELL_BARCODE, UMI,
            RT_SAMPLE_BARCODE, role("targeted gDNA or cDNA insert")))
        workflow = _sdr
    else:
        final = role_library("smol-seq final library", (ANALYTE_BARCODE,), index=SAMPLE_INDEX)
        workflow = _smol
    spec = SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT="", FINAL_LIBRARY=final, SEQ_PRIMERS=(),
        FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.endpoint,
        SEQUENCING_UNAVAILABLE=("The defining source establishes adapter and identifier roles, but the page does not promote panel-specific or supplementary oligo sequences into universal protocol constants."),
        workflow=workflow, oligos=lambda: (),
    )
    if key == "tirtl-seq":
        spec.protocol_graph = lambda: _tirtl_graph(d, final)
    elif key == "sdr-seq":
        spec.protocol_graph = lambda: _sdr_graph(d, final)
    else:
        spec.protocol_graph = lambda: _smol_graph(d, final)
    return spec
