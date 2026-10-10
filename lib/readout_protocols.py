"""Shared role-level models for targeted immune, DNA–RNA and sensor readouts."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from chemdraw import Construct, Segment, feature
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
    return SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT="", FINAL_LIBRARY=final, SEQ_PRIMERS=(),
        FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.endpoint,
        SEQUENCING_UNAVAILABLE=("The defining source establishes adapter and identifier roles, but the page does not promote panel-specific or supplementary oligo sequences into universal protocol constants."),
        workflow=workflow, oligos=lambda: (),
    )
