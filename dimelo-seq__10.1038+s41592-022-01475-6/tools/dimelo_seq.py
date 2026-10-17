"""Altemose et al. DiMeLo-seq directed-methylation model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from adapter_ligation import AdapterLigationGeometry
from chromatin_epigenetics import dimelo_mark_rows
from protocol_graph import (AnnotationSlot, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, StepKind)

TITLE = "DiMeLo-seq — antibody-directed methylation on long DNA"
NOTES = "01_dimelo-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41592-022-01475-6">Altemose et al., <i>Nature Methods</i> (2022)</a>.'
SUMMARY = "An antibody recruits Protein A–Hia5 to a chromatin target; SAM activation deposits nearby m6A marks that are read directly on unamplified long DNA molecules."
CAVEAT = "The scientific construct is the marked native DNA molecule. Oxford Nanopore library-adapter bases are kit-specific and are intentionally represented by role rather than invented sequence."
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_MECHANISM = "nanopore_motor"
MOTOR_ADAPTER_NAME = "ligation-sequencing motor adapter"
SEQUENCING_ENDING = "Extract high-molecular-weight native DNA, attach the chosen Oxford Nanopore ligation-sequencing motor adapter, and feed one strand through the pore. Modified-base calling reports exogenous m6A together with endogenous CpG methylation on the same molecule; no sequencing primer binds."


METHODS = "Altemose et al. 2022 Methods (PMC9189060 JATS 'DiMeLo-seq', 'Nanopore library preparation and sequencing')"
PROTOCOLSIO = "author protocols.io DiMeLo-seq v2 (dx.doi.org/10.17504/protocols.io.b2u8qezw) materials list and steps 10-14"


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    cells=graph.input("cells","1-5 million live, frozen or lightly fixed cells per condition",
                      MaterialKind.BIOLOGICAL,primary=True)
    dig_wash=graph.input("dig_wash","Dig-Wash: 0.02% digitonin, 20 mM HEPES-KOH pH 7.5, 150 mM NaCl, 0.5 mM spermidine, protease inhibitor, 0.1% BSA; 5 min on ice",
                         MaterialKind.CONDITION,source=METHODS)
    nuclei=graph.product("nuclei","permeabilized nuclei with native chromatin",MaterialKind.BIOLOGICAL)
    graph.react("permeabilize","Permeabilize nuclei in digitonin Dig-Wash",StepKind.PERMEABILIZATION,
                ((cells,InputRole.TEMPLATE),(dig_wash,InputRole.CONDITION)),(nuclei,),
                source=METHODS,changes_molecule=False)
    antibody=graph.input("primary_antibody","Protein A-compatible primary antibody, 1:50 in Tween-Wash, 4 °C 2 h",
                         MaterialKind.MOLECULAR,source=METHODS)
    bound=graph.product("antibody_bound","antibody-bound nuclei",MaterialKind.BIOLOGICAL)
    free_ab=graph.product("unbound_antibody","unbound antibody removed in two Tween-Wash washes",MaterialKind.WASTE)
    graph.react("bind_antibody","Bind primary antibody, then wash twice",StepKind.BINDING,
                ((nuclei,InputRole.TEMPLATE),(antibody,InputRole.STRUCTURAL)),(bound,free_ab),
                source=METHODS,changes_molecule=False)
    graph.discard(free_ab,"supernatant of the Tween-Wash washes")
    hia5=graph.input("pa_hia5","200 nM Protein A-Hia5 in Tween-Wash (v1: RT 1 h; v2: 4 °C 2 h)",
                     MaterialKind.CATALYST,source=METHODS)
    tethered=graph.product("tethered","antibody-targeted pA-Hia5 chromatin",
                           MaterialKind.MOLECULAR)
    free_hia5=graph.product("unbound_pa_hia5","unbound pA-Hia5 removed in two Tween-Wash washes",MaterialKind.WASTE)
    graph.react("target_hia5","Recruit pA-Hia5 to antibody-bound chromatin, then wash twice",
                StepKind.GENERIC,
                ((bound,InputRole.TEMPLATE),(hia5,InputRole.CATALYST)),(tethered,free_hia5),
                source=METHODS,changes_molecule=False)
    graph.discard(free_hia5,"supernatant of the Tween-Wash washes")
    sam=graph.input("sam","Activation Buffer with 800 µM SAM, 37 °C (v1: 30 min; v2: 2 h, SAM replenished halfway, 0.05 mM spermidine)",
                    MaterialKind.CONDITION,source=METHODS)
    marked=graph.product("marked_chromatin","m6A-marked native chromatin",
                         MaterialKind.MOLECULAR,state=tuple(dimelo_mark_rows()))
    graph.react("directed_methylation","Deposit m6A near the antibody-bound target",
                StepKind.CHEMICAL_CONVERSION,
                ((tethered,InputRole.TEMPLATE),(hia5,InputRole.CATALYST),
                 (sam,InputRole.CONDITION)),(marked,),source=METHODS)
    extraction=graph.input("extraction_condition","NEB Monarch Genomic DNA Purification (T3010S) or Monarch HMW DNA Extraction (T3050L) kit",
                           MaterialKind.CONDITION,source=METHODS)
    native_dna=graph.product("native_dna","long native m6A-marked DNA",
                             MaterialKind.MOLECULAR)
    graph.react("extract_dna","Extract marked genomic DNA without amplification",StepKind.EXTRACTION,
                ((marked,InputRole.TEMPLATE),(extraction,InputRole.CONDITION)),
                (native_dna,),source=METHODS)
    end_prep=graph.input("end_prep_module","NEBNext End Repair/dA-tailing module (E7546S); 3 µg DNA input",
                         MaterialKind.CATALYST,source=PROTOCOLSIO+"; "+METHODS)
    prepared=graph.product("end_prepared","end-prepared native DNA (ends per the ONT ligation-kit documentation)",
                           MaterialKind.MOLECULAR)
    graph.react("end_prep","End-prepare DNA (LSK110 branch: 1 h plus 30 min deactivation)",StepKind.END_REPAIR,
                ((native_dna,InputRole.TEMPLATE),(end_prep,InputRole.CATALYST)),(prepared,),
                source=METHODS+"; "+PROTOCOLSIO)
    sre=graph.input("sre_buffer","Circulomics SRE buffer pelleting, 70% ethanol washes",MaterialKind.SELECTION,source=METHODS)
    prepared_clean=graph.product("end_prepared_clean","pelleted end-prepared DNA",MaterialKind.MOLECULAR)
    sre_waste=graph.product("end_prep_supernatant","end-prep reaction supernatant and ethanol washes",MaterialKind.WASTE)
    graph.react("end_prep_cleanup","Pellet DNA with SRE buffer (LSK110 branch)",StepKind.CLEANUP,
                ((prepared,InputRole.TEMPLATE),(sre,InputRole.SELECTION)),(prepared_clean,sre_waste),source=METHODS)
    graph.discard(sre_waste,"supernatant after centrifugation at 10,000 g")
    adapter=graph.input("ont_adapter","motor-loaded Oxford Nanopore ligation adapter (SQK-LSK110; SQK-LSK109 for barcoded runs)",
                        MaterialKind.ADAPTER,inferred=True,
                        source=SOURCE+" (kit-specific terminal bases proprietary)")
    ligase=graph.input("adapter_ligase","NEBNext Quick Ligation Module (E6056S), half-volume 50 µl reaction, 1 h",
                       MaterialKind.CATALYST,source=PROTOCOLSIO+"; "+METHODS)
    library=graph.product("library","motor-adapted marked native-DNA library",
                          MaterialKind.LIBRARY)
    geometry=AdapterLigationGeometry.unavailable(
        "The chosen Oxford Nanopore kit's motor-adapter bases and ligation-facing end are proprietary and version-specific; the ONT kit documentation was not read in this audit.")
    graph.react("attach_motor_adapter","Ligate the motor-loaded Nanopore adapter",
                StepKind.LIGATION,
                ((prepared_clean,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(library,),geometry=geometry)
    lfb=graph.input("lfb_wash","centrifugal pelleting with two LFB washes",MaterialKind.SELECTION,source=METHODS)
    library_clean=graph.product("library_clean","pelleted adapted library, resuspended in EB (48 h at 4 °C)",MaterialKind.LIBRARY)
    lfb_waste=graph.product("ligation_supernatant","ligation supernatant and LFB washes",MaterialKind.WASTE)
    graph.react("library_cleanup","Pellet and wash the adapted library with LFB",StepKind.CLEANUP,
                ((library,InputRole.TEMPLATE),(lfb,InputRole.SELECTION)),(library_clean,lfb_waste),source=METHODS)
    graph.discard(lfb_waste,"supernatant after centrifugation at 10,000 g")
    graph.sequence(library_clean,SequencingPlan(
        limitation="MinION v9.4 flow cell (FLO-MIN106.1), 500 ng loaded with EXP-WSH004 wash and reload every 24 h (LSK110 branch); nanopore translocation uses a motor-loaded adapter and no sequencing primer; the proprietary adapter sequence is unavailable."))
    graph.annotate(Highlight("end_prep","Library branches named in the Methods",
        "The drawn path is the stated LSK110 method (2) with DNA pelleting cleanups, used for CTCF, H3K9me3 and their controls. "
        "Method (1), used for optimization and CENP-A/AlphaHOR-RES samples, follows SQK-LSK109 with EXP-NBD104/114 native barcode ligation "
        "(1 µg end-repaired DNA per barcode ligation, ligations at least 20 min, 3 µg pooled into the final ligation, LFB cleanup, 1 µg loaded). "
        "Native barcode bases and the PacBio SMRTbell Express 2.0/Sequel II comparison branch are not drawn.",
        slot=AnnotationSlot.BEFORE_WORKFLOW))
    graph.declare_reviewed_scope((
        "cells","dig_wash","permeabilize","primary_antibody","bind_antibody","pa_hia5",
        "target_hia5","sam","directed_methylation","extraction_condition","extract_dna",
        "end_prep_module","end_prep","sre_buffer","end_prep_cleanup","ont_adapter",
        "adapter_ligase","attach_motor_adapter","lfb_wash","library_cleanup","library_clean"),SOURCE)
    return graph.require_publishable()
