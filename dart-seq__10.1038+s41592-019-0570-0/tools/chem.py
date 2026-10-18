"""Source-owned protocol definition and reaction graph.

Authority: Meyer, Nat. Methods 16, 1275 (2019), doi:10.1038/s41592-019-0570-0
(PMC6884681 Methods "Next-generation sequencing" and "Data analysis"), and its
Supplementary Protocol (Springer MOESM10 = PMC NIHMS1537259-supplement-2, same md5).
The DART-seq chemistry is in-cell APOBEC1-YTH deamination; library construction is a
commercial NEBNext kit whose internal oligos and steps these sources do not print.
"""
from __future__ import annotations

from chemdraw import Construct, Scene, Segment
from types import SimpleNamespace
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role
from protocol_graph import Disposition, Terminal

DEFINITION = SimpleNamespace(
    title='DART-seq — antibody-free m6A mapping',
    doi='10.1038/s41592-019-0570-0',
    source_label='Meyer, Nature Methods (2019)',
    summary='APOBEC1 fused to the YTH m6A-binding domain of YTHDF2 is expressed in cells and writes nearby C-to-U marks; total RNA is then made into a standard directional RNA-seq library, where the marks appear as C-to-T mismatches.',
    sequencing='Typically single-read 50-bp Illumina HiSeq 4000 sequencing (Supplementary Protocol C); reads are reverse-complemented before C-to-U calling (Methods, Data analysis).',
)

KIT = "NEBNext Ultra II Directional RNA Library Prep Kit for Illumina"
LOW_INPUT_KIT = "NEBNext Single Cell/Low Input RNA Library Prep Kit"
READ_LENGTHS = {"Read 1": 50}

TITLE=DEFINITION.title
NOTES="01_dart_seq.html"
SOURCE=f'<a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.source_label}</a>.'
SUMMARY=DEFINITION.summary
CAVEAT=("The selected route uses1µg totalRNA and the source-named NEBNext UltraII Directional kit. "
        "Native fragmentation/RT/second-strand/end-prep/ligation/PCR steps and historical adapter/index species "
        "are not reconstructed; a generic TruSeq dual-index shell is not supplied. Low-input10/100ng and "
        "PacBio Iso-Seq are separate source variants outside this selected graph. RNA C→U chemistry is "
        "distinct from strand-dependent cDNA A/T and reference-oriented C→T variant reporting.")
FINAL_LIBRARY=Construct([role("NEBNext directional RNA-seq endpoint; native arms/index geometry unassigned",inferred=True),
                         role("cDNA encoding edited and unedited RNA; strand-specific base orientation unassigned")],
                        name="DART-seq source-limited directional library")
SEQ_PRIMERS=()
FINAL_CAPTION="NEBNext directional library, not a source-exact adapter reconstruction. RNA U is read as a reference-oriented T substitution after strand-aware processing."
SEQUENCING_INTRO=DEFINITION.sequencing
SEQUENCING_ENDING="Source also permits paired-end sequencing; the graph selects typical SR50."
SEQUENCING_UNAVAILABLE="Native sequencing primer/index geometry is unassigned; SR50 acquisition is source-declared."

def oligos():
    return ()

def extraction_mixture():
    sc=Scene()
    sc.strand("RNA",[role("edited and unedited total RNA")],label="RNA population")
    sc.strand("DNA",[role("contaminating genomic DNA")],label="independent genomic DNA")
    return sc


def _state(label, *names):
    return molecule(label, tuple(role(n) for n in names))


def protocol_graph():
    d = DEFINITION
    source = f"{d.source_label}; doi:{d.doi}"
    g = ProtocolGraph(d.title, authority=source)
    cells = g.input("cells", "HEK293T cells (plated ~1:8, transfected 12–24 h later)",
                    MaterialKind.BIOLOGICAL, primary=True,
                    state=_state("HEK293T cell RNA", "cellular RNA containing m6A"))
    plasmid = g.input("fusion_plasmid",
                      "pCMV APOBEC1–YTH plasmid (rat APOBEC1, 15-aa linker, YTHDF2 aa 385–579, HA)",
                      MaterialKind.SELECTION)
    fugene = g.input("fugene", "FuGENE HD transfection reagent", MaterialKind.CONDITION)
    expressing = g.product("expressing", "cells expressing APOBEC1–YTH",
                           MaterialKind.MOLECULAR,
                           state=_state("cellular RNA", "cellular RNA containing m6A"))
    g.react("transfect", "Transfect the APOBEC1–YTH expression plasmid",
            StepKind.GENERIC,
            ((cells, InputRole.TEMPLATE), (plasmid, InputRole.SELECTION),
             (fugene, InputRole.CONDITION)), (expressing,), changes_molecule=False)
    growth = g.input("growth", "24 h at 37 °C, 5% CO2 (DMEM, 10% FBS)", MaterialKind.CONDITION)
    marked = g.product("marked", "total cellular RNA with C-to-U marks beside YTH-bound m6A",
                       MaterialKind.MOLECULAR,
                       state=_state("marked RNA", "RNA", "C-to-U mark near m6A"))
    g.react("deaminate", "APOBEC1–YTH binds m6A and deaminates nearby cytidines in cells",
            StepKind.CHEMICAL_CONVERSION,
            ((expressing, InputRole.TEMPLATE), (growth, InputRole.CONDITION)), (marked,),
            reason="Editing is not a selection: RNA without YTH-bound m6A stays in the total "
                   "RNA and is sequenced; background APOBEC1 sites are removed computationally "
                   "using APOBEC1-alone and APOBEC1–YTHmut controls (Methods, Data analysis).")
    trizol = g.input("trizol", "TRIzol (or RNeasy Plus Mini; Supplementary Protocol B)",
                     MaterialKind.CONDITION)
    extracted = g.product("extracted", "total RNA with genomic DNA", MaterialKind.MOLECULAR,
                          state=extraction_mixture())
    g.react("extract", "Isolate total RNA", StepKind.EXTRACTION,
            ((marked, InputRole.TEMPLATE), (trizol, InputRole.CONDITION)), (extracted,),
            changes_molecule=False)
    dnase = g.input("dnase", "RNase-free DNase I, 20 min at 37 °C", MaterialKind.CATALYST)
    rna = g.product("total_rna", "DNase-treated total RNA (1 µg into library prep)",
                    MaterialKind.MOLECULAR, state=marked.state)
    dna_waste = g.product("dna_digest", "hydrolysis products of contaminating genomic DNA", MaterialKind.MOLECULAR)
    g.react("dnase_step", "Digest contaminating DNA", StepKind.DIGESTION,
            ((extracted, InputRole.TEMPLATE), (dnase, InputRole.CATALYST)), (rna, dna_waste))
    g.mark_terminal(Terminal(dna_waste.id,Disposition.UNRESOLVED,"DNase hydrolysis changes DNA into fragments/nucleotides, not a physical discard. Handling through source extraction/library preparation is not followed here."))
    kit = g.input("kit", KIT+";1µg totalRNA,no rRNA removal required", MaterialKind.CONDITION)
    indexes=g.input("multiplex_oligos","NEBNext Multiplex Oligos for sample barcoding; historical set/species unassigned",MaterialKind.OLIGO)
    library=g.product("library","DART-seq directional indexed library",MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True)
    g.react("kit_prep","Prepare directional RNA-seq library and barcode with source-named NEBNext reagents",StepKind.MAGIC,
        ((rna,InputRole.TEMPLATE),(kit,InputRole.CONDITION),(indexes,InputRole.STRUCTURAL)),(library,),
        reason="The article and Supplementary Protocol delegate library preparation/barcoding to named NEBNext kits. Native molecular intermediates, physical oligos and operation boundaries are not reconstructed here; no terminal PCR is invented.")
    g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,acquired_roles=("Read 1",),read_lengths=READ_LENGTHS,
        limitation=SEQUENCING_UNAVAILABLE+" "+SEQUENCING_ENDING))
    reviewed = ["cells", "fusion_plasmid", "fugene", "transfect", "growth", "deaminate",
                "trizol", "extract", "dnase", "dnase_step", "kit", "kit_prep",
                "multiplex_oligos", "library"]
    g.declare_reviewed_scope(reviewed, source)
    return g.require_publishable()
