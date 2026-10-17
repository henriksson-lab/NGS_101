"""Source-owned protocol definition and reaction graph.

Authority: Meyer, Nat. Methods 16, 1275 (2019), doi:10.1038/s41592-019-0570-0
(PMC6884681 Methods "Next-generation sequencing" and "Data analysis"), and its
Supplementary Protocol (Springer MOESM10 = PMC NIHMS1537259-supplement-2, same md5).
The DART-seq chemistry is in-cell APOBEC1-YTH deamination; library construction is a
commercial NEBNext kit whose internal oligos and steps these sources do not print.
"""
from __future__ import annotations

from batch_ngs import truseq_library
from chemdraw import Segment
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role
from protocol_specs import Definition, get

DEFINITION = Definition(
    "DART-seq — antibody-free m6A mapping", "10.1038/s41592-019-0570-0",
    "Meyer, Nature Methods (2019)",
    "APOBEC1 fused to the YTH m6A-binding domain of YTHDF2 is expressed in cells and "
    "writes nearby C-to-U marks; total RNA is then made into a standard directional "
    "RNA-seq library, where the marks appear as C-to-T mismatches.",
    "truseq", ("cDNA carrying m6A-adjacent C-to-T mark",),
    "Typically single-read 50-bp Illumina HiSeq 4000 sequencing (Supplementary Protocol C); "
    "reads are reverse-complemented before C-to-U calling (Methods, Data analysis).",
    "RNA")

KIT = "NEBNext Ultra II Directional RNA Library Prep Kit for Illumina"
LOW_INPUT_KIT = "NEBNext Single Cell/Low Input RNA Library Prep Kit"
READ_LENGTHS = {"Read 1": 50}

CAVEAT = (
    "<b>Library adapters are inferred.</b> Methods and the Supplementary Protocol name the "
    f"{KIT} (1 µg total RNA, no rRNA removal) and NEBNext Multiplex Oligos for Illumina, "
    "but print no adapter, index or primer sequence and do not say which index set was used. "
    "The Illumina arms, the i5/i7 windows and their 8-nt lengths are therefore a conventional "
    "TruSeq-compatible placeholder, not a source-specified layout; kit-internal steps "
    "(fragmentation, reverse transcription, second-strand synthesis, end prep, adaptor "
    "ligation) are collapsed into one unresolved step. The low-input branch "
    f"({LOW_INPUT_KIT}, 10 or 100 ng total RNA) and the PacBio Iso-Seq long-read branch "
    "are not drawn. Only Read 1 (SR50, the authors' typical run) is stated; index-read "
    "cycles are not.")


def _final():
    insert = [Segment("cDNA carrying m6A-adjacent C-to-T mark", "X" * 28, placeholder=True)]
    return truseq_library(insert, DEFINITION.title + " final library", inferred_adapters=True)


FINAL_LIBRARY, SEQ_PRIMERS = _final()
SPEC = get("dart-seq", definition=DEFINITION, final_override=(FINAL_LIBRARY, SEQ_PRIMERS))
SPEC.CAVEAT = CAVEAT
SPEC.READ_LENGTHS = READ_LENGTHS


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
                          state=marked.state)
    g.react("extract", "Isolate total RNA", StepKind.EXTRACTION,
            ((marked, InputRole.TEMPLATE), (trizol, InputRole.CONDITION)), (extracted,),
            changes_molecule=False)
    dnase = g.input("dnase", "RNase-free DNase I, 20 min at 37 °C", MaterialKind.CATALYST)
    rna = g.product("total_rna", "DNase-treated total RNA (1 µg into library prep)",
                    MaterialKind.MOLECULAR, state=marked.state)
    dna_waste = g.product("dna_digest", "digested genomic DNA", MaterialKind.WASTE)
    g.react("dnase_step", "Digest contaminating DNA", StepKind.DIGESTION,
            ((extracted, InputRole.TEMPLATE), (dnase, InputRole.CATALYST)), (rna, dna_waste))
    g.discard(dna_waste, "DNase digestion products; not carried as library template")
    kit = g.input("kit", KIT, MaterialKind.CONDITION)
    ligated = g.product("adaptor_cdna", "adaptor-ligated directional cDNA", MaterialKind.MOLECULAR,
                        inferred=True,
                        state=_state("adaptor-ligated cDNA", "NEBNext adaptor end",
                                     "cDNA carrying m6A-adjacent C-to-T mark",
                                     "NEBNext adaptor end"))
    g.react("kit_prep", "INFERRED — kit-internal directional cDNA and adaptor-ligation steps",
            StepKind.MAGIC, ((rna, InputRole.TEMPLATE), (kit, InputRole.CONDITION)), (ligated,),
            reason=f"Methods and the Supplementary Protocol name only the {KIT}; its "
                   "fragmentation, RT, second-strand, end-prep and adaptor reagents and "
                   "geometry are not printed and the kit manual was not inspected in this pass.")
    fwd = g.input("pcr_universal", "NEBNext Universal PCR primer (sequence not printed)",
                  MaterialKind.OLIGO, inferred=True)
    rev = g.input("pcr_index", "NEBNext Multiplex Oligos index primer (set not stated)",
                  MaterialKind.OLIGO, inferred=True)
    pol = g.input("pcr_polymerase", "kit PCR master mix", MaterialKind.CATALYST)
    library = g.product("library", "DART-seq indexed RNA-seq library", MaterialKind.LIBRARY,
                        state=FINAL_LIBRARY, inferred=True)
    g.react("index_pcr", "INFERRED — barcode samples by PCR with NEBNext Multiplex Oligos", StepKind.PCR,
            ((ligated, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL),
             (rev, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (library,),
            geometry=PCRGeometry.unavailable(
                FINAL_LIBRARY, reason="The sources name NEBNext Multiplex Oligos but print no "
                                      "primer sequence or index set."))
    g.sequence(library, SequencingPlan(library=FINAL_LIBRARY, primers=SEQ_PRIMERS,
                                       acquired_roles=("Read 1",),
                                       read_lengths=READ_LENGTHS))
    reviewed = ["cells", "fusion_plasmid", "fugene", "transfect", "growth", "deaminate",
                "trizol", "extract", "dnase", "dnase_step", "kit", "kit_prep",
                "pcr_universal", "pcr_index", "pcr_polymerase", "index_pcr", "library"]
    g.declare_reviewed_scope(reviewed, source)
    return g.require_publishable()


SPEC.protocol_graph = protocol_graph
