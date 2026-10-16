"""Falconer et al. Strand-seq parental-template and historical PE library model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, revcomp
from chromatin_epigenetics import strand_selection_rows
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


def seg(name, top, tag=None, **kw):
    return Segment(name, top, tag, **kw)


TITLE = "Strand-seq — parental template-strand sequencing"
NOTES = "01_strand-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.2206">Falconer et al., <i>Nature Methods</i> (2012)</a>.'
SUMMARY = "Cells replicate once in BrdU; after single-cell isolation, Hoechst/UV selectively nicks the BrdU-containing nascent strands so only the parental template strands amplify."
CAVEAT = ("Falconer et al. Online Methods, ‘DNA template strand library construction,’ prints seven N positions in the custom PCR primer but calls the assigned index a fault-tolerant hexamer; ‘Illumina sequencing’ specifies a seven-base index read. The model uses six variable bases and one adjacent fixed base, but the printed N run does not itself resolve that assignment. The paper names Illumina PE adapters and PE&nbsp;1.0; the adapter and standard read-primer bases below come from Illumina’s obsolete paired-end oligo table.")

# Falconer Online Methods. The apparent seven-N run is explicitly described as a
# fault-tolerant hexamer. The seven-cycle index read therefore continues into P7'.
INDEX_NT = 6
PE_ADAPTER_BOTTOM = il.TRUSEQ_READ1
INDEX_READ = "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG"
# The printed index-read primer spans the adapter oligo and four adjacent bases.
# The old Read 2 primer is its reverse complement with the dA-junction base included.
PE_ADAPTER_TOP = INDEX_READ[:-4]
PE_READ2 = revcomp("A" + INDEX_READ)


def adapter_scene() -> Scene:
    sc = Scene()
    sc.strand("bottom", [
        seg("Read 1 fork arm", PE_ADAPTER_BOTTOM[:-13], "r1"),
        seg("PE stem complement", PE_ADAPTER_BOTTOM[-13:-1], "r1"),
        seg("3′ T", PE_ADAPTER_BOTTOM[-1], "r1"),
    ], label="PE adapter oligo")
    sc.anneal("top", [
        seg("PE stem", PE_ADAPTER_TOP[:12], "r2"),
        seg("PE fork arm", PE_ADAPTER_TOP[12:], "r2"),
    ], to="bottom", pair=("PE stem", "PE stem complement"),
              label="PE adapter oligo", mod5="p", unpaired=("PE fork arm",))
    sc.mark("bottom", "3′ T", "ligation overhang")
    return sc


def adapter_geometry() -> AdapterLigationGeometry:
    adapter = adapter_scene()
    adapter.ligation_end("bottom", "3′ T", kind="cohesive", polarity="3-prime",
                         owner="adapter", text="adapter 3′ dT")
    substrate = Scene()
    substrate.strand("fragment", [seg("repaired DNA end", "X" * 18,
                                     placeholder=True), seg("3′ dA", "A")],
                     label="A-tailed fragment")
    substrate.overhang("fragment", "3′ dA", polarity="3-prime",
                       owner="substrate", text="fragment 3′ dA")
    return AdapterLigationGeometry.known(
        "cohesive", adapter, substrate,
        statement="Falconer Online Methods specify fragment end repair and A tailing before Illumina PE adapter ligation; the obsolete Illumina paired-end oligo table prints the adapter 3′ dT.")


def pe1_segments():
    return [seg("P5", il.P5, "p5"),
            seg("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1")]


def indexed_primer_segments():
    return [seg("P7", il.P7, "p7"),
            seg("hexamer index", "N" * INDEX_NT, "cbc", placeholder=True),
            seg("PE Read 2 arm", PE_READ2, "r2")]


SEQ_PRIMERS = (
    sp.custom("Read 1", "obsolete PE Read 1", il.TRUSEQ_READ1,
              'Illumina "Oligonucleotide Sequences for Paired End DNA" (obsolete)'),
    sp.custom("Index 1 (i7)", "Strand-seq custom index primer", INDEX_READ,
              "Falconer et al. 2012 Online Methods",
              "Seven cycles report the six-base index and the adjacent first base of P7'."),
    sp.custom("Read 2", "obsolete PE Read 2", PE_READ2,
              'Illumina "Oligonucleotide Sequences for Paired End DNA" (obsolete)'),
)


def final_library() -> Construct:
    lib = Construct([
        *pe1_segments(),
        seg("retained parental template fragment", "X" * 36, placeholder=True),
        seg("PE Read 2 site", revcomp(PE_READ2), "r2"),
        seg("i7 reverse complement", "N" * INDEX_NT, "cbc", placeholder=True,
            feature=feature("cell_i7", "cell_barcode", "whitelist", whitelist="published hexamer index set")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="Strand-seq library")
    problems = sp.verify(lib, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid Strand-seq library: " + "; ".join(problems))
    return lib


FINAL_LIBRARY = final_library()
FINAL_CAPTION = "PCR-completed, single-index historical paired-end library. The six-base i7 identifies the cell; there is no i5 index."
SEQUENCING_INTRO = "Paired 76-nt reads used the obsolete Illumina PE primers. A third seven-cycle read is reported. The model assigns six cycles to the barcode and one to an adjacent fixed base; see the source discrepancy above."


def sections():
    return [
        ("Label nascent DNA during one S phase", strand_selection_rows()[:2],
         "Semiconservative replication leaves each chromatid with one parental template and one BrdU-containing newly synthesized strand."),
        ("Repair, A-tail and ligate the historical Illumina PE adapter", adapter_geometry().rows(),
         "The paper repairs and A-tails fragmented DNA before ligating PE adapters. The adapter's 3′ dT joins the fragment's 3′ dA overhang; ligation precedes strand destruction."),
        ("Destroy the nascent strand before PCR", strand_selection_rows()[2:],
         "Hoechst sensitization and UV nick BrdU DNA. PE 1.0 and the custom indexed primer then amplify only intact parental templates."),
    ]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    cells=graph.input("cells","cells entering one S phase",MaterialKind.BIOLOGICAL,
                      primary=True)
    brdu=graph.input("brdu","BrdU nucleotide label",MaterialKind.MOLECULAR)
    replication=graph.input("replication_polymerase","cellular DNA replication machinery",
                            MaterialKind.CATALYST)
    replicated=graph.product("replicated","chromatids with BrdU-labelled nascent strands",
                             MaterialKind.MOLECULAR,
                             state=tuple(strand_selection_rows()[:2]))
    graph.react("label_nascent","Replicate once in BrdU",StepKind.EXTENSION,
                ((cells,InputRole.TEMPLATE),(brdu,InputRole.STRUCTURAL),
                 (replication,InputRole.CATALYST)),(replicated,))
    fragmented=graph.product("fragmented","fragmented BrdU-labelled DNA",
                             MaterialKind.MOLECULAR)
    mnase=graph.input("mnase","micrococcal nuclease",MaterialKind.CATALYST)
    graph.react("fragment_dna","Fragment genomic DNA",StepKind.FRAGMENTATION,
                ((replicated,InputRole.TEMPLATE),(mnase,InputRole.CATALYST)),
                (fragmented,))
    repair_enzymes=graph.input("repair_enzymes","DNA end-repair enzymes",
                               MaterialKind.CATALYST)
    repaired=graph.product("repaired","end-repaired BrdU-labelled fragments",
                           MaterialKind.MOLECULAR)
    graph.react("end_repair","Repair fragment ends",StepKind.END_REPAIR,
                ((fragmented,InputRole.TEMPLATE),(repair_enzymes,InputRole.CATALYST)),
                (repaired,))
    tailing_enzymes=graph.input("tailing_enzymes","DNA dA-tailing enzymes",
                                MaterialKind.CATALYST)
    tailed=graph.product("tailed","3′-dA-tailed BrdU-labelled fragments",
                         MaterialKind.MOLECULAR)
    graph.react("a_tail","A-tail fragment ends",StepKind.TAILING,
                ((repaired,InputRole.TEMPLATE),(tailing_enzymes,InputRole.CATALYST)),
                (tailed,))
    adapter=graph.input("pe_adapter","historical Illumina PE adapter",
                        MaterialKind.ADAPTER,state=tuple(adapter_scene().rows()))
    ligase=graph.input("adapter_ligase","PE adapter ligase",MaterialKind.CATALYST)
    adapted=graph.product("adapted","PE-adapter-ligated BrdU-labelled DNA",
                          MaterialKind.MOLECULAR)
    graph.react("ligate_adapter","Ligate the historical Illumina PE adapter",
                StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(adapted,),
                geometry=adapter_geometry())
    damage=graph.input("hoechst_uv","Hoechst and UV strand-destruction condition",
                       MaterialKind.CONDITION)
    parental=graph.product("parental_templates","intact parental template strands",
                           MaterialKind.MOLECULAR)
    nascent=graph.product("damaged_nascent","nicked BrdU-containing nascent strands",
                          MaterialKind.WASTE)
    graph.react("destroy_nascent","Nick BrdU-containing nascent strands selectively",
                StepKind.CHEMICAL_CONVERSION,
                ((adapted,InputRole.TEMPLATE),(damage,InputRole.CONDITION)),
                (parental,nascent))
    graph.discard(nascent,"UV-nicked BrdU strands do not support intact library amplification")
    pe1=graph.input("pe1_primer","Illumina PE 1.0 PCR primer",MaterialKind.OLIGO)
    indexed=graph.input("indexed_primer","custom P7/index/PE2 PCR primer",
                        MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    graph.react("library_pcr","Amplify only intact parental-template molecules",
                StepKind.PCR,
                ((parental,InputRole.TEMPLATE),(pe1,InputRole.STRUCTURAL),
                 (indexed,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The source gives the custom indexed primer and PE 1.0 by name; the historical PE1 PCR-primer base string is not independently transcribed here."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "cells","brdu","replication_polymerase","label_nascent","mnase","fragment_dna",
        "repair_enzymes","end_repair","tailing_enzymes","a_tail","pe_adapter",
        "adapter_ligase","ligate_adapter","hoechst_uv","destroy_nascent",
        "pe1_primer","indexed_primer","pcr_polymerase","library_pcr","library"),SOURCE)
    return graph.require_publishable()
