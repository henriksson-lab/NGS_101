"""Ingolia ribosome-profiling library with the 2012 indexed run geometry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Row, Scene, Workflow, feature, revcomp
from circular import CircularizationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from rna_special import adapter_ligation_scene, circular_cdna_rows

TITLE = "Ribo-seq — ribosome profiling"
NOTES = "01_ribo-seq.html"
SOURCE = 'Defining paper: <a href="https://doi.org/10.1126/science.1168978">Ingolia et al. (2009)</a>; library details: <a href="https://doi.org/10.1038/nprot.2012.086">Ingolia et al. (2012)</a>.'
SUMMARY = "Nuclease leaves ribosome-protected RNA footprints. Size-selected footprints receive a 3′ adapter, are reverse-transcribed and circularized, then PCR adds the Illumina run structure."
CAVEAT = "The six-base sample index is selected from the published indexed reverse-primer table; the model leaves only that position variable."

LINKER = "CTGTAGGCACCATCAAT"
RT_PRIMER = ("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC"
             "CACTCA"
             "TTCAGACGTGTGCTCTTCCGATCTATTGATGGTGCCTACAG")
FORWARD_PCR = il.P5

# The P5 PCR primer and Read 1 site share ACAC physically.  Store the overlap once.
FINAL_LIBRARY = Construct([
    seg("P5-exclusive", il.P5[:-4], "p5"),
    seg("P5 / Read 1 shared ACAC", il.P5[-4:], "r1"),
    seg("Read 1 site after shared ACAC", il.TRUSEQ_READ1[4:], "r1"),
    seg("ribosome footprint cDNA", "X" * 30, placeholder=True),
    seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
    seg("six-base i7 index reverse complement", "I" * 6, "cbc", placeholder=True,
        feature=feature("sample_i7", "sample_index", "whitelist", whitelist="published indexed reverse-primer set")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="single-index ribosome-profiling library")
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"])
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)")):
    raise ValueError("Ribo-seq final library: " + "; ".join(problems))

FINAL_CAPTION = "PCR-linearized product from the circular cDNA: single-end, six-base single-index, with the physical P5/Read 1 ACAC overlap represented once."
SEQUENCING_INTRO = "The protocol explicitly uses the standard Illumina genomic Read 1 and indexing primers. Read 1 enters the footprint; one six-base i7 read identifies the sample."

FOOTPRINT = Construct([seg("ribosome-protected RNA footprint","X"*30,placeholder=True)],
                      name="28–30-nt ribosome footprint")
ADAPTED_RNA = Construct([
    *FOOTPRINT.segments, seg("3-prime linker",LINKER,"r2")],
    name="3-prime-adapter-ligated ribosome footprint")
LINEAR_CDNA = Construct([
    seg("ribosome footprint cDNA","X"*30,placeholder=True),
    seg("RT-primer-derived handles","X"*28,"r2",placeholder=True)],
    name="linear ribosome-footprint cDNA")


def linker_geometry():
    end=Scene()
    end.strand("linker",[seg("5-prime rApp ligating end",LINKER),
                           seg("3-prime ddC block","[3ddC]",placeholder=True)],
               label="published footprint linker")
    end.ligation_end("linker","5-prime rApp ligating end",kind="single-strand",
                     polarity="5-prime",owner="adapter",
                     text="5-prime-rApp single-strand end")
    return AdapterLigationGeometry.known(
        "single-strand",end,
        statement=("The detailed protocol prints /5rApp/-CTGTAGGCACCATCAAT-/3ddC/: "
                   "its activated 5-prime end ligates to the footprint 3-prime hydroxyl."))

def workflow():
    w=Workflow(MolecularState("translated mRNA",(Row(chunks=[("mRNA — [80S-protected footprint] — mRNA",None,False)]),)))
    w.react("Digest unprotected RNA",[Row(chunks=[("RNase I ↓ retained 28–30-nt footprint","umi",False)])],name="ribosome footprint",note="RNase digestion and monosome isolation select protected fragments.")
    w.ligate_adapter("Ligate the footprint 3′ adapter",linker_geometry(),adapter_ligation_scene(fragment_name="ribosome footprint").rows(),name="adapter-ligated footprint",note="The 3′ block prevents linker concatemerization.")
    w.react("Reverse-transcribe and circularize cDNA",circular_cdna_rows([seg("footprint cDNA","X"*30,placeholder=True),seg("published RT-primer handles","X"*24,"r2",placeholder=True)]),name="circular cDNA",note="CircLigase supplies the second PCR boundary without a second RNA-adapter ligation.")
    return w


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    mrna=graph.input("translated_mrna","translated mRNA in ribosome complexes",
                     MaterialKind.BIOLOGICAL,primary=True)
    rnase=graph.input("rnase_i","RNase I",MaterialKind.CATALYST)
    protected=graph.product("protected_fragments","ribosome-protected RNA fragments",
                            MaterialKind.MOLECULAR,state=FOOTPRINT)
    digested=graph.product("digested_rna","unprotected digested RNA",MaterialKind.WASTE)
    graph.react("digest_unprotected","Digest RNA outside ribosome footprints",
                StepKind.DIGESTION,
                ((mrna,InputRole.TEMPLATE),(rnase,InputRole.CATALYST)),
                (protected,digested))
    graph.discard(digested,"RNA not protected by a ribosome is digested")
    selection=graph.input("footprint_selection","monosome recovery and 28–30-nt size selection",
                          MaterialKind.SELECTION)
    footprint=graph.product("footprint","selected 28–30-nt ribosome footprint",
                            MaterialKind.MOLECULAR,state=FOOTPRINT)
    rejected=graph.product("rejected_fragments","non-footprint RNA fragments",MaterialKind.WASTE)
    graph.react("select_footprint","Recover monosomes and select footprint-sized RNA",
                StepKind.SELECTION,
                ((protected,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),
                (footprint,rejected))
    graph.discard(rejected,"does not pass monosome recovery and footprint-size selection")
    linker=graph.input("linker","5-prime-rApp, 3-prime-ddC footprint linker",
                       MaterialKind.ADAPTER,
                       state=Construct([seg("footprint linker",LINKER,"r2")],name="footprint linker"))
    ligase=graph.input("rna_ligase","RNA 3-prime-adapter ligase",MaterialKind.CATALYST)
    adapted=graph.product("adapted_footprint","3-prime-adapter-ligated footprint",
                          MaterialKind.MOLECULAR,state=ADAPTED_RNA)
    graph.react("ligate_linker","Ligate the activated linker to the footprint 3-prime end",
                StepKind.LIGATION,
                ((footprint,InputRole.TEMPLATE),(linker,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(adapted,),geometry=linker_geometry())
    rt_primer=graph.input("rt_primer","published phosphorylated RT primer",
                          MaterialKind.OLIGO,
                          state=Construct([seg("RT primer",RT_PRIMER)],name="Ribo-seq RT primer"))
    reverse_transcriptase=graph.input("reverse_transcriptase","reverse transcriptase",
                                      MaterialKind.CATALYST)
    cdna=graph.product("linear_cdna",LINEAR_CDNA.name,MaterialKind.MOLECULAR,state=LINEAR_CDNA)
    graph.react("reverse_transcribe","Reverse-transcribe the adapter-ligated footprint",
                StepKind.REVERSE_TRANSCRIPTION,
                ((adapted,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (reverse_transcriptase,InputRole.CATALYST)),(cdna,))
    circligase=graph.input("circligase","CircLigase",MaterialKind.CATALYST)
    circle_rows=tuple(circular_cdna_rows(list(LINEAR_CDNA)))
    circle=graph.product("circular_cdna","circular ribosome-footprint cDNA",
                         MaterialKind.MOLECULAR,state=circle_rows)
    graph.react("circularize","Circularize first-strand cDNA",StepKind.CIRCULARIZATION,
                ((cdna,InputRole.TEMPLATE),(circligase,InputRole.CATALYST)),(circle,),
                geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="The detailed protocol circularizes first-strand cDNA with CircLigase."))
    reverse_sequence=il.P7+"N"*6+revcomp(il.INDEX1_PRIMER)
    forward=graph.input("forward_primer","P5 library PCR primer",MaterialKind.OLIGO,
                        state=Construct([seg("P5 PCR primer",FORWARD_PCR)],name="P5 PCR primer"))
    reverse=graph.input("reverse_primer","six-base indexed P7 PCR primer",MaterialKind.OLIGO,
                        state=Construct([seg("indexed P7 PCR primer",reverse_sequence)],name="indexed P7 PCR primer"))
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("library_pcr","Linearize and index circular cDNA by PCR",StepKind.PCR,
                ((circle,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.exact(
                    FINAL_LIBRARY,PCRPrimer("P5 PCR primer",FORWARD_PCR),
                    PCRPrimer("indexed P7 PCR primer",reverse_sequence),
                    statement="The detailed protocol prints the forward and six-base indexed reverse PCR series."))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),
        acquired_roles=("Read 1","Index 1 (i7)")))
    graph.declare_reviewed_scope((
        "translated_mrna","rnase_i","digest_unprotected","footprint_selection",
        "select_footprint","linker","rna_ligase","ligate_linker","rt_primer",
        "reverse_transcriptase","reverse_transcribe","circligase","circularize",
        "forward_primer","reverse_primer","pcr_polymerase","library_pcr","library"),SOURCE)
    return graph.require_publishable()
