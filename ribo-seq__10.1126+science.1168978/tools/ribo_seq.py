"""Ingolia ribosome-profiling library with the 2012 indexed run geometry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Row, Scene, Workflow, feature
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

def workflow():
    w=Workflow(MolecularState("translated mRNA",(Row(chunks=[("mRNA — [80S-protected footprint] — mRNA",None,False)]),)))
    w.react("Digest unprotected RNA",[Row(chunks=[("RNase I ↓ retained 28–30-nt footprint","umi",False)])],name="ribosome footprint",note="RNase digestion and monosome isolation select protected fragments.")
    end=Scene(); end.strand("linker",[seg("5′-rApp ligating end",LINKER),seg("3′ ddC block","[3ddC]",placeholder=True)],label="published footprint linker"); end.ligation_end("linker","5′-rApp ligating end",kind="single-strand",polarity="5-prime",owner="adapter",text="5′-rApp single-strand end")
    geometry=AdapterLigationGeometry.known("single-strand",end,statement="The detailed protocol prints /5rApp/-CTGTAGGCACCATCAAT-/3ddC/: its activated 5′ end ligates to the footprint 3′-OH; this is not a duplex overhang.")
    w.ligate_adapter("Ligate the footprint 3′ adapter",geometry,adapter_ligation_scene(fragment_name="ribosome footprint").rows(),name="adapter-ligated footprint",note="The 3′ block prevents linker concatemerization.")
    w.react("Reverse-transcribe and circularize cDNA",circular_cdna_rows([seg("footprint cDNA","X"*30,placeholder=True),seg("published RT-primer handles","X"*24,"r2",placeholder=True)]),name="circular cDNA",note="CircLigase supplies the second PCR boundary without a second RNA-adapter ligation.")
    return w
