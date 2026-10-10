"""Li et al. GRID-seq bivalent-linker library."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Scene, Workflow, feature
from rna_dna_bridge import BivalentRNADNABridge

TITLE = "GRID-seq — RNA–chromatin contacts"
NOTES = "01_grid-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nbt.3968">Li et al. (2017)</a>.'
SUMMARY = "A biotinylated bivalent RNA–DNA linker captures a chromatin RNA and a nearby AluI-cut genomic end. MmeI releases a short RNA tag–linker–DNA tag molecule for Illumina sequencing."
CAVEAT = "The tag bases are molecule-specific placeholders; linker bases, modifications and orientation are transcribed from the paper. The five-base PCR barcode is published, but its index-read priming recipe is not."

GRID_R1 = il.TRUSEQ_READ1
GRID_P7_ARM = "AGATCGGAAGAGCACACGTCT"  # printed phosphorylated Y-adapter strand
BRIDGE = BivalentRNADNABridge(
    "GRID-seq bridge", "GTTGGAGTTCGGTGTGTGGGAGTGAGCTGTGTC",
    "GUUGGAUUC", 3, "G", "GACACAGCTCACTCCCACACACCGAACTCCAAC",
    biotin_offset=8, preadenylated_rna_5=True, barcode_id="grid_linker_umi")
LINKER_DNA = BRIDGE.dna_strand

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("five-base library index", "IIIII", "cbc", placeholder=True,
        feature=feature("sample_index_i5", "sample_index", "unknown")),
    seg("Read 1 site", GRID_R1, "r1"),
    seg("RNA-derived tag", "R" * 20, placeholder=True, length_bp=20),
    seg("RNA-side MmeI boundary", "X" * 6, "me", placeholder=True),
    seg("linker RNA-arm fixed bases", BRIDGE.rna_fixed_5.replace("U", "T")),
    seg("three-base molecular barcode", "N" * BRIDGE.barcode_length, "umi", placeholder=True,
        feature=feature("grid_linker_umi", "umi", "random")),
    seg("linker RNA-arm terminal G", BRIDGE.rna_terminal),
    seg("linker duplex fixed bases", BRIDGE.duplex_arm),
    seg("DNA-side MmeI boundary", "X" * 6, "me", placeholder=True),
    seg("genomic-DNA tag", "D" * 20, placeholder=True, length_bp=20),
    seg("P7-side Y-adapter arm", GRID_P7_ARM, "r2"),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="GRID-seq single-end library")
SEQ_PRIMERS = (sp.TRUSEQ["R1"],)
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=("Read 1",)):
    raise ValueError("GRID-seq final library: " + "; ".join(problems))
READ_LENGTHS = {"Read 1": 100}
FINAL_CAPTION = "The 100-nt read crosses an ~20-nt RNA tag, the orientation-defining bivalent linker and an ~20-nt genomic-DNA tag. PCR Primer #1 places a five-base multiplex barcode between P5 and the Read 1 site."
SEQUENCING_INTRO = "The published HiSeq 2500 run prints the Read 1 primer. The linker orientation identifies which flanking MmeI tag came from RNA and which from DNA."
SEQUENCING_UNAVAILABLE = "The paper prints the five-base barcode inside PCR Primer #1 but does not state the index-read recipe or a separate index sequencing primer; its binding geometry is therefore not asserted."

_initial = Scene()
_initial.strand("RNA", [seg("chromatin-associated RNA", "R" * 28, placeholder=True)], label="crosslinked chromatin RNA")
_initial.strand("DNA", [seg("nearby genomic DNA", "D" * 34, placeholder=True)], label="nearby genomic DNA")
_initial.mark("RNA", "chromatin-associated RNA", "fixed in spatial proximity")
INITIAL_ROWS = tuple(_initial.rows())
INITIAL_NAME = "Fixed RNA–chromatin contact"


def workflow():
    reagent = BRIDGE.unligated_contact_scene(rna_length=28, dna_length=34)
    joined = BRIDGE.contact_scene()
    captured = Scene.duplex([seg("RNA-derived cDNA", "R" * 20, placeholder=True),
                             seg("biotin linker", "X" * 36, "me", placeholder=True),
                             seg("genomic DNA", "D" * 20, placeholder=True)],
                            label="captured, second-strand-completed contact")
    captured.mark("top", "biotin linker", "streptavidin capture")
    tags = Scene.duplex([seg("RNA tag", "R" * 20, placeholder=True),
                         seg("RNA MmeI boundary", "X" * 4, "me", placeholder=True),
                         seg("linker", "X" * 28, placeholder=True),
                         seg("DNA MmeI boundary", "X" * 4, "me", placeholder=True),
                         seg("DNA tag", "D" * 20, placeholder=True)], label="85-bp GRID tag pair")
    tags.mark("top", "RNA MmeI boundary", "MmeI cut")
    tags.mark("top", "DNA MmeI boundary", "MmeI cut")
    final = Scene.duplex(list(FINAL_LIBRARY), label="194-bp PCR library")
    w = Workflow(MolecularState(INITIAL_NAME, INITIAL_ROWS))
    w.react("Add the pre-annealed bivalent linker", reagent.rows(), name="Fixed contact plus free bridge",
            note="RNA, genomic DNA and both bridge strands are all retained; the bridge carries rNNN, internal biotin and a pre-adenylated RNA end.")
    w.react("Ligate RNA, reverse-transcribe, then ligate nearby DNA", joined.rows(), name="RNA–linker–DNA contact",
            note="T4 Rnl2(tr) joins RNA; SuperScript III copies it; T4 DNA ligase joins the AluI end.")
    w.react("Capture and make the second strand", captured.rows(), name="Captured duplex",
            note="Biotin capture, alkaline strand release and random-hexamer Klenow make dsDNA.")
    w.react("Cut both tags with MmeI", tags.rows(), name="MmeI tag pair",
            note="MmeI releases about 20 nt on each side; the desired band is 85 bp.")
    w.react("Ligate Y adapters and PCR", final.rows(), name="Sequencing library",
            note="Printed adapter and PCR-primer sequences produce the 194-bp library.")
    return w
