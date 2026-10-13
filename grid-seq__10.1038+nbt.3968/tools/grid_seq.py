"""Li et al. GRID-seq bivalent-linker library."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Scene, Workflow, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from rna_dna_bridge import BivalentRNADNABridge

TITLE = "GRID-seq — RNA–chromatin contacts"
NOTES = "01_grid-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nbt.3968">Li et al. (2017)</a>.'
SUMMARY = "A biotinylated bivalent RNA–DNA linker captures a chromatin RNA and a nearby AluI-cut genomic end. MmeI releases a short RNA tag–linker–DNA tag molecule for Illumina sequencing."
CAVEAT = "The tag bases are molecule-specific placeholders; linker bases, modifications and orientation are transcribed from the paper. The five-base PCR barcode is published, but its index-read priming recipe is not."

GRID_R1 = il.TRUSEQ_READ1
GRID_P7_ARM = "AGATCGGAAGAGCACACGTCT"  # printed phosphorylated Y-adapter strand
GRID_Y_BOTTOM = "ACACTCTTTCCCTACACGACGCTCTTCCGATCTNN"
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
    seg("RNA-side MmeI boundary", "X" * 6, "me", placeholder=True, length_bp=6),
    seg("linker RNA-arm fixed bases", BRIDGE.rna_fixed_5.replace("U", "T")),
    seg("three-base molecular barcode", "N" * BRIDGE.barcode_length, "umi", placeholder=True,
        feature=feature("grid_linker_umi", "umi", "random")),
    seg("linker RNA-arm terminal G", BRIDGE.rna_terminal),
    seg("linker duplex fixed bases", BRIDGE.duplex_arm),
    seg("DNA-side MmeI boundary", "X" * 6, "me", placeholder=True, length_bp=6),
    seg("genomic-DNA tag", "D" * 20, placeholder=True, length_bp=20),
    seg("P7-side Y-adapter arm", GRID_P7_ARM, "r2"),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="GRID-seq single-end library")
SEQ_PRIMERS = (sp.TRUSEQ["R1"],)
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=("Read 1",)):
    raise ValueError("GRID-seq final library: " + "; ".join(problems))
READ_LENGTHS = {"Read 1": 100}
SHOW_READ_CONTENT = True
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
    reagent.ligation_end("bridge hybrid", "RNA-arm fixed bases", kind="single-strand", polarity="5-prime", owner="adapter", text="5′-App RNA-ligation end")
    reagent.ligation_end("chromatin RNA", "chromatin-associated RNA", kind="single-strand", polarity="3-prime", owner="substrate", text="RNA 3′-OH")
    rna_geometry=AdapterLigationGeometry.known("single-strand",reagent,statement="The paper prints the mixed RNA/DNA bridge strand and states that its RNA-facing 5′ end is pre-adenylated for ATP-free ligation to chromatin RNA 3′ ends.")
    joined = BRIDGE.contact_scene()
    dna_ready=BRIDGE.unligated_contact_scene(rna_length=28,dna_length=34)
    dna_ready.ligation_end("bridge DNA","bridge DNA strand",kind="blunt",owner="adapter",text="5′-phosphorylated blunt bridge end")
    dna_ready.ligation_end("genomic DNA","nearby genomic DNA",kind="blunt",owner="substrate",text="blunt AluI genomic end")
    dna_geometry=AdapterLigationGeometry.known("blunt",dna_ready,statement="AluI produces a blunt genomic end, and the source prints the 5′-phosphorylated DNA bridge strand used for T4 DNA ligation.")
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
    y=Scene();y.strand("Y bottom",[seg("Y bottom duplex arm",GRID_Y_BOTTOM[:-2],"r1"),seg("two-base Y overhang",GRID_Y_BOTTOM[-2:],placeholder=True)],label="printed Y-adapter bottom")
    y.strand("Y phosphorylated strand",[seg("phosphorylated Y strand",GRID_P7_ARM,"r2")],label="5′-phosphorylated Y strand",mod5="p")
    y.overhang("Y bottom","two-base Y overhang",polarity="3-prime",owner="adapter",text="3′-NN MmeI-compatible overhang")
    y.strand("MmeI tag",[seg("two-base MmeI end","N"*2,placeholder=True),seg("85-bp tag pair","X"*30,placeholder=True)],label="MmeI-cut tag pair")
    y.overhang("MmeI tag","two-base MmeI end",polarity="5-prime",owner="substrate",text="two-base MmeI end")
    y_geometry=AdapterLigationGeometry.known("cohesive",y,statement="The GRID-seq Methods print a Y adapter whose bottom strand terminates in NN; MmeI creates the matching two-base tag end before ligation.")
    w = Workflow(MolecularState(INITIAL_NAME, INITIAL_ROWS))
    w.react("Add the pre-annealed bivalent linker", reagent.rows(), name="Fixed contact plus free bridge",
            note="RNA, genomic DNA and both bridge strands are all retained; the bridge carries rNNN, internal biotin and a pre-adenylated RNA end.")
    w.ligate_adapter("Ligate chromatin RNA to the bridge",rna_geometry,joined.rows(),name="RNA-ligated bridge",note="T4 Rnl2(tr) joins the pre-adenylated bridge to the RNA 3′ end.")
    w.react("Reverse-transcribe the bridge-bound RNA",joined.rows(),name="RNA-copied bridge",note="SuperScript III copies the ligated RNA.")
    w.ligate_adapter("Ligate the bridge to nearby DNA",dna_geometry,joined.rows(),name="RNA–linker–DNA contact",note="T4 DNA ligase joins the blunt AluI end.")
    w.react("Capture and make the second strand", captured.rows(), name="Captured duplex",
            note="Biotin capture, alkaline strand release and random-hexamer Klenow make dsDNA.")
    w.react("Cut both tags with MmeI", tags.rows(), name="MmeI tag pair",
            note="MmeI releases about 20 nt on each side; the desired band is 85 bp.")
    w.ligate_adapter("Ligate Y adapters and PCR",y_geometry,final.rows(),name="Sequencing library",
            note="Printed adapter and PCR-primer sequences produce the 194-bp library.")
    return w


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    contact = graph.input("rna_chromatin_contact", "fixed RNA–chromatin contact",
                          MaterialKind.BIOLOGICAL, state=INITIAL_ROWS, primary=True)
    bridge = graph.input("bivalent_bridge", "pre-annealed biotinylated RNA–DNA bridge",
                         MaterialKind.ADAPTER)
    rna_ligase = graph.input("rna_ligase", "T4 RNA ligase 2 truncated",
                             MaterialKind.CATALYST)
    rna_joined = graph.product("rna_joined", "RNA-ligated bivalent bridge",
                               MaterialKind.MOLECULAR)
    reagent = BRIDGE.unligated_contact_scene(rna_length=28, dna_length=34)
    reagent.ligation_end("bridge hybrid", "RNA-arm fixed bases", kind="single-strand",
                         polarity="5-prime", owner="adapter", text="5′-App bridge end")
    reagent.ligation_end("chromatin RNA", "chromatin-associated RNA", kind="single-strand",
                         polarity="3-prime", owner="substrate", text="RNA 3′-OH")
    rna_geometry = AdapterLigationGeometry.known(
        "single-strand", reagent,
        statement="The printed bridge has a pre-adenylated RNA-facing 5′ end.")
    graph.react("rna_ligation", "Ligate chromatin RNA to the bivalent bridge",
                StepKind.LIGATION,
                ((contact, InputRole.TEMPLATE), (bridge, InputRole.STRUCTURAL),
                 (rna_ligase, InputRole.CATALYST)), (rna_joined,), geometry=rna_geometry)

    rt_primer = graph.input("rt_primer", "bridge-bound reverse-transcription primer",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "SuperScript III reverse transcriptase",
                     MaterialKind.CATALYST)
    copied = graph.product("rna_copied", "RNA-copied bridge contact",
                           MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe", "Copy bridge-bound chromatin RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna_joined, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (copied,))

    alui = graph.input("alui", "AluI restriction enzyme", MaterialKind.CATALYST)
    cut_contact = graph.product("alui_cut_contact", "contact with a nearby blunt AluI end",
                                MaterialKind.MOLECULAR)
    digest_waste = graph.product("alui_waste", "distal AluI restriction fragments",
                                 MaterialKind.WASTE)
    graph.react("alui_digest", "Expose a nearby blunt genomic end with AluI",
                StepKind.DIGESTION, ((copied, InputRole.TEMPLATE),
                                     (alui, InputRole.CATALYST)),
                (cut_contact, digest_waste))
    graph.discard(digest_waste, "distal restriction fragments are not bridge-linked")

    dna_ligase = graph.input("dna_ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    chimera = graph.product("contact_chimera", "RNA–bridge–DNA contact chimera",
                            MaterialKind.MOLECULAR,
                            state=tuple(BRIDGE.contact_scene().rows()))
    ready = BRIDGE.unligated_contact_scene(rna_length=28, dna_length=34)
    ready.ligation_end("bridge DNA", "bridge DNA strand", kind="blunt",
                       owner="adapter", text="5′-phosphorylated blunt bridge end")
    ready.ligation_end("genomic DNA", "nearby genomic DNA", kind="blunt",
                       owner="substrate", text="blunt AluI genomic end")
    dna_geometry = AdapterLigationGeometry.known(
        "blunt", ready,
        statement="AluI leaves a blunt DNA end for ligation to the printed phosphorylated bridge strand.")
    graph.react("dna_ligation", "Join the bridge to the nearby AluI-cut DNA end",
                StepKind.LIGATION,
                ((cut_contact, InputRole.TEMPLATE), (bridge, InputRole.STRUCTURAL),
                 (dna_ligase, InputRole.CATALYST)), (chimera,), geometry=dna_geometry)

    beads = graph.input("streptavidin_beads", "streptavidin capture beads",
                        MaterialKind.SELECTION)
    captured = graph.product("captured_chimera", "captured biotinylated contact chimera",
                             MaterialKind.MOLECULAR)
    uncaptured = graph.product("uncaptured", "uncaptured contact material", MaterialKind.WASTE)
    graph.react("capture", "Capture the biotinylated bridge on streptavidin",
                StepKind.AFFINITY_CAPTURE,
                ((chimera, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, uncaptured))
    graph.discard(uncaptured, "material without a recoverable biotinylated bridge is washed away")

    random_primer = graph.input("random_primer", "random hexamer second-strand primer",
                                MaterialKind.OLIGO)
    klenow = graph.input("klenow", "Klenow fragment", MaterialKind.CATALYST)
    duplex = graph.product("contact_duplex", "second-strand-completed contact duplex",
                           MaterialKind.MOLECULAR)
    graph.react("second_strand", "Make the second DNA strand",
                StepKind.EXTENSION,
                ((captured, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (klenow, InputRole.CATALYST)), (duplex,))

    mmei = graph.input("mmei", "MmeI restriction enzyme", MaterialKind.CATALYST)
    tag_pair = graph.product("tag_pair", "MmeI-released RNA–linker–DNA tag pair",
                             MaterialKind.MOLECULAR)
    mmei_waste = graph.product("mmei_waste", "flanking contact DNA outside the MmeI tags",
                               MaterialKind.WASTE)
    graph.react("mmei_digest", "Release the paired contact tags with MmeI",
                StepKind.DIGESTION,
                ((duplex, InputRole.TEMPLATE), (mmei, InputRole.CATALYST)),
                (tag_pair, mmei_waste))
    graph.discard(mmei_waste, "flanking sequence outside the released tag pair is discarded")

    y_adapter = graph.input("y_adapter", "two-base-overhang Y adapter", MaterialKind.ADAPTER)
    y_ligated = graph.product("y_ligated", "Y-adapter-ligated GRID tag pair",
                              MaterialKind.MOLECULAR)
    y = Scene()
    y.strand("Y bottom", [seg("Y bottom duplex arm", GRID_Y_BOTTOM[:-2], "r1"),
                          seg("two-base Y overhang", GRID_Y_BOTTOM[-2:], placeholder=True)])
    y.strand("Y phosphorylated strand", [seg("phosphorylated Y strand", GRID_P7_ARM, "r2")],
             mod5="p")
    y.overhang("Y bottom", "two-base Y overhang", polarity="3-prime", owner="adapter")
    y.strand("MmeI tag", [seg("two-base MmeI end", "N" * 2, placeholder=True),
                          seg("85-bp tag pair", "X" * 30, placeholder=True)])
    y.overhang("MmeI tag", "two-base MmeI end", polarity="5-prime", owner="substrate")
    y_geometry = AdapterLigationGeometry.known(
        "cohesive", y,
        statement="The printed Y-adapter NN overhang matches the two-base MmeI tag end.")
    graph.react("y_ligation", "Ligate the MmeI tag pair to the Y adapter",
                StepKind.LIGATION,
                ((tag_pair, InputRole.TEMPLATE), (y_adapter, InputRole.STRUCTURAL),
                 (dna_ligase, InputRole.CATALYST)), (y_ligated,), geometry=y_geometry)

    pcr1 = graph.input("pcr_primer_1", "GRID PCR Primer #1 with five-base index",
                       MaterialKind.OLIGO)
    pcr2 = graph.input("pcr_primer_2", "GRID PCR Primer #2", MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY)
    graph.react("library_pcr", "PCR-complete the GRID-seq library", StepKind.PCR,
                ((y_ligated, InputRole.TEMPLATE), (pcr1, InputRole.STRUCTURAL),
                 (pcr2, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="The paper does not state a complete two-primer binding geometry."))
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY,
                                           primers=tuple(SEQ_PRIMERS),
                                           read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope((
        "rna_chromatin_contact", "bivalent_bridge", "rna_ligase", "rna_ligation",
        "rt_primer", "reverse_transcriptase", "reverse_transcribe", "alui", "alui_digest",
        "dna_ligase", "dna_ligation", "streptavidin_beads", "capture", "random_primer",
        "klenow", "second_strand", "mmei", "mmei_digest", "y_adapter", "y_ligation",
        "pcr_primer_1", "pcr_primer_2", "pcr_polymerase", "library_pcr", "library"), SOURCE)
    return graph.require_publishable()
