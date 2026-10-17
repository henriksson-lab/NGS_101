"""Li et al. GRID-seq bivalent-linker library."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Scene, Workflow, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from rna_dna_bridge import BivalentRNADNABridge

TITLE = "GRID-seq — RNA–chromatin contacts"
NOTES = "01_grid-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nbt.3968">Li et al. (2017)</a>.'
SUMMARY = "A biotinylated bivalent RNA–DNA linker captures a chromatin RNA and a nearby AluI-cut genomic end. MmeI releases a short RNA tag–linker–DNA tag molecule for Illumina sequencing."
CAVEAT = "The tag bases are molecule-specific placeholders. The 2017 paper, Methods 'Construction of GRID-seq library', prints linker L1 as /5Phos/rGrUrUrGrGrArUrUrCrNrNrNrGrACACAGC/iBiodT/CACTCCCACACACCGAACTCCAAC; the 2019 same-group Nature Protocols Table 1 prints /5Phos/rGrUrUrGrGrArUrUrCrNrNrNrGrACACAGC/iBiodT/CACTCCACACACCGAACTCCAAC, missing one C after CACT. The model uses the 2017 strand because it pairs fully with the L2 strand printed in both sources. Whether the shorter 2019 sequence was synthesized remains uncertain. The printed rG rA after rNNN are the first two bases of the L2-paired arm and are drawn as ribonucleotides within the duplex. The five-base PCR barcode is printed in both inspected methods, while a separate index-read priming recipe remains unverified."

GRID_R1 = il.TRUSEQ_READ1
GRID_P7_ARM = "AGATCGGAAGAGCACACGTCT"  # printed phosphorylated Y-adapter strand
GRID_Y_BOTTOM = "ACACTCTTTCCCTACACGACGCTCTTCCGATCTNN"
# 2017 Methods "Construction of GRID-seq library", printed PCR primers (B = barcode).
GRID_PCR1 = "AATGATACGGCGACCACCGAGATCTACACBBBBBACACTCTTTCCCTACACGACGCTCTTCCGATCT"
GRID_PCR2 = "CAAGCAGAAGACGGCATACGAGACGTGTGCTCTTCCGATCT"
GRID_BARCODE_LENGTH = GRID_PCR1.count("B")
if GRID_PCR1 != il.P5 + "B" * GRID_BARCODE_LENGTH + GRID_R1:
    raise ValueError("GRID PCR Primer #1 is not P5 + barcode + the printed Read 1 site")
if not GRID_PCR2.endswith(revcomp(GRID_P7_ARM)):
    raise ValueError("GRID PCR Primer #2 must end in the complement of the Y-adapter strand")
# Primer #2's 5' tail beyond the adapter complement; it is not the canonical
# P7 adapter string, so it is derived from the printed primer rather than il.P7_RC.
GRID_P7_TAIL = revcomp(GRID_PCR2[:-len(GRID_P7_ARM)])
# Printed L1: rG rU rU rG rG rA rU rU rC rN rN rN rG rA CACAGC /iBiodT/ CAC...
# The printed rG rA after rNNN are the first two bases of the L2-paired arm
# (revcomp(L2) begins GACACAGC), so no separate unpaired terminal base exists;
# the 2019 protocol likewise calls the linker 45 bp (9 + 3 + 33). Those two paired
# bases are declared as duplex ribonucleotide positions 0 and 1.
BRIDGE = BivalentRNADNABridge(
    "GRID-seq bridge", "GTTGGAGTTCGGTGTGTGGGAGTGAGCTGTGTC",
    "GUUGGAUUC", 3, "", "GACACAGCTCACTCCCACACACCGAACTCCAAC",
    biotin_offset=8, preadenylated_rna_5=True, barcode_id="grid_linker_umi",
    duplex_rna_positions=(0, 1))
LINKER_DNA = BRIDGE.dna_strand
LINKER_LENGTH = (len(BRIDGE.rna_fixed_5) + BRIDGE.barcode_length + len(BRIDGE.rna_terminal)
                 + len(BRIDGE.duplex_arm))
TAG_LENGTH = 20  # MmeI cuts ~20 nt beyond each linker-end site; 2019 protocol: 19-23 bp

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("five-base library index", "I" * GRID_BARCODE_LENGTH, "cbc", placeholder=True,
        feature=feature("sample_index_i5", "sample_index", "unknown")),
    seg("Read 1 site", GRID_R1, "r1"),
    seg("RNA-derived tag", "R" * TAG_LENGTH, placeholder=True, length_bp=TAG_LENGTH),
    seg("linker RNA-arm fixed bases", BRIDGE.rna_fixed_5.replace("U", "T")),
    seg("three-base molecular barcode", "N" * BRIDGE.barcode_length, "umi", placeholder=True,
        feature=feature("grid_linker_umi", "umi", "random")),
    seg("linker duplex fixed bases (first two printed as rG rA)", BRIDGE.duplex_arm),
    seg("genomic-DNA tag", "D" * TAG_LENGTH, placeholder=True, length_bp=TAG_LENGTH),
    seg("P7-side Y-adapter arm", GRID_P7_ARM, "r2"),
    seg("P7 tail from PCR Primer #2", GRID_P7_TAIL, "p7"),
], name="GRID-seq single-end library")
SEQ_PRIMERS = (sp.TRUSEQ["R1"],)
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=("Read 1",)):
    raise ValueError("GRID-seq final library: " + "; ".join(problems))
READ_LENGTHS = {"Read 1": 100}
SHOW_READ_CONTENT = True
FINAL_CAPTION = "The 100-nt read crosses an ~20-nt RNA tag, the orientation-defining bivalent linker and an ~20-nt genomic-DNA tag. The Y adapter is ligated to both MmeI ends, so either tag can lie next to Read 1; one orientation is drawn. PCR Primer #1 places a five-base multiplex barcode between P5 and the Read 1 site; the P7 end is the printed Primer #2 tail, not the canonical TruSeq P7 string. With 20-nt tags the modelled molecule is " + str(len(FINAL_LIBRARY)) + " bp; the 2019 protocol gives 20–21-nt tags, so the published 194-bp PCR band lies within the modelled range."
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
    tags = Scene.duplex([seg("RNA tag", "R" * TAG_LENGTH, placeholder=True),
                         seg("bivalent linker, MmeI site at each end", "X" * LINKER_LENGTH, "me", placeholder=True),
                         seg("DNA tag", "D" * TAG_LENGTH, placeholder=True)], label="~85-bp GRID tag pair")
    tags.mark("top", "RNA tag", "MmeI-cut end")
    tags.mark("top", "DNA tag", "MmeI-cut end")
    final = Scene.duplex(list(FINAL_LIBRARY), label=f"{len(FINAL_LIBRARY)}-bp PCR library (20-nt tags)")
    y=Scene();y.strand("Y bottom",[seg("Y bottom duplex arm",GRID_Y_BOTTOM[:-2],"r1"),seg("two-base Y overhang",GRID_Y_BOTTOM[-2:])],label="printed Y-adapter bottom")
    y.strand("Y phosphorylated strand",[seg("phosphorylated Y strand",GRID_P7_ARM,"r2")],label="5′-phosphorylated Y strand",mod5="p")
    y.overhang("Y bottom","two-base Y overhang",polarity="3-prime",owner="adapter",text="3′-NN MmeI-compatible overhang")
    y.strand("MmeI tag",[seg("85-bp tag pair","X"*30,placeholder=True),seg("two-base MmeI end","N"*2,placeholder=True)],label="MmeI-cut tag pair")
    y.overhang("MmeI tag","two-base MmeI end",polarity="3-prime",owner="substrate",text="two-base MmeI 3′ overhang")
    y_geometry=AdapterLigationGeometry.known("cohesive",y,statement="The 2017 Methods print a Y adapter whose unphosphorylated strand ends 3′-NN; the 2019 protocol states that MmeI leaves a two-base 3′ overhang and that the random NN adapter end preserves those bases. Shrimp alkaline phosphatase first removes tag 5′ phosphates to prevent self-ligation.")
    w = Workflow(MolecularState(INITIAL_NAME, INITIAL_ROWS))
    w.react("Cut chromatin with AluI, then treat with T4 PNK", INITIAL_ROWS, name="AluI-cut, PNK-treated contact",
            note="In fixed, SDS-permeabilized nuclei AluI leaves blunt genomic ends. The 2019 protocol states that T4 PNK phosphorylates 5′ ends and removes 3′ phosphates of DNA and RNA.")
    w.react("Add the pre-annealed bivalent linker", reagent.rows(), name="Fixed contact plus free bridge",
            note="RNA, genomic DNA and both bridge strands are all retained; the bridge carries rNNN, internal biotin and a pre-adenylated RNA end.")
    w.ligate_adapter("Ligate chromatin RNA to the bridge",rna_geometry,joined.rows(),name="RNA-ligated bridge",note="ATP-free T4 RNA ligase 2 truncated KQ joins the pre-adenylated bridge to the RNA 3′ end.")
    w.react("Reverse-transcribe the bridge-bound RNA",joined.rows(),name="RNA-copied bridge",note="No separate primer is added: SuperScript III extends the bridge DNA strand's 3′ end across rNNN and into the ligated RNA.")
    w.ligate_adapter("Ligate the bridge to nearby DNA",dna_geometry,joined.rows(),name="RNA–linker–DNA contact",note="After washing away free linker, T4 DNA ligase joins the blunt duplex bridge end to a nearby AluI end (16 °C overnight).")
    w.react("Capture and make the second strand", captured.rows(), name="Captured duplex",
            note="After proteinase K and DNA extraction, streptavidin captures the biotinylated strand; 150 mM NaOH releases the non-biotinylated cDNA-containing strand, which random hexamers and Klenow (3′→5′ exo−) make double-stranded.")
    w.react("Cut both tags with MmeI", tags.rows(), name="MmeI tag pair",
            note="MmeI cuts about 20 nt beyond each outward-facing linker site. A 12% native gel selects the 85-bp RNA–linker–DNA band; a 65-bp band diagnoses linkers joined to only one partner.")
    w.ligate_adapter("Ligate Y adapters and PCR",y_geometry,final.rows(),name="Sequencing library",
            note="After dephosphorylation and Y-adapter ligation, T4 PNK and ligase seal the remaining nick; a 10% gel band at 165–185 bp is amplified for 16 Phusion cycles and the 194-bp product is gel-recovered.")
    return w




def _y_adapter_scene():
    y = Scene()
    y.strand("Y bottom", [seg("Y bottom duplex arm", GRID_Y_BOTTOM[:-2], "r1"),
                          seg("two-base Y overhang", GRID_Y_BOTTOM[-2:])])
    y.strand("Y phosphorylated strand", [seg("phosphorylated Y strand", GRID_P7_ARM, "r2")],
             mod5="p")
    return y


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    contact = graph.input("rna_chromatin_contact", "fixed, permeabilized RNA–chromatin contact",
                          MaterialKind.BIOLOGICAL, state=INITIAL_ROWS, primary=True)

    # 2017 Methods: AluI in nuclei, then T4 PNK, both before linker ligation to RNA.
    alui = graph.input("alui", "AluI restriction enzyme", MaterialKind.CATALYST)
    cut_contact = graph.product("alui_cut_contact", "contact with blunt AluI genomic ends",
                                MaterialKind.MOLECULAR)
    graph.react("alui_digest", "Fragment chromatin in situ with AluI",
                StepKind.DIGESTION, ((contact, InputRole.TEMPLATE),
                                     (alui, InputRole.CATALYST)), (cut_contact,))
    pnk = graph.input("t4_pnk", "T4 polynucleotide kinase with ATP", MaterialKind.CATALYST)
    pnk_contact = graph.product("pnk_contact", "end-treated RNA and DNA in nuclei",
                                MaterialKind.MOLECULAR)
    graph.react("pnk_end_treatment", "Treat RNA and DNA ends with T4 PNK",
                StepKind.PHOSPHORYLATION, ((cut_contact, InputRole.TEMPLATE),
                                           (pnk, InputRole.CATALYST)), (pnk_contact,),
                reason="2019 protocol: PNK phosphorylates 5′ ends and dephosphorylates 3′ ends of DNA and RNA.")

    bridge = graph.input("bivalent_bridge", "pre-annealed biotinylated RNA–DNA bridge",
                         MaterialKind.ADAPTER, state=BRIDGE.scene())
    rna_ligase = graph.input("rna_ligase", "T4 RNA ligase 2 truncated KQ (ATP-free)",
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
                ((pnk_contact, InputRole.TEMPLATE), (bridge, InputRole.STRUCTURAL),
                 (rna_ligase, InputRole.CATALYST)), (rna_joined,), geometry=rna_geometry)

    # No separate RT primer: the annealed bridge DNA strand is extended into the RNA.
    rt = graph.input("reverse_transcriptase", "SuperScript III reverse transcriptase",
                     MaterialKind.CATALYST)
    copied = graph.product("rna_copied", "RNA-copied bridge contact",
                           MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe", "Extend the bridge DNA strand into the ligated RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna_joined, InputRole.TEMPLATE), (rna_joined, InputRole.SELF_PRIMER),
                 (rt, InputRole.CATALYST)), (copied,))

    wash = graph.product("free_linker", "unligated free linker washed from nuclei",
                         MaterialKind.WASTE)
    washed = graph.product("washed_nuclei", "washed nuclei with RNA-ligated bridges",
                           MaterialKind.MOLECULAR)
    graph.react("wash_free_linker", "Wash nuclei to remove free linker", StepKind.CLEANUP,
                ((copied, InputRole.TEMPLATE),), (washed, wash), changes_molecule=False)
    graph.discard(wash, "free linker is removed by washing nuclei in DNA ligase buffer")

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
    graph.react("dna_ligation", "Join the RNA-ligated bridge to a nearby AluI-cut DNA end in situ",
                StepKind.PROXIMITY_LIGATION,
                ((washed, InputRole.TEMPLATE), (dna_ligase, InputRole.CATALYST)),
                (chimera,), geometry=dna_geometry)

    proteinase = graph.input("proteinase_k", "proteinase K/SDS at 65 °C", MaterialKind.CATALYST)
    dna = graph.product("extracted_dna", "extracted total DNA with contact chimeras",
                        MaterialKind.MOLECULAR)
    graph.react("reverse_crosslink", "Digest protein, reverse crosslinks and extract DNA",
                StepKind.DEPROTEINIZATION,
                ((chimera, InputRole.TEMPLATE), (proteinase, InputRole.CATALYST)), (dna,),
                changes_molecule=False)

    beads = graph.input("streptavidin_beads", "streptavidin magnetic beads",
                        MaterialKind.SELECTION)
    captured = graph.product("captured_chimera", "bead-bound biotinylated contact chimera",
                             MaterialKind.MOLECULAR)
    uncaptured = graph.product("uncaptured", "DNA without a biotinylated bridge",
                               MaterialKind.WASTE)
    graph.react("capture", "Capture the biotinylated bridge on streptavidin",
                StepKind.AFFINITY_CAPTURE,
                ((dna, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, uncaptured))
    graph.discard(uncaptured, "removed by five B&W washes")

    naoh = graph.input("naoh", "150 mM NaOH", MaterialKind.CONDITION)
    released = graph.product("released_ssdna", "released cDNA-containing single strand",
                             MaterialKind.MOLECULAR)
    bead_strand = graph.product("bead_strand", "bead-retained biotinylated hybrid strand",
                                MaterialKind.WASTE)
    graph.react("alkaline_release", "Release the non-biotinylated strand with NaOH",
                StepKind.STRAND_SEPARATION,
                ((captured, InputRole.TEMPLATE), (naoh, InputRole.CONDITION)),
                (released, bead_strand))
    graph.discard(bead_strand, "stays on the beads; only the cleared supernatant is kept")

    random_primer = graph.input("random_primer", "random hexamer primers",
                                MaterialKind.OLIGO,
                                state=Construct([seg("random hexamer", "N" * 6)],
                                                name="random hexamer primer"))
    klenow = graph.input("klenow", "Klenow fragment (3′→5′ exo−)", MaterialKind.CATALYST)
    duplex = graph.product("contact_duplex", "random-primed double-stranded contact DNA",
                           MaterialKind.MOLECULAR)
    graph.react("second_strand", "Make the second DNA strand by random priming",
                StepKind.EXTENSION,
                ((released, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (klenow, InputRole.CATALYST)), (duplex,))

    mmei = graph.input("mmei", "MmeI restriction enzyme with SAM", MaterialKind.CATALYST)
    digested = graph.product("mmei_digest_mix", "MmeI-digested contact DNA",
                             MaterialKind.MOLECULAR)
    graph.react("mmei_digest", "Cut ~20 nt beyond both outward-facing linker MmeI sites",
                StepKind.DIGESTION,
                ((duplex, InputRole.TEMPLATE), (mmei, InputRole.CATALYST)), (digested,))

    gel85 = graph.input("gel_12pct", "12% native polyacrylamide gel", MaterialKind.SELECTION)
    tag_pair = graph.product("tag_pair", "85-bp RNA tag–linker–DNA tag pair",
                             MaterialKind.MOLECULAR)
    one_sided = graph.product("one_sided_65bp", "65-bp product with linker joined to only RNA or DNA",
                              MaterialKind.WASTE)
    other_dna = graph.product("mmei_waste", "flanking contact DNA and background smear",
                              MaterialKind.WASTE)
    graph.react("select_85bp", "Excise the 85-bp band", StepKind.SELECTION,
                ((digested, InputRole.TEMPLATE), (gel85, InputRole.SELECTION)),
                (tag_pair, one_sided, other_dna), changes_molecule=False)
    graph.discard(one_sided, "diagnostic 65-bp band left in the gel")
    graph.discard(other_dna, "outside the excised 85-bp band")

    sap = graph.input("sap", "shrimp alkaline phosphatase", MaterialKind.CATALYST)
    dephos = graph.product("dephosphorylated_tags", "5′-dephosphorylated tag pair",
                           MaterialKind.MOLECULAR)
    graph.react("dephosphorylate", "Dephosphorylate tag ends to prevent self-ligation",
                StepKind.DEPHOSPHORYLATION,
                ((tag_pair, InputRole.TEMPLATE), (sap, InputRole.CATALYST)), (dephos,))

    y_adapter = graph.input("y_adapter", "two-base 3′-overhang Y adapter", MaterialKind.ADAPTER,
                            state=_y_adapter_scene())
    y_ligated = graph.product("y_ligated", "Y-adapter-ligated tag pair with nicks",
                              MaterialKind.MOLECULAR)
    y = _y_adapter_scene()
    y.overhang("Y bottom", "two-base Y overhang", polarity="3-prime", owner="adapter")
    y.strand("MmeI tag", [seg("85-bp tag pair", "X" * 30, placeholder=True),
                          seg("two-base MmeI end", "N" * 2, placeholder=True)])
    y.overhang("MmeI tag", "two-base MmeI end", polarity="3-prime", owner="substrate")
    y_geometry = AdapterLigationGeometry.known(
        "cohesive", y,
        statement="The printed Y-adapter 3′-NN end accepts the two-base MmeI 3′ overhang (2019 protocol).")
    graph.react("y_ligation", "Ligate the tag pair to the Y adapter",
                StepKind.LIGATION,
                ((dephos, InputRole.TEMPLATE), (y_adapter, InputRole.STRUCTURAL),
                 (dna_ligase, InputRole.CATALYST)), (y_ligated,), geometry=y_geometry)
    sealed = graph.product("y_sealed", "nick-sealed adapter-ligated tag pair",
                           MaterialKind.MOLECULAR)
    graph.react("seal_nicks", "Phosphorylate and seal the unligated nick",
                StepKind.NICK_SEALING,
                ((y_ligated, InputRole.TEMPLATE), (pnk, InputRole.CATALYST),
                 (dna_ligase, InputRole.CATALYST)), (sealed,))
    gel10 = graph.input("gel_10pct", "10% native polyacrylamide gel", MaterialKind.SELECTION)
    selected = graph.product("y_selected", "165–185-bp adapter-ligated band",
                             MaterialKind.MOLECULAR)
    free_adapter = graph.product("excess_adapter", "excess adapters and other bands",
                                 MaterialKind.WASTE)
    graph.react("select_ligated", "Gel-select the adapter-ligated band", StepKind.SELECTION,
                ((sealed, InputRole.TEMPLATE), (gel10, InputRole.SELECTION)),
                (selected, free_adapter), changes_molecule=False)
    graph.discard(free_adapter, "outside the excised 165–185-bp band")

    pcr1 = graph.input("pcr_primer_1", "GRID PCR Primer #1 with five-base index",
                       MaterialKind.OLIGO, state=Construct([
                           seg("P5", il.P5, "p5"),
                           seg("five-base library index", "B" * GRID_BARCODE_LENGTH, "cbc",
                               placeholder=True,
                               feature=feature("sample_index_i5", "sample_index", "unknown")),
                           seg("Read 1 site", GRID_R1, "r1")], name="GRID PCR Primer #1"))
    pcr2 = graph.input("pcr_primer_2", "GRID PCR Primer #2", MaterialKind.OLIGO,
                       state=Construct([seg("PCR Primer #2", GRID_PCR2, "p7")],
                                       name="GRID PCR Primer #2"))
    polymerase = graph.input("pcr_polymerase", "Phusion High-Fidelity DNA polymerase",
                             MaterialKind.CATALYST)
    amplified = graph.product("pcr_product", "16-cycle PCR product", MaterialKind.MOLECULAR)
    graph.react("library_pcr", "PCR-amplify the GRID-seq library", StepKind.PCR,
                ((selected, InputRole.TEMPLATE), (pcr1, InputRole.STRUCTURAL),
                 (pcr2, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="The paper prints both primers, but the primer-binding geometry on the Y-ligated template is not drawn in this model."))
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY)
    pcr_waste = graph.product("pcr_gel_waste", "primers and off-size PCR products",
                              MaterialKind.WASTE)
    graph.react("select_194bp", "Gel-recover the 194-bp library", StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (gel10, InputRole.SELECTION)),
                (library, pcr_waste), changes_molecule=False)
    graph.discard(pcr_waste, "outside the recovered 194-bp band")
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY,
                                           primers=tuple(SEQ_PRIMERS),
                                           read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope(tuple(graph.materials) + tuple(graph.steps), SOURCE)
    return graph.require_publishable()
