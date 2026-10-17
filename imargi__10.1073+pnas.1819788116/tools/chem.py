"""iMARGI in-situ RNA–DNA bridge ligation."""
from multibranch import illumina_branch
from interaction_capture import BivalentRNAChromatinBridge
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState,Row,Scene,Segment,Workflow
from chemdraw import Construct
from endprep import dA_tailed_scene, repair_and_dA_tail
from pcr import PCRGeometry
from circular import CircularizationGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)
TITLE="iMARGI — in situ RNA–genome interaction mapping"
NOTES="01_imargi.html"
SOURCE='Defining source: <a href="https://doi.org/10.1073/pnas.1819788116">Yan et al., <i>PNAS</i> (2019)</a>; detailed protocol <a href="https://doi.org/10.1038/s41596-019-0229-4">Nat. Protocols (2019)</a>.'
SUMMARY="A bivalent linker is ligated in situ first to fragmented chromatin-associated RNA and then to nearby genomic DNA; reverse transcription and Illumina preparation preserve the RNA–linker–DNA order."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("iMARGI","RNA–linker–DNA contact chimera",inferred_adapters=True)
FINAL_CAPTION="INFERRED — paired reads recover the RNA-derived and DNA-derived ends; the exact PCR landing geometry is not represented by this compact library view. The linker has two ligation-bias bases, not a UMI."
SEQUENCING_INTRO="Declared Illumina primers bind the completed iMARGI library."
CAVEAT=(
    "The author protocol prints the modified bivalent linker, cut oligo and PCR primers. "
    "The workflow follows Yan et al. 2019 Methods and the Wu et al. 2019 Procedure (Steps 73–114): exonuclease biotin removal, streptavidin capture, on-bead RT, 5′ phosphorylation, NaOH release, CircLigase circularization, Cut_oligo/BamHI relinearization, PCR and gel size selection. "
    "Molecular states after cross-link reversal, and the strand-resolved BamHI split and PCR landing junctions, are not drawn; the compact final-library view is inferred. "
    "Published citation discrepancy: Yan et al. 2019 reference 10 cites "
    "<a href=\"https://doi.org/10.1016/j.cub.2017.01.068\">DOI 10.1016/j.cub.2017.01.068</a> "
    "for Sridhar et al. MARGI, but that DOI is the erratum; "
    "the <a href=\"https://doi.org/10.1016/j.cub.2017.01.011\">original article</a> "
    "has DOI 10.1016/j.cub.2017.01.011. The model uses the "
    "iMARGI detailed-protocol oligos. Original Sridhar 2017 MARGI Document S1 "
    "(Supplemental Experimental Procedures, lines 204–215) prints a top linker "
    "ending ...TCCGA/BiodT/CT without a C* mark; later Wu 2019 iMARGI Nature "
    "Protocols (JATS lines 275–285) prints ...TCC*GA/iBiodT/CT and defines C* "
    "as phosphorothioated C. The bases agree, but the printed modifications differ. "
    "This model uses the later iMARGI-specific C* oligo; the 2019 PNAS paper's "
    "statement that it reused the MARGI linker sequence does not establish "
    "which modification was in its experimental oligo."
    " Apparent reagent typo in Wu et al. 2019 detailed protocol: stock-preparation text "
    "(extracted JATS lines 395–397) says ATP for a reagent labeled 10 mM dATP, "
    "while its dA-tailing reaction table (lines 673–685) and Yan et al. 2019 "
    "Methods (line 203) specify dATP. This model uses a 3′-dA tail from the "
    "reaction descriptions; the reagent actually used in the experiment remains unverified."
    " Wu et al. 2019 also names the PCR master mix inconsistently: Step 104 (pilot PCR, JATS line 1159) says "
    "“2× NEBNext HC PCR Master Mix”, whereas Step 106 (line 1173 table) and the reagent list (line 259) give "
    "NEBNext High-Fidelity 2× PCR Master Mix (M0541). The model uses the High-Fidelity mix; “HC” appears to be a typo."
)
READ_LENGTHS={"Read 1":100,"Read 2":100}
def oligos():
 return [
  "<b>bivalent linker top</b> — 5′-/Phos/NNAGATCGGAAGAGCGTCGTGTAGGGAGGATCCGTTCAGACGTGTGCTCTTCC*GA/iBiodT/CT-3′ (C* is phosphorothioated)",
  "<b>bivalent linker bottom</b> — 5′-/Phos/GATCGGAAGAGCACACGTCTGAACGGATCCTCCCTACACGACGCTCT-3′",
  "<b>BamHI cut oligo</b> — 5′-TCGTGTAGGGAGGATCCGTTCAGACGTGTGCTCT/3InvdT/-3′",
  "<b>NEBNext universal PCR</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′",
  "<b>NEBNext index PCR</b> — 5′-CAAGCAGAAGACGGCATACGAGAT[6 bp barcode]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-S-T-3′",
 ]
BRIDGE=BivalentRNAChromatinBridge("iMARGI bivalent linker")
def rna_ligation_geometry():
    rna=Scene(); rna.strand("RNA arm",[Segment("preadenylated RNA-ligation end","[RNA-arm end]",placeholder=True),Segment("annealed bivalent-linker body","[bivalent linker]",placeholder=True)],label="iMARGI linker")
    rna.ligation_end("RNA arm","preadenylated RNA-ligation end",kind="single-strand",polarity="5-prime",owner="adapter",text="preadenylated 5′ RNA-ligation end")
    return AdapterLigationGeometry.known("single-strand",rna,statement="The detailed protocol identifies a preadenylated RNA-compatible end on the annealed bivalent linker; this is a single-strand end ligation, not a duplex overhang.")
def dna_ligation_geometry():
    linker=Scene()
    linker.strand("linker top",[
        Segment("5-prime RNA ligation tail", "NNAGATCGGA"),
        Segment("paired linker core", "AGAGCGTCGTGTAGGGAGGATCCGTTCAGACGTGTGCTCTTCCGATC"),
        Segment("3-prime T overhang", "T")],label="iMARGI bivalent linker")
    linker.anneal("linker bottom",[
        Segment("paired bottom core", "GATCGGAAGAGCACACGTCTGAACGGATCCTCCCTACACGACGCTCT")],
        to="linker top",pair=("paired bottom core","paired linker core"),
        label="5′-phosphorylated bottom strand")
    linker.ligation_end("linker top","3-prime T overhang",kind="cohesive",
                      polarity="3-prime",owner="adapter",text="3′-T linker overhang")
    substrate=dA_tailed_scene(repair_and_dA_tail(Construct([
        Segment("AluI-cut genomic fragment", "X"*24,placeholder=True)])),
        label="dA-tailed genomic DNA end")
    return AdapterLigationGeometry.known("cohesive",substrate,linker,
        statement="The detailed iMARGI protocol describes dA tailing and sticky-end DNA ligation (Wu et al. 2019, lines 101–103); the original MARGI Fig. 1 identifies a linker 3′-T overhang paired with the genomic 3′-A tail. The 2019 linker strands print the same paired core and terminal T (lines 273–285).")
def sections():
    wf=Workflow(MolecularState("Crosslinked and fragmented chromatin",tuple(BRIDGE.unjoined().rows())))
    wf.ligate_adapter("Ligate the RNA end of the bivalent linker in situ",rna_ligation_geometry(),BRIDGE.rna_joined().rows(),note="The linker end joins fragmented RNA; its DNA-reactive end remains free.")
    wf.ligate_adapter("Ligate the bridge DNA end to nearby genomic DNA",dna_ligation_geometry(),BRIDGE.chimera().rows(),note="The second ligation makes the diagnostic RNA–linker–DNA chimera in the nucleus.")
    wf.react("Reverse-transcribe, purify and make the paired-end library",[Row(chunks=[("RNA-derived cDNA—linker—genomic DNA → paired-end adapters","r1",False)])],note="The linker orientation distinguishes RNA and DNA read ends.")
    return [(r.action,r.after.rows,r.note) for r in wf.reactions]

def protocol_graph():
 source="Defining source: https://doi.org/10.1073/pnas.1819788116; detailed protocol: https://doi.org/10.1038/s41596-019-0229-4"
 graph=ProtocolGraph("iMARGI",authority=source)
 chromatin=graph.input("chromatin","formaldehyde-cross-linked, SDS-permeabilized nuclei",
                       MaterialKind.BIOLOGICAL,primary=True)
 nuclease=graph.input("restriction_enzyme","AluI (with RNasinPlus, CutSmart)",MaterialKind.CATALYST)
 cut=graph.product("alui_cut_chromatin","AluI-digested chromatin in nuclei",MaterialKind.MOLECULAR)
 graph.react("fragmentation","Digest chromatin with AluI overnight at 37 °C",StepKind.DIGESTION,
             ((chromatin,InputRole.TEMPLATE),(nuclease,InputRole.CATALYST)),(cut,))
 rnase=graph.input("rnase_i","RNase I (1:10, 3 min, 37 °C)",MaterialKind.CATALYST)
 fragmented=graph.product("fragmented_contacts","fragmented RNA and AluI-cut chromatin",
                          MaterialKind.MOLECULAR,state=BRIDGE.unjoined())
 graph.react("rna_fragmentation","Briefly fragment chromatin-associated RNA with RNase I",
             StepKind.FRAGMENTATION,((cut,InputRole.TEMPLATE),(rnase,InputRole.CATALYST)),(fragmented,))
 pnk=graph.input("t4_pnk_phosphatase","T4 PNK in pH 6.5 phosphatase buffer",MaterialKind.CATALYST)
 rna_oh=graph.product("rna_3oh","RNA fragments with 3′-OH ends",MaterialKind.MOLECULAR,
                      state=BRIDGE.unjoined())
 graph.react("rna_3_dephosphorylation","Remove RNA 3′ phosphates in situ",StepKind.DEPHOSPHORYLATION,
             ((fragmented,InputRole.TEMPLATE),(pnk,InputRole.CATALYST)),(rna_oh,))
 klenow=graph.input("klenow_exo_minus","Klenow fragment (3′→5′ exo−) with 0.1 mM dATP",MaterialKind.CATALYST)
 tailed=graph.product("da_tailed_contacts","RNA 3′-OH ends and dA-tailed genomic DNA ends",
                      MaterialKind.MOLECULAR,state=BRIDGE.unjoined())
 graph.react("dna_a_tailing","Add a 3′-dA tail to AluI-cut genomic DNA in situ",StepKind.TAILING,
             ((rna_oh,InputRole.TEMPLATE),(klenow,InputRole.CATALYST)),(tailed,))
 linker=graph.input("bivalent_linker","5′-adenylated, annealed iMARGI bivalent linker",MaterialKind.ADAPTER)
 rna_ligase=graph.input("rna_ligase","T4 RNA ligase 2, truncated KQ (20% PEG 8000)",MaterialKind.CATALYST)
 rna_joined=graph.product("rna_joined","RNA-end-ligated bivalent linker contacts",
                         MaterialKind.MOLECULAR,state=BRIDGE.rna_joined())
 graph.react("rna_ligation","Ligate the RNA end of the bivalent linker in situ",
             StepKind.LIGATION,
             ((tailed,InputRole.TEMPLATE),(linker,InputRole.STRUCTURAL),
              (rna_ligase,InputRole.CATALYST)),(rna_joined,),geometry=rna_ligation_geometry())
 washed=graph.product("washed_rna_joined","washed nuclei with RNA-ligated linker",
                      MaterialKind.MOLECULAR,state=BRIDGE.rna_joined())
 free_linker=graph.product("free_linker","excess unligated linker",MaterialKind.WASTE)
 graph.react("free_linker_wash","Wash nuclei five times to remove excess free linker",StepKind.CLEANUP,
             ((rna_joined,InputRole.TEMPLATE),),(washed,free_linker))
 graph.discard(free_linker,"removed in PNK-buffer washes of intact nuclei")
 dna_ligase=graph.input("dna_ligase","T4 DNA ligase (2 mL proximity-ligation volume)",MaterialKind.CATALYST)
 chimera=graph.product("contact_chimera","RNA–linker–DNA chimeras plus incomplete RNA–linker products",
                       MaterialKind.MOLECULAR,
                       state=BRIDGE.chimera())
 graph.react("dna_ligation","Ligate the bridge DNA end to nearby genomic DNA",
             StepKind.LIGATION,
             ((washed,InputRole.TEMPLATE),(linker,InputRole.STRUCTURAL),
              (dna_ligase,InputRole.CATALYST)),(chimera,),geometry=dna_ligation_geometry())
 protease=graph.input("proteinase_k","Proteinase K extraction buffer (1% SDS, 65 °C, 3 h)",MaterialKind.CATALYST)
 reversed_=graph.product("decrosslinked","cross-link-reversed nucleic acids",MaterialKind.MOLECULAR)
 graph.react("reverse_crosslinks","Lyse nuclei and reverse cross-links",StepKind.DEPROTEINIZATION,
             ((chimera,InputRole.TEMPLATE),(protease,InputRole.CATALYST)),(reversed_,))
 phenol=graph.input("phenol_chloroform","phenol:chloroform:isoamyl alcohol, ethanol precipitation",MaterialKind.CONDITION)
 nucleic=graph.product("purified_nucleic_acids","purified DNA and RNA",MaterialKind.MOLECULAR)
 graph.react("extract_nucleic_acids","Extract and precipitate DNA and RNA",StepKind.EXTRACTION,
             ((reversed_,InputRole.TEMPLATE),(phenol,InputRole.CONDITION)),(nucleic,))
 exo=graph.input("exo_i_t4_dna_polymerase","Exonuclease I, then T4 DNA polymerase with dATP/dGTP at 12 °C",MaterialKind.CATALYST)
 debiotinylated=graph.product("biotin_removed","chimeras keep biotin; biotin cleaved from linkers not ligated to DNA",
                              MaterialKind.MOLECULAR)
 graph.react("remove_unligated_biotin","Remove biotin from linkers not ligated to DNA",StepKind.DIGESTION,
             ((nucleic,InputRole.TEMPLATE),(exo,InputRole.CATALYST)),(debiotinylated,))
 streptavidin=graph.input("streptavidin_c1","Dynabeads MyOne Streptavidin C1, high-salt biotin washes",MaterialKind.SELECTION)
 captured=graph.product("captured_chimera","bead-bound biotinylated RNA–linker–DNA chimeras",MaterialKind.MOLECULAR)
 unbound=graph.product("unbound_nucleic_acids","de-biotinylated RNA–linker, unligated RNA/DNA and washed-off nucleic acids",MaterialKind.WASTE)
 graph.react("streptavidin_pulldown","Pull down biotinylated chimeras",StepKind.AFFINITY_CAPTURE,
             ((debiotinylated,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),(captured,unbound))
 graph.discard(unbound,"supernatant and high-salt washes of the streptavidin capture")
 rt=graph.input("reverse_transcriptase","SuperScript III, 50 °C, on beads",MaterialKind.CATALYST)
 cdna=graph.product("cdna_chimera","bead-bound duplex with cDNA extended from the linker bottom strand",MaterialKind.MOLECULAR)
 graph.react("reverse_transcription","Reverse-transcribe the RNA from the linker bottom strand",
             StepKind.REVERSE_TRANSCRIPTION,
             ((captured,InputRole.TEMPLATE),(linker,InputRole.STRUCTURAL),
              (rt,InputRole.CATALYST)),(cdna,))
 kinase=graph.input("t4_pnk_kinase","T4 PNK with 1 mM ATP (pH 7.6)",MaterialKind.CATALYST)
 phos=graph.product("phosphorylated_cdna","bead-bound duplex with 5′-phosphorylated DNA ends",MaterialKind.MOLECULAR)
 graph.react("phosphorylate_5_ends","Phosphorylate DNA 5′ ends on beads",StepKind.PHOSPHORYLATION,
             ((cdna,InputRole.TEMPLATE),(kinase,InputRole.CATALYST)),(phos,))
 naoh=graph.input("denaturing_buffer","0.1 M NaOH, 0.1 mM EDTA",MaterialKind.CONDITION)
 released=graph.product("released_ssdna","released ssDNA: genomic DNA–linker bottom strand–cDNA",MaterialKind.MOLECULAR)
 bead_strand=graph.product("bead_bound_strand","biotinylated top-strand/RNA side retained on beads",MaterialKind.WASTE)
 graph.react("release_ssdna","Release the cDNA-containing strand with NaOH",StepKind.STRAND_SEPARATION,
             ((phos,InputRole.TEMPLATE),(naoh,InputRole.CONDITION)),(released,bead_strand))
 graph.discard(bead_strand,"beads magnetically removed; supernatant kept")
 ssdna=graph.product("purified_ssdna","purified single-stranded chimera DNA",MaterialKind.MOLECULAR)
 hydrolysate=graph.product("hydrolysed_rna","alkali-hydrolysed RNA and SILANE flow-through",MaterialKind.WASTE)
 graph.react("hydrolyse_and_purify","Hydrolyse residual RNA at 98 °C, neutralize and purify on SILANE beads",
             StepKind.CLEANUP,((released,InputRole.TEMPLATE),),(ssdna,hydrolysate))
 graph.discard(hydrolysate,"not retained by SILANE bead purification")
 circligase=graph.input("circligase","CircLigase ssDNA ligase, 60 °C 4 h",MaterialKind.CATALYST)
 circle=graph.product("circular_ssdna","circularized single-stranded chimera",MaterialKind.MOLECULAR)
 graph.react("circularize","Circularize the released ssDNA",StepKind.CIRCULARIZATION,
             ((ssdna,InputRole.TEMPLATE),(circligase,InputRole.CATALYST)),(circle,),
             geometry=CircularizationGeometry.known("single-strand intramolecular ligation",
                 statement="Wu et al. 2019 Steps 89–98: the released ssDNA is 5′-phosphorylated on beads and circularized with CircLigase ssDNA ligase."))
 cut_oligo=graph.input("cut_oligo","BamHI Cut_oligo with 3′ inverted dT",MaterialKind.OLIGO)
 annealed=graph.product("cut_oligo_annealed","circle with a local Cut_oligo duplex over the linker BamHI site",MaterialKind.MOLECULAR)
 graph.react("anneal_cut_oligo","Anneal the Cut_oligo across the linker BamHI site",StepKind.HYBRIDIZATION,
             ((circle,InputRole.TEMPLATE),(cut_oligo,InputRole.STRUCTURAL)),(annealed,))
 bamhi=graph.input("bamhi","BamHI, 37 °C 1 h; SILANE bead purification",MaterialKind.CATALYST)
 linear=graph.product("relinearized_ssdna","left half-linker–DNA–cDNA–right half-linker ssDNA",MaterialKind.MOLECULAR)
 graph.react("bamhi_relinearize","Relinearize the circle with BamHI, splitting the linker",StepKind.DIGESTION,
             ((annealed,InputRole.TEMPLATE),(bamhi,InputRole.CATALYST)),(linear,))
 p5=graph.input("p5_primer","NEBNext Universal PCR primer",MaterialKind.OLIGO)
 p7=graph.input("p7_primer","NEBNext Index PCR primer (6-nt index unassigned)",MaterialKind.OLIGO)
 polymerase=graph.input("polymerase","NEBNext High-Fidelity 2× PCR Master Mix (12 cycles)",MaterialKind.CATALYST)
 amplified=graph.product("amplified_library",FINAL_LIBRARY.name+" (before size selection)",MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","PCR-amplify the relinearized ssDNA from its half-linker ends",StepKind.PCR,
             ((linear,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
              (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(amplified,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The compact Illumina branch records the completed library rather than separate PCR landing objects."))
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 out_of_range=graph.product("size_rejected","primer dimers and fragments outside 250–1000 bp",MaterialKind.WASTE)
 graph.react("size_selection","MinElute purify and cut 250–1000 bp from a 2% agarose E-gel",StepKind.CLEANUP,
             ((amplified,InputRole.TEMPLATE),),(library,out_of_range))
 graph.discard(out_of_range,"left in the gel outside the excised 250–1000 bp range")
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
     acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
     read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in SEQ_PRIMERS}))
 graph.declare_reviewed_scope((
     *graph.materials,
     *graph.steps),source)
 return graph.require_publishable()
