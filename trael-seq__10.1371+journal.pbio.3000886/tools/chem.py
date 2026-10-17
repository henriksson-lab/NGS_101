"""TrAEL-seq native 3-prime-end capture."""
import illumina as il
from batch_ngs import UNASSIGNED_LENGTH, truseq_library
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState, Row, Scene, Workflow, feature
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from topology import HairpinTopology
TITLE="TrAEL-seq — native DNA 3-prime-end sequencing"
NOTES="01_trael-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1371/journal.pbio.3000886">Kara et al., <i>PLOS Biology</i> (2021)</a>.'
SUMMARY="A pre-adenylated, biotinylated hairpin adapter ligates selectively to native DNA 3-prime OH ends; fill-in, fragmentation and affinity selection retain the captured end for directional Illumina sequencing."
UMI_FEATURE=feature("umi","umi","random")
# Read-strand order follows the Methods text and Fig 1A legend: Bst extends the
# hairpin's 3' arm across the 8-nt UMI, the TdT rA tail and the genomic strand, so
# Read 1 (primed on adapter 1) reads UMI, T(1-4) and then the genomic complement.
FINAL_LIBRARY,SEQ_PRIMERS=truseq_library(
    [seg("UMI (adapter 1 N8 copy)","N"*8,"umi",placeholder=True,feature=UMI_FEATURE),
     seg("T tract copied from 1-4 nt rA tail","TTT",placeholder=True),
     seg("genomic DNA complementary to the captured 3-prime end", "X" * 30, placeholder=True)],
    # E7335S (six-base i7) is only the paper's example set; the index length used is
    # not printed, so the i7 window carries no length.
    "TrAEL-seq", dual_index=False, i7_nt=UNASSIGNED_LENGTH,
    # The paper says Ultra II blunting; the dA follows from adapter 2's 3-prime T.
    inferred_junction=True)
ACQUIRED_ROLES=("Read 1","Index 1 (i7)")
FINAL_CAPTION=("Read 1 starts on adapter 1: 8-nt UMI, then the T tract copied from the 1-4 nt TdT rA tail, "
               "then genomic sequence complementary to the strand whose 3-prime end was captured. "
               "The dA junction reflects the published adapter 2 single 3-prime T overhang; the paper itself says the ends were blunted and polished with the Ultra II kit.")
SEQUENCING_INTRO=("The source reports 50-base single-read MiSeq or 75-base single-end "
                  "NextSeq 500 acquisition, so Read 2 is not acquired; the index-read length is not printed. "
                  "Amplification uses the NEBNext Universal primer and 'any NEBNext Index' primer from a "
                  "NEBNext Multiplex Oligos set, with E7335S given only as an example (a six-base i7 set); "
                  "the specific index set and primer used for each library are not named.")
CAVEAT=("Source wording: Fig 1A legend and Methods describe 'blunting' of sonicated ends with the NEBNext Ultra II kit, "
        "but the published adapter 2 ends in a single 3-prime T outside its 31-bp stem (see folded hairpin). The model "
        "draws the TA-compatible junction implied by that adapter and labels the dA tail as inferred; the Ultra II End Prep "
        "tailing activity was not checked in a vendor original in this pass.")
ADAPTER2 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCUUUUGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
ADAPTER1 = ("NNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTU"
            "GCGCAGGCCATTGGCCTGCGCUACACTCTTTCCCTACACGACGCT")
ADAPTER1_TOPOLOGY = HairpinTopology(
    "TrAEL adapter 1", ADAPTER1, first_stem_start=18,
    second_stem_start=64, stem_length=23, colour="var(--c-r1)")
ADAPTER2_TOPOLOGY = HairpinTopology(
    "TrAEL adapter 2", ADAPTER2, first_stem_start=0,
    second_stem_start=35, stem_length=31, colour="var(--c-r2)")
def oligos():
    return ["<b>TrAEL adapter 1</b> — 5′-[rApp]NNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTUGCGCAGGCCATTGGCC[BtndT]GCGCUACACTCTTTCCCTACACGACGCT-3′",f"<b>Adapter 2</b> — 5′-[Phos]{ADAPTER2[:-1]}*{ADAPTER2[-1]}-3′"]
def _adapter1_geometry():
    sc=Scene(); sc.strand("adapter 1",[
        seg("UMI","N"*8,"umi",placeholder=True,feature=feature("umi","umi","random")),
        seg("5′ adapter prefix","AGATCGGAAG"),
        seg("hairpin stem arm 1","AGCGTCGTGTAGGGAAAGAGTGT","r1"),
        seg("uracil/biotin loop","UGCGCAGGCCATTGGCCTGCGCU","w1"),
        seg("hairpin stem arm 2","ACACTCTTTCCCTACACGACGCT","r1")],label="published TrAEL adapter 1")
    sc.ligation_end("adapter 1","UMI",kind="single-strand",polarity="5-prime",owner="adapter",text="5′-App single-strand end")
    return AdapterLigationGeometry.known("single-strand",sc,statement="The protocol prints adapter 1 as a 5′-adenylated hairpin and joins that activated 5′ end to the TdT rA-tailed DNA 3′-OH with T4 RNA ligase 2 truncated KQ.")

def _adapter2_geometry():
    sc=Scene(); sc.strand("adapter 2",[
        seg("stem arm 1",ADAPTER2[:31],"r2"),
        seg("dU loop",ADAPTER2[31:35],"w1"),
        seg("stem arm 2",ADAPTER2[35:66],"r2"),
        seg("3′ T",ADAPTER2[66:])],label="published adapter 2 hairpin (5′-P; *T phosphorothioate)")
    sc.ligation_end("adapter 2","3′ T",kind="cohesive",polarity="3-prime",owner="adapter",text="single 3′ T overhang")
    return AdapterLigationGeometry.known("cohesive",sc,statement="The paper prints 5′-phosphorylated adapter 2 whose 31-bp stem leaves one 3′ T (phosphorothioate-linked) unpaired; it is ligated on beads with NEBNext Ultra II ligation mix after Ultra II End Prep.")

def _native_end():
    sc=Scene(); sc.strand("genomic strand",[seg("genomic DNA","X"*28,placeholder=True),seg("native 3′-OH","E",placeholder=True)],label="native DNA 3′ end in agarose")
    sc.labels("genomic strand"); return sc

def _tailed_end():
    sc=Scene(); sc.strand("genomic strand",[seg("genomic DNA","X"*28,placeholder=True),seg("1–4 nt rA tail","AAA",placeholder=True)],label="TdT-tailed 3′ end (ATP substrate)")
    sc.junction("genomic strand","genomic DNA","1–4 nt rA tail","TdT addition"); sc.labels("genomic strand"); return sc

def workflow():
    w=Workflow(MolecularState("native DNA 3′ end",tuple(_native_end().rows())))
    w.react("Tail native 3′ ends with terminal transferase and ATP (in agarose plugs)",_tailed_end().rows(),name="rA-tailed 3′ end",note="The paper reports 1–4 adenosines added to single-stranded 3′ ends (1–3 nt on an in-vitro substrate).")
    first=Scene(); first.strand("captured molecule",[seg("genomic DNA","X"*22,placeholder=True),seg("rA tail","AAA",placeholder=True),seg("adapter 1","X"*28,placeholder=True)],label="adapter-1-ligated end"); first.junction("captured molecule","rA tail","adapter 1","RNA-ligase junction")
    w.ligate_adapter("Ligate pre-adenylated adapter 1 with T4 RNA ligase 2 truncated KQ (PEG, 25 °C overnight), then wash plugs",_adapter1_geometry(),first.rows(),name="adapter-1-tagged end",note="Extensive plug washes remove unligated adapter 1; its 8 random bases are the UMI and internal biotin enables capture.")
    w.react("Digest agarose (β-agarase) and extend the hairpin with Bst 2.0 WarmStart (65 °C 30 min)",[Row(chunks=[("hairpin 3′ arm primes copying of UMI, rA tail and genomic strand → unnicked duplex","r1",False)])],name="tagged duplex",note="The paper cites strand-displacement and reverse-transcriptase activity of the thermophilic polymerase.")
    w.react("Sonicate (Covaris E220) and capture on MyOne C1 streptavidin beads",[Row(chunks=[("streptavidin ● biotin-tagged fragment; washes remove untagged fragments","w1",False)])],name="bead-bound tagged fragment",note="Capture precedes end polishing and adapter-2 ligation.")
    second=Scene.duplex([seg("tagged genomic fragment","X"*30,placeholder=True),seg("inferred dA","A",placeholder=True,inferred=True),seg("adapter 2","X"*20,placeholder=True)],label="two-adapter fragment on beads"); second.junction("top","inferred dA","adapter 2","adapter-2 ligation")
    w.ligate_adapter("INFERRED — on-bead Ultra II End Prep then adapter 2 ligation (paper says 'blunting'; adapter 2 has a 3′ T)",_adapter2_geometry(),second.rows(),name="two-adapter TrAEL molecule",note="Bead washes remove adapter-2 homodimers before elution.")
    w.react("Elute with USER (2 × 15 min, 37 °C) and amplify with NEBNext Universal + Index primers",[Row(chunks=[("USER cuts the dU sites flanking biotin and the adapter-2 loop → released fragment → indexed PCR","w1",False)])],name="indexed TrAEL library",note="Test PCR sets the cycle number; the final PCR is cleaned twice with AMPure XP.")
    return w

def protocol_graph():
    source="Kara et al., PLOS Biology 19:e3000886 (2021); doi:10.1371/journal.pbio.3000886"
    methods=source+"; Methods 'TrAEL-seq library preparation and sequencing' and S2 File (pbio.3000886.s009)"
    graph=ProtocolGraph("TrAEL-seq",authority=source)
    dna=graph.input("native_dna","agarose-embedded DNA bearing native 3-prime-OH ends",MaterialKind.BIOLOGICAL,
                    state=_native_end(),primary=True,source=methods)
    tdt=graph.input("tdt","Terminal Transferase (NEB M0315) with 0.4 mM ATP",MaterialKind.CATALYST,source=methods)
    tailed=graph.product("ra_tailed","3-prime ends carrying a 1-4 nt rA tail",MaterialKind.MOLECULAR,state=_tailed_end())
    graph.react("tailing","Tail 3-prime ends with TdT and ATP in the plug (37 C, 2 h)",StepKind.TAILING,
                ((dna,InputRole.TEMPLATE),(tdt,InputRole.CATALYST)),(tailed,),source=methods)
    adapter1=graph.input("adapter1","pre-adenylated biotinylated hairpin adapter 1",
                         MaterialKind.ADAPTER,source=methods)
    rna_ligase=graph.input("rna_ligase","T4 RNA ligase 2 truncated KQ (NEB M0373) with PEG 8000",MaterialKind.CATALYST,source=methods)
    tagged=graph.product("adapter1_tagged","adapter-1-tagged native DNA end",MaterialKind.MOLECULAR)
    graph.react("adapter1_ligation","Ligate pre-adenylated adapter 1 to rA-tailed 3-prime ends (25 C overnight)",
                StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adapter1,InputRole.STRUCTURAL),
                 (rna_ligase,InputRole.CATALYST)),(tagged,),geometry=_adapter1_geometry(),source=methods)
    plug_wash=graph.input("plug_washes","Tris buffer plug washes (4-5 x 10 ml, incl. overnight)",MaterialKind.SELECTION,source=methods)
    washed=graph.product("washed_plug","washed adapter-1-tagged DNA in agarose",MaterialKind.MOLECULAR)
    free_adapter=graph.product("free_adapter1","unligated adapter 1",MaterialKind.WASTE)
    graph.react("plug_wash","Wash plugs to remove unligated adapter 1",StepKind.CLEANUP,
                ((tagged,InputRole.TEMPLATE),(plug_wash,InputRole.SELECTION)),(washed,free_adapter),
                changes_molecule=False,source=methods)
    graph.discard(free_adapter,"plug wash buffer")
    agarase=graph.input("agarase","beta-agarase (NEB M0392) after melting at 65 C",MaterialKind.CATALYST,source=methods)
    melt=graph.input("plug_melt","agarase buffer, 65 C melt then 42 C",MaterialKind.CONDITION,source=methods)
    released=graph.product("released_dna","agarose-free adapter-1-tagged DNA",MaterialKind.MOLECULAR)
    graph.react("agarose_removal","Melt and digest agarose, then ethanol-precipitate",StepKind.EXTRACTION,
                ((washed,InputRole.TEMPLATE),(melt,InputRole.CONDITION),(agarase,InputRole.CATALYST)),(released,),
                changes_molecule=False,source=methods)
    bst=graph.input("bst","Bst 2.0 WarmStart DNA polymerase (NEB M0538)",MaterialKind.CATALYST,source=methods)
    filled=graph.product("filled_hairpin","hairpin-extended unnicked tagged duplex",MaterialKind.MOLECULAR)
    graph.react("hairpin_fill","Extend the adapter-1 hairpin across UMI, rA tail and genomic strand (65 C, 30 min)",
                StepKind.EXTENSION,
                ((released,InputRole.TEMPLATE),(bst,InputRole.CATALYST)),(filled,),source=methods)
    sonication=graph.input("sonication","Covaris E220 sonication (DF10, PIP175, 200 cycles)",MaterialKind.CONDITION,source=methods)
    fragmented=graph.product("fragmented","sonicated tagged duplex",MaterialKind.MOLECULAR)
    graph.react("sonicate","Sonicate the tagged duplex",StepKind.FRAGMENTATION,
                ((filled,InputRole.TEMPLATE),(sonication,InputRole.CONDITION)),(fragmented,),source=methods)
    beads=graph.input("streptavidin_beads","Dynabeads MyOne streptavidin C1",MaterialKind.SELECTION,source=methods)
    captured=graph.product("captured","bead-bound biotin-tagged fragments",MaterialKind.MOLECULAR)
    waste=graph.product("uncaptured","untagged DNA fragments",MaterialKind.WASTE)
    graph.react("capture","Capture adapter-1-tagged fragments on streptavidin beads",
                StepKind.AFFINITY_CAPTURE,
                ((fragmented,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(captured,waste),source=methods)
    graph.discard(waste,"unbound molecules are removed with the binding supernatant and 1 M NaCl / 0.1x TE washes")
    repair=graph.input("end_repair_mix","NEBNext Ultra II End Prep buffer/enzyme (+1 ng salmon sperm carrier)",MaterialKind.CATALYST,source=methods)
    blunt=graph.product("end_prepped","bead-bound end-prepped fragment (dA tail inferred)",MaterialKind.MOLECULAR)
    graph.react("end_repair","On-bead Ultra II End Prep (20 C 30 min, 65 C 30 min)",StepKind.END_REPAIR,
                ((captured,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(blunt,),source=methods)
    adapter2=graph.input("adapter2","annealed 5-prime-phosphorylated hairpin adapter 2",MaterialKind.ADAPTER,source=methods)
    ligase=graph.input("dna_ligase","NEBNext Ultra II ligation mix + ligation enhancer",MaterialKind.CATALYST,source=methods)
    dual=graph.product("dual_adapter_fragment","bead-bound two-adapter TrAEL molecule",MaterialKind.MOLECULAR)
    graph.react("adapter2_ligation","Ligate adapter 2 on beads (20 C, 30 min)",StepKind.LIGATION,
                ((blunt,InputRole.TEMPLATE),(adapter2,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(dual,),geometry=_adapter2_geometry(),source=methods)
    bead_wash=graph.input("bead_washes","1 M NaCl wash buffer and 0.1x TE bead washes",MaterialKind.SELECTION,source=methods)
    washed_beads=graph.product("washed_beads","washed bead-bound two-adapter molecules",MaterialKind.MOLECULAR)
    dimers=graph.product("adapter2_dimers","ligation mix including adapter-2 homodimers",MaterialKind.WASTE)
    graph.react("bead_wash","Discard ligation mix and wash beads",StepKind.CLEANUP,
                ((dual,InputRole.TEMPLATE),(bead_wash,InputRole.SELECTION)),(washed_beads,dimers),
                changes_molecule=False,source=methods)
    graph.discard(dimers,"ligation supernatant and bead washes")
    user=graph.input("user_enzyme","USER enzyme (NEB), 2 x 15 min at 37 C",MaterialKind.CATALYST,source=methods)
    opened=graph.product("user_opened","USER-eluted fragments",MaterialKind.MOLECULAR)
    graph.react("user_open","Elute by cleaving dU sites flanking biotin and in the adapter-2 loop",StepKind.DIGESTION,
                ((washed_beads,InputRole.TEMPLATE),(user,InputRole.CATALYST)),(opened,),source=methods)
    forward=graph.input("pcr_forward","NEBNext Universal PCR primer",MaterialKind.OLIGO,
                        source=methods+"; NEBNext Multiplex Oligos set, e.g. E7335S")
    reverse=graph.input("pcr_reverse","NEBNext Index primer (set and number not named)",MaterialKind.OLIGO,
                        source=methods+"; NEBNext Multiplex Oligos set, e.g. E7335S")
    polymerase=graph.input("pcr_polymerase","NEBNext Ultra II Q5 PCR master mix",MaterialKind.CATALYST,source=methods)
    amplified=graph.product("amplified",FINAL_LIBRARY.name+" (PCR product)",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    graph.react("indexed_pcr","Amplify (98 C 30 s; n x 98 C 10 s / 65 C 75 s; cycle number from a test PCR)",
                StepKind.PCR,
                ((opened,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(amplified,),
                geometry=PCRGeometry.exact(
                    FINAL_LIBRARY,
                    PCRPrimer("NEBNext Universal",il.NEBNEXT_UNIVERSAL_PRIMER),
                    PCRPrimer("NEBNext indexed primer",il.nebnext_i7_primer("N"*6)),
                    statement="The paper uses the NEBNext Universal primer and any NEBNext Index primer, giving E7335S as an example set; the six-base i7 window follows that example (E7335 manual). The set and index number used per library are not reported."),
                source=methods)
    ampure=graph.input("ampure","AMPure XP, two sequential cleanups (40 ul then 20 ul beads)",MaterialKind.SELECTION,source=methods)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    pcr_waste=graph.product("pcr_cleanup_waste","primers and short products",MaterialKind.WASTE)
    graph.react("pcr_cleanup","Clean the amplified library twice with AMPure XP",StepKind.CLEANUP,
                ((amplified,InputRole.TEMPLATE),(ampure,InputRole.SELECTION)),(library,pcr_waste),
                changes_molecule=False,source=methods)
    graph.discard(pcr_waste,"AMPure XP supernatants")
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
        acquired_roles=ACQUIRED_ROLES,
        read_uses={"Read 1":"include","Index 1 (i7)":"demultiplex"}))
    graph.declare_reviewed_scope((
        "native_dna","tdt","tailing","adapter1","rna_ligase","adapter1_ligation","plug_washes","plug_wash",
        "plug_melt","agarase","agarose_removal","bst","hairpin_fill",
        "sonication","sonicate","streptavidin_beads","capture","end_repair_mix","end_repair","adapter2","dna_ligase",
        "adapter2_ligation","bead_washes","bead_wash","user_enzyme","user_open",
        "pcr_forward","pcr_reverse","pcr_polymerase","indexed_pcr","ampure","pcr_cleanup","library"),source)
    return graph.require_publishable()
