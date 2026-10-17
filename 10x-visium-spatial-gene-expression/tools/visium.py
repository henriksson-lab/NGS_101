"""Visium Spatial Gene Expression (fresh-frozen workflow) architecture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg, truseq_library
from chemdraw import Construct, Row, Scene, Segment, feature, oligo
from adapter_ligation import AdapterLigationGeometry
from illumina import TRUSEQ_READ1
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE = "10x Visium Spatial Gene Expression — array capture"
NOTES = "01_visium.html"
SOURCE = ('Commercial protocol: 10x Genomics <a href="https://cdn.10xgenomics.com/image/upload/v1660261286/support-documents/CG000239_Visium_Spatial_Gene_Expression_User_Guide_Rev_F.pdf">CG000239 Rev F</a>. '
          'The <a href="https://cdn.10xgenomics.com/raw/upload/v1655151898/support/in-line%20documents/Dual_Index_Kit_TT_Set_A.json">Dual Index TT Set A table</a> publishes the 96 well-specific i7/i5 pairs.')
SUMMARY = "A tissue section is permeabilized over an array of spatially indexed oligo-dT spots. Released polyadenylated RNA is copied on the slide so each cDNA receives the spot barcode and a UMI before standard short-read library construction."
CAVEAT = "The Rev F appendix prints the capture-oligo pattern and library primer/adapter strands. The TT Set A table publishes each index pair; the selected well and spot barcode depend on the physical slide and run. The UMI is random and slide-attachment chemistry remains unverified."
# CG000239 Rev F Appendix p. 69 prints the slide primer and cDNA forward primer as the
# 22-nt partial Read 1 CTACACGACGCTCTTCCGATCT (an earlier model slice took 23 nt).
SLIDE_READ1 = TRUSEQ_READ1[-22:]
TSO = "AAGCAGTGGTATCAACGCAGAGTACATrGrGrG"
SECOND_STRAND_PRIMER = "AAGCAGTGGTATCAACGCAGAG"
CDNA_PCR_FORWARD = SLIDE_READ1
CDNA_PCR_REVERSE = SECOND_STRAND_PRIMER
ADAPTER_TOP = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_BOTTOM = "TCTAGCCTTCTCG"[::-1]  # source prints this strand left-to-right as 3′→5′
INDEX_I5_PREFIX = "AATGATACGGCGACCACCGAGATCTACAC"
INDEX_I5_SUFFIX = "ACACTCTTTCCCTACACGACGCTC"
INDEX_I7_PREFIX = "CAAGCAGAAGACGGCATACGAGAT"
INDEX_I7_SUFFIX = "GTGACTGGAGTTCAGACGTGT"
SPATIAL=feature("spatial_barcode","spatial_barcode","whitelist",whitelist="10x-visium-spot-barcode-set")
UMI_FEATURE=feature("umi","umi","random")

FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([
    seg("spatial barcode", "B" * 16, "cbc", placeholder=True, feature=SPATIAL),
    seg("UMI", "U" * 12, "umi", placeholder=True, feature=UMI_FEATURE),
    seg("poly(dT)30", "T" * 30),
    seg("anchoring V", "V", placeholder=True), seg("anchoring N", "N", placeholder=True),
    seg("captured cDNA", "X" * 38, placeholder=True)],
    "Visium gene-expression library", dual_index=True, i5_nt=10, i7_nt=10)
FINAL_CAPTION = "Read 1 contains the 16-nt spatial barcode followed by the 12-nt UMI; Read 2 enters the captured transcript. The fixed adapter and primer arms follow the Rev F oligonucleotide appendix; spot and sample-index windows remain variable."
SEQUENCING_INTRO = "The 28-cycle Read 1 reports spot barcode plus UMI. Read 2 reports transcript sequence; both 10-cycle i7 and i5 reads identify the sample library."
READ_LENGTHS = {"Read 1": 28, "Index 1 (i7)": 10, "Index 2 (i5)": 10, "Read 2": 90}

def printed(name: str, *parts: tuple[str, str]) -> Construct:
    """One CG000239 Rev F appendix oligo as an exact whole strand (N windows as printed)."""
    return Construct([Segment(part, bases) for part, bases in parts], name=name)
def tso_state() -> Construct:
    # The appendix prints 5'-...ACATrGrGrG-3'; the three 3' riboguanosines are drawn as GGG.
    if not TSO.endswith("rGrGrG"):
        raise ValueError("Visium TSO no longer ends in the printed rGrGrG")
    return printed("template-switch oligo", ("DNA part", TSO[:-6]), ("rGrGrG", "GGG"))
def capture_scene() -> Scene:
    sc=Scene()
    sc.strand("mRNA", [seg("transcript", "X"*28, placeholder=True), seg("poly(A)", "A"*30)], label="tissue mRNA")
    primer=[seg("partial Read 1", SLIDE_READ1, "r1"),
            seg("spatial barcode", "B"*16, "cbc", placeholder=True, feature=SPATIAL),
            seg("UMI", "U"*12, "umi", placeholder=True, feature=UMI_FEATURE),
            seg("poly(dT)30", "T"*30), seg("V", "V", placeholder=True),
            seg("N", "N", placeholder=True)]
    sc.anneal("capture primer", primer, to="mRNA", pair=("poly(dT)30","poly(A)"),
              label="published slide-primer pattern", unpaired=("V","N"))
    sc.arrow("capture primer", "reverse transcription")
    sc.labels("capture primer")
    return sc

def adapter_scene() -> Scene:
    sc=Scene()
    sc.strand("bottom", [seg("paired stem", ADAPTER_BOTTOM[:-1], "r2"),
                          seg("3′ T", ADAPTER_BOTTOM[-1], "r2")], label="Visium adapter bottom, 5′→3′")
    sc.anneal("top", [seg("stem complement", ADAPTER_TOP[:12], "r2"),
                       seg("adapter top fork arm", ADAPTER_TOP[12:], "r2")],
              to="bottom", pair=("stem complement","paired stem"),
              label="Visium adapter top", unpaired=("adapter top fork arm",))
    sc.ligation_end("bottom", "3′ T", kind="cohesive", polarity="3-prime",
                    owner="adapter", text="published 3′-T ligation overhang")
    return sc

def oligos():
    return [oligo("slide capture primer pattern", [seg("partial Read 1",SLIDE_READ1),
           seg("spatial barcode","N"*16,placeholder=True),seg("UMI","N"*12,placeholder=True),
           seg("poly(dT)30","T"*30),seg("V","V",placeholder=True),seg("N","N",placeholder=True)]),
           f"<p>template-switch oligo: 5′- {TSO} -3′</p>",
           oligo("second-strand / cDNA PCR reverse primer",[seg("primer",SECOND_STRAND_PRIMER)]),
           oligo("cDNA PCR forward primer",[seg("primer",CDNA_PCR_FORWARD)]),
           oligo("ligation adapter top",[seg("top",ADAPTER_TOP)]),
           oligo("ligation adapter bottom (ordering orientation)",[seg("bottom",ADAPTER_BOTTOM)]),
           oligo("i5 indexing primer pattern",[seg("P5 arm",INDEX_I5_PREFIX),
                 seg("well-specific i5", "N"*10, placeholder=True),seg("3′ footprint",INDEX_I5_SUFFIX)]),
           oligo("i7 indexing primer pattern",[seg("P7 arm",INDEX_I7_PREFIX),
                 seg("well-specific i7", "N"*10, placeholder=True),seg("3′ footprint",INDEX_I7_SUFFIX)])]

def sections():
    return [
        ("Place, image and permeabilize tissue over the capture array",
         [Row(chunks=[("tissue section:       [cell][cell][cell]", None, False)]),
          Row(chunks=[("                         ↓ released poly(A) RNA", None, False)]),
          Row(chunks=[("capture area:    ● ● ● ● ● spatially barcoded spots", "cbc", False)])],
         "Each printed spot carries many copies of one spatial-barcode species. Imaging preserves the mapping from spot coordinates to tissue morphology."),
        ("Capture RNA and reverse-transcribe on the slide",
         capture_scene().rows(),
         "The surface oligo contributes spatial barcode, UMI and oligo-dT. Reverse transcription copies the locally captured RNA while it remains associated with its spot."),
        ("Release cDNA and construct the sequencing library",
         [Row(chunks=[("on-slide second strand → 0.08 M KOH denaturation → transfer → cDNA PCR", None, False)]),
          Row(chunks=[("→ fragmentation / end repair / dA-tail → adapter ligation → sample-index PCR", "me", False)])],
         "Second-strand synthesis happens on the slide; KOH denaturation then transfers the cDNA to a tube (CG000239 Rev F Steps 2.1–2.2). The spatial identity is already inside the cDNA before bulk recovery; downstream fragmentation must retain the barcode-bearing end used for Read 1."),
    ]

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    tissue=g.input("tissue","imaged tissue section containing polyadenylated RNA",MaterialKind.BIOLOGICAL,state=molecule("tissue section",(role("cells"),role("poly(A) RNA"))),primary=True)
    permeabilization=g.input("permeabilization","tissue permeabilization condition",MaterialKind.CONDITION)
    released=g.product("released_rna","locally released poly(A) RNA",MaterialKind.MOLECULAR,state=molecule("released RNA",(role("transcript"),role("poly(A)"))))
    g.react("permeabilize","Permeabilize tissue over the capture array",StepKind.EXTRACTION,((tissue,InputRole.TEMPLATE),(permeabilization,InputRole.CONDITION)),(released,))
    surface=g.input("surface_oligo","surface-bound spatial-barcode/UMI oligo(dT) with published sequence pattern",MaterialKind.OLIGO,state=tuple(capture_scene().rows()),inferred=True,source=SOURCE+" (Appendix; spot-specific barcode and surface chemistry remain variable)")
    tso=g.input("tso","template-switch oligo: "+TSO,MaterialKind.OLIGO,state=tso_state(),source=SOURCE+" (Appendix, PN-3000228)")
    rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
    surface_cdna=g.product("surface_cdna","spatial-barcode/UMI-tagged surface cDNA",MaterialKind.MOLECULAR,state=molecule("surface cDNA",(role("spatial barcode"),role("UMI"),role("captured cDNA"))))
    g.react("reverse_transcribe","Capture RNA and reverse-transcribe on the slide",StepKind.REVERSE_TRANSCRIPTION,((released,InputRole.TEMPLATE),(surface,InputRole.STRUCTURAL),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(surface_cdna,))
    second_primer=g.input("second_primer","second-strand synthesis primer: "+SECOND_STRAND_PRIMER,MaterialKind.OLIGO,state=printed("second-strand primer",("primer",SECOND_STRAND_PRIMER)),source=SOURCE+" (Appendix p. 69, Step 2.1, PN-2000217)")
    second_pol=g.input("second_polymerase","Second Strand Enzyme",MaterialKind.CATALYST,source=SOURCE+" (Step 2.1, PN-2000218/2000183)")
    surface_duplex=g.product("surface_duplex","double-stranded spatial cDNA on the capture area",MaterialKind.MOLECULAR,state=surface_cdna.state)
    g.react("second_strand","Synthesize the second cDNA strand on the slide",StepKind.EXTENSION,((surface_cdna,InputRole.TEMPLATE),(second_primer,InputRole.STRUCTURAL),(second_pol,InputRole.CATALYST)),(surface_duplex,),
            source=SOURCE+" (Step 2.1: Second Strand Mix added to the wells on the slide)")
    koh=g.input("koh","0.08 M KOH denaturation, 10 min room temperature, then 1 M Tris pH 7.0",MaterialKind.CONDITION,source=SOURCE+" (Step 2.2 a–g)")
    cdna=g.product("released_cdna","transferred spatially barcoded cDNA strand",MaterialKind.MOLECULAR,
                   state=molecule("transferred cDNA strand",(role("TSO"),role("captured cDNA"),role("poly(dA)BN"),role("UMI complement"),role("spatial barcode complement"),role("partial Read 1 complement"))),
                   inferred=True,source=SOURCE+" (Step 2.2 transfers the denatured eluate; Appendix p. 69 prints this second strand. The guide does not state slide attachment chemistry, so the identity of the released strand is inferred from the surface-primer origin of the first strand.)")
    slide=g.product("slide_residue","capture slide with surface-retained first strand",MaterialKind.WASTE)
    g.react("denature","Denature with KOH and transfer cDNA from each capture area",StepKind.STRAND_SEPARATION,((surface_duplex,InputRole.TEMPLATE),(koh,InputRole.CONDITION)),(cdna,slide),
            source=SOURCE+" (Step 2.2)")
    g.discard(slide,"Step 2.2h: the Visium Slide Cassette and slide may be discarded after transfer")
    wf=g.input("wta_forward","cDNA preamplification forward primer: "+CDNA_PCR_FORWARD,MaterialKind.OLIGO,state=printed("cDNA forward primer",("primer",CDNA_PCR_FORWARD)),source=SOURCE+" (Appendix, PN-2000089)")
    wr=g.input("wta_reverse","cDNA preamplification reverse primer: "+CDNA_PCR_REVERSE,MaterialKind.OLIGO,state=printed("cDNA reverse primer",("primer",CDNA_PCR_REVERSE)),source=SOURCE+" (Appendix, PN-2000089)")
    wp=g.input("wta_polymerase","cDNA preamplification polymerase",MaterialKind.CATALYST)
    amplified=g.product("amplified","amplified spatial cDNA",MaterialKind.MOLECULAR,state=surface_cdna.state)
    g.react("wta","Amplify spatially barcoded cDNA",StepKind.PCR,((cdna,InputRole.TEMPLATE),(wf,InputRole.STRUCTURAL),(wr,InputRole.STRUCTURAL),(wp,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(amplified.state,reason="The appendix prints both primers; exact binding coordinates vary with the captured transcript and are not represented by this role-level cDNA state."))
    spri_cdna=g.input("spri_cdna","SPRIselect cDNA cleanup",MaterialKind.SELECTION,source=SOURCE+" (Step 3.3)")
    cleaned=g.product("cleaned_cdna","SPRIselect-purified amplified cDNA",MaterialKind.MOLECULAR,state=amplified.state)
    cdna_waste=g.product("cdna_cleanup_waste","primers and short products removed by SPRIselect",MaterialKind.WASTE)
    g.react("cdna_cleanup","Clean up amplified cDNA with SPRIselect",StepKind.CLEANUP,((amplified,InputRole.TEMPLATE),(spri_cdna,InputRole.SELECTION)),(cleaned,cdna_waste),changes_molecule=False,source=SOURCE+" (Step 3.3)")
    g.discard(cdna_waste,"SPRIselect supernatant")
    frag=g.input("fragmentation","cDNA fragmentation condition",MaterialKind.CONDITION)
    fragments=g.product("fragments","fragmented barcode-bearing cDNA",MaterialKind.MOLECULAR,state=amplified.state)
    g.react("fragment","Fragment amplified cDNA while retaining barcode-bearing ends",StepKind.FRAGMENTATION,((cleaned,InputRole.TEMPLATE),(frag,InputRole.CONDITION)),(fragments,),source=SOURCE+" (Step 4.1 Fragmentation, End Repair & A-tailing)")
    endmix=g.input("endprep","end-repair and dA-tailing enzyme mix",MaterialKind.CATALYST)
    repaired=g.product("repaired","end-repaired cDNA fragments",MaterialKind.MOLECULAR,state=fragments.state)
    g.react("end_repair","Repair fragmented cDNA ends",StepKind.END_REPAIR,((fragments,InputRole.TEMPLATE),(endmix,InputRole.CATALYST)),(repaired,))
    tailed=g.product("tailed","dA-tailed cDNA fragments",MaterialKind.MOLECULAR,state=repaired.state)
    g.react("da_tail","Add one 3-prime dA to each repaired end",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(endmix,InputRole.CATALYST)),(tailed,))
    spri_frag=g.input("spri_frag","post-fragmentation double-sided SPRIselect size selection",MaterialKind.SELECTION,source=SOURCE+" (Step 4.2)")
    sized=g.product("sized","size-selected dA-tailed cDNA fragments",MaterialKind.MOLECULAR,state=tailed.state)
    frag_waste=g.product("frag_size_waste","fragments outside the double-sided SPRIselect window",MaterialKind.WASTE)
    g.react("frag_size_select","Double-sided SPRIselect size selection",StepKind.SELECTION,((tailed,InputRole.TEMPLATE),(spri_frag,InputRole.SELECTION)),(sized,frag_waste),changes_molecule=False,source=SOURCE+" (Step 4.2)")
    g.discard(frag_waste,"removed by the two SPRIselect cuts")
    adapter=g.input("adapter","Visium library adapter with published top and bottom strands",MaterialKind.ADAPTER,state=adapter_scene(),source=SOURCE+" (Appendix, PN-2000094)")
    ligase=g.input("ligase","library adapter ligase",MaterialKind.CATALYST)
    adapted=g.product("adapted","adapter-ligated spatial cDNA",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    g.react("ligate","Ligate sequencing adapters",StepKind.LIGATION,((sized,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),geometry=AdapterLigationGeometry.known("cohesive",adapter_scene(),statement="CG000239 Rev F appendix prints both adapter strands and the 3′-T end paired with a dA-tailed insert."))
    spri_lig=g.input("spri_lig","post-ligation SPRIselect cleanup",MaterialKind.SELECTION,source=SOURCE+" (Step 4.4)")
    ligated_clean=g.product("ligated_clean","purified adapter-ligated spatial cDNA",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    lig_waste=g.product("ligation_cleanup_waste","free adapter and ligation reagents removed by SPRIselect",MaterialKind.WASTE)
    g.react("ligation_cleanup","Clean up the ligation with SPRIselect",StepKind.CLEANUP,((adapted,InputRole.TEMPLATE),(spri_lig,InputRole.SELECTION)),(ligated_clean,lig_waste),changes_molecule=False,source=SOURCE+" (Step 4.4)")
    g.discard(lig_waste,"SPRIselect supernatant")
    reviewed=["tissue","permeabilization","permeabilize","surface_oligo","tso","rt","reverse_transcribe","second_primer","second_polymerase","second_strand","koh","denature","wta_forward","wta_reverse","wta_polymerase","wta","spri_cdna","cdna_cleanup","fragmentation","fragment","endprep","end_repair","da_tail","spri_frag","frag_size_select","adapter","ligase","ligate","spri_lig","ligation_cleanup"]
    finish_pcr_branch(g,"visium",ligated_clean,FINAL_LIBRARY,SEQ_PRIMERS,reviewed,label="Visium gene-expression",read_lengths=READ_LENGTHS,inferred_primers=False,limitation="The appendix prints both indexed primer patterns; the exact 10-nt indices depend on the selected Dual Index TT well.",
        # Appendix Step 4.5 prints each primer with an N10 well-specific index window.
        forward_state=printed("i5 indexing primer",("P5 arm",INDEX_I5_PREFIX),("i5 window","N"*10),("3' footprint",INDEX_I5_SUFFIX)),
        reverse_state=printed("i7 indexing primer",("P7 arm",INDEX_I7_PREFIX),("i7 window","N"*10),("3' footprint",INDEX_I7_SUFFIX)),
        primer_source=SOURCE+" (Appendix p. 71, Protocol Step 4.5, Dual Index TT Set A PN-1000215)")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
