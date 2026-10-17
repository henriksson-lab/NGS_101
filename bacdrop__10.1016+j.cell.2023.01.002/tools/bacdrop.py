"""BacDrop molecular model.

Source: Ma et al., Cell 2023 (PMC10014032 author manuscript JATS): STAR Methods
lines 177-217, Figure 1B image (nihms-1869349-f0002.jpg) and Key Resources Table
lines 1018-1083; 10x CG000209 Rev G line 3766 for the gel-bead oligo.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Row, Scene, Segment, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


# Key Resources Table (JATS lines 1018-1022), printed with spaces; spaces removed.
SMRT_DT = "AAGCAGTGGTATCAACGCAGAGTTTTTTTTTTTTTTTVN"
SMRT_PCR = "AAGCAGTGGTATCAACGCAGAGT"
P5_PARTIAL = "AATGATACGGCGACCACCGAGA"
# 10x CG000209 Rev G (line 3766): P5 - 16-nt barcode - partial Nextera Read 1 (S5).
BEAD_P5 = il.P5
BEAD_PARTIAL_READ1 = nx.S5
PLATE_BC_LEN, UMI_LEN, RANDOM_LEN, BEAD_BC_LEN = 13, 8, 6, 16
PBS1_REMAINDER = nx.READ1_PRIMER[len(nx.S5):]
PBS1_NOTE = ("Remainder of the Read 1 primer site assumed on the RT primer 5' handle; "
             "Figure 1B names it PBS1 and the bead primer's 3' S5 must anneal to it, "
             "but Table S1 (the RT primer sequences) is unavailable")


def _features(method):
    return (feature(method + "_plate", "cell_barcode", "combinatorial",
                    group=method + "_cell", part="plate well"),
            feature(method + "_droplet", "cell_barcode", "combinatorial",
                    group=method + "_cell", part="droplet"),
            feature(method + "_umi", "umi", "random"))


def indexed_random_rt_scene(*, method: str, barcode_len: int = PLATE_BC_LEN,
                            umi_len: int = UMI_LEN) -> Scene:
    """In-cell random-primed RT; Figure 1B order PBS1 - UMI - CB1 - random primer."""
    plate, _, umi = _features(method)
    rna = [seg("bacterial RNA", "X" * 34, placeholder=True)]
    primer = [
        seg("PBS1 (partial Read 1 site)", "H" * 16, "s5", placeholder=True),
        seg("UMI", "U" * umi_len, "umi", placeholder=True, feature=umi),
        seg("plate barcode CB1", "B" * barcode_len, "cbc", placeholder=True, feature=plate),
        seg("random primer", "N" * RANDOM_LEN, placeholder=True),
    ]
    sc = Scene(); sc.strand("RNA", rna, label="fixed, permeabilized bacterium")
    sc.anneal("RT primer", primer, to="RNA", pair=("random primer", "bacterial RNA"),
              label="well-indexed random RT primer (Table S1)",
              unpaired=("PBS1 (partial Read 1 site)", "UMI", "plate barcode CB1",
                        "random primer"))
    sc.arrow("RT primer", "in situ reverse transcription")
    sc.labels("RT primer")
    return sc


def droplet_barcode_rows(method: str) -> list[Row]:
    return [Row(chunks=[
        ("P5 — [CB2] — PBS1 — [UMI] — [CB1] — cDNA — A(n) — SMRT", "cbc", False)]),
        Row(chunks=[
            ("plate (CB1) + droplet (CB2) barcodes jointly identify the bacterium", None, False)])]


def bacterial_rna_library(method: str, name: str):
    """Figure 1B 'For Sequencing': P5-CB2-PBS1-UMI-CB1-cDNA-PBS2-Index1-P7."""
    plate, droplet, umi = _features(method)
    lib = Construct([
        seg("P5", BEAD_P5, "p5"),
        seg("droplet barcode CB2 (i5 position)", "J" * BEAD_BC_LEN, "cbc", placeholder=True,
            feature=droplet),
        seg("bead partial Read 1 (S5)", BEAD_PARTIAL_READ1, "s5"),
        seg("PBS1 remainder", PBS1_REMAINDER, "me", inferred=True, note=PBS1_NOTE),
        seg("UMI", "U" * UMI_LEN, "umi", placeholder=True, feature=umi),
        seg("plate barcode CB1", "B" * PLATE_BC_LEN, "cbc", placeholder=True, feature=plate),
        seg("random-primer copy", "N" * RANDOM_LEN, placeholder=True),
        seg("bacterial cDNA", "X" * 40, placeholder=True),
        seg("PBS2: Tn5 mosaic end (rc)", nx.ME_RC, "me",
            note="Nextera XT tagmentation end; kit adapter sequence not printed by the article"),
        seg("PBS2: S7 (rc)", nx.S7_RC, "s7",
            note="Nextera XT tagmentation end; kit adapter sequence not printed by the article"),
        seg("i7 reverse complement (N7 index)", "I" * 8, "cbc", placeholder=True,
            feature=feature("sample_i7", "sample_index", "unknown")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name=name)
    primers = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["I2"], sp.NEXTERA["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError(f"{name}: " + "; ".join(problems))
    return lib, primers



TITLE="BacDrop — droplet bacterial single-cell RNA-seq"
NOTES="01_bacdrop.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2023.01.002">Ma et al., <i>Cell</i> (2023)</a>.'
SUMMARY="Fixed bacteria receive a plate-specific RT index, are pooled at several cells per 10x droplet, and acquire a second bead barcode after TdT poly(A)-tailing of cDNA and in-droplet oligo-dT second-strand synthesis; Nextera XT tagmentation completes the library."
CAVEAT="The two-index cell identity, Figure 1B library order and the 10x bead primer's constant arms are published; the RT primer's PBS1 handle, individual bead barcode members and the 384 plate-indexed RT primer sequences remain unverified because the cited Table S1 is unavailable. The PBS1 remainder shown is inferred from bead-primer complementarity and standard Read 1 sequencing."
FINAL_LIBRARY,SEQ_PRIMERS=bacterial_rna_library("bacdrop","BacDrop RNA library")
FINAL_CAPTION="Figure 1B order: the droplet barcode sits in the i5 position, and Read 1 reads UMI, then plate barcode CB1, then the random-primed cDNA; Read 2 enters cDNA from the Nextera XT tagmentation end."
SEQUENCING_INTRO="NovaSeq 6000 with standard primers: Read 1 60, Index 1 8, Index 2 16, Read 2 39 cycles (STAR Methods)."
READ_LENGTHS={"Read 1":60,"Read 2":39,"Index 1 (i7)":8,"Index 2 (i5)":16}
SHOW_READ_CONTENT=True

def oligos():
 return [
  "<b>Well-indexed random RT primers (Table S1, 384)</b> — PBS1 – UMI (8) – CB1 (13) – random 6-mer (Figure 1B and legend); member sequences not available here.",
  "<b>SMRT_dT second-strand primer</b> — 5′-AAGCAGTGGTATCAACGCAGAGTTTTTTTTTTTTTTTVN-3′",
  "<b>SMRT_PCR cDNA enrichment primer</b> — 5′-AAGCAGTGGTATCAACGCAGAGT-3′",
  "<b>Partial P5 cDNA/library primer</b> — 5′-AATGATACGGCGACCACCGAGA-3′",
  '<b>10x PN-2000210 gel-bead oligo</b> — 5′-AATGATACGGCGACCACCGAGATCTACAC-[16-nt barcode]-TCGTCGGCAGCGTC-3′ (<a href="https://cdn.10xgenomics.com/image/upload/v1666737555/support-documents/CG000209_Chromium_NextGEM_SingleCell_ATAC_ReagentKits_v1.1_UserGuide_RevG.pdf">CG000209 Rev G, Appendix p. 51</a>; the barcode member is a placeholder).',
  "<b>Nextera index 1 primer (N7**, FC-131-2001)</b> — Illumina kit primer; sequence not printed by the article.",
 ]
INITIAL_ROWS=(Row(chunks=[("formaldehyde-fixed bacterial cells containing total RNA",None,False)]),)
INITIAL_NAME="Fixed bacterial cells"

def sections(): return [
 ("Permeabilize, deplete rRNA and remove gDNA in cells",[Row(chunks=[("bacterium: Tween-20/lysozyme — rRNA probes + thermostable RNase H — DNase I",None,False)])],"NEBNext bacterial rRNA probes and thermostable RNase H (50 °C, 30 min) reduce ribosomal RNA; DNase I then removes genomic DNA before indexing."),
 ("Install the first barcode by random-primed in-situ RT",indexed_random_rt_scene(method="bacdrop").rows(),"Each plate well contributes one CB1 and each copied molecule a UMI (Maxima H Minus RT; Figure 1B order PBS1–UMI–CB1–N6)."),
 ("Pool cells and TdT-tail their cDNA",[Row(chunks=[("pooled plate-indexed cells → TdT + dATP poly(A) tails on cDNA",None,False)])],"The poly(A) tail is added in cells (37 °C, 1 h; EDTA stop) before droplet generation."),
 ("Partition several bacteria per 10x droplet",[Row(chunks=[("cells + KAPA HiFi + SMRT_dT, 10x ATAC gel beads, Chip H", "cbc",False)])],"Combinatorial plate and droplet indexes permit super-loading."),
 ("Second-strand synthesis and droplet barcode capture",droplet_barcode_rows("bacdrop"),"SMRT_dT primes the poly(A) tail (39 °C 5 min, 65 °C 10 min); four cycles let the bead primer's 3′ S5 anneal to PBS1 and copy CB2 onto the cDNA."),
 ("Purify, enrich, tagment and index",[Row(chunks=[("Dynabeads + 0.6X AMPure → P5/SMRT_PCR enrichment → Nextera XT → P5/N7 PCR → 0.6X AMPure",None,False)])],"Only fragments retaining the bead-derived P5 end and a Nextera tagmentation end form the P5/N7 library."),
]


def _oligo(name, seq, **kw):
    return Construct([seg(name, seq, **kw)])


def protocol_graph():
 graph=ProtocolGraph(TITLE,authority=SOURCE)
 cells=graph.input("fixed_cells","bacteria fixed in 4% formaldehyde overnight at 4 °C, then 50% ethanol",
                   MaterialKind.BIOLOGICAL,state=INITIAL_ROWS,primary=True)
 permeabilizer=graph.input("permeabilizer","0.04% Tween-20 on ice, then 2.5 mg/mL lysozyme in Tris/EDTA, 37 °C 15 min",
                           MaterialKind.CONDITION)
 permeable=graph.product("permeable_cells","permeabilized fixed bacteria",MaterialKind.MOLECULAR)
 graph.react("permeabilize","Permeabilize fixed bacterial cells",StepKind.PERMEABILIZATION,
             ((cells,InputRole.TEMPLATE),(permeabilizer,InputRole.CONDITION)),(permeable,))
 probes=graph.input("rrna_probes","NEBNext Bacterial rRNA depletion solution (E7850)",MaterialKind.OLIGO)
 hybrid=graph.product("rrna_hybrid","rRNA–depletion-probe hybrids in cells",MaterialKind.MOLECULAR)
 graph.react("hybridize_rrna","Hybridize universal probes to abundant rRNA (50 °C 2 min, ramp to 22 °C)",
             StepKind.HYBRIDIZATION,
             ((permeable,InputRole.TEMPLATE),(probes,InputRole.STRUCTURAL)),(hybrid,))
 rnaseh=graph.input("rnase_h","Thermostable RNase H (E7850)",MaterialKind.CATALYST)
 depleted=graph.product("depleted_cells","rRNA-depleted permeabilized bacteria",
                        MaterialKind.MOLECULAR)
 rrna_waste=graph.product("rrna_waste","cleaved rRNA fragments",MaterialKind.WASTE)
 graph.react("deplete_rrna","Cleave probe-bound rRNA (50 °C, 30 min)",StepKind.DIGESTION,
             ((hybrid,InputRole.TEMPLATE),(rnaseh,InputRole.CATALYST)),
             (depleted,rrna_waste))
 graph.discard(rrna_waste,"cells are pelleted into the DNase reaction; released fragments leave with the supernatant")
 dnase=graph.input("dnase_i","DNase I (Sigma AMPD1), room temperature 30 min, EDTA stop at 50 °C",MaterialKind.CATALYST)
 dna_free=graph.product("dna_depleted_cells","rRNA- and gDNA-depleted bacteria",MaterialKind.MOLECULAR)
 gdna_waste=graph.product("gdna_waste","DNase I digestion products",MaterialKind.WASTE)
 graph.react("remove_gdna","Remove genomic DNA with DNase I",StepKind.DIGESTION,
             ((depleted,InputRole.TEMPLATE),(dnase,InputRole.CATALYST)),(dna_free,gdna_waste))
 graph.discard(gdna_waste,"cells are washed twice in PBS-RI")
 rt_primer=graph.input("indexed_rt_primer","Table S1 well-indexed random RT primers (PBS1–UMI–CB1–N6)",
                       MaterialKind.OLIGO)
 rt=graph.input("reverse_transcriptase","Maxima H Minus reverse transcriptase",MaterialKind.CATALYST)
 indexed_a=graph.product("indexed_cells_a","one group of plate-indexed cells",MaterialKind.MOLECULAR,
                         state=tuple(indexed_random_rt_scene(method="bacdrop").rows()))
 indexed_b=graph.product("indexed_cells_b","other plate-indexed cell groups",MaterialKind.MOLECULAR)
 graph.react("in_situ_rt","Install UMI and plate barcode by in-situ RT (multi-step 8–50 °C program)",
             StepKind.REVERSE_TRANSCRIPTION,
             ((dna_free,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
              (rt,InputRole.CATALYST)),(indexed_a,indexed_b))
 pooled=graph.product("pooled_cells","pooled plate-indexed bacteria",MaterialKind.MOLECULAR)
 graph.react("pool_wells","Pool all plate-indexed cell groups",StepKind.POOL,
             ((indexed_a,InputRole.TEMPLATE),(indexed_b,InputRole.TEMPLATE)),(pooled,),
             changes_molecule=False)
 tdt=graph.input("tdt","terminal transferase (NEB M0315L) with dATP",MaterialKind.CATALYST)
 tailed=graph.product("tailed_cdna","poly(A)-tailed indexed cDNA in pooled cells",MaterialKind.MOLECULAR)
 graph.react("tail_cdna","Add poly(A) tails to first-strand cDNA before droplet loading (37 °C 1 h)",StepKind.TAILING,
             ((pooled,InputRole.TEMPLATE),(tdt,InputRole.CATALYST)),(tailed,))
 droplets=graph.input("droplet_condition","10x Chip H with Single Cell ATAC gel beads and partitioning oil",MaterialKind.CONDITION)
 occupied=graph.product("occupied_droplets","droplets containing indexed bacteria and a bead",
                        MaterialKind.MOLECULAR)
 empty=graph.product("excluded_droplets","empty or unusable droplets",MaterialKind.WASTE)
 graph.react("partition","Partition pooled bacteria with barcode beads",StepKind.PARTITION,
             ((tailed,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,empty),
             changes_molecule=False)
 graph.discard(empty,"empty and otherwise unusable droplets yield no retained cell library")
 smrt_dt=graph.input("smrt_dt","SMRT_dT second-strand primer",MaterialKind.OLIGO,
                     state=Construct([seg("SMRT handle + dT",SMRT_DT[:-2]),
                                       seg("VN anchor","VN",placeholder=True,note="degenerate V then N")]),source=SOURCE+"; Key Resources Table")
 kapa=graph.input("kapa_hifi","KAPA HiFi HotStart ReadyMix",MaterialKind.CATALYST)
 second=graph.product("second_strand_cdna","oligo-dT-primed second strands in droplets",MaterialKind.MOLECULAR)
 graph.react("second_strand","Prime the poly(A) tail with SMRT_dT and synthesize the second strand (39 °C 5 min, 65 °C 10 min)",
             StepKind.EXTENSION,
             ((occupied,InputRole.TEMPLATE),(smrt_dt,InputRole.STRUCTURAL),
              (kapa,InputRole.CATALYST)),(second,))
 bead_primer=graph.input("bead_primer","10x PN-2000210 gel-bead primer",
                         MaterialKind.OLIGO,inferred=True,
                         source=SOURCE+"; 10x CG000209 Rev G, Appendix p. 51 (constant arms printed; barcode member unresolved)")
 barcoded=graph.product("dual_barcoded_cdna","plate- and droplet-barcoded cDNA",
                        MaterialKind.MOLECULAR,state=tuple(droplet_barcode_rows("bacdrop")))
 graph.react("transfer_droplet_barcode","Capture with the bead primer over four cycles (98/62/72 °C)",
             StepKind.EXTENSION,
             ((second,InputRole.TEMPLATE),(bead_primer,InputRole.STRUCTURAL),
              (kapa,InputRole.CATALYST)),(barcoded,))
 silane=graph.input("silane_ampure","Recovery Agent, Dynabeads MyOne SILANE, then 0.6X AMPure XP",MaterialKind.SELECTION)
 purified=graph.product("purified_cdna","purified droplet-barcoded cDNA",MaterialKind.MOLECULAR)
 cleanup_waste=graph.product("emulsion_waste","oils, recovery agent, bead supernatants and washes",MaterialKind.WASTE)
 graph.react("cdna_cleanup","Break emulsions and purify cDNA",StepKind.CLEANUP,
             ((barcoded,InputRole.TEMPLATE),(silane,InputRole.SELECTION)),(purified,cleanup_waste))
 graph.discard(cleanup_waste,"removed with the oil phases and bead supernatants")
 p5=graph.input("p5_partial","P5 partial primer",MaterialKind.OLIGO,
                state=_oligo("P5 partial",P5_PARTIAL,tag="p5"),source=SOURCE+"; Key Resources Table")
 smrt_pcr=graph.input("smrt_pcr","SMRT_PCR primer",MaterialKind.OLIGO,
                      state=_oligo("SMRT_PCR",SMRT_PCR),source=SOURCE+"; Key Resources Table")
 enriched=graph.product("enriched_cdna","P5/SMRT_PCR-enriched cDNA",MaterialKind.MOLECULAR)
 graph.react("cdna_enrichment","Enrich cDNA with P5 and SMRT_PCR (qPCR-determined cycles; 0.6X AMPure XP)",StepKind.PCR,
             ((purified,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),(smrt_pcr,InputRole.STRUCTURAL),
              (kapa,InputRole.CATALYST)),(enriched,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The enrichment amplicon's PBS1 handle and barcode members are not established; Table S1 is unavailable."))
 transposome=graph.input("nextera_xt_atm","Nextera XT Amplicon Tagment Mix adapters",MaterialKind.ADAPTER,
                         inferred=True,source=SOURCE+" (Nextera XT kit FC-131-1096; adapter strands not printed)")
 tn5=graph.input("tn5","Nextera XT transposase (ATM)",MaterialKind.CATALYST)
 tagmented=graph.product("tagmented","tagmented enriched cDNA (~2 ng per 50 µL reaction)",MaterialKind.MOLECULAR)
 graph.react("tagment","Tagment enriched cDNA with Nextera XT",StepKind.TAGMENTATION,
             ((enriched,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagmented,))
 n7=graph.input("pcr_reverse","Nextera Index 1 primer (N7**, FC-131-2001)",MaterialKind.OLIGO,inferred=True,
                source=SOURCE+"; Illumina kit primer, sequence not printed by the article")
 npm=graph.input("pcr_polymerase","Nextera PCR Master Mix (NPM)",MaterialKind.CATALYST)
 amplified=graph.product("indexed_pcr","P5/N7 library PCR product",MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","Index with P5 partial and N7 primers (5 cycles plus qPCR-determined cycles)",StepKind.PCR,
             ((tagmented,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
              (n7,InputRole.STRUCTURAL),(npm,InputRole.CATALYST)),(amplified,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The P5 partial primer is printed, but the N7 primer, Nextera XT adapter strands and the RT primer PBS1 handle are not printed by the inspected sources."))
 ampure=graph.input("ampure_xp","0.6X AMPure XP beads",MaterialKind.SELECTION)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 lib_waste=graph.product("library_supernatant","0.6X bead supernatant and washes",MaterialKind.WASTE)
 graph.react("library_cleanup","Purify the library with 0.6X AMPure XP",StepKind.CLEANUP,
             ((amplified,InputRole.TEMPLATE),(ampure,InputRole.SELECTION)),(library,lib_waste))
 graph.discard(lib_waste,"bead supernatant and washes are removed")
 graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),
                                       read_lengths=READ_LENGTHS))
 graph.declare_reviewed_scope(("fixed_cells","permeabilizer","permeabilize","rrna_probes",
   "hybridize_rrna","rnase_h","deplete_rrna","dnase_i","remove_gdna","indexed_rt_primer",
   "reverse_transcriptase","in_situ_rt","pool_wells","droplet_condition","partition","tdt",
   "tail_cdna","smrt_dt","kapa_hifi","second_strand","bead_primer","transfer_droplet_barcode",
   "silane_ampure","cdna_cleanup","p5_partial","smrt_pcr","cdna_enrichment","nextera_xt_atm",
   "tn5","tagment","pcr_reverse","pcr_polymerase","library_pcr","ampure_xp","library_cleanup",
   "library"),SOURCE)
 return graph.require_publishable()
