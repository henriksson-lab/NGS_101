"""Xi Chen laboratory Modified SMART-Seq chemistry (Zheng et al. 2022, Methods Mol Biol 2520)."""
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind, Terminal)

TITLE="Modified SMART-Seq — 3′ counting with well barcodes and UMIs"
NOTES="01_modified-smart-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1007/7651_2021_435">Zheng et al., <i>Methods in Molecular Biology</i> (2022)</a>.'
SUMMARY=("A barcoded oligo(dT) installs an 8-nt well barcode and 10-nt UMI during reverse transcription in each well; "
         "template switching and cDNA PCR follow. Wells of one plate are then pooled, column-purified and Exonuclease I-treated "
         "before tagmentation and a plate-indexed library PCR.")
CAVEAT=("The Tn5 transposome is a choice of Illumina TDE1, other vendors or in-house Tn5 (Note 10); its adapter oligos are not printed in the protocol. "
        "The model draws the Nextera Read 2-type Tn5 end that the printed Library Index PCR Primer 3′ end requires. Fates of tagmentation "
        "fragments lacking the RT-primer handle are not described by the source and remain unresolved.")

# Table 1 oligos, 5'->3' (Zheng et al. 2022, p. 237).
RT="CTACACGACGCTCTTCCGATCT"                 # 5' handle of the Barcoded RT Primer = cDNA PCR Fwd Primer
TSO_DNA="AAGCAGTGGTATCAACGCAGAGTACAT"       # TSO DNA part; the protocol adds rGrGrG
CDNA_REV="AAGCAGTGGTATCAACGCAGAG"
UNIVERSAL="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC"
INDEX_5="CAAGCAGAAGACGGCATACGAGAT"          # before [8-bp plate barcode]
INDEX_3="GTCTCGTGGGCTCGG"                   # after [8-bp plate barcode]


def _overlap(primer, handle):
    """Longest suffix of the primer that equals a prefix of the handle."""
    for n in range(min(len(primer),len(handle)),0,-1):
        if primer[-n:]==handle[:n]:
            return n
    raise ValueError("library universal primer does not overlap the RT handle")


if not UNIVERSAL.startswith(il.P5):
    raise ValueError("Library PCR Universal Primer does not begin with P5")
_OV=_overlap(UNIVERSAL,RT)
PRIMER_ADDED=UNIVERSAL[len(il.P5):len(UNIVERSAL)-_OV]   # tail added 5' of the RT handle
if revcomp(INDEX_3)!=nx.S7_RC or revcomp(INDEX_5)!=il.P7_RC:
    raise ValueError("Library Index PCR Primer arms are not the Nextera i7 adapter / P7")

FINAL_LIBRARY=Construct([
    seg("P5",il.P5,"p5"),
    seg("universal-primer tail",PRIMER_ADDED,"r1"),
    seg("RT-primer handle",RT,"r1"),
    seg("well barcode","B"*8,"cbc",placeholder=True,feature=feature(
        "well_barcode","cell_barcode","whitelist",
        whitelist="Zheng et al. 2022, Table 2 well-barcode set")),
    seg("UMI","U"*10,"umi",placeholder=True,feature=feature("umi","umi","random")),
    seg("oligo(dT)","T"*25),
    seg("VN anchor","VN",placeholder=True),
    seg("transcript 3′ cDNA","X"*36,placeholder=True),
    seg("mosaic end reverse complement",nx.ME_RC,"me"),
    seg("Nextera i7 adapter (from index primer 3′ arm)",revcomp(INDEX_3),"s7"),
    seg("plate index reverse complement","I"*8,"cbc",placeholder=True,feature=feature(
        "plate_index","sample_index","whitelist",
        whitelist="Zheng et al. 2022, Table 3 plate-index set")),
    seg("P7 reverse complement",revcomp(INDEX_5),"p7")],name="Modified SMART-Seq library")
SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.NEXTERA["I1"],sp.NEXTERA["R2"])
FINAL_CAPTION=("Read 1 begins with the well barcode and UMI, then the oligo(dT); Read 2 reads transcript sequence from the tagmented end. "
               "The TSO-derived 5′ cDNA end is not part of this 3′ library structure.")
SEQUENCING_INTRO=("Section 3.8 step 11 lists Read 1 ≥18 cycles (well barcode and UMIs), Index 1 8 cycles (plate barcode), "
                  "Index 2 8 cycles (identity N/A) and Read 2 ≥50 cycles (cDNA). Note 13 calls Index 2 optional because the left-hand "
                  "side carries no index, and says Read 2 is normally 75 cycles; Index 2 is therefore not drawn as an informative read.")
READ_LENGTHS={"Read 1":18,"Index 1 (i7)":8,"Read 2":75}
SHOW_READ_CONTENT=True


def oligos(): return [
 "<b>Barcoded RT primer</b> — 5′-CTACACGACGCTCTTCCGATCT[8-bp well barcode]NNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTVN-3′",
 "<b>Template-switching oligo</b> — 5′-AAGCAGTGGTATCAACGCAGAGTACATrGrGrG-3′ (rG = riboguanosine)",
 f"<b>cDNA PCR forward / reverse</b> — 5′-{RT}-3′ / 5′-{CDNA_REV}-3′",
 f"<b>Library PCR universal</b> — 5′-{UNIVERSAL}-3′",
 f"<b>Library PCR index</b> — 5′-{INDEX_5}[8-bp plate barcode]{INDEX_3}-3′"]


def _oligo(name,seq):
    return Construct([seg(name,seq)],name=name)


def protocol_graph():
 source="Defining source: https://doi.org/10.1007/7651_2021_435"
 graph=ProtocolGraph("Modified SMART-Seq",authority=source)
 cell=graph.input("cell","single sorted cell in 2 µl PEG-8000/Triton X-100/RNase inhibitor/dNTP lysis buffer",MaterialKind.BIOLOGICAL,primary=True)
 rt_primer=graph.input("barcoded_rt_primer","Barcoded RT Primer, 1 µl of 2 µM per well (one well barcode per well)",MaterialKind.OLIGO)
 heat=graph.input("lysis_denature","72 °C 10 min",MaterialKind.CONDITION)
 primed=graph.product("primed_lysate","lysed cell RNA with barcoded oligo(dT)",MaterialKind.MOLECULAR)
 graph.react("lyse_denature","Lyse/denature with the barcoded RT primer in the well",StepKind.THERMAL_TREATMENT,
             ((cell,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(heat,InputRole.CONDITION)),(primed,))
 rt=graph.input("reverse_transcriptase","Maxima H-minus RT (8 U/µl) in Tris/NaCl/MgCl2/GTP/DTT/RNase inhibitor mix",MaterialKind.CATALYST)
 tso=graph.input("tso","Template Switching Oligo with 3′ rGrGrG (8 µM in RT mix)",MaterialKind.OLIGO)
 program=graph.input("rt_program","42 °C 90 min; 10 × (50 °C 2 min, 42 °C 2 min); 85 °C 5 min",MaterialKind.CONDITION)
 switched=graph.product("switched_cdna","well-barcode/UMI-tagged, template-switched first-strand cDNA",MaterialKind.MOLECULAR)
 graph.react("rt_template_switch","Reverse-transcribe and template-switch in one RT program",StepKind.REVERSE_TRANSCRIPTION,
             ((primed,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST),(program,InputRole.CONDITION)),(switched,))
 cdna_forward=graph.input("cdna_pcr_forward","Table 1 cDNA PCR Fwd Primer",MaterialKind.OLIGO,state=_oligo("cDNA PCR Fwd Primer",RT))
 cdna_reverse=graph.input("cdna_pcr_reverse","Table 1 cDNA PCR Rev. Primer",MaterialKind.OLIGO,state=_oligo("cDNA PCR Rev. Primer",CDNA_REV))
 polymerase=graph.input("pcr_polymerase","KAPA HiFi HotStart ReadyMix",MaterialKind.CATALYST)
 amplified=graph.product("amplified_cdna","well-tagged amplified full-length cDNA (per well)",MaterialKind.MOLECULAR)
 graph.react("cdna_pcr","Amplify cDNA in each well (98 °C 5 min; 16 × 98/65/72 °C; 72 °C 5 min)",StepKind.PCR,
             ((switched,InputRole.TEMPLATE),(cdna_forward,InputRole.STRUCTURAL),
              (cdna_reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(amplified,),
             geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                 reason="The cDNA-PCR intermediate is not represented base-for-base.",
                 symmetric_species=True))
 column=graph.input("qiaquick","5 vol Buffer PB, one QIAquick column per plate, 80% ethanol wash, 50 µl elution",MaterialKind.SELECTION)
 eluted=graph.product("purified_pool","column-purified plate pool of well-tagged cDNA",MaterialKind.MOLECULAR)
 flow=graph.product("column_flowthrough","column flow-through and wash",MaterialKind.WASTE)
 graph.react("pool_column_purify","Pool all wells of one plate (plates kept separate; Note 5) and purify on one QIAquick column",StepKind.CLEANUP,
             ((amplified,InputRole.TEMPLATE),(column,InputRole.SELECTION)),(eluted,flow),changes_molecule=False)
 graph.discard(flow,"Column flow-through and washes are discarded.")
 exo=graph.input("exonuclease_i","Exonuclease I (1 U/µl, NEBuffer 3.1), 37 °C 30 min, 80 °C 2 min",MaterialKind.CATALYST)
 digested=graph.product("exo_treated","cDNA with residual single-stranded primers digested",MaterialKind.MOLECULAR)
 graph.react("exo_digest","Digest leftover primers with Exonuclease I after pooling",StepKind.DIGESTION,
             ((eluted,InputRole.TEMPLATE),(exo,InputRole.CATALYST)),(digested,))
 spri1=graph.input("spri_07","0.7 vol SPRI beads, three 80% ethanol washes",MaterialKind.SELECTION)
 cleaned=graph.product("cleaned_cdna","purified full-length cDNA of one plate",MaterialKind.MOLECULAR)
 waste=graph.product("spri1_supernatant","bead supernatant with primer digestion products and short material",MaterialKind.WASTE)
 graph.react("spri_cdna","Purify full-length cDNA with 0.7× SPRI beads",StepKind.CLEANUP,
             ((digested,InputRole.TEMPLATE),(spri1,InputRole.SELECTION)),(cleaned,waste),changes_molecule=False)
 graph.discard(waste,"The 0.7× bead supernatant is removed.")
 transposome=graph.input("transposome","Tn5 transposome (TDE1, other vendor or in-house; Note 10) with 2× TD buffer; adapter oligos not printed",MaterialKind.ADAPTER)
 tn5=graph.input("tn5","Tn5 transposase, 50 ng cDNA, 55 °C 5 min, stopped with 0.2% SDS at 55 °C 5 min",MaterialKind.CATALYST)
 tagged=graph.product("tagmented_cdna","tagmented well-barcoded cDNA",MaterialKind.MOLECULAR)
 offtarget=graph.product("tagmentation_off_products","Tn5–Tn5 and Tn5–TSO-end fragments lacking the RT-primer handle",MaterialKind.MOLECULAR)
 graph.react("tagment","Tagment pooled cDNA",StepKind.TAGMENTATION,
             ((cleaned,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),
             (tagged,offtarget))
 graph.mark_terminal(Terminal(offtarget.id,Disposition.UNRESOLVED,
     "These fragments lack the RT-primer handle bound by the universal primer and so cannot form the declared P5 library end; "
     "their amplification by the index primer alone, carryover or removal is not described by the source."))
 spri2=graph.input("spri_08_tag","0.8 vol SPRI beads after tagmentation",MaterialKind.SELECTION)
 tagged_clean=graph.product("tagmented_clean","purified tagmented cDNA",MaterialKind.MOLECULAR)
 waste2=graph.product("spri2_supernatant","post-tagmentation bead supernatant",MaterialKind.WASTE)
 graph.react("spri_tagmented","Purify tagmented cDNA with 0.8× SPRI beads",StepKind.CLEANUP,
             ((tagged,InputRole.TEMPLATE),(spri2,InputRole.SELECTION)),(tagged_clean,waste2),changes_molecule=False)
 graph.discard(waste2,"The 0.8× bead supernatant is removed.")
 forward=graph.input("library_forward","Library PCR Universal Primer (0.5 µM)",MaterialKind.OLIGO,state=_oligo("Library PCR Universal Primer",UNIVERSAL))
 reverse=graph.input("indexed_reverse","Library PCR Index Primer (0.5 µM; Nextera XT-type, one plate index per plate)",MaterialKind.OLIGO)
 raw_library=graph.product("raw_library",FINAL_LIBRARY.name+" (unpurified PCR)",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
 graph.react("library_pcr","Gap-fill and amplify RT-handle/Tn5 fragments, adding the plate index (72 °C 5 min; 98 °C 2 min; 8 cycles)",StepKind.PCR,
             ((tagged_clean,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
              (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(raw_library,),
             geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                 reason="The tagmented pre-PCR fragment is variable; the final primer sites are canonical."))
 spri3=graph.input("spri_08_lib","0.8 vol SPRI beads after library PCR (repeat if primer leftover persists; Fig. 7)",MaterialKind.SELECTION)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 waste3=graph.product("spri3_supernatant","library bead supernatant with primer leftover",MaterialKind.WASTE)
 graph.react("spri_library","Purify the library with 0.8× SPRI beads",StepKind.CLEANUP,
             ((raw_library,InputRole.TEMPLATE),(spri3,InputRole.SELECTION)),(library,waste3),changes_molecule=False)
 graph.discard(waste3,"The 0.8× bead supernatant is removed.")
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
     acquired_roles=tuple(p.role for p in SEQ_PRIMERS),read_lengths=READ_LENGTHS,
     read_uses={"Read 1":"extract","Index 1 (i7)":"demultiplex","Read 2":"include"}))
 graph.declare_reviewed_scope(tuple(graph.materials)+tuple(graph.steps),source)
 return graph.require_publishable()
