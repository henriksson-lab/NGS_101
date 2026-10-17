"""ICICLE-seq molecular model.

Authority: Swanson et al., eLife 2021;10:e63632 (PMC8034981 JATS): Methods "ICICLE-seq"
(oligonucleotide sources, Tn5 complexing, antibody staining, library preparation,
sequencing), Figure 3-figure supplement 1 legend (a-e), and Supplementary file 5
("Table 5 - ICICLE Oligos").  The 10x 3' v3.1 gel-bead oligo is a vendor reagent whose
sequence is not printed in those sources; it is drawn as an inferred scaffold.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Row, Segment, feature, oligo, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

SOURCE='<a href="https://doi.org/10.7554/eLife.63632">Swanson et al., <i>eLife</i> (2021)</a>.'
TABLE5=SOURCE+" Supplementary file 5 (Table 5, ICICLE Oligos)"
FIG_S1=SOURCE+" Figure 3—figure supplement 1 legend"

# Supplementary file 5 oligos.  Canonical sequences are imported and the selftest checks
# each complete base sequence against the extracted original table.
MOSAIC_BOT=nx.ME_RC                       # /5phos/, standard desalting
TN5ME_S7_TOP=nx.ADAPTOR_S7                 # "Compatible with Nextera i7 primers"
POLY_T_TOP="T"*30                         # followed by VN and the mosaic end
SI_P5_22=il.TRUSEQ_P5_FULL                # "TruSeq P5 primer"
F_BC_PRIMER=il.TRUSEQ_READ1[11:]          # "Forward"
ADT_REV_AMP="CCTTGGCACCCGAGAATTCC"        # 3' part of the ADT-i7 arm (selftest checks)
ADT_I7_ARM="GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"  # ...CC*A: final two bases phosphorothioate-linked
TRUSEQ_R1_SEQ=il.TRUSEQ_READ1             # spiked in as Custom Read 1 primer
if not ADT_I7_ARM.startswith("GTGACTGGAGTTC") or ADT_REV_AMP not in ADT_I7_ARM:
    raise ValueError("ADT-Rev-AMP must lie within the ADT-i7 binding arm")

CELL=feature("icicle_cell","cell_barcode","whitelist",whitelist="10x v3 3' GEX whitelist 3M-february-2018")
UMI=feature("icicle_umi","umi","random")
ATAC_I7=feature("icicle_atac_i7","sample_index","unknown")
ADT_I7=feature("icicle_adt_i7","sample_index","unknown")
ANTIBODY=feature("icicle_antibody","feature_barcode","whitelist",whitelist="TotalSeq-A antibody-tag panel (Supplementary file 6)")


def bead_scaffold():
    """10x 3' v3.1 bead oligo: inferred vendor scaffold, not printed in the article.

    Read 1 = 28 bp (cell barcode and UMI) and the v3 16-base whitelist fix the
    16 + 12 split; the partial Read 1 handle is implied by the F BC Primer."""
    return [Segment("partial TruSeq Read 1",F_BC_PRIMER,"r1",inferred=True),
            Segment("cell barcode","N"*16,"cbc",placeholder=True,inferred=True,feature=CELL),
            Segment("UMI","N"*12,"umi",placeholder=True,inferred=True,feature=UMI),
            Segment("oligo(dT)","T"*30,inferred=True),
            Segment("VN anchor","VN",placeholder=True,inferred=True,length_bp=2)]


def atac_precursor():
    # Fig. 3-S1b-c: bead extension across a Poly-T/R2N fragment, then F BC Primer copies.
    return Construct([*bead_scaffold(),
        Segment("mosaic end",nx.ME,"me"),
        Segment("accessible genomic DNA","X"*40,placeholder=True),
        Segment("mosaic end reverse complement",nx.ME_RC,"me"),
        Segment("Tn5ME-s7 R2N reverse complement",revcomp(TN5ME_S7_TOP[:-len(nx.ME)]),"s7")],
        name="ICICLE-seq barcoded ATAC precursor")


def adt_precursor():
    # Fig. 3-S1b-c: bead extension across a poly(A) TotalSeq-A tag, then ADT-Rev-AMP PCR.
    return Construct([*bead_scaffold(),
        Segment("antibody barcode reverse complement","N"*15,"cbc",placeholder=True,inferred=True,feature=ANTIBODY),
        Segment("copied ADT handle",revcomp(ADT_I7_ARM[len(ADT_I7_ARM)-len(ADT_REV_AMP)-1:]),"r2")],
        name="ICICLE-seq barcoded ADT precursor")


def atac_library():
    return Construct([Segment("P5",il.P5[:-4],"p5"),Segment("TruSeq Read 1",il.TRUSEQ_READ1,"r1"),
        *atac_precursor().segments[1:-1],
        Segment("Nextera Read 2 arm",revcomp(TN5ME_S7_TOP[:-len(nx.ME)]),"s7"),
        Segment("i7 reverse complement","N"*8,"cbc",placeholder=True,feature=ATAC_I7),
        Segment("P7 reverse complement",il.P7_RC,"p7")],name="ICICLE-seq ATAC library")


def adt_library():
    return Construct([Segment("P5",il.P5[:-4],"p5"),Segment("TruSeq Read 1",il.TRUSEQ_READ1,"r1"),
        *adt_precursor().segments[1:-1],
        Segment("small-RNA Read 2 arm",revcomp(ADT_I7_ARM),"r2"),
        Segment("i7 reverse complement","N"*8,"cbc",placeholder=True,feature=ADT_I7),
        Segment("P7 reverse complement",il.P7_RC,"p7")],name="ICICLE-seq antibody-tag library")


ATAC=atac_library()
ADT=adt_library()
R1_CUSTOM=sp.custom("Read 1","TruSeq R1 Seq Primer (Supplementary file 5, Custom Read 1)",TRUSEQ_R1_SEQ,TABLE5,
    "Spiked in at 0.3 µM to avoid off-target Nextera Read 1 priming on the partial R1N in the Poly-T mosaic end.")
ATAC_P=(R1_CUSTOM,sp.NEXTERA["I1"],sp.NEXTERA["R2"])
ADT_P=(R1_CUSTOM,)
for _lib,_primers in ((ATAC,ATAC_P),(ADT,ADT_P)):
    _problems=sp.verify(_lib,_primers,required_roles=tuple(p.role for p in _primers))
    if _problems:raise ValueError("; ".join(_problems))

TITLE="ICICLE-seq — chromatin accessibility and epitopes"
NOTES="01_icicle-seq.html"
SUMMARY=("Tn5 loaded with a 1:1 mix of Poly-T and Nextera R2N (s7) mosaic-end duplexes tagments permeabilized, "
         "TotalSeq-A-stained cells. In 10x 3′ v3.1 GEMs the poly(T) end and polyadenylated antibody tags anneal "
         "to the bead oligo; ATAC is copied linearly with F BC Primer and ADTs exponentially with ADT-Rev-AMP, "
         "then split 3:1 into separately indexed ATAC and ADT libraries.")
CAVEAT=("<strong>INFERRED regions:</strong> the 10x 3′ v3.1 gel-bead oligo (partial Read 1, 16-base cell barcode, "
        "12-base UMI, oligo(dT)30VN) is not printed in the article; its 16 + 12 split is derived from the 28-cycle "
        "Read 1 and the v3 16-base whitelist. TotalSeq-A barcode length and the Chromium i7 Multiplex Kit N primer "
        "sequence are vendor values not printed here. Figure 3—figure supplement 1a states that only about half of "
        "the tagmented fragments carry one Poly-T and one R2N end; same-end fragments stay in the mixture and are "
        "not drawn as a separated product. The ADT Read 2/Index 1 entry is named only as a "
        "\"TruSeq small RNA/i7 Read 2 priming site\"; the loaded primer sequence is not asserted.")
FINAL_LIBRARIES=(("ATAC library",ATAC,ATAC_P,"Single-end ATAC insert read from the R2N end; Read 1 gives cell barcode and UMI.",
                  "Custom TruSeq Read 1 reads cell barcode + UMI (28); i7 via the Nextera Read 2 site (8); Read 2 reads gDNA (100)."),
                 ("Antibody-tag library",ADT,ADT_P,"Antibody barcode amplified with SI-P5-22 and ADT-i7.",
                  "Custom TruSeq Read 1 reads cell barcode + UMI; Read 2 identifies the antibody."))
READ_LENGTHS={"ATAC library":{"Read 1":28,"Index 1 (i7)":8,"Read 2":100},
              "Antibody-tag library":{"Read 1":28,"Index 1 (i7)":8,"Read 2":100}}
ADT_LIMITATION=("Figure 3—figure supplement 1e names a TruSeq small RNA/i7 Read 2 priming site for ADT Index 1 and Read 2, "
                "but the loaded sequencing-primer sequence is not printed; those binding sites and cycle coordinates are not asserted.")
READ_SPECS={"Antibody-tag library":(
    sp.read_spec("Read 1",28,"16-base cell barcode and 12-base UMI","extract"),
    sp.read_spec("Index 1 (i7)",8,"ADT-i7 sample index","demultiplex"),
    sp.read_spec("Read 2",100,"antibody barcode; exact primer entry unresolved","extract"))}
SEQUENCING_UNAVAILABLE={"Antibody-tag library":ADT_LIMITATION}


def oligos():
    """Supplementary file 5 oligos, verbatim with modifications."""
    yield oligo("Poly-A Top-L",[Segment("poly(T)",POLY_T_TOP),Segment("VN","VN",placeholder=True),Segment("mosaic end",nx.ME,"me")])
    yield oligo("Tn5ME-s7_Top",[Segment("Nextera R2N (s7)",TN5ME_S7_TOP[:-len(nx.ME)],"s7"),Segment("mosaic end",nx.ME,"me")])
    yield oligo("MOSAIC_Bot",[Segment("mosaic end reverse complement",MOSAIC_BOT,"me")],mods="/5phos/")
    yield oligo("F BC Primer",[Segment("partial TruSeq Read 1",F_BC_PRIMER,"r1")])
    yield oligo("ADT-Rev-AMP",[Segment("ADT handle primer",ADT_REV_AMP,"r2")])
    yield oligo("SI-P5-22",[Segment("P5",il.P5[:-4],"p5"),Segment("TruSeq Read 1",il.TRUSEQ_READ1,"r1")])
    yield oligo("ADT-i7 Primer",[Segment("P7",il.P7,"p7"),Segment("8-base i7 index","N"*8,placeholder=True),
                                 Segment("ADT-binding arm",ADT_I7_ARM[:-2],"r2"),Segment("phosphorothioate-linked 3′ end","C*A",placeholder=True)])
    yield oligo("TruSeq R1 Seq Primer",[Segment("TruSeq Read 1",TRUSEQ_R1_SEQ,"r1")])


def sections():return []


def _cleanup(g,ident,label,template,condition,reason,kept_label,kind=None,state=None):
    sel=g.input(ident+"_beads",condition,MaterialKind.SELECTION)
    kept=g.product(ident+"_kept",kept_label,kind or template.kind,state=template.state if state is None else state)
    waste=g.product(ident+"_waste","removed fraction ("+condition+")",MaterialKind.WASTE)
    g.react(ident,label,StepKind.CLEANUP,((template,InputRole.TEMPLATE),(sel,InputRole.SELECTION)),(kept,waste),changes_molecule=False)
    g.discard(waste,reason)
    return kept,[sel.id,kept.id,waste.id,ident]


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 reviewed=[]
 # Tn5 complexing (Methods; Mulqueen et al. 2019 procedure is cited, not re-read here).
 polyt=g.input("poly_t_top","Poly-A Top-L",MaterialKind.OLIGO,state=Construct([Segment("poly(T)",POLY_T_TOP),Segment("VN","VN",placeholder=True,length_bp=2),Segment("mosaic end",nx.ME)]),source=TABLE5)
 s7top=g.input("s7_top","Tn5ME-s7_Top",MaterialKind.OLIGO,state=Construct([Segment("Tn5ME-s7_Top",TN5ME_S7_TOP)]),source=TABLE5)
 bot=g.input("mosaic_bot","MOSAIC_Bot (5′ phosphate)",MaterialKind.OLIGO,state=Construct([Segment("MOSAIC_Bot",MOSAIC_BOT)]),source=TABLE5)
 buffer=g.input("anneal_buffer","2× Dialysis Buffer, 95°C 5 min then 5°C steps every 2 min to 20°C",MaterialKind.CONDITION)
 duplexes=g.product("me_duplexes","Poly-T/MOSAIC_Bot and Tn5ME-s7/MOSAIC_Bot duplexes mixed 1:1",MaterialKind.ADAPTER,
     state=molecule("mosaic-end duplex mix",(role("Poly-T overhang duplex"),role("R2N overhang duplex"))))
 g.react("anneal","Anneal top oligos to MOSAIC_Bot and mix complexes 1:1",StepKind.HYBRIDIZATION,
     ((bot,InputRole.TEMPLATE),(polyt,InputRole.STRUCTURAL),(s7top,InputRole.STRUCTURAL),(buffer,InputRole.CONDITION)),(duplexes,))
 tn5=g.input("tn5","Tn5 transposase (Beta Lifescience TN5-BL01) + NaCl, complexed with the duplex mix at 25°C 60 min",MaterialKind.CATALYST)
 reviewed+=["poly_t_top","s7_top","mosaic_bot","anneal_buffer","anneal","me_duplexes","tn5"]
 cells=g.input("cells","FACS-depleted PBMCs",MaterialKind.BIOLOGICAL,state=molecule("ICICLE-seq cells",(role("accessible chromatin"),role("surface epitopes"))),primary=True)
 antibodies=g.input("antibodies","TruStain FcX then 46-antibody TotalSeq-A panel (Supplementary file 6)",MaterialKind.SELECTION,
     inferred=True,source=SOURCE+" (selected TotalSeq-A barcode sequences not printed)")
 stained=g.product("stained","cells retaining antibody-derived tags",MaterialKind.BIOLOGICAL,state=molecule("stained cells",(role("accessible chromatin"),role("antibody tags"))))
 unbound=g.product("unbound","unbound antibodies",MaterialKind.WASTE)
 g.react("stain","Stain on ice 30 min; wash three times in DPBS + 2% BSA",StepKind.BINDING,((cells,InputRole.TEMPLATE),(antibodies,InputRole.SELECTION)),(stained,unbound))
 g.discard(unbound,"removed by three 4 mL washes")
 digitonin=g.input("digitonin","0.01% digitonin permeabilization",MaterialKind.CONDITION)
 perm=g.product("permeabilized","permeabilized stained cells",MaterialKind.BIOLOGICAL,state=stained.state)
 g.react("permeabilize","Permeabilize with 0.01% digitonin",StepKind.PERMEABILIZATION,((stained,InputRole.TEMPLATE),(digitonin,InputRole.CONDITION)),(perm,),changes_molecule=False)
 td1=g.input("td1","TD1 buffer, 20,000 cells, 37°C 120 min",MaterialKind.CONDITION)
 tagged=g.product("tagged","cells with Poly-T/R2N-flanked fragments and antibody tags",MaterialKind.MOLECULAR,
     state=molecule("tagmented cells",(role("Poly-T/R2N fragments (~half; barcodable and indexable)"),
                                       role("Poly-T/Poly-T or R2N/R2N fragments (remain; not barcoded and indexed)"),role("antibody tags"))))
 g.react("tagment","Tagment with 6 µL Poly-T overhang Tn5 complexes",StepKind.TAGMENTATION,((perm,InputRole.TEMPLATE),(duplexes,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST),(td1,InputRole.CONDITION)),(tagged,))
 reviewed+=["cells","antibodies","stain","stained","unbound","digitonin","permeabilize","permeabilized","td1","tagment","tagged"]
 chip=g.input("gems","Chip G, 3′ v3.1 Gel Beads, Partitioning Oil, Chromium Controller",MaterialKind.CONDITION)
 bead=g.input("bead_oligo","10x 3′ v3.1 gel-bead oligo (INFERRED scaffold)",MaterialKind.OLIGO,state=Construct(bead_scaffold()),inferred=True,
     source=SOURCE+" names Chromium Single Cell 3′ v3.1 Gel Beads; oligo sequence not printed")
 fbc=g.input("f_bc","F BC Primer (0.2 µM)",MaterialKind.OLIGO,state=Construct([Segment("F BC Primer",F_BC_PRIMER)]),source=TABLE5)
 revamp=g.input("adt_rev_amp","ADT-Rev-AMP (0.2 µM)",MaterialKind.OLIGO,state=Construct([Segment("ADT-Rev-AMP",ADT_REV_AMP)]),source=TABLE5)
 q5=g.input("q5","NEBNext Ultra II Q5 Master Mix + Reducing Agent B",MaterialKind.CATALYST)
 gems=g.product("gem_contents","partitioned cells with bead oligo and barcoding mix",MaterialKind.MOLECULAR,state=tagged.state)
 empty=g.product("unresolved_gems","empty or multiply occupied GEMs",MaterialKind.WASTE)
 g.react("partition","Partition cells with barcoding master mix and gel beads",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(chip,InputRole.CONDITION)),(gems,empty),changes_molecule=False)
 g.discard(empty,"reads lacking a whitelist-correctable cell barcode are excluded computationally")
 filled=g.product("gap_filled","gap-filled fragments with poly(A) complement at the Poly-T end",MaterialKind.MOLECULAR,state=gems.state)
 g.react("gap_fill","Gap-fill at 72°C 5 min (Fig. 3—S1b)",StepKind.EXTENSION,((gems,InputRole.TEMPLATE),(q5,InputRole.CATALYST)),(filled,))
 barcoded=g.product("barcoded","bead-barcoded ATAC (linear) and ADT (exponential) products",MaterialKind.MOLECULAR,
     state=molecule("barcoded mixture",(role("cell-barcoded ATAC copies"),role("cell-barcoded ADT amplicons"))))
 g.react("barcode_pcr","Denature, anneal bead oligo, extend; 12 cycles 98/42/65°C",StepKind.PCR,
     ((filled,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(fbc,InputRole.STRUCTURAL),(revamp,InputRole.STRUCTURAL),(q5,InputRole.CATALYST)),(barcoded,),
     geometry=PCRGeometry.unavailable(adt_precursor(),reason="One in-GEM reaction: F BC Primer copies ATAC linearly while F BC Primer + ADT-Rev-AMP amplify ADTs exponentially (Fig. 3—S1c); the vendor bead oligo is inferred."))
 reviewed+=["gems","bead_oligo","f_bc","adt_rev_amp","q5","partition","gem_contents","unresolved_gems","gap_fill","gap_filled","barcode_pcr","barcoded"]
 rec,ids=_cleanup(g,"recover","Break GEMs with Recovery Agent; Dynabeads MyOne SILANE cleanup",barcoded,"Recovery Agent and SILANE beads","oil phase and barcoding reagents removed","recovered barcoded products");reviewed+=ids
 sized,ids=_cleanup(g,"size_select","Dual-sided 0.6×/2.0× SPRIselect",rec,"0.6×/2.0× SPRIselect","large DNA fragments and unused primers removed","size-selected barcoded products");reviewed+=ids
 atac_t=g.product("atac_branch","ATAC aliquot (3 parts)",MaterialKind.MOLECULAR,state=atac_precursor())
 adt_t=g.product("adt_branch","ADT aliquot (1 part)",MaterialKind.MOLECULAR,state=adt_precursor())
 g.react("split","Split 3:1 ATAC:ADT",StepKind.SPLIT,((sized,InputRole.TEMPLATE),),(atac_t,adt_t),changes_molecule=False)
 reviewed+=["split","atac_branch","adt_branch"]
 # ATAC indexing: SI-P5-22 + Chromium i7 Multiplex Kit N Set A, Amp Mix, 7 cycles.
 p5=PCRPrimer("SI-P5-22",SI_P5_22,anneal_nt=len(F_BC_PRIMER))
 atac_f=g.input("atac_forward","SI-P5-22 (20 µM)",MaterialKind.OLIGO,state=Construct([Segment("SI-P5-22",SI_P5_22)]),source=TABLE5)
 atac_r=g.input("atac_reverse","Chromium i7 Multiplex Kit N Set A primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (vendor index primer sequence not printed)")
 atac_pol=g.input("atac_polymerase","10x Amp Mix",MaterialKind.CATALYST)
 atac_amp=g.product("atac_amplified","indexed ATAC amplicons",MaterialKind.MOLECULAR,state=ATAC)
 g.react("atac_pcr","ATAC index PCR: 98°C 45 s; 7 × (98/54/72°C)",StepKind.PCR,((atac_t,InputRole.TEMPLATE),(atac_f,InputRole.STRUCTURAL),(atac_r,InputRole.STRUCTURAL),(atac_pol,InputRole.CATALYST)),(atac_amp,),
     geometry=PCRGeometry.unavailable(atac_t.state,reason="SI-P5-22 is printed in Table 5; the Chromium i7 Multiplex Kit N primer is a vendor reagent whose sequence is not printed (Table 5 calls Tn5ME-s7 'compatible with Nextera i7 primers')."))
 atac_lib_state,ids=_cleanup(g,"atac_spri","Dual-sided 0.4×/1.2× SPRIselect",atac_amp,"0.4×/1.2× SPRIselect","outside the retained size window",ATAC.name,MaterialKind.LIBRARY,ATAC)
 reviewed+=["atac_forward","atac_reverse","atac_polymerase","atac_pcr","atac_amplified",*ids]
 g.sequence(atac_lib_state,SequencingPlan(library=ATAC,primers=ATAC_P,acquired_roles=tuple(READ_LENGTHS["ATAC library"]),read_lengths=READ_LENGTHS["ATAC library"]))
 # ADT indexing: SI-P5-22 + ADT-i7, KAPA HiFi, 15 cycles; exact primer landing.
 i7=PCRPrimer("ADT-i7 Primer",il.P7+"N"*8+ADT_I7_ARM,anneal_nt=len(ADT_REV_AMP))
 adt_f=g.input("adt_forward","SI-P5-22 (10 µM)",MaterialKind.OLIGO,state=Construct([Segment("SI-P5-22",SI_P5_22)]),source=TABLE5)
 adt_r=g.input("adt_reverse","ADT-i7 Primer (10 µM; 3′ phosphorothioate)",MaterialKind.OLIGO,
     state=Construct([Segment("P7",il.P7),Segment("i7 index","N"*8,placeholder=True,feature=ADT_I7),Segment("ADT binding arm",ADT_I7_ARM)]),source=TABLE5)
 adt_pol=g.input("adt_polymerase","KAPA HiFi HotStart ReadyMix",MaterialKind.CATALYST)
 adt_amp=g.product("adt_amplified","indexed ADT amplicons",MaterialKind.MOLECULAR,state=ADT)
 g.react("adt_pcr","ADT index PCR: 95°C 3 min; 15 × (95/60/72°C)",StepKind.PCR,((adt_t,InputRole.TEMPLATE),(adt_f,InputRole.STRUCTURAL),(adt_r,InputRole.STRUCTURAL),(adt_pol,InputRole.CATALYST)),(adt_amp,),
     geometry=PCRGeometry.exact(adt_t.state,p5,i7,statement="Table 5 prints SI-P5-22 and ADT-i7; the precursor uses the inferred bead scaffold."))
 adt_lib_state,ids=_cleanup(g,"adt_spri","1.2× SPRIselect",adt_amp,"1.2× SPRIselect","primers and small products removed",ADT.name,MaterialKind.LIBRARY,ADT)
 reviewed+=["adt_forward","adt_reverse","adt_polymerase","adt_pcr","adt_amplified",*ids]
 g.sequence(adt_lib_state,SequencingPlan(library=ADT,primers=ADT_P,acquired_roles=("Read 1",),
     read_lengths={"Read 1":READ_LENGTHS["Antibody-tag library"]["Read 1"]},limitation=ADT_LIMITATION))
 g.declare_reviewed_scope([r for r in dict.fromkeys(reviewed) if r in g.materials or r in g.steps],SOURCE)
 return g.require_publishable()
