"""MUSIC molecular model.

Source scope (round 2026-10-17-breadth-first-13): published Zhao et al. 2024,
PMC11023937 JATS Methods "Critical reagents" (lines 169-203) and the library
procedure from "Ligation of the RNA linker with RNA" to "Sequencing" (lines
221-258); Supplementary Note 2 (MOESM1 text lines 100-200); Supplementary Table 2
(authors' public sheet, export.xlsx.txt). The graph follows the published order:
RNA 5'-phosphorylation and RNA-linker ligation precede chromatin digestion,
dA-tailing and DNA-linker ligation.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct,Scene,Segment,feature,oligo,revcomp
from interaction_capture import FixedMolecularComplex
from pcr import PCRGeometry
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind
from protocol_model import molecule,role

RNA_LINKER="CGAGGAGCGCTTNNNNN"+"ATAGCATTGC"
DNA_TOP="CTAGACACTGTGCGTATCTNB"+"A"*30  # NB + A30 (Methods l.177; Supp. Table 2); was A31 before round 13
DNA_BOTTOM="CGAGGAGNNNNNACAACGCACAGTGTCTAGT"
BC1_TOP="CTCCTCGGCAACCATTACCTT"; BC1_BOTTOM="GTGCTCAAAGGTAATGGTTGC"
BC2_TOP="TGAGCACTCTAGGTGTCGCTG"; BC2_BOTTOM="GCTCGTACAGCGACACCTAGA"
BC3_TOP="TACGAGCCCTTCACTGCAATAAGATCGGAAGA"
BC3_BOTTOM="CAGACGTGTGCTCTTCCGATCTTATTGCAGTGAAGG"
UNIVERSAL="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
INDEX_A1="CAAGCAGAAGACGGCATACGAGATGTCGGTAAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"


def _gem():
 return [Segment("10x GEM complex barcode","G"*16,"cbc",placeholder=True,feature=feature("music_gem_complex","molecular_complex_barcode","whitelist",group="music_complex",part="GEM",whitelist="10x 3-prime whitelist")),
         Segment("10x UMI","U"*12,"umi",placeholder=True,feature=feature("music_complex_umi","umi","random")),
         Segment("poly(dT)VN","T"*30+"VN",placeholder=True)]


def _cell_and_p7():
 return [Segment("cell barcode round 1","A"*14,"cbc",placeholder=True,feature=feature("music_cell_1","cell_barcode","combinatorial",group="music_cell",part="round 1",whitelist="96 barcodes")),
         Segment("cell barcode round 2","B"*14,"cbc",placeholder=True,feature=feature("music_cell_2","cell_barcode","combinatorial",group="music_cell",part="round 2",whitelist="96 barcodes")),
         Segment("cell barcode round 3","C"*14,"cbc",placeholder=True,feature=feature("music_cell_3","cell_barcode","combinatorial",group="music_cell",part="round 3",whitelist="96 barcodes")),
         Segment("Read 2 arm reverse complement",revcomp(il.TRUSEQ_READ2),"r2"),
         Segment("i7 complex barcode reverse complement","I"*8,"cbc",placeholder=True,feature=feature("music_i7_complex","molecular_complex_barcode","whitelist",group="music_complex",part="i7",whitelist="eight indexes per library")),
         Segment("P7 reverse complement",il.P7_RC,"p7")]


_P5=[Segment("P5",il.P5,"p5"),Segment("Read 1 remainder",il.TRUSEQ_READ1[4:],"r1")]
# Supplementary Note 2: cDNA = cell barcode + RNA linker + RNA insert + 10X barcode;
# dsDNA = cell barcode + DNA-linker(bt) + DNA insert + DNA-linker(top) + 10X barcode.
# Linker regions are drawn as role placeholders sized from the printed oligos.
FINAL_RNA=Construct(_P5+_gem()+[
 Segment("RNA insert","X"*48,placeholder=True),
 Segment("RNA-linker RNA portion copy","R"*10,placeholder=True),
 Segment("RNA-linker UMI","N"*5,"umi",placeholder=True,feature=feature("music_rna_linker_umi","umi","random")),
 Segment("RNA-linker DNA arm copy","L"*12,placeholder=True)]+_cell_and_p7(),name="MUSIC RNA-derived library molecule")
FINAL_DNA=Construct(_P5+_gem()+[
 Segment("DNA-linker top-strand copy","L"*19,placeholder=True),
 Segment("DNA insert","X"*48,placeholder=True),
 Segment("DNA-linker bottom duplex region + T","L"*19,placeholder=True),
 Segment("DNA-linker UMI","N"*5,"umi",placeholder=True,feature=feature("music_dna_linker_umi","umi","random")),
 Segment("DNA-linker bottom arm","L"*7,placeholder=True)]+_cell_and_p7(),name="MUSIC DNA-derived library molecule")
PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
for _lib in (FINAL_RNA,FINAL_DNA):
 _e=sp.verify(_lib,PRIMERS,required_roles=("Read 1","Index 1 (i7)","Read 2"))
 if _e:raise ValueError("MUSIC library: "+"; ".join(_e))
FINAL=FINAL_RNA
TITLE="MUSIC — multinucleic-acid interaction mapping in single cells"
NOTES="01_music.html"
SOURCE='<a href="https://doi.org/10.1038/s41586-024-07239-w">Zhao et al., <i>Nature</i> (2024)</a>.'
SUMMARY="RNA and fragmented DNA receive distinct linkers, three rounds of cell barcodes and a two-part molecular-complex barcode. Reads sharing both identities reveal multiway DNA–DNA and RNA–DNA complexes."
CAVEAT=("The public oligo sheet resolves the linkers, representative barcode duplexes and final PCR primers. "
        "The page shows one representative barcode from each 96-member round; the full barcode sets remain in the source sheet. "
        "Linker regions in the final molecules are role placeholders sized from the printed oligos; the barcode-ligation "
        "overhang geometry and the purpose of the Bst 3.0 incubation are not modelled. "
        "Published discrepancy: the paper's Critical reagents DNA-linker paragraph (Methods, extracted PDF lines 1811–1819; "
        "JATS line 177) calls the linker duplex 15 bp, but its DNA-linker ligation paragraph (PDF lines 1991–1995; JATS "
        "line 231) calls the same printed strands 14 bp. The model uses the printed top and bottom oligos from "
        "Supplementary Table 2, whose complementary overlap is 15 bp with a one-base 3′-dT overhang; the source of the "
        "14-bp count is unresolved.")
_INTRO_RNA="Read 1 gives the 16-bp GEM complex barcode and 12-bp UMI; i7 completes complex identity; Read 2 traverses three cell barcodes, the RNA linker and the RNA insert."
_INTRO_DNA="Read 1 gives the 16-bp GEM complex barcode and 12-bp UMI; i7 completes complex identity; Read 2 traverses three cell barcodes, the DNA-linker bottom-strand region and the DNA insert."
FINAL_LIBRARIES=(("RNA-derived library molecule",FINAL_RNA,PRIMERS,"Copied from the poly(A)-tailed RNA by the 10x poly(dT) barcode oligo; the RNA linker sits between insert and cell barcodes.",_INTRO_RNA),
                 ("DNA-derived library molecule",FINAL_DNA,PRIMERS,"Primed by the 10x poly(dT) on the DNA-linker top-strand poly(A); the DNA linker appears on both sides of the insert (Supplementary Note 2).",_INTRO_DNA))
_READS={"Read 1":28,"Index 1 (i7)":8,"Read 2":150}
READ_LENGTHS={h:_READS for h,*_ in FINAL_LIBRARIES}
COMPLEX=FixedMolecularComplex(("RNA molecule","DNA fragment A","DNA fragment B"),"DSG + 3% formaldehyde",label="multiway nucleic-acid complex")


def oligos():
 yield oligo("RNA linker",[Segment("DNA portion + 5-nt UMI",RNA_LINKER[:17],"umi",placeholder=True),Segment("RNA portion",RNA_LINKER[17:])],mods="/5OH/",three="/3OH/-3′ (final 10 bases RNA)")
 yield oligo("DNA linker top",[Segment("phosphorylated top",DNA_TOP,placeholder=True)],mods="/5Phos/",three="/3OH/-3′")
 yield oligo("DNA linker bottom",[Segment("bottom with 5-nt UMI",DNA_BOTTOM,"umi",placeholder=True)],mods="/5OH/",three="/3OH/-3′")
 yield oligo("barcode 1 top, member 1",[Segment("BC1 top",BC1_TOP)],mods="/5Phos/")
 yield oligo("barcode 1 bottom, member 1",[Segment("BC1 bottom",BC1_BOTTOM)],mods="/5Phos/")
 yield oligo("barcode 2 top, member 1",[Segment("BC2 top",BC2_TOP)],mods="/5Phos/")
 yield oligo("barcode 2 bottom, member 1",[Segment("BC2 bottom",BC2_BOTTOM)],mods="/5Phos/")
 yield oligo("barcode 3 top, member 1",[Segment("BC3 top",BC3_TOP)],mods="/5Phos/")
 yield oligo("barcode 3 bottom, member 1",[Segment("BC3 bottom",BC3_BOTTOM)],mods="/5OH/",three="/3OH/-3′")
 yield oligo("universal PCR adapter",[Segment("universal PCR adapter",UNIVERSAL[:-1]),Segment("terminal phosphorothioate T","*T",placeholder=True)],mods="/5OH/",three="/3OH/-3′")
 yield oligo("index PCR adapter A1",[Segment("index PCR adapter A1",INDEX_A1[:-1]),Segment("terminal phosphorothioate T","*T",placeholder=True)],mods="/5OH/",three="/3OH/-3′")


def dna_linker_geometry():
 sc=Scene(); sc.strand("linker",[Segment("Y-linker body","[published Y-linker duplex]",placeholder=True),Segment("single-base dT","T")],label="DNA linker")
 sc.ligation_end("linker","single-base dT",kind="cohesive",polarity="3-prime",owner="adapter",text="published 3′-dT sticky end")
 sub=Scene(); sub.strand("DNA",[Segment("dA-tailed genomic DNA","X"*24,placeholder=True),Segment("genomic 3-prime dA","A")],label="DNA")
 sub.ligation_end("DNA","genomic 3-prime dA",kind="cohesive",polarity="3-prime",owner="substrate",text="3′-dA")
 return AdapterLigationGeometry.known("cohesive",sc,sub,statement="The published MUSIC Y linker has a single-base sticky end complementary to dA-tailed genomic DNA; its full top and bottom oligos remain copyable above.")


def rna_linker_geometry():
 """Methods l.171/223: 5'-OH DNA / 3'-OH RNA chimeric linker; RNA 5' ends are first phosphorylated by T4 PNK."""
 sc=Scene(); sc.strand("RNA linker",[Segment("DNA arm + 5-nt UMI",RNA_LINKER[:17],"umi",placeholder=True),Segment("RNA portion (rArUrArGrCrArUrUrGrC)",RNA_LINKER[17:])],label="chimeric RNA linker",mod5="OH",mod3="OH")
 sc.ligation_end("RNA linker","RNA portion (rArUrArGrCrArUrUrGrC)",kind="single-strand",polarity="3-prime",owner="adapter",text="linker RNA 3′-OH")
 sub=Scene(); sub.strand("RNA",[Segment("nuclear RNA","R"*24,placeholder=True)],label="5′-phosphorylated nuclear RNA",mod5="P")
 sub.ligation_end("RNA","nuclear RNA",kind="single-strand",polarity="5-prime",owner="substrate",text="RNA 5′-phosphate (T4 PNK)")
 return AdapterLigationGeometry.known("single-strand",sc,sub,statement="T4 RNA ligase 1 joins the linker's 3′-OH RNA end to the PNK-phosphorylated 5′ end of nuclear RNA; the linker's 5′-OH DNA end cannot act as a donor until the later PNK step.")


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE+" PMC11023937 Methods; Supplementary Note 2; Supplementary Table 2")
 reviewed=[]
 def inp(i,label,kind,**kw):
  reviewed.append(i); return g.input(i,label,kind,**kw)
 def step(i,label,kind,ins,outs,**kw):
  reviewed.append(i); return g.react(i,label,kind,ins,outs,**kw)
 def nuc(i,label,state=None,**kw):
  return g.product(i,label,MaterialKind.MOLECULAR,state=state if state is not None else complex_state,**kw)
 def wash(i,label,reason):
  w=g.product(i,label,MaterialKind.WASTE); return w,reason
 complex_state=COMPLEX.complex()
 nuclei=g.input("nuclei","DSG + 3% formaldehyde crosslinked nuclei, lysed and SDS/Triton-permeabilized in rCutSmart",MaterialKind.BIOLOGICAL,state=complex_state,primary=True)
 reviewed.append("nuclei")

 pnk=inp("pnk_rna","T4 PNK + 1 mM ATP, 37 °C 1 h (5′ phosphorylation master mix)",MaterialKind.CATALYST)
 rna5p=nuc("rna_5p","nuclei with 5′-phosphorylated RNA")
 step("phosphorylate_rna","Phosphorylate RNA 5′ ends",StepKind.PHOSPHORYLATION,((nuclei,InputRole.TEMPLATE),(pnk,InputRole.CATALYST)),(rna5p,))
 rna_linker=inp("rna_linker","chimeric RNA linker, 5′-OH DNA / 3′-OH RNA, 4 µM",MaterialKind.ADAPTER,state=rna_linker_geometry().scenes[0])
 rnl1=inp("t4_rnl1","T4 RNA ligase 1, 15% PEG 8000, 1 mM ATP; 25 °C 2 h then 16 °C overnight",MaterialKind.CATALYST)
 rna_linked=nuc("rna_linked","nuclei with RNA-linker–RNA junctions")
 rl_waste,rl_reason=wash("rna_linker_wash","unligated RNA linker","removed by PBS wash buffer 1 and three PBS wash buffer 2 washes")
 step("ligate_rna_linker","Ligate the RNA linker to nuclear RNA 5′ ends",StepKind.LIGATION,((rna5p,InputRole.TEMPLATE),(rna_linker,InputRole.STRUCTURAL),(rnl1,InputRole.CATALYST)),(rna_linked,rl_waste),geometry=rna_linker_geometry())
 g.discard(rl_waste,rl_reason)

 hpy=inp("hpych4v","HpyCH4V in rCutSmart with RNase inhibitor, 37 °C 3 h",MaterialKind.CATALYST)
 digested=nuc("digested","nuclei with HpyCH4V-digested chromatin")
 dig_waste,dig_reason=wash("digestion_supernatant","digestion supernatant","removed with the supernatant and four washes")
 step("digest","Digest chromatin in nuclei",StepKind.DIGESTION,((rna_linked,InputRole.TEMPLATE),(hpy,InputRole.CATALYST)),(digested,dig_waste))
 g.discard(dig_waste,dig_reason)
 klenow=inp("klenow","Klenow fragment in NEBNext dA-tailing buffer, 37 °C 1.5 h",MaterialKind.CATALYST)
 tailed=nuc("da_tailed","nuclei with dA-tailed DNA fragments")
 step("da_tail","dA-tail DNA fragment 3′ ends",StepKind.TAILING,((digested,InputRole.TEMPLATE),(klenow,InputRole.CATALYST)),(tailed,))
 dna_linker=inp("dna_linker","published UMI-bearing Y DNA linker, 4.5 µM",MaterialKind.ADAPTER,state=dna_linker_geometry().scenes[0])
 sticky=inp("sticky_ligase","Instant Sticky-end Ligase Master Mix with Quick Ligase Buffer and 6% 1,2-propanediol, 20 °C 6 h",MaterialKind.CATALYST)
 dna_linked=nuc("dna_linked","nuclei with RNA and DNA linkers attached")
 step("ligate_dna_linker","Ligate the Y linker to dA-tailed genomic DNA",StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(dna_linker,InputRole.STRUCTURAL),(sticky,InputRole.CATALYST)),(dna_linked,),geometry=dna_linker_geometry())
 pnk2=inp("pnk_linker","T4 PNK + ATP, 37 °C 1 h",MaterialKind.CATALYST)
 ready=nuc("linker_5p","nuclei with 5′-phosphorylated linker arms; strained and counted")
 step("phosphorylate_linkers","Phosphorylate the linker 5′ ends for barcode ligation",StepKind.PHOSPHORYLATION,((dna_linked,InputRole.TEMPLATE),(pnk2,InputRole.CATALYST)),(ready,))

 current=ready
 for n in (1,2,3):
  adapter=inp(f"cell_barcode_{n}",f"barcode set {n}: 96 duplexes with 7-nt top overhang, 14-bp unique duplex, {'11' if n==3 else '7'}-nt bottom overhang",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (strand-resolved overhang pairing not transcribed in this model)")
  lig=inp(f"barcode_ligase_{n}",f"ligation master mix {n}: 5.4 µM barcodes, Instant Sticky-end Ligase, 1,2-propanediol; overnight 20 °C; EDTA/EGTA quench",MaterialKind.CATALYST)
  out=nuc(f"cell_indexed_{n}",f"pooled nuclei after barcode round {n}",inferred=True,source=SOURCE+" (barcode identity known; overhang geometry not transcribed)")
  w,r=wash(f"barcode_wash_{n}",f"unligated round-{n} barcodes","removed by quench, pooling and three PBS wash buffer 2 washes")
  step(f"ligate_cell_barcode_{n}",f"Split nuclei into 96 wells, ligate barcode set {n}, quench and pool",StepKind.LIGATION,((current,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(lig,InputRole.CATALYST)),(out,w),geometry=AdapterLigationGeometry.unavailable(f"Set {n} overhang complementarity is published (Extended Data Fig. 2d,i; Supplementary Table 2) but its strand-resolved junction is not drawn in this model."))
  g.discard(w,r); current=out

 dephos=inp("pnk_dephos","T4 PNK 0.5 U/µl in PNK buffer (no ATP listed), 37 °C 1 h",MaterialKind.CATALYST)
 deph=nuc("rna_3oh","nuclei with RNA 2′,3′-cyclic phosphates converted to 3′-OH")
 step("dephosphorylate","Convert RNA 3′ cyclic phosphates to 3′-OH",StepKind.DEPHOSPHORYLATION,((current,InputRole.TEMPLATE),(dephos,InputRole.CATALYST)),(deph,))
 pap=inp("polya_polymerase","E. coli poly(A) polymerase + 1 mM ATP, 37 °C 10 min",MaterialKind.CATALYST)
 polya=nuc("polya","nuclei with poly(A)-tailed RNA")
 step("polya_tail","Poly(A)-tail all RNA 3′ ends",StepKind.TAILING,((deph,InputRole.TEMPLATE),(pap,InputRole.CATALYST)),(polya,))
 count=inp("nuclei_count","10-µm strainer, ethidium homodimer count; 5,000 nuclei taken",MaterialKind.CONDITION)
 sampled=nuc("sampled_nuclei","5,000 barcoded nuclei")
 surplus,s_reason=wash("surplus_nuclei","nuclei not taken forward","not used for the 5,000-nucleus GEM loading")
 step("sample_nuclei","Strain, count and take 5,000 nuclei",StepKind.SAMPLING,((polya,InputRole.TEMPLATE),(count,InputRole.CONDITION)),(sampled,surplus))
 g.discard(surplus,s_reason)
 covaris=inp("sonication","Covaris M220, 50 W, duty factor 5, 6 °C, 5 min in 15 µl",MaterialKind.CONDITION)
 released=nuc("released_complexes","released barcoded molecular complexes")
 step("sonicate","Sonicate nuclei to release molecular complexes",StepKind.FRAGMENTATION,((sampled,InputRole.TEMPLATE),(covaris,InputRole.CONDITION)),(released,))
 chromium=inp("gem_condition","Chromium Next GEM Single Cell 3′ kit (PN-1000268), steps 1.1–1.5 with RT reagent B, reducing agent B, RT enzyme C",MaterialKind.CONDITION)
 partitioned=nuc("partitioned","molecular complexes in GEMs")
 gem_loss,gem_reason=wash("gem_unrecovered","complexes not encapsulated with a gel bead","not barcoded in a GEM")
 step("partition","Partition complexes into 10x GEMs",StepKind.PARTITION,((released,InputRole.TEMPLATE),(chromium,InputRole.CONDITION)),(partitioned,gem_loss),changes_molecule=False)
 g.discard(gem_loss,gem_reason)
 gem_oligo=inp("gem_oligo","10x barcode oligo: partial Read 1, 16-nt GEM barcode, 12-nt UMI, dT30VN",MaterialKind.OLIGO)
 rt=inp("rt_enzyme","10x RT enzyme C",MaterialKind.CATALYST)
 first=nuc("first_strand","10x-barcoded copies: cDNA of poly(A) RNA and linker-poly(A)-primed DNA copies",state=molecule("10x-barcoded first strand",tuple(_gem())+(role("RNA or DNA copy with linker and cell barcodes"),)))
 step("gem_rt","Prime poly(A) RNA and the DNA-linker top-strand poly(A) with 10x poly(dT) oligos; copy in GEMs",StepKind.REVERSE_TRANSCRIPTION,((partitioned,InputRole.TEMPLATE),(gem_oligo,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
 revx=inp("reverse_crosslink_buffer","break emulsion; eight aliquots; NaCl/SDS/EDTA/EGTA + proteinase K, 50 °C 2 h then 55 °C overnight",MaterialKind.CONDITION)
 decross=nuc("decrosslinked","reverse-crosslinked nucleic acids",state=first.state)
 step("reverse_crosslink","Reverse crosslinks",StepKind.DEPROTEINIZATION,((first,InputRole.TEMPLATE),(revx,InputRole.CONDITION)),(decross,))
 monarch=inp("monarch_rna","Monarch RNA purification column, elute 21 µl",MaterialKind.SELECTION)
 purified=nuc("purified","purified RNA and DNA",state=first.state)
 m_waste,m_reason=wash("column_flowthrough","column flow-through","not retained by the column")
 step("column","Purify nucleic acids",StepKind.CLEANUP,((decross,InputRole.TEMPLATE),(monarch,InputRole.SELECTION)),(purified,m_waste))
 g.discard(m_waste,m_reason)
 bst=inp("bst3","Bst 3.0 DNA polymerase, isothermal amplification buffer II, extra 6 mM MgSO4, dNTPs, 55 °C 15 min",MaterialKind.CATALYST)
 extended=nuc("bst_extended","Bst 3.0-treated products",state=first.state)
 step("bst","Bst 3.0 incubation (purpose not stated in the inspected Methods)",StepKind.EXTENSION,((purified,InputRole.TEMPLATE),(bst,InputRole.CATALYST)),(extended,))
 rnaclean=inp("rnaclean","1.8× RNAClean XP beads, elute 20 µl",MaterialKind.SELECTION)
 cleaned=nuc("bst_cleaned","bead-purified PCR templates",state=first.state)
 b_waste,b_reason=wash("rnaclean_supernatant","bead supernatant","not retained by 1.8× beads")
 step("rnaclean_cleanup","Bead-purify",StepKind.CLEANUP,((extended,InputRole.TEMPLATE),(rnaclean,InputRole.SELECTION)),(cleaned,b_waste))
 g.discard(b_waste,b_reason)

 universal=inp("universal_adaptor","Universal Adaptor (P5 + Read 1), 3′ phosphorothioate",MaterialKind.OLIGO,state=UNIVERSAL)
 index=inp("index_adaptor","index adaptor, one of eight per aliquot (example A1), 3′ phosphorothioate",MaterialKind.OLIGO,state=INDEX_A1)
 q5=inp("q5","NEBNext Ultra II Q5 Master Mix, 13–14 cycles",MaterialKind.CATALYST)
 rna_amp=g.product("rna_amplicons","RNA-derived amplicons",MaterialKind.MOLECULAR,state=FINAL_RNA)
 dna_amp=g.product("dna_amplicons","DNA-derived amplicons",MaterialKind.MOLECULAR,state=FINAL_DNA)
 step("pcr","Amplify each aliquot with the Universal Adaptor and its index adaptor",StepKind.PCR,((cleaned,InputRole.TEMPLATE),(universal,InputRole.STRUCTURAL),(index,InputRole.STRUCTURAL),(q5,InputRole.CATALYST)),(rna_amp,dna_amp),geometry=PCRGeometry.unavailable(FINAL_RNA,reason="Universal Adaptor hybridizes the 10x partial Read 1 and the index adaptor the Read 2 part of barcode set 3 (Methods l.187, 197–203); the pre-PCR 10x Read 1 bases are in the unread CG000315 guide, so the primer sites are not drawn on a pre-PCR template."))
 ampure=inp("ampure","1.2× AMPure beads per aliquot; eight aliquots combined",MaterialKind.SELECTION)
 rna_c=g.product("rna_clean","purified RNA-derived amplicons",MaterialKind.MOLECULAR,state=FINAL_RNA)
 dna_c=g.product("dna_clean","purified DNA-derived amplicons",MaterialKind.MOLECULAR,state=FINAL_DNA)
 a_waste,a_reason=wash("ampure_supernatant","AMPure supernatant","not retained by 1.2× beads")
 step("ampure_cleanup","Bead-purify and pool the eight index libraries",StepKind.CLEANUP,((rna_amp,InputRole.TEMPLATE),(dna_amp,InputRole.TEMPLATE),(ampure,InputRole.SELECTION)),(rna_c,dna_c,a_waste))
 g.discard(a_waste,a_reason)
 gel=inp("gel","4% E-gel, excise 300 bp–1.2 kb, Monarch gel extraction",MaterialKind.SELECTION)
 rna_lib=g.product("rna_library","RNA-derived MUSIC library molecules",MaterialKind.LIBRARY,state=FINAL_RNA)
 dna_lib=g.product("dna_library","DNA-derived MUSIC library molecules",MaterialKind.LIBRARY,state=FINAL_DNA)
 g_waste,g_reason=wash("gel_rejected","fragments outside 300 bp–1.2 kb","left in the gel")
 step("gel_select","Size-select by gel",StepKind.SELECTION,((rna_c,InputRole.TEMPLATE),(dna_c,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(rna_lib,dna_lib,g_waste))
 g.discard(g_waste,g_reason)
 for lib,final in ((rna_lib,FINAL_RNA),(dna_lib,FINAL_DNA)):
  g.sequence(lib,SequencingPlan(library=final,primers=PRIMERS,acquired_roles=tuple(p.role for p in PRIMERS),read_lengths=_READS))
  reviewed.append(lib.id)
 g.declare_reviewed_scope(reviewed,g.authority)
 return g.require_publishable()
