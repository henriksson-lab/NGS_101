"""Cao et al. sci-CAR RNA/ATAC co-assay."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
seg=Segment
import illumina as il
import nextera as nx
import seqprimers as sp
from protocol_graph import Disposition, Terminal
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from pcr import PCRGeometry

TITLE="sci-CAR-seq — joint combinatorial RNA and ATAC profiling"
NOTES="01_sci-car-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1126/science.aau0730">Cao et al., <i>Science</i> (2018)</a>.'
SUMMARY="One intact nucleus receives a barcoded poly(dT) RT primer and a separately barcoded Tn5 adapter in the first plate; FACS redistribution supplies a second well index while RNA and ATAC aliquots become separate libraries."
CAVEAT=("<b>Source scope and inferred adapter sequences.</b> Step order, reagents, conditions, the RT primer (Read 1 partial arm, 8-N UMI, 10-bp RT index, T30VN), the RNA and ATAC P5/P7 PCR primers, the cycle numbers and the RNA read recipe come from the published Supplementary Materials (Materials and Methods, &#8216;Sample processing for cultured cells&#8217;, pp. 16&#8211;18); the article&#8217;s Fig. 1 legend assigns each identifier to a read. The supplementary table workbook (Tables S1&#8211;S13) holds analysis results, not oligos. Dotted segments remain inferred: the mosaic end of the &#8216;i7 only TDE1&#8217; adaptor on the RNA cDNA end (the P7 primer establishes only S7), the interior of the &#8216;indexed TDE1&#8217; ATAC transposons, which follows the Amini 2014 Table 4 architecture (connector, 8-nt index, extended mosaic end) because the article attributes transposase-complex preparation to that reference; the ATAC P5/P7 primers establish only the S5/S7 ends. ATAC libraries were read 51 + 51 with a custom recipe and &#8216;custom primers as previously described&#8217; (Cusanovich 2018, not read here); the ATAC sequencing primers drawn are therefore inferred and ATAC index-read cycles are not stated. "
        "RNA and ATAC first-round indexes are chemically distinct but map to the same first well; their second-round PCR indexes (carried by both PCR primers of each branch) likewise map to the same FACS well. "
        "<b>Variants not drawn:</b> PFA-fixed mouse-kidney nuclei are reverse-crosslinked (SDS, proteinase K, 65 &#176;C 16 h) after sorting and split before second-strand synthesis, which is then done on the AMPure-purified RNA portion only; &#8216;sci-RNA-seq only&#8217; and &#8216;sci-ATAC-seq only&#8217; control branches. "
        "<b>Published predecessor correction:</b> <a href=\"https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fng.3119/MediaObjects/41588_2014_BFng3119_MOESM2_ESM.xlsx\">Amini 2014 Supplementary Table 4</a>, Correction Notice and corrected rows 1–5, says the first five indexed Tn5 oligo rows were omitted from the original online table and restored on 20 May 2015. Their corrected 8-base i5 values are TATAGCCT, ATAGAGGC, CCTATCCT, GGCTCTGA and AGGCGAAG. This model uses generic index placeholders and does not assign those Amini members to sci-CAR; the inspected Cao 2018 supplementary materials print no barcode lists.")
RNA1=feature("rna_rt_well","cell_barcode","combinatorial",group="cell_id",part="RT well",whitelist="published sci-CAR RT primer plate")
ATAC1=feature("atac_tn5_well","cell_barcode","combinatorial",group="cell_id",part="Tn5 well (P5 side)",whitelist="published sci-CAR Tn5 plate")
UMI=feature("rna_umi","umi","random")
# Fig. 1 legend (PMC6571013 JATS line 475): RNA -- index 2 and read 1 cover the i5
# index, UMI and RT barcode; index 1 and read 2 cover the i7 index and cDNA.
# ATAC -- reads 1 and 2 cover genomic DNA; indices 1 and 2 cover the Tn5 and PCR
# barcodes.  Line 78: both PCR primers of each branch carry the second (PCR-well)
# index.  The published Supplementary Materials (NIHMS1033869-supplement-1.pdf,
# 'Sample processing for cultured cells', pp. 16-18; text twin lines 62-243) print
# the same primers, enzymes and cycles as the earlier-read GEO text, plus the steps
# GEO truncated (RNA/ATAC cleanups, ATAC sequencing recipe).  GEO GSE117089
# sample extraction protocol (author-deposited; identical in GSM3271040-45,
# geo-GSE117089-all-full.txt line 344488) prints the RT primer
# ACGACGCTCTTCCGATCT-N8-[10bp index]-T30VN, RNA P5 primer
# P5-[i5]-ACACTCTTTCCCTACACGACGCTCTTCCGATCT, ATAC P5 primer P5-[i5]-TCGTCGGCAGCGTC,
# P7 primer CAAGCAGAAGACGGCATACGAGAT-[i7]-GTCTCGTGGGCTCGG (both branches),
# 'i7 only TDE1' for cDNA tagmentation and the cycle numbers below.
GEO="GEO GSE117089 sample extraction protocol (author-deposited; GSM3271040-GSM3271045)"
SM="Cao 2018 Supplementary Materials, Materials and Methods, 'Sample processing for cultured cells' (NIHMS1033869-supplement-1.pdf pp. 16-18)"
AMINI = "Amini 2014 Supplementary Table 4 (doi:10.1038/ng.3119), cited by Cao 2018 line 140 for transposase-complex preparation"
# Amini Table 4 transposon architecture: S5 + connector A + 8-nt index + extension A + ME;
# S7 + connector B + 8-nt index + extension B + ME.
CONNECTOR_A="TCCACGC"; EXT_A="GCGATCGAGGACGGC"
CONNECTOR_B="CTGTCCCTGTCC"; EXT_B="CACCGTCTCCGCCTC"
ATAC_R1=EXT_A+nx.ME
ATAC_R2=EXT_B+nx.ME
ATAC_PRIMERS=(
 sp.custom("Read 1","extended mosaic end A (inferred)",ATAC_R1,AMINI),
 sp.custom("Index 1 (i7)","reverse complement of extended mosaic end B (inferred)",revcomp(ATAC_R2),AMINI),
 sp.custom("Index 2 (i5)","reverse complement of extended mosaic end A (derived, inferred)",revcomp(ATAC_R1),AMINI+"; derived as the reverse complement of Read 1"),
 sp.custom("Read 2","extended mosaic end B (inferred)",ATAC_R2,AMINI),
)

def _pcr_index(name,feat):
 return seg(name,"[PCR index; length unassigned]","cbc",placeholder=True,feature=feat)

def rna_library():
 i5=feature("rna_pcr_i5","cell_barcode","combinatorial",group="cell_id",part="FACS/PCR well (i5)",whitelist="published sci-CAR PCR plate")
 i7=feature("rna_pcr_i7","cell_barcode","combinatorial",group="cell_id",part="FACS/PCR well (i7)",whitelist="published sci-CAR PCR plate")
 lib=Construct([seg("P5",il.P5,"p5"),_pcr_index("RNA PCR i5",i5),
  seg("TruSeq Read 1 arm (P5 primer tail + RT primer 18 nt)",il.TRUSEQ_READ1,"r1"),
  seg("UMI","U"*8,"umi",placeholder=True,feature=UMI),
  seg("RT well index","A"*10,"cbc",placeholder=True,feature=RNA1),
  seg("oligo(dT) and 3-prime cDNA tag","X"*38,placeholder=True),
  seg("unindexed Tn5 mosaic end",nx.ME_RC,"me",inferred=True),
  seg("Tn5 adaptor S7 (P7 primer 3' end)",nx.S7_RC,"s7"),
  Segment("RNA PCR i7'","[PCR i7; length unassigned]","cbc",True,feature=i7),
  seg("P7'",il.P7_RC,"p7")],name="sci-CAR RNA library")
 primers=(sp.TRUSEQ["R1"],sp.NEXTERA["I1"],sp.TRUSEQ["I2"],sp.NEXTERA["R2"])
 errors=sp.verify(lib,primers)
 if errors: raise ValueError("invalid sci-CAR RNA library: "+"; ".join(errors))
 return lib,primers

def atac_library():
 i5=feature("atac_pcr_i5","cell_barcode","combinatorial",group="cell_id",part="FACS/PCR well (i5)",whitelist="published sci-CAR PCR plate")
 i7=feature("atac_pcr_i7","cell_barcode","combinatorial",group="cell_id",part="FACS/PCR well (i7)",whitelist="published sci-CAR PCR plate")
 tn5_b=feature("atac_tn5_well_b","cell_barcode","combinatorial",group="cell_id",part="Tn5 well (P7 side)",whitelist="published sci-CAR Tn5 plate")
 lib=Construct([seg("P5",il.P5,"p5"),_pcr_index("ATAC PCR i5",i5),
  seg("S5 (ATAC P5 primer 3' end)",nx.S5,"s5"),seg("connector A",CONNECTOR_A,"r2",inferred=True),
  seg("Tn5 well index (P5 side)","A"*8,"cbc",placeholder=True,inferred=True,feature=ATAC1),
  seg("extension A",EXT_A,"r1",inferred=True),seg("ME",nx.ME,"me",inferred=True),
  seg("accessible DNA","X"*38,placeholder=True),
  seg("ME'",nx.ME_RC,"me",inferred=True),seg("extension B'",revcomp(EXT_B),"r1",inferred=True),
  Segment("Tn5 well index (P7 side)'","B"*8,"cbc",True,inferred=True,feature=tn5_b),
  seg("connector B'",revcomp(CONNECTOR_B),"r3",inferred=True),seg("S7' (P7 primer 3' end)",nx.S7_RC,"s7"),
  Segment("ATAC PCR i7'","[PCR i7; length unassigned]","cbc",True,feature=i7),
  seg("P7'",il.P7_RC,"p7")],name="sci-CAR ATAC library")
 errors=sp.verify(lib,ATAC_PRIMERS)
 if errors: raise ValueError("invalid sci-CAR ATAC library: "+"; ".join(errors))
 return lib,ATAC_PRIMERS

RNA,RNA_P=rna_library()
ATAC,ATAC_P=atac_library()
FINAL_LIBRARIES=(("RNA library",RNA,RNA_P,"INFERRED — the cDNA-end mosaic end (vendor 'i7 only TDE1' adaptor) ; physical PCR-index lengths remain unassigned; P5/P7 primers, Read 1 arm, 8-N UMI and 10-bp RT index follow the published Supplementary Materials (pp. 16-17).","Fig. 1 legend: Read 1 covers the UMI and RT barcode, Index 2 the i5 PCR index, Index 1 the i7 PCR index and Read 2 the cDNA."),("ATAC library",ATAC,ATAC_P,"INFERRED — transposon interior (connector, Tn5 index, extended mosaic end) follows Amini 2014 Table 4 (the cited transposase-complex method); the supplementary materials call it 'indexed TDE1 enzyme (Illumina)' without bases. S5/S7 ends follow the printed PCR primers (Supplementary Materials p. 18).","Fig. 1 legend: Reads 1 and 2 cover genomic DNA; Index 1 and Index 2 each cover a Tn5 barcode and a PCR barcode."))
# Supplementary Materials p. 17 (and GEO GSE117089): RNA "Read 1: 18 cycles, Read 2: 52 cycles, Index 1: 10
# cycles, Index 2: 10 cycles" (NextSeq V2 75-cycle kit); ATAC p. 18 "paired end 51 bp
# reads with index reads that covered both the Tn5 barcode and the library amplification
# barcode", custom recipe and primers from ref 29 (index cycles not stated).
READ_LENGTHS={"RNA library":{"Read 1":18,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":52},"ATAC library":{"Read 1":51,"Read 2":51}}

def _rt_parts():
 return [seg("RT handle","ACGACGCTCTTCCGATCT"),seg("UMI","N"*8,"umi",placeholder=True,feature=UMI),seg("RT-well barcode","N"*10,"cbc",placeholder=True,feature=RNA1),seg("oligo(dT)VN","T"*30+"VN",placeholder=True),seg("native cDNA extent unassigned","[cDNA]",placeholder=True)]

def population(*, copied=False, tagged=False, ds=False, rna_tagged=False):
 sc=Scene()
 dna=[seg("genomic DNA","[chromatin DNA]",placeholder=True)]
 if tagged:
  dna=[seg("first-well Tn5 barcode; native adapter length unassigned","[Tn5 barcode]","cbc",placeholder=True,feature=ATAC1),*dna]
 sc.strand("dna",dna,label="independent ATAC population; native strand ends/gaps unassigned")
 sc.strand("rna",[seg("native RNA; fate through lysis/cleanup unassigned","[RNA]",placeholder=True)],label="RNA population")
 if copied:
  parts=_rt_parts()
  if rna_tagged: parts += [seg("cDNA i7-only TDE1 end; position and gap unassigned","[TDE1 end]",placeholder=True)]
  sc.strand("cdna",parts,label="cDNA RT-primer strand")
  if ds:
   sc.anneal("second",complement_segments(parts),"cdna",pair=(parts[0].name+"'",parts[0].name),label="second strand; native terminal extent unassigned",unpaired=tuple(x.name+"'" for x in parts if x.is_role_token()))
 return sc

SHOW_READ_CONTENT=False
READ_SPECS={
 "RNA library":(sp.read_spec("Read 1",18,"N8 UMI + RT-well N10 barcode","extract"),sp.read_spec("Index 1 (i7)",10,"PCR-well i7; physical length unassigned","extract"),sp.read_spec("Index 2 (i5)",10,"PCR-well i5; physical length unassigned","extract"),sp.read_spec("Read 2",52,"RNA-derived cDNA","include")),
 "ATAC library":(sp.read_spec("Read 1",51,"genomic DNA; primer geometry conditional","include"),sp.read_spec("Index 1 (i7)","unassigned","Tn5 + PCR-well barcodes; custom recipe unassigned","undetermined"),sp.read_spec("Index 2 (i5)","unassigned","Tn5 + PCR-well barcodes; custom recipe unassigned","undetermined"),sp.read_spec("Read 2",51,"genomic DNA; primer geometry conditional","include"))}
READ_SPECS_NOTE="RNA cycle counts are printed in the published Supplementary Materials. ATAC index acquisition is positive but cycle counts/custom recipe are unassigned; neither run implies physical PCR-index widths."


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 nucleus=g.input("nucleus","intact nucleus containing accessible chromatin and poly(A) RNA",MaterialKind.BIOLOGICAL,state=population(),primary=True)
 rtprimer=g.input("rt_primer","well-indexed anchored oligo(dT) RT primer (Read 1 partial arm, 8-N UMI, 10-bp index, T30VN)",MaterialKind.OLIGO,source=SM,
  state=Construct([seg("TruSeq Read 1 partial arm","ACGACGCTCTTCCGATCT","r1"),seg("UMI","N"*8,"umi",placeholder=True,feature=UMI),
   seg("RT well index","A"*10,"cbc",placeholder=True,feature=RNA1),seg("oligo(dT)30","T"*30),seg("VN anchor","VN",placeholder=True)],name="sci-CAR RT primer"))
 rt=g.input("rt","SuperScript IV reverse transcriptase with RNaseOUT (primer/dNTP 55 °C 5 min, then 4→55 °C gradient)",MaterialKind.CATALYST,source=SM)
 copied=g.product("copied","nucleus retaining first-well-indexed cDNA and chromatin",MaterialKind.MOLECULAR,state=population(copied=True))
 g.react("reverse_transcribe","Index RNA by in-situ reverse transcription",StepKind.REVERSE_TRANSCRIPTION,((nucleus,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))
 transposome=g.input("atac_transposome","well-indexed TDE1 (Illumina) in Nextera TD buffer, 55 °C 30 min, stopped with EDTA/spermidine; transposon bases not printed (Amini 2014 architecture, inferred)",MaterialKind.ADAPTER,inferred=True,source=AMINI)
 tn5=g.input("atac_tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=g.product("tagged","nucleus carrying indexed cDNA and tagmented chromatin",MaterialKind.MOLECULAR,state=population(copied=True,tagged=True))
 twaste=g.product("tn5_waste","candidate non-library tagmentation products",MaterialKind.MOLECULAR)
 g.react("tagment_atac","Index accessible DNA by in-situ tagmentation",StepKind.TAGMENTATION,((copied,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,twaste)); g.mark_terminal(Terminal(twaste.id, Disposition.UNRESOLVED, 'Candidate non-library end classes remain with the nuclei and enter ATAC PCR. They lack the intended P5/P7 end combination; amplification and physical recovery are unresolved.'))
 facs=g.input("facs_wells","pool, DAPI stain (3 µM), singlet-gated FACS of 25 nuclei per well into 5 µL EB",MaterialKind.CONDITION,source=SM)
 redistributed=g.product("redistributed","uniquely redistributed indexed nuclei",MaterialKind.MOLECULAR,state=tagged.state)
 excluded=g.product("excluded","DAPI-gated doublets and unsorted nuclei",MaterialKind.WASTE)
 g.react("redistribute","Pool nuclei and redistribute them by FACS",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(facs,InputRole.CONDITION)),(redistributed,excluded),changes_molecule=False); g.discard(excluded,"physically excluded by the sorter gate; barcode collisions between separately sorted nuclei remain in the wells and can only be recognised computationally")
 second_pol=g.input("second_polymerase","mRNA Second Strand Synthesis buffer and enzyme (NEB), 16 °C 180 min; no added primer is printed",MaterialKind.CATALYST,source=SM)
 ds=g.product("double_stranded","nuclei containing double-stranded cDNA and ATAC fragments",MaterialKind.MOLECULAR,state=population(copied=True,tagged=True,ds=True))
 g.react("second_strand","Make second-strand cDNA before lysis",StepKind.EXTENSION,((redistributed,InputRole.TEMPLATE),(second_pol,InputRole.CATALYST)),(ds,))
 cleanup=g.input("cleanup","6 µL Zymo DNA binding buffer (5 min), 24 µL AMPure XP, elute 12.5 µL EB",MaterialKind.SELECTION,source=SM)
 purified=g.product("purified","12.5 µL bead-recovered mixed cDNA/ATAC material",MaterialKind.MOLECULAR,state=ds.state)
 cwaste=g.product("cleanup_waste","AMPure supernatant and washes",MaterialKind.WASTE)
 g.react("purify","Stop with DNA binding buffer and AMPure-purify",StepKind.CLEANUP,((ds,InputRole.TEMPLATE),(cleanup,InputRole.SELECTION)),(purified,cwaste),changes_molecule=False)
 g.discard(cwaste,"physical bead supernatant and washes; molecular-class depletion unassigned")
 rna=g.product("rna_branch","6 µL mixed aliquot for RNA-directed reactions",MaterialKind.MOLECULAR,state=ds.state)
 atac=g.product("atac_branch","6 µL mixed aliquot for ATAC-directed reactions",MaterialKind.MOLECULAR,state=ds.state)
 remainder=g.product("aliquot_remainder","0.5 µL nominal eluate remainder",MaterialKind.MOLECULAR,state=ds.state)
 g.react("split","Aliquot 6 µL each to RNA and ATAC plates",StepKind.SPLIT,((purified,InputRole.TEMPLATE),),(rna,atac,remainder),changes_molecule=False)
 g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,"Nominal 12.5 µL eluate minus two 6 µL aliquots; fate is not assigned."))
 rna_transposome=g.input("rna_transposome","unindexed cDNA transposome ('i7 only TDE1', 25 nM, Illumina) in Nextera TD buffer, 55 °C 5 min; adaptor bases not in inspected sources",MaterialKind.ADAPTER,source=SM)
 rna_tn5=g.input("rna_tn5","cDNA Tn5 transposase",MaterialKind.CATALYST)
 tagged_rna=g.product("tagged_rna","tagmented cell-indexed cDNA",MaterialKind.MOLECULAR,state=population(copied=True,tagged=True,ds=True,rna_tagged=True))
 rwaste=g.product("rna_tn5_waste","candidate non-library cDNA tagmentation products",MaterialKind.MOLECULAR)
 g.react("tagment_rna","Tagment second-strand cDNA",StepKind.TAGMENTATION,((rna,InputRole.TEMPLATE),(rna_transposome,InputRole.STRUCTURAL),(rna_tn5,InputRole.CATALYST)),(tagged_rna,rwaste)); g.mark_terminal(Terminal(rwaste.id, Disposition.UNRESOLVED, 'Candidate RNA-side products lacking the RT-primer end remain in the RNA PCR pool. The end-class calculation does not establish absence of amplification or physical removal.'))
 # Only the RNA-directed aliquot receives the source SDS/BSA/P7 incubation.
 p7=g.input("p7","P7-[i7]-S7, printed PCR reverse primer",MaterialKind.OLIGO,state=Construct([seg("P7",il.P7),seg("PCR i7; length unassigned","[i7]","cbc",placeholder=True,feature=feature("pcr_i7","cell_barcode","combinatorial",group="cell_id",part="FACS well i7")),seg("S7",nx.S7)],name="sci-CAR P7 primer"),source=SM)
 release=g.input("rna_release","0.4 µL 1% SDS, 0.4 µL BSA; 55 °C 15 min with P7 primer, then Tween-20",MaterialKind.CONDITION,source=SM)
 released=g.product("rna_released","RNA-directed mixed aliquot after transposase release/quench",MaterialKind.MOLECULAR,state=tagged_rna.state)
 g.react("release_rna_tn5","Release Tn5 with SDS/BSA and P7; quench SDS with Tween-20",StepKind.THERMAL_TREATMENT,((tagged_rna,InputRole.TEMPLATE),(p7,InputRole.STRUCTURAL),(release,InputRole.CONDITION)),(released,),changes_molecule=False)
 for prefix,template,library,primers,arm,ratio,anneal,gap in (("rna",released,RNA,RNA_P,il.TRUSEQ_READ1,"0.8×",66,5),("atac",atac,ATAC,ATAC_P,nx.S5,"1×",63,3)):
  f=g.input(prefix+"_p5","printed indexed P5 primer",MaterialKind.OLIGO,state=Construct([seg("P5",il.P5),seg("PCR i5; length unassigned","[i5]","cbc",placeholder=True,feature=feature(prefix+"_pcr_i5","cell_barcode","combinatorial",group="cell_id",part="FACS/PCR well (i5)")),seg("source forward primer arm",arm)],name=prefix+" indexed P5 primer"),source=SM)
  pol=g.input(prefix+"_polymerase",f"NEBNext High-Fidelity 2× PCR Master Mix; 72 °C {gap} min, 18–22 cycles, {anneal} °C anneal",MaterialKind.CATALYST,source=SM)
  amp=g.product(prefix+"_amplified","indexed target products within mixed aliquot",MaterialKind.MOLECULAR,state=library)
  residual=g.product(prefix+"_pcr_residual","co-carried other modality/RNA and off-target amplification products",MaterialKind.MOLECULAR,state=template.state)
  g.react(prefix+"_pcr","Gap-fill/extend and index-PCR mixed "+prefix+" aliquot",StepKind.PCR,((template,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(p7,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amp,residual),geometry=PCRGeometry.unavailable(template.state,reason="Actual incoming mixed native Tn5 product is shown. Primer installation, native gap/end geometry and conditional reference assembly are not proven by an expected final library."))
  g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,"Both physical aliquots retain the mixture; modality-directed PCR does not prove class-specific removal or no amplification at downstream pooling/cleanup."))
  beads=g.input(prefix+"_beads",ratio+" AMPure XP after pooling",MaterialKind.SELECTION,source=SM)
  recovered=g.product(prefix+"_recovered","pooled bead-recovered indexed material",MaterialKind.LIBRARY if prefix=="rna" else MaterialKind.MOLECULAR,state=library)
  waste=g.product(prefix+"_cleanup_waste","bead supernatant/washes",MaterialKind.WASTE)
  g.react(prefix+"_cleanup","Pool and purify with "+ratio+" AMPure",StepKind.CLEANUP,((amp,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(recovered,waste),changes_molecule=False)
  g.discard(waste,"physical bead supernatant/washes; off-target polymer class recovery unassigned")
  if prefix=="atac":
   gel=g.input("atac_gel","gel extraction excluding fragments below 200 bp",MaterialKind.SELECTION,source=SM)
   final=g.product("atac_library","gel-recovered ATAC library",MaterialKind.LIBRARY,state=library)
   rejected=g.product("atac_gel_waste","gel size fraction not recovered",MaterialKind.WASTE)
   g.react("atac_size_select","Gel-extract ATAC fragments above 200 bp",StepKind.SELECTION,((recovered,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(final,rejected),changes_molecule=False)
   g.discard(rejected,"physically excluded gel size fraction")
  else:
   final=recovered
  lengths=READ_LENGTHS["RNA library" if prefix=="rna" else "ATAC library"]
  acquired=tuple(lengths) if prefix=="rna" else ("Read 1","Read 2","Index 1 (i7)","Index 2 (i5)")
  g.sequence(final,SequencingPlan(library=library,primers=tuple(primers),acquired_roles=acquired,read_lengths=lengths))
 g.declare_reviewed_scope(REVIEWED_TARGETS,SOURCE)
 return g.require_publishable()

REVIEWED_TARGETS=("nucleus","rt_primer","rt","reverse_transcribe","atac_transposome","atac_tn5","tagment_atac","facs_wells","redistribute","second_polymerase","second_strand","cleanup","purify","split","rna_transposome","rna_tn5","tagment_rna","p7","rna_release","release_rna_tn5","rna_p5","rna_polymerase","rna_pcr","rna_beads","rna_cleanup","rna_recovered","atac_p5","atac_polymerase","atac_pcr","atac_beads","atac_cleanup","atac_gel","atac_size_select","atac_library")
