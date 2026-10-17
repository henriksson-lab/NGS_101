"""Drop-BS molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import UNASSIGNED_LENGTH,nextera_library,seg
from chemdraw import Construct,Scene,feature,oligo,revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)

CELL_LEN=15
RANDOM_PRIMER="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"+"N"*9
P5_TRANSPOSON_PRIMER="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
P7_TRANSPOSON_BS_PRIMER="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTAAGTAGAAGATGGTATATGAGAT"
# Supplementary Fig. 1 (MOESM1 text lines 29-38): bead oligo
# 5'-[photocleavable linker]-PCR priming site-D15-ligation site-3'.
BEAD_PCR_SITE="CAAGCAGAAGACGGCATACGAGAT"
BEAD_LIGATION_END="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
BLOCKED_SPLINT="CTGTCTCTTATACACATCT"  # 5′→3′; source figure prints this strand 3′→5′
CELL_BARCODE=feature("dropbs_cell","cell_barcode","whitelist",whitelist="Drop-BS bead barcode set")

def bisulfite_unmethylated(seq):
 """Synthetic bead oligo: every C is unmethylated and reads as T after conversion."""
 return seq.replace("C","T")

# The barcode alphabet is A/G/T (Methods; Fig. S1), so conversion leaves it unchanged.
CONVERTED_ARM_RC=revcomp(bisulfite_unmethylated(BEAD_LIGATION_END))
CONVERTED_PCR_SITE_RC=revcomp(bisulfite_unmethylated(BEAD_PCR_SITE))
# The P7 transposon BS primer = Nextera Read 2 arm + the converted bead-site complement.
if P7_TRANSPOSON_BS_PRIMER!=BEAD_LIGATION_END+revcomp(CONVERTED_PCR_SITE_RC):
 raise AssertionError("P7 transposon BS primer does not prime the converted bead PCR site")

def bead_oligo():
 return Construct([seg("bead PCR priming site",BEAD_PCR_SITE,"p7"),
                   seg("15-nt cell barcode (A/G/T)","D"*CELL_LEN,"cbc",placeholder=True,feature=CELL_BARCODE),
                   seg("ligation site",BEAD_LIGATION_END,"me")],name="ChemGenes barcode-bead oligo (lot L142464-1)")

def blocked_complement():
 return Construct([seg("3′-end blocked complement oligo",BLOCKED_SPLINT,"me")],name="3′-SpC3 / 5′-Phos complement oligo")

def p5_transposon_primer():
 return Construct([seg("P5 transposon primer",P5_TRANSPOSON_PRIMER,"me")],name="P5 transposon primer")

def p7_transposon_bs_primer():
 return Construct([seg("Nextera Read 2 arm",BEAD_LIGATION_END,"me"),
                   seg("converted bead-site complement",revcomp(CONVERTED_PCR_SITE_RC),"p7")],name="P7 transposon BS primer")

def random_primer():
 return Construct([seg("S5+ME",P5_TRANSPOSON_PRIMER,"me"),seg("N9","N"*9,placeholder=True)],name="S5–ME–N9 random primer")

FINAL_LIBRARY,SEQ_PRIMERS=nextera_library([
 seg("bisulfite-converted genomic insert","X"*38,placeholder=True),
 seg("converted bead ligation-site complement",CONVERTED_ARM_RC,"me"),
 seg("15-nt cell barcode complement (C/A/T)","H"*CELL_LEN,"cbc",placeholder=True,feature=CELL_BARCODE),
 seg("converted bead PCR-site complement",CONVERTED_PCR_SITE_RC,"p7"),
],"Drop-BS library",
 # Supplementary Fig. 1 prints the N5XX/Ad2.XX primers with "[i5 index]"/"[i7 index]"
 # only, so neither index window is given a length.
 i5_nt=UNASSIGNED_LENGTH,i7_nt=UNASSIGNED_LENGTH)
TITLE="Drop-BS — droplet barcoding before bisulfite conversion"
NOTES="01_drop-bs.html"
SOURCE='<a href="https://doi.org/10.1038/s41467-023-40411-w">Zhang et al., <i>Nature Communications</i> (2023)</a>.'
SUMMARY="MNase-fragmented DNA and a photocleavable barcode bead meet in one droplet. Ligation attaches cell identity before a second droplet conversion, and random priming plus indexed PCR completes the library."
FINAL_CAPTION="Nextera-ended Drop-BS library carrying the bead-derived 15-nt cell barcode on the converted genomic molecule."
SEQUENCING_INTRO=("NovaSeq 6000 S4, 150 bp paired-end (Methods). Read 2 starts in the converted bead PCR-site complement, "
                  "then the 15-nt cell barcode and the converted ligation arm before genomic sequence; i5/i7 come from the N5XX/Ad2.XX indexing PCR, "
                  "whose selected index members and index lengths are not printed, so the index windows carry no length.")
READ_LENGTHS={"Read 1":150,"Read 2":150}

def bead_ligation_geometry():
 sc=Scene();sc.strand("released bead oligo",[seg("bead barcode and PCR arm","CAAGCAGAAGACGGCATACGAGAT"+"D"*15,"cbc",placeholder=True),seg("ligation-end arm",BEAD_LIGATION_END[:-len(BLOCKED_SPLINT)],"me"),seg("splint-paired terminus",BEAD_LIGATION_END[-len(BLOCKED_SPLINT):],"me")],label="photocleaved bead oligo")
 sc.anneal("blocked complement",[seg("3′-blocked ligation splint",BLOCKED_SPLINT,"me")],to="released bead oligo",pair=("3′-blocked ligation splint","splint-paired terminus"),label="3′-SpC3; 5′-phosphorylated")
 sc.ligation_end("released bead oligo","splint-paired terminus",kind="splinted",polarity="3-prime",owner="adapter",text="published splinted bead-oligo 3′ end")
 sc.strand("insert DNA",[seg("prepared insert end","X"*20,placeholder=True)],label="MNase fragment after end repair")
 sc.ligation_end("insert DNA","prepared insert end",kind="splinted",polarity="5-prime",owner="substrate",text="insert end presented to splint")
 return AdapterLigationGeometry.known("splinted",sc,statement="Drop-BS Supplementary Figure 1 prints the bead ligation sequence and the 3′-SpC3-blocked, 5′-phosphorylated complementary oligo that presents the ligation site.")
def oligos():
 yield oligo("Barcode-bead oligo (5′ photocleavable linker)",list(bead_oligo()))
 yield oligo("3′-end blocked complement oligo",list(blocked_complement()),mods="/5Phos/",three="/3SpC3/-3′")
 yield oligo("S5–ME–N9 random primer",list(random_primer()))
 yield oligo("P5 transposon pre-PCR primer",list(p5_transposon_primer()))
 yield oligo("P7 transposon BS pre-PCR primer",list(p7_transposon_bs_primer()))
def protocol_graph():
 source="Defining source: https://doi.org/10.1038/s41467-023-40411-w"
 graph=ProtocolGraph("Drop-BS",authority=source)
 nucleus=graph.input("nucleus","single nucleus",MaterialKind.BIOLOGICAL,primary=True)
 mnase=graph.input("mnase","lysis buffer with MNase (Thermo 88216) and Ca2+",MaterialKind.CATALYST)
 fragments=graph.product("mnase_fragments","MNase-fragmented genomic DNA",MaterialKind.MOLECULAR)
 graph.react("mnase_fragmentation","Lyse and fragment one nucleus with MNase in a droplet (ice, RT 15 min, 65 °C 30 min)",
             StepKind.FRAGMENTATION,
             ((nucleus,InputRole.TEMPLATE),(mnase,InputRole.CATALYST)),(fragments,))
 bead=graph.input("barcode_bead","photocleavable barcode-bead oligo",MaterialKind.ADAPTER,state=bead_oligo())
 complement=graph.input("blocked_complement","3′-SpC3-blocked, 5′-phosphorylated 19-nt complement oligo",MaterialKind.OLIGO,state=blocked_complement())
 ready_bead=graph.product("ligation_ready_bead","bead oligo with double-stranded ligation tail",MaterialKind.MOLECULAR)
 graph.react("anneal_complement","Anneal the blocked complement to the bead oligos (85 °C 5 min, slow cool)",
             StepKind.HYBRIDIZATION,((bead,InputRole.TEMPLATE),(complement,InputRole.STRUCTURAL)),(ready_bead,))
 fusion=graph.input("fusion_condition","AC-field droplet fusion (300 Vp-p, 200 kHz)",MaterialKind.CONDITION)
 paired=graph.product("paired_droplet","fragment-and-barcode-bead droplet",MaterialKind.MOLECULAR)
 excluded=graph.product("fusion_loss","unpaired droplet material",MaterialKind.WASTE)
 graph.react("droplet_fusion","Fuse the fragment droplet with a bead/enzyme/ligation-mix droplet",
             StepKind.SAMPLING,
             ((fragments,InputRole.TEMPLATE),(fusion,InputRole.CONDITION)),(paired,excluded),changes_molecule=False)
 graph.discard(excluded,"unpaired or unusable fused droplets are not retained")
 repair=graph.input("repair_enzymes","End-It enzyme mix (Lucigen) with dNTPs and ATP",MaterialKind.CATALYST)
 repaired=graph.product("repaired_fragments","end-repaired genomic fragments in droplets",
                        MaterialKind.MOLECULAR)
 graph.react("end_repair","End-repair the MNase fragments",StepKind.END_REPAIR,
             ((paired,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
 uv=graph.input("uv_release","365 nm UV, 160 mJ/cm2, 8 min (photocleaves the bead linker)",MaterialKind.CONDITION)
 ligase=graph.input("ligase","Fast-Link DNA ligase (Lucigen), 16 °C overnight",MaterialKind.CATALYST)
 barcoded=graph.product("barcoded_dna","bead-barcode-ligated genomic DNA",MaterialKind.MOLECULAR)
 graph.react("barcode_ligation","Photorelease and ligate the bead barcode",
             StepKind.LIGATION,
             ((repaired,InputRole.TEMPLATE),(ready_bead,InputRole.STRUCTURAL),(uv,InputRole.CONDITION),
              (ligase,InputRole.CATALYST)),(barcoded,),geometry=bead_ligation_geometry())
 recovery=graph.input("recovery","70 °C inactivation, perfluorooctanol break, proteinase K, 1.4× then 1× SPRIselect",MaterialKind.SELECTION)
 recovered=graph.product("recovered_dna","purified barcoded DNA",MaterialKind.MOLECULAR)
 unbound=graph.product("ligation_cleanup_loss","oil, proteins and SPRI-excluded material",MaterialKind.WASTE)
 graph.react("recover_barcoded","Break droplets and purify barcoded DNA",StepKind.CLEANUP,
             ((barcoded,InputRole.TEMPLATE),(recovery,InputRole.SELECTION)),(recovered,unbound),changes_molecule=False)
 graph.discard(unbound,"removed by droplet breaking, proteinase K digestion and SPRIselect cleanup")
 bisulfite=graph.input("bisulfite","MethylEdge bisulfite reagent in droplets (95 °C 5 min, 54 °C 1 h) after 95 °C denaturation",MaterialKind.CONDITION)
 converted=graph.product("converted_dna","bisulfite-converted barcoded DNA",MaterialKind.MOLECULAR)
 graph.react("bisulfite_conversion","Bisulfite-convert the barcoded DNA in droplets",
             StepKind.CHEMICAL_CONVERSION,
             ((recovered,InputRole.TEMPLATE),(bisulfite,InputRole.CONDITION)),(converted,))
 desulph=graph.input("desulphonation","MethylEdge desulphonation/column cleanup, 1.1× SPRIselect",MaterialKind.SELECTION)
 clean_converted=graph.product("clean_converted_dna","desulphonated converted single-stranded DNA",MaterialKind.MOLECULAR)
 bs_waste=graph.product("bisulfite_waste","bisulfite reagent and column/SPRI flow-through",MaterialKind.WASTE)
 graph.react("desulphonate","Desulphonate and clean up converted DNA",StepKind.CLEANUP,
             ((converted,InputRole.TEMPLATE),(desulph,InputRole.SELECTION)),(clean_converted,bs_waste))
 graph.discard(bs_waste,"removed in the kit column and SPRIselect cleanup")
 rp=graph.input("random_primer","S5–ME–N9 random primer",MaterialKind.OLIGO,state=random_primer())
 klenow=graph.input("klenow","Klenow exo− (Enzymatics), slow ramp 4→37 °C, 37 °C 60 min",MaterialKind.CATALYST)
 extended=graph.product("extended_dna","random-primer-extended converted DNA",MaterialKind.MOLECULAR)
 graph.react("random_extension","Random-prime converted DNA to install the S5–ME handle",
             StepKind.EXTENSION,
             ((clean_converted,InputRole.TEMPLATE),(rp,InputRole.STRUCTURAL),
              (klenow,InputRole.CATALYST)),(extended,))
 rp_clean=graph.input("rp_cleanup","1.1× SPRIselect, 95 °C denaturation, second 1.1× SPRIselect",MaterialKind.SELECTION)
 ss_copy=graph.product("single_stranded_copy","denatured random-primed copy strands",MaterialKind.MOLECULAR)
 rp_waste=graph.product("random_priming_waste","free primer and SPRI-excluded material",MaterialKind.WASTE)
 graph.react("denature_cleanup","Denature and purify the random-primed products",StepKind.CLEANUP,
             ((extended,InputRole.TEMPLATE),(rp_clean,InputRole.SELECTION)),(ss_copy,rp_waste))
 graph.discard(rp_waste,"removed by the two 1.1× SPRIselect cleanups")
 p5t=graph.input("p5_transposon_primer","P5 transposon primer",MaterialKind.OLIGO,state=p5_transposon_primer())
 p7t=graph.input("p7_transposon_bs_primer","P7 transposon BS primer (primes the converted bead site)",MaterialKind.OLIGO,state=p7_transposon_bs_primer())
 kapa=graph.input("kapa_uracil","KAPA HiFi HotStart Uracil+ ReadyMix with DMSO",MaterialKind.CATALYST)
 pre=graph.product("pre_pcr_product","pre-PCR product (S5–ME … converted bead site–Read 2 arm)",MaterialKind.MOLECULAR)
 graph.react("pre_pcr","Pre-amplify with P5/P7 transposon primers (5 cycles), 1× SPRIselect",StepKind.PCR,
             ((ss_copy,InputRole.TEMPLATE),(p5t,InputRole.STRUCTURAL),(p7t,InputRole.STRUCTURAL),
              (kapa,InputRole.CATALYST)),(pre,),
             geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason="The pre-PCR template is not drawn as an exact molecular state; the P7 transposon BS primer's landing on the converted bead site is derived in the module."))
 i5=graph.input("i5_primer","N5XX primer family (i5 member not printed)",MaterialKind.OLIGO)
 i7=graph.input("i7_primer","Ad2.XX primer family (i7 member not printed)",MaterialKind.OLIGO)
 kapa2=graph.input("kapa_uracil_index","KAPA HiFi HotStart Uracil+ ReadyMix with EvaGreen",MaterialKind.CATALYST)
 amplified=graph.product("indexed_library","indexed Drop-BS library before cleanup",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
 graph.react("library_pcr","Index PCR with N5XX/Ad2.XX (qPCR-monitored)",StepKind.PCR,
             ((pre,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
              (i7,InputRole.STRUCTURAL),(kapa2,InputRole.CATALYST)),(amplified,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="Index-primer members are published only as N5XX/Ad2.XX families."))
 final_clean=graph.input("final_spri","0.85× SPRIselect",MaterialKind.SELECTION)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 small=graph.product("final_spri_loss","primers and short fragments not retained at 0.85×",MaterialKind.WASTE)
 graph.react("final_cleanup","Purify the indexed library",StepKind.CLEANUP,
             ((amplified,InputRole.TEMPLATE),(final_clean,InputRole.SELECTION)),(library,small),changes_molecule=False)
 graph.discard(small,"left in the 0.85× SPRIselect supernatant")
 graph.annotate(SourceDiscrepancy("random_extension","Copied cell-barcode symbol in Supplementary Fig. 1",
     (SourceClaim("Methods 'Design and preparation of barcode beads' (PMC10400590 text line193) and Fig. S1 label","15-nt barcode of A, T, G bases, drawn D15"),
      SourceClaim("Supplementary Fig. 1 random-priming and pre-PCR strands (MOESM1 text lines84–100)","copy strand drawn V15 (IUPAC A/C/G)")),
     "Use the complement of the stated A/G/T alphabet (C/A/T, IUPAC H) for the copied barcode.",
     "Both symbols are printed in the published supplement; the complement of D is H, not V.",
     "Likely a labelling slip that does not change any printed base; not confirmed by the authors."))
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
     acquired_roles=tuple(p.role for p in SEQ_PRIMERS),read_lengths=READ_LENGTHS,
     read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in SEQ_PRIMERS}))
 graph.declare_reviewed_scope((
     "nucleus","mnase","mnase_fragmentation","barcode_bead","blocked_complement","ligation_ready_bead",
     "anneal_complement","fusion_condition","droplet_fusion","paired_droplet","fusion_loss",
     "repair_enzymes","end_repair","repaired_fragments","uv_release","ligase","barcode_ligation","barcoded_dna",
     "recovery","recover_barcoded","recovered_dna","ligation_cleanup_loss","bisulfite",
     "bisulfite_conversion","converted_dna","desulphonation","desulphonate","clean_converted_dna","bisulfite_waste",
     "random_primer","klenow","random_extension","extended_dna","rp_cleanup","denature_cleanup","single_stranded_copy",
     "random_priming_waste","p5_transposon_primer","p7_transposon_bs_primer","kapa_uracil","pre_pcr","pre_pcr_product",
     "i5_primer","i7_primer","kapa_uracil_index","library_pcr","indexed_library","final_spri","final_cleanup",
     "library","final_spri_loss"),source)
 return graph.require_publishable()
