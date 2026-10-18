"""scSNV-seq transcript, barcode and genotype branches."""
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from protocol_model import role
import seqprimers as sp
import illumina as il
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind

TITLE = "scSNV-seq — genotype-linked single-cell transcriptomes"
NOTES = "01_scsnv-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1186/s13059-024-03169-y">Cooper et al., <i>Genome Biology</i> (2024)</a>.'
SUMMARY = "A guide-associated 6-nt iBAR and an independent 12-nt PuroR barcode connect targeted DNA genotypes to 10x 5′ transcriptomes."
CAVEAT = ("Supplementary Table S2 prints the vector-construction primers, both custom "
          "10x spike-in RT primers, the shared 10-base-indexed nested primer, the JAK1 "
          "amplicon panel and the integrated pKLV2 region (guide 19N, scaffold, iBAR 6N, "
          "PuroR 5′ UTR barcode 12N). Those sequences are exact below, and both spike-in "
          "RT arms are checked against that vector sequence. INFERRED: the paired "
          "commercial 10x 5′ HT strand assembly and reverse nested primer, and the Mission "
          "Bio Tapestri native endpoint. Cell-barcode/UMI widths, native index geometry and acquired cycles remain unassigned here. Genotyping uses separately sampled expanded cells, with clone linkage computational through vector barcodes. Sequencing was on a "
          "NovaSeq 6000; study read allocation remains unassigned in this capped model.")

ULTRAMER_5UTR = "CTTGTCTCGCTCCGGGGAACGCTCGGAAACTCCCGGCCGCCGCCACCCGCGTCTGTTCTGTTACACAAGGGAAGAAAAGCCGCTGCCGCACTCCGAGTGT"
UTR_BC_F1 = "CCGGGCCTTTCGACCTGCATCCATCTTGTCTCGCTCCGGGGAACG"
UTR_BC_R1 = "CTCGGTCATGGTAAGCTTCAGCTGCNNNNNNNNNNNNACACTCGGAGTGCGGCAG"
IBAR_PCR1_F = "TGTGGAAAGGACGAAACACC"
IBAR_PCR1_R = "CAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTNNNNNNTCCAGCATAGCTCTTAAAC"
IBAR_GIBSON_F = "GGCTTTATATATCTTGTGGAAAGGACGAAACACC"
IBAR_GIBSON_R = "GAAAAGCGCCTCCCCTACCCGGTAGAATTGGATCCAAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCC"
GUIDE_RT = "AAGCAGTGGTATCAACGCAGAGTACCAAGTTGATAACGGACTAGCC"
PURO_RT = "AAGCAGTGGTATCAACGCAGAGTACGGCTTGTACTCGGTCATGGT"
I5_LONG = "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT"
# Supplementary Table S2 publishes the integrated pKLV2 region (guide 19N, scaffold,
# iBAR 6N, PuroR 5'UTR barcode 12N); these are verbatim stretches of it.
SCAFFOLD_LOOP1 = "GTTTAAGAGCTATGCTGGA"
SCAFFOLD_AFTER_IBAR = "AACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTG"
PURO_UTR_TAIL = "CTTGTCTCGCTCCGGGGAACGCTCGGAAACTCCCGGCCGCCGCCACCCGCGTCTGTTCTGTTACACAAGGGAAGAAAAGCCGCTGCCGCACTCCGAGTGT"
PURO_CDS_START = "GCAGCTGAAGCTTACCATGACCGAGTACAAGCC"
# Each spike-in RT primer is a SMART handle plus a gene-specific 3' arm; the arm must be
# the reverse complement of the transcribed vector sequence it primes on.
SMART_HANDLE = "AAGCAGTGGTATCAACGCAGAGTAC"
GUIDE_RT_ARM = GUIDE_RT[len(SMART_HANDLE):]
PURO_RT_ARM = PURO_RT[len(SMART_HANDLE):].upper()
if revcomp(GUIDE_RT_ARM) not in SCAFFOLD_AFTER_IBAR:
    raise AssertionError("oJR160 does not prime on the published scaffold")
if revcomp(PURO_RT_ARM) not in PURO_CDS_START:
    raise AssertionError("RT_Puro_Singlecell does not prime on the published PuroR start")
# Both nested forward primers are P5 + a 10-base library index + the TruSeq Read 1 arm.
if I5_LONG != il.P5 + "N" * 10 + il.TRUSEQ_READ1:
    raise AssertionError("printed i5long primer is not P5 + 10N index + TruSeq Read 1")
I5_INDEX_NT = 10
NESTED_FWD = PCRPrimer("Kosuke_Guide_i5longFwd / Puro_i5longFwd", I5_LONG,
                       anneal_nt=len(il.TRUSEQ_READ1))

CELL=feature("chromium_cell_barcode","cell_barcode","unknown",note="10x5′HT identity; native width unassigned")
UMI_FEATURE=feature("chromium_umi","umi","random",note="native width unassigned")
CB=lambda:role("Chromium cell barcode; native length/position unassigned",molecular_feature=CELL)
UMI=lambda:role("Chromium UMI; native length/position unassigned",molecular_feature=UMI_FEATURE)
GUIDE_INSERT = [CB(), UMI(),
    Segment("vector transcription-start G", "G"),
    Segment("guide spacer", "N" * 19, "cbc", placeholder=True,
            feature=feature("guide_barcode", "guide_barcode", "unknown")),
    Segment("scaffold loop 1", SCAFFOLD_LOOP1),
    Segment("iBAR", "I" * 6, "cbc", placeholder=True,
            feature=feature("ibar", "guide_barcode", "random")),
    Segment("scaffold to the oJR160 priming site", SCAFFOLD_AFTER_IBAR)]
PURO_INSERT = [CB(), UMI(),
    Segment("PuroR 5' UTR (NeoUTR3 ultramer)", PURO_UTR_TAIL),
    Segment("PuroR barcode", "B" * 12, "cbc", placeholder=True,
            feature=feature("puror_barcode", "feature_barcode", "random")),
    Segment("PuroR start to the RT-primer site", PURO_CDS_START)]
TAPESTRI_CELL=feature("tapestri_cell_barcode","cell_barcode","unknown",note="MissionBio DNA-cell identity; distinct from Chromium cell barcode")
SAMPLE=feature("scsnv_nested_i5","sample_index","random",note="TableS2 prints N10 indexing barcode; selected member unassigned")
GENO_INSERT=[role("MissionBio DNA-cell barcode; native length/position unassigned",molecular_feature=TAPESTRI_CELL),role("target-dependent JAK1 or vector-barcode amplicon; separate target populations")]
RNA=Construct([role("native10x end A"),CB(),UMI(),role("transcript-derived insert"),role("native10x end B")],name="scSNV transcriptome native boundary")
def _nested_endpoint(payload,name):
 return Construct([Segment("source nested P5",il.P5),Segment("source N10 library index","N"*10,"cbc",placeholder=True,feature=SAMPLE),Segment("source nested forward arm",il.TRUSEQ_READ1),role("native assembly between nested arm and target; unassigned"),*payload,role("reverse native end; primer geometry unassigned")],name=name)
GUIDE=_nested_endpoint(GUIDE_INSERT,"scSNV guide/iBAR conditional native endpoint")
PURO=_nested_endpoint(PURO_INSERT,"scSNV PuroR conditional native endpoint")
GENO=Construct([role("nativeTapestri end A"),*GENO_INSERT,role("nativeTapestri end B")],name="scSNV DNA-genotype native endpoint")
TENX_PRIMERS=GUIDE_PRIMERS=PURO_PRIMERS=GENO_PRIMERS=()
FINAL_LIBRARIES=(("Transcriptome library",RNA,(),"Native10x5′HT reference; complete strands and run allocation unassigned.","Chromium identity links RNA-derived modalities."),("Guide and iBAR library",GUIDE,(),"Conditional ordered target and printed nested forward primer; intervening native assembly unassigned.","Guide19N and iBAR6N retain vector lineage."),("PuroR-barcode library",PURO,(),"Conditional PuroR target and printed nested forward primer; native reverse end unassigned.","PuroR12N links clonal descendants."),("Targeted genotype library",GENO,(),"Distinct Tapestri cells sampled from the expanded population; no Chromium barcode length or endpoint is inherited.","Clone linkage through vector identifiers is computational."))
SHOW_READ_CONTENT=False
READ_SPECS={name:tuple(sp.read_spec(r,"unassigned","NovaSeq6000 reported; native read/index acquisition unassigned in this source scope","undetermined") for r in ("Read 1","Read 2","Index 1 (i7)","Index 2 (i5)")) for name,*_ in FINAL_LIBRARIES}


def oligos():
    return (
        f"<b>Ultramer_5UTR</b> — 5′-{ULTRAMER_5UTR}-3′",
        f"<b>5UTR_BC_F1_Ext</b> — 5′-{UTR_BC_F1}-3′",
        f"<b>5UTR_BC_R1_Ext</b> — 5′-{UTR_BC_R1}-3′",
        f"<b>iBAR_PCR1_F</b> — 5′-{IBAR_PCR1_F}-3′",
        f"<b>iBAR_PCR1_R</b> — 5′-{IBAR_PCR1_R}-3′",
        f"<b>iBar_Gibson_Nested_F</b> — 5′-{IBAR_GIBSON_F}-3′",
        f"<b>iBar_Gibson_Nested_R</b> — 5′-{IBAR_GIBSON_R}-3′",
        f"<b>oJR160 guide RT spike-in</b> — 5′-{GUIDE_RT}-3′",
        f"<b>RT_Puro_Singlecell spike-in</b> — 5′-{PURO_RT}-3′",
        f"<b>Kosuke_Guide_i5longFwd / Puro_i5longFwd</b> — 5′-{I5_LONG}-3′",
    )

def population(*,copied=False,amplified=False,bead=False):
 sc=Scene()
 sc.strand("gDNA",[role("genomic DNA; separate polymer")],label="retained DNA; no Tapestri routing from Chromium droplets")
 for kind,payload,primer,arm in (("guide",GUIDE_INSERT[2:],GUIDE_RT,GUIDE_RT_ARM),("puro",PURO_INSERT[2:],PURO_RT,PURO_RT_ARM)):
  sc.strand(kind+"RNA",payload,label=kind+" transcribed source vector population")
  if copied:
   # RT-primer arm already supplies the priming-site sequence; copy only the upstream target.
   last=payload[-1]
   if not last.top.endswith(revcomp(arm)):raise ValueError("source RT arm is not terminal target complement")
   before=[*payload[:-1],Segment(last.name+" upstream of RT site",last.top[:-len(arm)])]
   sc.strand(kind+"cDNA",[Segment("source "+kind+" RT primer",primer),*complement_segments(before),role("native barcoding/strand continuation; topology unassigned")],label="native amplified target reference" if amplified else "ss RT-primer-derived cDNA; barcoding extent unassigned")
   if amplified:sc.strand(kind+" complementary copies",[role("complementary "+kind+" PCR-copy population; native strand boundaries unassigned")],label="separate amplified opposite-strand class")
 sc.strand("mRNA",[role("other poly(A) transcripts")],label="other RNA population")
 if copied:sc.strand("mRNAcDNA",[CB(),UMI(),role("native transcript cDNA; complete assembly unassigned")],label="separate native10x cDNA class")
 if bead or copied:sc.strand("bead",[CB(),UMI(),role("native10x bead oligo arrangement unassigned")],label="bead barcoding component; no assumed16/10widths")
 return sc


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 expanded=g.input("expanded","expanded edited population carrying source iBAR6N/PuroR12N guide cassette; source lentiviral integration and doxycycline-induced base editing precede these library branches",MaterialKind.BIOLOGICAL,state=population(),primary=True)
 chromium=g.product("chromium_cells","cells sampled for Chromium5′HT transcriptomics",MaterialKind.BIOLOGICAL,state=population())
 tapestri=g.product("tapestri_cells","separately sampled/expanded cells for Tapestri genotyping",MaterialKind.BIOLOGICAL,state=population())
 g.react("sample_populations","Sample expanded population for separate Chromium and Tapestri experiments",StepKind.SPLIT,((expanded,InputRole.TEMPLATE),),(chromium,tapestri),changes_molecule=False)
 partition=g.input("partition_condition","ChromiumX5′HT droplet partitioning",MaterialKind.CONDITION)
 gem=g.product("gem","Chromium reaction containing independent cellular DNA/RNAs and bead oligos",MaterialKind.MOLECULAR,state=population(bead=True))
 other_gems=g.product("other_gems","other native droplets/material outside representative cellular reaction",MaterialKind.MOLECULAR)
 g.react("partition","Partition transcriptomic cells with native5′HT beads",StepKind.PARTITION,((chromium,InputRole.TEMPLATE),(partition,InputRole.CONDITION)),(gem,other_gems),changes_molecule=False)
 g.mark_terminal(Terminal(other_gems.id,Disposition.UNRESOLVED,"Other partition contents/recovery not assigned; computational cell calling is not physical discard."))
 gp=g.input("guide_rt_primer","oJR160 guide RT spike-in",MaterialKind.OLIGO,state=Construct([Segment("oJR160",GUIDE_RT)],name="oJR160"))
 pp=g.input("puro_rt_primer","RT_Puro_Singlecell spike-in",MaterialKind.OLIGO,state=Construct([Segment("RT_Puro_Singlecell",PURO_RT)],name="RT_Puro_Singlecell"))
 rt=g.input("rt","native10x5′HT RT reagents; component identity unassigned",MaterialKind.CATALYST)
 bead=g.input("bead","native10x5′HT bead barcoding oligos; strands/widths unassigned",MaterialKind.OLIGO)
 copied=g.product("copied","native RT/barcoding mixture; ss target copies distinct from final libraries",MaterialKind.MOLECULAR,state=population(copied=True))
 g.react("rt_capture","Reverse-transcribe native10x and custom guide/PuroR templates",StepKind.REVERSE_TRANSCRIPTION,((gem,InputRole.TEMPLATE),(gp,InputRole.STRUCTURAL),(pp,InputRole.STRUCTURAL),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))
 f=g.input("cdna_forward","native10x cDNA amplification primer A",MaterialKind.OLIGO)
 r=g.input("cdna_reverse","native10x cDNA amplification primer B; species relationship unassigned",MaterialKind.OLIGO)
 pol=g.input("cdna_polymerase","native10x cDNA amplification polymerase",MaterialKind.CATALYST)
 amplified=g.product("amplified","native amplified cDNA population",MaterialKind.MOLECULAR,state=population(copied=True,amplified=True))
 g.react("cdna_pcr","Amplify cDNA according to named5′HT workflow",StepKind.PCR,((copied,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(copied.state,reason="Actual native RT mixture; complete kit primer/strand assembly unassigned."))
 selection=g.input("spri","source SPRI fractionation after cDNA amplification; exact ratio unassigned",MaterialKind.SELECTION)
 large=g.product("large","long cDNA fraction for transcriptome construction",MaterialKind.MOLECULAR,state=amplified.state)
 small=g.product("small","small SPRI fraction for guide/iBAR and PuroR nested enrichment",MaterialKind.MOLECULAR,state=amplified.state)
 g.react("fractionate","Retain transcriptome and small SPRI fractions",StepKind.SELECTION,((amplified,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),(large,small),changes_molecule=False)
 kit=g.input("rna_kit","named10x5′HT transcriptome construction reagents",MaterialKind.CONDITION)
 rna=g.product("rna_library","native transcriptome endpoint",MaterialKind.LIBRARY,state=RNA)
 g.react("rna_construct","Native5′HT transcriptome construction boundary",StepKind.MAGIC,((large,InputRole.TEMPLATE),(kit,InputRole.CONDITION)),(rna,),reason="Published Methods delegates transcriptome construction to manufacturer; exact historical operations/oligos remain outside this two-original cap.")
 guides=g.product("guide_aliquot","mixed small-fraction aliquot for guide/iBAR",MaterialKind.MOLECULAR,state=small.state)
 puros=g.product("puro_aliquot","mixed small-fraction aliquot for PuroR",MaterialKind.MOLECULAR,state=small.state)
 g.react("split_small","Aliquot small fraction for analogous nested enrichments",StepKind.SPLIT,((small,InputRole.TEMPLATE),),(guides,puros),changes_molecule=False)
 for key,template,endpoint in (("guide",guides,GUIDE),("puro",puros,PURO)):
  first_f=g.input(key+"_first_forward","target-specific first enrichment primer; selected strand unassigned",MaterialKind.OLIGO)
  first_r=g.input(key+"_first_reverse","native first enrichment opposite primer; unassigned",MaterialKind.OLIGO)
  enzyme=g.input(key+"_polymerase","nested enrichment PCR polymerase; identity unassigned",MaterialKind.CATALYST)
  first=g.product(key+"_first","first enrichment products within mixed aliquot",MaterialKind.MOLECULAR,state=template.state)
  g.react(key+"_first_pcr","First "+key+" enrichment stage implied by source nested PCR",StepKind.PCR,((template,InputRole.TEMPLATE),(first_f,InputRole.STRUCTURAL),(first_r,InputRole.STRUCTURAL),(enzyme,InputRole.CATALYST)),(first,),geometry=PCRGeometry.unavailable(template.state,reason="Mixed small-fraction precursor; first-round primer species, cycle counts and native products remain unassigned."))
  nested=g.input(key+"_nested_forward","TableS2 shared N10 i5long nested forward",MaterialKind.OLIGO,state=Construct([Segment("P5",il.P5),Segment("N10index","N"*10,"cbc",placeholder=True,feature=SAMPLE),Segment("nested3′arm",il.TRUSEQ_READ1)],name=key+" nested forward"))
  reverse=g.input(key+"_nested_reverse","native nested opposite primer; selected strand unassigned",MaterialKind.OLIGO)
  lib=g.product(key+"_library","conditional source-target nested endpoint",MaterialKind.LIBRARY,state=endpoint)
  other=g.product(key+"_residual","co-carried other transcripts/DNA and non-target PCR products",MaterialKind.MOLECULAR,state=first.state)
  g.react(key+"_nested_pcr","Nested "+key+" PCR using printed i5long forward primer",StepKind.PCR,((first,InputRole.TEMPLATE),(nested,InputRole.STRUCTURAL),(reverse,InputRole.STRUCTURAL),(enzyme,InputRole.CATALYST)),(lib,other),geometry=PCRGeometry.unavailable(first.state,reason="Actual first-enrichment mixture; reverse primer/native assembly unassigned. Final reference is not a PCR input."))
  g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,"Nested primer selectivity is not physical modality purification; native off-target/RNA/DNA recovery remains unassigned."))
  g.sequence(lib,SequencingPlan(library=endpoint,limitation="NovaSeq6000; native acquired reads/index cycles and primer landing unassigned."))
 dna_kit=g.input("tapestri_kit","Tapestri native custom amplicon workflow; TableS2 JAK1 and vector target panel",MaterialKind.CONDITION)
 dna=g.product("genotype_library","native Tapestri genotype endpoint from separate cells",MaterialKind.LIBRARY,state=GENO)
 g.react("tapestri_construct","Native Tapestri targeted genotyping boundary",StepKind.MAGIC,((tapestri,InputRole.TEMPLATE),(dna_kit,InputRole.CONDITION)),(dna,),reason="Source names manufacturer custom-amplicon workflow; native cell barcode/primer assembly and internal reaction order remain deferred.")
 g.sequence(dna,SequencingPlan(library=GENO,limitation="Distinct MissionBio DNA-cell identities, computationally linked to Chromium clonal descendants via vector identifiers; read/index geometry unassigned."))
 g.sequence(rna,SequencingPlan(library=RNA,limitation="Named10x5′HT endpoint; native run geometry/cycles unassigned."))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()

REVIEWED=("expanded","sample_populations","partition_condition","partition","guide_rt_primer","puro_rt_primer","rt","bead","rt_capture","cdna_forward","cdna_reverse","cdna_polymerase","cdna_pcr","spri","fractionate","rna_kit","rna_construct","split_small","guide_first_forward","guide_first_reverse","guide_polymerase","guide_first_pcr","guide_nested_forward","guide_nested_reverse","guide_nested_pcr","guide_library","puro_first_forward","puro_first_reverse","puro_polymerase","puro_first_pcr","puro_nested_forward","puro_nested_reverse","puro_nested_pcr","puro_library","tapestri_kit","tapestri_construct","genotype_library","rna_library")
