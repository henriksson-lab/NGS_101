"""scSNV-seq transcript, barcode and genotype branches."""
from batch_ngs import truseq_amplicon_library
from chemdraw import Row, Segment, feature, revcomp
import illumina as il
from multibranch import illumina_branch
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind

TITLE = "scSNV-seq — genotype-linked single-cell transcriptomes"
NOTES = "01_scsnv-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1186/s13059-024-03169-y">Cooper et al., <i>Genome Biology</i> (2024)</a>.'
SUMMARY = "A guide-associated 6-nt iBAR and an independent 12-nt PuroR barcode connect targeted DNA genotypes to 10x 5′ transcriptomes."
CAVEAT = ("Supplementary Table S2 prints the vector-construction primers, both custom "
          "10x spike-in RT primers, the shared 10-base-indexed nested primer, the JAK1 "
          "amplicon panel and the integrated pKLV2 region (guide 19N, scaffold, iBAR 6N, "
          "PuroR 5′ UTR barcode 12N). Those sequences are exact below, and both spike-in "
          "RT arms are checked against that vector sequence. INFERRED: the paired "
          "commercial 10x 5′ HT bead construct and i7-side nested primer, and the Mission "
          "Bio Tapestri shell including its cell-barcode length. Sequencing was on a "
          "NovaSeq 6000; read lengths are not published.")

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

CB = lambda: Segment("cell barcode", "B" * 16, "cbc", placeholder=True,
                     feature=feature("cell_barcode", "cell_barcode", "whitelist", whitelist="10x Chromium 5-prime HT cell-barcode whitelist"))
UMI = lambda: Segment("UMI", "U" * 10, "umi", placeholder=True,
                      feature=feature("umi", "umi", "random"))
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
GENO_INSERT = [Segment("Tapestri cell barcode", "B" * 16, "cbc", placeholder=True,
                       inferred=True,
                       feature=feature("cell_barcode", "cell_barcode", "whitelist", whitelist="10x Chromium 5-prime HT cell-barcode whitelist")),
               Segment("JAK1, iBAR or PuroR amplicon", "X" * 34, placeholder=True)]

RNA, TENX_PRIMERS = illumina_branch("scSNV-seq transcriptome library", "10x 5-prime cDNA", cell=16, umi=10, inferred_adapters=True)
GUIDE, GUIDE_PRIMERS = truseq_amplicon_library(GUIDE_INSERT, "scSNV-seq guide/iBAR library", i5_nt=I5_INDEX_NT, inferred_adapters=True)
PURO, PURO_PRIMERS = truseq_amplicon_library(PURO_INSERT, "scSNV-seq PuroR-barcode library", i5_nt=I5_INDEX_NT, inferred_adapters=True)
GENO, GENO_PRIMERS = truseq_amplicon_library(GENO_INSERT, "scSNV-seq targeted genotype library", inferred_adapters=True)
FINAL_LIBRARIES = (
    ("Transcriptome library", RNA, TENX_PRIMERS, "10x gene-expression reads carry cell barcode and UMI.", "The transcript branch uses the declared 10x-compatible run primers."),
    ("Guide and iBAR library", GUIDE, GUIDE_PRIMERS, "The nested enrichment reads the droplet cell barcode and UMI, then the guide spacer, the published scaffold loop and its 6-nt iBAR, up to the oJR160 priming site. The P5 arm is the printed nested primer with its 10-base library index; the P7 arm is INFERRED.", "The guide branch retains the same droplet cell identity."),
    ("PuroR-barcode library", PURO, PURO_PRIMERS, "The analogous PuroR enrichment reads the published 5' UTR, the 12-nt PuroR barcode and the PuroR start up to the RT-primer site. The P5 arm is the printed nested primer; the P7 arm is INFERRED.", "The PuroR branch retains the same droplet cell identity."),
    ("Targeted genotype library", GENO, GENO_PRIMERS, "INFERRED — Tapestri amplicons (JAK1 panel plus the iBAR and PuroR amplicons) call the edited allele; the Mission Bio cell-barcode length and sequencing shell are not published here.", "The DNA branch supplies the genotype joined computationally through the barcode pair."),
)


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

def sections():
    return [
        ("Construct the doubly barcoded guide library", [Row(chunks=[("PuroR 5′ UTR—[12-nt barcode]    guide—[6-nt iBAR]", "cbc", False)])], "Both transcribed identifiers are integrated with the guide cassette."),
        ("Base-edit and expand cells", [Row(chunks=[("guide + base editor → edited clone descendants", "r1", False)])], "Descendants of one integration inherit its guide and barcode pair."),
        ("Partition cells for 10x 5′ RNA capture", [Row(chunks=[("cell barcode + UMI → transcript, guide/iBAR and PuroR cDNA", "umi", False)])], "Scaffold- and PuroR-specific RT primers are spiked into the reaction."),
        ("Split cDNA into three RNA-derived libraries", [Row(chunks=[("transcriptome | nested guide/iBAR PCR | nested PuroR-barcode PCR", "r2", False)])], "All three products preserve the droplet cell barcode."),
        ("INFERRED — complete compatible sequencing shells", [Row(chunks=[("four source-established payloads → complete Illumina-compatible endpoints", "me", False)])], "The drawn outer adapter completion is inferred; genotype linkage and barcode payloads are source-established."),
    ]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input("cells","cells before guide-library integration",MaterialKind.BIOLOGICAL,primary=True)
    cassette=g.input("guide_cassette","doubly barcoded guide cassette with 6-nt iBAR and 12-nt PuroR barcode",MaterialKind.MOLECULAR,state=sections()[0][1])
    lentiviral=g.input("integration_machinery","lentiviral integration machinery",MaterialKind.CATALYST)
    integrated=g.product("integrated","cells carrying a guide, iBAR and PuroR barcode pair",MaterialKind.MOLECULAR,state=sections()[0][1])
    g.react("integrate","Integrate the doubly barcoded guide cassette",StepKind.INTEGRATION,((cassette,InputRole.TEMPLATE),(cells,InputRole.ANALYTE),(lentiviral,InputRole.CATALYST)),(integrated,))

    guide=g.input("guide_rna","integrated guide RNA",MaterialKind.OLIGO)
    editor=g.input("base_editor","base-editing enzyme",MaterialKind.CATALYST)
    edited=g.product("edited","base-edited cells retaining the barcode pair",MaterialKind.MOLECULAR,state=sections()[1][1])
    g.react("edit","Base-edit the guide-directed genomic locus",StepKind.GENOME_EDITING,((integrated,InputRole.TEMPLATE),(guide,InputRole.STRUCTURAL),(editor,InputRole.CATALYST)),(edited,))
    expanded=g.product("expanded","expanded edited clone descendants",MaterialKind.MOLECULAR,state=sections()[1][1])
    g.react("expand","Expand edited clone descendants",StepKind.GENERIC,((edited,InputRole.TEMPLATE),),(expanded,),changes_molecule=False)

    droplet_condition=g.input("droplet_condition","10x 5-prime droplet partitioning condition",MaterialKind.CONDITION)
    occupied=g.product("occupied","occupied 10x droplets containing edited cells",MaterialKind.MOLECULAR,state=sections()[2][1])
    empty=g.product("droplet_waste","empty or unusable droplets",MaterialKind.WASTE)
    g.react("partition","Partition edited cells for 10x 5-prime capture",StepKind.PARTITION,((expanded,InputRole.TEMPLATE),(droplet_condition,InputRole.CONDITION)),(occupied,empty),changes_molecule=False)
    g.discard(empty,"does not contain a usable edited-cell template")

    bead=g.input("bead_oligo","proprietary 10x 5-prime cell-barcode and UMI bead oligo",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (commercial oligo sequence not printed as a complete endpoint)")
    guide_rt=g.input("guide_rt_primer","source-printed oJR160 guide-scaffold spike-in RT primer",MaterialKind.OLIGO,state=GUIDE_RT)
    puro_rt=g.input("puro_rt_primer","source-printed RT_Puro_Singlecell spike-in RT primer",MaterialKind.OLIGO,state=PURO_RT.upper())
    rt=g.input("reverse_transcriptase","10x 5-prime reverse transcriptase",MaterialKind.CATALYST)
    captured=g.product("captured","droplet cell retaining gDNA and three RNA-derived cDNA classes",MaterialKind.MOLECULAR,state=sections()[2][1],inferred=True,source=SOURCE+" (complete commercial endpoint oligos not printed)")
    g.react("reverse_transcribe","Capture transcript, guide/iBAR and PuroR-barcode cDNA",StepKind.REVERSE_TRANSCRIPTION,((occupied,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(guide_rt,InputRole.STRUCTURAL),(puro_rt,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(captured,))

    transcript=g.product("transcript_template","10x transcriptome cDNA",MaterialKind.MOLECULAR,state=RNA,inferred=True,source=SOURCE+" (outer completion inferred)")
    guide_template=g.product("guide_template","guide and 6-nt iBAR cDNA",MaterialKind.MOLECULAR,state=GUIDE,inferred=True,source=SOURCE+" (outer completion inferred)")
    puro_template=g.product("puro_template","12-nt PuroR-barcode cDNA",MaterialKind.MOLECULAR,state=PURO,inferred=True,source=SOURCE+" (outer completion inferred)")
    genotype=g.product("genotype_template","targeted edited-locus genomic DNA",MaterialKind.MOLECULAR,state=GENO,inferred=True,source=SOURCE+" (outer completion inferred)")
    g.react("separate_readouts","Separate transcript, guide/iBAR, PuroR-barcode and genotype readouts",StepKind.SPLIT,((captured,InputRole.TEMPLATE),),(transcript,guide_template,puro_template,genotype),changes_molecule=False)

    reviewed=["cells","guide_cassette","integration_machinery","integrate","guide_rna","base_editor","edit","expand","droplet_condition","partition","bead_oligo","guide_rt_primer","puro_rt_primer","reverse_transcriptase","reverse_transcribe","separate_readouts"]
    finish_pcr_branch(g,"rna",transcript,RNA,TENX_PRIMERS,reviewed,label="scSNV-seq transcriptome",inferred_primers=True,limitation="The complete commercial 10x 5-prime endpoint primer pair is not printed as one construct.")
    finish_pcr_branch(g,"guide",guide_template,GUIDE,GUIDE_PRIMERS,reviewed,label="scSNV-seq guide and iBAR",inferred_primers=True,geometry=PCRGeometry.regional(GUIDE,NESTED_FWD,reverse_pool="10x 5' HT library-construction primer (not printed)",reverse_region="Index 1 / Read 2 arm",statement="Methods: nested PCR from the small SPRI fraction after cDNA amplification. Supplementary Table S2 prints the forward nested primer (P5 + 10-base library index + TruSeq Read 1); its reverse partner is not printed."),limitation="The complete nested guide/iBAR enrichment and outer primer pair is not printed as one construct.")
    finish_pcr_branch(g,"puro",puro_template,PURO,PURO_PRIMERS,reviewed,label="scSNV-seq PuroR barcode",inferred_primers=True,geometry=PCRGeometry.regional(PURO,NESTED_FWD,reverse_pool="10x 5' HT library-construction primer (not printed)",reverse_region="Index 1 / Read 2 arm",statement="The analogous PuroR nested PCR uses the same printed forward primer (Puro_i5longFwd); its reverse partner is not printed."),limitation="The complete nested PuroR-barcode enrichment and outer primer pair is not printed as one construct.")
    finish_pcr_branch(g,"genotype",genotype,GENO,GENO_PRIMERS,reviewed,label="scSNV-seq targeted genotype",inferred_primers=True,limitation="The target-dependent genotype primer panel is not one universal published endpoint pair.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
