"""scSNV-seq transcript, barcode and genotype branches."""
from chemdraw import Row
from multibranch import illumina_branch
from multimodal_spatial import finish_pcr_branch
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind

TITLE = "scSNV-seq — genotype-linked single-cell transcriptomes"
NOTES = "01_scsnv-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1186/s13059-024-03169-y">Cooper et al., <i>Genome Biology</i> (2024)</a>.'
SUMMARY = "A guide-associated 6-nt iBAR and an independent 12-nt PuroR barcode connect targeted DNA genotypes to 10x 5′ transcriptomes."
CAVEAT = "The paper establishes the four readout branches, barcode lengths and 10x 5′ capture. The exact complete endpoint oligos are not printed as one final construct; their standard-compatible sequencing shells are therefore marked inferred."
RNA, TENX_PRIMERS = illumina_branch("scSNV-seq transcriptome library", "10x 5-prime cDNA", cell=16, umi=10, inferred_adapters=True)
GUIDE, GUIDE_PRIMERS = illumina_branch("scSNV-seq guide/iBAR library", "guide plus 6-nt iBAR", cell=16, umi=10, guide=26, inferred_adapters=True)
PURO, PURO_PRIMERS = illumina_branch("scSNV-seq PuroR-barcode library", "12-nt PuroR barcode", cell=16, umi=10, target=12, inferred_adapters=True)
GENO, GENO_PRIMERS = illumina_branch("scSNV-seq targeted genotype library", "edited genomic locus", cell=16, inferred_adapters=True)
FINAL_LIBRARIES = (
    ("Transcriptome library", RNA, TENX_PRIMERS, "10x gene-expression reads carry cell barcode and UMI.", "The transcript branch uses the declared 10x-compatible run primers."),
    ("Guide and iBAR library", GUIDE, GUIDE_PRIMERS, "A nested enrichment reads the guide together with its 6-nt iBAR.", "The guide branch retains the same droplet cell identity."),
    ("PuroR-barcode library", PURO, PURO_PRIMERS, "A separate nested enrichment reads the 12-nt transcribed PuroR barcode.", "The PuroR branch retains the same droplet cell identity."),
    ("Targeted genotype library", GENO, GENO_PRIMERS, "Targeted single-cell DNA amplicons call the edited allele.", "The DNA branch supplies the genotype joined computationally through the barcode pair."),
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
    guide_rt=g.input("guide_rt_primer","guide-scaffold spike-in RT primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete endpoint oligo not printed)")
    puro_rt=g.input("puro_rt_primer","PuroR-barcode spike-in RT primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete endpoint oligo not printed)")
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
    finish_pcr_branch(g,"guide",guide_template,GUIDE,GUIDE_PRIMERS,reviewed,label="scSNV-seq guide and iBAR",inferred_primers=True,limitation="The complete nested guide/iBAR enrichment and outer primer pair is not printed as one construct.")
    finish_pcr_branch(g,"puro",puro_template,PURO,PURO_PRIMERS,reviewed,label="scSNV-seq PuroR barcode",inferred_primers=True,limitation="The complete nested PuroR-barcode enrichment and outer primer pair is not printed as one construct.")
    finish_pcr_branch(g,"genotype",genotype,GENO,GENO_PRIMERS,reviewed,label="scSNV-seq targeted genotype",inferred_primers=True,limitation="The target-dependent genotype primer panel is not one universal published endpoint pair.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
