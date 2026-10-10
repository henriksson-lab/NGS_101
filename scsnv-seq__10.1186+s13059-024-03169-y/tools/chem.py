"""scSNV-seq transcript, barcode and genotype branches."""
from chemdraw import Row
from multibranch import illumina_branch

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
