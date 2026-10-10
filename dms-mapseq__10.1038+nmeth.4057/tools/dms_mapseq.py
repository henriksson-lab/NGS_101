"""Zubradt et al. genome-wide DMS-MaPseq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Scene, Segment, circle_rows, revcomp
from circular import circularize_ssdna
from rna_special import adapter_ligation_scene, mutational_read_rows

TITLE = "DMS-MaPseq — mutations report RNA structure"
NOTES = "01_dms-mapseq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.4057">Zubradt et al., <i>Nature Methods</i> (2017)</a>.'
SUMMARY = "DMS modifies exposed A and C bases in living cells; TGIRT reads through those lesions and records them as substitutions before cDNA circularization and custom-primer Illumina sequencing."
CAVEAT = "This page follows the paper's genome-wide branch. Its targeted RT-PCR and UMI-targeted branches are related protocols with different library conversion."
LINKER = "CACTCGGGCACCAAGGA"
RT_HANDLE = "GATCGTCGGACTGTAGAACTCTGAACCTGTCG"
READ1 = "GCAGCGACAGGTTCAGAGTTCTACAGTCCGACGATC"
CDNA = Construct([Segment("RT-primer handle", RT_HANDLE, "r1"),
                  Segment("mutation-bearing cDNA", "X"*42, placeholder=True)],
                 name="DMS-MaPseq cDNA")
CIRCLE = circularize_ssdna(CDNA, five_prime_phosphate=True)
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = "Custom Read 1 primer oNTI202 anneals to the PCR-added handle and reads into the mutation-bearing cDNA; the paper used 50-nt single-end HiSeq 4000 reads."

def final_primer_scene():
    template = [Segment("oNTI202 binding site", revcomp(READ1), "r1"),
                Segment("mutation-bearing cDNA", "X"*42, placeholder=True)]
    primer = [Segment("oNTI202", READ1, "r1")]
    sc=Scene(); sc.strand("library",template,label="PCR-completed library")
    sc.anneal("Read 1 primer",primer,to="library",pair=("oNTI202","oNTI202 binding site"),label="oNTI202")
    sc.arrow("Read 1 primer","Read 1 enters cDNA")
    sc.labels("library")
    return sc

def sections():
    return [
        ("Modify exposed A and C bases with DMS", mutational_read_rows(),
         "DMS reacts preferentially at unpaired Watson-Crick faces; the genome-wide treatment averages about one modification per 50 nt."),
        ("Fragment RNA, dephosphorylate and ligate the 3-prime linker", adapter_ligation_scene(adapter_name="5-prime-rApp linker-2").rows(),
         "Zn2+ fragmentation is followed by size selection and rSAP treatment; truncated T4 RNA ligase 2 joins the preadenylated, 3-prime-blocked linker."),
        ("Reverse-transcribe through lesions with TGIRT", [*mutational_read_rows(), Row(chunks=[("TGIRT cDNA carries substitutions at DMS adducts", "r1", False)])],
         "The long 57-degree reaction favors processive read-through so multiple modifications can be encoded on one cDNA."),
        ("Remove RNA and circularize the phosphorylated cDNA", circle_rows(CDNA,"CircLigase closure"),
         "The RT primer supplies the 5-prime phosphate required for CircLigase; the cDNA becomes a single-stranded circle."),
        ("PCR-add Illumina arms and bind the custom sequencing primer", final_primer_scene().rows(),
         "Nine to thirteen PCR cycles add the indexed flow-cell arms. oNTI202—not a stock TruSeq Read 1 primer—starts the single-end read."),
    ]
