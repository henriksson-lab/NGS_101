#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from storage import data_root
import chem
import seqprimers as sp
from checks import Check, Source, have

check = Check()
source_dir = data_root() / "sources" / HERE.parent.name
article = Source(source_dir / "gad.350434.123_PMC10270195.xml.txt",
                 "https://pmc.ncbi.nlm.nih.gov/articles/PMC10270195/")
workbook = Source(source_dir / "gad.350434.123_supp_supp_gad.350434.123_Supplemental_Primers.xlsx.txt",
                  "PMC10270195 Supplemental_Primers.xlsx")
if have(check, article, workbook, label="Butt-seq source oligo transcription"):
    methods = article.path.read_text(encoding="utf-8")
    table = workbook.path.read_text(encoding="utf-8")
    check("SCR2 blocked RNA is printed in article Methods", chem.SCR2_RNA in methods)
    check("SCR2 blocked RNA is printed in supplemental workbook", chem.SCR2_RNA in table)
    check("full-length R1MER blocked RNA is printed in workbook", chem.R1ME_RNA in table)
    for name, pool in (("SCR2_DNA", chem.SCR2_DNA_POOL), ("R1ME", chem.R1ME_DNA_POOL)):
        for base, sequence in zip("ATCG", pool):
            check(f"{name}_{base} occurs as a whole workbook row",
                  f"{name}_{base}\t{sequence}\t" in table)
    check("SCP7 universal PCR oligo is printed in workbook", f"SCP7\t{chem.SCP7}\t" in table)
    check("representative Ad1.60 index PCR oligo is printed in workbook",
          f"v2_Ad1.60_TACTGACC\t{chem.AD1_60}\t" in table)
    check("article's shorter R1ME DNA is distinct from workbook sequence",
          "GCAGCGTCAGATGTGTATAAGAGACAGN" in methods
          and "GCAGCGTCAGATGTGTATAAGAGACAGN" not in table)
chem.protocol_graph()
spans = sp.feature_spans(chem.FINAL_LIBRARY, chem.SEQ_PRIMERS, chem.READ_LENGTHS)
check("eight-cycle UMI and i5 sample index land in different index acquisitions",
      [(span.read, span.feature.role, span.cycles) for span in spans],
      [("Index 1 (i7)", "umi", "1–8"),
       ("Index 2 (i5)", "sample_index", "1–8")])
read2 = sp.locate(chem.FINAL_LIBRARY, sp.TRUSEQ["R2"])
check("Read 2 begins with the first-starter fixed T before the variable overhang",
      (read2.strand, read2.reads[:2]), ("top", "Tn"))
check.report()
