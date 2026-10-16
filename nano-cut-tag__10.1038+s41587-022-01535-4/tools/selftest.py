#!/usr/bin/env python3
"""Transcription checks against the original publisher workbook and kit guide."""
from pathlib import Path
import re
import sys
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1] / "lib")]
from storage import data_root
import nano_cut_tag as P
from checks import Check, Source, have

c = Check()
c.section("nano-CT original oligo transcription")
data = data_root() / "sources"
workbook = Source(data / H.parent.name /
    "s41587-022-01535-4_supp_41587_2022_1535_MOESM3_ESM.xlsx.txt",
    "doi:10.1038/s41587-022-01535-4, Supplementary Tables 1–4 workbook; tools/get_sources.py")
if have(c, workbook):
    rows = {}
    for line in workbook.path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) >= 2:
            rows[fields[0]] = fields[1]
    for name, barcode in (("A", P.BARCODE_A), ("B", P.BARCODE_B), ("C", P.BARCODE_C)):
        for stagger in range(4):
            key = f"Tn5_P5_MeA_Bcd{name}_{stagger}N"
            c(key, "".join(s.top for s in P.mea_parts(barcode, name, stagger)), rows[key])
    c("MeB", "".join(s.top for s in P.meb_parts()), rows["Tn5_P7_MeB"])
    c("Tn5_Rev with source phosphorylation", "5Phos/" + P.ME_REV, rows["Tn5_Rev"])
    c("custom R1", P.CUSTOM_R1, rows["R1_seq"])
    c("custom I2", P.CUSTOM_I2, rows["I2_seq"])
guide = Source(data / "scifi-atac-seq__10.1101+2023.09.17.558155" /
    "10x_CG000209_Chromium_NextGEM_SingleCell_ATAC_ReagentKits_v1.1_UserGuide_RevG.pdf.txt",
    "10x CG000209 Rev G Appendix A1; source guide registered for scifi-atac-seq")
if have(c, guide):
    text = guide.path.read_text().replace(" ", "").replace("-", "")
    c("source forward SI-PCR scaffold", P.si_pcr_forward().top() in text)
    c("source reverse SI-PCR scaffold", P.si_pcr_reverse().top() in text)
config = Source(data / H.parent.name / "gh_config.yaml",
    "Author study config/config.yaml at commit 7657a3a9f15e798b0a8d2ba01ab7f3006a238d21; catalogue/extra_sources.tsv")
if have(c, config):
    mappings = dict(re.findall(r"^\s+(ATAC|H3K27ac|H3K27me3): ([ACGT]+)\s*$",
                               config.path.read_text(), re.MULTILINE))
    c("study modality barcode assignments",
      (P.BARCODE_A, P.BARCODE_B, P.BARCODE_C),
      tuple(mappings[k] for k in ("ATAC", "H3K27ac", "H3K27me3")))
    c("representative antibody-targeted barcode matches study H3K27ac",
      next(s.top for s in P.mea_parts() if s.feature), mappings["H3K27ac"])
P.workflow()
P.protocol_graph()
sys.exit(c.report())
