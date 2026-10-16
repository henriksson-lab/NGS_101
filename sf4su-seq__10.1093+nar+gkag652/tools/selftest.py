#!/usr/bin/env python3
"""Source oligos, anti-dimer handling and sf4sU-specific boundary checks."""
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from storage import data_root
from checks import Check, Source, have
from protocol_graph import Disposition
import chem

check = Check()
base = data_root() / "sources"
bundle = base / HERE.parent.name
manual = Source(bundle / "manuale7300_e7330_e7560_e7580.pdf.txt", "Original NEB E7300 manual Version9.0_07/25")
if have(check, manual, label="Dated original E7300 guide"):
    text = manual.path.read_text()
    check("guide activated/blocked3′ adapter order retained",
          "rApp" + chem.ADAPTER_3P + "-NH2" in text and "rApp" + chem.ADAPTER_3P + "-NH2" in chem.kit_oligos()[0])
    check("guide whole RNA5′ adapter remains source-transcribed",
          "".join("r" + base for base in chem.ADAPTER_5P) in text)
    check("guide whole RT and protected PCR oligos agree with the displayed reference",
          chem.RT_PRIMER in text and chem.PCR_FORWARD[:-1] + "*" + chem.PCR_FORWARD[-1] in text and
          chem.PCR_REVERSE_3P[:-1] + "*" + chem.PCR_REVERSE_3P[-1] in text)
    steps = list(chem.protocol_graph().steps)
    check("source anti-dimer RT-primer hybridization remains before5′ ligation",
          "Hybridize the Reverse Transcription Primer" in text and
          steps.index("ligate_3p") < steps.index("block_adapter") < steps.index("ligate_5p"))
methods = Source(bundle / "PMC13309788_supp_gkag652_supplemental_file.pdf.txt", "Original2026 supplementary capture/sf4sU Methods and Figure4C")
if have(check, methods, label="Original sf4sU source"):
    text = methods.path.read_text()
    graph = chem.protocol_graph()
    check("source-observed adapter dimer has an explicit unresolved physical fate",
          "second-highest bands are a result of adaptor dimerisation" in text and
          graph.terminals["adapter_dimer"].disposition is Disposition.UNRESOLVED)
    check("source DTT release and decapping recovery survive graph assembly",
          "dithiothreitol (DTT)" in text and "heat-inactivated" in text and
          all(step in graph.steps for step in ("release", "stop_decap", "decap_recovery")))
    notice = next(a for a in graph.annotations if a.title == "Published wash-buffer difference")
    check("actual2026 wash is preserved as a source value without a corrected model concentration",
          "1 mM NaCl" in text and "1 mM NaCl" in notice.claims[0].claim and
          "unassigned" in notice.model_choice)
archive = Source(bundle / "efetch.fcgi", "Original NCBI SRA experiment and run metadata")
if have(check, archive, label="Deposited size-selected library metadata"):
    root = ET.parse(archive.path).getroot()
    runs = root.findall(".//RUN")
    check("sample-specific single genomic read length follows all four original run records",
          root.find(".//EXPERIMENT").get("accession") == "SRX11439177" and
          root.find(".//LIBRARY_LAYOUT/SINGLE") is not None and len(runs) == 4 and
          all(run.find("./Statistics/Read[@index='0']").get("average") == str(chem.SPEC.READ_LENGTHS["Read 1"]) and
              run.find("./Statistics/Read[@index='0']").get("stdev") == "0" and
              run.find("./Statistics/Read[@index='0']").get("count") == run.get("total_spots") and
              run.find("./Statistics/Read[@index='1']").get("count") == "0" for run in runs))
sample = Source(bundle / "geo-GSE179306-gsm-full.txt", "Original GEO full sample metadata")
series = Source(base / "snu-seq__10.1093+nar+gkag652" / "geo-GSE179306-self-full.txt", "Original GEO series metadata")
if have(check, sample, series, label="Original size-window records"):
    sample_body = sample.path.read_text().split("^SAMPLE = GSM5452295", 1)[1].split("^SAMPLE =", 1)[0]
    notice = next(a for a in chem.protocol_graph().annotations if a.title == "Size-selection windows differ across source records")
    check("conflicting deposited size windows stay visible without selecting a universal interval",
          "70-100" in sample_body and "50-100" in series.path.read_text() and "unassigned" in notice.model_choice)
check.report()
