#!/usr/bin/env python3
"""Guard the xpress endpoint against its original author-linked order workbook."""
from pathlib import Path
import csv
import os
import build_page as model
from checks import Check, Source, have
from chemdraw import Construct, Segment

check = Check()
data = Path(os.environ.get("CHEM_DATA", model.HERE.parents[1] / "_data"))
source = Source(data / "sources" / model.HERE.parent.name / "gh7u2zw6.xls.txt",
                "Original protocol DOI yxmvmk1yng3p/v1; 384set_Nextera_BC_order.xls attachment https://content.protocols.io/files/gh7u2zw6.xls")
if have(check, source, label="Original custom index order workbook"):
    sheets = {}
    for line in source.path.read_text().splitlines():
        if line.startswith("## sheet "):
            sheet = line.removeprefix("## sheet ")
            sheets[sheet] = {}
        else:
            fields = line.split("\t")
            if len(fields) == 3 and fields[0] != "Well Position":
                sheets[sheet][fields[0]] = fields[1:]
    check("complete ordered i5 A1 oligo", model.INDEX_FORWARD,
          sheets["UD_i5_indices"]["A1"][1])
    check("complete ordered i7 A1 oligo", model.INDEX_REVERSE,
          sheets["UD_i7_indices"]["A1"][1])
    for title, prefix, suffix in (("UD_i5_indices", model.il.P5, model.nx.S5),
                                  ("UD_i7_indices", model.il.P7, model.nx.S7)):
        rows = list(sheets[title].values())
        check(title + " all order rows agree with local source architecture",
              len(rows) == 384 and all(len(index) == len(model.I5_A1)
                  and sequence == prefix + index + suffix for index, sequence in rows))

library = model.final_library()
shortened = Construct([Segment(s.name, s.top[:8], s.tag, feature=s.feature)
                       if s.name == "illustrative i5 A1" else s for s in library.segments])
check.raises("full ordered primers reject inherited eight-base index endpoint",
             lambda: model.require_index_sites(shortened))
check.raises("amplified product cannot omit the forward-primer overhang",
             lambda: model.require_preamp_sites(model.switched_cdna()))
annotation = Source(data / "sources" / model.HERE.parent.name /
                    "Smartseq3xpress_SeqAmp_TSOs_RTtest_v6.barcode_annotation.txt",
                    "Original E-MTAB-11467 source annotation, named in its SDRF and BioStudies file index")
if have(check, annotation, label="Original TSO design annotation"):
    with annotation.path.open() as handle:
        row = next(r for r in csv.DictReader(handle, delimiter="\t")
                   if r["TSO"] == "M462_TSO")
    check("source-labelled identifier includes both eight-N region and WW",
          "".join(s.top for s in library if s.feature == model.UMI), row["UMIsequence"])
    spans = model.sp.feature_spans(library, model.SEQ_PRIMERS, {"Read 1": 24})
    check("derived read spans preserve the source identifier length",
          sum(s.cycle_end - s.cycle_start + 1 for s in spans if s.feature == model.UMI),
          int(row["UMIlength"]))
check.report()
