#!/usr/bin/env python3
"""Source transcription checks for SUM-seq's published oligo table."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from storage import data_root

import chem as P
from checks import Check, Source, have


check = Check()
check("sample and droplet identity parts", P.CELL_PARTS,
      (("sample index", 11), ("droplet barcode", 16)))

source_dir = data_root() / "sources" / HERE.parent.name
table = Source(source_dir / "s41592-025-02700-8_supp_41592_2025_2700_MOESM3_ESM.xlsx.txt",
               "Llorens-Rico et al. 2025 Supplementary Table 1a–d, MOESM3 workbook; "
               "fetch with tools/get_sources.py and convert with tools/doctext.py")
if have(check, table, label="SUM-seq Supplementary Table 1 oligo transcription"):
    rows = {}
    index_rows = {}
    with table.path.open(encoding="utf-8") as source:
        for line in source:
            if line.startswith("## sheet Supp. Table 2"):
                break
            cells = line.rstrip("\n").split("\t")
            if len(cells) > 1:
                rows[cells[0]] = cells[1].strip().replace("\u00a0", "")
                if cells[0].startswith(("SUMseq-ATAC-", "SUMseq-RNA-")) and len(cells) > 2:
                    index_rows[cells[0]] = cells[2]

    for source_name, model_sequence in (
        ("SUMseq-ATAC-1", P.ATAC_INDEX_1),
        ("SUMseq-RNA-1", P.RNA_INDEX_1),
        ("Mosaic-end oligo (ME-ddC)", P.TN5_ME_DDC),
        ("Blocking invdT oligo", P.TN5_BLOCK),
        ("Nextera Read 1 (R1N)", P.NEXTERA_READ_1),
        ("Partial-P5-PTO", P.PARTIAL_P5_PTO),
        ("TruseqR2", P.TRUSEQ_R2),
        ("ATAC-spacer", P.ATAC_SPACER),
        ("P7-i7-spacer", P.P7_I7_SPACER),
        ("Partial-P7-PTO", P.PARTIAL_P7_PTO),
        ("P7_i7_TruSeqR2", P.P7_I7_TRUSEQ_R2),
    ):
        check(f"published {source_name}", rows.get(source_name), model_sequence)

    families = []
    for name in ("SUMseq-ATAC-", "SUMseq-RNA-"):
        indices = [index for key, index in index_rows.items() if key.startswith(name)]
        families.append((name, len(indices), len(set(indices)), all(len(index) == 11 for index in indices)))
    check("both published pre-index families are 96 distinct 11-base members", families,
          [("SUMseq-ATAC-", 96, 96, True), ("SUMseq-RNA-", 96, 96, True)])

check.report()
