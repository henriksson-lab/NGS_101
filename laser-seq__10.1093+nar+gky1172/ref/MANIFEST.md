# Source manifest: LASER-Seq

DOI https://doi.org/10.1093/nar/gky1172 . Originals remain external under `/husky/henriksson/NGS101_papers/sources/laser-seq__10.1093+nar+gky1172`.

Re-obtain with `python3 tools/get_sources.py --slugs laser-seq__10.1093+nar+gky1172`; then `python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/laser-seq__10.1093+nar+gky1172`. The external MANIFEST.tsv maps originals/URLs/MD5s; a filename is not source identity.

Inspected originals:

- `nar_gky1172_PMC6326810.xml`: Published VOR, PMC licence/identity explicitly not author manuscript. Introduction38–46/Fig1; Methods52–80 including treatment58–68 and complete library76, acquisition76 and processing80.
- `nar_gky1172_supp_gky1172_supplemental_files_180423_LASERseq_SuppInfo_Table_S1_Oligos.xlsx`: Original published SI workbook from original supplement ZIP. Sheet1 rows2–5 full modifications and PCR family rows6–29; workbook-derived text inspected, exactRT/linker source transcription guarded.
- `nar_gky1172_supp_gky1172_supplemental_files_180831_LASER_SEQ_SuppInfo_.docx`: Original published SI from same ZIP. Title and FigureS1–S8 captions; no additional preparation recipe inferred from plots.

Pending directly linked originals/source versions:

- MainMethods52/76 references8,22,23,31: Actual directly cited NAz synthesis/library predecessors unreviewed at three-original cap; circle-to-PCR spacer traversal/primer geometry remains unassembled.
- MainMethods64 ZymoDirectZol original manual: Historical edition unreviewed; K562 omits DNase per main, no contemporary default substituted.

No third-party source material is distributed here. Source-dependent tests skip when originals are absent. Complete protocol/source/absence review is not claimed.
