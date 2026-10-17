# Source manifest

Run `python3 tools/get_sources.py "Spatial C&T"`. The defining paper, original hsrChST-seq
preprint and detailed Nature Protocols article are linked in the research note. No
third-party source files are committed.

The defining Science Supplementary Materials DOCX is behind the PMC download gate: fetch
`https://pmc.ncbi.nlm.nih.gov/articles/instance/7612972/bin/EMS145971-supplement-Supplementary_Material.docx`
in a browser and save it as
`science.abg7216_supp_EMS145971-supplement-Supplementary_Material.docx` in
`$CHEM_DATA/sources/spatial-cut-tag__10.1126+science.abg7216/` (done 2026-10-17), then run
`python3 tools/doctext.py` on it.
