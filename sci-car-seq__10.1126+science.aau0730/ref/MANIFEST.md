# Source manifest

- Cao et al. 2018, DOI `10.1126/science.aau0730`, PMC6571013; obtain with
  `python3 tools/get_sources.py --doi 10.1126/science.aau0730`.
- Cao 2018 Supplementary Materials PDF and supplementary table workbook (Tables S1-S13):
  PMC instance URLs
  `https://pmc.ncbi.nlm.nih.gov/articles/instance/6571013/bin/NIHMS1033869-supplement-1.pdf`
  and `.../NIHMS1033869-supplement-Supplimentary_table.xlsx`. The PMC download gate blocks
  scripted retrieval; download by hand in a browser and place them in the source bundle
  (`$CHEM_DATA/sources/sci-car-seq__10.1126+science.aau0730/`), then run
  `get_sources.py --revalidate`.
