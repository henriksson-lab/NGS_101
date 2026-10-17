# Source manifest

- Defining paper: https://doi.org/10.1073/pnas.1105422108
- Open full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC3111315/

Fetch with `python3 tools/get_sources.py "Safe-SeqS"`.

SI Table S4 (image-only) is read from the original PNAS SI facsimile in a USPTO PTAB
exhibit, fetched by the `extra_sources.tsv` row for this slug as `download-documents.pdf`
(https://ptacts.uspto.gov/ptacts/public-informations/petitions/1549805/download-documents?artifactId=kofo6g92QN3i2XDa55gLrZUDekd75G4UPgP6LPFbcqDOjexY8FWVqaI&format=pdf;
SI on exhibit pages 0007–0016, Table S4 on 0015–0016). Render those pages to read the table.
