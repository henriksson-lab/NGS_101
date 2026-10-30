# Source retrieval

Use `python3 tools/get_sources.py --slugs spatial-vdj__10.1126+science.adf8486`. Sources and text conversions belong in the external source store, never this checkout. Exact URLs are also declared in `catalogue/extra_sources.tsv`.

## Retrieved author manuscript

[Zenodo 13148082](https://zenodo.org/records/13148082) explicitly identifies the deposited PDF as the **accepted author version**, with a separate definitive Science publication. Downloaded and revalidated on 30 October 2026 as `zenodo13148082_adf8486_accepted_reduced.pdf` (24,024,126 bytes); its extracted text is beside it in the external store. It includes consolidated supplementary Methods, but lists Data S1 (primer and probe sequences) as a separate file. This is not a preprint or evidence of equivalence to the final publication. See the external `MANIFEST.tsv` for the exact download URL and checksum.

## Final published article and user-supplied supplements retrieved

The final full-access article was saved through the authorized visible browser as `science_adf8486_final_browser.html` (943,856bytes), validated and converted. Its final Methods1623–1669 and Results1580 support the selected SR-TRB preparation.

The user subsequently supplied the actual published downloads in `/tmp/camilla`. Both were moved with checksum verification to the external source store, source-manifested and revalidated:

- [Supplementary PDF](https://www.science.org/doi/suppl/10.1126/science.adf8486/suppl_file/science.adf8486_sm.pdf):`science.adf8486_sm.pdf`,6,201,464bytes. Title/DOI and Methods1–72 inspected, particularly distinct nonspatial validation and enrichment-test branches.
- [Data S1–S5 ZIP](https://www.science.org/doi/suppl/10.1126/science.adf8486/suppl_file/science.adf8486_data_s1_to_s5.zip):126,270bytes, intact ZIP verified. It contains `science.adf8486_data_s1.xlsx` through `_data_s5.xlsx`; all extracted originals and text twins remain external. Data S1 spatial PartRead1/TRBV/probe/blocker worksheets were inspected. Data S2–S5 are available but uninspected, not chemistry authorities in this pass.

The external `MANIFEST.tsv` records original publisher URLs, archive-member identity, sizes and MD5 checksums. No further download is requested for these successfully obtained files. Earlier browser verification failures are retrieval history, not their current availability. Historical vendor-source questions remain separate from obtained files.

Related 2022 preprint entries in the external manifest remain superseded and excluded. The 25.6 GB Zenodo data/scripts archive was not downloaded or substituted for chemistry sources. No fingerprint spoofing, headless challenge solving or proxy rotation was used; automated fetcher cool-downs were preserved.

See [selected scope and follow-up](../01_spatial_vdj.md).
