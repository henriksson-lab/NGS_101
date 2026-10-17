# Telomere Profiling — third-party reference material

No source files are committed. Re-obtain Karimian et al. 2024 from
[doi:10.1126/science.ado0431](https://doi.org/10.1126/science.ado0431) (author manuscript
PMC12425158). Its Supplementary Material PDF (Materials and Methods, figs. S1–S6, tables S1
and S3) must be fetched by hand:
<https://pmc.ncbi.nlm.nih.gov/articles/instance/12425158/bin/NIHMS2109074-supplement-Supplemental_Material.pdf>
(PMC download gate). The bioRxiv v2 preprint
([doi:10.1101/2023.12.21.572870](https://doi.org/10.1101/2023.12.21.572870)) is fetched by
`tools/get_sources.py`. Fetched copies belong under `$CHEM_DATA/sources/`.

Round 2026-10-17-deferred-followup-02 found no other legitimate route. The PMC manuscript PDF
<https://pmc.ncbi.nlm.nih.gov/articles/PMC12425158/pdf/nihms-2109074.pdf> is also gated, and Europe PMC
and <https://www.science.org/doi/suppl/10.1126/science.ado0431/suppl_file/science.ado0431_sm.pdf>
return 403. Dryad `doi:10.5061/dryad.dz08kps5d` holds data only.

Round 2026-10-17-deferred-followup-03: the maintainer obtained the published Supplementary
Materials by hand on the institutional network from the Science URL above
(`science.ado0431_sm.pdf`, md5 ff94c1ae81f5e1195f35a5fa39be705b; MANIFEST.tsv row). It is
the source for the model. It contains Materials and Methods, figs. S1–S6 and Tables S1–S3,
with the TeloTag sequences in Table S3. Re-obtain it the same way and run `tools/doctext.py`
on it.
