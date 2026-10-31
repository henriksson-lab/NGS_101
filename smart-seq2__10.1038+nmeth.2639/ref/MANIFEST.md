# SMART-seq2 reference originals

Third-party originals and text twins remain outside the checkout under
`/husky/henriksson/NGS101_papers/sources/smart-seq2__10.1038+nmeth.2639/`.
No paper, vendor PDF, converted text or rendered source page is checked in here.
The external `MANIFEST.tsv` records original URLs, bytes, checksums and source identity.

## Re-obtain

Both source-era vendor originals are registered in `catalogue/extra_sources.tsv`.
Run from the repository root:

```sh
CHEM_DATA=/husky/henriksson/NGS101_papers python3 tools/get_sources.py --slugs smart-seq2__10.1038+nmeth.2639
CHEM_DATA=/husky/henriksson/NGS101_papers python3 tools/get_sources.py --revalidate --slugs smart-seq2__10.1038+nmeth.2639
python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/smart-seq2__10.1038+nmeth.2639/
```

| Original URL | External filename | Identity and inspected source scope |
|---|---|---|
| [Picelli2013 article](https://link.springer.com/content/pdf/10.1038/nmeth.2639.pdf) | `nmeth.2639.pdf` | Published version of record; Online Methods RT/ISPCR and named Nextera DNA back end. See main research note for existing inspected locations. |
| [Picelli2013 Supplementary Table5](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.2639/MediaObjects/41592_2013_BFnmeth2639_MOESM304_ESM.xlsx) | `nmeth.2639_41592_2013_BFnmeth2639_MOESM304_ESM.xlsx` | Text19–37, especially27–29: actual source-paper kit catalogue IDs. |
| [Illumina Nextera DNA dual-indexing insert](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexteradna/15029964_C_Nextera_DNA_Sample_Preparation_Insert.pdf) | `illumina_15029964_C_Nextera_DNA_Sample_Preparation_Insert.pdf` | Printed part15029964 RevC, copyright2012,2pp. Both pages read; text34–68 kit identities and123–136 standard index programme. Renderedp2 inspected. |
| [Illumina customer sequence letter, public vendor-original mirror](https://www.seqanswers.com/filedata/fetch?id=304385) | `source_pdf_fbf4f1d1ac38acc3ee4b.pdf` | Printed24May2012,16pp; filename/metadataMay23 is not the printed date. Printedpp1–3 read. Text86–183 contains kit IDs, transferred adapter/PCR orders and adapter-versus-sample-sheet index tables. Renderedpp1/3 inspected to verify identity and table orientation. |

The opaque PDF URL is registered with the `pdf:` handler; its query parameter is part
of source identity. Do not rename it to a generic cached `fetch.pdf` or replace it with
a response for another attachment ID. The two vendor originals were fetched, validated
and converted by the source tools on31October2026. The2012 letter's components/widths
are now positive source support for the named kit family, not modern-family inference.
They do not identify the authors' actual study-well index choices or disclose all
reagent modifications, NPM/PPC contents or the exact full preparation-guide revision.

Individual source-letter sequence reproduction carries:
**Oligonucleotide sequences ©2007–2012 Illumina, Inc. All rights reserved.** The whole
letter remains external and is not republished by this project.
