# Source manifest

No third-party source files are committed here. Re-obtain the article, JATS XML,
PDF and supplements with:

```sh
python3 tools/get_sources.py --doi 10.1016/j.mcpro.2021.100168
```

The transcription uses PMC8633680, especially Figure 1 and the PCR1, PCR2,
sequencing and data-processing sections.

Public metadata/documentation consulted to look for the undisclosed sequence dictionary:

- https://patents.google.com/patent/WO2021191442A1/en
- https://olink.com/software/download
- https://7074596.fs1.hubspotusercontent-na1.net/hubfs/7074596/01-User%20Manuals%20for%20website/1344-Olink-explore-preprocessing.pdf
- https://github.com/crukci-bioinformatics/combine-olink

These were inspected as public documentation, not copied into the repository. The patent
publishes its P5 and P7 adapter portions as SEQ ID 1 and 2, but does not disclose the
complete primers or barcode dictionary. The other sources describe proprietary
preprocessing and symbolic forward/reverse barcode identifiers without the corresponding
bases. The Olink download page says the software requires a purchased licence and
restricts relaying intermediate BCL/FASTQ data; no software or panel-data archive is
vendored here.
