# NEBNext Ultra II USER-adaptor workflow — third-party reference material

Third-party files are deliberately not committed. Re-obtain the sources from:

- [E7103/E7645 Ultra II DNA manual v7.0_9/22](https://www.neb.com/en-us/-/media/nebus/files/manuals/manuale7103-e7645.pdf).
- [E7600 Dual Index Primers Set 1 manual v6.0_6/24](https://www.neb.com/en/-/media/nebus/files/manuals/manuale7600.pdf).
- [E7780 Dual Index Primers Set 2 manual](https://www.neb.com/en-us/-/media/nebus/files/manuals/manuale7780.pdf), p. 3 independently prints the same modified hairpin.
- [Original third-party adaptor FAQ](https://www.neb.com/faqs/can-i-use-this-nebnext-kit-with-adaptors-and-primers-from-other-vendors-than-neb). Dated and current canonical URLs return403 to the downloader in this environment; the original official page is readable through browser retrieval.
- [FAQ-linked E7395 DNA UDI/UMI manual](https://www.neb.com/en/-/media/nebus/files/manuals/manuale7395.pdf), v2.0_2/21, separate full-length indexed-adaptor branch.
- [FAQ-linked E7416 RNA UDI/UMI manual](https://www.neb.com/en/-/media/nebus/files/manuals/manuale7416.pdf), v5.0_9/25, separate RNA adaptor branch.
- [NEB USER Enzyme documentation](https://www.neb.com/products/m5505-user-enzyme).
- [Illumina Adapter Sequences](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/CDIndexes.htm), document 1000000002694 v22.

Run `python3 tools/get_sources.py "NEBNext Ultra II DNA Library Prep (USER adaptor)"`.
Files belong under `$CHEM_DATA/sources/nebnext-ultra-ii-dna-library-prep-user-adaptor/`,
never here.

Local transcription checks read original E7600/E7780 text twins, and skip when those originals are absent. They cover the complete modified hairpin and the selected i501/i701 PCR oligos.
