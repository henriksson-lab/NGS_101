# Source manifest

- Defining article and Supplementary Table S2: https://doi.org/10.1186/s13059-016-0938-8
- Detailed author protocol: Additional file 4 of that article. It specifies Illumina RP1 and RPIX for library PCR but refers their sequences to Illumina.
- Original vendor oligo table: https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-SmallRNA.htm (Illumina Adapter Sequences, document 1000000002694 v22). It prints RP1, the RPI family and RPI1's index in separate rows.
- Retrieve originals with `python3 tools/get_sources.py --slugs cel-seq2__10.1186+s13059-016-0938-8`. Files and extracted text belong in `$CHEM_DATA/sources/`, not this directory.
