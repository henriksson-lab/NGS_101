# Sources

- PacBio procedure: <https://www.pacb.com/wp-content/uploads/Procedure-checklist-Preparing-Kinnex-libraries-using-the-Kinnex-full-length-RNA-kit.pdf>
- PacBio Kinnex overview: <https://www.pacb.com/technology/kinnex/>
- PacBio SMRT Link v13.1 user guide (MAS8 segmentation reference and Iso-Seq v2 primer FASTA example): <https://www.pacb.com/wp-content/uploads/SMRT-Link-User-Guide-v13.1.pdf>
- PacBio Kinnex full-length RNA public dataset, including `REF-MAS_adapters/MAS-Seq_Adapter_v3/mas8_primers.fasta`, adjacent BarcodeSet XML and the README: <https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/>
- PacBio Kinnex barcoded adapter index FASTA archive, linked from the official multiplexing resource page: <https://www.pacb.com/wp-content/uploads/Kinnex_barcoded_adapter_indexes.fasta_.zip> and <https://www.pacb.com/multiplexing/>
- PacBio Kinnex full-length RNA application note, 102-326-591 Rev03, linked from the dataset README: <https://www.pacb.com/wp-content/uploads/Application-note-Kinnex-full-length-RNA-kit-for-isoform-sequencing.pdf>

No third-party source files are committed.


The application note explicitly cites the MAS-Seq research predecessor:

- [Original article](https://doi.org/10.1038/s41587-023-01815-7), available as PMC12236355 via the registered `pmccloud:` source route.
- [Original Supplementary Information](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41587-023-01815-7/MediaObjects/41587_2023_1815_MOESM1_ESM.pdf), Table 1 PDF pp. 25–26, and [Reporting Summary](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41587-023-01815-7/MediaObjects/41587_2023_1815_MOESM2_ESM.pdf).
- [Longbow original README.rst, pinned snapshot](https://raw.githubusercontent.com/broadinstitute/longbow/248a5ce9ac4c0d2bc803f2b7edcec208a8f0828c/README.rst).
- [MAS-Seq paper-data original README, pinned snapshot](https://github.com/broadinstitute/mas-seq-paper-data/raw/ae9c079f8d9a54f0ae268f3c2edacc7a515f7bbb/README.md).

These are research-method and analysis sources. Their physical primers and software versions are not silently assigned to the commercial eight-member Kinnex kit. The SMRT Link v13.1 primer FASTA is explicitly for read trimming, including a reverse-complemented, poly(T)-excluded 3′ reference.


The MAS reporting summary names exact historical analysis versions. Registered source
rows retrieve the bounded originals used for source-scope checks:

- Longbow `v0.6.4` (`9844452445cc803f814613e203eaf971814babbc`): original `README.rst`, `src/longbow/utils/model.py`, `model_utils.py`, and `wdl/MasSeqSingleFlowcell.wdl`.
- long-read-pipelines `82bfe460778e56ac642942e62638d594f22bd6bf`: `wdl/PB10xMasSeqArraySingleFlowcellv5.wdl`, `docker/lr-10x/tag_mas_sirv_umi_positions.py`, and the research `resources/PB10xMasSeqArraySingleFlowcell/cDNA_array_8x.gb` plus its alignment-reference FASTA.
- mas-seq-paper-data `0c166c63ce4d8ee37c835c89c9693bb65becb51c`: original `README.md` and `scripts/downsampling_analysis/mb_02_long_adata_seurat_sct.R`.

Sources remain under `$CHEM_DATA/sources/pacbio-kinnex-full-length-rna/`; these historical
analysis/reference files do not select commercial Kinnex reagent sequences. No datasets
or third-party source files are committed.
