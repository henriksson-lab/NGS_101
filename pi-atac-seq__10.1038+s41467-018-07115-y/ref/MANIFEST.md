# Pi-ATAC source attachments

Fetch the defining article and its supplements with
`python3 tools/get_sources.py "Pi-ATAC-seq"`; originals and extracted text live only
under `$CHEM_DATA/sources/pi-atac-seq__10.1038+s41467-018-07115-y/`.

The Nature Communications MOESM12 attachment is an Adobe LiveCycle XFA PDF, served
from the [publisher supplement](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-018-07115-y/MediaObjects/41467_2018_7115_MOESM12_ESM.pdf)
and also fetched by the PMC Cloud route in `MANIFEST.tsv`. Ordinary PDF extraction
shows a “Please wait” placeholder. `tools/doctext.py` extracts its filled XFA dataset
into the ignored `.pdf.txt` twin. The useful methods field is
`Maincontent/General/FlowCytometry/MethodologicalDetails/sampleprep` (line 63 of
the current extraction). It is a reporting summary, not a primer table; it names
fixation, lysis, tagmentation, EDTA quenching and sorting but gives no Ad1/Ad2 base
strings or 96 × 90 plate map.


Source-audit follow-up: registered routes in `catalogue/extra_sources.tsv` recover
public protocols.io API integer 13224 (API metadata is current, not a frozen 2018 snapshot),
its actual overview `https://content.protocols.io/z43fw26.png` and separate current
header image, plus GEO GSE112091 family SOFT metadata via the official FTP host.
The NCBI series HTML route returned a challenge; family SOFT is a validated original
metadata alternative. No raw reads or result tables were downloaded for this check.


The current Index Adapters Pooling Guide is recovered through its public MadCap
`Data/HelpSystem.xml`: the default topic is `Content/IndexAdapterPooling.htm`.
Exact registered topic routes retrieve the entry plus Overview, ColorBalance,
SequencingChemistry, NexteraPoolingGuidelines, NexteraIndexKit and RevisionHistory.
The inspected body identifies document 1000000041074 v17, July 2026; it is not
silently substituted for a historical Pi-ATAC well-assignment source. The blank
HTML5 loader does not establish an access barrier. Named external adapter/readout
references remain listed as source-scope work in the review ledger.


Further registered originals cover the current Adapter Sequences guide (1000000002694
v22), Indexed Sequencing guide (15057455 v10, August 2024) and the actual named
NextSeq/MiniSeq bulletin. Topic URLs are identified from each public HelpSystem/TOC;
project-qualified local names keep generic Overview/RevisionHistory pages distinct.
The old bulletin URL now serves current Illumina Knowledge Article 1241, not a frozen
2016 document. Current/historical orientation statements and obsolete Epicentre
Nextera arms are not silently assigned to Pi-ATAC's unknown custom well pairs.
