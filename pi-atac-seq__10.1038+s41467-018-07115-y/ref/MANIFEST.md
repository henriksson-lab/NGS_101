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
