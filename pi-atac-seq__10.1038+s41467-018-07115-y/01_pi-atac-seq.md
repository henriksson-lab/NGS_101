# Pi-ATAC-seq — research-note pointer

Pi-ATAC-seq is documented from Chen *et al.*, 2018,
doi:[10.1038/s41467-018-07115-y](https://doi.org/10.1038/s41467-018-07115-y).

The detailed source record was originally written jointly with Plate_scATAC-seq and is
retained as provenance in the
[combined research note](../plate-scatac-seq__10.1038+s41467-018-07771-0/01_plate-scatac-seq.md#3b-step-by-step-pi-atac-differences).
The checked Pi-ATAC article and supplements do not give a well-to-primer map;
the schematic marks its reconstructed library regions as inferred. The recovered
MOESM12 LiveCycle PDF is a reporting summary, not a primer table. Its filled form
field `Maincontent/General/FlowCytometry/MethodologicalDetails/sampleprep` says
the cells were fixed in 1% formaldehyde for 10 minutes, quenched with 0.125 M
glycine for 5 minutes, permeabilized, stained, tagmented in bulk for 30 minutes,
quenched with 40 mM EDTA and index sorted. This was read from the original XFA
data with `tools/doctext.py`; it does not settle the indexed oligo assignments.

The authors' public [Pi-ATAC protocol](https://www.protocols.io/view/pi-atac-q6gdzbw.pdf)
(doi:10.17504/protocols.io.q6gdzbw, 23 October 2018), pp. 6 and 10, was also
checked. It specifies Nextera kit FC-121-1030 and adds `Ad1` and `Ad2` as distinct
25 µM primer stocks during per-well PCR, but prints neither primer sequence nor an
index-plate table. This confirms the primer roles without resolving their bases or
modifications. It cites Buenrostro *et al.* 2015 for Nextera adapters; those
cannot establish which indexed `Ad1`/`Ad2` oligos occupied each Pi-ATAC well.

The linked source was followed in this audit: Buenrostro *et al.*, *Nature* 523,
486–490 (2015), doi:[10.1038/nature14590](https://doi.org/10.1038/nature14590),
Supplementary Table 1 (`nature14590_41586_2015_BFnature14590_MOESM36_ESM.xlsx.txt`).
It prints unindexed `v2_Ad1_NoMX` and `v2_Ad2_NoMX` at lines 3–4, 92 indexed
`v2_Ad1` entries at lines 6–97 and 96 indexed `v2_Ad2` entries at lines 99–194.
The oligo scraper found these original table rows; the Pi-ATAC protocols.io PDF itself
has no sequence hits. 🟡 These predecessor primers support the inferred Nextera arms,
but no checked Pi-ATAC source maps its reported 96 × 90 combinations to those entries.


## Published conditions and source boundaries

The article Methods, **Index sorting and library preparation for Pi-ATAC**, states
20 µl reverse-crosslink buffer, 65°C overnight and 25 mM barcoding-primer stock.
The original author protocol dated **23 October 2018**, steps 26, 32 and 34, states
25 µl buffer, 60°C and 25 µM Ad1/Ad2 stocks. These published versions differ; this
schematic does not choose among their temperatures, volumes or stock concentrations.

The protocol's step 25 literally prints **200 µg/µl** proteinase K and a **1:100**
dilution of **20 mg/ml** stock. The stated dilution computes to **200 µg/ml**, an
apparent internal unit discrepancy. The article's reverse-crosslink buffer description
instead prints **5 ng/ml**. No disputed dose is assigned to the molecular graph.
These original locations were checked independently by two readers; they do not prove
which conditions were actually used.

The original public protocols.io API for historical integer ID **13224** exposes
an overview figure and current header-image metadata. The GEO accession **GSE112091**
is a SuperSeries with eight SubSeries. Its family SOFT metadata contains sample
protocols and sample identifiers; bulk peak-file names include some literal indexes,
which do not establish the single-cell well-to-primer map. Direct source/link inspection
continues; no publication-wide absence claim is accepted.

## Bounded curation round 3 (2026-10-16)

Inspected: Published PMC6214962 JATS202–226.

Added overnight65°C reverse crosslinking,80°C/10min proteinase-K inactivation and20-cycle PCR labels. Source confirms75-base paired genomic reads. Candidate Tn5 products no longer presumed physically discarded.

Deferred: Published supplementary oligos, precise index-read mask, gel recovery stages;25mM primer stock wording not assigned a corrected value. This capped pass does not establish a complete source review or publication-wide absence.
