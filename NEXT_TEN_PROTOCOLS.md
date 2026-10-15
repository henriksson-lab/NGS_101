# Second ten-protocol pass

This pass applies the local-ownership and source/chemistry audits from
`PILOT_PROTOCOL_OWNERSHIP.md` to ten further documented entries in
`catalogue/ours.tsv`. These methods add bead capture, IVT, bacterial split-pool
ligation, contact indexing, historical SOLiD adapters, bisulfite conversion,
and two library formats with a published oligo conflict. The first 15 audited
protocols are excluded.

| Protocol | Distinct chemistry | Initial graph/source boundary |
|---|---|---|
| MARS-seq | Plate barcode, pooled IVT, ligation | Parameterized graph in its page module; 2014 author oligo tables recovered, local source bundle absent |
| SCRB-seq | UMI RT, template switching, selective 3′ amplification | Parameterized graph in its page module; author Table S2 reviewed, local source bundle absent |
| Seq-Well | Barcoded bead capture and selective Nextera PCR | Local graph; source bundle absent; exact structural inputs incomplete |
| sci-Hi-C | DpnII contact bridge plus two combinatorial barcode rounds | Source-specific molecular states local, no typed graph; original supplement includes bridge/Y-adapter table |
| PETRI-seq | Bacterial RT and two split-pool ligations | Local graph, but split-pool oligos and linkers lack exact graph states |
| Tang 2009 | First single-cell RNA-seq with SOLiD terminal adapters | Local graph; original oligo table reviewed, adapter input state incomplete |
| VASA-seq | Total-RNA IVT with plate/drop branches | Local graph; source prints an incompatible drop RA3/RTP and i5-primer combination |
| inDrop v2 | Droplet hydrogel bead and flipped read orientation | Page delegates to v1 family module; no local graph |
| Plate_scATAC-seq | Bulk tagmentation then indexed single-nucleus PCR | Local graph; transposome has no exact input state; local source bundle absent |
| scBS-seq | Bisulfite DNA conversion and historical indexed PCR | Local graph; inherited primer evidence and final chemical states require audit |

For each method, verify defining and supplementary sources and named vendor originals;
put source-owned claims and the ordered graph in the protocol directory; declare exact
published reagent states and mark unknowns; run the source-oligo guard and scoped
chemistry checker; inspect source conflicts and show them on the public page; then run
its self-tests, pairing lint, documentation and full site build. A published-source
discrepancy is a page notice, while a prior error in this repository is recorded as
a model correction.

Initial discovery found six graphs through the generic source guard, ten exact
inputs checked with no mismatches, eighteen graph inputs without exact states, and
one source-bundle skip. The generic chemistry checker stops at MARS-seq because its
parameterized graph is not discoverable under the local no-argument interface. These
are baseline coverage observations, not evidence of chemical completeness.

## Source intake and public discrepancies

The defining-source fetch has now populated local bundles for all ten. MARS-seq's
2014 supplemental figure file, Tang 2009's main text, inDrop v2's main text and
the later scBS-seq companion main text remain publisher or PMC-gated; their fetched
supplements and author materials are retained. The VASA-seq cited-source fetch
stopped after filling the workspace, but its defining paper and supplements and
several predecessor sources are present. Large Seq-Well video attachments fetched
in this pass were removed to restore disk space; its oligo table and written
Master Protocol remain. With the available bundles, the generic source guard
checks 13 exact inputs, reports no mismatches or skipped bundles, and still lists
18 inputs without exact whole-strand states across six discoverable graphs.

VASA-seq's public page now states the Supplementary Table 12 RA3/RTP versus
`VD_ILMN8_i5` incompatibility, the alternate Table 14 attribution, and why its
drop amplicon stays unresolved. Plate_scATAC-seq's public page now distinguishes
the main-text qPCR optimization's 72 °C for 5 min from the Supplementary Methods
library PCR's 72 °C for 10 min and follows the latter for the drawn workflow.

MARS-seq's ordered graph now lives in its own `tools/mars_graph.py`; the shared
page renderer calls it for both published versions, while the no-argument
`protocol_graph()` selects the 2014 named method for independent audits. The
chemistry checker reports one established, three conditional and eight unresolved
claims with no contradiction. The source guard checks its exact RT2 primer. The
primary 2014 Supplementary Table S7 prints both exact library PCR primers, but
PDF extraction interleaves the role labels between each primer's two sequence
lines; two reasoned source-guard exceptions record this extraction limit rather
than treating an absent whole-string match as a publication gap.

After that migration, the ten-directory source guard discovers seven graphs,
checks fourteen exact inputs with no failures or source-bundle skips, and reports
twenty inputs without exact states. This remains a partial chemistry audit until
SCRB-seq, sci-Hi-C and inDrop v2 expose their own no-argument graphs.

The affected MARS-seq, VASA-seq and Plate_scATAC-seq pages build and pairing lint
finds no bad duplexes. The full site build generated all 405 diagram pages,
accepted all 400 documented NGS graph contracts and resolved every relative
link. These are integration checks for the current partial pass, not completion
of the ten-protocol objective.

## Local graph and source audit across all ten

SCRB-seq's ordered graph and RT scenes now live in its own `tools/scrb_graph.py`;
the shared page module renders the original and mcSCRB variants. The source guard
found that its transposome had `ME` and `S7` in the wrong order. Illumina's Nextera
sequence table supports `S7` then `ME`; the graph was corrected. This was a
repository model error, not an error in the publication. Its three exact inputs
now match the author table or the declared vendor source.

sci-Hi-C now declares DpnII digestion, both 96-well barcode rounds, the intervening
pool and proximity ligation, then second-well reverse crosslinking, AMPure cleanup,
MseI/HaeIII digestion, end repair, dA-tailing, Y-adaptor ligation, second pooling,
streptavidin capture and PCR in `tools/scihic.py`. The detailed Supplementary
Methods orders the four-cutter digest before Y ligation and capture, whereas the
main-text overview lists it after capture; the page now states both locations and
follows the detailed method. The A1 bridge and Y-adapter duplexes are
source-exact; the other wells remain explicit unresolved families. The page uses
this local graph. Source parsing was extended to count `/iBiodT/` as a modified T
base, which lets the A1 bridge pass the exact source guard.

inDrop v2 now assembles its own `tools/indrop_v2_graph.py`. It still imports the
published family's molecular constructors and drawing renderer, but its source
items, reaction order and sequencing plan are local. The author supplementary
oligo table supports its bead seed, PE2-N6, PE2 PCR and representative PE1 index
primer; all four pass exact source checks. The varying barcode plates remain
declared as a set rather than an invented single oligo.
The v2 supplementary bead-QC method disagrees with itself about the intermediate
band (102 nt in the gel instruction, 101 nt in Anticipated Results). Both the v2
page and the v1 page that inherits the bead drawing now state the conflict and
leave the barcode-dependent intermediate length unchosen.

PETRI-seq's graph now carries two named members of each source-printed barcode
set and both exact splint linker sequences. Its linker scenes pair the published
IUPAC S/W ends to source-derived barcode termini; the eight pairing claims are
conditional because S and W admit several base assignments. The other 95 family
members are not individually assembled. scBS-seq's graph now includes the
published Exonuclease I removal of unused oligo1 and AMPure cleanup before
streptavidin capture; its two N9 random primers and PE1.0 are checked against
Smallwood's Online Methods. Its notes now identify the cited Quail supplement,
not the scg drawing, as the source of the historical iPCRTag design. Tang's graph
now carries all four SOLiD adapter strands and classifies dATP as a tailing
condition. The source guard was corrected to compare strands printed 3′→5′ in
their molecular orientation; all ten Tang oligos now pass.

The MARS-seq 2014 RT1 and pool adapter are now source-exact graph inputs, as is
SCRB's E3V6NEXT oligo-dT. The generic scanner now recognizes adjacent
`(N)10(T)30` shorthand in the original SCRB source. The cited Quail 2012
`iPCRtagT1` is a representative exact input for scBS-seq's variable indexed
primer family; the public schematic still shows the family placeholder. Together
these changes bring the ten-graph source guard to 43 checked whole strands,
zero failures, zero missing source bundles and 12 inputs deliberately left
without an exact whole-strand state. These remaining inputs are unspecified
PCR primer pairs, variable barcode families, iso-base modified strands or
vendor kit components whose exact oligo was not selected by the paper.
The cited Quail supplement was already present in the separate
`_data/sources/doi_10.1038+nmeth.1814/` bundle and was read there. A direct
Springer supplement URL is now declared for reproducible fetching into the
scBS bundle; this run's fetch returned a DNS error, recorded as a manual row
in that local manifest. The matching Quail T1 sequence is also printed in
the fetched Clark 2017 supplement held with scBS, so the local exact guard
does not depend on that retry.

All ten selected graphs are independently discoverable. The scoped chemistry
checker currently reports 16 established, 17 conditional, 116 unresolved and
zero contradictory local claims across them. The large unresolved count is a
visible limit of the declarations, not a successful chemical proof. The
source-backed public discrepancies are visible on the affected VASA,
Plate_scATAC, scBS, sci-Hi-C, inDrop v2 and inherited inDrop v1 pages.

Final local gates: the source guard checked 43 exact whole strands in all ten
graphs with no mismatch or missing bundle; the 12 other oligo/adapter inputs
remain explicit. The shared suite passed 1,150 checks and the source-fetcher
suite passed 80. Pairing lint found no bad duplexes in the affected pages.
`build_site.py` built all 405 diagram pages, accepted all 400 documented NGS
graph contracts, and link-checked the assembled site successfully.
