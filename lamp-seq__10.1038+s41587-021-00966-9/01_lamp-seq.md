# LAMP-Seq — source notes

Defining source: Ludwig *et al.*, *Nature Biotechnology* (2021),
[doi:10.1038/s41587-021-00966-9](https://doi.org/10.1038/s41587-021-00966-9).

## Sources read

- 🟢 Main article and Methods, PMC8678193.
- 🟢 Supplementary Data 1, including the primer list, 192 validated Tier-1 FIP
  barcodes, the larger Tier-2 set and pool-index primer table.
- 🟢 Supplementary Notes and application note.

## Molecular interpretation

🟢 Six target primers drive RT-LAMP: F3, B3, FIP, BIP, LF and LB. The 10-nt
sample identifier is physically embedded between the F1c and F2 portions of FIP.
The published open-source implementation uses LNA substitutions in F3 and B3 and
a strand-displacing Bst polymerase at 65 °C for 80 minutes.

🟢 The paper gives two PCR branches. The two-stage branch shown on the page pools
and dilutes LAMP products 1:100,000 before an 18-cycle target PCR, then uses
another 18-cycle pool-index PCR. Adding 5 µl of diluted pool to the 50 µl first
PCR makes the effective dilution 1:1,000,000. The open-source one-step branch
instead dilutes 1:40,000 before a 25-cycle indexed PCR (effective 1:400,000).
The dilution limits barcode swapping caused by relabelling across pooled wells.

🟢 The PCR primers install two additional 10-nt pool identifiers together with P5,
P7 and standard TruSeq sequencing arms. Supplementary Data 1 publishes every
primer; the schematic uses the first validated FIP barcode (`AGCAAACTAC`) and
Pool_fwd_1/Pool_rev_1 as representatives.

🟢 The pool-barcoding primers carry only partial sequencing-primer tails
(Supplementary Data 1 `Pool barcoding primers` row 1; Methods two-step PCR, PMC8678193 text line 277–279):
Pool_fwd_1 is P5 + `TCGTGGAGCG` + `{{= lamp_seq.LAMP_POOL_FWD_TAIL }}` ({{= len(lamp_seq.LAMP_POOL_FWD) }} nt) and Pool_rev_1 is
P7 + `CGCTCAGTTC` + `{{= lamp_seq.LAMP_POOL_REV_TAIL }}` ({{= len(lamp_seq.LAMP_POOL_REV) }} nt). The first PCR supplies the rest of
both arms. PCR-N-rev-11 prints its tail as `TGACTGG…` in both the Methods and
Supplementary Data 1; the 5' G of the pool reverse tail completes canonical `GTGACT`,
so this is a design, not a source error. *Model correction (2026-10-17):* earlier
versions appended full 33/34-nt TruSeq arms to the pool primers (72/68 nt).
The full-arm forms belong to the open-source one-step Ultramer primers (Methods, line 353–355).

🟢 The two-step branch drawn is the commercial-reagent protocol (Methods lines 205–289):
QuickExtract lysis with pX330 carrier DNA, ≥10 min RT then 95 °C 5 min; 30 min
activated-carbon/ion-exchange bead treatment and HCl pH adjustment; NEB LAMP master mix
E1700L plus Bst 3.0; C-F3-LNA, C-B3-LNA, C-BIP, C-LF, C-LB and barcoded C-FIP; 65 °C
80 min; pooling on ice; 1:100,000 dilution; two 18-cycle NEBNext 2× PCRs (98 °C 20 s,
65 °C 20 s, 72 °C 30 s); pooling, 2% E-gel band excision, two silica spin columns,
NanoDrop. The Methods call the target primers PCR-C-fwd/PCR-C-rev; Supplementary
Data 1 prints the identical sequences as PCR-N-fwd-10/PCR-N-rev-11. The commercial
master mix also contains an ACTB control primer set and barcoded ACTB-FIP. Their native primer/control-product roles and distinct structured barcode identity are retained, while ordered strands and barcode width/member remain unassigned.

🟡 LAMP does not yield a unique linear intermediate. It yields stem-looped,
branched concatemers containing repeated and inverted target segments. The final
PCR selects compatible units from that population. Consequently, the page preserves
the exact primer landing sequences and embedded barcode but does not assert one
universal order or length for the intervening LAMP-derived sequence.

## Sequencing

🟢 The MiSeq Nano and iSeq runs use Read 1 100 cycles, Index 1 8 cycles and Index 2
8 cycles (Methods lines 365–385, open-source section); Read 2 is not acquired. Only
the first 8 nt of each 10-nt pool index are therefore read. 🟡 The commercial branch
names MiSeq or iSeq without uniquely linking its preparation to that cycle recipe; the graph keeps actual acquisition unspecified. Read 1 is retained after extracting the embedded
FIP barcode and confirming non-primer viral sequence. The two index reads identify
the post-pooling PCR group.

## Native graph migration (2026-10-18, batch 12)

Read two published originals: author-manuscript `s41587-021-00966-9_PMC8678193.xml.txt`, commercial Methods lines 205–297 (lysis, adsorbent transfer, independent ACTB control, native LAMP primers, physical pooling/dilution, two PCR aliquots, library pooling/gel/two columns and diagnostic read analysis), plus named sequencing sections 365–387; and Supplementary Data1 `MOESM3_ESM.xlsx.txt`, Primer List lines 1–18, Tier1 representative barcode rows 64–66 and pool-primer pair1 row 4102. Mandatory full-bundle scraper run; relevant sequence and modification hits were checked at those original locations.

Removed the shared LAMP product/topology recipes and canonical amplicon shell. The local final component reference is assembled from the printed primer tails and representative pool1 indexes, preserving actual ten-nucleotide i5/i7 identities. Both incoming PCR states are now actual native mixtures/aliquots, not the proposed selected linear interval or completed final library. First-cycle selection among branched/inverted LAMP products remains unresolved.

Original viral/human RNA, pX330 DNA from lysate and commercial master mix, independent ACTB control products, original LNA-modified outer primers, ordinary copied DNA, heterogeneous LAMP DNA, PCR daughters and possible residual original templates/oligos remain distinct. LNA sugars on original primers are not copied into daughter DNA. This graph does not assign diagnostic sensitivity, decision thresholds or clinical performance. Correct genome sequence excluding primer-covered portions and barcode co-occurrence are computational diagnostic analysis.

The source lysate/pre-PCR/second-PCR aliquots have explicit physical remainders with unresolved fate. Adsorbent treatment, fresh-tube transfer and pH adjustment are separated. Independently barcoded wells and independently indexed PCR libraries enter actual pool operations; the companion members/ratios are unassigned and do not inherit representative pair1 identities. Dilution does not physically remove molecular classes. Gel-band selection and two silica columns each expose physical retained/nonretained fractions without proving class purity.

The commercial branch names MiSeq/iSeq without uniquely assigning its acquired recipe. The separately printed named sequencing recipe is 100 Read1 plus 8 i7 plus 8 i5, covering only eight of each physical ten-base pool index; no Read2 is acquired in that recipe. The graph keeps actual acquisition unspecified while retaining source-supported primer-site and physical index geometry. Open-source single-PCR, ACTB detailed primer architecture, N100/PhiX admixture and other diagnostic-validation branches remain outside the selected graph.

Deferred: native LAMP stem-loop/isomer assemblies and first-cycle PCR selection; detailed ACTB/control product architecture and actual pool membership; modified-outer-primer individual lineage fate; clinical readout analysis/thresholds; exact selected-library linkage to run recipe and N100/PhiX admixtures; cited predecessors and complete linked-source review. Existing missing-information status is unchanged.

## Independent implementation review and bookkeeping corrections (2026-10-18)

No original documents were reread. The implementation review checked native incoming PCR mixtures, independent residual and modified-primer lineages, actual aliquots and pools, structured identifiers, and source-qualified acquisition against the existing migration record and research note.

Corrected the stale sequencing paragraph that applied the separately named 100/8/8 recipe to the commercial branch. Its actual run linkage remains unspecified. The graph labels now retain the commercial conditions already recorded above rather than generic “source condition” wording. The barcoded ACTB-FIP primer and its independent control DNA now share a distinct structured control-barcode feature with unassigned width/member; the viral representative barcode is not assigned to that control. These are documentation and identifier corrections, not new source findings or an extension of the source investigation. Existing missing-information dispositions remain unchanged.
