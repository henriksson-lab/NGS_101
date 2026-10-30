# SiMSen-Seq

DefiningDOI[10.1093/nar/gkw224](https://doi.org/10.1093/nar/gkw224). This is a capped published-source pass and private incomplete model, not a finished schematic or accepted absence audit.

## Selected sources

PublishedPMC4914102 Methods65–75,Results79–93/Fig1; publishedTableS1File008workbook universalprimers/selectedTP2row; File006SIcaptionsS1–S3. Mandatory bundle oligo scraper run; selected source strings/modifications checked at original locations.

## Known preparation boundary

Hairpin-protectedN12 molecularbarcoding inthreePCRcycles; TE/protease65°C then95°C,physical10µLtransferintoQ5libraryPCR18–30cycles;36µLcleaned0.83–1xAMPure,FragmentAnalyzerQC,normalize/pool andMiSeq/HiSeqSEorPE150.

## Source limits

Maintext reports hairpinTm74.0±0.3°C butFig1caption70.01±0.24; retainboth anddonotfit/predictfoldingtemperature. Full molecular products, unreviewed linked originals and source-specific read geometry remain explicit; no perfect yield/affinity/purity or absence claim. Detailed model and verification record follow this source-first note.

## Source sequence construction and physical graph

TableS1 specifies a universal forward component beginning with canonical Read1 sequence plus5′GG, thenN12, thenAT and the14nt closing stem. The left14nt is `GGACACTCTTTCCC`; the closing stem is computed by reverse complement as `GGGAAAGAGTGTCC`. The table's locus primer is appended as a3′overhang. The local constructor validates stem complementarity with neutral `hairpin.Hairpin`, rejects a wrong barcode width, and preserves the source target arm outside the stem. It does not predict folding occupancy or mispriming at a temperature. The generic partial input represents an unassembled target family, not an invented full manufactured primer. The TP2 target row is used only in a sequence-construction/source-transcription check.

The graph begins with source gDNA or cfDNA; extraction and any input-specific shearing remain outside this chosen barcoding entry. Three AccuPrime cycles generate new barcoded DNA while original input/unused primers remain. Adding20µLTE/protease digests Taq at65°C15min;95°C15min follows. Neither step physically removes DNA or proteolysis remnants.10µL is allocated to the40µLQ5 reaction, with an explicit unresolved remainder. PCR2 uses its distinct slow80→72→76°C ramp, adds library arms/sample index and generates new copies; molecular consensus remains computational.36µL goes to0.83–1.0×AMPure depending amplicon length, with other PCR material kept separate. A FragmentAnalyzer analytical aliquot is separate from normalized production libraries and pooling.

Source acquisition includes low-depth MiSeqNanoV2SE150 and full MiSeq/HiSeq2500SEorPE150. Exact cohort allocation, full target/library products, primer landing, historical kit details and directly linked barcoding predecessors remain unreviewed. N12 source molecule identifiers and6nt sample indexes have separate structured features; no perfect one-barcode-per-original efficiency is predicted.

The originalFig1p4 was rendered: the example curve shows74.1°C and plotD clustersnear74°C, while its caption prints70.01±0.24°C; Results text prints74.0±0.3°C. Both locations are retained on the private preview. No source value is silently repaired or used to fit a folding calculation. Safe-SeqS appears as a comparison/control universal primer in TableS1; similarity alone does not authorize a full inherited recipe or documented ancestry edge.

TableS1 also explicitly notes noGGstabilizer forSMA_1/SMA_4 and threeGs forSMA_2. The displayed standardGGcomponent and `StandardProtectedPrimer` constructor do not stand in for those variant orders. Their full bodies and pool allocation remain unassembled; no universalGGrequirement is imposed on every study target.
