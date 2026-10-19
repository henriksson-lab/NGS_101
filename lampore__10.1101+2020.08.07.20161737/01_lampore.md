# LamPORE

## Defining source

🟢 James *et al.*, “LamPORE: rapid, accurate and highly scalable molecular screening for SARS-CoV-2 infection, based on nanopore sequencing”, medRxiv (2020), [doi:10.1101/2020.08.07.20161737](https://doi.org/10.1101/2020.08.07.20161737).

🟢 Read 2026-10-17: medRxiv **v3** (posted 30 Oct 2020; latest of three versions, unpublished per the medRxiv API), obtained from its open medRxiv location (`2020.08.07.20161737v3.full.pdf.txt`): Methods 1 "Amplification and library preparation" lines 196–256, Methods 2 data analysis lines 257–378, Results "Assay design" lines 380–398, "Barcode demultiplexing" lines 948–968, "FIP barcode optimisation" lines 1002–1022, references lines 1690–1790. Its Supplementary Table 1 (clinical results) was not obtained; it is not a chemistry source. Workflow claims were earlier checked against the later clinical paper, Peto *et al.*, *J Clin Microbiol* (2021), [doi:10.1128/JCM.03271-20](https://doi.org/10.1128/JCM.03271-20) (PMC8315953, JATS Methods "LamPORE", `.xml.txt` lines 218–220), which used the ONT LamPORE protocol dated 1 July 2020 (version 1, revision 4). That vendor protocol revision is not in the bundle. The JCM paper is a separate clinical study (its own preprint is medRxiv 10.1101/2020.09.18.20195370, cited as reference 20 of the defining preprint, lines 1820–1824, with the JCM title and first author), not a published version of the defining preprint; no published version of the defining preprint is recorded in the source bundle or manifest (checked 2026-10-19).

## Molecular workflow

🟢 Four multiplexed RT-LAMP primer sets target SARS-CoV-2 ORF1a, E and N plus human ACTB. One of eight sample barcodes is embedded in each forward inner primer (FIP), so it is copied repeatedly into the branched LAMP products.

🟢 (preprint lines 198–222) Primer sets come from NEB and cited work; short barcodes are added to the FIPs as in LAMP-Seq (Schmid-Burgk) and Yamagishi; primers are IDT HPLC-purified and actin primers are deliberately at lower concentration. One 10x primer pool per FIP barcode; each plate row receives one barcoded FIP mix. Each reaction holds 25 µl 2x LAMP Master Mix (NEB E1700), 5 µl 10x primer pool and 20 µl RNA (50 µl total, computed), incubated at 65 °C for 35 min then 80 °C for 5 min. The JCM summary (65–80 °C, 40 min, line 220) is consistent. Eight of an initial twelve barcoded FIPs were kept after sensitivity screening (lines 1002–1022).

🟢 (preprint lines 222–256) Reactions are pooled by column (12 pools of 8). Each pool gets a 10 µl rapid-barcoding reaction (6.5 µl water, 1 µl pooled LAMP product, 2.5 µl SQK-RBK004 rapid barcode), 30 °C 2 min then 80 °C 2 min. All reactions are pooled into one tube, purified with 0.8x AMPure beads, washed in fresh 80% ethanol and eluted in 15 µl EB. 11 µl eluate receives 1 µl rapid adapter (RAP), 5 min at room temperature, and is sequenced on one MinION flow cell for 1 h. After RBK the fragments have a modal length of about 500 bp, from LAMP strands up to about 5 kb (lines 950–956).

🟢 (JCM lines 218–220) Eight differently FIP-barcoded reactions are pooled by plate column. A rapid transposase-based step adds one of twelve second barcodes, yielding 8 × 12 combinatorial sample identities. All samples are then pooled into a single library without normalization, the pool receives a magnetic bead cleanup, and it is loaded on a MinION flow cell (GridION instrument) and sequenced for 60 min. In that study FIP7 was not used (lower β-actin read counts), so plates held 80 samples plus 2 no-template and 2 synthetic-RNA positive controls; this is study-specific, not a protocol change.

🟢 The preprint resolves two earlier gaps in the JCM Methods: RAP is added after the AMPure cleanup, and the beads are 0.8x AMPure. The adapter attachment chemistry and bases remain commercial and undisclosed, so the motor-adapter input stays an inferred molecular state.

🔴 The commercial FIP-barcode, rapid-barcode, attachment-tag and motor-adapter sequences are not disclosed by the defining sources. The page therefore represents their molecular roles and does not invent bases or one linear consensus for heterogeneous LAMP concatemers.

🟡 The defining preprint (Assay design, lines 386–394; references 10 and 14) cites the ORF1a primer set "AS1" to Rabe & Cepko (medRxiv 10.1101/2020.04.23.2007687, not read) and E1, N2 and actin to the Zhang *et al.* bioRxiv preprint. Whether Zhang's "As1e" set is identical to the AS1 set used is not checked here.

🟢 The later clinical LamPORE paper (Peto *et al.*, PMC8315953, Methods) says its ORF1a, E and N target sets use previously published primers and cites Zhang *et al.*, [doi:10.2144/btn-2020-0078](https://doi.org/10.2144/btn-2020-0078). The published BioTechniques full text was not obtained (manual row), so the primer bases come from Table 1 of the [bioRxiv version](https://doi.org/10.1101/2020.06.03.132894) prints six unbarcoded primers each for As1e (ORF1a), E1, N2 and ACTB. As examples, its F3 primers are As1e `CGGTGGACAAATTGTCAC`, E1 `TGAGTACGAACTTATGTACTCAT`, N2 `ACCAGGAACTAATCAGACAAG`, and ACTB `AGTACCCCATCGAGCACG` (bioRxiv v1 `pdf.txt` lines 275–330, checked 2026-10-17; preprint provenance, to be replaced by the published table when obtained).

🟡 The cited set gives defensible target-binding **precursor** oligos, but the LamPORE product adds a sample barcode into each FIP and the clinical paper does not identify the complete resulting FIP molecules. The cited unbarcoded sequences should not be copied into a LamPORE ordering list as finished FIPs.

## Read interpretation

🟢 (JCM line 220; preprint lines 335–378, same thresholds chosen to maximise F1 for 1 h runs) The report counts reads assigned to each target: invalid < 50 classified SARS-CoV-2 + β-actin reads; positive ≥ 50 SARS-CoV-2 reads; inconclusive 20–49; negative < 20.

🟢 (preprint lines 257–333) Classification uses the rapid barcode (guppy_barcoder 4.0.11, SQK-RBK004, mid-strand barcodes detected), the repeatedly copied FIP barcode, and target-derived sequence assigning a read to ORF1a, E, N, ACTB or unclassified. Nanopore acquisition is primer-free and variable-length.

🟢 (preprint lines 962–968) Every LamPORE read carries an ONT rapid barcode at one end; because not every repeat unit carries the FIP barcode, reads containing a LAMP barcode are selected, typically retaining about 70% of reads. This is computational read selection.

🟡 The physical fate of fragments without a rapid-barcode attachment end is not stated (all reactions are pooled before the 0.8x AMPure cleanup), so the model marks them as an unresolved terminal rather than discarded; the LAMP-barcode read selection above concerns a different population and does not assign this fate.
