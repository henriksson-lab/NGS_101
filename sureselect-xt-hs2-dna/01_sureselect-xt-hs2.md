# SureSelect XT HS2 DNA target enrichment

## 1. What it is

Agilent SureSelect XT HS2 converts fragmented genomic DNA into dual-indexed Illumina libraries, optionally adds inline molecular barcodes (MBCs) at both insert ends, pools indexed libraries before capture, hybridizes them to biotinylated RNA baits, retrieves bait–library duplexes on streptavidin beads and performs post-capture PCR.

## 2. Authoritative source

🟢 Agilent, *SureSelect XT HS2 DNA Kits: Library Preparation (+/−MBC), Overnight-Hyb Target Enrichment, Pre-capture Pooling Workflow*, G9957-90000 version B0, May 2026. The manual explicitly specifies mechanical or enzymatic fragmentation, combined end repair/dA-tailing, adaptor ligation, pre-capture indexing PCR, blocker-plus-probe hybridization, streptavidin capture, post-capture PCR and standard Illumina paired-end sequencing.

🟢 The manual states that MBC-tagged libraries carry inline molecular barcodes at both insert ends. Step 4 "Process and analyze the reads" (`G9957-90000.pdf.txt` lines 3528–3567) masks or trims the first 5 bases of Read 1 and Read 2 (base mask `N5Y*,I8,I8,N5Y*`) and gives one approach: "trim off 5 bases at the beginning of each read, take the first 3 of the 5 bases as MBC, and discard the remaining 2". The model therefore draws a 3-nt MBC plus 2 non-MBC bases at each end (previously a 5-nt MBC).

Source bundle (hand-downloaded by the maintainer 2026-10-17; identities checked on the printed title pages): `G9957-90000.pdf` (Version B0, May 2026, 78 pp., the applicable manual workflow), `G9985-90011.pdf` (Bravo automation Option B, Version A2, October 2023; spot-checked only), `XTHS2_Index_Sequence_Resource.xlsx` (PN G9983-90600 version A.0; pairs 1, 25 and 384 match manual Tables 61–68). The scraper finds no oligo in any of them (indexes are 8-nt, below its length threshold); `grep` for runs of ≥12 A/C/G/T finds no adaptor or primer bases in either manual.

## 3. Construct logic

All 🟢 from G9957-90000 B0 unless marked.

1. Fragment 10–200 ng gDNA in 50 µL by Covaris shearing (duration by read length; Table 12) or SureSelect enzymatic fragmentation (Chapter 2).
2. End Repair-A Tailing Buffer/Enzyme Mix: 20 °C 15 min, 72 °C 15 min (Table 19).
3. Add ligation master mix (Ligation Buffer + T4 DNA Ligase), then 5 µL XT HS2 Adaptor Oligo Mix (MBC-tagged, white cap) or SureSelect MBC-Free Adaptor Oligo Mix (black cap); 20 °C 30 min (Table 21). AMPure XP 80 µL, elute 35 µL (Table 22).
4. Pre-capture PCR with Herculase II and one SureSelect XT HS2 Index Primer Pair: 98 °C 2 min; 8–14 cycles of 98 °C 30 s / 60 °C 30 s / 72 °C 1 min by input amount and quality (Tables 23–24); 72 °C 5 min. AMPure XP 50 µL (Table 26).
5. Pool 8 or 16 indexed libraries to 4 µg per capture, dry to 4 µL (Table 30); add XT HS2 Blocker Mix, 95 °C 5 min, hold at 65 °C or 67.5 °C (Tables 32–33); add probe, RNase Block and hybridization buffer; hybridize 16–24 h. 🟡 The probe is RNA per the Bravo guide's troubleshooting text ("The RNA Probe", `G9985-90011.pdf.txt` line 20077) and the RNase Block; biotinylation is implied by streptavidin capture.
6. Capture on SureSelect Streptavidin Beads, 30 min room temperature; Wash Buffer 1, then six 70 °C Wash Buffer 2 washes; resuspend in 25 µL water.
7. Post-capture PCR on the beads with Herculase II and the SureSelect Post-Capture Primer Mix: 10–16 cycles by probe size (Table 40); remove the beads on a magnet; AMPure XP 50 µL (Table 42). The mix's two end-specific members are a model inference; their bases are not published.

## 4. Uncertainty and scope

🟡 The B0 guide prints 384 paired eight-base indexes and the five-base MBC read layout, but not complete HS2 adaptor or indexing-primer arms. The schematic uses inferred canonical TruSeq arms only for orientation and primer landing; these bases are not established as the Agilent kit oligos. No specific index pair is chosen. Capture baits are design-specific and are shown as a complementary target interval.

## 5. Read layout

Read 1 and Read 2 each begin with five adaptor-tag bases (3-nt MBC + 2 discarded, per one suggested approach) before target DNA. Index 1 and Index 2 are fixed 8-cycle reads of the paired eight-base sample indexes. Read 1/Read 2 length is the user's choice (Table 46: 2 × 100 or 2 × 150 bp; Table 47 example 151 + 8 + 8 + 151 cycles). No custom primers are used. Primer placement is computed on the inferred canonical Illumina shell; it does not verify the exact proprietary HS2 arm bases.

## 6. Bounded pass 2026-10-17 (round deferred-followup-03)

First located reading of the hand-downloaded originals (the 2026-10-16 passes used web-tool extractions, later withdrawn). Read in `G9957-90000.pdf.txt`: title/edition (1–40), What's New in B0 (144–148), Chapter 3 library preparation (1413–2090), Chapter 4 pooling/hybridization/capture/post-capture PCR (2225–3060), Chapter 5 sequencing Tables 46–47 and Step 4 read processing (3405–3567), Table 61 header. Corrections: 3-nt MBC + 2 non-MBC read-start bases instead of a 5-nt MBC, with the Read 2 site's terminal base marked as an inferred arm base rather than the insert dA; named T4 ligase, Herculase II, Blocker Mix, Post-Capture Primer Mix; added AMPure clean-ups with discard, pre-capture pooling, bead removal after on-bead PCR; reaction temperatures and cycle ranges; 8-cycle index reads. Still open: HS2 adaptor, index-primer and post-capture primer bases; adaptor end geometry and whether the 2 discarded bases include the ligation T; MBC-free branch not drawn; Fast-Hyb (G9985-90000) workflow not read.


## Critical re-review, 2026-10-17

The source-supported two non-MBC read-start bases now have explicit length_bp=2; unknown base identities no longer make their known read span variable.

Inspected scope: G9957 manual .txt 3565–3567; existing source-qualified molecular-barcode/read-start declaration. This bounded correction does not complete the linked-source or two-reader absence audit.
