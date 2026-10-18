# scSPLAT — research notes

Defining source (published version, used as authority): Raine *et al.*, "scSPLAT, a scalable
plate-based protocol for single cell WGBS library preparation", *Scientific Reports* (2022),
[10.1038/s41598-022-09798-2](https://doi.org/10.1038/s41598-022-09798-2) (PMC8986790), Methods
"scWGBS library preparation and sequencing" and Supplementary Table S1 (MOESM1). The directory is
named for the superseded bioRxiv preprint
[10.1101/2021.10.14.464375](https://doi.org/10.1101/2021.10.14.464375); it is not used for active claims.

## Sources read

- Published JATS (PMC Cloud `PMC8986790.1.xml`, text twin lines 59, 63, 131–135): Results first
  section, Fig. 1 legend, Methods library preparation, adapter preparation and demultiplexing.
- Published Supplementary Table S1 (MOESM1 text twin lines 1–150): SSO_H_1–32, ss1a, ss1b.
- Not yet read: the original SPLAT method (Raine *et al.*, *NAR* 2017, doi:10.1093/nar/gkw1110,
  reference 14), which the Methods cite for the adapter sequences; the snmC-seq2 low-volume
  conversion protocol (reference 5); Zymo EZ Methylation Direct and NEBNext primer documents.

## Oligos (verbatim, published Table S1 / Methods)

- 🟢 `SSO_H_1`: `5AmMC6/TTCCCTACACGACGCTCTTCCGATCTTAGCTTHHHHHHHH`. SSO_H_2–32 differ only in the
  six-base barcode (positions 27–32), e.g. SSO_H_2 `…ACTTGA…`, SSO_H_3 `…GCCAAT…`.
- 🟢 `ss1a`: `5'AmMC6/GACGTGTGCTCTTCCGATCTNNNNNN/3'AmMo`.
- 🟢 `ss1b`: `5'Phos/AGATCGGAAGAGCACACGTC/3'AmMo`.
- 🟢 PCR: "NEBNext dual index primers" (sequences not printed) with KAPA HiFi 2× mastermix.

## How they interlock (computed by the model)

- 🟡 The 26-nt SSO_H tail is the 3′ end of TruSeq Read 1, and ss1b is the 5′ start of the
  TruSeq Read 2 / Index 1 arm; both are asserted when `tools/chem.py` is imported.
- 🟡 ss1a's 20-nt arm pairs with ss1b over its full length (checked by `Scene`), leaving the
  3′ N6 as the splint that captures the single-stranded 3′ end beside the ss1b 5′ phosphate.
- 🔴 The NEBNext primer bases (P5, i5, first seven Read 1 bases, Read 2 remainder, i7, P7) are
  not assigned as native components in this capped reconstruction; no completed canonical shell or default index width is substituted.

## Step by step (Methods, JATS line 131)

1. FACS into 2 µl Zymo M-digestion buffer with proteinase K; 15 µl EZ Methylation Direct reagent.
2. Reduced-volume 384 spin-plate binding, M-wash, M-desulphonation (15 min RT) and two washes;
   elution in 7 µl buffer containing 500 nM SSO_H_X.
3. 95 °C 2 min, ice; Klenow exo− (2.5 U, 2× Blue buffer, 1 mM dNTP) 4 °C 5 min, 25 °C 5 min,
   37 °C 60 min.
4. Exonuclease I (0.5 µl), 37 °C 30 min, digests excess oligos.
5. Pool up to 32 wells; 1.2× AMPure XP; elute 10 µl Tris pH 8.5.
6. 95 °C 2 min with 50 ng ET-SSB, ice — separates the second strand from the original strand.
7. SPLAT: adapter 5 µM, T4 ligase buffer, PEG4000, 15 U T4 DNA ligase HC, 10 U T4 PNK, 15 µl,
   20 °C 1 h; 1:1 SPRI.
8. 2.5 µl NEBNext dual index primers + 12.5 µl KAPA HiFi, 12–16 cycles; two 0.8× SPRI.
9. PE150 on HiSeqX (liver batch 1) or NovaSeq 6000.

The released original converted strand lacks the SSO_H tail; the model leaves its fate in the
bulk ligation unresolved rather than claiming physical removal.

## Final library and reads

Native outer-primer boundary–SSO tail–cell barcode (6)–H8 random-primer copy–copied insert–ss1b–opposite native boundary. 🟢 Read 1 begins with the six-base cell barcode (demultiplexed with je). 🟢 Analysis
clips 15 bases from the 5′ and 10 from the 3′ end of reads (computational, not a separation).
🟡 Index-read cycles are not stated; native acquisition remains unassigned, separately from physical pool indexing.

## Published discrepancy

🟡 The Methods say the SPLAT adapter was prepared "from oligo ss1a and ss2b (Supplementary
Table S1)", but the same paragraph and Table S1 print only `ss1b`. The model uses the printed
ss1b sequence (it pairs with ss1a); whether ss2b is a typo or a distinct reagent is open until
the 2017 SPLAT paper is checked. The discrepancy is shown on the schematic page.

## Protocol boundary

The page retains the cell-first barcode/pool/splint order that distinguishes scSPLAT from PBAT and snmC-seq2.


## Native graph migration (2026-10-18, batch06)

Two originals were reread: defining published JATS Methods131–135 and published
MOESM1 TableS1 opening rows11–65 and closing rows120–149; middle family rows were not separately reread. The mandatory scraper returned five representative
hits, including the32-member SSO family and blocked ss1a/ss1b. Preprint material
was excluded; no absence follows from extraction counts.

The graph and endpoint declarations are protocol-owned. Actual SPLAT-ligated ssDNA
enters index PCR, retaining the source5′AmMC6 and ss1b3′AmMo chemistry in intermediate
states. The ss1a3′random splint also retains its source3′AmMo block; a splint is not
a ligating hydroxyl. Physical bisulfite U originals remain distinct from Klenow
copies. Exonuclease hydrolysis is not physical discard; original converted-strand
ligation/recovery and off-target amplification remain unresolved. Native NEBNext
outer ends, selected members/index widths, complete PCR geometry and sequencing
primer sites are explicit roles. Source PE150 and six-base inline Read1 barcode
remain positive facts; computational clipping is not a physical selection.

Deferred originals: cited2017SPLAT, linked low-volume conversion protocol and
applicable Zymo/NEBNext documents. Historical missing-information status unchanged.

The two source post-PCR0.8×SPRI recoveries are separate successive graph steps. Native blocked ss-template/off-target classes are unresolved PCR outputs; copying does not erase them or prove their later removal.
