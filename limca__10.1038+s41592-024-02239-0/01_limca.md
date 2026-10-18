# LiMCA

**Defining source:** Wu, Zhang, … Xie, “Simultaneous single-cell three-dimensional genome and
gene expression profiling uncovers dynamic enhancer connectivity underlying olfactory receptor
choice”, *Nature Methods* (2024), doi:10.1038/s41592-024-02239-0; open full text PMC11166570
(PMC Cloud JATS/PDF, CC BY). The Research Square preprint (rs.3.rs-3210240/v1) is superseded and
was not used.

## Sources read (bounded pass 2026-10-17)

- PMC11166570 JATS text, Methods: “LiMCA protocol”, “Single-cell nucleus–cytoplasm
  separation”, “Cytoplasm Smart-seq2 procedure”, “Single-nucleus chromosome conformation
  capture”, “Single-cell WGA”, “Library construction”, “Sequencing” (xml.txt lines 179–203).
- Supplementary files: MOESM1 (Reporting Summary, no chemistry), MOESM2 (Peer Review File;
  only NlaIII vs MboI resolution remark, line 720), MOESM3 (cell/OR tables; no oligos).
- Oligo scraper on the bundle: only `dT30VN` shorthand and the Nextera mosaic-end trimming
  strings used in cutadapt (analysis, not library oligos).

## Molecular boundary

🟢 One picked cell is lysed in 7 µl soft lysis buffer (0.45% IGEPAL-CA630, SUPERaseIn), 30 min
on ice, vortexed and spun at 500g; 5 µl supernatant (cytoplasm) goes to a new tube, the
nucleus stays in the remaining pellet. The design follows Trio-seq physical separation.

🟢 Cytoplasm, Smart-seq2: 1 µM oligo-dT30VN + dNTP, 72 °C 5 min; RT with SuperScript II,
betaine, 6 mM MgCl2, 1 mM GTP and 1 µM template-switch oligo, 42 °C 90 min then 10 cycles of
50 °C/42 °C; KAPA HiFi HotStart with ISPCR primer, 21 cycles; 0.7× AMPure XP. 1–5 ng cDNA is
tagmented with Nextera-type Tn5 (Vazyme TD502) for library preparation; cells are pooled and
purified with 0.6× and then 0.2× SPRI.

🟢 Nucleus, single-nucleus 3C: 2.5% PFA 10 min RT, glycine quench (with magnetic beads as a
pellet carrier), Hi-C lysis, 0.3% SDS 62 °C then Triton X-100, NlaIII digestion (4 U/µl,
rCutSmart, 37 °C 2 h; the article prints “NIAIII”), T4 ligase 10 U/µl RT 2.5 h, then a
Qiagen-protease lysis buffer at 50 °C 1 h, 65 °C 1 h, 70 °C 15 min; stored at −80 °C. No biotin
fill-in or junction pull-down is described.

🟢 The lysate is amplified by META (Tn5 transposase Vazyme S111-01) “as previously described”
(refs 18, 26) with **20 META tags**.

🔴 The article prints no oligo bases: not the full oligo-dT30VN, TSO or ISPCR sequences, not
the twenty META tags, and not the library index primers. META amplification conditions and how
the META amplicons become a sequencing library are deferred to the cited earlier papers.

## Final libraries and reads

🟢 Both the gDNA (3C) and cDNA libraries were sequenced paired-end 2 × 150 bp on HiSeq X10 or
NovaSeq (3–6 Gb gDNA and 0.2–0.6 Gb cDNA per cell).

🟡 Because cells are pooled after library preparation, each cell's library must carry a
distinguishing index; the page draws it as an inferred `per-cell library index` role. gDNA
and RNA data are joined by experimental labels, not by a shared molecular barcode.

## Open questions

- Author protocol Wu, Zhang & Xie, Research Square doi:10.21203/rs.3.pex-2423/v1 (article ref.
  44; step-by-step protocol and oligos): behind Cloudflare, not retrieved (manifest `(manual)`).
- META tag sequences, amplification and gDNA library construction: Tan et al. 2018/2021
  (refs 18, 26), not inspected in this pass.
- Exact Smart-seq2 oligo sequences/modifications used here (ref. 24 Picelli 2014 not
  inspected); Vazyme TD502 kit adapters and index primers.
- Which fraction is kept in the 0.6× then 0.2× SPRI steps is not stated.

This is a capped bounded pass, not a confirmation that any of this information is absent.


## Native graph migration — 2026-10-18 batch10

Sources read (three originals, capped pass): Published PMC11166570 Methods181–213; directly cited author protocol10.21203/rs.3.pex-2423/v1 originalPDF120–246; original attached LiMCAOligos workbook3b3ac9309196696067b340eb.xlsx, META20/META40/MErows1–64, UDPexamples65onward and Smartseq2 oligos. Attachment manifest identity checked; workbookMD5 31b5de54798ff75a5f17f0c022b4bd86.

Native source-owned compartment graphs replace endpoint/branch recipes. Physical5µl cytoplasmic supernatant and2µl pellet are explicit without assigning purity. RNA firstSS and switchedSS are actual ISPCR substrates; native TD502 input and initial extension feed source9-cycle UDP PCR. Nuclear NlaIII3′CATG/contact DNA feeds actual META20/native tagging, protease,13-cycle PCR,ExoI,META40 two cycles,ExoI,DM two cycles and source recoveries. META technical family is not promoted to UMI/cell barcode. Printed META20/ME pair has13bp overlap and6nt companion overhang; no missing bases silently added. Primary72°C preheat/0.6×+0.2× selection versus v1 70°C/0.55×+0.25× are displayed as version differences.

Remaining scope: Complete native META20 loaded-transposase insertion/first-cycle geometry and DM5xx/DM7xx fullstrands/member/width remain unassigned. RNA UDP example physical10nt does not assign nuclear DMwidth or acquired cycles. Primary/article versus directly cited protocol ordinary version differences remain explicit; no complete revision equivalence claimed.

Verification: 7 original-workbook transcription checks;34graphsteps;43unresolved/0contradictory; page/pairing/schema pass. No missing-information acceptance follows from this pass.
