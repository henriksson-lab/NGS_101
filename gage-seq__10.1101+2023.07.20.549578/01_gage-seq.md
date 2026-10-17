# GAGE-seq — source notes

Defining preprint: Zhou et al. (2023), [doi:10.1101/2023.07.20.549578](https://doi.org/10.1101/2023.07.20.549578). Published version: Zhou et al. (2024), [doi:10.1038/s41588-024-01745-3](https://doi.org/10.1038/s41588-024-01745-3).

- 🟢 The full preprint, published paper, supplementary protocol and the published Supplementary Table 1 workbook were read.
- 🟢 RNA is reverse-transcribed with biotinylated poly(dT) or random-hexamer primers; genomic DNA undergoes in-situ contact ligation.
- 🟢 DNA and cDNA receive two rounds of 96-well ligation barcodes using analyte-specific adapter sets. Streptavidin beads then separate cDNA from DNA.
- 🟢 Both branches are preamplified, Nextera-tagmented, indexed and sequenced PE150.
- 🟢 The version-of-record workbook supplies the material that was missing from
  the preprint bundle: every 96-well RNA and Hi-C adapter plate, strand-specific
  modifications, RNA splints and blockers, the post-TdT primers, Hi-C pre-amplification
  primers, and terminal index-primer sets.
- 🟢 The two RNA ligation rounds can be reconstructed exactly. Linker 1 bridges
  the RT-primer handle to Round1_01; the RNA linker 2 bridges Round1_01 to Round2_01.
  Both complementary arms are checked by the shared splinted-ligation constructor.
- 🟢 The terminal shells differ by branch, as Supplementary Figure 1 indicates:
  chromatin-contact molecules use Nextera P5 with TruSeq P7, whereas RNA molecules use
  TruSeq P5 with Nextera P7. Representative A1/index primers are modeled exactly;
  alternate wells change only their printed whitelist sequences.

## Bounded source pass 2026-10-17 (round breadth-first-10)

Read: published PMC11323187 JATS Methods lines 141–171 (lysis, RT, 1st/2nd-round
chromatin fragmentation, combinatorial barcoding, reverse crosslinking/separation,
sequencing) and data-processing line 177; publisher Supplementary Methods
(MOESM1) lines 91–131 (scHi-C and scRNA-seq library construction) and Figure S1
caption lines 535–551; Supplementary Table 1a–1i workbook (MOESM3) A1 rows.
The bioRxiv preprint was not re-read (superseded).

- 🟢 Order (Methods): high-salt lysis and SDS → in-situ RT (Maxima H Minus,
  biotinylated RT primers) → MseI + CviQI digestion (25 °C, 2 h) → T4 DNA ligase
  proximity ligation (16 °C overnight) → DdeI re-digestion (37 °C, 1 h) → round-1
  ligation in 96 wells (25 °C, 3 h) → pool, ~7,500 cells per 96-well plate, 60 °C
  SDS lysis ≥24 h → round-2 ligation with Hi-C-AD2, RNA-AD2, 5′-P TNA/TA
  Nextera-P5-AD and RNA ligation-1 block (25 °C, 24 h) → proteinase K, 60 °C
  reverse crosslinking → isopropanol precipitation → MyOne C1 streptavidin
  separation. The Methods describe no end fill-in: MseI/CviQI (TA) and DdeI (TNA)
  5′ overhangs are what the AD-1 TNA/TA ends are designed to match (Figure S1).
  The model previously had a fill-in step and RT after ligation; corrected.
- 🟢 Hi-C supernatant (Supplementary Methods): precipitation, ligation of
  5′-P-TNA/TA-Nextera-P5-AD with Hi-C-AD1-Block (22 °C, ≥20 h), AMPure, 9-cycle
  Q5 pre-PCR with Nextera-P5-pre-primer / Trueseq-P7-pre-P-S, TDE1 tagmentation,
  7-cycle indexed Nextera P5 / TruSeq P7 PCR.
- 🟢 RNA beads: TdT with CoCl2 and 0.95 µM dGTP / 0.05 µM ddGTP, so the tail is
  poly(G) (the earlier model said poly(C)); Post-TdT PCR 13–18 cycles with
  Post-TdT-PolyC12-S and Post-TdT-P5(T7); Nextera XT tagmentation; 12-cycle
  indexed TruSeq P5 / Nextera P7 PCR.
- 🟢 Table 1a prints `/5Phos/AGA CAG GAC ACT GAC /iBiodT/TC GTT…VN`: the biotin-dT
  is a base, so the RT primer is {{= gage_seq.RT_DT }}; the earlier model dropped it.
- 🟢 Table 1e A1 `AD-1_F_1` = `/5Phos/TNAA AGTACACG CTG TCT CTT ATA CAC ATC TGA` and
  Table 1g A1 `AD-1_R_1` = `/5Phos/TGACTTG TAGTG CAG CGTGTACT T`. The earlier model
  used the RNA Round1 barcode and dropped the 3′ G of the forward strand; corrected
  and drawn as a Scene duplex with the Y arms unpaired. AD-2 (Tables 1h/1i) pairs
  through `TATCAGCGAGATCGG` and leaves `CAAGTCA` to anneal to the AD-1_R `TGACTTG`.
- 🟢 Read structure from Methods line 177: DNA read 2 = BC2 (1–8), 15-nt joint
  (9–23), BC1, then 35 nt trimmed; RNA read 1 = N10 UMI, BC2, 30-nt linker
  (19–48), BC1 (49–56), round-1 arm (57–71). The final constructs were reordered
  to this layout (the earlier DNA construct placed both barcodes on the P5 side).
- 🟡 The UMI is the N10 of the Round2 oligo (Table 1d), not the RT primer.
- Open: the Methods trim 12 nt from DNA read 1, matching an AD-1 BC1 + A + TNA end
  rather than a Tn5 end. The model draws the Tn5 end after tagmentation and does
  not decide which P5 end read 1 starts from.
- Open: Hi-C-AD1-Block sequence and the annealed Nextera-P5-AD duplex are not in
  the inspected tables; that ligation is drawn as an unavailable-geometry state.
- Open: the published Methods call the round-1 plate "CARE-seq 1st-round adaptors"
  and print "500 mM dNTP, 10 mM Biotinylated RT primers"; neither changes the
  modelled molecules and the concentrations are not modelled.
