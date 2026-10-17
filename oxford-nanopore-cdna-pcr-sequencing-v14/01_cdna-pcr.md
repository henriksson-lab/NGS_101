# Oxford Nanopore cDNA-PCR Sequencing V14 (SQK-PCS114)

## 1. What it is

SQK-PCS114 converts RNA to full-length strand-switched cDNA, incorporates a UMI, PCR-amplifies molecules with rapid-attachment tags, and attaches a motor-loaded Rapid Adapter for nanopore sequencing.

## 2. Sources and evidence

- 🟢 Oxford Nanopore, *cDNA-PCR Sequencing V14 (SQK-PCS114)*, document version
  PCS_9200_v114_revF_31Mar2026 (local text line 249):
  https://nanoporetech.com/document/pcr-cdna-sequencing-v14-sqk-pcs114. Rev F says
  the method was “updated extensively” in March 2026 with new third-party reagents
  and PCR conditions; earlier revisions are a different recipe and are not modelled.
  Sections read: overview lines 254–288, kit/third-party reagent lists 312–339 and
  488–507, RT and strand switching 640–1000, PCR 1040–1283, adapter addition 1336–1392.
- 🟢 Oxford Nanopore, *Chemistry Technical Document*, CHTD_500_v1_revAV_24July2026
  (line 164): cDNA-PCR Sequencing Kit V14 workflow lines 467–474; section 15
  adapter sequences lines 2069–2112 (CRTA, RTP, SSPII, cPRM forward and reverse,
  Rapid Adapter / RAP T top strand); section 16 reagent roles from line 2119.
- The oligo scraper was run on both documents; `--find` locates the SSPII handle in
  CHTD lines 2095, 2110 and 2115.

## 3. Oligos (verbatim, CHTD section 15)

- 🟢 CRTA: `5'-CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG-3'`. The protocol calls
  CRTA “a double stranded adapter with a poly(T) overhang” (line 676); the overhang
  strand is not printed.
- 🟢 RTP: `5'-CTTGCCTGTCGCTCTATCTTCAGAGGAG-3'`
- 🟢 SSPII: `5'-TTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTmGmGmG-3'`
- 🟢 cPRM forward: `5'-ATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGACTTGCCTGTCGCTCTATCTTC-3'`
- 🟢 cPRM reverse: `5'-ATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGC-3'`
- 🟢 Rapid Adapter (RA) top strand / RAP T: `5'-CCTGTACTTCGTTCAGTTACGTATTGCT-3'`

Model correction (2026-10-17): the shared `nanopore_library.CDNA_PRIMER_REVERSE`
omits the 7-nt `CTTTGTG` after the common 32-nt block. That is a repository
transcription error, not a vendor discrepancy; the protocol module now declares both
cPRM sequences from CHTD line 2108/2110, and the source guard checks them.

🟡 Derived interlocks (computed in `cdna_pcr.py`): RTP is fully complementary to the
CRTA 3′ {{= cdna_pcr.CRTA_PAIRED }} nt; the cPRM forward 3′ {{= cdna_pcr.FORWARD_ANNEAL }} nt
equal the RTP 5′ end and the cPRM reverse 3′ {{= cdna_pcr.REVERSE_ANNEAL }} nt equal the
SSPII handle. Both cPRM sequences share the same {{= len(cdna_pcr.CPRM_TAG) }}-nt 5′ tag,
so both PCR ends carry the same attachment tag.

## 4. Molecular path (Rev F)

1. 🟢 10 ng poly(A)+ RNA or 500 ng total RNA + CRTA + Annealing Buffer, 60 °C 5 min,
   then 5 min at room temperature (lines 642–687).
2. 🟢 T3 DNA Ligase in StickTogether buffer with RNaseOUT, 10 min RT.
3. 🟢 Lambda Exonuclease + USER, 37 °C 5 min: “digest the bottom strand of the ligated
   CRTA so that the RT Primer (RTP) can bind” (line 747). 🟡 That implies a uracil-
   containing, 5′-phosphorylated bottom strand; neither is printed.
4. 🟢 RNAClean XP (42 µl) with two Short Fragment Buffer washes; elute 12 µl.
5. 🟢 RTP + dNTPs, 5 min RT; Maxima H Minus 5× RT buffer, RNaseOUT, SSPII; 42 °C 2 min;
   add Maxima H Minus RT; 42 °C 30 min, 85 °C 5 min (lines 835–935).
6. 🟢 AMPure XP (40 µl), two 80% ethanol washes; elute 22.5 µl.
7. 🟢 Four 50 µl PCRs, each 5 µl cDNA with cPRM and TaKaRa Ex Premier: 5 cycles of
   98 °C 10 s / 50 °C 30 s / 68 °C 3 min, then 17 cycles of 98 °C 10 s / 68 °C 3 min.
8. 🟢 Thermolabile Exonuclease I (1 µl per PCR) removes excess primers; 37 °C 5 min,
   80 °C 2 min. Pool the four reactions (204 µl); AMPure XP 140 µl; elute 15 µl EB.
9. 🟢 Take 50 fmol, add 1 µl RA diluted in ADB (1.5 + 3.5 µl), 5 min RT (lines 1283–1392).

🟡 The drawn strand-switched cDNA places RTP, the copied CRTA 5′ region, poly(T), the
transcript, three non-templated C's, the complemented UMI and the SSPII-handle
complement in that order. The vendor writes that SSPII pairs to deoxycytidines “at the
5′ end of the first cDNA strand” (line 888); the model places them at the growing 3′
end opposite the RNA 5′ end, and the page shows this as a highlight.

This protocol is kept separate from Direct cDNA because PCR is both the enrichment/amplification step and the means of installing the Rapid Adapter attachment boundary.

## 5. Open questions

🔴 Not established by the inspected documents: the CRTA overhang strand (bases,
uracils, phosphorylation); the rapid-attachment chemistry carried by the cPRM 5′ tags;
the complete motor-bearing Rapid Adapter assembly and which or how many ends are
adapted. Still to inspect: an ONT component-level source for RA/cPRM attachment, and the
pre-March-2026 SQK-PCS114 revisions if a dataset used them. No absence claim is accepted.
