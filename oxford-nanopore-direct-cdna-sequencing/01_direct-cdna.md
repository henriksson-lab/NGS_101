# Oxford Nanopore Direct cDNA Sequencing (SQK-DCS109)

## 1. What it is

The historical SQK-DCS109 workflow makes strand-switched duplex cDNA without PCR, dA-tails it, then ligates a motor-loaded dT-tailed sequencing adapter.

## 2. Sources and evidence

- 🟢 Oxford Nanopore, *Chemistry Technical Document* CHTD_500_v1_revAV_24July2026: https://nanoporetech.com/document/chemistry-technical-document. Section 12 "RNA and cDNA sequencing kits" lists SQK-DCS109 among "previously available" kits; the "Direct cDNA Sequencing Kit" paragraph (local `.txt` lines 1072–1074) is the only kit-specific passage read.
- 🟢 Oxford Nanopore product brochure BR_1009(EN)_V5_30Nov2020, "Which RNA kit?" table (printed p15; PDF page 8 rendered): SQK-DCS109 for "Full-length transcripts without PCR bias", 275 min prep, 100 ng poly(A)+ RNA input, PCR not required, RT required, multiplexing with the Native Barcoding Expansion pack. https://store.nanoporetech.jp/sites/default/files/s3/literature/product-brochure.pdf
- 🔴 The SQK-DCS109 kit protocol (reagent names, volumes, conditions) was not obtained: guessed nanoporetech.com document URLs returned 404 on 2026-10-17, and the legacy protocol is behind the Nanopore Community login.

## 3. What the vendor states (verbatim facts)

- 🟢 PCR-free; input 100 ng poly(A)-tailed RNA, or 70–200 ng of already prepared cDNA (tech doc line 1073). This is far more than the 1 ng poly(A)+ RNA listed for SQK-PCS109 (brochure p15).
- 🟢 Order (tech doc line 1074): anneal the reverse transcription primer; add the strand-switching primer; reverse transcription with strand switching; remove the RNA template with an RNase cocktail; second-strand synthesis; the double-stranded cDNA "is treated similarly to gDNA"; the fragment ends are dA-tailed; a dT-tailed sequencing adapter is ligated.
- 🟡 Poly(A) priming of the RT primer is inferred from the poly(A)-tailed RNA input; the DCS109 paragraph does not describe the primer.
- 🟡 End repair alongside dA tailing is inferred from "treated similarly to gDNA"; the paragraph names only dA tailing.
- 🔴 Kit RT-primer, strand-switching-primer, sequencing-adapter and motor sequences are not given in the inspected passages. The tech doc's "Sequencing adapter sequences" table (lines 2070–2105) is introduced as "the adapter sequences for our current kits" (LA, NA, RA, RTP, SSPII, CRTA, RAP T, cPRM, RBA). None is assigned to SQK-DCS109, so the model uses none of them.

## 4. Molecular path

The cDNA duplex is treated like a genomic-DNA ligation library after second-strand synthesis. The A:T overhang junction installs the sequencing adapter, rather than a PCR-added rapid tag. The ligated adapter end enters the pore first.

## 5. Open questions

- The kit protocol: RT primer and strand-switching primer identity (bases, modifications, any rG 3′ end), the RT enzyme, the RNase cocktail, the second-strand polymerase, end-prep and ligation reagents, cleanups.
- Sequencing-adapter identity and bases for this kit, and which strand carries the motor.

## 6. Follow-up pass 2026-10-17 (round deferred-followup-02)

Routes tried for the SQK-DCS109 kit protocol, none successful 🔴:

- nanoporetech.com document sitemap (`/__sitemap__/document.xml`, 380 entries): has no
  DCS109 document. Legacy kits such as SQK-RNA002 are still listed, but SQK-DCS109 is not.
- Guessed `nanoporetech.com/document/…dcs109…` and product/store pages return 404.
  `community.nanoporetech.com/…/direct-cdna-sequencing-sqk-dcs109` (with and without a
  revision) returns "page doesn't exist any longer".
- Internet Archive Wayback/CDX was temporarily offline on 2026-10-17. Retry it with
  `community.nanoporetech.com/protocols/direct-cdna-sequencing-sqk-dcs109`.
- protocols.io API search requires a token; not pursued.

Lead read but **not used**: the successor protocol *Ligation sequencing V14 – Direct cDNA
sequencing (SQK-LSK114)*, DCS_9187_v114_revK_24Feb2026
(https://nanoporetech.com/document/ligation-sequencing-v14-direct-cdna-sequencing; local
`successor_…html.txt` lines 164–345, 600–720, 912, 964–994). It names Maxima H Minus RT,
RNaseOUT, the RNase Cocktail, LongAmp Taq for the second strand and NEBNext Ultra II End
Prep. Its user-supplied VN primer (`/5phos/ACTTGCCTGTCGCTCTATCTTC…T20VN`, line 267),
strand-switching primer and PR2 primer belong to a different kit, and the document never
mentions SQK-DCS109. A current protocol for another kit does not identify the DCS109
reagents or sequences, so the model and page are unchanged. Scraper run on it: 3 hits
(VN, strand-switching, PR2 primers), all successor-kit-specific.
