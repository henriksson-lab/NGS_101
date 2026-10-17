# TeloPCR-seq

## Published source discrepancy

🟢 Within the published FEBS Letters article (PMC5561429 author manuscript), the Methods (`1873-3468.12444_PMC5561429.xml.txt` line 48) remove free primers with an **Amicon Ultra 0.5 mL centrifugal filter (Millipore)** followed by the Qiagen PCR purification kit, whereas the Figure 1A legend (line 374) says the products were purified "using a **Microcon concentrator**". The model and page follow the Methods; the filter cut-off is not stated, so which device (or both) was used remains unresolved.

Version note: the superseded author bioRxiv v1 Methods (`047837_v1.full.pdf.txt` lines 205–209) wrote "0.25 pM" anchor where the published Methods say 0.25 picomoles per 20 µL. That is an ordinary preprint-to-publication revision, not a published discrepancy; the published value governs (round 2026-10-17-breadth-first-20 correction of earlier framing).

## 1. What it is

TeloPCR-seq attaches an anchor to native *S. pombe* chromosome termini, amplifies each end with a subtelomeric primer, and reads the resulting telomere amplicons as PacBio circular libraries.

## 2. Sources read

- 🟢 Bennett et al. (2016), [doi:10.1002/1873-3468.12444](https://doi.org/10.1002/1873-3468.12444), PMC5561429 full text, `/husky/henriksson/NGS101_papers/sources/telopcr-seq__10.1002+1873-3468.12444/1873-3468.12444_PMC5561429.xml.txt` lines 48 (Telo-PCR) and 56 (sequencing pipeline), Materials and Methods; Figure 1 legend line 374; Table 1 lines 392–415. The article itself prints the anchor, its terminal modifications, both subtelomeric primers and the anchor PCR primer.
- Superseded (not used for values): the same four oligos in the authors' bioRxiv preprint, [doi:10.1101/047837](https://doi.org/10.1101/047837), `/husky/henriksson/NGS101_papers/sources/telopcr-seq__10.1002+1873-3468.12444/047837_v1.full.pdf.txt` lines 205 and 213–214 — an independently extracted second copy of the Methods paragraph.
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), method overview.

## 3. Construct and reaction order

1. 🟢 CircLigase II ssDNA ligase (Epicentre) attaches the anchor to telomere ends: 20 µL with 10 µg phenol/chloroform genomic DNA, 0.25 picomoles anchor, kit buffers plus betaine; 60 °C 2 h, 80 °C 10 min (Methods line 48).
2. 🟢 PCR from 2 µL ligation per 50 µL with subtelomeric primer 1 or 2 and the anchor reverse-complement primer: Taq (Roche, 5% DMSO; 95 °C 3 min, 45 × [95 °C 30 s, 52 °C 30 s, 65 °C 40 s], 65 °C 5 min) for primer 1, iProof GC buffer (98 °C 3 min, 45 × [98 °C 10 s, 60.5 °C 15 s, 72 °C 45 s], 72 °C 5 min) for primer 2; LongAmp Taq 2X is also named. Table 1 libraries used Taq/1 or iProof/2.
3. 🟢 Six 50 µL PCRs are pooled, concentrated and buffer-exchanged on an Amicon Ultra 0.5 mL filter to remove free primers, diluted with 100 µL water and purified with the Qiagen PCR purification kit.
4. 🟢 Libraries were generated and sequenced per the manufacturer's ~250 bp protocol (PacBio RS) at the Yale Center for Genome Analysis (line 56). 🔴 The kit revision, end processing and hairpin ligating ends are not stated; the model draws this as an unresolved boundary rather than assuming dA-tailing/T-overhang hairpins (round 2026-10-17-breadth-first-20 correction).
5. 🟢 Circular consensus reads are filtered computationally: both amplifying primers present (≤2 mismatches) and at least one GGTTAC; the subtelomeric primer is trimmed before length measurement.

The defining article gives these ordered oligos 🟢 (5′→3′):

| Role | Oligo as printed |
|---|---|
| Ligated anchor | `5′-[Phos]AGGGTTAATAAGCGGCCGCGTCG[AmC3]-3′` |
| Subtelomeric primer 1 | `5′-CGAGGCTGCGGGTTACAAGGTTACG-3′` |
| Subtelomeric primer 2 | `5′-GGGAATTTAGGAAGTGCGGTAAGTTG-3′` |
| Anchor PCR primer | `5′-CGACGCGGCCGCTTATTAACCCT-3′` |

The Methods call the anchor PCR primer "the reverse complement of the anchor sequence" and print it; the model therefore derives it as `revcomp(ANCHOR)` rather than retyping it, and the printed string is what that computation yields 🟢. `[Phos]` is the ligating 5′ phosphate; `[AmC3]` blocks the anchor 3′ end. The model shows primer 1 as one representative subtelomeric boundary. The PacBio hairpin oligos and sequencing-primer site remain undisclosed 🔴; a modern SMRTbell kit sequence cannot establish which historical kit revision was used.

## 4. Final library and sequencing

The selected telomere amplicon lies between two hairpin adapters. A PacBio sequencing primer binds a hairpin primer site and the polymerase makes repeated passes around the closed template.
