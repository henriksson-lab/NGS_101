# ChAIR

**Defining paper:** Chai et al., "Tri-omic single-cell mapping of the 3D epigenome and
transcriptome in whole mouse brains throughout the lifespan", *Nature Methods* (2025),
doi:10.1038/s41592-025-02658-7 (author manuscript PMC12955729). The bioRxiv preprint
doi:10.1101/2024.05.03.592322 is superseded and not used.

> 🟢 verbatim from a source read for this note · 🟡 derived/inferred · 🔴 not published or
> not read. The model is `tools/chem.py`.

## Sources read (2026-10-17 bounded pass)

All in `/husky/henriksson/NGS101_papers/sources/chair__10.1038+s41592-025-02658-7/`.

| File | What | Sections read |
|---|---|---|
| `s41592-025-02658-7_PMC12955729.xml.txt` | published article (author manuscript JATS) | Methods "ChAIR" line 179; "ChAIR data processing" line 187 |
| `s41592-025-02658-7_41592_2025_2658_MOESM1_ESM.pdf.txt` | Springer Supplementary Information, "A step-by-step ChAIR protocol" (dated 9.30.2024; fetched 2026-10-17) | Materials list; steps 1–82; Addendum "Bridge linker preparation" |
| `10x-chromium-single-cell-multiome-atac-gene-expression/10x_CG000338_…_RevG.pdf.pdf.txt` | the 10x guide ChAIR names (no revision cited) | Appendix "Sequences", p. 80 (ATAC and GEX library structures) |

The PMC NIHMS supplement rows in `MANIFEST.tsv` remain behind the PMC download gate; the
Springer Supplementary Information was read instead, and their identity with it was not
compared 🔴.

## Oligos

🟢 Bridge linker, printed identically in Methods (line 179) and the SI Addendum:

```
top     5'-/5Phos/CGTGATATT/iBIOdT/CACGACTCT-3'
bottom  5'-/5Phos/GAGTCGTGAAATATCACGT-3'
```

`/iBIOdT/` is the biotinylated T at position 10 of the top strand, as printed. Annealed (SI: 95 °C 2 min then a stepped 0.1 °C/s ramp to 4 °C, top:bottom
ratio chosen on a TBE gel), the strands form an 18-bp duplex with a single 3′-T overhang
on each strand 🟡 (computed by `chem.bridge_linker_ends()`, which raises otherwise). That
matches the preceding Klenow (3′→5′ exo−) dA-tailing: the linker joins two dA-tailed
restriction ends, and the sealed top strand reads dA of fragment A + the top oligo + fragment B in one
of the two linker orientations (`chem.contact_junction()` checks that this equals the
reverse complement of dA + the bottom oligo) 🟡.

All other oligos are 10x kit reagents (ATAC Enzyme B adapters, Multiome Gel Beads,
Pre-Amp Primers, SI-PCR Primer B, Sample Index N, cDNA Primers, GEX adaptors) whose
sequences ChAIR does not print 🔴.

## Steps (Methods line 179 and SI steps; conditions 🟢)

1. FA (2%) then EGS (2 mM) crosslinking; lysis in 0.1% SDS FA buffer (omitted for nuclei);
   permeabilization in 0.1% SDS; quench with Triton X-100.
2. In-situ digestion with AluI, or AluI + HpyCH4V (SI step 8 uses both), CutSmart, 37 °C ≥2 h.
3. dA-tailing: Klenow (3′→5′ exo−), dATP, BSA, 37 °C 1 h; stop at 65 °C 20 min.
4. In-situ ligation of 3 µl 2 ng/µl bridge linker with T4 DNA ligase in NEB Quick ligation
   buffer, RT ≥4 h or 16 °C overnight.
5. Wash (ATAC-RSB with RNase inhibitor), resuspend in 10x Diluted Nuclei Buffer, brief
   sonication, 20 µm strainer; 12,000 nuclei tagmented with ATAC Buffer B + ATAC Enzyme B,
   37 °C 1 h.
6. Chip J / Chromium Controller GEMs "strictly following" CG000338; GEM incubation 37 °C
   45 min, 25 °C 30 min.
7. Recovery Agent; proteinase K 65 °C 2 h (decrosslinking); Dynabeads cleanup.
8. Streptavidin capture on M-280 beads (pre-blocked with iBlock and 500 ng genomic DNA):
   beads = **DNA component** (washed 5× 0.5% SDS/2× SSC); supernatant = **cDNA component**
   (1.8× SPRI).
9. DNA: pre-amplification on beads (NEBNext HF 2×, Pre-Amp Primers, 72 °C 5 min + 7 cycles),
   beads discarded, 1.6× SPRI; sample-index PCR (SI-PCR Primer B + Sample Index N, 7
   cycles); 1× SPRI, then 0.55×/0.8× double-sided selection → **ChAIR-DNA** library.
10. cDNA: same pre-amplification, 1.6× SPRI; cDNA amplification (cDNA Primers, 7 cycles),
    0.6× SPRI; GEX fragmentation, end repair/A-tailing, adaptor ligation, sample-index PCR
    per CG000338 → **ChAIR-RNA** library.

The graph keeps the 10x GEM step and the GEX library construction as warned unmodelled
boundaries.

## Libraries and readouts

There are **two physical libraries**, ChAIR-DNA and ChAIR-RNA 🟢. ChAIR-ATAC and
ChAIR-PET are computational partitions of the ChAIR-DNA reads: ChAIR-PIPE starts with
ChIA-PIPE linker trimming and mapping, then assigns cell barcodes from Cell Ranger ARC
(Methods line 187) 🟢. An earlier revision of this model drew three separately amplified
libraries (RNA, ATAC and PET) from a generic shared recipe; that was a repository
modelling error and has been replaced.

The drawn ChAIR-DNA molecule uses the CG000338 Rev G appendix ATAC shell
(P5 · 16-nt cell barcode · `CGCGTCTG` spacer · Read 1N · insert · Read 2N · 8-nt
Sample Index N · P7) 🟢 for Rev G; ChAIR cites no revision 🔴. Read lengths: "prepared and
sequenced in accordance with the manufacturer's guidelines" 🟢; no ChAIR cycle counts are
assigned 🔴.

## Published discrepancies (shown on the page)

- Permeabilization: Methods "10 µl 0.1% SDS at RT for 2 h" vs SI steps 4–5 "1.5 h"
  (optimize per cell type).
- SDS quench: Methods 37 °C 20 min vs SI step 7 15 min.
- Declumping sonication: Methods "1 second" vs SI step 27 "43% for 1.5 s".

The model records none of these durations.

## Open questions

- 🔴 Kit revision of CG000338 used; kit oligo sequences and the selected sample indices.
- 🔴 Fate of linker-free (non-biotinylated) DNA left in the M-280 supernatant and of
  half-ligated linker products; not stated in the inspected passages.
- 🔴 Whether the gated PMC NIHMS supplements differ from the Springer Supplementary
  Information.
- 🔴 Read lengths actually acquired for each dataset.
- 🔴 Predecessor ChIATAC (ref. 11) and the Current Protocols ChIA-PET linker source
  (Wang et al. 2021) were not read in this pass.
