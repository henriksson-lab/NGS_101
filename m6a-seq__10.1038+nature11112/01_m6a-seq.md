# m6A-seq / MeRIP-seq — research notes

Defining source: Dominissini et al., 2012, “Topology of the human and mouse m6A RNA methylomes revealed by m6A-seq”, *Nature* 485:201–206, [doi:10.1038/nature11112](https://doi.org/10.1038/nature11112). Closed access (no PMC copy); read from the author-hosted PDF of the published article.

## Sources read

- 🟢 Author-hosted PDF of the published *Nature* article (<https://www.weizmann.ac.il/molgen/Sorek/files/publications/Dominissini_Nature_2012.pdf>, text twin `Dominissini_Nature_2012.pdf.txt`): Results first section (lines 69–79), Methods Summary (lines 768–779), Methods “RNA preparation” and “m6A-seq” (lines 957–980), “Calculation of coverage” (lines 1031–1040). Bounded pass 2026-10-17.
- 🟢 Publisher Supplementary Information PDF (MOESM81): Supplementary Note 1a (line 1225), elution only with excess free m6A nucleoside. The supplement figures and tables do not supply a library oligo table.
- Oligo scraper run on both texts: the only hits are the METTL3 siRNA and the biotinylated RSV RNA-affinity baits (Methods lines 978–991), neither of which is part of the library.

## Boundary

🟢 m6A-seq/MeRIP-seq enriches a population of RNA fragments; the modification does not become a barcode or a deterministic base change. As the Results note, fragmentation before IP makes the method blind to single-transcript patterns; peak position is broader than single-nucleotide resolution.

## Step by step (Methods)

1. 🟢 Total RNA (HepG2 or mouse liver PerfectPure kits, 5 Prime; human brain RNA from Biochain) is DNase treated; poly(A)+ RNA is enriched with one round of the GenElute mRNA miniprep kit (Sigma-Aldrich).
2. 🟢 Chemical fragmentation to ~100 nt: 5 min at 94 °C in 10 mM ZnCl2, 10 mM Tris-HCl pH 7; stopped with 0.05 M EDTA, ethanol precipitated, resuspended in H2O. Bioanalyzer fragment range 90–150 nt (“Calculation of coverage”).
3. 🟢 Fragmented RNA (mRNA or total RNA) is incubated 2 h at 4 °C with affinity-purified anti-m6A polyclonal antibody (Synaptic Systems) in IPP buffer (150 mM NaCl, 0.1% NP-40, 10 mM Tris-HCl pH 7.4), then 2 h at 4 °C with protein-A beads (Repligen), and washed extensively. Input control fragments are sequenced alongside.
4. 🟢 Bound RNA is eluted competitively with free N6-methyladenosine (Sigma-Aldrich) in IPP buffer and ethanol precipitated.
5. 🟢 IP and input RNA are used “for library generation with mRNA sequencing kit (Illumina)”.
6. 🟢 Sequencing on Illumina GAIIx, 10 pM template per sample for cluster generation, genomic sequencing kit V2. The coverage analysis extends reads to 100 nt from a “sequenced 36-nucleotide stretch”.

The extracted text renders the µ sign as “m”, so quantities (RNA input, antibody, elution concentration) are left out of the model until checked against the rendered PDF.

## Final library and read layout

🔴 The Methods give no kit catalogue number and no reverse-transcription primer, adapter or PCR oligo sequence. The model therefore treats library construction as a warned kit boundary and draws the final library ends as role tokens. A 2026-10-17 model correction removed the dual-indexed TruSeq shell and 2 × 50 reads drawn earlier: neither the source nor a kit original supports them.

🟡 Read 1 is 36 nt (from the coverage Methods). A second read is not described, and the read's orientation relative to the RNA is not stated. Sample indexing is not described (“10 pM template per sample”).

## Open questions

- 🔴 The exact 2012 Illumina mRNA sequencing kit generation, its library steps (the kit normally includes its own fragmentation, which the Methods do not discuss for pre-fragmented IP RNA) and its adapter/primer sequences. Next source to read: the historical Illumina kit guide.
- The later author protocol, [Dominissini et al. 2013](https://doi.org/10.1038/nprot.2012.148), describes a detailed implementation. It is a later protocol and should not silently determine the 2012 library construct.
- RNA input, antibody amount and m6A elution concentration need checking against the rendered PDF.
