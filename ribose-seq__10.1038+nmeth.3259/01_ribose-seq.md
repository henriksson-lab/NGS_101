# Ribose-seq: capture a ribonucleotide-bearing circular fragment

Koh et al. (2015), [Ribose-seq: global mapping of ribonucleotides embedded in genomic DNA](https://doi.org/10.1038/nmeth.3259), and [2019 Addendum](https://doi.org/10.1038/s41592-019-0505-9).

## Capped sources and revision identity

One reader,28October2026; three chemistry originals: article-associated accepted manuscript PMC4686381 (Results and Methods, especially “Ribose-seq library construction”, “DNA sequencing”, analysis); full published2019 Addendum; its corrected publisher Supplementary Information (Figure1/chemistry captions and Table8, physically rendered PDF page20). The article manifest identifies an author manuscript, not the published PDF. Corrected supplement is under the addendum's canonical DOI bundle; no copied source is kept in the checkout. Mandatory scraper was run on both protocol and addendum bundles. The2015 supplement remains historical, superseded for the corrected figures/tables by the2019 PDF.

🟢 **Addendum:** upstream/downstream coordinates in sequence-context plots were reversed. It replaces Fig3 context plots, Supplementary Figs5/6/8 and Table4, adjusts directional text, and changes the reverse-strand offset handling in `nuc_frequencies.py`. Use corrected positions; do not reproduce the original directional claims. This is an analysis correction, not a reversal of the physical capture strand. No whole-version equivalence beyond the correction's stated scope is assumed.

## Ordered preparation

1. 🟢 Extract yeast genomic DNA. Digest overnight37°C with SspI, DraI and EcoRV to500–3000bp, average about1.5kb. Klenow exo− and dATP A-tail fragments30min37°C; spin-column purify.
2. 🟢 T4DNA ligase attaches preannealed Adaptor.L/Adaptor.S overnight15°C. The adapter has one3′T overhang, eight random UMI bases and phosphorylated5′ ends. Adaptor.S is20nt with a terminal3′amine, not21nt: `Am` is a non-nucleotide modification, as Table8's rendered row and footnote establish. The scraper's21-base reading is rejected. AMPure purification follows this and subsequent recovery steps.
3. 🟢 Treat adapter-ligated DNA with0.3MNaOH2h55°C, neutralize and purify. This denatures DNA and cleaves at embedded ribonucleotides. The source describes2′,3′cyclic-phosphate and2′phosphate capture substrates and nonproductive3′phosphate ends; do not treat every cleaved fragment as ligatable or substitute a free3′OH.
4. 🟢 AtRNL with ATP self-ligates compatible rNMP-bearing fragments at30°C1h. The circle retains the rNMP and a residual2′phosphate at the newly formed junction. An AMPure purification follows AtRNL; remaining unclosed fragments are T5exonuclease substrates at37°C2h. Digestion products remain in the tube through Tpt1 treatment, not magically removed by catalysis. Tpt1 with NAD+ removes the junction2′phosphate30°C1h; source-described subsequent AMPure purification supplies the physical separation, without a quantitative removal guarantee.
5. 🟢 PCR uses PCR.2 and one of PCR.1.Index1–4 with Phusion or EconoTaq. These are source alternatives, not interchangeable proof of rNMP bypass. Check by6%PAGE, purify, pool and sequence on MiSeq,50-cycle single-end.

## Exact Table8 library oligos,5′→3′

`P` means phosphate; `Am` means terminal amino group, neither contributes a base. Bases and modifications were checked against the rendered corrected table.

| Name | Sequence |
|---|---|
| Adaptor.L | P-NNNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGGGAGTTCAGACGTGTGCTCTTCCGATCTAGCCAGCGCAGACCGTGAGGT |
| Adaptor.S | P-CCTCACGGTCTGCGCTGGCT-Am |
| PCR.1.Index1 | CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC |
| PCR.1.Index2 | CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC |
| PCR.1.Index3 | CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC |
| PCR.1.Index4 | CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC |
| PCR.2 | AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT |

## Physical chemistry versus sequencing

🟢 Reads begin with eight UMI cycles, then genome-derived sequence. The source removes UMI before alignment and identifies the captured rNMP as the reverse complement of the first aligned read base. A retained ribonucleotide is not chemically replaced by a deoxynucleotide merely because PCR reports its complementary base. For a contiguous rNMP stretch the source states the5′most rNMP is captured; it does not claim independent detection of every base in that stretch.

🟡 Local states declare end chemistry, topology and junction2′phosphate separately. The shared free3′OH ssDNA-circle helper is inappropriate for this AtRNL pathway; no false3′OH parameter is supplied. Exact circular-PCR strand assembly and run-primer acquisition remain unassembled in this capped model. Displaying source-known primers does not itself validate PCR product geometry.

## Ancestry and deferred work

🟢 The article explicitly cites Schutz, Hesselberth and Fields2010, [Capture and sequence analysis of RNAs with terminal2′,3′cyclic phosphates](https://doi.org/10.1261/rna.1934910), for AtRNL/Tpt1 and related cyclic-phosphate capture. This supports enzyme/terminal-processing ancestry; its complete RNA library preparation is not imported. The predecessor original remains unread within this cap.

Deferred: corrected published article Methods; complete linked-original checklist and second reader; source-specific strand-preservation/AtRNL chemistry validation; exact circular-PCR endpoint and sequencing/index primer mapping; study index assignments. No sequence or reaction detail is claimed absent from the publication.
