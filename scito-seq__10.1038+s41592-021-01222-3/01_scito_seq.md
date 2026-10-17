# SCITO-seq

**Defining paper:** Hwang et al., *Nature Methods* (2021), doi:10.1038/s41592-021-01222-3. The open author manuscript and methods were read from PMCID PMC8643207.

## Molecular boundary

🟢 Each antibody carries a reusable 20-nt handle. A hybridized splint supplies its complement, a TruSeq Read 2 handle, antibody identity, a 10-nt staining-pool barcode and the 10x feature-capture sequence.

🟢 Cells are stained in separate pools and combined before high-loading droplet capture. Gel beads add a droplet barcode and UMI; antibody-plus-pool identity and droplet identity resolve cells combinatorially.

## Evidence boundary

🟢 The source publishes multiple splint *architectures*, including
commercial-antibody and RNA-comodality variants. For the original design it specifies,
in order, a 20-nt antibody-handle complement, 34-nt TruSeq Read 2 region, 10-nt
compound antibody/pool barcode and 22-nt 10x feature-capture region.

🟡 We rechecked the defining article and Methods, the open author manuscript,
Supplementary Tables 1–4, the bioRxiv supplement and the authors' linked SCITO-seq
analysis repository. They do not print a universal set of antibody-handle, antibody-
barcode or pool-barcode bases. Those fields are deliberately antibody/pool specific;
the final capture/indexing chemistry is inherited from the commercial 10x Chromium 3′
v3 Feature Barcoding workflow. The page therefore keeps those regions typed by role
instead of inventing a representative sequence.

## Bounded source pass (2026-10-17, round breadth-first-17)

Read: PMC8643207 author-manuscript JATS (`s41592-021-01222-3_PMC8643207.xml.txt`)
Results lines 113, 135 and 143–145; Methods lines 189–205 (conjugation, hybridization
design, saturation), 207–209 (10x run and sequencing), 260–267 (TS-C and comodality
variants); Fig. 1b legend line 608; MOESM1 table headings (Supplementary Tables 1–4:
antibody panels and cost only). The bundle scraper finds no oligo sequences.

🟢 Each clone is ThunderLink-conjugated (1 antibody : 3 oligos) to its own 20-bp
3′-amine-modified, HPLC-purified handle oligo; one conjugation per clone is reused
across pools.

🟢 1 µg of each antibody is pooled with 1 µl of 1 µM secondary (splint) oligo, 15 min at
room temperature, then purified on an Amicon 50K MWCO column to remove excess free
oligonucleotides. The model now represents this separation explicitly.

🟢 The libraries follow the 10x 3′ v3 Feature Barcoding (Cell Surface Protein) guide;
sequencing was NovaSeq 6000, Read 1 28, i7 8 and Read 2 98 cycles. Read 2 can be
shortened according to the antibody-plus-pool barcode length.

🟡 Methods line 195 labels the 10-bp field "pool barcode", while the Results (line 113)
and the preceding Methods sentence describe a compound antibody-plus-pool barcode. The
model keeps the compound Ab+PBC reading; this is wording, not a sequence conflict.

Open: the 10x guide revision for the drawn 3′ v3 route (the paper cites CG000185 Rev D,
page 52, only for the TS-C variant; local copies are Rev A/B only); Fig. 1b and Extended
Data Fig. 9b artwork (not in the source bundle); TS-C 15-bp, scATAC splint and scifi
comodality bridge-oligo variants are not drawn.
