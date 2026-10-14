# Double-stranded ancient DNA with partial UDG

## 1. What it is

This double-stranded ancient-DNA workflow removes most internal deamination damage while deliberately retaining terminal damage, then installs short molecular barcodes before target capture or sequencing.

## 2. Authoritative source

🟢 Rohland et al., *Partial uracil–DNA–glycosylase treatment for screening of ancient DNA*, Philosophical Transactions B (2015), DOI `10.1098/rstb.2013.0624`.

## 3. Construct-changing operations

🟢 Limited UDG and Endonuclease VIII treatment cleaves internal uracil sites while terminal uracils survive. End repair and ligation attach a unique pair of short barcoded adapters; later PCR completes the platform arms.

## 4. Adapter and read layout

🟢 USER acts for 3 hours before T4 polymerase and T4 polynucleotide kinase are added. This ordering removes internal uracils while inefficient removal at the termini preserves a reduced damage signal.

🟢 Two partially double-stranded adapters carry independent 7-nt molecular barcodes. Table 2 prints `ATCGATT` / `GACTTAT` as one representative P5/P7 pair; the barcodes are the first seven bases of Read 1 and Read 2 and are trimmed before alignment.

🟢 PreHyb PCR leaves short adapter arms for capture; a later indexing PCR completes i5 and i7 arms for sequencing.

🔴 The authors place the exact adapter, PreHyb and index-primer oligos in Supplementary
Table S1. The PMC article, PMC Cloud bundle and Europe PMC record do not contain that
supplement, and the Royal Society supplementary-file endpoint currently returns an anti-bot
page rather than the PDF. Later papers reproduce parts of this chemistry, but those secondary
transcriptions are not promoted to exact constants here. The unresolved bases therefore mean
"original source file presently unavailable", not "the paper omitted the sequences".
