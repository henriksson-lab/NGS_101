# icSHAPE — research notes

Primary protocol: Flynn *et al.*, *Nature Protocols* (2016), [10.1038/nprot.2016.011](https://doi.org/10.1038/nprot.2016.011) ([open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC4896316/)). Defining application: Spitale *et al.*, *Nature* (2015), [10.1038/nature14263](https://doi.org/10.1038/nature14263).

## 1. What it is

In vivo click selective 2′-hydroxyl acylation and profiling experiment (icSHAPE) turns RNA flexibility into enriched reverse-transcription stops. NAI-N3 acylates accessible ribose 2′-hydroxyls and also supplies an azide for copper-free attachment of biotin.

## Sources read

- 🟢 The full protocol was read through reagent lists, steps 18–68, anticipated results and sequencing instructions.
- The source fetch command is `python3 tools/get_sources.py --doi 10.1038/nprot.2016.011`; publisher material remains external.

## Key oligos

- 🟢 Treated-sample 3′ adapter as printed: `/5rApp/AGATCGGAAGAGCGGTTCAG/3ddC/`. T4 RNL2 truncated K227Q ligates this blocked, preadenylated adapter to repaired RNA 3′ ends.
- 🟢 RT-primer-1 as printed: `/5phos/DDDNNAACCNNNNAGATCGGAAGAGCGTCGTGGA/iSp18/GGATCC/iSp18/TACTGAACCGC`.
- 🟢 The two random regions total 9 nt; `AACC` is the experimental barcode. The primer architecture supports later intramolecular cDNA circularization.

## Construct-changing steps

1. 🟢 Modify living-cell RNA with NAI-N3; process a DMSO background in parallel.
2. 🟢 Click DIBO-biotin to the azide, fragment RNA to roughly 50–100 nt and repair its 3′ ends with T4 PNK.
3. 🟢 Ligate the 3′ RNA adapter, size-select, then reverse-transcribe with the barcoded RT primer.
4. 🟢 Streptavidin selection occurs after RT. RNases release cDNA from selected RNA–cDNA hybrids.
5. 🟢 Size-select cDNA, circularize it with CircLigase II, and PCR-amplify with Solexa primers.

## Final library and read layout

- 🟢 Standard single-end Illumina conditions are used; 75 bases suffice, and the normal external 6/8-base index read is omitted because the experimental barcode is inline.
- 🟢 The first 13 read bases are the 4-base experimental barcode plus 9 random bases. After trimming, the first mapped transcript base marks the RT stop; analysis shifts one nucleotide to identify the modified base.
- 🟡 Exact Solexa primer binding geometry is left role-labelled on the public page because the accessible protocol prints the PCR primer names separately from a fully assembled top-strand library.
