# BID-seq

Defining original: https://doi.org/10.1038/s41587-022-01505-w

Dai et al., *Quantitative sequencing using BID-seq uncovers abundant pseudouridines in mammalian mRNA at base resolution* (online2022; issue2023).

🟢 Published original `s41587-022-01505-w_PMC10017504.xml.txt`, Results81–113, limitations209, Methods231–237 and processing281; publisher SI `MOESM1_ESM.pdf.txt` FigS1 caption264–266. PMC identifies non-manuscript content. The matching SI title/DOI and original caption were checked; full graphical chemistry panels and table workbooks were not exhaustively reviewed. Full-bundle scraper run; adapter/RT-primer hits verified in actual Methods233–237.

## Molecular distinction

🟢 The optimized bisulfite treatment produces a Ψ-BS RNA adduct. SuperScriptIV can skip a complementary base while copying this adduct. The Ψ residue is still present in RNA: cDNA deletion is a copying event, not a physical RNA deletion. A nondeletion copy of the same adduct is permitted. Results101 discusses motif-dependent deletion frequencies, including high-background contexts; the model assigns no rate or complete conversion probability. Adjacent uridines can make deletion position ambiguous (limitations209).

## Source oligos and order

🟢 Methods233: 3′adapter `5′App-NNNNNATCACGAGATCGGAAGAGCACACGTCT-3SpC3`, with `ATCACG` explicitly an inline barcode. The 3′C3 block is a non-nucleotide spacer, unlike m6ACE's terminal ddC. Methods235: 5′RNA adapter `GUUCAGAGUUCUACAGUCCGACGAUCNNNNN`. Methods237: SR RT primer `AGACGTGTGCTCTTCCGATCT`. Canonical adapter motifs are imported; activation, random windows, barcode and spacer identity remain local. The two random windows receive different identifiers. Their UMI/deduplication interpretation is inferred from the source's random-barcode removal/deduplication analysis, rather than an assumed consequence of N5 alone.

🟢 Poly(A) RNA10–20ng; fragmentation70°C14min; sequential PNK end repair/phosphorylation; 3′ligation with truncatedKQ RNA ligase2,25°C2h then16°C10h; 5′deadenylase/RecJf removal of excess adapters and column purification; 5′RNA ligase1 ligation25°C8h; column cleanup; split1.5µl input/8.5µl treated. Treat with2.4MNa2SO3/0.36MNaHSO3 at70°C3h; column desulfonation1h and washes/elution; SSIV50°C1h; RNaseH37°C20min; heat70°C5min; DNA column;15cycle PCR using SR Primer/NEBNext indexed primers;3.5%agarose purification; NovaSeq6000SE100.

## Local scope and open work

The graph starts after source adapter ligations and split, with separate treated/input constructors. It covers treatment/desulfonation, RT, RNaseH and actual DNA purification. RNaseH products remain molecular material until column separation. Adapters have an explicit inventory but ligation-ready geometry, exact historical NEBNext PCR/index set and acquired read-cycle plan are not assembled. Notes status; no finished schematic or confirmed missing-information claim.

Ancestry suggestion: Results87 explicitly derives the Ψ-BS-adduct/deletion principle from RBS-seq (refs15–16), while modifying conditions and removing the conventional cytosine-deamination tradeoff. This establishes component strategy, not inherited library adapters or quantitative equivalence. Detailed2024 procedure DOI10.1038/s41596-023-00917-5 is a deferred successor source, not a replacement for the defining2022 publication.
