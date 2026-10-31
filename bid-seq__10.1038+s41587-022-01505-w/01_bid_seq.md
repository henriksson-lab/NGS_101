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

The graph now traces one shared RNA preparation through both adapter ligations and the actual input/treated split, then through both RT, DNA cleanup, production PCR and gel-purified library paths. RNaseH products remain molecular material until column separation. Adapters have an explicit inventory but ligation-ready geometry, exact historical NEBNext PCR/index set and acquired read-cycle plan are not assembled. The public page is a source-qualified partial published instance; no confirmed missing-information claim is accepted.

Ancestry suggestion: Results87 explicitly derives the Ψ-BS-adduct/deletion principle from RBS-seq (refs15–16), while modifying conditions and removing the conventional cytosine-deamination tradeoff. This establishes component strategy, not inherited library adapters or quantitative equivalence. Detailed2024 procedure DOI10.1038/s41596-023-00917-5 is a deferred successor source, not a replacement for the defining2022 publication.

## Bounded complete-preparation model

The same closed defining publication Methods 231–237 supports the added physical
producer path; no successor protocol or historical kit recipe was imported.
Fragmentation and two actual PNK additions precede the 3′ and 5′ adapter
ligations. Activated-adapter AMP/coproducts are not mature DNA backbone, and
unreacted/background classes are not assumed removed before physical cleanup.
The two N5 windows stay distinct source random windows; unsupported molecule-
UMI semantics are not imposed on the actual adapter inventories.

The actual shared 10 µl RNA eluate is split into 1.5 µl input and 8.5 µl
treated material. Only the treated branch receives molecular bisulfite reagents,
column binding and desulfonation. Both routes reach 10.5 µl before actual RT
supplies. RNA Ψ-BS remains distinct from a possible complementary-base omission
in cDNA. RNase H digestion products remain until physical DNA cleanup.

Each source 20 µl cDNA stock contributes 4 µl to one representative actual
15-cycle PCR allocation, leaving 16 µl nominal untransferred stock. Additional
PCR count is unassigned. Actual historical SR/NEBNext primer functional classes
and aggregate enzyme/nucleotide/buffer supplies are represented without
claiming commercial premix packaging or exact primer bodies. Gel purification,
source SE100 loading contributions and stock remainders complete both branches;
no extra QC aliquot or pool between controls is invented.

Normal publication checks pass. Full occupied transcripts, ligation junctions,
mutation fractions, whole primer orders, acquired index mask and read landing
remain unassembled. The intended end/acceptor view is explicitly qualified; the
page retains the complete workflow and source/inference limits.
