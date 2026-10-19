# PEM-seq — research notes

Yin et al., *Cell Discovery* (2019), [doi:10.1038/s41421-019-0088-8](https://doi.org/10.1038/s41421-019-0088-8).

## Assay target

The sequenced material is human genomic DNA around a bait double-strand break, quantifying repair outcomes. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (a locus bait primer, a random-barcode bridge duplex and Illumina tagging primers). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 The defining open-access article and Supplementary Methods.
- 🟢 Defining-paper Supplementary Table S11, which prints the bridge strands, P5/P7 tagging primers, I7 series, and locus-specific biotin/nested primers with modifications.
- 🟢 Liu et al., *STAR Protocols* (2022), [doi:10.1016/j.xpro.2021.101088](https://doi.org/10.1016/j.xpro.2021.101088), the authors' stepwise PEM-seq procedure. The running text extracts normally; **Table 2 (the oligo table) does not**. It is a vector object: `pdftotext` returns only its caption, and the separately registered `fx6.gif` is a 100×160-px thumbnail. Read it by rendering page 11 of `PMC9019705_PMC9019705.pdf` at ≥200 dpi, e.g.
  `pdftoppm -r 400 -f 11 -l 11 -png $CHEM_DATA/sources/pem-seq__10.1038+s41421-019-0088-8/PMC9019705_PMC9019705.pdf out`.
  Table 1 (bait DSBs, page 10) has the same problem. The colour coding in Table 2 is load-bearing: the index bases are
  printed in red and the RMB tract orange and underlined, so a text-only reading cannot recover either boundary.

## Construct logic

1. 🟢 Cells are allowed to repair a recurrent bait DSB. Genomic DNA therefore already contains germline/perfect-rejoin, insertion, deletion, inversion, and translocation outcomes before library preparation.
2. 🟢 Extracted DNA is sonicated to 0.3–2 kb. A 5′-biotinylated locus primer then makes one extension product from each recoverable template. The single extension avoids the 80-cycle linear amplification used by LAM-HTGTS.
3. 🟢 Free biotin primer is removed, strands are denatured, and biotinylated extension products are immobilized on streptavidin beads.
4. 🟢 A 5′-phosphorylated bridge upper strand (BA-up, "link ssDNA with RMB") is ligated to the free 3′ end with T4 DNA ligase in 15% PEG8000 at room temperature overnight. Both bridge strands carry 3′ C7 amino blocks. The upper strand contains the pattern `NNNNTNNNANNNTNNNN`: fourteen random bases separated by fixed T/A/T bases. 🟢 BA-lower ("common bridge oligo") is printed as `TGT AGA GCA CGC GTG GNN NNN N-NH2(C7)`; 🟡 computed from the printed strands, its first 16 bases pair the fixed anchor and its six degenerate bases lie past BA-up's 5′ phosphate. 🟢 Liu et al. 2022 STAR Protocols (steps 45–47) anneal the two strands into the 50 µM bridge adapter. 🟡 Those six bases are modelled as the splint over the captured strand's own 3′ end, whose bases the source does not disclose. There is no end repair and no dA-tailing anywhere in this protocol: the ligation substrate is a bead-bound single strand.
5. 🟢 On-bead nested PCR uses a bait-side I5 primer and an I7-series bridge-side primer. The inline sample barcode is the **four** bases between the Read 1 tail and the locus-specific nested primer: `AGCG` (I5-R1), `GCCT` (R2), `AGGA` (R3), `TCAG` (R4). For the RAG1A/B example the nested primer is therefore `AGGATCTCACCCGGAACAGC` (20 nt), and the oligo reads `ACTCTTTCCCTACACGACGCTCTTCCGATCT` + `AGCG` + `AGGATCTCACCCGGAACAGC`. I7-B01 carries the six-base index `CGTGAT` in the ordering oligo; its reverse complement is present on the displayed library strand.
6. 🟢 A tagging PCR with the published P5-I5 and P7 oligos completes the flow-cell ends. The authors sequence 2×150 bp on HiSeq.

## Why the inline I5 barcode is 4 nt, not 5

The concatenated oligo cannot distinguish the parses on its own — `…CCGATCT` + `AGCGA` +
`GGATCTCACC…` and `…CCGATCT` + `AGCG` + `AGGATCTCACC…` are the same string — so this was
resolved two independent ways.

- 🟢 **Within Table S11.** The table prints I5-R1…R4 for three loci. Under the 4-nt parse
  each locus keeps one shared nested primer across all four barcodes: RAG1A/B
  `AGGATCTCACCCGGAACAGC`, MYC1 `AGAGCCATTCTCTGGCTCAG`, RAG1C
  `TATCAATATCCCACTGATGTATC`. Under the 5-nt parse the barcode set would have to be
  `AGCGA/GCCTA/AGGAA/TCAGA` at RAG1A/B and MYC1 but `AGCGT/GCCTT/AGGAT/TCAGT` at RAG1C —
  i.e. a sample barcode whose last base depends on which locus is being amplified. The
  4-nt parse is the only one that makes the barcode a property of the sample.
- 🟢 **STAR Protocols Table 2.** The rendered table (page 11, see above) prints the index
  bases in red for rows I5-Nested-1…4: exactly `AGCG`, `GCCT`, `AGGA`, `TCAG`, each
  following the complete `…CTTCCGATCT` tail. Row I5-Nested-1-Myc6 shows the same red
  `AGCG` followed by `GGAAACCAGAGGGAATCCTC`, the Myc6 nested primer, in black.

The library's bases are unchanged by this correction; only the feature boundary, the
segment length and the Read 1 cycle span move by one base.

## Cross-source conflict: the RMB length

🟡 The two sources disagree about the bridge upper strand, and the model deliberately
follows the defining paper.

| Source | Bridge-adapter-upper random region | Random bases |
|---|---|---|
| Defining paper, Suppl. Table S11 | `NNNNTNNNANNNTNNNN` | 14 |
| STAR Protocols Table 2 | `NNNTNNNNTNNN` | 10 |

The STAR table contradicts its own prose, which says the adapter contains "a 14-bp RMB"
(and its Step 7 narrative repeats it); the defining paper also compares 14-bp with 7-bp
RMB libraries in its supplementary figures. Table S11 is internally consistent with both
papers' text, so `PEM_BA_UP` keeps the 14-base, T/A/T-separated form. The STAR Table 2
printing is recorded here as an error in the later paper, not as a variant to model; a
comment in `tools/chem.py` says so at the constant.

## Represented branch

The page shows the RAG1A/B bait with I5-R1 and I7-B01. Supplementary Table S11 provides alternative bait primers, four I5 barcode/nested variants for each locus, and B/R I7 series. These alternatives change target or sample identity without changing the molecular workflow.

## Read interpretation

Read 1 enters from the nested bait primer and covers the known bait sequence before reaching the repair junction. Read 2 enters from the bridge side and encounters the structured random molecular barcode. The model records the four random tracts as parts of one RMB group, excluding the fixed T/A/T separators from the 14-base UMI length.

PEM-seq quantifies outcomes recovered from the chosen bait side. A repair product that deletes either required primer-binding site is not recoverable, and translocation abundance also reflects three-dimensional contact probability; the library should not be read as a direct frequency measurement for every off-target locus.

## Source audit correction (2026-10-16)

🟢 The defining 2019 article, Methods → Primer extension, specifies sonicated genomic DNA of **0.3–2 kb**. The previous 300–700-bp description was a repository transcription error and has been corrected in the page workflow and graph. The subsequent 2022 STAR Protocols procedure is a separately dated implementation.

## Bounded pass 2026-10-17 (round breadth-first-15)

🟢 Reagents and conditions read in the defining article's Materials and methods
(`s41421-019-0088-8_PMC6434046.xml.txt` lines 131–173) and now named in the graph:
sonication to 0.3–2 kb from 20 µg genomic DNA; a biotinylated bait primer placed within
200 bp of the cut site, annealed over five anneal/denature cycles at a length-dependent
annealing temperature (50–60 °C for 20–50-nt primers); extension with Bst polymerase 3.0
(65 °C 10 min, 80 °C 5 min); 1.2× AxyPrep Mag bead removal of excess bait primer followed
by 95 °C denaturation on ice; capture on Dynabeads MyOne Streptavidin C1 with B&W washes;
on-bead bridge-adapter ligation with T4 DNA ligase in 15% PEG8000 at room temperature
overnight; on-bead nested PCR with Taq for 16 cycles; size-selection bead recovery; and a
tagging PCR with FastPfu that adds the P5 and P7 ends, sequenced 2 × 150 bp on HiSeq.

🟢 The graph now draws the primer-depletion/denaturation step and the post-nested size
selection explicitly, each with its discarded material, and declares the tagging PCR
locally (P7 as an exact source Construct) instead of using the shared generic finisher.

🟡 Table S11 prints the P5-I5 tagging primer in one cell wrapped over two lines, so no contiguous
source line carries the whole oligo; it stays a named input whose transcription is
checked in the selftest. The bridge-ligation substrate geometry and the nested-PCR
amplicon are still not reconstructed as Scenes, and the 2019 vs 2022 random-barcode
conflict remains displayed on the page unchanged.

## Bounded migration check 2026-10-19 (round 2026-10-19-migration-final)

🟢 The module no longer imports the shared `batch_ngs` library recipe or the
`protocol_model` schematic shortcuts `da_tailed` / `closed_circle` (they were imported but
never called; the graph was already native), nor the unused `nebnext`, `nextera`,
`circular` and `reported_adapter_end` helpers. No generic TruSeq library assembly, no
dual-index default and no dA junction is inherited: the P5/P7 ends, the single six-base
i7 index and the four-base inline I5 index are the module's own transcriptions of
Supplementary Table S11, and the canonical Illumina arms come from `lib/illumina.py`.

🟢 The bridge-adapter ligation now carries an explicit
`AdapterLigationGeometry.known("splinted", …)` built from two Scenes: the bead-bound
captured single strand with its 3′-OH and an undisclosed terminal sequence, and the
BA-up/BA-lower duplex drawn from the printed oligos — BA-up 5′-phosphorylated and
3′-NH2(C7)-blocked, BA-lower pairing the 16-base anchor with its six degenerate bases
reaching past that phosphate. 🟡 The `splinted` classification is inferred from that
printed structure (the hexamer annealing over the captured strand's undisclosed 3′ end);
the Methods state only the ligation. The ligated intermediate is a separate state whose
junction is marked. Because BA-up wraps over three printed lines of one cell and the extraction
interleaves the name cells, the adapter input is listed in `SOURCE_OLIGO_EXCEPTIONS` with
the printed lines quoted; both strands remain transcribed verbatim in the selftest.

🟢 The final library's Read 2 side is now assembled from the bridge adapter's own 3′
31 bases (`PEM_BA_UP[-31:]`) plus the three bases that only the I7 primer contributes,
instead of a whole `illumina.INDEX1_PRIMER` segment preceded by a separate "A"; the
library bases are unchanged.

🟡 Still deferred: the nested-PCR amplicon as an exact Construct (both nested primers and
the P5-I5 tagging primer stay named inputs, since Table S11 prints them across broken
cells), and the dead `workflow()` prose chain that the graph-driven page no longer
renders.
