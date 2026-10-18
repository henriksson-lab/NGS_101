# Illumina DNA Prep shotgun metagenomics

## 1. What it is

Illumina DNA Prep applies bead-linked transposomes to total DNA extracted from a microbial community. On-bead tagmentation fragments many genomes without locus selection, and a limited-cycle PCR adds dual indexes and complete flow-cell arms. The resulting shotgun library samples both taxonomic and functional sequence content.

## 2. Authoritative sources

🟢 Illumina, *DNA Prep Reference Guide*, document 1000000025416 v12 (April 2026), defines bead-linked tagmentation, post-tagmentation wash and index PCR.
🟢 The guide's downloadable PDF was checked directly. Its consumables section directs readers to Illumina's separate *Adapter Sequences* reference for index-adapter bases; the guide itself names compatible index sets but does not print an ordered index-primer pair. The index values in a particular library depend on the set and well chosen.
🟢 The official *Illumina Prep and Nextera Kits* adapter table prints the Read 1 transposase adapter as `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`, the Read 2 adapter as `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG`, and index-PCR primer arms as `CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG` and `AATGATACGGCGACCACCGAGATCTACAC[i5]TCGTCGGCAGCGTC`. The official *IDT for Illumina Unique Dual Indexes* table gives the well-specific i7/i5 bases for the Set A kit used in the application note. These public sequences resolve the adapter arms; a specific library still requires its chosen well.
🟢 The Set A i7 and i5 index regions are each ten bases long. The schematic uses that source-known width while keeping the unspecified well's bases symbolic.

🟢 Illumina, *Shotgun metagenomics with Illumina DNA Prep on the NovaSeq 6000 System*, M-NA-00010 v2.0, identifies DNA extraction followed by Illumina DNA Prep as the demonstrated shotgun workflow. Library-preparation paragraph (text lines 215–221): eleven libraries were made with Illumina DNA Prep, (M) Tagmentation (96 Samples, IPB), catalog 20060059, from 10 ng total DNA normalized to 30 µl with nuclease-free water before on-bead tagmentation, and amplified with IDT for Illumina DNA/RNA UD Indexes Set A, catalog 20091654. Table 2 lists these libraries as 10 ng input, 350 bp insert. Sequencing (lines 239–242): pooled libraries ran 2 × 251 bp on HiSeq 2500 (Rapid Run) and NovaSeq 6000 (SP flow cell). Extraction used the QIAGEN DNeasy PowerSoil or Revvity chemagic DNA Stool kits.

🟢 Guide reaction steps used here (pp. 23–29): TB1 + BLT, TAG 55 °C 15 min; TSB, PTC 37 °C 15 min; supernatant discarded, two TWB washes plus a third TWB kept until PCR; EPM + pre-paired i7/i5 index adapters, BLT PCR 68 °C 3 min, 98 °C 3 min, X × (98 °C 45 s, 62 °C 30 s, 68 °C 2 min), 68 °C 1 min; double-sided IPB cleanup and RSB elution. The guide (p. 26) requires configuration for 10-base UD codes; p. 32 positively recommends 2 × 151 cycles plus 10 cycles per index read.

🟡 The ©2024 application note and April 2026 guide are separate source versions. The native graph selects the v12 reference implementation at 10 ng: eight PCR cycles from its input table, and its recommended 151 × 10 × 10 × 151 run. These guide choices are not proof of the recipe or index acquisition executed in the 2024 application-note experiment, which reports 2 × 251 biological reads. No full revision equivalence has been established.

## 3. Construct logic

1. Extract total DNA without selecting a marker locus.
2. Bead-linked Tn5 complexes fragment and tag the DNA (TAG 55 °C 15 min). TSB stops the reaction; the supernatant and TWB washes physically remove unbound DNA and reaction components.
3. Index PCR with EPM and a Set A pre-paired well completes the P5/P7 arms and introduces 10-base i5/i7 sample indexes (8 cycles derived for 10 ng).
4. A double-sided IPB cleanup discards the larger first-bead fraction and the short second supernatant.
5. The selected 2026 guide configuration uses 2 × 151 biological cycles plus 10-cycle indexes. The 2024 application note separately reports 2 × 251 biological reads; historical index acquisition remains unassigned.

At 10 ng the library is below the guide's 100 ng normalization range, so per the guide its yield is not normalized and it is quantified individually.

## 4. What distinguishes it from amplicon profiling

There is no locus-specific primer and therefore no 16S or ITS boundary in the molecule. “Shotgun metagenomics” describes the biological sampling strategy; the library molecule itself is an Illumina DNA Prep/Nextera-style tagmentation product.

## 5. Proprietary regions

🟡 Current reagent mixes are commercial, but the final Nextera-compatible primer architecture is public. The schematic models the public adapter/read structure and treats the bead linkage as reaction geometry rather than DNA sequence.

## 6. Open questions

- The current reference recipe does not establish historical 2024 PCR or index acquisition; that comparison remains deferred.
- The specific Set A wells used are not reported, so index bases stay symbolic.
- The umbrella *Illumina Adapter Sequences* landing page (`adapter-sequences.htm`) is a manual row in the source manifest (no article text returned); its child pages were read.

## Bounded native-graph migration, 2026-10-18

Three originals were inspected: v12 guide identity and pp. 23–32 (text lines 1576–2230), the ©2024 application note library preparation/sequencing (211–242 and identity/footer), and the official *Illumina Prep and Nextera* adapter table (182–205). The mandatory bundle scraper ran; relevant bases and method passages were checked directly. The Set A well table was not reread this pass.

The local final assembly and index primers use the source-table arms, 10-base guide-defined UD regions and a variable-length genomic role. The preceding native BLT states retain the S5/ME and S7/ME supplied arms on separate strand classes; native gaps, pairing extent and bead occupancy remain unassigned. Index PCR consumes that preceding state rather than a completed expected duplex. Non-amplification is not called physical discard; wash/IPB removal remains tied to actual source physical separations.

This pass does not establish complete historical BLT composition, precise pre-PCR strand geometry, an actual Set A well, or equivalence of the 2024 experiment and the 2026 guide. It accepts no missing-information claim.
