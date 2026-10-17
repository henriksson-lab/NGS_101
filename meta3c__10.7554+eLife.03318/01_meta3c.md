# meta3C — research notes

## Citation and sources read

- 🟢 Marbouty M, Cournac A, Flot J-F, Marie-Nelly H, Mozziconacci J, Koszul R. Metagenomic chromosome conformation capture (meta3C) unveils the diversity of chromosome organization in microorganisms. *eLife* 2014;3:e03318. DOI `10.7554/eLife.03318` (published version; no preprint used).
- Original read: PMC4381813 JATS (`eLife.03318_PMC4381813.xml.txt`), with the eLife v1 XML checked for the same passages. Sections inspected in the 2026-10-17 capped pass: Results first meta3C paragraph (line 129) and no-enrichment statement (197); Materials and methods — bacterial library (356), yeast library (360–362), environmental library (364–366), Illumina sequencing (368–371), Table 1 (372–391), size selection/PCR/sequencing (393), processing of PE reads (395–397); river assembly read length (423). Table 2 (lines 427–446) lists PCR primers for Southern-blot probes of the PFGE experiment, not library oligos.
- The oligo scraper was run on the source bundle; its MM70/MM71/MM76/MM77 hits were checked against Table 1.

## Construct path

🟢 Mixed microbial cells are cross-linked with formaldehyde before lysis. The three libraries use different restriction enzymes: the bacterial mix HpaII (`C^CGG`, NEBuffer 1, 3 h 37 °C), the yeast mix DpnII (`G^ATC`, 3 h plus overnight with fresh enzyme) and the enriched river sediment HaeIII (`GG^CC`). Digested chromatin is diluted into 8 ml ligation buffer and ligated with T4 DNA ligase (16 °C, 4 h). Cross-links are reversed overnight at 65 °C with proteinase K in EDTA; DNA is isopropanol-precipitated, RNase-treated, phenol–chloroform extracted twice and ethanol-precipitated. 🟢 Unlike Hi-C, meta3C omits junction enrichment, so ~80 % of read pairs are ordinary paired-end shotgun sequence used for assembly, while contact chimeras group and scaffold contigs.

🟢 Library: 5 µg of each 3C library is sheared on a Covaris S220, purified on QIAquick columns and "processed using a commercial kit (Paired-End DNA sample Prep Kit — Illumina — PE-930-1001)". The DNA is then ligated to custom versions of the Illumina PE adapters (Table 1; 3 h RT, NEB T4 DNA ligase 400 U/µl, 4 µl of 10 µM adapter solution), followed by 65 °C for 20 min. Fragments of 400–800 bp are purified with a PippinPrep; test PCRs set the cycle number for eight large-scale PCRs; products are purified on Qiagen MinElute columns and sequenced on a HiSeq2000.

🔴 Which kit steps were used before the custom-adapter ligation (end repair, A-tailing, cleanups) and the kit PCR primers are not printed in the article. The graph therefore uses a warned boundary for the kit processing step and role tokens for the final library ends. The PCR polymerase is not named.

The model selects the bacterial (HpaII, MM76/MM77) route. The yeast and river branches are named on the page but not drawn.

## Sequencing-library oligos

🟢 Table 1 ("custom-made Illumina adapters"), 5′→3′, no modifications printed:

- MM70 (11-yeast mix; river sediment): `GTANNNNNNAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG`
- MM71 (11-yeast mix; river sediment): `ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNTACT`
- MM76 (3-bacteria mix): `CAGNNNNNNAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG`
- MM77 (3-bacteria mix): `ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNCTGT`

🟢 Methods (line 397): reads are demultiplexed by "the small tag present at the beginning of each sequence (contained in the custom-made adapters)", and PCR duplicates are collapsed using "the six Ns present on each adapter". The model attaches an `inline_barcode` feature to the tag and a `umi` feature to each N6.

🟡 How the pair interlocks (computed, not stated in the paper): the reverse complement of the first nine bases of MM76 is `{{= revcomp("CAGNNNNNN") }}`, which matches MM77 bases 34–42 except at the random positions. MM77 then has one extra 3′ T, consistent with a T-overhang adapter for an A-tailed insert. The reverse complement of the rest of MM76 is `{{= revcomp("AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG") }}`, the 3′ end of the Illumina PE Read 2 primer. MM77 begins with the 3′ part of the Read 1 primer site. Because the two N6 stretches are synthesized independently, they cannot base-pair except by chance. The paper does not say how the adapter solution was annealed, whether the oligos were phosphorylated, or whether the insert was A-tailed. The graph keeps them as separate ordering strands and does not draw an annealed adapter.

🟡 Read layout implied by that polarity: each read would start with N6, then the tag (`CTG` plus the T/A junction base for the bacterial pair), then the insert. The paper only says reads begin with the tag. This order has not been checked against the kit primer sequences, so the page gives read content as text instead of computed cycle spans.

## Published discrepancy shown on the page

- Results line 129: the bacterial library was sequenced "Hiseq2000 − Paired End − 2 × 104 bp".
- Methods line 393: libraries were sequenced on a "HiSeq2000; PE 2 × 100".
- The eLife v1 XML gives the same two values. The model records both and does not assign one acquired cycle count. Line 423 reports "91 bp useful" for the river library, so effective lengths differ further after trimming. It is unresolved whether 104 counts extra cycles or 100 is a generic Methods value.

## Model corrections (2026-10-17 pass)

- Removed the former canonical Illumina endpoint, which drew two sample-index reads. The article describes PE-kit reads with inline adapter tags and gives no index reads. The final ends are now role tokens and no sequencing primers are declared.
- Added the source-stated QIAquick, PE-930-1001 kit boundary, PippinPrep 400–800 bp and MinElute steps, with their discarded column fractions. Named both T4 DNA ligases, the cross-link reversal and recovery chemistry, the adapter heat step and the test/large-scale PCR design.
- Selected the bacterial MM76/MM77 adapter pair as an explicit graph input, so the source guard checks both strands.

## Open questions (deferred, not confirmed absent)

- PE-930-1001 kit original: pre-ligation end processing, PCR primers PE 1.0/2.0 and their polymerase. Not obtained in this pass.
- Dekker 2002 and Oza 2009 3C predecessors cited for the protocol base. Not read.
- Adapter annealing, phosphorylation and exact read cycle boundaries; which read length was acquired for the bacterial library.
- The yeast/river branches (DpnII/HaeIII digestion conditions, Precellys lysis) are not drawn.
