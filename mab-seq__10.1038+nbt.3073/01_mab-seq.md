# MAB-seq: protect unmodified CpG cytosines before bisulfite

Wu et al. (2014), [Single-base resolution analysis of active DNA demethylation using methylase-assisted bisulfite sequencing](https://doi.org/10.1038/nbt.3073).

## Capped sources

One reader,28October2026: article-associated accepted manuscript PMC4269366, Results “Experimental strategy” and Methods “Library generation for genome-scale MAB-seq/BS-seq”, conversion controls, locus-specific and caMAB branches; publisher Supplementary Information, Figure1/2/3 chemistry and controls and Table2 locus-primer inventory. Both source manifests and the mandatory oligo scraper were inspected. PMC article is an author manuscript, not a verified copy of the complete version-of-record Methods. This is a capped first pass, not a completed source or absence audit.

The model selects **genome-scale MAB-seq**. It does not replace its order with the locus-specific protocol, whose M.SssI treatment precedes locus PCR and later NGS library generation. The source's caMAB and redBS experiments remain separate branches in the notes rather than a second registration or a shared default.

## Ordered genome-scale preparation

1. 🟢 Start with1µg unenriched genomic DNA, or200–500ng H3K4me1/H3K27me3 ChIP DNA. Add unmethylated lambda DNA at1:400. The schematic's entry is whole-genome DNA; ChIP extraction is a separately described optional front end, not a universal MAB step.
2. 🟢 Repair DNA and ligate **5mC-containing custom adapters** using NEBNext UltraDNA LibraryPrep. Purify with1.2×AMPure. Printed adapter identities,5′phosphate and terminal phosphorothioate are retained. The enzyme product's exact genomic ends are unassigned; no acquired index member is fabricated.
3. 🟢 Treat adapter-ligated DNA for two rounds with M.SssI and SAM. In each round,1U/µL M.SssI plus640µMSAM incubates4h in25µL; add25µL fresh same-concentration reaction and incubate another8h. Round1 uses Mg-free buffer with EDTA; round2 uses NEB buffer2 containing10mMMgCl2. PCI extraction/ethanol precipitation follow **each** round.
4. 🟢 EpiTect59104 bisulfite conversion follows methylase treatment; the source repeats the kit thermal cycle twice. The matched BS-seq control bypasses M.SssI. These are different analytical treatments.
5. 🟢 Pre-amplify five cycles with KAPA HiFi Uracil+ and NEBNext Multiplex Oligos Index Primers Set1 (indexed and universal primers). Determine final cycle count by qPCR, perform scaled-up PCR, and1.2×AMPure purify. HiSeq2500 single-end60or76bp. Do not turn the five pre-amplification cycles into the final cycle count.

## Source-exact adapters,5′→3′

| Component | Printed bases/modifications |
|---|---|
| Forward | ACACTCTTTCCCTACACGACGCTCTTCCGATC*T |
| Reverse | /5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTC |

The source says the custom adapters are methylated with5mC; the model preserves that source qualification without independently assigning a per-position modification map. `*` marks the final phosphorothioate bond. Ordinary printed `C` letters alone are not evidence that these adapters are unmodified. The adapter pair forms a Y structure with a3′T end; exact pairing is computed from the printed strands. Historical kit end-prep implementation and the exact Index Primers Set1 members remain open source/model questions.

## Molecular products and reported bases

| Starting CpG cytosine | After M.SssI | After bisulfite, as specified by this source | Sequenced base |
|---|---|---|---|
| C | 5mC | protected5mC | C |
| 5mC | 5mC | protected5mC | C |
| 5hmC | 5hmC | source describes bisulfite-resistant behavior; exact derivative not assigned here | C |
| 5fC | 5fC | U | T |
| 5caC | 5caC | 5caU | T |

🟢 Source Results explicitly distinguishes chemical5caU from reportedT. The model preserves that distinction. Non-CpG C is poorly protected by M.SssI, so this table is **CpG-qualified** and must not be reused as a general genome-wide C-to-5mC rule. Conversion is not perfect: incompletely protected C can contribute falseT; measured controls/background correction are analytical necessities, not chemical guarantees. Neither a T call nor a diagram by itself establishes a pure5fC or pure5caC identity: the primary signal combines them.

## Related branches and ancestry

Locus-specific MAB uses four rounds of M.SssI with different additions, then bisulfite, locus-specific KAPA PCR and NGS/Sanger alternatives. Source caMAB treats methylase-protected DNA with NaBH4 before bisulfite to distinguish5caC. These source-local branches are not folded into the genome-scale graph.

The paper identifies standard BS-seq as the analytical baseline; it supplies its own custom adapters and protection order. Similar NEBNext products or comparison with TAB/redBS do not establish whole-protocol ancestry. The later2016 detailed protocol DOI10.1038/nprot.2016.069 is a follow-up source, not a predecessor supplying defaults to the2014 original.

Deferred: complete version-of-record Methods, later detailed protocol, source-linked/historical vendor identities, exact indexed PCR members, full final strand architecture and run primers, and independent second reader. No missing sequence or reagent is confirmed absent.


## Public preparation view, 2026-10-30

The public page declares a **published preparation instance**, qualified to the source
and selected molecular scope described above. Role-level molecular classes are not
particular genomic sequences or complete final-library assemblies. The existing
normal typed graph gate is preserved; computed unresolved chemistry stays visible
in the full workflow.

SAM is now a molecular contribution in each source protection round. Each round remains a clearly labelled compound schedule: initial25µL4h, then fresh25µL same-concentration M.SssI/SAM8h. A single diagram step does not mean simultaneous50µL treatment. Catalyst, nucleotide and buffer classes are exposed from the actual named KAPA Uracil+ ReadyMix supply; this is not three additional fresh stocks. Exact formulation and allocations between preamplification/final PCR remain unassigned. The retained library and analytical qPCR aliquot remain distinct; sequencing now selects a loading aliquot and leaves an explicit stock remainder. No unreported sample pool, historical primer member or whole methylated-adapter modification map is invented.

The source investigation was closed before implementation. This pass rechecked two originals: the associated2014 manuscript genome-scale Methods110 and published supplement original FigS2p2/caption. No accepted/final-version equivalence, later protocol inheritance or complete historical vendor recipe is asserted.

Mandatory bundle scrapers ran in the read-only preflight. The explicit source inventory
was updated with literal curator-selected targets, not graph enumeration. Source and
material-fate checks were added. This is not a full linked-source audit or acceptance
of any missing-information claim. Earlier capped-review statements remain historical.
