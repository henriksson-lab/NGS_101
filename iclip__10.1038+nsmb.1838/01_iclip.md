# iCLIP: source-qualified first pass

König et al., *Nature Structural & Molecular Biology* (2010), [doi:10.1038/nsmb.1838](https://doi.org/10.1038/nsmb.1838).

This is a capped first source pass, not a finished protocol review. The local graph stops at the replicate-1 cDNA circle and intentionally fails the sequencing publishability gate. The full original protocol continues through BamHI opening and PCR. Their exact reconstruction is deferred, not declared unavailable.

## Sources actually inspected

- Article, published author manuscript `nsmb.1838_PMC3000544.pdf.txt`, Introduction lines 71–83; Results 87–110; Discussion 288–296. These establish UV CLIP, crosslink-induced RT truncation, circular cDNA recovery and random-barcode deduplication.
- Published Springer supplement `nsmb.1838_41594_2010_BFnsmb1838_MOESM21_ESM.pdf.txt`, Supplementary Methods 2048–2195; especially L3 ligation 2080–2084, RT 2106–2126, circularization/opening 2141–2158 and PCR/read mapping 2159–2195. Original PDF pages 15–16 (printed Supplementary Methods pp. 2–3) were visually inspected to verify wrapped primers and bold experiment identifiers.
- The separately cached NIHMS supplement is the author-manuscript counterpart, not an extra protocol. Its corresponding Methods 2024–2155 were compared with the publisher original. The publisher original is the transcription authority below.

The complete bundle oligo scraper was run before reading. Its primer hits were checked in the original passages; it also returns many unrelated validation primers, which are not library oligos. No source-absence conclusion is made.

## What is new

🟢 Crosslink-bound protein remnants often stop reverse transcriptase before the RNA 5′ end. Instead of requiring traversal of both RNA ligation adapters, iCLIP captures the truncated cDNA through self-circularization. The phosphorylated RT primer supplies two inverted PCR handles, a random barcode, an experiment identifier and a BamHI site. The circle is reopened at this internal site and amplified.

## Verbatim oligo inventory

🟢 Publisher Supplementary Methods pp. 2–4 prints:

| Oligo | As written, 5′ to 3′ |
|---|---|
| L3 linker | `5′phosphate-UGAGAUCGGAAGAGCGGTTCAG-3′-Puromycin` |
| RT replicate 1 | `5′-phosphate-NNNCAAGATCGGAAGAGCGTCGTGGATCCTGAACCGCTC-3′` |
| RT replicate 2 | `5′-phosphate-NNNGAAGATCGGAAGAGCGTCGTGGATCCTGAACCGC-3′` |
| RT replicate 3 | `5′-phosphate-NNNTGAGATCGGAAGAGCGTCGTGGATCCTGAACCGCTC-3′` |
| BamHI linearization oligo | `5′-GTTCAGGATCCACGACGCTCTTCAAAA-3′` |
| PCR, P7 side | `5′-CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT-3′` |
| PCR, P5 side | `5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′` |

L3 uses mixed U/T notation in the original. It is preserved, without assuming that the notation identifies every sugar chemistry. The local DNA-copy arithmetic normalizes U to T only for complement calculation. The 5′ phosphate and 3′ Puromycin are retained explicitly in the L3 drawing. PCR oligos reuse canonical historical Illumina constants after comparison with the printed sequences.

🟢 `NNN` is a three-base random barcode; the bold two bases after it identify the experiment. 🟡 The replicate-1 terminal `GAACCGCTC` pairs exactly to an internal L3 interval, leaving two L3 3′ nucleotides beyond that interval. The extension adds the reverse complement of the upstream L3 portion and then the variable RNA. This is computed in `tools/iclip.py`, not eyeballed. The initial manual suspicion of an RT1 mismatch was rejected by this calculation; there is no such discrepancy.

The printed replicate-2 primer is shorter and carries `GA` where mapping text reports the set `TG`, `TC`, `CA`. These observations need a dedicated orientation/variant comparison before any discrepancy interpretation. No repaired primer or whole-family equivalence is used here; the partial model scopes only replicate 1.

## Reaction order

1. 🟢 UV-C irradiation of HeLa cells (100 mJ/cm²), RNase I partial digestion and hnRNP C immunoprecipitation precede the model's starting post-IP material.
2. 🟢 Shrimp alkaline phosphatase removes RNA 3′ phosphate (37 °C, 10 min). RNA ligase attaches the phosphorylated L3 linker overnight at 16 °C. L3's Puromycin blocks its other end.
3. 🟢 Gel electrophoresis and nitrocellulose transfer isolate the complexes. Proteinase K followed by urea/proteinase K extraction recovers peptide-bearing RNA.
4. 🟢 SuperScript III extends the phosphorylated RT primer. Heating and RNase A remove RNA. Denaturing-gel fractions contain 100–175 or 175–350 nt cDNA. The compact graph groups these source operations; it does not claim nucleotide-level peptide chemistry or quantitative recovery.
5. 🟢 CircLigase II joins the cDNA ends at 60 °C for 1 h with MnCl₂. The local circle preserves the complete primer and copied RNA order, with a structured circle-closure relation.
6. 🟢 A complementary oligo exposes the internal BamHI site; FastDigest BamHI acts for 30 min at 37 °C. Illumina paired-end PCR primers and Immomix amplify the linearized template (35 cycles), followed by G-25 desalting. These downstream states are read but not yet modelled.

## Reported readout and open work

🟢 The original reports one Illumina GA2 lane for three hnRNP C replicates, 54-nt sequence reads, experiment/random identifiers and genomic mapping after adapter removal. Calling the PCR primers “paired-end” does not establish that both reads were acquired. The local graph declares no sequencing terminal, index read, run-primer identity or invented read cycles.

Deferred: explicit upstream biological/IP graph, complete paired RT scene, chemical RNA-removal states, gel kept/rejected branches, linearization oligo's nonmatching terminal tail, BamHI cut coordinates on the circular substrate, actual PCR binding/product calculation, replicate-2 reconciliation, and exact acquisition geometry. The source checklist in `catalogue/source_reviews/` records this first pass; a second reader and linked-source review are still required before any missing-information assertion.

## Ancestry candidates

The original explicitly describes extending earlier CLIP by recovering crosslink-truncated cDNA with circularization (article Introduction 71–83, Discussion 288–296). Its cited earlier CLIP originals are uninspected in this cap and have no established parent directory in this note. This supports a conceptual CLIP → iCLIP link with scope limited to crosslink/IP/RNA recovery, not inheritance of exact adapter chemistry.

Existing eCLIP notes name iCLIP predecessors but explicitly say they were not re-read. Candidate `iclip__10.1038+nsmb.1838` → `eclip__10.1038+nmeth.3810` therefore remains a lineage-review suggestion until the eCLIP defining source's specific inherited component is checked. Shared primers alone supply no ancestry claim.
