# MethylRAD: standard FspEI preparation

Wang et al. (2015), **MethylRAD: a simple and scalable method for genome-wide DNA methylation profiling using methylation-dependent restriction enzymes**, [Open Biology 5, 150130](https://doi.org/10.1098/rsob.150130).

## Inspected originals and scope

🟢 The published article is available as the PMC4680569 JATS XML and printed PDF. Inspected overview §2.1, Figure 1 and its printed caption (PDF p3), Methods §§5.1–5.3 (XML text lines368–382), and the introduction's explicit 2b-RAD development link (line60). The published supplement `rsob150130supp1.pdf` was read at Table S5, **rendered PDF p6**, alongside its extracted text. These are two chemistry originals; the article's XML/PDF are representations of one work.

The mandatory whole-bundle scraper also reports result-tag sequences from the separate results workbooks. Those are pointers, not adopted oligos or another library-preparation recipe. The remaining original identity is counted conservatively as an incidental scraper-only workbook encounter; neither results workbook was inspected as chemistry authority. All investigations closed before implementation. No confirmed missing-information claim is made.

The selected graph starts with prepared Arabidopsis Col-0 aerial-tissue CTAB-extracted genomic DNA and represents the **standard NNNN adapter branch**. It does not reconstruct CTAB extraction, introduce a glucosyltransferase step, or merge the selective RTR branch. Other Mrr-family enzymes named in the paper, such as LpnPI, are alternatives within this method, not extra protocol counts.

## Whole source orderings

🟢 Table S5 prints all oligos 5′→3′. The independent whole-ordering fixtures are retained in `tools/chem.py`; canonical parts are imported and the actual source adapter pairing is computed locally. `AminoC6` is a physical **3′ amine block**, not a nucleotide. No 5′ phosphate is silently added. `XXXXXX` denotes six sample-index positions; a particular index sequence is unassigned.

| Oligo | Verbatim ordering envelope |
|---|---|
| Standard/RTR Adap-1 sense | `ACACTCTTTCCCTACACGACGCTCTTCCGATCT` |
| Standard/RTR Adap-2 sense | `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` |
| Standard antisense, both adapters | `NNNNAGATCGGAAGAGC(AminoC6)` |
| RTR Adap-1 antisense | `TNNNAGATCGGAAGAGC(AminoC6)` |
| RTR Adap-2 antisense | `CNNNAGATCGGAAGAGC(AminoC6)` |
| P1 | `ACACTCTTTCCCTACACGACGCT` |
| P2 | `GTGACTGGAGTTCAGACGTGTGCT` |
| P3 | `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT` |
| Index primer | `CAAGCAGAAGACGGCATACGAGATXXXXXXGTGACTGGAGTTCAGACGTGT` |

🟡 Antiparallel terminal pairing of the source-written strands is computed with `Scene.anneal`. The terminal NNNN is a degenerate cohesive synthesis pattern, **not a UMI**. Pairing establishes that declared molecular relation; it does not establish the preferred fold, annealing conditions, complete ligation product or end phosphorylation. The RTR TNNN/CNNN source classes are preserved as ordering fixtures but do not enter the selected standard-library graph.

## Source preparation and physical fates

1. Select genomic DNA within the source's 1–200 ng input range. Digest with 4 U FspEI in 15 µl for 4 h at 37°C. The article describes methylation-dependent N12/N16 cleavage and eligible short tags from symmetrically modified sites. Complete genomic cut products, site occupancy, historical activator/buffer allocation and recovery are unassembled. Modified-cytosine enrichment does not independently distinguish 5mC from 5hmC.
2. Remove **5 µl for the 1% agarose digestion QC**; the remaining digest supplies production. Add the source's actual 10 µl ligation master containing the two adapters, ATP and T4 DNA ligase. Its printed 0.2 µM-per-adapter/1 mM ATP/800 U ligase values are retained without deciding whether concentration wording refers to master or final mixture. Adapter, ATP and ligase nodes are components of this same actual addition, not extra doses. Incubate 6–12 h at 4°C.
3. Select **7 µl of ligated material** for a 20 µl P1/P2 PCR with separately named Phusion polymerase, dNTPs and HF buffer. The source specifies 16–22 cycles of 98°C 5 s/60°C 20 s/72°C 10 s and a final 72°C 5 min. Keep untransferred ligation stock separate.
4. Load an unassigned portion of PCR product on 8% PAGE, excise the approximately 100 bp band and diffuse DNA in nuclease-free water for 6–12 h at 4°C. Unloaded PCR stock, gel/workup fractions and retained DNA are distinct. A gel-selected fraction is not assumed perfectly pure.
5. Select **25 ng of recovered DNA** for a separate 20 µl P3/index-primer PCR, with the source's separately named Phusion/dNTP/HF-buffer supplies and 5–7 cycles using the same profile. Purify with QIAquick. The graph retains unused gel-recovered and indexed-library stocks.
6. Figure 1 states physical pooling after sample indexing. Separate sample contributions, the prepared pool, a loading aliquot and the untransferred pooled stock. Methods specify **HiSeq2000 single-end 36 nt**. The six index positions do not establish an acquisition-cycle mask or orientation; those are not inferred from the platform name.

The model does not generate a complete genomic library from matching Illumina handles. Primer envelopes, ordering modifications and adapter-end geometry can be known while whole ligation/PCR products and actual read landing remain unassembled.

## Published overhang discrepancy

🟢 The published overview §2.1 (XML70) says **four-base 5′ overhangs**. The published Figure 1 caption (XML74 and printed PDF p3) instead says **four-base 3′ overhangs**. Both originals' values remain explicit here and in a graph `SourceDiscrepancy` annotation. Table S5 places the NNNN at the antisense **5′ end**, and the computed adapter scene retains that ordering and its terminal paired core. This supports the declared adapter relation; the graph leaves whole enzymatic products and junction chemistry unresolved. It does not silently correct the published caption or claim an independently verified complete digest simulation.

The article explicitly describes a 2b-RAD-inspired method-development route. That documented relationship does not authorize importing the predecessor's SOLiD adapters, enzyme conditions, molecular products or readout.

## Limits of validation

Whole-ordering fixtures cover all nine unique source envelopes, including the unselected RTR alternatives. The external-input source guard covers the four PCR primer envelopes; the adapter molecules produced from the actual master are **outside that guard's external-input scope**, and retain local whole-ordering fixtures instead. A match to the index envelope does not assign its X positions.

The partial graph records source-supported preparation, actual stocks/aliquots and narrow adapter pairing. Its public preparation page retains these open model questions: genomic end distributions, complete occupied substrates/products, historical reagent formulations, ligase end chemistry, exact index assignment/acquisition and full primer landing. The normal source-inventory contract is satisfied for the independently inspected selected preparation; this does not establish complete molecular validation. No source absence, complete reaction validation or full source-review claim is accepted.

## Public preparation scope

The published main overview/Figure1 and Methods5.1–5.3, plus the original rendered TableS5p6, were reread for publication wiring. The existing positive `main`/`s5` source items remain and are mapped to the corresponding digestion/adapter components. A literal independently reviewed operational inventory covers the selected preparation. The5µl digest aliquot is typed as analytical QC; it never enters production again. Whole genomic molecules, exact end chemistry, index assignments and primer/read geometry remain unresolved on the public page. No precursor recipe, current vendor formulation, or confirmed missing-information claim is introduced.
