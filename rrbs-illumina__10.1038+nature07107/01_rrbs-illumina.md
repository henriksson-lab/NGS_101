# RRBS (Illumina)

Meissner et al. (2008), Genome-scale DNA methylation maps of pluripotent and differentiated cells. DOI: https://doi.org/10.1038/nature07107

The catalogue entry represents the 2008 Illumina implementation of reduced representation bisulfite sequencing, not a claim that RRBS originated in 2008. MspI selects a reproducible CpG-rich genomic fraction before methylation-sensitive chemical conversion.

## Sources read

🟢 The defining article's associated accepted author manuscript, PMC2896277, Methods Summary (text lines 112–120), assay rationale lines 68–70. This is an author manuscript, not a verified publisher version of record.

🟢 Publisher Supplementary Information, `nature07107_41586_2008_BFnature07107_MOESM261_ESM.pdf`, Supplementary Figure S2 and its caption (text lines 1010–1186), Figure S1 caption and Table S2 sequencing context. The illustration supplies C→U chemical conversion and subsequent copy-strand sequencing; its illustrative short genome strings are not library oligos. Other fetched presentation slides were not inspected. The subsequently fetched supplementary-data workbook was sampled only at the HCP-state header (text lines1–5: methylation, expression and histone-state data); its full contents remain deferred. A final whole-bundle scraper rerun still gave only the same illustrative Figure S2 strings.

## Preparation and molecular interpretation

1. 🟢 Digest 1–10 µg mouse genomic DNA with 10–100 U MspI. Isolate 40–120 bp and 120–220 bp fractions; the selected genomic roles remain unbounded in the model, because a size interval is not an exact insert sequence.
2. 🟢 Fill in and add terminal 3′ A. Phenol extraction and ethanol precipitation follow this end processing. DNA and reaction components remain together until that physical workup.
3. 🟢 Ligate pre-annealed Illumina adapters containing 5-methylcytosine in place of cytosine using the specified Illumina DNA preparation kit/protocol. QIAquick cleanup follows. The present source pass does not assign the historical complete adapter polymers or infer them from later TruSeq libraries.
4. 🟢 Treat adapter-bearing DNA using the Qiagen EpiTect Bisulphite kit. The Figure S2 idealization shows unmethylated C becoming U and methylated C retained. Chemical U precedes PCR; the model does not relabel it as T before copying. Conversion efficiency and noncanonical cytosine species are not established by this class diagram.
5. 🟢 Preparative PCR, QIAquick cleanup, then a final 4% NuSieve agarose selection: 130–210 bp or 210–310 bp adapter-ligated products; QIAquick gel recovery.
6. 🟢 Illumina 1G Genome Analyser, 36 bp end reads. Align against a bisulfite-converted reduced representation reference; C/T evidence is interpreted relative to strand and reference. An individual C call alone does not establish complete methylation specificity.

## Evidence boundaries

🟡 A partial graph preserves the source order and separates chemical conversion from copying and readout. Its historical kit, exact adapter ligation geometry, PCR primer members and complete sequencing-primer landing sites remain unresolved. No adapter or index polymer is borrowed from another RRBS implementation. Unknown index/run settings are not a claim of absent indexes.

The PMC supplementary attachment is currently gated, while a publisher Supplementary Information PDF is available and inspected above. Exact equivalence of those two exports was not asserted. The complete publisher article, original Illumina kit revision/protocol, EpiTect workup detail, and cited 2005 RRBS precursor remain follow-up targets. This capped first pass is not a two-reader missing-information review.

## Lineage

The article explicitly cites Meissner et al. (2005), DOI 10.1093/nar/gki901, for RRBS/reduced representation; source scope is the restriction-selected bisulfite strategy. The 2008 Illumina adapter chemistry is owned here and is not inherited from that older method. Later RRBS-based methods may cite this paper; each child needs its own component evidence.
