# PIP-seq (protein interaction profile sequencing)

**Defining source:** Silverman *et al.*, *Genome Biology* (2014),
doi:10.1186/gb-2014-15-1-r3. The open article and four-branch workflow were read.

## Evidence used

- 🟢 Cross-linked lysate is divided between protein-retaining footprint samples and
  proteinase-K/SDS-treated structure controls.
- 🟢 Each protein state is divided again between single- and double-strand-specific RNases,
  producing four libraries.
- 🟢 Strand-specific sequencing compares footprint and protein-depleted controls to call
  protein-protected sites while the two RNases retain RNA-structure information.
- 🟢 The defining paper's “PIP-seq library preparation” section says strand-specific library construction followed earlier nuclease-structure studies by Li *et al.* (2012) and Zheng *et al.* (2010), then added Illumina DSN normalization. It does **not** identify a commercial library-preparation kit. The four library branches used standard 50-bp single-end HiSeq 2000 reads.
- 🟢 The cited Zheng *et al.* method's supplementary Text S1 specifies Illumina Small RNA Sample Prep v1.5 kit after fragment end repair (PMC2947979, `pgen.1001141.s017.doc`, “dsRNA-seq library preparation”). Thus a commercial kit is part of the *cited predecessor chemistry*, even though PIP-seq itself describes the library step by reference rather than giving a kit name.
- 🔴 The defining article and checked predecessor text do not print the complete library adapter, RT, PCR or sequencing-primer oligos. Li *et al.* (2012), the other cited predecessor, remains to be checked because its full text was not retrieved by the source fetcher.

## Naming

This is unrelated to the particle-templated single-cell RNA-seq method also called
PIP-seq. The catalogue and URL qualify this page as “protein interaction profile.”
