# scMT-seq — transcriptome and DNA methylome

Hu et al., *Genome Biology* (2016), [doi:10.1186/s13059-016-0950-z](https://doi.org/10.1186/s13059-016-0950-z).

🟢 scMT-seq physically separates a nucleus from cytosol before constructing matched RNA and DNA-methylome libraries. It is distinct from scM&T-seq.

## Active source scope

The 2026-10-18 bounded pass inspected three originals: published PMC4858893 Methods205,215,219 and the predecessor identities in references18/19; the study's AdditionalFile1 PDF contents/figuresS1–S5; and the cited 2015 scRRBS protocol's MOESM430 supplement, figures1–5. The supplements inspected contain performance/methylome-coverage material; this is not a publication-wide search for absent sequences. The mandatory bundle scraper was run; its five source-data files exceeding10MB were not scanned. Superseded preprint hits were excluded.

🟢 Cell lysis uses2% Triton,20mM NaCl/Tris and RNase inhibitor with ERCC RNA. The nucleus is physically pipetted in0.2µl into4µl RRBS lysis buffer; cytosol and ERCC RNA remain in the RNA branch.

🟢 The cytosolic branch follows Smart2-seq, substituting SuperScriptIII forII. Amplification and purification precede the selection of0.1ng cDNA for Nextera tagmentation/library construction. The native RT-primer, TSO and PCR strands remain unassigned in this model. First-strand cDNA, switched single-strand cDNA and recovered amplified cDNA are distinct states; RNA templates are independent molecules.

🟢 Released nuclear DNA is digested with MspI with a lambda-DNA spike-in. End repair/dA tailing precedes adapter ligation, bisulfite conversion and purification. Lambda DNA remains an independent population, first introduced at digestion. The ligation-ready drawing shows opposed3′ A ends derived from the supplied duplex, while5′ phosphorylation and adapter geometry remain unassigned. A canonical Y-adapter or TruSeq index length is not substituted.

🟢 Two methylome PCR rounds use20 then12 cycles, followed by180–500bp PAGE selection and HiSeq2500 sequencing. The first PCR receives actual purified converted strands, where converted unprotected cytosines are physical U. T/A pairs arise on copying; final sequencing signal is not substituted for that substrate. The second PCR receives the first-round product. Primer identity relations, nestedness, aliquot choice and any inter-round purification remain unassigned; two cycle counts do not establish those details.

🟡 Both native endpoints and sequencing-primer landing remain roles. Acquired read/index cycles and index identities are unassigned. Nonamplified RNA/DNA or candidate tagmentation outcomes are not automatically physical discard fractions. Explicit cleanup/PAGE fractions state only their described physical separation.

## Predecessor boundaries and retired model choices

The Methods cites Picelli2014, *Nature Protocols*, doi:10.1038/nprot.2014.006, and Guo2015, doi:10.1038/nprot.2015.039. Their exact inherited molecular components still need version-specific source verification. The previous model selected a Smart-seq2 TSO linker through a shared constructor and used an older generic RRBS endpoint with inferred canonical TruSeq arms. Those choices have been removed. Earlier review records remain historical; they do not establish native oligo identity for this study.

Deferred targets: the directly cited2014 Smart2-seq detailed original and2015 scRRBS detailed original, exact study adapter/primer members and modifications, native gap-fill/PCR geometry, inter-round methylome handling, and read/index acquisition. This capped pass does not close linked-source review or accept missing-information claims.
