# ClickSeq

Routh et al.(2015), *ClickSeq: Fragmentation-Free Next-Generation Sequencing via Click Ligation of Adaptors to Stochastically Terminated3′-Azido cDNAs*, [defining article](https://doi.org/10.1016/j.jmb.2015.06.011).

## Source scope

🟢 Associated accepted manuscript PMC4523409/NIHMS705616, main30–34 for the core preparation,38–58 for experiments/limitations and Fig1–2 captions. The source explicitly differs from conventional RNA fragmentation plus enzymatic ligation. Final published-method applicability remains unverified. This is the RNA library method, not the unrelated DNA-lesion Click-code-seq.

The mandatory actual-bundle scraper returned no oligos in the acquired main. Main JATS explicitly identifies `NIHMS705616-supplement.docx` containing TableS1; the generic fetch did not acquire it. A declared standard attempt at that exact source link reached the PMC preparing-download gate. No absence claim follows. Exact full oligos, modifications beyond those explicitly named in main, reaction formulation and acquisition parameters remain unreviewed. No later protocol is substituted. One original read; the named chemistry supplement is blocked. No complete independent original/linked-material audit.

## Chemical and physical sequence

🟢 The selected source100ngFHV RNA preparation is reverse-transcribed from a semirandom primer containing an Illumina3′-adapter region and six random priming nucleotides (Oligo1, source TableS1). The random priming site is not a UMI. dNTPs are supplemented with3′-azido-2′,3′-dideoxynucleotides, which can terminate newly synthesized DNA. The FHV and independent CrPV libraries use a1:20AzNTP:dNTP ratio (main40), distinct from the broader targeted-Primer5 optimization series (main38). This produces a range of copied cDNA lengths without cutting the original RNA. Site choices, ratio-dependent distributions and complete-runoff background are not predicted by the model; the original RNA is not fragmented.

Source DNA cleaning/desalting removes unincorporated AzNTPs. RNA, free-primer and other DNA-class recoveries are not inferred to be perfect. Azido-cDNA is then mixed with excess5′-hexynyl-modified Illumina5′-adapter (Oligo2), active Cu(I)-TBTA and ascorbate; the stock formulation remains unreviewed. Successful joining makes a triazole-linked DNA backbone. Unclicked cDNA, possible nonterminated DNA, free adapter and other reaction material remain classes in the mixture; no universal conversion fraction is assumed. The validation ligation percentage and targeted viral-RNA oligos are not silently imposed on every RNAseq library.

The local selected-template constructor preserves the supplied terminal nucleotide and both DNA bodies while consuming the selected azide/alkyne termini into a nonphosphodiester junction. It adds no DNA nucleotide for the triazole bridge. Explicit finite demonstration sequences default to inferred; default source states retain unassembled sequences. The selected relation predicts neither successful joining of all molecules nor a preferred RNA termination position.

🟢 Main reports direct PCR from the click reaction using Primers3/4. It names common DNA polymerases capable of inefficient read-through but does not establish full library-specific enzyme/thermal details in the acquired source. No invented pre-PCR cleanup is inserted. PCR produces new DNA alongside original modified templates and remnants; the template's chemical triazole is not copied onto every newly synthesized strand. Exact copied sequence, errors/indels near the junction, length and yield remain unassembled. Source AMPure size selection and gel extraction follow PCR; main40 specifies final300–500nt library selection for the FHV/CrPV context.

## Source branches and readout limits

Main40 reports separate libraries from100ngFHV and100ngCrPV RNA,1:20AzNTP:dNTP,300–500nt selection and IlluminaHiSeq single-end100nt acquisition. The graph selects the FHV preparation; host mRNA mapping in main42 comes from that FHV dataset, not an independent cellular-mRNA preparation. Targeted-Primer5/gel validation remains a separate experiment. Upstream extraction is outside this library graph. Index acquisition, primer landing and full end geometry remain unassembled; TableS1 and precise chemistry methods are required. No contemporary kit recipe or index value is borrowed.

The source reports low observed recombination and inefficient read-through in its experiments. These observations are not source-independent accuracy guarantees or constructor invariants. The model remains a private incomplete operational/chemical-class graph, with no finished schematic or accepted missing-information claim.
