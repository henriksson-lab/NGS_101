# Start-seq

[Nechaev et al., Science (2010)](https://doi.org/10.1126/science.1181421) sequences short capped nuclear RNA. Its 5′-end libraries map initiation; separate custom libraries read RNA 3′ termini associated with promoter-proximal pausing. These are distinct library preparations, not the two reads of one asserted paired-end construct.

## Original deposited library routes

The defining [GSE18643 series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18643) links nine samples. The displayed examples are the untreated [5′-read sample GSM463298](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM463298) and LacZ-depleted [3′-read sample GSM463301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM463301). Their independently read original extraction fields specify:

1. Select **25–100-nt nuclear RNA** on a gel.
2. Treat sequentially with 5′ polyphosphatase and Terminator 5′ exonuclease.
3. Ligate the branch-specific **3′ linker**, then gel-purify.
4. Treat sequentially with alkaline phosphatase and tobacco acid pyrophosphatase.
5. Ligate the branch-specific **5′ linker**, then reverse-transcribe and amplify with the named Illumina small-RNA kit reagents.

The 5′-read sample uses standard Illumina small-RNA linkers. The 3′-read sample specifies two custom linkers and a custom RT primer. No standard or later TruSeq oligos are substituted for that branch. The new protocol-owned graph expands this order explicitly; shared constructors no longer choose its adapters, RT/PCR formulation or endpoint.

The two displayed deposits have different sample conditions. The graph starts from two independent nuclear-RNA inputs and does not claim they were physically split from one recovered sample. A common enzyme-type input expresses the shared named polyphosphatase reagent, not shared biological identity or a known reagent lot.

## Molecule fates and read processing

Terminator treatment creates a logical class of degraded RNA products. Their eventual recovery/removal is unresolved; enzymatic degradation is not an immediate physical discard of intact RNA. No complete depletion efficiency is asserted. Excising selected gel material is a physical separation, with the precise composition of excluded post-ligation material unassigned.

The original samples report **Illumina Genome Analyzer II**. Their reads were **trimmed to 26 nt** before alignment: this is a computational rule, not proof of acquisition cycles. Mappable reads from biological replicates were merged during analysis. Exact acquired masks, index reads and molecular primer sites remain unassigned.

The article's observed promoter-proximal pause-position interval is not an RNA gel-selection interval. The earlier repository used the former as a selection rule; the original GEO records resolve the RNA interval as 25–100 nt.

## Historical kit and successor boundaries

The recovered [Illumina Small RNA Sample Preparation Guide 1004239 Rev. D](https://sapac.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_legacy/SmallRNA_SamplePrep_Guide_1004239_D.pdf) is dated **August 2010**, after the experiment. It identifies SRA adapters, RT primer and GX1/GX2 PCR components by part number and describes a standard 5′-then-3′ ligation workflow. The author-deposited Start-seq order is explicitly modified. This positive manufacturer original does not identify the exact experimental kit edition, custom oligos or selected enzyme formulation.

[Henriques et al., Molecular Cell (2013)](https://doi.org/10.1016/j.molcel.2013.10.001) is a later study: its scRNA-seq section delegates the 3′ library method to Nechaev 2010 and adds synthetic capped spike-ins for selected experiments. Its short-RNA qPCR uses a distinct poly(A)/RT procedure; that is not assigned to the 2010 sequencing library. Its extended-methods original remains access-blocked.

## Original supplement and remaining work

The exact defining original is `NIHMS401007-supplement-supp_material.pdf`, identified by the PMC JATS. Current PMC download routes return a gate; the TDM Cloud record provides no attachment. The original published article also identifies the Science DC1 landing; its registered retrieval returned HTTP 403. The identified author ResearchGate SOM record currently offers request-only access, not the original file. No speculative PDF filename is counted as another original document.

The missing source remains an access blocker, not confirmed scientific absence. The remaining seven GEO samples, linked processing README, experiment metadata and other chemistry-relevant source children remain explicit obligations in the review. Ordered 2010 library reagents and run assignments stay unassembled model data; no missing-information claim is accepted.
