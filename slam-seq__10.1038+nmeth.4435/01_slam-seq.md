# SLAM-seq

## 1. What it is

SLAM-seq labels newly synthesized RNA with 4-thiouridine and converts the label into a sequencing-detectable base substitution through thiol alkylation.

## 2. Authoritative source

🟢 Herzog et al., *Thiol-linked alkylation of RNA to assess expression dynamics*, Nature Methods (2017), DOI `10.1038/nmeth.4435`.

🟢 Ameres, Herzog and Reichholf, [*Thiol-linked alkylation for the metabolic sequencing of RNA (SLAMseq)*](https://doi.org/10.1038/protex.2017.105), Protocol Exchange (2017), original authored six-page PDF in the [Research Square archive](https://assets-eu.researchsquare.com/files/nprot-6257/v1/8bc3034f-b1da-42a1-a3bb-5446016a9efc.pdf), pp. 2–5. Its separately linked Tables 1–4 are not in that PDF and remain unavailable.

## 3. Construct-changing operations

🟢 Iodoacetamide alkylates 4-thiouridine. Reverse transcription across the modified base generates characteristic T-to-C mismatches; the downstream RNA library method is otherwise conventional.

## 4. Chemistry and library endpoint

🟢 The paper reports standard alkylation as 10 mM iodoacetamide in 50% DMSO and 50 mM sodium phosphate, pH 8, for 15 min at 50 °C, followed by DTT quenching. IAA attaches a carboxyamidomethyl group to s4U.

🟢 Reverse transcription of alkylated s4U produces the diagnostic T-to-C readout. The defining mRNA experiment used the Lexogen QuantSeq 3′ mRNA library kit and Illumina single-read sequencing.

🟢 The authors' detailed Protocol Exchange PDF specifies Lexogen QuantSeq 3′ mRNA-Seq Library Kit for Illumina, catalog 015 (p. 2), identifies the **FWD** kit (p. 4), and recommends HiSeq 2500 single-read 100 (SR100) mode (pp. 4–5). It describes 100 µM s4U for the tested mESC setting, DTT during RNA isolation, and IAA treatment for 15 min at 50 °C followed by DTT quenching. Table 4, the full IAA reaction setup, is linked separately and has not been checked; the paper's own 10 mM IAA recipe above remains the source for that concentration.

🟡 The paper names the commercial kit but does not print its adapter oligos. The schematic therefore marks the Illumina-compatible single-index adapter skeleton as inferred rather than presenting kit bases as published facts.

🟢 Supplementary Table 3 does print the oligos for the **separate synthetic let-7 primer-extension control**: RT primer `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`, Solexa_PCR_fwd `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`, and Solexa_IDX_rev `CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` (`N6` is the index in the table). The same table identifies adapter-flanked synthetic let-7 and its s4U-containing counterpart for that control. These oligos substantiate the control-library construction; they do **not** resolve the QuantSeq kit oligos used for the defining transcriptome experiment.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC5712218 JATS100–138, carboxyamidomethylation, SLAM-seq in mESCs and RNA library preparation.

Added DTT quench to alkylation transition. Corrected acquisition explanation to study SR50 with SR100 for Fig3 transcriptional-output measurements, distinct from later protocol recommendations; no universal run silently assigned.

Deferred: Published linked ProtocolExchange version, original QuantSeq kit primers and end geometry, RNA precipitation/total RNA mixture states and selected library run. This pass does not establish a complete source review or publication-wide absence.
