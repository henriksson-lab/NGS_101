# Primary template-directed amplification (PTA) — source notes

## Sources read

- 🟢 Gonzalez-Pena et al., *PNAS* (2021),
  [doi:10.1073/pnas.2024176118](https://doi.org/10.1073/pnas.2024176118), is the defining
  PTA paper. The full Methods and SI were read from PMCID PMC8214697.
- 🟢 The linked bioRxiv version,
  [doi:10.1101/2020.11.20.391961](https://doi.org/10.1101/2020.11.20.391961), contains the
  same protocol and a higher-resolution reaction schematic.

## Defining amplification chemistry

🟢 A single cell is freeze-thawed, then lysed for 10 min at 40 °C after adding 2 µl of a
1:1 mixture of 5% Triton X-100 and 20 mg/ml Proteinase K. Four microlitres of REPLI-g D2
buffer and 1 µl 500 µM exonuclease-resistant random primer denature the DNA for 15 min at
65 °C; 4 µl Stop solution neutralizes it.

🟢 The 56-µl phi29 amplification mix contains an equal mixture of alpha-thio-ddNTPs at
1,200 µM total final concentration. Amplification runs 8 h at 30 °C and is stopped for
3 min at 65 °C. The alpha-thio-ddN creates an exonuclease-resistant phosphorothioate bond
at the terminated 3′ end. Unlike an ordinary ddN that phi29 proofreading can remove, this
end cannot be rescued, so products stay shorter and undergo less secondary amplification.

🔴 The paper names an exonuclease-resistant random primer but does not print its sequence,
length or modifications. The schematic therefore keeps that end as an undisclosed role;
it does not import the commonly used phosphorothioate random hexamer from another MDA
protocol.

## Defining sequencing-library route

🟢 The entire PTA reaction is purified with 2:1 AMPure XP beads and enters KAPA HyperPlus
**without fragmentation**. IDT unique-dual-index adapters are ligated at 2.5 µM, followed
by 15 PCR cycles. A 400–700-bp gel fraction is recovered and sequenced on NovaSeq 6000.

🔴 The article identifies the kit and adapter supplier but does not print the IDT UDI
adapter or PCR-primer sequences, index values or read lengths.

🟡 KAPA HyperPlus and IDT UDI adapters create a standard dual-index TruSeq-shaped library.
The page uses canonical TruSeq arm bases, visibly marked inferred, to show how the four
standard Illumina primers bind. This does not claim a particular IDT index set or index
value. No genomic read-cycle lengths are asserted.

## Bounded curation round 2 (2026-10-16)

Inspected: Original PMC8214697 JATS172–179.

Added source400–700bp gel selection after15-cycle library PCR; retained no-fragmentation PTA branch distinct from SCMDA. Removed unverified8-cycle acquired index mask and publication-wide primer-absence wording; canonical geometry remains inferred.

Deferred: Original SI, random-primer/UDI vendor documents and exact acquisition mask; source lysis/neutralization and post-PTA cleanup represented only at summarized scope. This capped pass does not establish a complete source review or publication-wide absence.
