# Primary template-directed amplification (PTA) — source notes

## Sources read

- 🟢 Gonzalez-Pena et al., *PNAS* (2021),
  [doi:10.1073/pnas.2024176118](https://doi.org/10.1073/pnas.2024176118), is the defining
  PTA paper. The full Methods and SI were read from PMCID PMC8214697.
- The linked bioRxiv version is superseded by the published article and SI. It was not used in the 2026-10-18 native-graph pass.

## Defining amplification chemistry

🟢 A single cell is freeze-thawed, then lysed for 10 min at 40 °C after adding 2 µl of a
1:1 mixture of 5% Triton X-100 and 20 mg/ml Proteinase K. Four microlitres of REPLI-g D2
buffer and 1 µl 500 µM exonuclease-resistant random primer denature the DNA for 15 min at
65 °C; 4 µl Stop solution neutralizes it.

🟢 The 56-µl phi29 amplification mix contained alpha-thio-ddNTPs (TriLink) at equal ratios
"at a concentration of 1,200 µM" (Methods173). 🟡 Whether 1,200 µM is per ddNTP or total, and whether it refers to the mix or the final reaction, is not stated. Amplification runs 8 h at 30 °C and is stopped for
3 min at 65 °C. The alpha-thio-ddN creates an exonuclease-resistant phosphorothioate bond
at the terminated 3′ end. Unlike an ordinary ddN that phi29 proofreading can remove, this
end cannot be rescued, so products stay shorter and undergo less secondary amplification.

The inspected published WGA Methods name an exonuclease-resistant random primer. Its exact bases, length and modification pattern remain unassigned in this bounded scope; no phosphorothioate random hexamer is imported from another MDA protocol.

## Defining sequencing-library route

🟢 The entire PTA reaction is purified with 2:1 AMPure XP beads and enters KAPA HyperPlus
**without fragmentation**. IDT unique-dual-index adapters are ligated at 2.5 µM, followed
by 15 PCR cycles. A 400–700-bp gel fraction is recovered and sequenced on NovaSeq 6000.

The historical KAPA/IDT native adapter and PCR-primer strands, index members and widths, end geometry and run mask remain unassigned in this bounded scope. A named kit and UDI supplier do not establish canonical TruSeq strands or eight-base indexes. The graph carries these native roles without sequence substitution.

The native KAPA preparation boundary is explicitly unresolved. The graph does not convert the alpha-thio-ddN-terminated single-strand class into a repairable duplex by changing its name. Component-specific processing and recovery of blocked ends require the applicable historical kit chemistry. PCR geometry receives the actual native ligation population, and ordinary daughter DNA does not acquire a template phosphorothioate bond through copying. Original genomes, terminated products, other possible PTA classes and later ordinary PCR daughters are independent polymers in the population views.

Source AMPure cleanup and gel excision create physical retained/rejected fractions. The specific allocation of original, terminated, unligated and daughter classes within those fractions remains unresolved; successful amplification or later computational filtering does not establish their physical removal.

## Bounded curation round 2 (2026-10-16)

Inspected: Original PMC8214697 JATS172–179.

Added source400–700bp gel selection after15-cycle library PCR; retained no-fragmentation PTA branch distinct from SCMDA. Removed unverified8-cycle acquired index mask and publication-wide primer-absence wording; canonical geometry remains inferred.

Deferred: Original SI, random-primer/UDI vendor documents and exact acquisition mask; source lysis/neutralization and post-PTA cleanup represented only at summarized scope. This capped pass does not establish a complete source review or publication-wide absence.

## Bounded native-graph migration (2026-10-18, batch 12)

Inspected two published originals: PMC8214697 JATS Results lines 83–85 and WGA/Library Preparation Methods lines 172–177; published SI Appendix title and Fig S1 caption lines 1–24 and 67–73. The mandatory oligo scraper returned zero hits; this does not establish absent information. No preprint chemistry was consulted.

Removed shared PTA source defaults and the canonical TruSeq/UDI endpoint. Native random-primer and kit roles are local declarations. The physical 2:1 AMPure cleanup precedes preparation of the entire PTA reaction without fragmentation, 2.5 µM IDT UDI ligation, 15 PCR cycles and 400–700-bp gel recovery. The initial graph boundary is the denatured genome after source lysis and neutralization; those preparatory operations and analytical quantification remain summarized, with no invented aliquot volumes.

Deferred: applicable historical KAPA HyperPlus guide and IDT UDI ordering/strand specifications; native blocked-end processing, exact random-primer formulation, first-cycle PCR and ligation geometry, sample-specific UDI members, class-specific cleanup recoveries, sequencing run mask, and other SI figures. This bounded migration does not change the blocked missing-information audit status or count as either independent absence review.
