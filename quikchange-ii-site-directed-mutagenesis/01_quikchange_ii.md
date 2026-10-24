# QuikChange II

## Source and scope

🟢 Agilent instruction manual for catalog200523/200524, **Revision E1**, own cover and printed pp1–5 read from the original official PDF through the browser. The ordinary validated local fetch returned HTTP403; that is recorded in the external manifest. The browser-readable original identifies Revision E1 (2023 footer), not the older Revision C mirrored elsewhere. No revision equivalence is assumed. This schematic represents the vendor's control-primer example with an unassigned plasmid body, not a target-specific design.

The original source is <https://www.agilent.com/cs/library/usermanuals/public/200523.pdf>. A local PDF/text twin remains desirable for repeatable transcription checks. Running the primer scraper on the external bundle finds no eligible original while that local acquisition is blocked; that is not evidence of missing oligo disclosure. The original's printed materials table discloses the two control oligos, which were read directly.

## Chemistry boundary

🟢 The primers contain the intended change, anneal at the same location on opposite template strands and are complementary to one another. PfuUltra performs primer-directed strand synthesis without primer displacement. The stated product is a duplex with staggered nicks. DpnI digests methylated/hemimethylated parental DNA. The nicked product enters competent cells, where nick repair occurs. There is no in-vitro kinase/ligase closure step in this selected source recipe.

The graph keeps this separate from back-to-back Q5/KLD inverse PCR. It does not interpret the two overlapping QuikChange oligos as an inward-facing conventional PCR pair on a fabricated linear template, or treat the pre-transformation product as a covalently closed strand.

🟢 Control primer1 is `CCATGATTACGCCAAGCGCGCAATTAACCCTCAC`; primer2 is `GTGAGGGTTAATTGCGCGCTTGGCGTAATCATGG`. Their complete complementarity is enforced in construction. A full control plasmid, exact nick positions and strand-by-strand parental methylation/carryover remain unassigned.

## Limits

The intended repaired endpoint is a partial construct class. Source primer transcription and complementarity do not prove nick positions, clone genotype, quantitative background removal or reaction yield. Digested parent fragments remain a mixture class with unresolved subsequent fate; a physical purification step is not invented. Local-source checks skip when the official original is unavailable. No absence claim or complete independent source audit is accepted.
