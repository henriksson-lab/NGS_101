# Nano3P-seq — research notes

Defining source: Begik *et al.*, [10.1038/s41592-022-01714-w](https://doi.org/10.1038/s41592-022-01714-w).

## Evidence used

- 🟢 Nano3P-seq is Nanopore 3′-end-capture cDNA sequencing and accepts total, ribodepleted or poly(A)-selected RNA.
- 🟢 Template switching captures RNA-derived molecules from the native 3′ end without RNA-adapter ligation and without PCR amplification.
- 🟢 Long reads support simultaneous abundance, tail-length and tail-composition analysis for polyadenylated and non-polyadenylated RNAs.
- 🟢 The Methods print the complete blocked `R_RNA`, 5′-phosphorylated
  `D_DNA`, and `CompA_DNA` oligos. TGIRT template-switches from the preannealed
  RNA/DNA starter duplex directly onto the native RNA 3′ end.
- 🟢 After RNA digestion, `CompA_DNA` is annealed to the cDNA handle before
  ligation of Oxford Nanopore Adapter Mix (AMX) from SQK-DCS109. The published
  Supplementary Protocol section 5 (extracted lines 294–322) explicitly says
  this annealing produces a double-stranded cDNA end with an **A overhang**.
- 🟢 The published Supplementary Protocol (MOESM1, oligo table at extracted lines
  99–122) independently prints the three custom oligos, and its preannealing
  instructions follow immediately. Supplementary Software (MOESM4 ZIP, `README.md`
  lines 209–216) supplies custom front/end primer recognition sequences for tail
  analysis; these are analysis parameters, not an established AMX motor-adapter sequence.
- 🔴 The exact SQK-DCS109 AMX motor-adapter sequence and its complementary
  ligating end remain unresolved pending review of the kit-specific ONT source.
  The cDNA-side A overhang is published and does not share that uncertainty.

## Protocol boundary

This is distinct from Oxford Nanopore Direct RNA Sequencing: Nano3P-seq sequences cDNA produced from the RNA 3′ end.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC9834059 JATS241–243, Nano3P-seq library preparation.

Added RNase Cocktail digestion (37°C/10min) and 0.8xAMPure cleanup before CompA annealing. Kept RNA-fragment recovery unresolved; confirmed modified starter/CompA oligos in Methods and PCR-free AMX finish.

Deferred: Supplementary step-by-step protocol and kit-specific motor-adapter end geometry, explicit reagent release and final WSB purification. This pass does not establish a complete source review or publication-wide absence.

## Native graph migration (2026-10-18, batch02)

Inspected: Published PMC9834059 JATS Methods241–243; published supplementary protocol MOESM1 oligo table86–125 and first Nano3P workflow235–355. Mandatory whole-bundle scraper reviewed; source oligos/modifications verified directly.

Native graph preserves modified RNA/DNA starter annealing, physical capture N, actual ss-cDNA and distinct RNA populations. RNase digestion precedes source recovery, followed by exact CompA annealing with its published cDNA-side A overhang and unavailable AMX-side ligation geometry. Source article94°C annealing/100µL ligation selected, with supplement90°C/50µL recipe kept as a distinct published condition. Final WSB recovery and historical SQK-DCS109/R9.4.1 motor boundary retained; comparator dd-cDNA and current ONT recipes are not substituted.

Deferred: Kit-specific AMX sequence/ligating-end/motor geometry and full molecular ligation proof; original historical kit guide not newly inspected. This capped pass does not establish publication-wide absence or complete molecular proof.
