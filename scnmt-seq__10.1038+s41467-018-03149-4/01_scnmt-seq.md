# scNMT-seq — RNA, methylation and accessibility

Clark et al., *Nature Communications* (2018), [doi:10.1038/s41467-018-03149-4](https://doi.org/10.1038/s41467-018-03149-4).

## Active source scope

The capped 2026-10-18 source pass inspected three originals: published PMC5823944, Methods136–144 and interpretation Results; original MOESM1 SupplementaryInformation figures, including context filtering; and the MOESM2 EditorialNote/reviewer-report head. The editorial material supplies no recipe facts. The mandatory source-bundle scraper returned24hits in8hit-files; one15MB source-data workbook was skipped. This is not a source-absence finding. Implementation used the supplied capped-source facts and opened no additional original.

## Native reaction branches

🟢 Cells enter a2.5µl M.CviPI/SAM reaction for15min at37°C. Accessible GpC can acquire5mC; endogenous CpG methylation is independently retained.5µl RLTplus stops the reaction. Oligo-dT magnetic beads capture RNA; the DNA-containing lysate and three15µl bead washes are transferred into the DNA branch. RNA and DNA remain independent molecules, sharing plate-well provenance rather than an invented inline cell barcode.

🟢 The RNA branch uses100U SuperScriptII with42°C60min,50°C30min and60°C10min incubation. A named TSO contributes a template-switched single-strand cDNA precursor. ISPCR/KAPA amplification uses18cycles;1:1AMPure purification recovers cDNA in20µl water.100–400pg cDNA enters a one-fifth-volume NexteraXT preparation. Its actual tagged cDNA precursor is distinct from the completed indexed endpoint; native adapter/gap/PCR geometry remains unassigned.

🟢 DNA is purified with0.8×AMPure into10µl water, then bisulfite converted with the half-volume EZ Methylation Direct MagBead kit. Original converted strands carry physical U at converted unprotected-C positions; final T readout is not substituted for those molecules.

🟢 The Methods names **6NF for both synthesis stages**. The first uses50U Klenow at37°C for30min and repeats four more times, for five total rounds.20U ExonucleaseI treatment at37°C precedes0.75×AMPure purification. Chemical hydrolysis and physical cleanup are distinct; the source scope does not assign the exact surviving original/copy strand spectrum.

🟢 The second stage again names6NF, with98°C2min/ice setup and Klenow at37°C for90min. Another0.75×AMPure purification recovers the actual native copied precursor. PE1.0/iTAG/KAPA PCR uses14cycles, followed by0.7×AMPure purification. No N9 tail, historical scBS-seq final scaffold or default8-base iTAG index is supplied.

🟡 Complete native6NF strands, handles and tail widths remain unassigned. Likewise, native PE1.0/iTAG species and sequencing-primer geometry remain unassigned. Identifier metadata records sample/plate-well meaning without supplying a barcode count, sequence or cycle span. No UMI or inline cell barcode is invented.

## TSO boundary

Earlier model notes compared the AC linker in Picelli2013 SupplementaryTable2, rGrG+G row, with ACAT in Picelli2014 ReagentSetup. The study cites both sources, which does not choose its ordered TSO. The current model retains an unassigned native TSO role and inherits neither sequence. This removes the previous guessed ACAT declaration; historical source-review records are retained without upgrading their scope.

## Cohorts and interpretation

🟢 The graph selects EL16 RNA and DNA libraries sequenced on HiSeq2500 with paired125-base reads. The initial20-library DNA QC cohort, including3controls, used MiSeq single-end50. The separate E14 embryoid-body cohort used NextSeq paired75 for DNA and single-end75 for RNA. Its RNA run has no acquired Read2. Index acquisition remains unassigned for all these recipes; separate runs are not merged into one cycle mask.

🟢 GCG contexts are ambiguous between GpC accessibility marking and endogenous CpG methylation; CCG can reflect off-target marking. Their exclusion is computational interpretation, not physical DNA removal. The GpC and CpG readouts therefore retain their context restrictions after sequencing.

Deferred targets: version-specific Smart-seq2 predecessor components and ordered TSO, native6NF and PE1.0/iTAG strand identities/handle lengths, native NexteraXT index/gap/PCR geometry, index acquisition and source-data workbook relevance. No linked-source closure or missing-information claim is accepted by this bounded migration.
