# TaDRIM-seq

[Published defining article](https://doi.org/10.1038/s41467-024-53534-5).

🟢 Two originals: published PMC11811046 main69–71/169–193 and published SI1 SupplementaryNote1 detailed protocol, source oligos/stock tables at originalPDF13/14 and downstream selected recipe through final200ng library input. Original modified oligos and tables were rendered/read. Mandatory scraper45candidates/four files; two large source-data files skipped its normal cap, not absence. FISH targets, comparator RT&Tag/ChRD-PET/CUT&Tag and result workbooks remain separate scopes, not imported recipes. No complete source/absence audit.

TaDRIM-seq combines antibody-targeted proteinG-Tn5 DNA tagging with RNA3′ bridge ligation, RT and in-nucleus proximal DNA ligation. The same preparation produces RNA-DNA and RNA-RNA proximity reporters; these are not proof of direct RNA-RNA binding or separate target-panel protocols. Selected model scope is mammalian cells; no plant extraction values are imported.

## Whole source stock orderings

pG-Tn5 A: `5′phos-CTGTCTCTTATACACATCTT-3′`; B: `5′phos-AAGATGTGTATAAGAGACAG-3′`. Both actual terminal bases and5′ phosphates are retained, without replacing them by generic Nextera ME sequences. SI's heading calls its reagent “circularization mixture”, but the actual operation mixes20µl pG-Tn5,3.5µl annealed A/B and14µl coupling buffer30°C1h to form TTE; no DNA circularization is inferred from that heading. A5µl selected TTE aliquot has an untransferred37.5µl-stock remainder.

Bridge forward precursor: `5′[5Phos]CCTCGATAC/iBIOdT/TGAGCTGAC-3′`; reverse: `5′[5Phos]GTCAGC/iBIOdT/CAAGTATCGAG-3′`. Each internal biotin-dT is a real T nucleotide, not an omitted token or extra free biotin molecule. SI adenylates the forward precursor using12µl10×buffer,12µl1mMATP,6µl100µM forward,12µlMthRNA ligase and78µlwater65°C1h/85°C5min. Add6µl reverse, source stepped anneal, column/elute50µl. The final forward5′App stated in main185 follows this preparation; it is not a contradictory phosphate claim. The bridge used for this method ends `GAC`, whereas the separate ChRD-PET comparator forward includes an extra terminalT: versions/components are not merged.

🟡 The local RNA3′OH acceptor class is an inferred intended substrate after PNK, with constructor inference enabled by default. Actual RNA end distribution, bridge priming occupancy and complete sealed products are unassembled. Source-known bridge bases and biotin-dT positions remain explicit.

## Selected sequence of operations

Crosslink cells1%formaldehyde10min/0.2Mglycine5min, wash. Isolate~one million cell nuclei with lysis buffer, bind10µg target antibody4°C2h, wash, add5µl preloaded TTE4°C1h and wash. Retained transposase/adapter components, rather than re-added free stocks, perform Mg-initiated tagging37°C1h. EDTA stop, wash, SDS62°C10min and Triton quench follow. RNaseIf37°C3min/ice5min precedes PNK/ATP37°C40min.

Use prepared biotinylated bridge in400µl truncatedRnl2 RNA-ligation mixture16°Covernight. Washfive times, then MMLVRNaseH-minus RT with0.5mMdNTP42°C1h; full priming/ligation products remain unassembled, and any primer component must come from the existing bridge rather than an invented new RT oligo. Washthree times, DNA ligase/ATP16°Covernight for proximal joining.

ProteinaseK/SDS65°C3h reverses crosslink/protein context; SI specifies organic aqueous-phase recovery, explicit2µlGlycoBlue molecular-carrier addition followed by acetate/isopropanol precipitation and ethanol washes,50µl recovery. Gap-fill with T4DNApolymerase/dNTP40min; main says1mMdNTP while SI's table gives8µl10mM stock in100µl mixture. Preserve the source table and label this concentration convention/difference unresolved rather than silently normalize it. SI thenAMPure/elute50µl, second-strandE6111S16°C1h,AMPure/elute50µlEB. Select200ng for historicalTruePrepTD501 library preparation, preserving other stock; M280 capture of biotin-bridge fragments precedes PCR, thenAMPure size selection and HiSeqX-Ten2×150. Whole TruePrep components, capture fractions, indexed PCR oligos and final read landing remain unassembled.

The model must preserve source molecular contributions of nucleotides/bridge strands and actual selected stocks, and retain noncaptured/residual material separately. Entire sealed chimera sequence, linkage efficiencies/physical recovery and complete source-linked predecessor/vendor applicability remain open. No complete-reaction or source-absence acceptance.
