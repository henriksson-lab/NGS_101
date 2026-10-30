# MOLseq — optically controlled ligation for spatial RNA sequencing

Venkatramani et al., *Multiplexed optical barcoding and sequencing for spatial omics*, [DOI10.1038/s41598-026-41186-y](https://doi.org/10.1038/s41598-026-41186-y).

## Sources and selected boundary

🟢 This capped pass read the published main PMC13136315: concept/Fig1–2 at54–70, co-culture100–110 and Methods124–168. The original MethodsPDFp9 was rendered for the literal stock-buffer unit. Published SI Figs1/2/4 and Table1 originalp6–8 were inspected, including rendered whole modified oligo rows. Both originals were selected before opening; no preprint or superseded article-in-press original was adopted. The mandatory primer scraper ran and its source pointers were checked. Source investigation is closed at two originals, not a full linked-source/absence audit.

🟢 The selected RNA-sequencing route uses polyT-primed in-situ RT in fixed, adjacent U2OS/MEF regions. Subsequent light-controlled additions assemble spatial barcode letters at the cDNA5′ end. Photolysis exposes a phosphate on the retained DNA letter, permitting a further splinted ligation. This differs from protein-reporter cell sorting, photocrosslinking a barcode to RNA, or joining staggered duplex ends. The locked-T size-measurement primers, no-RT TapeStation controls and fluorescent detection letter variants are separate source components, not production RNA-seq primers.

## Whole source orderings and non-nucleotide linker

🟢 The RT primer is source5P–handle–N20–T30; N20 purpose remains unassigned, with no automatic UMI designation. The model stores the whole modified ordering and checks the original Table1 row after removing only its interleaved column label. Table1 supplies the TSO with two rG residues and terminal LNA-G, four letters, four splints, SmartSeq3 tag, LD reverse primer and both5′-enrichment primers. The model imports the exact canonical Nextera segments used by the latter; it does not import a Nextera or SmartSeq reaction constructor.

🟢 Each non-fluorescent source letter contains a DNA cap, `/iSpPC/` and a retained DNA body. The local typed `PhotoLetter` preserves their order and counts only DNA nucleotides. The linker is not replaced by a string of Ns or counted as bases. Its cleaved-class descriptor retains the downstream body with5P and records the separated cap DNA; detailed linker-fragment chemistry, intact/damaged populations and quantitative conversion remain unassembled. Spatial-identifier metadata is preserved on the retained letter. Full linker-bearing backbone SVG geometry is not claimed.

🟡 Exact local source-component interfaces are computed for L2→the initial capture handle with Splint21 and L4→L2 with Splint42. Scene rejects a supplied splint that does not pair to the incoming3′ five bases and prior5′ ten bases. These are conditional local interfaces using supplied source components. They do not establish the complete chronological co-culture schedule, pool occupancy, ROI pattern or final whole cDNA. The graph therefore retains the full spatial barcoding phase as an explicitly opaque operation. Unknown occupied splint members or complete schedules are open modelling questions, not invented unavailable disclosure documents.

## Physical preparation and source values

🟢 Source plated cells receive4%PFA/DPBS, the polyT primer1µM/2×SSC overnight37°C, then actual formamide/SSC, DPBS, water and SSC washes. RT supplies MaximaHMinus, TSO, molecular dNTP substrate and two distinct RNase inhibitors. The original dNTP concentration/stock wording is retained without selecting a per-base-versus-total interpretation. Primer components retained after washing are not supplied as a second fresh primer addition. Actual source wash fractions remain physical workups, not computational specificity guarantees.

🟢 Per-round letter/splint doses, QuickLigase, dark incubation, patterned365nm irradiation and washes are source-qualified, while detailed round scheduling and full products remain opaque. Actual published mainXML and originalPDFp9 both print **132nM Tris-HCl in the2×ligation stock**. The model retains that literal value; its intended concentration remains unresolved. It is not silently repaired from a familiar ligation-buffer recipe.

🟢 Source SmartSeq3-tag ligation precedes heated Triton/water extraction. The aspirated DNA-containing solution and unaspirated cell/surface residue have distinct physical fates. LD-PCR uses the source tag/TSO-related primers and SeqAmp. It is followed by0.8×AMPure cleanup, historical NexteraXT tagmentation, custom5′-enrichment PCR, another0.8×cleanup, index-adding PCR and0.8×cleanup. Each represented PCR uses components from one preparation stock, including molecular nucleotides; no separate fresh enzyme is invented from a buffer name. Source kit versions, reaction programmes, complete adapters, sealed molecules and exact primer landing remain unassembled where not derived.

## Sequencing and remaining scope

🟢 The main reports MiSeqReagentKitv3/MS-102-3001 with a150-cycle label and also says150paired-endreads. Both statements remain explicit and the exact run mask is unassigned; no current kit specification silently replaces the source. The source filters Read1 for barcodes and aligns corresponding Read2 separately. Alignment filters and count generation are computation, not removal of molecules from the experimental preparation. No perfect ligation, photolysis, targeting or species-assignment rate is claimed.

All selected originals were accessible. Remaining full molecular schedules, occupied splints, original source interpretations and directly linked predecessor/vendor details are unreviewed or unassembled, not confirmed absent. This is a source-positive private partial graph, with no finished public schematic or full reaction-validation claim. A narrow independent implementation recheck is recorded separately from the two-reader missing-information gate.
