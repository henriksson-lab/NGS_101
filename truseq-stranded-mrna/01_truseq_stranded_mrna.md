# TruSeq Stranded mRNA

## Sources read and version scope

🟢 Illumina, *TruSeq Stranded mRNA Reference Guide*, document **1000000040498 v00, October 2017**. The source defines the workflow, optional controls and both HS/LS handling branches. Located readings include overview 174–258; capture/fragmentation and first strand 595–1004; SMM and cleanup 1005–1240; A-tailing and ligation 1245–1677; PCR and final cleanup 1775–1885; separately purchased components 2099–2135; inline controls 2490–2567. Line references are to the registered PDF's text twin.

🟢 Illumina, *Illumina Adapter Sequences*, document **1000000002694 v22, September 2025**, CD section 21090–21244. It prints both complete adapter templates and D-index members. This is positive source information, not a claim that an historical reagent lacked sequences. Its later revision does not establish the exact adapter revision in a 2017 preparation.

🟢 Illumina, *Index Adapters Pooling Guide*, **1000000041074 v17**, [TruSeq Kits](https://support-docs.illumina.com/SHARE/IndexAdaptersPooling/Content/TruSeqPoolingGuidelines.htm), local body 135–259 and footer 299 onward. This current page identifies compatible UD, CD and single-index component families, plate versus tube formats, and CD single- or dual-index use. Compatibility does not select a particular component, D member or acquisition recipe.

The full source-bundle oligo scraper was run before these readings. It is a reading aid; the sequence templates were checked in the original adapter-guide section.

## Selected workflow

The schematic selects the standard **HS handling branch**, **SuperScript II** from the listed reverse-transcriptase options, and omits optional inline DNA controls. Protoscript II is another option in the source. HS/LS handling and the separately purchased adapter component are distinct choices: the overview illustrates LS single indexing and describes HS dual indexing, while Appendix A names RNA UD, CD and single-index components.

1. Oligo-dT magnetic beads capture poly(A)-containing RNA twice. Supernatants and washes are discarded. Their composition is not equated to all and only non-polyadenylated RNA.
2. Elevated-temperature divalent-cation fragmentation and random priming occur during the second elution. The source does not establish a universal fragment length; Appendix B offers alternate fragmentation conditions.
3. Random-hexamer reverse transcription with FSA/Actinomycin D makes antisense first-strand cDNA.
4. SMM removes the RNA template and uses RNase H/DNA polymerase I with dUTP to make the marked second strand. The source states that the resulting cDNA is blunt-ended. **There is no additional standalone end-repair reaction between this step and A-tailing.** The inline-control table assigns end repair and 5′ phosphorylation to SMM.
5. AMPure purification retains cDNA, discards nonretained supernatant/washes and elutes the recovered material.
6. ATL adds one 3′ dA to the blunt ends. The source describes complementary single 3′-T adapter ends; the model shows that known polarity/base without assigning an exact historical adapter-core sequence.
7. LIG ligates the separately purchased **indexed adapters**. STL stops the ligation. Two sequential AMPure cleanup rounds precede PCR.
8. PPC anneals to adapter ends and PMM selectively enriches adapter-bearing fragments using the source's standard PCR program. Polymerase does not copy past dUTP. Indexes were installed by ligation; this is not an additional generic index-PCR recipe.
9. A final AMPure purification recovers the library. Quantification, normalization and pooling follow outside the construct-changing schematic.

The RNA degradation products, residual single-strand/under-adapted material and individual cleanup losses are not assigned complete physical removal. The guide says zero-adapter fragments do not hybridize to flow-cell primers and one-adapter fragments cannot form clusters; that functional exclusion is not a demonstrated physical discard at ligation or PCR. Marking the second strand quenches its amplification; it does not prove immediate destruction of the complete strand.

## Conditional CD reference and readout

🟢 The v22 source gives these full CD templates, with placeholder indexes as printed:

- i5 adapter: `AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT`
- i7 adapter: `GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7]ATCTCGTATGCCGTCTTCTGCTTG`

The preceding insert-derived dA is shown separately in the source's parenthesized example. D701–D712 and D501–D508 give 8-nt families. The molecular features identify i5/i7 as sample indexes; no individual D member, pairing or uniqueness claim is selected.

🟡 The schematic's complete **RNA CD reference duplex** uses those later source-positive templates as a compatible endpoint. Historical reagent-version assignment remains conditional, and constructor defaults preserve inferred styling. The biological interval is variable, not a literal 42-nt insert.

🟢 The workflow guide states that Read 1 maps antisense to the original RNA and Read 2 maps sense. It delegates cluster generation/sequencing to the chosen platform. The model has not selected acquired read lengths, single- versus paired-end acquisition, index masks or a physical index-2 priming mechanism. It therefore does not turn adapter compatibility into an acquisition claim or assign one universal index-2 primer.

The exact ordered PPC strands and their historical binding sites remain unassigned model data. This is not an accepted claim that no source prints them.

## Linked-source obligations

The review remains **in progress**, with no accepted missing-information claim. Fresh inspection of actual PDF annotations and text identifies:

- Pooling text names document 1000000041074, but the PDF's annotation also targets the older `truseq-library-prep-pooling-guide-15042173` landing page. The current v17 topic is read; comparison with that actual older linked target remains pending.
- The guide names a NextSeq/MiniSeq pooling bulletin, the sequencing-library qPCR guide 11322363, IEM guide 15031335 and quick card 15037152, Custom Protocol Selector, Local Run Manager/BaseSpace documentation, and a customer-sequence-letter endpoint. These targets are not counted as read merely from their titles or links.
- The checklist link literally contains `150XXXXX`; no exact completed checklist identifier or unavailable document is invented from that placeholder.
- Supporting platform dilution/denaturation/system guides depend on a platform selection that has not been made. That selection question is not an identified unavailable historical primer document.
- Remaining supporting sections, optional controls and alternative index families retain source-specific review boundaries. A read of the CD section does not exclude an entire vendor manual or support site.

The source checklist preserves earlier independent readings and records fresh bounded checks. Neither an unselected index nor a source-compatible reference endpoint closes the linked-source gate.
