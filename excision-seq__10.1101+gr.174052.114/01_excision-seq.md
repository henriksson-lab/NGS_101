# Excision-seq

Defining DOI: https://doi.org/10.1101/gr.174052.114

Bryan et al. (2014), High resolution mapping of modified DNA nucleobases using excision repair enzymes.

## Source scope

🟢 Published open-access PMC4158761 main article (manifest author manuscript=False), Results Development of Excision-seq/Figure1 (text50–62), uracil post-digestion interpretation70, UV branches/Figure4 (78–88), Methods110–122. Source explicitly refers detailed recipes and oligonucleotides to Supplemental Methods. The fetched EuropePMC supplementary ZIP was inspected for original filenames; it contains only eight article figure JPG/GIFs, not that Methods document. This does not establish that the Methods are absent from publication or unreachable elsewhere.

## Two distinct uracil branches

🟢 Predigestion: UDG removes uracil, generating abasic sites; the main text names “T4 endonuclease IV” for gap production. That printed enzyme name is preserved as source wording and not silently substituted with another enzyme or a mechanistic recipe. At sufficiently high modification density, double-stranded fragments can be recovered and converted by end repair, A tailing, adapter ligation and PCR. Positions immediately upstream of mapped5′ read ends report removed uracil (referenceT), rather than sequencing a retained U base.

🟢 Post-digestion: mechanically shear genomic DNA and prepare adapter-ligated molecules before excision. Enzyme treatment interrupts strands containing the modified base; intact strands remain PCR-competent. Broken DNA and excision products are physically still in the tube unless removed by an actual documented workup. Lower sequenced coverage is an amplification-selection effect, not evidence that broken molecules were physically discarded at digestion. One surviving intact strand may still seed amplification. Relative coverage is inversely related to uracil content at fragment scale, not an exact basewise lesion call.

## Other source branches

🟢 For UV products, UVDE cleaves UV-damaged DNA at lesion-associated positions and separate CPD or6-4PP photolyases restore the corresponding damaged bases before ordinary library operations. Those enzyme conditions and lesion-end geometries require Supplemental Methods; the partial graph does not merge this route into uracil cleavage or borrow CPD-seq's later T4endoV/APE1 procedure.

The source sequences libraries on Illumina MiSeq or HiSeq2000. Read length, exact sequencing/adapter polymers, index architecture and reaction conditions remain unassembled in this pass. The main text shows chemical U→abasic→gap classes without fixing genomic coordinates or assuming complete conversion. A source-qualified private graph represents pre/post uracil branches and leaves the full library recipe unresolved. The preview has no finished marker.

## Deferred sources and lineage

The actual published Supplemental Methods original is the primary follow-up target; its attachment identity is not resolved from the current PMC export. Desired adapter/end specifications remain open questions, not fictional unavailable vendor documents. Defining main text cites enzyme predecessors and LM-PCR as context; neither a familiar base-excision enzyme nor similar library ends establishes parent chemistry. No component-specific existing-catalogue parent is accepted from this pass. No second-reader missing-information audit occurred.
