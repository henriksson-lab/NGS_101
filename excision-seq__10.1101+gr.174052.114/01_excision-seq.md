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


## Public-page source-identity hold, 2026-10-30

The actual normal graph still requires `supplemental_methods`, which remains
unaccounted for, as well as a sequencing terminal. This source obligation was
neither excluded nor marked covered by the main article's coarse libraries.
No model code or publication gate changed.

The inspected defining original is the published article
[PMC4158761](https://pmc.ncbi.nlm.nih.gov/articles/PMC4158761/)
([defining DOI](https://doi.org/10.1101/gr.174052.114)),
`gr.174052.114_PMC4158761.xml(.txt)` and its originalFigure1, available in the
article PDF and alternate article-figure ZIP. The actual main version is
CCBY/non-author-manuscript, not a substitute preprint. Results50–62/Figure1
and Methods110–122 were read; the originalFig1JPG was inspected directly.

The main explicitly assigns strains, oligonucleotides and detailed library
recipes to **Supplemental Methods**. The identity of that actual cited
attachment remains unresolved in the local export. The particular downloaded
EuropePMCZIP contains only these article-figure representations:
`1534fig1.jpg`, `1534fig1.gif`, `1534fig2.jpg`, `1534fig2.gif`,
`1534fig3.jpg`, `1534fig3.gif`, `1534fig4.jpg`, `1534fig4.gif`. It is not the
Supplemental Methods. Local mainJATSlink attributes provide the defining
article DOI link but no exact methods attachment URI. This does not establish
that the methods are absent or unavailable elsewhere. No PDF/DOCXattachment
path is guessed, and the article DOI is not reported as a failed known
methods-download URL. Later authorized retrieval must resolve the actual
source-file identity.

The main supports actual predigestion and post-digestion DNA libraries,
including end polishing/adaptor ligation/PCR, sequenced on MiSeq or HiSeq2000.
ActualFig1predigestion also depicts chew-back/fill-in before polishing and
adaptor ligation. Specific enzymes/end states, oligos and reaction/workup
fates still require the named methods original; raw run masks are not
inferred from platform names. The main repeatedly prints “T4endonucleaseIV”,
while the diagram labels “UDG/EndoIV”. That literal name and unresolved
identity are retained. Neither E.coliEndoIV nor T4EndoV is silently substituted,
and no confirmed publication-error conclusion is made.

The physical input U denotes **deoxyuridine in DNA**, not RNA. UDG releases
free uracil nucleobase and leaves abasic DNA; subsequent cleavage/chew-back/
fill-in is not a direct U→Tconversion. In predigestion, the mapped referenceT
immediately upstream of a readstart represents the former uracil site; it is
not a retainedUbase in the sequence read. Post-digestion coverage reflects
relative survival/amplification: broken strands remain physical material
even when unable to seed the represented PCR. Species/sample and pre/post
preparations, rather than a fictitious universal allocation, remain open
where their detailed source identity is unresolved.

The bounded preflight ran the mandatory bundle scraper
(`/tmp/unlinked-excision9-scraper.log`), inspected actual source/figure/ZIP
identities and confirmed the actual normal gate failure. Zero oligo hits
is not an absence conclusion. Two existing source-boundary/no-marker checks
passed; chemistry8unresolved/0contradictory. No new original was fetched, and
tests/builds were not repeated for this documentation-only adjudication.
