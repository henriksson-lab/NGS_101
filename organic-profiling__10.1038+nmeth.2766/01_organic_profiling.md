# ORGANIC profiling — research notes

Kasinathan *et al.*, *Nature Methods* (2014), [10.1038/nmeth.2766](https://doi.org/10.1038/nmeth.2766).

## Sources and selected scope

🟢 The defining article's author manuscript (PMC3929178), published Supplementary Figure 1, and original deposited protocol in [GSE45672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45672) describe native MNase digestion and affinity purification. The graph selects the yeast 80 mM salt branch. The sources also describe yeast 150/600 mM conditions, fly profiling and mixing experiments; those are separate conditions.

🟢 The study metadata actually links `GSE45672_ORGANIC_profiling_protocol.rtf.gz`. The recovered author RTF names the study and authors, gives the yeast/fly/mixed-nuclei scope, and inherits yeast native immunoprecipitation from [Zentner *et al.* 2013, PLoS Genetics](https://doi.org/10.1371/journal.pgen.1003317), with the stated salt changes. It is a nuclei/IP procedure, not a library-preparation SOP. Word's `ADDIN EN.CITE` field artifact in the text extraction is not a biological claim.

🟢 The actual inherited 2013 original (PMC3585014, Methods 150–168) gives native calcium/MNase digestion, EDTA quenching, soluble S1/S2 extraction and pooling, an input aliquot, FLAG affinity purification, RNase A treatment, SDS/Proteinase K treatment, phenol extraction and ethanol DNA recovery. Its separate formaldehyde X-ChIP branch is not transferred to this native model.

🟢 Published Supplementary Figure 1 (MOESM254 PDF p. 3, visually inspected) independently shows the yeast and fly workflows. Its panel c gives the common input/IP library path: end repair → phenol/chloroform/S-300 cleanup → 3′ A addition → phenol/chloroform/S-300 cleanup → TruSeq adapter ligation → AMPure cleanup → PCR → AMPure cleanup → 25-cycle paired-end HiSeq 2000 sequencing.

## Historical library sources

🟢 The actual [Henikoff 2011 predecessor](https://doi.org/10.1073/pnas.1110731108), PMC3215028, is now recovered. Its Methods names Illumina **PE-930-1001**, omitted gel/Qiagen steps, phenol/chloroform plus S-300 cleanup, and AMPure steps after ligation and PCR. Its Results discuss the short-fragment recovery and older adapter-length contribution. Its directly linked detailed SI still returns a download gate/403 through the attempted original routes. Recovery of the main article does not establish that its older adapter configuration equals the later TruSeq example.

🟢 [Orsi *et al.* 2015](https://doi.org/10.1002/0471142727.mb2131s110), PMC4410783, Basic Protocol 3, supplies a later same-author ORGANIC library implementation. It was not cited by the 2014 article. The exact oligos and downstream enzyme choices illustrated on this page are attributed to this later procedure, not silently assigned to each 2014 sample.

The later source prints:

- Universal Adapter: `5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′`.
- Indexed Adapter: `5′-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT*G-3′`.
- P5 primer: `5′-AATGATACGGCGACCACCGAG-3′`.
- P7 primer: `5′-CAAGCAGAAGACGGCATACGAG-3′`.

🟢 The source marks the indexed strand's 5′ phosphate and both terminal phosphorothioate bonds; N6 denotes its multiplexing index. The exact sample-index member remains unassigned.

🟡 The complete 58-base Universal Adapter already includes the overlapping P5/Read 1 sequence. The local endpoint includes it once, followed by insert, fragment-derived dA and the complete 63-base Indexed Adapter. This corrects the repository's previous duplication of the four-base `ACAC` overlap. The correction concerns the model, not a publication error.

The Orsi illustration locally declares combined T4 DNA polymerase/Klenow/PNK end polishing, phenol/S-300 cleanups, M0212M Klenow exo− A-tailing, Rapid DNA ligase, positive bead cleanups and KAPA amplification. It does not introduce a separate kinase incubation. The initial Y-adapter duplex and its unpaired terminal dT are drawn from the source strands. The complete pre-PCR ligated Y molecule is not constructed, so its PCR geometry remains explicitly unresolved rather than substituting the amplified endpoint as a pre-PCR duplex.

## Readout and physical fates

🟢 The defining figure and all 24 instrument fields in the original GEO metadata identify HiSeq 2000; the defining figure specifies 25-cycle paired-end sequencing. The six-position Index 1 view belongs to the illustrated Orsi variant and is not a recovered 2014 instrument mask or sample assignment.

Both repaired fragment boundaries and insert span are useful observations. Soluble fraction selection, affinity capture/washing, DNA recovery and positive library cleanup have explicit physical fates in the graph. RNase/protease action is represented before DNA purification; the enzyme treatment itself is not called a physical discard. Mapping, fragment-length filters and occupancy analyses do not become molecular-removal steps.

## Remaining source review

The local PMC author manuscript lacks the referenced **Online Methods** section: its text 120–160 contains supplementary captions, not Methods. Earlier ledger entries that claimed those lines as Methods reads have been withdrawn and preserved as history. The original publisher entry currently supplies a subscription preview; the deposited RTF resolves the native-protocol inheritance but does not replace the complete published Online Methods or the cited 2011 SI.

The exact 2014 adapter/index assignment remains a model question. The 2011 SI, Online Methods and relevant historical vendor procedures remain outstanding originals; nuclei/fly/predecessor-linked sources have their actual identities and pending scope recorded in the review ledger. No claim that information was unpublished has passed the two-reader gate.

## Bounded curation round 2 (2026-10-16)

Inspected: Published original PMC3929178 JATS106–126; original MOESM254 PDF text174–280, SupplementaryFigure1.

No correction identified in inspected scope: native extraction/fraction pooling, RNase/proteinase digestion and purification, end repair/dA-tail/TruSeq ligation and25-cycle paired-end readout agree with existing qualified model.

Deferred: Rendered original SupplementaryFigure1 for garbled text; missing published Online Methods, exact historical adapter/index assignment and cited procedures; no preprint used. This capped pass does not establish a complete source review or publication-wide absence.

The rendered original Supplementary Figure1 (PDF page3) was also inspected because text extraction is partly garbled; it confirms the yeast/fly fraction handling and end repair → A-tail → TruSeq ligation → cleanup → PCR → cleanup →25-cycle paired sequencing order. Rendering this figure is no longer deferred.
