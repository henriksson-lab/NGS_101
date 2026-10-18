# G&T-seq — research notes

## Sources read

- 🟢 Macaulay et al. 2015, DOI [10.1038/nmeth.3370](https://doi.org/10.1038/nmeth.3370), [accepted author manuscript with Online Methods (Zenodo record 895966)](https://zenodo.org/records/895966), linked from Europe PMC as the free full text for PMID 25915121. Online Methods "Cell lysis, cDNA isolation and amplification", "Genomic DNA precipitation and amplification" and "Library preparation and sequencing" read at extracted lines 309–367. The supplement workbooks and figures were also checked; the large MOESM220 workbook is Fig. 1b/1c measurement data, not a laboratory protocol. An earlier pass read a ResearchGate author upload through a web viewer (its line numbers are no longer used); that copy returns HTTP 403 to the fetcher.
- 🟢 Macaulay et al. 2016, [authors' detailed G&T-seq protocol](https://www.pure.ed.ac.uk/ws/files/181359881/Separation_and_parallel_sequencing_of_the_genomes_and_transcriptomes_of_single_cells_using_G_T_seq.pdf), DOI 10.1038/nprot.2016.138: oligonucleotide list at extracted lines 644–658 and capture-bead preparation at lines 718–734.

## Molecular order

🟢 The original 2015 Online Methods specify this order: a lysed single cell is mixed with a biotinylated oligo(dT) primer. Poly(A) RNA is
captured on streptavidin beads and magnetically separated from the genomic-DNA-containing
supernatant. The bead fraction undergoes Smart-seq2-like reverse transcription,
template-switching and WTA. The DNA fraction undergoes either MDA or PicoPLEX WGA. The
two amplified fractions are converted into separate indexed sequencing libraries.

## Protocol boundary

MDA and PicoPLEX are alternative genome-amplification choices tested by the paper. The
transcriptome and genome products never become one chimeric construct. Their common cell
identity is inherited from physical separation of a single-cell lysate.

🟢 The 2015 Online Methods print the capture/RT primer as
`5′-biotin-triethyleneglycol-AAGCAGTGGTATCAACGCAGAGTAC-T30-VN-3′` (accepted-manuscript lines 319–320),
the TSO `5′-AAGCAGTGGTATCAACGCAGAGTACrGrG+G-3′` (lines 336–337), and the
PCR primer `5′-AAGCAGTGGTATCAACGCAGAGT-3′` (lines 341–342). The method names
MDA/PicoPLEX as alternative WGA branches (lines 348–351) and Nextera XT for
the cDNA and gDNA libraries (lines 358–360).

🟢 The authors' 2016 detailed protocol separately prints the capture/RT primer verbatim as
`5′-Biotin-TEG-AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN-3′`, with a 30-base dT tract and
3′ `VN` anchor (PDF text lines 644–651). It also prints the TSO
`5′-AAGCAGTGGTATCAACGCAGAGTACATrGrG+G-3′` and ISPCR
`5′-AAGCAGTGGTATCAACGCAGAGT-3′` (lines 652–658). The same protocol names
Nextera XT library preparation and PicoPLEX or MDA as alternative WGA branches
(lines 171–185, 421–440, 1373–1387).

The 2015 TSO ends `…AGTACrGrG+G`; the 2016 author protocol adds `AT` before
`rGrG+G`, giving `…AGTACATrGrG+G`. This is a **source-version difference**: the
public 2015 page uses the defining paper's shorter TSO. The 2016 recipe is a separately scoped later version and does not replace the defining2015 oligo. Kit-specific Nextera XT index
choices and WGA reagent oligos remain unresolved for a selected branch; the
public page retains explicit historical kit boundaries for those branches rather than assigning generic final-library geometry.

## Bounded pass 2026-10-17 (round breadth-first-10)

The 2015 author-uploaded Online Methods PDF returned HTTP 403 to `tools/get_sources.py`
and Europe PMC lists no open copy (PMID 25915121); its 2015 oligo values above rest on
the earlier web-viewer reading and were not re-read in this pass. Re-read here: the 2016
author protocol oligo list (lines 644–658), separation steps 17–28 (lines 782–884, washes
pooled into the gDNA plate), RT on beads (step 29, 42/50 °C then 60 °C inactivation),
ISPCR/KAPA HiFi WTA (steps 31–34) and AMPure cleanups (steps 35–54); and the 2015
figure supplement MOESM222 (legends only, no chemistry).

🟡 Model correction: the template-switch panel previously used a generic shared scene
with a 14-nt symbolic handle, dT14 and an `rGrGrG` TSO end. Both originals print an
`rGrG+G` end (two RNA G and one LNA G). The local panel now draws the printed capture
oligo priming first-strand cDNA and the 2015 TSO; the untemplated C additions that
pair with `rGrG+G` remain an inferred mechanism of unprinted length.

Deferred: re-reading the 2015 Online Methods (lines 554–623) once a copy is in the
bundle; 2016-only details (wash pooling, AMPure ratios, KAPA HiFi) are not imported into
the 2015 graph; Nextera XT index members and WGA kit primers; Smart-seq2 predecessor.

## Bounded pass 2026-10-17 (round deferred-followup-02)

The defining 2015 text is now in the local bundle: an accepted author manuscript
deposited on Zenodo (record 895966, DOI-linked; reached from the Europe PMC full-text
link). It is the peer-reviewed accepted version, not the typeset article; its title
line reads "G&T-seq: Separation and parallel sequencing of the genomes and
transcriptomes of single cells". Read: Online Methods lines 309–367.

🟢 Verified there: biotin-TEG oligo(dT)30VN on streptavidin Dynabeads, 20 min capture
(lines 318–324); four washes in 50 mM Tris-HCl pH 8.3, 75 mM KCl, 3 mM MgCl2, 10 mM
DTT, 0.5% Tween 20, 0.2× SUPERase-In, each wash pooled with the gDNA supernatant
(lines 325–332); on-bead RT with SuperScript II, betaine, MgCl2 and the TSO
`AAGCAGTGGTATCAACGCAGAGTACrGrG+G`, 42 °C 60 min, 50 °C 30 min, 60 °C 10 min
(lines 333–340); KAPA HiFi HotStart ReadyMix with the 23-nt PCR primer, 18 cycles
(lines 341–345); cDNA AMPure 1:1; gDNA precipitated on AMPure 0.6× and eluted
directly into GenomiPhi V2 (MDA) or PicoPlex; WGA cleanup AMPure 1:1 (lines 345–356);
1–5 ng into Nextera XT (lines 358–360). The 2015 TSO lacks the `AT` printed in the
2016 protocol, as already shown on the page. The wash, RT and PCR conditions above
are therefore 2015 facts, no longer only 2016 details.

Corrections: source locators in the module and note now point to the bundled
manuscript; the RT, WTA PCR and WGA enzyme inputs name the 2015 reagents; the TSO
source-oligo exception was removed because the 2015 TSO is now verified verbatim in
the bundle.

Deferred: typeset 2015 article not compared with the manuscript; Nextera XT index
members and historical guide; GenomiPhi V2 / PicoPlex kit primers; Picelli 2014
Smart-seq2 predecessor; WTA PCR binding geometry still unplaced in the generic cDNA state.

## Native graph migration — 2026-10-18 batch02

This pass re-read one original: the defining2015 accepted manuscript, Online Methods310–368. The2016 oligo/reagent comparison remains historical evidence from earlier bounded checks; it is not substituted. Source capture includes the pooled genomic supernatant, all four bead washes and10µl tipwash. RT/WTA has actual single-strand first/switched cDNA and one physical ISPCR species; the staged graph is a mechanistic decomposition of the source reaction, not a claim that TSO was added in a separate tube. Untemplated C additions/copy extent remain inferred. GenomiPhiV2 MDA and PicoPLEX are alternatives. Historical NexteraXT native adapter/index/PCR geometry remains unassigned; no modern kit or generic shell is substituted. Run platform is source-named HiSeq2500 fast mode or MiSeq; acquired read roles, cycles and selected physical index members remain unassigned.

The local graph owns the named reagents, source order, physical branch decisions, primers, endpoint and read declarations. Independent molecular populations occupy separate Scene rows; subsequent RNA-template fate and complete native end geometry remain explicit model boundaries. Curator-written review targets are local and are not derived from graph contents. This is a capped migration pass, not a completed missing-information audit; no missing-information claim is accepted.
