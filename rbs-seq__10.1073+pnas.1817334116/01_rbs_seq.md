# RBS-seq

Khoddami et al.(2019), *Transcriptome-wide profiling of multiple RNA modifications simultaneously at single-base resolution*. DOI https://doi.org/10.1073/pnas.1817334116 .

## Source identity and correction boundary

🟢 Published non-manuscript PMC6452723 main, provenance57, Fig1/chemical rationale63–89, Methods107–131; published SI original pnas.1817334116.sapp.pdf, methodspp3–6 and mechanistic-assay discussionpp9–11. Full-bundle scraper run; adapter-trimming strings and validation primers are not substituted for actual library oligos. Other data workbooks and DatasetS10 original primer/oligo document remain deferred under the three-original cap.

The actual linked correction is *Correction for Khoddami et al., Transcriptome-wide profiling of multiple RNA modifications simultaneously at single-base resolution*, PNAS2019;116(18):9136, DOI10.1073/pnas.1905628116, PMID31010919/PMC6500115. Identity is verified from PubMed's erratum record only; PRINTED CORRECTION CONTENTS WERE NOT OBTAINED OR READ. Standard DOI retrieval gave EBI500 and no article content from PMC; normal PMC browsing showed recaptcha, and the PNAS DOI page returned403. No bypass attempted. All chemical/library statements below are a provisional main/SI-qualified model pending that correction; its contents are not assumed immaterial. Explicit uncovered correction/source inventory prevents finished publication.

## Distinct physical and copying paths

🟢 Complete BS/desulfonation workflow converts selected surviving RNA C to U, retains selected m5C, converts m1A to m6A by Dimroth rearrangement, and creates a stable pseudouridine-monobisulfite adduct. These changes do not replace RNA with its reported DNA base. A typed multi-site RNA helper preserves original residue identity and chemical state, allowing multiple modifications on one explicitly supplied synthetic molecule.

Copying RNA U gives complementary DNA A and reference-oriented T; copying m5C gives DNA G/referenceC; m6A derived from m1A gives DNA T/referenceA. Native m1A copy choices remain caller-supplied, not a deterministic mismatch. The selected single-Psi bypass path removes a position from cDNA, while the modified RNA residue and RNA length remain intact. Source also reports1–2base deletion signatures; the helper's single-site case does not predict every multi-base/ambiguous alignment outcome.

Source experiments show BS treatment alone does not create the characteristic Psi deletion signal: subsequent heated/Mg-containing library handling matters. The helper explicitly qualifies a selected heat/Mg-conditioned bypass path. Exact RNA sugar conformation and the Mg-assisted reorientation model remain mechanistic inference; source nucleoside structural data are not silently treated as a complete structure for every polymeric RNA site.

## Actual source workflow through sequencing

Three RNA fractions (small RNA, fragmented rRNA-depleted total RNA, fragmented poly(A) RNA) were prepared/processed separately; collective analysis is not a physical pooling step. About one-third of each was retained as NBS control. BS uses formamide95°C5min, then312µl5Mbisulfite stock plus3µl100mMhydroquinone stock in600µl total,50°C16h; the stock is not called a5M final reaction. NAP-10 cleanup is followed by mixing0.5ml sample with0.5ml2MTris stock pH9,37°C2h, then source ethanol recovery. The exact SI mixing instruction is retained rather than calling the diluted solution2M.

Antarctic phosphatase37°C30min/65°C5min precedes PNK+ATP37°C60min and RNeasy cleanup. Source v1.5 small-RNA3′adapter preannealing70°C2min precedes truncatedT4RNA-ligase2 reaction22°C1h with addedMgCl2. Separately heated5′adapter is ligated with T4RNA ligase20°C1h. SSII with source barcodedRTprimers, Phusion15cycles, AMPure cleanup, NaOH denaturation and HiSeq2000PE101 follow. The graph carries these known operations, including both BS/NBS branches, through reported sequencing metadata. Actual historical adapter/polymer/end/primer/index/run geometry and allocation of each relevant heat/Mg exposure remain unassembled; no later BID-seq buffer recipe is borrowed.

## Deferrals and lineage

Actual correction targets: https://doi.org/10.1073/pnas.1905628116 , https://pmc.ncbi.nlm.nih.gov/articles/PMC6500115/ and https://www.pnas.org/doi/10.1073/pnas.1905628116 . Correction relevance remains blocked until its original contents are read. DatasetS10, linked historical RNA-bisulfite/small-RNA originals and independent source review remain open. BID-seq explicitly cites RBS-seq for its Psi-adduct/deletion principle but changes chemistry and cytosine-conversion behavior; that supports a scoped RBS→BID component-strategy edge, not whole-library equivalence. No confirmed absent-information claim or finished schematic.

## Bounded physical-array correction with correction/geometry holds retained

The same closed published main and SI support actual molecular formamide,
sodium bisulfite and hydroquinone inputs, rather than condition-only labels.
The NAP-10 1 ml eluate is now physically divided into two 0.5 ml portions
before separate Tris desulfonation. Each actual preparation receives acetate
and GlycoBlue and is divided into three 377 µl precipitation portions. Six
independent precipitation/recovery operations each produce a 10 µl RNA
contribution, which are physically pooled afterwards into the source 60 µl
stock. This does not clone one completed reaction into supposed replicates.
Curator-written tube target inventories are declared independently of graph
node enumeration.

Thirteen focused checks pass. These physical corrections do not waive the
actual unread correction or complete library-geometry obligations. Both
`unread_correction` and `unreviewed_library_geometry` remain uncovered. The
previous actual correction HTTP/reCAPTCHA records remain history; no retry,
new failure or judgment that the correction is irrelevant was made. Remaining
RNA-fraction/control, adapter, RT/PCR stock, QA and full read-assembly scope
is still incomplete. There is no PAGE_FILE builder or finished-page claim.

Physical RNA residues and adducts retain their origins. A selected complementary
cDNA omission does not delete a Ψ-BS residue from RNA, and a native m1A copying
outcome is not silently predicted. No source absence claim follows from the
bounded pass.
