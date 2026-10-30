# 2b-RAD: published SOLiD/AlfI scope

Wang, Meyer, McKay and Matz2012, *2b-RAD: a simple and flexible method for genome-wide genotyping*, DOI10.1038/nmeth.2023.

## Original evidence

🟢 One published original read: SupplementaryFile MOESM591, title/DOI verified. Selected scopes: Table6 actual SOLiD/AlfI libraries; Table8 imperfect RTR enrichment; Table10 lines528–580 and originalp11 rendered; DetailedProtocol848–1029; Appendix1031–1067. Mandatory scraper found87 candidates in this SI; relevant oligos were checked against the actual text/table. Defining main DOI route found no open full original; softwareZIP is not treated as chemical evidence. This is a private incomplete model, not a finished schematic or missing-information audit.

## Oligos and noncanonical cap

🟢 AlfI long adapters are `CTCTCTATGGGCAGTCGGTGATNN` and `CTGCTGTACGGCCAAGGCGAANN`. Their terminalNN is the two-base **3′ cohesive-end** family, **not a UMI**. The ligation geometry now explicitly preserves each reported long-strand body and3′NN end via the neutral end constructor; this partial known end does not assemble the short cap, genomic junction or full sealed duplex. BsaXI NNN and RTR NNG are separate source branches. Table8 reports55.6% target specificity for selective adapters; the graph does not predict perfect restriction/ligation specificity or represent all non-target DNA as discarded.

🟢 Table10 short1 is written `3'-(invT)ACCCGTCAGCCACTA-5'`; short2 `3'-(invT)GCCGGTTCCGCTT-5'`. Reversing only the ordinary cores gives `ATCACCGACTGCCCA` and `TTCGCCTTGGCCG`, respectively. Each is followed by **one invertedT nucleotide cap**, making16 and14 nucleotide residues. Appendix writes `5'-ATCACCGACTGCCCAT(invT)-3'` and `5'-TTCGCCTTGGCCGT(invT)-3'`. Its postfix notation may annotate the terminalT or indicate addition; **full chemical equivalence and applicability to the Table10 form are unreviewed**. The model chooses the Table10 inventory of ordinary15/13nt cores plus one invertedT nucleotide, preserving both raw source forms. It does not add an invented extraT or erase the inverted base. Exact cap bond/polarity/duplex geometry remains unassembled; the selected Table10 class has blocked extendability and is not drawn as ordinary3′OH DNA. This source-form ambiguity is not asserted to be a publication error.

🟢 DetailedProtocol calls the four-primer mixture sld-p5/sld-p6 plus sld-p3/sld-p4; Appendix names them sld-5/6 and3/4. The selected complete primers are:

- sld-3: `CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT`
- sld-4: `CTGCCCCGGGTTCCTCATTCTCT[bc]CTGCTGTACGGCCAAGGCG`
- sld-5: `CCACTACGCCTCCGCTTTCCTCTCTATG`
- sld-6: `CTGCCCCGGGTTCCTCATTCT`

Table10 separately lists sld-p1/p2 as the shorter adapter-body primers. This is a source naming/set scope distinction, not an asserted erratum. Appendix/Table10 permit6–10nt sample barcodes; the study used6nt. The model represents an unassigned6nt sample-index family and does not instantiate a barcode or infer barcode acquisition cycles. It is not a UMI. The strict source matcher cannot instantiate the width of generic `[bc]`/`(bc)`. The coordinator independently checked the whole Table10/Appendix prefix, suffix and actual-study6nt width and authorized one narrowly documented p4-family source-transcription exemption. All five other exact source strands remain checked; capped short-strand bond geometry is deliberately unassembled. No source text or chemistry was altered, and this independent selected-source transcription review is not a missing-information audit.

## Source operation order and fractions

🟢 Prepare1µg intactDNA in4µl, add2µl bufferR/SAM/AlfI master mix,37°C1h, then65°C20min and ice. Split1µl for digest gelQC; the retained5µl goes to ligation. Independently combine the long and capped short oligos for each4µM adapter stock. Add20µl ligationmix to the5µl digest,16°C1h forAlfI, then hold onice. The source warns not to expose ligated material above65°C beforePCR. No cleanup occurs merely because restriction digestion produces short tags: flanking/uncutDNA and all reaction companions remain physical material.

🟢 Allocate4µl ligation to20µl testPCR and20µl to100µl preparatoryPCR; the unallocated physical remainder has unresolved fate. The testPCR's four5µl samples at6/8/10/12cycles go to analytical gels. Its products are not mixed into the productionPCR. Source Phusion four-primer PCR uses the minimum cycle count yielding a visible~130bp product, not a universal12cycle rule. Native templates, cap-bearing strands, unextended primers and PCR remnants are not assumed absent from the PCR mixture.

🟢 Load the entire prepPCR on2%gel and recover~130bp band while avoiding high-molecular-weight material and≤90bp dimers. The selected water-elution branch adds40µl water to the gel slice, spins2min, holds4°C overnight and transfers~50µl eluent next day to a new tube. Gel material is physically separated at collection. Commercial gel extraction is an explicit alternative, not the chosen route. Table6 verifies actual SOLiD/AlfI acquisition; exact run cycles/bead enrichment/read-primer geometry remain unreviewed.

## Conditional reference and boundaries

🟡 The optional two-sided-cut helper requires explicit complementary top/bottomDNA and **four caller-supplied cut coordinates** in each strand's5′→3′coordinate system. It returns left/tag/right fragments for both strands and retains the complete nucleotide inventory. The reference is inferred by default and rejects coordinate choices whose selected central strands have no shared duplex interval after mapping the bottom coordinates onto the top axis. It is a selected cleavage reference, not an AlfI recognition/offset/tag-length predictor. Finite test strings are illustrative inferred references, not a published genome sequence or evidence for cleavage yield. Native genomic bases, digestion fragments, PCR copies and reported genotype calls remain separate.

The model covers the source-known SOLiD/AlfI operational route while complete cap/end/ligation/PCR and acquisition geometry is unresolved. It imports no BsaXI/Illumina/RTR or modern kit recipe. No source-supported parent/component inheritance edge has been accepted from this one-original scope. Defining full main, full figure scope and exact historical SOLiD chemistry remain later original-reading targets, not confirmed absent information.

## Scoped verification

The16 focused checks pass, including strand conservation, overlapping duplex support, cap inventory/blocking and separate analytical/preparative fates. Source guard checks five complete literal source strands, with one independently authorized unassigned p4-family macro exemption and two capped strands left without invented full canonical geometry. Chemistry13unresolved/0contradictory; source-review schema validates. The transferred, source-ready eluate is a library endpoint, while its complete molecular/run geometry remains unresolved. This verification is an implementation review, not a full-source or missing-information audit.
