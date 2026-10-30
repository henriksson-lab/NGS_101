# Mint-ChIP — bounded source notes

Defining DOI: https://doi.org/10.1016/j.molcel.2015.11.003

Associated author manuscript PMC4707994 main70–86/134–170; one original. Published-final equivalence not established.

🟢 Native cells are lysed/digested at separately calibrated MNase concentrations; low- and high-MNase samples form separate pools. Production uses20ul cellsPBS+20ul2xlysis, ice20min/37°C10min. Carrier cells are digested concurrently, up to400000 perChIP, and receive water/EGTA/butyrate instead of T7adapter ligation mix.

🟢 Double-stranded T7 adapters contain T7promoter/SBS3/sample index#1 and a5′C3 spacer preventing self-ligation (main72). C3 is a non-nucleotide block, not three cytosines. Main146 cites actual TableS2 for whole sequences; gated supplement is not an inspected oligo table. Stock500uM/annealed250uM/frozen−20°C; source ligation descriptor4ulannealed adapter(12.5uM), 30mMEGTA/5mMbutyrate/1xkitbuffer/1.5mMATP/1xNTPmix/1.6ulEndRepairmix/1.6ulFastLink per80ul; RT2h then80ul dilutionbuffer. The printed NTPmix is retained without guessing sugar composition.

🟢 Pool indexed chromatin at the same MNase concentration, add unlabeled carrier, split four ChIPs: H3(5ulAbcam1791), H3K4me3(1ulMillipore07-473), H3K27ac(2ulActiveMotif39133), H3K27me3(5ulMillipore07-449). Overnight4°C then25ulProteinG2–4h4°C, five200ul washes(2RIPA/1highsalt360mMNaCl/1LiCl250mM/1TE), elute62°C1h with.25mgmlPK;1xSPRI>100bp.

🟢 HiScribeT7Quick IVT37°C2–16h;1ulDNase37°C15min thenSilaneRNAcleanup. RNA amplification does not imply perfect loss of every untagged carrier molecule. SSIII/50uMRTprimer(randomhexamer+SBS12) then1xSPRI>100bp. Randomhexamer is a priming class, not a UMI.

🟢 PCR1 forward1/reverse1 .2uM each/PfuUltraIIFusionHS,up to12cycles(usually6);1xSPRI>100bp. Select2–40ng intoPCR2 forward2/reverse2 .2uM, up to12cycles(usually8), reverse2 has assay index#2. DoubleSPRI.6x/.95x200–700bp; analyticalQubit/BioAnalyzer/KAPAquant aliquots; paired-endNextSeq500/HiSeq2500. Sample index#1 and antibody-assay index#2 differ in meaning. Exact runmask/whole primer landing unassembled.

🟡 Local private partial graph preserves actual independent carrier factory, source same-MNase pool, four assay allocations, IVT→DNase→RNAcleanup→RT order, actualPCR2/QC allocations, per-assay selected pool contributions and untransferred library-stock remainders, then actual assay-library pool. Whole duplex geometry, closed nucleotide/catalyst inventories, exact adapter junctions and complete final molecules remain unassembled. No predecessor recipe or confirmed absence is asserted.

Source investigation closed within the selected three-original ceiling before model construction. This capped positive read is not a complete source/reaction/absence audit. Missing-information claims are empty.

Deferred originals and applicability questions:

- Historical EndIt/FastLink/HiScribe/SSIII/PfuUltra original vendor recipes unread at bounded limit; kit names do not authorize whole component inheritance.
