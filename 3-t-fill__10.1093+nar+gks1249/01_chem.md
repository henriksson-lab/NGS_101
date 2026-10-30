# 3′T-fill — capped source-qualified notes

Wilkening et al. (2013), An efficient method for genome-wide polyadenylation site mapping and RNA quantification. DOI10.1093/nar/gks1249; corrigendum DOI10.1093/nar/gkt364.

🟢 Published main PMC3597643, Figure1/Methods42–59, Supplementary Table1 and formal corrigendum PMC3695521 were inspected. This consumes three originals; the detailed supplement containing Methods/TableS3 remains available but unreviewed. Mandatory bundle primer scraping completed; hits in that unread document do not become adopted sequences.

🟢 The formal correction DOI10.1093/nar/gkt364 changes YPD glucose and YPGal galactose from1% to2% in the culture Methods. It does not correct a library oligo or T-fill reagent. The selected molecular graph starts from isolated DNase-treated RNA and does not perpetuate the erroneous growth composition.

🟢 Source path: 10µg total RNA with independently prepared polyadenylated Bacillus IVT controls → Mg-acetate fragmentation80°C5min →1.5×AMPure/12.8µl →11.2µl into biotinylated P5_dT16VN-primed SSII RT42°C50min →1.5×AMPure →RNaseH/DNApolI second strand16°C2.5h →0.8×AMPure →M280 streptavidin capture and washes →NEBNext end repair/wash →Klenow exo-minus dA tail/wash →P7_T1_Mpx/T4DNA ligase adapter ligation/wash →50µl bead suspension →24µl for18cycle Phusion PCR →1.8×AMPure →300bp eGel size selection. The1.6µl pre-RT remainder and26µl pre-PCR bead remainder remain unresolved physical material.

🟢 The cBot priming buffer is replaced with unlabeled dTTP, genomic sequencing primerV2, non-hot-start Phusion and TaqE; 60°C5min primes and fills through the poly(A) tract before labeled acquisition. Both polymerases remain declared; the main says Taq is crucial for complete filling. This is a preparative on-cluster extension, not a computational trimming trick or a relabeling of an ordinary poly(A)-containing read. Final acquisition is HiSeq2000; the barcode is read from the opposite paired-end direction. Exact run lengths/index coordinates and historical primer sequences are unassembled in this capped pass.

🟡 P5_dT16VN, P7_T1_Mpx, PE1.0/PE2.0 and genomicV2 are retained as source-named inputs. No whole sequence is transcribed from the unreviewed TableS3 or filled from a later kit by name. No derived construction/Scene claims perfect T-fill, exact polymerase stop positions or capture purity. Mapping-quality/poly(A)/internal-priming filters occur after sequencing and do not physically remove library molecules. Private incomplete graph only; no accepted absence claim.

Open source scope:

- PublishedFile006.doc detailed Methods/TableS3: available, unreviewed because main/Table1/correction used the three-original cap.
- Historical P5_dT16VN/P7_T1_Mpx/PE1.0/PE2.0/genomicV2 complete sequence and endpoint geometry.
- Source-listed alternative3primeBpmI/Internal preparations are comparison branches, not extra counted protocols here.
