# Spatial-ATAC-Hi-C — bounded source notes

Defining DOI: https://doi.org/10.1038/s41592-026-03217-4

Published main PMC13645645 is CC-BY-NC-ND, open access, not an author manuscript; title/DOI/body verified. Published SI1 contains actual Tables1–2 (originalPDFp11–13 rendered/inspected). These are two originals. Source investigation closed before implementation; no third original needed for this selected partial scope. Preprint and other supplementary workbooks are excluded/unread, not absence authorities.

🟢 Sourcemain243–281: fresh-frozen7–10um tissue sections are thawedRT10min, fixed2%FA20min, quenched.2Mfinalglycine5min,2x1mlDPBSwashed.30ullysisbuffer4°C20min; source then literally adds36ul10%SDS62°C10min and30ul1%Triton37°C15min. Same-original publishedPDFp15 confirms36ul10%SDS; it is retained, not silently changed to a familiar recipe. Pre-tagmentation wash lists1MTris7.4/5MNaCl/1MMgCl2/10%Tween20/5%BSA/H2O without component volumes; these are kept as printed descriptors with final formulation unassembled, not treated as modelled final concentrations or a proven publication error.

🟢 Restrictionmix20.5ul:11ul10xCutSmart/1ul10%Triton/8ulDpnII/1.5ulHinfI37°C1h;2x200ulwater. Endrepair50ul:20ul.4mMbiotin-14-dATP,1ul10mM eachdCTP/dGTP/dTTP,5ul5UulKlenowM0210L,5ul10xCutSmart/17ulwater;RT1h. Contactligation100ul:10ul10xT4buffer/1ul20mgmlBSA/1ul10%Triton/10ulT4ligase/78ulwater16°Covernight. This physical contactligation precedes Tn5 fragmentation, unlike canonical accessibility-first HiCAR constructions; no predecessor recipe is inherited.

🟢 Tn5source orderings(main251–255/SI1Table1):

- Tn5ME-A: `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`
- Tn5MErev: `/5Phos/CTGTCTCTTATACACATCT`
- Tn5ME-B: `/5Phos/CATCGGCGTACGACTAGATGTGTATAAGAGACAG`

🟢 Source transposome assembled from unloadedDiagenodeC01070010 by vendor instructions; fullhistoricalloading stoichiometry remains unassembled. Main100ultranspositionmix=7.5ulhomemadetransposome/50ul2xtagmentbuffer/33ulDPBS/1ul10%Tween/3.5ulwater;37°C1h;100ul40mMEDTART5min;2x500ulNEBuffer3.1 wash.

🟢 Two separate orthogonal50-channel devices apply annealedA1–50 thenB1–50,37°C30min each. Sourcepreanneal10ul100uMbarcode+10ul100uMlinker1or2+20ul2xannealbuffer(20mMTris7.5–8/100mMNaCl/2mMEDTA), finalbarcode25uM. Perchannel5ul=2ulligationmix+2ulNEBuffer3.1+1ul25uMannealedbarcode;ligationmix27ulT4buffer/72.4ulwater/5.4ul5%Triton/11ulligase. Wholelinker1/2 strands and occupiedligation products remain unassembled. Images guide the two device alignments; coordinates are row/column **spatial identities**, not molecularUMIs.

🟢 SI1Table2 actual representative whole orderings:

- A1: `/5Phos/AGGCCAGAGCATTCGAACGTGATGTGGCCGATGTTTCG`
- B1: `CAAGCGTTGGCTTCTCGCATCTAACGTGATATCCACGTGCTTGAG`

🟡 Source family comparison gives the shared8nt variable identity windowAACGTGAT in these two representative orderings; A and B identify different orthogonal channel axes even when their8nt bases match. Source-knownA5P is visible; B is shown as source-written unmodifiedDNA without inventing phosphorylation. Wholearray geometry/ligation occupancy are unassembled. Applying all channel reagents is not pooling the50barcode solutions into one vessel or physically splitting2500pixels into tubes.

🟢 WholeROIlysis100ul2mgmlPK/TE8/1%SDS65°C1h onslide, collectlysate65°Covernight, ZymoChIPDNA/elute50ulEB. ActualT1beads50ul10mgmlwash2x200ul1xB&W/resuspend50ul2xB&W;50ulDNA+50ulbeads1hrotator.2x200ulTWB55°C2minwash; retainbeads/resuspend30ulEB. Select10ul of this actual stock forPCR; remaining20ul is a separate untransferredstock, not consumed by assumption.

🟢 PCRsource50ul=25ul2xQ5/2.5ul25uMcustomindexedi7/2.5ul25uMcustomP5/10ulbeads/10ulwater.The one actual25ulQ5master-mix stock exposes copresent DNApolymerase, dNTP and buffer components; these are not fresh independent additions or an asserted complete historical formulation.72°C5min/98°C30s then10–13cycles98°C10s/63°C30s/72°C1min. SI1Table1 listsN501andN701–707 whole source primers. The local graph retains the firstN501/N701 source candidates, not an asserted actualsample pairing, samplemask or primerlanding.1xAMPure/elute50ulwater, double-size-selection0.55xthen0.15x as printed; exactintermediate transfer schedule remains unassembled. IlluminaXplus2x150PE.

🟢 Main287 preprocessing splits originalR2 (genomicDNA+spatialbarcodes) computationally into newR1DNA/newR2barcodes and keeps originalR1 as newR3DNA. The same sourcecontact-library readout yields accessibility and contact analyses; there are no physically separate ATAC and Hi-C libraries invented here.2500pixels are computational barcode classes, not aliquots.

🟡 Private partial graph preserves sourcephysicalorder/stocks/wholeROI/T1selection/PCRaliquot/readclass distinction. LocalTn5duplex pairing is computed from source wholeoligos; it is not a complete transposition, contactligation or final-library proof. Wholelinker/vendorloading/occupiedDNA/enddistributions/PCRproducts remain unassembled. Actualcontrols(in-situHiC and no-spatialbarcodeHiC), microfabrication and other specimen branches remain outside selected scope; no absence conclusion.
