# S1-END-seq

Matos-Rodrigues et al. (2022), *S1-END-seq Reveals DNA Secondary Structures In Human Cells*, [defining article](https://doi.org/10.1016/j.molcel.2022.08.007).

## Selected evidence

🟢 Associated author manuscriptPMC9547894 mainMethods203/209, Fig1caption636 and KeyResources821/823 read. One original; mandatory scraper run. The exact associated SupplementaryFiguresPDF route is gated; it has not been read. Published final/AM recipe equivalence is unverified. No superseded preprint, inferred final replacement or upstream diagram is authority.

## Preparation

🟢 Source10–20millioncells embedded in1%agaroseplugs. For S1 branch, ProteinaseK1h50°C then7h37°C; oneWB/threeTE washes, RNaseA1h37°C; threeWB, EB, twoS1buffer equilibration washes, each15min. S1buffer40mMsodiumacetatepH4.5/300mMNaCl/2mMZnSO4. Add200UThermoFisherS1 in100µlbuffer perplug,37°C30min; quench10mMEDTA. This source explicitly distinguishes its200U from1.8Uusedfor neuronalSSB detection; do not merge those applications or infer equal cleavage behavior.

🟢 Source then delegates standardEND-seq. Fig1caption outlines biotinylatedadapterligation toDSBends, sonication, streptavidincapture and addition ofIlluminasequencingadapters. The historical complete END recipe remains unassembled here; known outline is not permission to import a predecessor constructor wholesale. Cleaved/unreacted structural classes remain mixed until actual workup; S1 sensitivity is not a guaranteed structure-specific conversion.

🟢 KeyResources explicitly lists both END hairpin whole-orderings with5′phosphate, dU, biotindT and terminalphosphorothioate. Inventory is retained verbatim; U is chemically dU, not silently replacedbyT/N in molecule inventory. Without inspecting linked END detailed originals, this inventory alone does not validate wholeproducts/USERcleavage/readlanding.

🟢 DefaultNextSeq550SE75; source also says150bp paired-end kits in'MiSeq200' (literal source label retained, not instrument correction). Treatment/control libraries run inparallel sameflowcell; exact loading quantities/pool composition unassigned. Neither independent PacBio WGS norcellcycleEdU/DAPI staining is incorporated into the S1library chemistry.

## Scope and open work

🟡 Notes/privatepartialgraph establish the source-specific200US1frontend and actual known outline; linkedEND operations, whole products, PCR/index/read geometry and originalSI remain unassembled. No publicationabsence or final/AM equivalence claim. P1 and baselineEND are same-paper branches, not automatically separate collected methods. No fabricated source document represents an unassigned value.

Computedwholemodified5′→3′orderings (dUretainedchemically; starsarePSbonds):

- Adapter1: {{= chem.RAW1 }}
- Adapter2: {{= chem.RAW2 }}
