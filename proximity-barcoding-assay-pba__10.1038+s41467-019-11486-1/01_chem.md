# Proximity barcoding assay: selected exosome preparation

Wu et al. (2019), **Profiling surface proteins on individual exosomes using a proximity barcoding assay**, DOI10.1038/s41467-019-11486-1.

## Sources and scope

🟢 Published PMC6710248, Fig1/Methods163–187 and tag-meaning paragraphs; published SI MOESM1 Note1/Table1, selected whole orderings, with original p11 checked for modified blocker and row assignments. Two originals; investigation closed. Mandatory scraper ran. No complete linked-source or absence audit. The private incomplete graph selects exosome PBA rather than the separate streptavidin validation or SP-PLA experiments.

## Distinct preparation and identifiers

🟢 Each antibody DNA carries an eight-base proteinTag and an eight-base random moleculeTag (UMI). A separate circular padlock contains a fifteen-base random complexTag, amplified into uracil-containing RCA bundles. Antibody-labelled exosomes bind CTB-coated wells; nearby reporters hybridize to one RCA bundle and copy its complexTag by T4-polymerase extension. The model does not import MPX, add a second pixel round, or infer that barcode collision/neighbor capture is impossible.

## Source orderings and chemical marks

🟢 Ab1: `/5ThioMC6-D/AAGCGTGTGATTCCTAGTAATGAACAATGGNNNNNNNNTGCTATTATGATGTCTCAGGT`.
Padlock: `5′/Phos/GCTATTATGATGTCTCAGGUNNNNNNNNNNNNNNNAGAUCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCTAGTGCTGGATGATCGTCC`.
Ligation template/RCA primer: `CTGAGACATCATAATAGCGGACGATCATCCAGCACT`.
Blocking oligo, original p11: `CTAGTGCTGGATGATCGTCCmUmUmUmUmU` (five 2′-O-methyl RNA bases).

🟢 Literal padlock U characters and phosphate are preserved; the terminal blocker mU modifications are not silently changed to DNA dU. The helper retains one complete Ab1 ordering and its UMI, but the fixed proteinTag boundary is not assigned locally. The full forty-conjugate pool and antibody-specific proteinTag mapping remain unassembled. Main PCR names pba-fwd/rev, fwd-index/rev-index and library-fwd/rev are not automatically equated to the table's differently named Universal-primer/D-index rows. Table1 starred orderings are explicitly for the streptavidin validation and are not used for exosome PBA.

## Order and physical provenance

🟢 Main163: SMCC antibody activation and separate thiol reduction/desalting precede overnight conjugation/dialysis. The graph starts from the resulting preconjugated pool; this upstream synthesis is described but not reconstructed.

🟢 Main167: padlock/template/T4 ligation37°C30min → ExoI/III37°C30min →85°C20min → RCA-primer adjustment/reannealing37°C20min → phi29 with d(A,U,G,C)TP37°C10min →65°C10min → store−20°C. Whole circular closure, repeat product and enzyme occupancy remain unassembled. Primer adjustment is explicit; the produced co-retained primer component enters RCA without an accidental duplicate use of the original template stock. Source adjustment remains a distinct addition of unassigned amount. Source splinted single-strand closure is typed as known; the complete paired/covalent product state remains unassembled.

🟢 Main175–187: overnight antibody/exosome binding → selected1µl dilution →20µl capture → washes → RCA/blocker37°C15min → washes → T4 extensionRT15min → wash → one Phusion/UNG reaction, first37°C15min then95°C/PCR →1µl allocation to two-cycle indexing → twenty-fold dilution into fresh PCR → fifteen-cycle library amplification → pool/purify → sequence. The model retains untransferred RCA, staining, dilution, amplification and indexed-library stocks; source-unspecified amounts stay unassigned. Washes are physical fates, not perfect-selectivity assertions. Other sample members of the final pool remain unassembled.

🟢 The first UNG/Phusion reaction is one bulk addition, represented by produced same-stock catalyst, primer and buffer components. Index and final-library stages consume co-transferred prior-enzyme components of the selected reaction aliquots; amounts and activity remain unassigned. Main187 names Phusion HF buffer in the final fresh mix, so the model does not infer a fresh polymerase addition from that buffer name. New source primer/buffer mixes likewise contribute produced same-stock components, rather than double-counting stock and free primers. These unresolved PCR boundaries use compound operations, not checked amplification products. Main183 and SI Note1 describe uracil-template inactivation/cleavage and T4 exonuclease suppression of unextended reporters. These are not labelled as physical removal before an actual workup.

## Endpoint and open questions

🟢 Source names MiSeq v2 300 cycles or NextSeq75SE. Exact per-run masks, full circle/extension/product/primer geometry and the mapping of source PCR names remain unassigned locally. Indexing and physical reagent transfers are source-positive; complete NGS readiness is not claimed. The source's singleton/read-count filtering occurs computationally, never as physical removal of DNA. No source-information absence is claimed; linked component originals and complete antibody/sequence mapping remain unreviewed.

## Source-qualified public preparation page

Source-qualified public page adds literal inventory with original selected SourceItem retained and covered. RCA ATP/nucleotides and extension dNTP are actual molecular contributions. Same PCR bulk stock exposes physical nucleotide/dye/DMSO inventory without a new enzyme addition. Actual source PCR pool consumes the selected contribution plus other prepared products, not cloned copies of one product. Prepared exosome/AOC input boundaries remain visible. ProteinTag/sample-index meanings are structured unpositioned metadata; index members/whole primer coordinates remain unresolved. Known complexTag/moleculeTag and five2primeO-methylU blocker remain. No STARRED-STV-only recipe or MPX secondcycle imported.

The public page includes a preparation overview, known molecular components, read scope and a collapsed complete workflow with unresolved-chemistry warnings. Full occupied molecules are not supplied as invented biological sequences. This page is a capped source-supported preparation declaration, not a full linked-original or missing-information audit.
