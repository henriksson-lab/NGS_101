# Mod-seq

Talkish et al. (2014), *Mod-seq: high-throughput sequencing for chemical probing of RNA structure*, DOI https://doi.org/10.1261/rna.042218.113 .

## Sources and official correction

🟢 Published archived PMC3988572 main, overview/Fig1 and Methods98–120 inspected. The actual downloaded EuropePMC supplementary ZIP contains only figure aliases (three figures and an equation), not the named chemistry/software Supplemental Files. Its inventory is not an inspected chemistry supplement. Mandatory scraper source pointers and modification strings were checked in main108–116. No sequence/reagent absence claim.

🟢 Official corrigendum (2019), DOI https://doi.org/10.1261/rna.072405.119 , primary publisher indexed full correction paragraph inspected at https://rnajournal.cshlp.org/content/25/9/1226_1.full . The archived original Methods116 prints reverse PCR arm `TGTACTGGAGTTCAGACGTGTGCTCTTCCG`; the correction explicitly replaces it with `GTGACTGGAGTTCAGACGTGTGCTCTTCCG`. The local model uses the corrected arm and preserves both written values. The correction says the online article was corrected; the cached PMC original still contains the earlier value, not a claim that today's publisher original still does. Local full correction retrieval failed (EuropePMC XML HTTP500, PMC HTML no articlecontent, DOI noopenfulltext); primary indexed correction reading is qualified rather than described as inspected original PDF/JATS. The source guard explicitly exempts only that corrected primer from matching the old local text. No private preview is generated; any future affected page must display this correction.

## Preparation and physical fates

Source branches are DMS-treated yeast, no-DMS control and STOP control (DMS added **after** quenching). DMS100mM30°C2min is followed by ice/β-mercaptoethanol/isoamyl-alcohol quench, cell wash, phenol extraction and ethanol RNA recovery. The source discusses RT-blocking m1A/m3C and weaker G-associated signals; the model does not make every DMS adduct a guaranteed stop or turn a reported DNA/reference base into a physical RNA substitution.

Source3µgRNA fragmentation uses Ambion reagents70°C3min, Zymo>17nt cleanup, PNK phosphatase37°C1h followed by separate ATP1mM/PNK5′ phosphorylation37°C1h and column cleanup. Full fragmented-end population is unassembled.3′ linker `rAppCTGTAGGCACCATCAAT–NH2` is adenylated and3′ blocked; NH2 is not a nucleotide. Rnl2truncated16°C overnight ligation precedes Zymo>200nt cleanup. The5′ End Acceptor17.93R is `ATCGTaggcaccugaaa`: uppercase DNA and lowercase RNA are distinct sugar blocks, not a homogeneous DNA adapter. Rnl1/QuickLigase buffer22°C overnight and a second Zymo>200nt cleanup follow.

The RNA capture probe `TTTCAGGTGCCTACGAT` carries **3′ Biotin-TEG** and hybridizes to the mixed-sugar5′ linker. Source90°C90s denaturation,37°C→RT annealing3°C/min and MyOneC1 capture precede80°C2min TE elution/immediate magnet recovery of the whole supernatant. Separated flowthrough/beads remain explicit unresolved companions; ideal specificity/release is not claimed.

SSII50°C RT primer is5′ phosphorylated and contains two `(SpC18)` non-nucleotide nodes: `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC-(SpC18)-CACTCA-(SpC18)-TTCAGACGTGTGCTCTTCCGATCTATTGATGGTGCCTACAG`. The actual RT mixture retains independent original RNA/unextended primer/reagent material before NaOH100mM98°C20min. Hydrolysis does not itself physically remove products; ethanol and8%ureaPAGE150–330bp recovery follow. CircLigase60°C1.5h closes selected cDNA while retaining both C18 nodes. The constructor refuses to flatten spacers into ordinary DNA bases. Circle-to-PCR polymerase traversal remains unassembled, not borrowed from similar irCLIP/LASER primers.

Subtractive cDNA probe `ATCGTAGGCACCTGAAA` has **5′ Biotin-TEG**, the opposite orientation to the RNA capture probe. MyOneC1 removes the bound full-runoff-linker cDNA fraction; remaining cDNA is ethanol recovered. This physical subtraction differs from later computational5′-adapter read filtering and does not guarantee every remaining molecule is a true DMS stop. Phusion/HF12cycles64°C, MinElute, TapeStation, pooling and HiSeq2000 follow. Corrected reverse PCR family is `CAAGCAGAAGACGGCATACGAGAT[index6]GTGACTGGAGTTCAGACGTGTGCTCTTCCG`; index meaning is structured, study member unassigned. Actual run cycles, sequencing primer landings and complete junction geometry remain unassembled.

## Component ancestry and limits

Main112 explicitly attributes the RT-primer component to Brar et al.2012; main108 attributes the5′ linker to Lau et al.2001. These are specific source-owned component leads, not evidence for importing a complete predecessor recipe. Direct originals and historical vendor guides remain pending. The linked genuine chemistry supplement is unreviewed; no two-reader absence audit, finished schematic or full-publication claim.
