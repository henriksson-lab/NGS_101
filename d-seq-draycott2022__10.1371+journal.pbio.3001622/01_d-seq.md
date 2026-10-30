# D-seq — reduced-dihydrouridine RT-stop libraries

Draycott et al. *Transcriptome-wide mapping reveals a diverse dihydrouridine landscape including mRNA*. PLOS Biology20:e3001622(2022). [Published article](https://doi.org/10.1371/journal.pbio.3001622).

## Source scope

Read the published main PMC9129914 Methods155/159–164, Results108 and actual Figure1 originalPDFp3, plus main supplemental-file identities. One original; PDF/XML are representations of that same article. Source investigation closed. Separate supplement originals, directly cited predecessors and vendor specifications remain unreviewed. Mandatory scraper ran; its18hits concern syntheticU/D RNAs and knockout-genotyping primers, which are not operational library adapters. Its default10MB limit skipped the32MB raw-image PDF; main identifies that file as uncropped Top7 SDS-PAGE/mRNA-gel data, but its original was not inspected. This is a positive capped preparation pass, not a full source or absence audit.

D-seq reduces D to tetrahydrouridine and selects RT-stop cDNA. It does not inherit Rho-seq's dye-derivatized substrate or recipe. WT/knockout conditions are controls, not separate protocol identities. The model covers one WT D-seq preparation; parallel RNA-seq controls omit NaBH4 and select full-length90–100nt cDNA, a separate source branch not merged here.

## Preparation

🟢 Source main155 describes two successive oligo-dT-cellulose poly(A) selections. The selected graph starts from that prepared mRNA stock; total-RNA extraction, culture/genotyping and initial8mg totalRNA stock are outside this graph. The >200nt cleanup in that paragraph applies to LC–MS analysis, not the modeled D-seq library.

🟢 Main159: fragment in10mMZnCl2,94°C1min, then precipitate. Add2µl10mg/ml NaBH4 in500mMTris pH7.5 stock to18µl RNA/water,0°C1h;500mM describes the added stock, not the final20µl reaction. Quench with4µl6N acetic acid and precipitate. D reduction and RT stopping are not assumed perfect or absolutely specific; source108 notes other modified nucleosides can also react.

🟢 CIP/PNK end processing precedes RNA size selection on8%urea-TBE PAGE. Their exact order, ATP conditions and full resulting end populations remain unassembled; two fresh serial additions and universal5′P/3′OH products are not invented. RNA is eluted overnight4°C in400µl300mMNaOAc pH5.5/1mMEDTA/100UmlRNasin and precipitated.

**Published-source wording difference:** Fig1D labels RNA size selection100–120, while Methods159 specifies RNA fragments70–80nt. The model uses explicit Methods70–80nt. Whether the figure uses a different counting convention or represents another window is unresolved; no correction is asserted. Both source locations remain recorded in the model notice. Truncated cDNA50–80nt and control full-length90–100nt are distinct later selections.

🟢 Ligate a preadenylated IDT RNA3′adapter with sourceT4RNA ligase in ATP-free50mMTris7.8/10mMMgCl2/10mMDTT,22°C2.5h, then precipitate. The local end view shows source preactivation and an explicitly inferred ligation-compatible RNA3′OH class. Whole adapter body/blocking, actual complete substrate and sealed junction remain unassembled. First-ligation ligase isoform is unassigned and not borrowed from the second reaction.

🟢 RNA/RTprimer95°C5min then ice5min; add1µl10mMdNTPs,1µl0.1MDTT,1µlRNasin,3µl5×FirstStrandBuffer,1µlSSIII;50°C1h. The same annealed primer component enters RT without a second stock addition. Source RT stops are observed1nt3′ of D; the model retains stopped/full-length products and does not turn every RT stop into a specific D call or physically cleave RNA there.

🟢 Select50–80nt cDNA on8%urea-TBE PAGE; source400µl300mMNaCl/10mMTris pH8 overnight room-temperature gel elution/precipitation workup is retained as a compound physical operation. Transfer5µl actual eluate for the next reaction and retain the untransferred stock; no perfect400µl recovery or exact395µl remainder is claimed.

🟢 Mix0.8µl80µM IDT linker,1µlDMSO and5µl cDNA;75°C2min then ice2min. Add2µlRNA-ligase buffer,0.2µl0.1MATP,6.5µlPEG8000,3.6µlwater and0.5µlT4RNA ligase2,overnight22°C. The literal “5′adapter” name is preserved without automatically assigning it to the chemical cDNA5′end. Whole linker/end/junction geometry and exact ligase variant remain unassembled. MyOneSilane cleanup precedes Phusion PCR; primer bodies/cycles/index identities remain unassigned. Gel-purify, precipitate and physically pool PCR products, then sequence on HiSeq2500. Run/index cycles and whole sequencing-primer geometry remain unassembled.

## Remaining scope

Main source identifies S1–S3 supplementary figures and S1tRNA/S2snRNA/S3mRNA data tables; these originals are available but unread, not access failures or proven irrelevant. Directly cited validation/chemistry methods and applicable historical vendor specifications require later depth-one review. Whole library oligos, physical RNA/DNA end populations, linker polarity/junctions and complete read geometry remain open model questions, not fabricated unavailable disclosure documents. No missing-information assertion or full-reaction validation is accepted.
