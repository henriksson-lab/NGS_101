# Slic-seq

[Wei et al., Transcriptome-wide profiling of A-to-I RNA editing by Slic-seq](https://doi.org/10.1093/nar/gkad604), Nucleic Acids Research2023.

🟢 Selected sources: published main PMC10484733 Methods119, Data Processing127, Results173–181 and actualPDFp2–4; published supplementary `revised_NAR-SI_4.pdf`, originalp13 TableS1, Figure1 caption and the distinct targeted-amplicon scope of TableS2. Mandatory whole-bundle scraper pointers were verified against these originals. The source pass is closed at two originals; the nested database archive was not opened. This is a private partial preparation model, not a complete reaction or source audit.

## Selected RNA-end preparation

🟢 The selected treated HEK293T branch begins with total RNA isolated using TRIzol, oligo(dT) enrichment using Dynabeads Oligo(dT)25 and TURBO DNase treatment. Purified poly(A) RNA is fragmented with Ambion reagent70°C10min. Detailed historical kit internals and actual RNA transcripts, amounts and end distributions remain unassembled. Other cell types, tissue/input/control preparations, model-RNA tests and targeted validation libraries are separate branches.

🟢 Main119 states T4 PNK end repair, followed for treated samples by periodate and EndoV. The selected primary description does not specify an ATP addition or complete end-processing reaction; the graph does not assert PNK-generated5′ phosphate or borrow the PNK mixture from reference29. RNA3′ termini are blocked using50mM periodate25°C30min in the dark. EndoV treatment uses10mMTris7.5/0.5mMMnCl2/50mMKCl/1mMDTT/5%glycerol,37°C60min.

🟢 Results173–181 describes periodate oxidation of terminal vicinal hydroxyls to aldehyde ends, followed by EndoV cleavage at the second phosphodiester3′ of inosine, producing3′OH/5′P termini. The local qualitative end descriptor distinguishes old oxidized ends from newly generated acceptors. It does not assert quantitative oxidation/cleavage, perfect specificity or exclusive selection; unreacted/background/other end classes remain possible. Inosine remains a physical RNA residue, whereas its apparent A-to-G signal in copied reads is a separate readout. Complete atom mapping, cleavage distributions and actual occupied transcripts are unassembled.

## Published adapter disagreement

🟢 Main Methods119, actualPDFp2, prints RNA3′ adapter `/5rApp/NNNNNNTGGAATTCTCGGGTGCCAAGG/3ddC/` and DNA5′ adapter `/5Phos/NNNNNNNNNNAGATCGGAAGAGCACACGTCTG/3ddC/`.

🟢 Published SI TableS1, originalp13, instead prints RNA3′ linker `5′rAPP-NNNNNNAGATCGGAAGAGCGTCGTG-3SpC3` and DNA5′ adapter `5′Phos-NNNNNNNNNNAGATCGGAAGAGCACACGTCTG-3SpC3`. Both sources give RT primer `ACACGACGCTCTTCCGATCT`.

The model retains all four whole modified orderings as separate source-qualified components. Main RNA3′ has a different fixed arm from SI RNA3′; both main adapters have an appended physical ddC nucleotide, while both SI C3 spacers are non-nucleotide blockers. These alternatives are not merged or called an author error. Actual study usage remains unresolved. The main RNA3′ backbone is{{= chem.adapter_length('main', 'RNA3') }}nt including terminal ddC; SI RNA3′ is{{= chem.adapter_length('si', 'RNA3') }}nt. Main DNA5′ is{{= chem.adapter_length('main', 'DNA5') }}nt including terminal ddC; SI DNA5′ is{{= chem.adapter_length('si', 'DNA5') }}nt.

🟡 For its limited local binding/ligation geometry, the graph explicitly selects the **SI-qualified** adapters. The RT primer pairs to the SI19nt arm with an unpaired5′A; it does not provide that binding relation for the main small-RNA arm. The partial RNA-bound template has inferred upstream RNA and an inferred junction, retains the SI3′C3 blocker and has already released the activating AMP. Narrow source-known adapter ends do not establish whole ligated RNA/cDNA products.

```
{{scene: chem.si_rt_binding_scene()}}
```

🟢 Main119 orders RNA3′ adapter ligation,5′ deadenylase and RecJf digestion, purification, SuperScriptIII RT, DNA5′ adapter ligation, purification, qPCR-guided PCR using named NEB primers, and cleanup. The model preserves that sequence without importing ligase identities, volumes, primer doses, RT conditions, RNA hydrolysis or other operations from reference29. Activated-adapter ligation releases AMP into the same reaction mixture; byproduct fates and cleanup recovery remain unassigned. Physical ddC/C3 blockers are not copied as sugar or spacer modifications onto newly synthesized DNA.

## Identifier and reporting limits

🟢 Main127 calls the ligated barcode a “random hexamer” while printing ten Ns, and says these random barcodes distinguish PCR duplicates from different fragments with identical sequences. Main181 separately says six random positions on the RNA3′ adapter improve sequencing quality. The physical model retains source N10 on DNA5′ and N6 on RNA3′. N10 has a structured random inline-barcode feature with the source's duplicate-discrimination purpose and an explicit mapping limitation; N6 is sequencing-diversity material, not assigned a deduplication, cell or sample role. The exact processing-to-window mapping and the hexamer/N10 wording remain unresolved. The targeted-amplicon section's explicit UMI processing is not imported into this branch.

🟢 A real qPCR aliquot, of unspecified volume, is separated from production stock. Whole qPCR reagent/primer identities and PCR reagent composition and physical stock packaging are unassigned in this selected main/SI scope. Production uses the remaining cDNA and source-named NEB primers, whose full bodies, index identities and acquired landing remain unassembled. The chosen final cleanup is NEBNext Sample Purification Beads; the source's low-melting-point agarose alternative remains a separate option. Purification retains uncertain molecular recovery/carryover rather than perfectly deleting background.

🟢 Main127 reports HiSeqXTen PE150. An actual loading allocation is separated from retained stock. Trimming, alignment, duplicate removal, background filtering and inosine calling are computational operations; they do not physically purify the library. Exact run primers/index cycles and complete molecular read geometry remain unassembled.

## Remaining scope

🟡 Whole RNA/cDNA sequences, end heterogeneity, both adapter variants' actual experimental use, complete ligation/PCR products, primer/index identities, processing-window interpretation and reagent/kit details remain open. Reference29 is explicitly cited as a modified predecessor but its recipe is not inherited. Vendor/predecessor originals remain ordinary unread targets, not confirmed missing disclosures. No absence claim or second absence-review pass is accepted.

The defining article's bibliography identifies reference29 as Zhou et al., *Evolution of a reverse transcriptase to map N1-methyladenosine in human messenger RNA*, [10.1038/s41592-019-0550-4](https://doi.org/10.1038/s41592-019-0550-4). This is a source-qualified library-preparation attribution, not evidence that Slic-seq inherits that paper's modified reverse transcriptase, reaction conditions or complete chemistry. That original was not opened.
