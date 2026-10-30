# APEX-RIP

Kaewsapsak et al.(2017), *Live-cell mapping of organelle-associated RNAs via proximity biotinylation combined with protein-RNA crosslinking*. DOI https://doi.org/10.7554/eLife.29224 .

🟢 Publisher eLife v1 VOR XML, MethodsPartsI–III source516–540 inspected. Supplementaryfile5 original construct/qRT-PCR primer sheets inspected; those primers are not assigned as NGS adapters. Mandatory bundle scraper run; two large data workbooks skipped by the default scan, explicitly not evidence of missing oligos. No chemistry reading of alternative archived manuscript/preprint versions. One bounded pass, no absence audit.

Source optimizedBP-first/quench/crosslink route is selected; alternative Mito and crosslink-first order experiments remain separate branches. BP500µM30min37°C; H2O2final1mM1min; azide-free ascorbate/Troloxquench1min;0.1%formaldehyde crosslink-quench1+9min; glycine125mM5min; washes/harvest. Local protein biotinylation and subsequent RNA–protein crosslinking differ physically from direct RNA biotinylation. No guarantee every recovered RNA directly binds a uniquely tagged protein.

RIPA lysis, sonication,15,000gclarification and NLBdilution precede20%input allocation and streptavidin2h4°C. Eight source washes precede ProteinaseK/sarcosine42°C1h then55°C1h. Proteolysis yields peptide/remnant-bearing RNA classes, not automatic complete reversal/removal. RNACleanXP→DNaseI37°C30min→RNACleanXP follows. Input/control routes remain separate physical libraries.

RiboZeroGold procedure is partly explicit in this source: removal mix68°C10min/25°C5min, magneticbeadbindingRT5min/50°C5min, recover RNA-containing SUPERNATANT. The graph may represent that source-established partition, unlike opaque unreviewed kits. RNACleanXP then elution into19.5µlTruSeqv2ElutePrimeFragmentmix, kit fromIncubateRFP handoff; exact internal steps/oligos remain opaque. Unique-indexlibraries pool equimolarly,≤10perpool, HiSeq2500rapidPE50plusindex. No UMI/index architecture or uninspected kit chemistry assumed.
