# APEX-RIP

Kaewsapsak et al.(2017), *Live-cell mapping of organelle-associated RNAs via proximity biotinylation combined with protein-RNA crosslinking*. DOI https://doi.org/10.7554/eLife.29224 .

🟢 Publisher eLife v1 VOR XML, MethodsPartsI–III source516–540 inspected. Supplementaryfile5 original construct/qRT-PCR primer sheets inspected; those primers are not assigned as NGS adapters. Mandatory bundle scraper run; two large data workbooks skipped by the default scan, explicitly not evidence of missing oligos. No chemistry reading of alternative archived manuscript/preprint versions. One bounded pass, no absence audit.

Source optimizedBP-first/quench/crosslink route is selected; alternative Mito and crosslink-first order experiments remain separate branches. BP500µM30min37°C; H2O2final1mM1min; azide-free ascorbate/Troloxquench1min;0.1%formaldehyde crosslink-quench1+9min; glycine125mM5min; washes/harvest. Local protein biotinylation and subsequent RNA–protein crosslinking differ physically from direct RNA biotinylation. No guarantee every recovered RNA directly binds a uniquely tagged protein.

RIPA lysis, sonication,15,000gclarification and NLBdilution precede20%input allocation and streptavidin2h4°C. Eight source washes precede ProteinaseK/sarcosine42°C1h then55°C1h. Proteolysis yields peptide/remnant-bearing RNA classes, not automatic complete reversal/removal. RNACleanXP→DNaseI37°C30min→RNACleanXP follows. Input/control routes remain separate physical libraries.

RiboZeroGold procedure is partly explicit in this source: removal mix68°C10min/25°C5min, magneticbeadbindingRT5min/50°C5min, recover RNA-containing SUPERNATANT. The graph may represent that source-established partition, unlike opaque unreviewed kits. RNACleanXP then elution into19.5µlTruSeqv2ElutePrimeFragmentmix, kit fromIncubateRFP handoff; exact internal steps/oligos remain opaque. Unique-indexlibraries pool equimolarly,≤10perpool, HiSeq2500rapidPE50plusindex. No UMI/index architecture or uninspected kit chemistry assumed.

## Bounded publication follow-up (2026-10-30)

The required `pending_library` item remains uncovered: actual historical TruSeqv2/RiboZero originals and full polymer/read assembly still require review. Published Methods516–540 supports the physical procedure, but its opaque kit handoff is not complete molecular coverage. No modern kit or alternate labeling order is substituted.

The limited fixes preserve the optimized BP-first route: labeling medium is aspirated, the azide-free quench is incubated1min and removed, then the first crosslink-quench solution is applied1min and removed before a fresh solution is applied9min. BP, oxidant, quench reagents, formaldehyde and glycine are molecular contributors distinct from conditions. This does not borrow the azide-containing crosslink-first alternative.

RNA-QA material is separated after DNase cleanup and before the selected RiboZero contribution; the source Bioanalyzer assay does not heat-cool the RNA. The32.5µl prepared-bead contribution retains an unassembled upstream90µl-slurry wash/resuspension interface, without inventing a57.5µl remaining stock. RNA remains in the actual recovered supernatant. The19.5µl ElutePrimeFragment mixture is an explicit molecular eluent after RNACleanXP retention, rather than an invented prior water elution or known fresh primer recipe. Residual cleanup material retains an unresolved fate.

DNA-HS library-QA material, the selected equimolar pool contribution, untransferred library stock and final two-lane loading contribution are now distinct. Other actual source library contributions remain external prepared boundaries; source-equimolar pooling does not consume the entire own library stock by assumption. Exact normalization volumes, index width and per-lane allocations remain unassigned. The20% input control remains separately unresolved and is not pooled as a molecular normalizer.

A literal reviewed subset replaces any implication of full historical kit reconstruction. The same one publisher-version1 original was used, without new acquisition or source expansion. Eight affected checks and a temporary private render pass; no public builder or completion marker. Notes status remains a source-scope hold, with no absence claim.
