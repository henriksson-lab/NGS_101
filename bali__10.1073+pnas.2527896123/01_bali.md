# BALI

Battistoni et al. (2026), **Spatially tunable multiomic sequencing using light-driven combinatorial barcoding of molecules in tissues**, [PNAS](https://doi.org/10.1073/pnas.2527896123), published18May2026. The2024preprint is superseded and not used as chemistry authority.

🟢 Read: published mainPMC13214022 root/index mechanism129–141 and published SI Appendix extendedMethods321–389, selected joint450–555, TableS2originalp32. Two originals, source investigation closed. Mandatory scraper46candidates inoneSI text; selected literal oligos and modifications checked in original table layout. DatasetS1 and other original components remain outside this pass, not absent.

🟢 BALI writes spatial codes through light-activated ligation of staggered double-stranded DNA indices. A root has a photocaged5′phosphate and6ntsingle-stranded overhang; subsequent index modules regenerate orthogonal6/7ntcaged overhangs for serial writing. Final modules carry a UMI and PCR handle. Published examples include RNA, chromatin accessibility and joint RNA/ATAC; these are one method with different branches. The selected graph selects the two-region, single-digit **RNA demonstration**. It does not claim to model the full multiomic branch or generalize the exact RNA recipe to all barcode lengths.

🟢 Selected RNA route: fresh-frozen E16.5coronal brain sections equilibrate5minRT/5min37°C, hydrophobic boundary dries5min;500µlPBS removesOCT. Fix500µl.5%PFA **exactly5minRT**, quench1.25Mglycine quick+5min wash, PBS twice5min. .5%Triton15minRT; PBS3×5min; .1NHCl5min; PBS2quick+5min; prewarmedPBS65°C5min thenicePBS10min. Subsequent preparation is light protected.

🟢 RTfirst14µlmix includes2.5µldNTP10mMeach/1.5µl100µMBL621/water,65°C5min→ice>2min **before enzyme addition**. Add36µlmix with10µl5×SSIVbuffer/2.5µlSSIV/1µlRNAsin/3µl100µMBL617/5µl5Mbetaine/.5µlETSSB/2.5µl100mMDTT/water; total50µl under22×22mmcoverslip,42°C~16h. This combined RT/template-switch transition is deliberately opaque until complete RNA/cDNA and switched-junction geometry is assembled. Whole source orderings5′→3′:

```
BL621: TAA/iSpPC/AATAACC CTAG NNNNNNNN ACTG TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN
BL617: tacactctttccctacacgacgctcttccgatct rGrGrG
BL599: GGTTATT CG ATCTAGGT GC
BL623: /5deSBioTEG/AGACGTGTGCTCTTCCGATCT NNNNNN GC ACCTAGAT CG
BL601: GGTTATT CG GGACTCGT GC
BL624: /5deSBioTEG/AGACGTGTGCTCTTCCGATCT NNNNNN GC ACGAGTCC CG
```

BL621 has30T and terminalVN; its8Nwindow remains a source variable window with unassigned functional identity in this scope, distinct from the terminal-index6NUMI. `/iSpPC/` occurs aftertheinitialTAA; no PC modification is replaced with an ordinary base. `/5deSBioTEG/` is desthiobiotinTEG, not an ordinary biotin annotation. The TSO includes the literal leadingT before the canonical TruSeqRead1body and three terminal RNA guanosines; no SMART-seq2 TSO is substituted.

🟢 PBSwash/2×SSCtwice5min; Draq5 1:1000in2×SSC10min thenwash twice5min. Index pairsBL599/BL623(SVZ) andBL601/BL624(cortex) are annealed45µMin2×SSC at95°C5min→65°C5min→4°C, ramps.1°C/sec. Each region receives405nmSP5confocal30mWmaximum-power5min illumination in200µl2×SSC. Washquick+5min2×SSC, equilibrate1×ligasebuffer5min. Ligate50µl with3µMannealedindex/T4ligase **printed100U/µl** in50mMTris7.4/10mMMgCl2/10mMDTT/1mMATP,25°C30min. Each ligation is followed by2×SSC+Draq510min,2×SSC+10%EC5min,2×SSC5min. These are sequential region exposures in the same section, not physical removal of all other regions. The full paired index strands are assembled as a checkedScene: the8ntcode and two flankingGC/CG segments form a12bpstem, with a source7nt5′GGTTATT cohesive overhang. The bottom5′desthiobiotin/PCRarm/6NUMItail remains on its source strand. Actual preannealed index stocks contribute selected ligation aliquots; untransferred stocks remain unresolved. Adapter-end geometry is declared known; full RNA/cDNA substrate and sealed product remain unassembled. The generic6ntroot description does not override these source RNA7ntinterfaces, which fall within the published6/7ntdesigns.

🟢 Collect200µlRIPA-basedlysate60°C~16h; RNA paragraph prints.8U/mlwithout assigning enzyme identity there. The joint branch later explicitly discusses proteinaseK inactivation, but no unreviewed enzyme specification is silently filled into the RNA paragraph. The RNA lysis boundary is therefore opaque with unassigned active-stock type; it is not typed as proteolysis. The known SSIV RT reaction retains its explicit catalyst role. MinElute/elute25µl. A2µl analytical qPCR allocation withQ5/universal/A012/EvaGreen determines the cycle count atone-thirdplateau; it is not returned to production. Use10µl for productionPCR, retaining13µl outside that reaction. Production uses sourceA0XX, not an invented chosen sample code. Both primer patterns preserve the exact published body:

```
universal: AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
A012: CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
A0XX: CAAGCAGAAGACGGCATACGAGAT[6nt index]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
```

The right primer body is one base shorter than full TruSeqRead2. SourcePCR98°C3min then98°C15s/68°C20s/72°C20s; actualcycles fromqPCR. .8×AMPure/elute20µlEB; separateTapestationD5000/KAPAlibraryQC; pool andNovaSeq. Main141 specifies100bpR1/200bpR2 in300cyclereagents; exact index mask/primer landings remain unassembled.

🟢 Alternative joint preparation, **not graphed**: adult brain .5%PFA/permeabilization→Tn5transposition35min37°C→EDTAwashes→phosphorylatedBL730RT/TSovernight. Three bulk index ligations followed by final region-specific modules use365nmbulkuncaging and405nmROIuncaging; firstligation45min, others30min. From200µllysate only150µl is PMSF-treated before streptavidin capture. Beads contain biotin-Tn5-enrichedDNA, supernatant cDNA-enrichedRNA; each has separate library finishing. RNA additionally undergoesRNaseH, preamplification andsingle-adapterTn5; DNA useson-beadgap-fill/PCR. No whole joint-library constructor is inherited from SHARE-seq or DBiT merely by conceptual similarity.

Open: complete switched/ligated products and read geometry; actual production sample code/cycle count; precise RNA-lysis catalyst assignment; whole joint route and all relevant linked component originals. No absence conclusion or complete chemistry is claimed.

## Public preparation scope

The public schematic covers the selected published preparation through its physical sequencing contribution. Its source inventory is a literal curator-written set of reviewed inputs and operations; adding graph targets cannot automatically extend review. Existing scoped source findings remain required and mapped locally, without exclusions or deletion. Full preparation, stock/analytical fractions and chemistry warnings remain visible in a collapsed workflow.

RNAparagraph0.8Uml/RIPAdoesnotnameactivecomponent; nojointproteaseinheritance. Generic6ntrootversusactualRNA7ntindexscopeexplained,notpublicationerrororrepair. Full occupied molecular products, primer landing and unspecified read masks remain unassembled rather than silently supplied from another method or current kit.

BALI source RT and Q5 supplies now expose their co-present enzyme, nucleotide, inhibitor/companion, dye and buffer functions. The actual 25 µl 2× Q5 master mix is distinct from existing primer inputs. These are views of the same actual source additions, not fresh doses or a proven aggregate premix. Both used and remaining annealed index portions retain the known source duplex geometry.
