# Circle-damage-seq

The major mechanism of melanoma mutations is based on deamination of cytosine in pyrimidine dimers as determined by circle damage sequencing. Science Advances (2021), [doi:10.1126/sciadv.abi6508](https://doi.org/10.1126/sciadv.abi6508).

Source scope: 🟢 published main article PMC8324051, Methods text139–185, including genomic CPD, deaminated-CPD and separate E. coli DpnI pilot; adapter spelling checked in original PDF p11. One original read. Mandatory scraper reviewed (7 candidates from PDF/XML representations of the same article). The EuropePMC supplementary ZIP contains main-article figure JPG/GIF files, not the cited supplementary text. No inference of supplementary absence follows. The private graph models the selected genomic CPD route; the two other routes remain explicitly separate. No complete linked-source or absence audit.

🟢 Genomic DNA is circularized before lesion-directed opening. This is not ordinary linear CPD mapping and not a general inverse-PCR library. Fibroblast CPD libraries use UVB1000J/m² and immediate extraction; deaminated-CPD libraries use4000J/m² and24or48h at37°C. QuickDNA MiniprepPlus supplies isolated DNA. CovarisE220 generates300–400bp fragments at4°C; DCC5 cleanup elutes42µl. A50µl repair reaction contains1µgDNA, RNaseH2U, T4DNApolymerase3U, T4PNK10U and dNTPs,24°C30min/70°C20min; DCC5 elutes100µl and NanoDrop provides separate QC. Blunt1µgDNA is circularized using1600UT4ligase in200µl16°Covernight. After DCC5 cleanup, PlasmidSafe20U digests linear DNA37°C30min/70°C20min, then1.8×AMPure elutes35µl. Digestion products remain until cleanup; circularity or recovery is not assumed perfect.

🟢 CPD route: T4PDG10U+APE110U in40µlNEBuffer4/BSA37°C20min;1.8×cleanup/elute25µl; E.coliphrBphotolyase0.75µg in50µl under365nmUVA at15cm for90minRT;1.8×cleanup/elute32µl; S1nuclease5U in40µl for4minRT. EDTA/Tris quenches,70°C10min follows; source states1.8×AMPure72µl/elute48µl. Full lesion-specific cuts and off-products remain unassembled rather than assigning genomic coordinates.

🟢 Deaminated-CPD route instead applies photolyase first, cleans/elutes35µl, then UDG5U+APE110U in40µl37°C20min, cleans/elutes32µl and applies S1. It does not inherit T4PDG from the CPD branch. The separate E.coli pilot uses RsaI/AluI fragmentation and a USER-cleavable hairpin before circularization/DpnI opening and NextSeq500PE75; it is not incorporated into the genomic graph.

🟢 After S1: Klenowexo-minus10U/dATP0.2mM in55µl37°C30min/70°C30min. UltraIILigationMasterMix30µl+Enhancer1µl+1.5µMadapter3µl is added directly and incubated20°C15min; USER3U37°C20min follows. Source states93µl ligation, although listed55+30+1+3=89µl. Preserve this apparent volume discrepancy, with4µl unassigned. Cleanup is1×AMPure96µl followed by0.5×300–1000bp selection/elute25µl.

🟢 Adapter orderings are literal, including both terminal variants:

```
5'-phos-GATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCTACACGACGCTCTTCCGATC*T-3'
5'-phos-GATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCTACACGACGCTCTTCCGATC*C-3'
```

The5′ phosphate, internal deoxyuridine and final phosphorothioate bond are retained. The printed right arm has `ACACTCTTTCCT`, differing from a canonical `ACACTCTTTCCCT` sequence; the original PDF confirms this printing. Neither replacing it with the canonical adapter nor declaring a publication error is justified by the selected scope. Complete folding, ligation and USER product geometry remains unassembled.

🟢 PCR lists25µllibrary+25µlQ5+2.5µli5+2.5µli7 but calls this50µl; listed sum55µl is a second apparent source volume discrepancy. Initial98°C30s,12cycles98°C10s/65°C75s,65°C5min;1×AMPure50µl/elute20µl. Qubit3/Bioanalyzer/KAPAlibraryquantification are distinct analytical QC allocations. HiSeq2500PE150 is acquisition; later2×80 hard clipping/divergent alignments are computational. Their selected cut interpretation is not an assertion of universal perfect lesion specificity.

Open scope: published supplementary original, historical NEB adapter/primer documents, i5/i7 identities, whole circular/opened-product geometry and other route implementations. These are unfinished checks, not confirmed missing disclosure. No curated ancestry established in this bounded pass.
