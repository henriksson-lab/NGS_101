# GLORI

**Defining source:** Liu *et al.*, *Nature Biotechnology* (2023),
doi:10.1038/s41587-022-01487-9. **Detailed protocol:** Shen *et al.*, *Nature
Protocols* (2024), doi:10.1038/s41596-023-00937-1. Both main texts were hand-downloaded by the
maintainer on 2026-10-17 (`s41587-022-01487-9.pdf`, `s41596-023-00937-1.pdf`; printed DOIs,
titles and authors checked) and read in round deferred-followup-03 (section below). Also read: the 2023 Supplementary Figures
(MOESM1), the 2023 Reporting Summary (MOESM2, image PDF, rendered p. 1), Supplementary
Data 8 (MOESM9) and the 2024 Supplementary Tables 1–2 (MOESM3 XLSX, in
`$CHEM_DATA/sources/doi_10.1038+s41596-023-00937-1/`).

## Evidence used

- 🟢 Glyoxal and nitrite catalysed deamination of A, with G protected by glyoxal and
  then deprotected (2023 Supplementary Fig. 1 legend, MOESM1 lines 76–90).
- 🟢 Liu 2023 Methods, "RNA protection, deamination and deprotection" (`.pdf.txt` lines
  2612–2633): fragmentation at 94 °C for 4 min (NEB E6150S); 200 ng RNA plus spike-ins with
  6 µl 8.8 M glyoxal and 20 µl DMSO, 50 °C 30 min; 10 µl saturated H3BO3, 50 °C 30 min;
  deamination in 750 mM NaNO2, 40 mM MES pH 6.0 and glyoxal for 8 h at 16 °C; ethanol
  precipitation; deprotection buffer (triethylammonium acetate/formamide) at 95 °C for
  10 min. (Model correction: an earlier note's additional "4 °C overnight" incubation is not
  in the 2023 Methods, nor in 2024 Step 23, which specifies 8 h at 16 °C.)
- 🟢 Liu 2023 Methods, "Library preparation and sequencing" (lines 2635–2658), eCLIP branch:
  T4 PNK (M0201S) dephosphorylates 3′ ends; 3′ linker `5′rAPP-AGATCGGAAGAGCGTCGTG-3SpC3`
  ligated with T4 RNA ligase 2 truncated KQ (M0373S); excess linker removed with 5′
  deadenylase (M0331S) and RecJf; RT with primer `ACACGACGCTCTTCCGATCT` and SuperScript III;
  excess primer digested by exonuclease I; silane beads; 5′ adapter
  `5′Phos-NNNNNNNNNNAGATCGGAAGAGCACACGTCTG-3SpC3` ligated on beads with T4 RNA ligase 1 high
  concentration (M0437M); PCR "with P5 and Illumina-indexed P7 primers"; 1.8× XP beads;
  165–250 bp from an 8% TBE gel; "Illumina Hiseq X-ten with paired-end (PE) run and read
  lengths of 2 × 150 base pairs". (The model's step "repair" was a phosphorylation; it is now
  the 3′ dephosphorylation the Methods describe.)
- 🟢 Liu 2023 Methods, data filtering (lines 2688–2696): "reads R2 yields sequences sense
  to the original RNA. Thus, only reads R2 was used"; duplicates removed on the "10-base-pair
  UMI at the 5′ end of reads R2", then the UMI is trimmed (`fastx_trimmer -f 11`).
- 🟢 Shen 2024 Nature Protocols: overview lines 360–391 (eCLIP Steps 36–122, alkaline
  hydrolysis of the RNA template after RT, HiSeq X Ten "paired-end 150 bp reads"); Reagents
  lines 1871 and 1891–1898 (NEBNext Index Primers Set 1 E7335S; the 3′ RNA adaptor is
  5′-adenylated with a 3′ three-carbon chain; the 3′ DNA adaptor is 5′-phosphorylated with a
  3′ C3 spacer and a 5′ N10); Steps 92–93 (lines 3070–3200): NEBNext Index Primer and
  NEBNext Universal PCR Primer for Illumina with NEBNext Ultra II Q5, 10–12 cycles (cycle
  table by input). The 2024 alkaline hydrolysis step is not in the 2023 Methods and is not
  drawn.
- 🟢 Modified eCLIP libraries were compared with NEB #E7300 and Takara #634487 kit
  libraries (2023 Supplementary Fig. 14 legend, MOESM1 lines 327–357).
- 🟢 The 2023 Reporting Summary (Data collection) states: "Next-generation sequencing data
  was collected and demultiplexed by Illumina HiSeq X-Ten and Nova." Paired-end mode
  and cycle counts are not stated there; the 2023 Methods give 2 × 150 paired-end (above).
- 🟢 The defining article's Supplementary Data 8 (`41587_2022_1487_MOESM9_ESM.xlsx`)
  publishes sequences for synthetic m6A spike-in **RNA templates** and other assay
  controls. These are not the adapters or PCR primers that form the sequencing library.
- 🟢 The detailed protocol's [Supplementary Tables 1–3](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41596-023-00937-1/MediaObjects/41596_2023_937_MOESM3_ESM.xlsx) give the **3′ RNA adapter, RT primer and 3′ cDNA adapter with a 10-nt UMI**. Table 3's RT/PCR/sequencing primers are specifically for conversion evaluation of synthetic Probe 5, not the complete RNA-seq library.
- 🟡 Derived from Supplementary Table 1 (computed with `lib/illumina`): the 3′ RNA adaptor
  is the reverse complement of the TruSeq Read 1 primer 3′ end and the RT primer is that
  3′ end; the 3′ DNA adaptor arm is the reverse complement of the TruSeq Read 2 primer 3′
  end. The cDNA strand is therefore Read 1 site – antisense insert – N10 UMI – Read 2
  site. Under the standard P5/Read 1 convention, Read 2 reads the UMI first and then the
  RNA-sense insert. The earlier model placed the UMI next to P5 at the start of reads; that
  was a repository modelling error and is corrected.
- 🟢 GEO GSE210563 sample records (`GSE210563_family.soft.txt`, all 58 samples): "mRNA was
  purified by oligo (dT) beads. 100-200 ng of mRNA was used for GLORI treatment and
  construction of sequencing libraries"; "RNA libraries were prepared using modified eCLIP
  library construction, NEBNext Small RNA library kits (NEB) and SMARTer Stranded Total
  RNA-Seq kits (v3) - Pico Input mammalian (Takara) ... following manufacturer's
  protocols"; instruments Illumina NovaSeq 6000 (32) and HiSeq X Ten (26); processing
  trims adapters (trim_galore 0.6.6), removes duplicates (SeqKit) and then UMIs
  (FASTX-Toolkit). GSE233875 (the 2024 protocol's two HEK293T samples) names modified
  eCLIP only, HiSeq X Ten.
- 🟢 SRA run info (PRJNA866163, 58 runs; PRJNA978429, 2 runs): every deposited run is
  `SINGLE`, 150 bases per spot. This is deposition metadata, not a stated run
  configuration.
- 🟢 The authors' GLORI-tools README (linked from the GEO records,
  github.com/liucongcas/GLORI-tools, lines 26–35) accepts only single-end, A-to-G-converted
  reads and trims a 5′ UMI (`fastx_trimmer -f {length(5′UMI)+1nt}`) after UMI-based
  duplicate removal.
- 🟡 Derived: the deposited single-end read begins with the UMI and is RNA-sense, which
  matches the read primed from the Read 2 site of the DNA adaptor in the model. Whether
  a Read 1 was also sequenced but not deposited was open before the 2023 Methods were read:
  they state 2 × 150 paired-end runs with only R2 analysed.
- 🔴 (not printed in the inspected texts; not confirmed absent) the bases of the 2023 "P5 and Illumina-indexed P7 primers" and the 2024 NEBNext Universal/Index primers, and the sequencing primers. The final outer arms and sample index remain inferred placeholders; the NEBNext E7335 vendor oligo document was not read.

## Oligo source boundary

🟢 Detailed protocol Supplementary Table 1 (re-read 2026-10-17 from the original XLSX)
prints, as ordered: 3′ RNA adaptor
`rAppAGATCGGAAGAGCGTCGTG`, RT primer `ACACGACGCTCTTCCGATCT`, and 3′ cDNA adaptor
`Phos-NNNNNNNNNNAGATCGGAAGAGCACACGTCTG`. The `N` span is the 10-nt UMI. The 2023 Methods
print the same three sequences, with `-3SpC3` on the linker and the UMI adapter. The page
now uses the 2023 strings with their modifications. The 2024 table omits the 3′ blocks, but
its Reagents text states them, so this is not a source conflict.

🔴 The complete modified-eCLIP library PCR oligo bases remain unprinted in the inspected texts (named only as P5 / Illumina-indexed P7, and NEBNext primers in 2024). The detailed
protocol names NEBNext Small RNA and SMARTer kit alternatives, but the defining branch is
modified eCLIP; those kit oligos cannot be substituted for it. Raw FASTQ from GSE210563
could support a read-visible adapter consensus and UMI boundary, but cannot establish an
ordered oligo's end modifications or unsequenced tail. Keep only the unverified outer
library arms inferred.

## Construct logic

GLORI is a negative conversion assay: an unmodified genomic A should be read as G, while
the retained-A fraction estimates m6A stoichiometry. The model therefore keeps chemical
state distinct from the downstream library branch and carries UMI identity explicitly.

## Bounded source pass 2026-10-17 (round breadth-first-10)

Inspected: 2023 MOESM1 Supplementary Fig. 1 and 14 legends; 2023 MOESM2 Reporting
Summary p. 1 (rendered); 2024 MOESM3 Supplementary Tables 1–2. Corrections: final role
order (UMI beside the Read 2 site, not P5); RNA adaptor and RT primer tagged as Read 1
side; instrument now HiSeq X Ten and NovaSeq per the Reporting Summary; unverified
reaction conditions downgraded to carried claims. Deferred: 2023 Methods and 2024 main
text (paywalled), GEO GSE210563 library descriptions, library PCR primers and read
configuration, ligase/RT enzyme identities, adaptor 3′ blocking groups (none printed in
Table 1). No missing-information claim is accepted by this capped pass.


## Bounded source pass 2026-10-17 (round deferred-followup-02)

Routes for the defining originals: OpenAlex and Unpaywall report all three DOIs
(2023 article, 2024 Nature Protocols) closed with no repository copy; Europe PMC has no
PMCID for either; no bioRxiv/Research Square preprint found; both main texts remain
publisher-only (manual URLs above). Read instead: GEO GSE210563 and GSE233875 SOFT
records, SRA run info for PRJNA866163/PRJNA978429 (all placed in the sources directory
with MANIFEST rows) and the GLORI-tools README. Corrections: sequencing text now states
the GEO library branches and instruments, single-end 150-nt deposited runs and the
README's 5′-UMI single-end read, labelled as deposition/derived evidence. Still
deferred: 2023 Methods and 2024 main text, library PCR primers, run cycle settings,
ligase/RT enzymes, reaction conditions. No missing-information claim is accepted.

## Bounded source pass 2026-10-17 (round deferred-followup-03)

The defining Liu 2023 article and the Shen 2024 Nature Protocols PDF were hand-placed in the
source bundle; their printed DOIs, titles and authors match. Read: 2023 Methods "RNA
protection, deamination and deprotection", "Library preparation and sequencing" and data
filtering (lines 2612–2696); 2024 overview (360–391), Reagents (1834–1898), Step 23 (2448)
and Steps 92–93 (3070–3200). Corrections: graph steps now carry the 2023 conditions and
enzymes; "repair" is a 3′ dephosphorylation; excess-linker and excess-primer digestions,
1.8× XP cleanup and the 165–250 bp gel are explicit, with waste; the linker, RT primer and
UMI adapter are exact graph inputs, with the 3′ C3 blocks printed. The read layout is
2 × 150 paired-end, with R2 analysed. The unsupported "4 °C overnight" claim was removed.
Still deferred: PCR primer bases (NEBNext E7335 vendor document), index read length, 2024
Procedure steps beyond those listed, and the kit-alternative branches.
