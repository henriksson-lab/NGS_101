# scEdU-seq — selected single-pulse preparation

Defining paper: [Quantifying DNA replication speeds in single cells by scEdU-seq](https://doi.org/10.1038/s41592-024-02308-4), 2024.

This is a private, incomplete source-qualified model of the single-pulse RPE-1 route. It is not a complete reaction review, historical kit reconstruction or accepted missing-information audit.

The bounded source pass inspected published main54/150–188. The first supplementary PDF was a reporting summary, and a file-selection glob also opened the peer-review PDF header. These administrative documents consume two originals; the main consumes the third. The glob did not satisfy the intended before-open accounting for the peer-review header, but it was immediately recorded and no fourth original was opened. Published Table2/3 workbooks were deferred; their scraper pointers do not constitute reading or authorize adopting their sequences. No peer-review chemistry is used.

🟢 The selected cells receive10µM EdU15min, followed by trypsinization/PBS washing and75% ice-cold ethanol fixation. Short storage is−20°C up to24h; the longer storage-buffer branch is separate. Cells are BSA-blocked before in-situ50µl click reactions using the EdUClick647 kit, replacing azide647 with2mM azide–PEG3–biotin and adding6mM THPTA. Copper(I) is the catalyst, whereas azide–biotin becomes part of the labeled product. Exact historical kit copper formulation is unassembled.

🟢 PBS washing and DAPI/BSA suspension/20µm filtering precede index sorting into384wells with5µl mineral oil. The selected well is explicit; other cells/wells are not silently discarded. ProteinaseK lysis55°C2h/80°C20min precedes NlaIII37°C4h/65°C30min.

🟢 Main176 specifies Klenow-based end repair37°C30min/75°C20min, AmpliTaq360/dATP A-tailing72°C15min, then T7-adapter ligation4°C20min→16°C16h→65°C20min. The adapter amount is retained literally as50nl of5mM. The cell-barcode and UMI meanings are known, but exact adapter bases/windows are unassembled in this capped pass.

**Apparent source formulation conflict:** main54 names Klenow and polynucleotide kinase; main176 lists both1.6nl and23.3nl entries as PNK reaction buffer, then refers to inactivation of both enzymes. The graph records known Klenow and leaves the actual PNK addition/formulation unresolved. It does not silently replace a buffer entry with an enzyme.

🟢 Indexed wells are physically pooled, the aqueous phase separated from mineral oil, and EdU-biotin DNA captured on MyOne StreptavidinC1. Noncaptured material has an unresolved physical fate. Heat95°C retrieves the non-EdU complementary strand; bead-associated EdU strands remain a separate unresolved fraction. This is not idealized quantitative label/capture specificity.

🟢 The retrieved complementary DNA is primed while ramping0.1°C/s to20°C, using:

```
ATGCCGGTAATACGACTCAC
```

Klenow extension25°C45min/75°C20min and1:1 sample:AMPure cleanup/elution7µl precede MEGAscriptT7 IVT37°C12h. Main180's50mM dNTP value is retained literally. TurboDNase37°C15min precedes RNA fragmentation94°C2min, ice cooling/50mMEDTA and1:1RNAcleanXP cleanup/elution12µl. DNA digestion products are present until physical workup, rather than disappearing at the nuclease reaction.

🟢 Only5µl of the12µl RNA eluate enters RT; the remainder is explicit. The complete RT-primer pattern is:

```
GCCTTGGCACCCGAGAATTCCANNNNNN
```

The N6 portion is a random priming window, not an assigned UMI. Annealing65°C5min/ice precedes the single SuperscriptII addition and25°C10min→42°C60min→70°C20min. Both this external primer and the earlier Klenow primer are tracked as produced, co-retained same-tube components between annealing and extension, consumed as STRUCTURAL rather than re-added or assigned an intramolecular SELF_PRIMER role. RNaseA37°C30min precedes NEBNextUltraIIQ5 PCR with source RP1/RPIx, whose full sequences/cycle programme await the deferred workbooks. Two separate0.8:1 sample:AMPure cleanups precede Qubit/Bioanalyzer QC.

🟢 Main188 literally specifies NextSeq500/NextSeq2000 with100cyclesR1(cell index/UMI) and100cyclesR2(sample index). This wording is retained as a limitation; no molecular read mask or index orientation is manufactured. The read spans/actual barcode widths need the directly linked published originals and full library assembly.

No method ancestry is asserted from matching T7/small-RNA handles. The cited adapter predecessor and kit originals are unreviewed. Source items remain open, not confirmed absent.
