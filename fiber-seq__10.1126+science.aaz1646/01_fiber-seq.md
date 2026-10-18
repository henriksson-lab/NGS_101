# Fiber-seq — research notes

## Active native migration, 2026-10-18

The graph selects original2020 Hia5 marking of isolated nuclei:200U,10min at37°C, source0.8mMSAM in BufferA. It preserves both original genomic DNA strands with independent, unassigned methylation patterns. There is no invented poly(A)-like accessibility tract or fixed methylation mask. Protein/RNA populations and their digestion products are distinct from the marked DNA.

The source order is explicit: SDS stop; RNaseA/ProteinaseK50°C1h; first phenol/chloroform extraction and upper aqueous transfer; first sodium-acetate/ethanol precipitation overnight, pellet/70%ethanol wash/resuspension; second RNaseA37°C1h; second phenol/chloroform extraction; separate chloroform/isoamyl extraction; second ethanol precipitation1h, pellet/70%wash/resuspension; AMPure repurification; Megaruptor10–15kb target. These are successive physical recoveries, not one catch-all extraction node. Digestion does not itself physically discard RNA/protein products, and individual species-specific recovery remains unassigned.

Original2020 library preparation identity remains unassigned in the inspected passages. The SequelII system/version designation is not a preparation-kit identity. A bounded native PacBio preparation stage retains original modified DNA without importing the living bench protocol's SMRTbellkit3.0 repair/tailing/hairpin/nuclease recipe. The later25°C/commercial-Hia5 implementation is a distinct version, not a source discrepancy or proof of2020 execution. No arbitrary18base adapter, A/T ligation geometry, exact circular circuit, methylation copying or polymerase-entry polarity is supplied.

PacBio kinetics interrogate original modified template; a sequencer's normal DNA synthesis/read reporting is not inherited methylation on a copied preparative library. System1.0/Sequel,20U and untreated/deproteinized controls are separate source experiments, outside the selected graph. The original unassigned kit/geometry is a model boundary, not an accepted publication-absence claim.

This capped pass read 2 originals within the selected round and reran the mandatory bundle primer scraper. All historical review statuses and the independent two-reader missing-information gate are unchanged.

## Historical observations

The earlier material below retains its original scope; current source-owned declarations above replace family-selected geometry and later-version execution defaults.


## Sources read

- 🟢 Stergachis et al. 2020, *Science* 368:1449–1454, DOI [10.1126/science.aaz1646](https://doi.org/10.1126/science.aaz1646),
  defining paper (author-hosted PDF of the published article; main text pp1–6).
- 🟢 The article's Supplementary Materials (aaz1646_stergachis_sm.pdf, 22 pages; printed
  title, authors and DOI 10.1126/science.aaz1646 checked), hand-downloaded 2026-10-17.
  Materials and Methods read: MTase cloning/production, in vitro assay, 'AdMTase-seq'
  (nuclei and MTase reaction, text lines 76–150) and 'Fiber-seq' (lines 185–231).
- 🟢 GEO SuperSeries GSE146942 sample records (all 32 GSM, SOFT text), the data accession named in the article's data-availability statement; author-deposited sample treatment and extract protocols.
- 🟢 Fiber-seq bench protocol, linked by the authors' Fiber-seq documentation; current Hia5 and PacBio workflow.
- 🟢 PacBio SMRTbell prep kit 3.0 procedure, for the current PCR-free dumbbell endpoint.

## Molecular order

🟢 Hia5 methylates solvent-accessible adenines in intact nuclei using SAM. Protein-bound
DNA is protected, so the resulting m6dA pattern records accessibility and footprints on
individual native molecules. High-molecular-weight DNA is extracted without PCR. The
current workflow shears to long inserts, repairs and dA-tails the ends, ligates PacBio
hairpins, removes open molecules enzymatically, then anneals sequencing primer and binds
polymerase for circular-consensus sequencing.

## What the 2020 article itself states

🟢 Five nonspecific m6A-MTases were titrated on *Drosophila* S2 nuclei and compared with
DNase I–seq (AdMTase-seq, an m6A-immunoprecipitation readout of ~110-bp fragments, is a
separate validation assay, not the Fiber-seq library); Hia5 showed the highest
efficiency. Fiber-seq then treated S2 nuclei, built libraries by "polymerase chain
reaction–free library construction on high–molecular-weight DNA", and sequenced them by
PacBio CCS (fibers up to ~30 kb, mean 10.9 kb); K562 nuclei were profiled the same way.
The main text names neither the reaction conditions nor the SMRTbell kit.

## 2020 sample protocols in GEO GSE146942

🟢 Every Fiber-seq sample record (e.g. GSM4411218 "S2-DSRC Fiber-seq (500 U Hia5)",
GSM4411222–GSM4411225 K562) states: nuclei extracted with 0.05% IGEPAL in Buffer A
(50 mM Tris-HCl pH 8, 100 mM NaCl, 1 mM DTT) for 10 min at 4 °C; m6A-MTase
conditions "500 U Hia5 for 10min at 37C" (or 20 U for the low-dose samples); genomic
DNA "extracted then subjected to shearing and then PacBio library construction";
"PCR-free PacBio libraries were constructed using standard protocol"; "Long-read
PacBio HiFi/CCS sequencing" on Sequel or Sequel II. The AdMTase-seq validation samples
use 0.6–200 U of each MTase for 20 min at 37 °C followed by m6A immunoprecipitation
(antibody Millipore Sigma ABE572) and Illumina libraries. The article acknowledgements
say Hia5 was FPLC-purified in house from a donated plasmid (author PDF lines 670–675).

🟡 Version difference: 2020 marking at 37 °C versus 25 °C in the current bench protocol
(commercial Hia5). The page caveat states both; neither is called a source error.

## 2020 Supplementary Materials: the original Fiber-seq reaction

🟢 Nuclei: S2 cells in Buffer A plus an equal volume of cold 2× lysis buffer (0.1% IGEPAL
CA-630 in Buffer A), 10 min on ice, pelleted 350 × g; K562 nuclei "as previously
described". Fiber-seq used nuclei from 2.5 million S2 or 300,000 K562 cells, resuspended
in 57.5 µl Buffer A at 37 °C, with 1 µl MTase and 1.5 µl SAM (0.8 mM final), "with the
exception of using 200 units (2.5 μl) of Hia5", for 10 min. Stop: 3 µl 20% SDS (1%
final), then 130 µl Buffer A and 7 µl 20% SDS. RNase A and proteinase K for 1 h at 50 °C;
phenol:chloroform:isoamyl alcohol extraction and ethanol precipitation; a second RNase A,
phenol-chloroform and chloroform extraction and precipitation; AMPure re-purification;
Fragment Analyzer; Megaruptor shearing to 10–15 kb; PacBio library "using the Sequel II
System 2.0" (half of one K562 sample on System 1.0; 20 U replicates on Sequel I). No
SMRTbell prep kit is named. The deproteinized S2 control was marked as purified DNA with
500 U Hia5. The untreated S2 control was Covaris-sheared to ~1 kb, prepared with "Sequel
Sequencing Kit 3.0" and BluePippin size-selected.

🟡 Published discrepancy, shown on the page: the SM gives 200 U Hia5 for the main
Fiber-seq samples, whereas the GEO sample titles and characteristics give 500 U
(GSM4411218, GSM4411222–GSM4411224). The model lists both and does not choose.

🟡 The SM supersedes the earlier current-protocol labels for stop, extraction and
shearing. The repair, dA-tailing, hairpin-ligation and nuclease steps still follow the
current bench protocol because the SM names no SMRTbell kit.

## Current bench protocol (later implementation)

🟢 The authors' bench protocol (fiberseq.github.io, Parts 2–5, read 2026-10-17)
isolates nuclei with IGEPAL, marks them with Hia5 and SAM (NEB B9003S) at 25 °C for
10 min, quenches with SDS to ~1%, extracts HMW DNA (Qiagen MagAttract or Promega Wizard
HMW), shears on a Megaruptor 3 (or g-TUBE), and follows PacBio SMRTbell prep kit 3.0
("Repair, Ligate, Nuclease"), then anneals sequencing primer and binds polymerase for
Revio. 🟡 These are current values, not the 2020 conditions.

## Interpretation boundary

🟡 The schematic shows representative accessible and protected intervals, not a fixed
nucleosome spacing. 🟡 The defining 2020 experiment predates SMRTbell prep kit 3.0; the
endpoint shown is the current standard Fiber-seq library workflow, linked as an additional
protocol source. Hairpin sequences are proprietary.

## Review trail

- 2026-10-17 bounded pass (round 2026-10-17-breadth-first-10): re-read the published
  article main text (author-hosted PDF text lines 90–146, 214–232, 425–445, 665–749);
  Supplementary Materials re-attempted and still HTTP 403; read the current bench
  protocol Parts 2–5. Added the SDS quench step and current-protocol extraction and
  shearing labels to the graph, and made the page caveat say which steps come from the
  later protocol.
- 2026-10-17 bounded pass (round 2026-10-17-deferred-followup-02): Science
  Supplementary Materials still unobtainable (science.org SM URL HTTP 403; Stergachis
  lab Science referrer link science.sciencemag.org/cgi/rapidpdf/368/6498/1449?ijkey=…
  HTTP 403; not in PMC or Europe PMC; OpenAlex closed with no repository copy; Harvard
  DASH no record; Churchman lab hosts only the article; no bioRxiv preprint found; the
  article's Zenodo citation 10.5281/zenodo.3743228 is analysis code). Read the GEO
  GSE146942 sample records instead; nuclei buffer, Hia5 units/time/37 °C, shearing and
  Sequel/Sequel II now labelled from them in the graph and caveat.
- 2026-10-17 bounded pass (round deferred-followup-03): Science Supplementary Materials
  now on disk (hand download, identity checked) and read (Materials and Methods,
  AdMTase-seq and Fiber-seq paragraphs). Graph inputs now carry the 2020 nuclei, Hia5,
  SAM, SDS stop, extraction and Megaruptor values; page caveat rewritten; added the
  200 U vs 500 U Hia5 discrepancy. SMRTbell kit and Sequel binding/sequencing kit for the
  2020 libraries remain unnamed.
