# Duplex Sequencing

## 1. What it is

Duplex Sequencing attaches complementary random tags to both strands of each original DNA duplex. Reads are grouped into single-strand consensuses (SSCS); reciprocal tag pairs associate the two original strands for duplex consensus (DCS).

## 2. Original sources and version boundary

🟢 Schmitt *et al.*, *Detection of ultra-rare mutations by next-generation sequencing*, PNAS 109:14508–14513 (2012), [10.1073/pnas.1208715109](https://doi.org/10.1073/pnas.1208715109).

🟢 A [public USPTO exhibit copy of the original article and SI](https://ptacts.uspto.gov/ptacts/public-informations/petitions/1558377/download-documents?artifactId=Xfv653zneGBqP8v2e6wy_WkvKbPgmB-RXmi2RhE1DK4x45zyth1sMRs) supplies nine pages: the six-page article and three-page SI. The SI explicitly says **“Supporting Information Corrected March 27, 2013.”** The model uses that version. Its first page was visually inspected as well as read in extracted text (1010–1107); the remaining processing overview, worked tag examples, figure and yield-table text were read too. The pre-correction SI has not been compared, so the correction's exact changes and identity with the initial 2012 synthesis are not asserted.

The old local claim that the adapter/PCR bases could not be recovered is superseded by this source. This was a repository source-recovery limitation. It does not establish a publication error.

## 3. Source-owned oligos and adapter construction

🟢 The corrected SI prints a 58-base adapter primer strand and a phosphorylated 50-base template strand. The page lists them verbatim, including `/5phos/`. The template is `ACTG` + N12 + the reverse complement of the canonical TruSeq Read 2 primer. The primer strand factors as P5 without its final four bases followed by the complete TruSeq Read 1 primer; the overlap of their common sequence restores the full P5 site in the library.

🟢 The two source oligos anneal over their shared 13-base core. Klenow exo− copies the random tag and the template's `ACTG` terminus. The SI explicitly reports a mixture of blunt and A-tailed products after this copying step, followed by a separate dATP-only Klenow treatment that A-tails residual blunt adapters.

🟡 Copying `ACTG` gives `CAGT`; adding the adapter dA gives `CAGTA`. The copied adapter has a 29-position duplex core (13 + 12 + 4) and an **unpaired** terminal dA. Its unequal outer Y arms remain unpaired. The model computes the pairing from the source oligos and draws the one-base tail separately, avoiding the previous incorrectly paired A-overhang drawing.

🟢 Corrected-SI PCR primers are a 21-base forward primer and P7 + **six fixed multiplexing-index positions** + the 23-base arm `GTGACTGGAGTTCAGACGTGTGC`. The actual index values are unassigned here. The page does not import an eight-position index default from a generic TruSeq-library constructor.

## 4. Molecular path and junctions

🟢 Corrected SI specifies shearing, NEB DNA End Repair kit treatment, two-stage AMPure size selection, fragment dT-tailing with Klenow exo−, T4 DNA ligase adapter attachment and KAPA HiFi PCR. The 2026-10-17 bounded pass re-read SI page 1 (text 1048–1093) and added the source's bead purifications to the graph: 2 volumes AMPure XP after Covaris shearing, the two-sided 0.7 + 0.65-volume selection (>500 bp bound and discarded; ~200 bp or greater retained), 1.2 volumes after dT-tailing and 0.8 volumes after the 15 min, 25 °C ligation. The SI does not itemize what each cleanup removes, so the graph names only bead supernatant/wash waste. PCR amplified 40 attomoles of adapter-ligated DNA for 18–20 cycles; the SI says input and cycle number must be titrated to give at least three duplicates per tag family. The two ethanol precipitations in adapter synthesis are not drawn as separate steps. Adapter dA pairs with fragment dT. The later [Kennedy 2014 protocol](https://doi.org/10.1038/nprot.2014.170) explicitly redesigns the polarity to adapter dT/fragment dA; its oligos are not assigned to this corrected-original model.

🟡 A top-strand library has the sequence roles:

`primer-derived left arm — α(N12) — CAGTA — genomic insert — TACTG — β(N12) — template-derived Read 2 arm — six-position index complement — P7 complement`.

`TACTG` is the reverse complement of `CAGTA`. Reading inward from either end therefore gives N12 followed by CAGTA and then genomic sequence. The original ligated molecule has asymmetric Y arms: its opposite strand is not a complete complement of both outer arms. The page shows this physical state separately. During PCR, reverse-primer extension first copies the left primer arm, making a forward-primer site available. The exact inward-site panel describes that copied duplex core, with the reverse primer's P7/index flap shown separately; it is not presented as the initial Y-arm duplex.

## 5. Readout and computational filtering

🟢 Corrected SI page 2, processing overview (text 1144 onward), explicitly prints **CAGTA after exactly 12 random nucleotides**. Page 1 places the five fixed bases at read positions 13–17. The arithmetic above therefore agrees with a separately printed source sequence.

🟢 The two end tags form a 24-base paired identifier. The worked example explains reciprocal αβ/βα families and Read 2's complementary orientation. Processing removes the tags, the five fixed bases and four additional genomic bases at each end before alignment/consensus construction. These are computational operations; they do not physically remove library molecules. The graph's bead-size-selection waste is a separate physical fate.

🟢 The corrected SI names HiSeq 2000 and the TruSeq system. No run-cycle count is assigned from those platform names. The model derives primer sites and read starts while retaining the source's six-position sample-index design.

## 6. Remaining source scope

The original SI names the NEB DNA End Repair kit and manufacturer protocols for bead handling without giving a historical kit revision. The graph keeps end repair at the named-kit level rather than assigning individual enzymes from a current guide. The initial pre-correction SI, historical referenced vendor procedures, and their bounded relevance still need original-source reconciliation before a full audit gate can close. No missing-information claim is accepted.
