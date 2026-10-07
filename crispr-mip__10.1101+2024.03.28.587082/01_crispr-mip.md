# CRISPR-MIP — sgRNA quantification by molecular inversion probe

> **Evidence marking.** 🟢 verbatim from the preprint · 🟡 derived · ✅ **computed by us from
> the real lentiCRISPRv2 map** and asserted in `tools/selftest.py` · 🔴 not published.

**Selinger M, Yakovenko I, Nazir I, Henriksson J.** "CRISPR-MIP replaces PCR and reveals GC
and oversampling bias in pooled CRISPR screens." *bioRxiv* 2024.03.28.587082.
doi:[10.1101/2024.03.28.587082](https://doi.org/10.1101/2024.03.28.587082)

---

## 1. What it replaces, and why

🟢 A pooled screen's readout is normally a PCR off genomic DNA, counted as reads. That
conflates two things: how many cells carried a guide, and how well that guide's amplicon
amplified. The paper's finding is that the second is not negligible — **GC% bias affects
PCR**, and sequencing depth introduces an oversampling bias.

CRISPR-MIP replaces the readout PCR with a **padlock capture carrying a UMI**. Each captured
genomic molecule is labelled *before* any amplification, so the readout counts molecules
rather than reads. 🟢 The paper proposes it as a new gold standard for sgRNA quantification,
"especially for genes that are not top ranked but still of broad interest".

See [`../ref/concepts/padlock-circularization.md`](../ref/concepts/padlock-circularization.md)
for the general chemistry.

## 2. The probe 🟢

Table S2 is in `ref/TableS2_primers_and_probes.xlsx` (supplementary `media-2.xlsx`).
One ssDNA molecule, written 5'→3' as `ligation arm — backbone — extension arm`:

```
{{oligo: crisprmip.probe_construct()}}
```

(generated from `crisprmip.probe_construct()`; `/5Phos/` on the 5' end of the ligation arm.
Lengths: ligation arm {{= len(crisprmip.LIG_ARM) }} · UMI {{= crisprmip.UMI_LEN }} ·
Read 2 site {{= len(crisprmip.READ2_SITE) }} · i7 {{= crisprmip.INDEX_LEN }} ·
Read 1 site {{= len(crisprmip.READ1_SITE) }} · extension arm {{= len(crisprmip.EXT_ARM) }};
backbone {{= crisprmip.BACKBONE_LEN }}, total **{{= crisprmip.PROBE_LEN }} nt**.
Arm Tm: ligation {{= round(tm(crisprmip.LIG_ARM), 1) }} °C, extension
{{= round(tm(crisprmip.EXT_ARM), 1) }} °C.)

✅ Verified: the **Read 2 site reverse-complements exactly to the TruSeq Read 2 primer**,
and the **Read 1 site is the TruSeq Read 1 primer less its 5' `ACAC`**
({{= len(crisprmip.READ1_SITE) }} nt). The i7 is **{{= crisprmip.INDEX_LEN }} nt**, and the
**nine probes differ only in that index** (all checked against Table S2 in the selftest):
`{{= " ".join(crisprmip.PROBE_INDICES) }}`.

> ⚠ **The probe is {{= crisprmip.PROBE_LEN }} nt, not the
> {{= crisprmip.PROBE_LEN_IN_TEXT }} the Methods text states.** ✅ Table S2's sequences are
> {{= crisprmip.PROBE_LEN }} and decompose exactly as above. The
> {{= crisprmip.PROBE_LEN_IN_TEXT - crisprmip.PROBE_LEN }} nt discrepancy is in the paper's
> prose, not the oligo. We follow Table S2.

> ✅ **P5 and P7 are not in the probe at all** — they arrive on separate primers (§6).

**The Tm difference is deliberate.** 🟢 The ligation arm is the hotter of the two so that it
stays bound while the polymerase approaches it — otherwise the polymerase displaces the
ligation arm and the circle never closes. ✅ Asserted.

**Where the arms land**, on the real vector:

- ✅ The extension arm ends **exactly at the U6 +1 position** — so synthesis begins at the
  first base of the guide and runs straight through it.
- ✅ The ligation arm sits past the scaffold.
- ✅ **Exactly one capture site in the whole 13 kb vector.** Specificity comes from requiring
  both arms, in order, on the same strand, a fixed distance apart.

## 3. The captured molecule ✅

| | computed | published |
|---|---|---|
| gap the polymerase fills | **112 nt** | 🟢 112 bp |
| captured span incl. arms | 154 bp | — |
| closed circle | **238 nt** | — |

The gap-fill is `+1 G` + the 20 nt spacer + the scaffold + terminator.

> ⚠ **The 112 assumes a non-G-initiated spacer.** ✅ A spacer already starting with G gives
> **111**, because the vector's `+1 G` is then supplied by the spacer itself rather than
> appended. Brunello spacers are not G-initiated, so 112 is the working value — the same
> off-by-one that shifts the GPP amplicon sizes in `../lenticrispr-gecko-screen__10.1126+science.1247005/01_lenticrispr_gecko.md`.

## 4. Protocol 🟢

| Step | Conditions |
|---|---|
| **1. denature & hybridise** | 1–10 µg gDNA, 0.2–0.002 µM probe, 1.5× Ampligase buffer. 94 °C 5 min → ramp to 60 °C at −0.1 °C/s → hybridise **overnight at 60 °C** |
| **2. extend & ligate** | add Ampligase + 4 U Phusion HF + 0.2 mM dNTPs, pre-heated to 60 °C. **60 °C, 1 h** |
| **3. exonuclease** | **37 °C, 45 min**: 10 U Exonuclease I (removes un-circularised probe) + 50 U Exonuclease III (removes genomic DNA). Inactivate 80 °C, 20 min |
| **4. PCR** | amplify the whole reaction |

## 5. Four independent selection steps — none of them a cleanup

This is the part worth internalising. The method's specificity is enzymological, not
physical:

1. both arms must anneal to the same molecule, a fixed distance apart;
2. the polymerase must cross the gap **and** the ligase must close the circle;
3. 🟢 **exonuclease I/III destroy everything still linear** — unreacted probe, genomic DNA,
   and any capture that annealed but failed to ligate ✅;
4. 🟢 the **P5/P7 primer sites lie within the captured sequence**, not the backbone, so only
   a probe that was both extended *and* ligated can amplify.

No gel, no size selection, no bead cleanup decides what gets sequenced.

## 6. Amplification off the circle ✅

```
P5_tracrRNA_fwd   AATGATACGGCGACCACCGAGATCTACAC GTGCTTTTTTGAATTCGC
P7_tracrRNA_rev   CAAGCAGAAGACGGCATACGAGAT      GTTGATAACGGACTAGCC
```

✅ **Both 3' tails anneal inside the captured sequence** (the tracrRNA/scaffold), and
neither is present anywhere in the probe. That is the fourth selection filter made
concrete: a probe that annealed but was never extended has no primer sites at all.

✅ The circle gives **exactly one product, 269 bp**, and it is **inverse PCR** — the product
wraps the circle's origin, running outward from the captured region through the backbone.
It starts at P5, ends at revcomp(P7), and carries both the UMI and the sgRNA spacer.

## 7. Sequencing 🟢

**Read 1 = 60 cycles** (reads the captured sgRNA), **Read 2 = 15 cycles** (enough for the
13 nt UMI) ✅. Deduplication on the UMI gives absolute molecule counts, which is what
removes the sequencing-depth bias.

## 8. Screen context 🟢

Brunello **kinome** library (Addgene #75314), backbone **lentiCRISPRv2** (#52961) — the map
in `../lenticrispr-gecko-screen__10.1126+science.1247005/ref/plasmids/`, which is what all the ✅ values above are computed
against. Library amplified per the Broad "Amplification of pDNA Libraries" protocol and QC'd
by the Broad "PCR of sgRNAs for Illumina Sequencing" protocol — i.e. the GPP ARGON/BEAKER
scheme documented in `../lenticrispr-gecko-screen__10.1126+science.1247005/01_lenticrispr_gecko.md`. GFP-targeting sgRNA from
LentiGuide-Puro-GFPg1 (BB09).

## 9. Open 🔴

Cycle numbers for each PCR are set per sample by qPCR, so there is no fixed value to record.

Everything else in the construct is now resolved: Table S2 supplies the full probe, the
nine indices and both amplification primers, and `tools/selftest.py` recomputes the probe
decomposition, the capture, the circle and the 269 bp final library from the real vector on
every run. `tools/show_probe.py` prints the whole thing to scale.
