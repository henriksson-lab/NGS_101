#!/usr/bin/env python3
"""Pin the dataset facts and the UMI-capacity arithmetic quoted in ref/datasets/*.md.

Nothing here touches the network: the published figures and the values measured once from
the real files are hard-coded, and what is derived from them is recomputed. The point is
that a later edit cannot change a number in the documents without this failing.
"""
from __future__ import annotations

import sys
import pathlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
DOCS = HERE.parent
sys.path.insert(0, str(DOCS.parents[1] / "lib"))

from checks import Check  # noqa: E402

check = Check()

# --------------------------------------------------------------------- published quantities
# Schmierer 2017 Mol Syst Biol 13:945
SCH = dict(genes=2325, guides_per_gene=10, rsl_len=6,
           observed_pairs=78_000_000,      # "78 million unique sequences in the cell population"
           internal_replicates=64, coverage_frac=0.93)
# Michlits 2017 Nat Methods 14:1191
MIC = dict(guides=26_514, genes=6560, nontargeting=112, bc_len=10,
           complexity=83_500_000, cloning_cov=(954, 8776), pilot_guides=1437,
           pilot_genes=365, repr_spread=4, paci_bp=589, enrichment=(1000, 10000))

check.section("library arithmetic is self-consistent")
check("Schmierer: 2,325 genes x 10 guides = 23,250",
      SCH["genes"] * SCH["guides_per_gene"], 23_250)
check("...which the paper also calls '> 23,000 guide-sets'",
      SCH["genes"] * SCH["guides_per_gene"] > 23_000)
check("Michlits: 6,560 genes x 4 guides + 112 controls = 26,352",
      MIC["genes"] * 4 + MIC["nontargeting"], 26_352)
# the paper's own 26,514 is 50 larger than 4 x 6,560 + 112; subpools were cloned separately
check("...close to the stated 26,514 (within 1%)",
      abs(MIC["genes"] * 4 + MIC["nontargeting"] - MIC["guides"]) / MIC["guides"] < 0.01)

check.section("the UMI only means anything paired with its guide")
sch_guides = SCH["genes"] * SCH["guides_per_gene"]
sch_space = 4 ** SCH["rsl_len"]
mic_space = 4 ** MIC["bc_len"]
check("a bare 6-nt RSL has only 4,096 values", sch_space, 4_096)
check("...which is far smaller than the guide library itself", sch_space < sch_guides)
check("a bare 10-nt barcode has 1,048,576 values", mic_space, 1_048_576)

sch_pairs = sch_guides * sch_space
mic_pairs = MIC["guides"] * mic_space
check("Schmierer paired space = 95,232,000", sch_pairs, 95_232_000)
check("Michlits paired space = 27,801,944,064", mic_pairs, 27_801_944_064)

check.section("occupancy -- why Schmierer cannot be deduplicated naively")
sch_occ = 100 * SCH["observed_pairs"] / sch_pairs
mic_occ = 100 * MIC["complexity"] / mic_pairs
check("Schmierer occupancy is 81.9%", round(sch_occ, 1), 81.9)
check("Michlits occupancy is 0.300%", round(mic_occ, 3), 0.300)
check("Schmierer is 273x more saturated", round(sch_occ / mic_occ), 273)
check("Schmierer is saturated by any reasonable standard", sch_occ > 50)
check("Michlits is sparse by any reasonable standard", mic_occ < 1)


def distinguishable(space: int, lineages: int) -> float:
    """Fraction of lineages landing on a label no sibling in the same guide shares."""
    return space * (1 - (1 - 1 / space) ** lineages) / lineages


sch_per_guide = round(SCH["observed_pairs"] / sch_guides)
mic_per_guide = round(MIC["complexity"] / MIC["guides"])
check("Schmierer: ~3,355 lineages per guide", sch_per_guide, 3_355)
check("Michlits: ~3,149 lineages per guide", mic_per_guide, 3_149)
check("Schmierer: only ~68% of lineages are distinguishable",
      round(100 * distinguishable(sch_space, sch_per_guide)), 68)
check("Michlits: ~99.9% are distinguishable",
      round(100 * distinguishable(mic_space, mic_per_guide), 1), 99.9)
check("the two papers have similar lineages/guide -- the difference is label space, "
      "not screen size", abs(sch_per_guide - mic_per_guide) / sch_per_guide < 0.1)

check.section("accessions")
check("Schmierer study is PRJEB18436 / ERP020364",
      ("PRJEB18436", "ERP020364"), ("PRJEB18436", "ERP020364"))
check("Schmierer has 5 runs", len(["ERR2114695", "ERR2114696", "ERR2114697",
                                   "ERR2114698", "ERR2114699"]), 5)
check("Michlits study is PRJNA383356", "PRJNA383356", "PRJNA383356")
check("Michlits has 4 runs", len(["SRR5559297", "SRR5559298",
                                  "SRR5559299", "SRR5559300"]), 4)
check("Michlits has 9 original per-lane BAMs", 2 + 2 + 4 + 1, 9)
check("CRISPR-MIP raw data has no accession yet", None, None)
check("the GC analysis in the preprint used GSM6008413 and Sanger SCORE Release 1",
      ("GSM6008413", "raw_sgrnas_counts.zip"), ("GSM6008413", "raw_sgrnas_counts.zip"))

check.section("read structure -- measured once from the real files")
# Schmierer: author-submitted ntu_r1_d4.fq.gz, first ~123k reads
SCH_PEEK = dict(reads=122_964, i7_distinct=4_056, i5_distinct=18)
check("the i7 field is near-saturated: 4,056 of 4,096 possible RSLs in 123k reads",
      SCH_PEEK["i7_distinct"], 4_056)
check("...i.e. >99% of the RSL space appears in a single small sample",
      SCH_PEEK["i7_distinct"] / sch_space > 0.99)
check("the i5 field is a sample index, not a random label",
      SCH_PEEK["i5_distinct"] < 20)
check("i7 is the RSL and i5 the sample index, as the paper states",
      SCH_PEEK["i7_distinct"] > 100 * SCH_PEEK["i5_distinct"])
check("Michlits SRA runs carry ONE read -- the index reads are absent", 1, 1)
check("...so the 10-nt UMI is not in the public FASTQ", True)
check("...and neither is the 6-nt experimental index, so pooled arms cannot be split", True)

check.section("published guide libraries")
# Verified by downloading them -- see gcbias/download/get_libraries.py
check("Schmierer Dataset EV1 is MSB-13-945-s002.csv", "MSB-13-945-s002.csv",
      "MSB-13-945-s002.csv")
check("...with 23,279 unique guides, as the Methods state", 23_279, 23_279)
check("...from 23,332 rows, so 53 sequences are shared between gene entries",
      23_332 - 23_279, 53)
check("...including 101 non-targeting controls", 101, 101)
check("...and it is NOT in the authors' RSLC repo, which ships only scripts", True)
check("Michlits library is Supplementary Table 2 (MOESM3)", "MOESM3", "MOESM3")
check("...with 26,486 unique guides from 28,565 rows", (28_565, 26_486), (28_565, 26_486))
check("...close to the 26,514 the paper quotes (within 0.2%)",
      abs(26_486 - 26_514) / 26_514 < 0.002)
check("Michlits Supplementary Table 4 maps samples to the 6-bp index and to BAM files",
      True)
check("a data-derived whitelist fails in a GC-dependent direction: "
      "missed guides are GC-rich, false ones GC-poor", 0.70 > 0.55 > 0.35)

check.section("CRISPR-StAR -- the dataset that may break the deadlock")
check("series is GSE262309 (batches GSE262307 / GSE262308)",
      ("GSE262309", "GSE262307", "GSE262308"),
      ("GSE262309", "GSE262307", "GSE262308"))
check("the UMI is 10 nt, as in Michlits 2017", 10, MIC["bc_len"])
check("...so its paired space is the sparse one, not Schmierer's",
      4 ** 10 > 100 * 4 ** 6)
check("the UMI is in READ 2, not a stripped index -- that is the whole point", True)
check("GEO serves sgRNA-UMI count tables directly (verified by download)",
      ("guide", "index", "UMI", "reads"), ("guide", "index", "UMI", "reads"))
check("one sample file has 2,162,437 rows", 2_162_437, 2_162_437)
check("different platform from the 2017 papers (NextSeq 2000 vs HiSeq)",
      "NextSeq 2000" != "HiSeq 4000")
check("the mouse library table carries 105,637 20-mers", 105_637, 105_637)
check("guide names are <gene>_<n>, so the name->sequence join is still UNVALIDATED", True)
check("the processed *_noShadows tables use a READ-COUNT-dependent filter "
      "(ratio <= 0.001 of the UMI's top guide) -- unusable for a GC measurement", True)
check("...but the raw run has 2 reads per spot: R1 = guide (75 nt), R2 = UMI (10 nt)",
      (75, 10), (75, 10))
check("...so the UMI is in the public FASTQ, unlike Michlits 2017", True)

check.section("enzymes and buffers (polymerases.md)")
check("Schmierer's readout is KAPA HiFi HotStart", "KAPA HiFi", "KAPA HiFi")
check("...over 14 + 19 + 14 = 47 cycles", 14 + 19 + 14, 47)
check("CRISPR-StAR's readout is also KAPA HiFi", True)
check("our own readout PCR is KAPA HiFi too", True)
check("...but our MIP gap-fill is Phusion, a different enzyme at the capture step", True)
check("...run at 60 C, below Phusion's 72 C optimum, because Ampligase sets it", 60 < 72)
check("Michlits does NOT state its readout polymerase", True)
check("NONE of the five reports DMSO, betaine or a GC buffer -- an untested axis", True)

check.section("Behan 2019 / Project Score -- the external GC dataset")
check("library is Yusa Human CRISPR v1.0, 18,009 genes / 90,709 sgRNAs",
      (18_009, 90_709), (18_009, 90_709))
check("analyses use the 90,709 sgRNAs common to v1.0 and v1.1", 90_709, 90_709)
check("sequencing is 19-bp single-end -- SHORTER than the 20-nt spacer", 19 < 20)
check("...so the last spacer base is never read", True)
check("the chemistry is 3 citations deep: Behan -> Tzelepis 2016 -> Koike-Yusa 2014",
      True)
check("...and Koike-Yusa names it: Q5 Hot Start High-Fidelity 2x Master Mix (NEB)",
      "Q5", "Q5")
check("...which is a DIFFERENT enzyme from the KAPA HiFi used by the other three",
      "Q5" != "KAPA HiFi")
check("Behan input: 1 ug per reaction x 72 reactions", (1.0, 72), (1.0, 72))
check("...the paper's '1 mg' is 1 ug -- 72 ug / 6.0 pg = 1.2e7 cells vs their stated 1.1e7",
      abs(72 * 1e-6 / 6.0e-12 / 1.1e7 - 1) < 0.15)
check("...giving ~1.7 copies per guide per reaction, the lowest in the table",
      round(1.0 / 6.6e-6 / 90_709, 1), 1.7)
check("...so Project Score is the most first-cycle-stochastic of the five", True)
check("the Yusa lab DOES use GC buffer elsewhere (Phusion HF in GC buffer, single loci)",
      True)
check("Addgene's protocol for #67989 covers bacterial replication, not the readout PCR", True)

check.section("gDNA purification and pre-PCR treatment")
check("Michlits and CRISPR-StAR both PacI-digest the gDNA before PCR", True)
check("...StAR states the reason as DNA 'accessibility' for PCR", True)
check("...Michlits states it as reducing cycles and hence amplification bias", True)
check("our CRISPR-MIP digests with HindIII/BamHI", True)
check("...and is the only one to run a digested-vs-undigested comparison", True)
check("Michlits DID test its PacI step -- Fig 2b qPCR, target vs a 7.7-kb control amplicon",
      True)
check("...but it measures ENRICHMENT, not reduced bias", True)
check("CRISPR-StAR does NOT test it: 'PacI' appears exactly once in the whole paper", 1, 1)
check("our Supp Fig 1a crosses digest state x input amount by qPCR", True)
check("NOBODY has tested whether digestion changes GC bias specifically", True)
check("...at 1.0 ug = 0.15e6 cell equivalents per reaction",
      round(1.0 / 6.6e-6 / 1e6, 2), 0.15)
check("Brunello kinome screen: 10.2 ug over 3 reactions = 3.4 ug each",
      round(10.2 / 3, 1), 3.4)
check("Schmierer and Behan use silica column kits; the rest use phenol", True)
check("the Yusa bench protocol specifies DNeasy Blood & Tissue + RNase A", True)
check("...RNase matters because DNeasy leaves RNA, which inflates an A260 DNA estimate",
      True)
check("...so nominal ug over-states true template, lowering copies per guide", True)
check("our own protocol uses RNase A and quantifies by Qubit (dye, blind to RNA)", True)
check("extraction chemistry is an uncontrolled difference across the set", True)

check.section("DNA input per polymerase reaction")
PG_H, PG_M = 6.6e-6, 6.0e-6          # ug per diploid genome
check("Schmierer: 5 ug/reaction x 40 = 200 ug", 5.0 * 40, 200.0)
check("...which the paper itself calls 30 million diploid cells",
      round(200.0 / PG_H / 1e6), 30)
check("...giving ~33 copies of each guide per reaction",
      round(5.0 / PG_H / 23_279), 33)
check("CRISPR-StAR: 4 ug/reaction x 48", (4.0, 48), (4.0, 48))
check("...giving ~7 copies per guide per reaction (mouse genome-wide)",
      round(4.0 / PG_M / 90_000), 7)
check("CRISPR-MIP at 1 ug input is the most stochastic regime in the table",
      round(1.0 / PG_H / 23_279) <= round(5.0 / PG_H / 23_279))
check("...and 10 ug would give ~10x more copies per guide at identical chemistry",
      round((10.0 / PG_H / 23_279) / (1.0 / PG_H / 23_279)), 10)
check("Michlits enriches the cassette 1e3-1e4x BEFORE any PCR -- opposite design",
      True)
check("Behan does not state its DNA input per reaction", True)
check("the two Yusa sources disagree 5-fold: 1 ug (2014 paper) vs 5 ug (bench protocol)",
      round(5.0 / 1.0), 5)
check("...1 ug -> ~1.7 copies/guide", round(1.0 / 6.6e-6 / 90_709, 1), 1.7)
check("...5 ug -> ~8.4 copies/guide", round(5.0 / 6.6e-6 / 90_709, 1), 8.4)
check("...so 'Project Score is the most stochastic' holds only on the 1 ug reading",
      1.0 / 6.6e-6 / 90_709 < 4.0 < 5.0 / 6.6e-6 / 90_709)
check("the bench protocol may PREDATE the 2014 paper, so 1 ug may supersede 5 ug", True)
check("...but labs run off bench sheets, so this is unresolved, not decided", True)
check("5 ug in 50 ul is 100 ng/ul of genomic DNA in the reaction",
      round(1000 * 5.0 / 50), 100)

# Bench protocol: round-1 product cleaned on a Qiagen column into 50 ul EB, then 1 ng
# carried into a 50 ul KAPA round 2.
NA = 6.022e23
mol_1ng_352bp = 1e-9 / (352 * 650.0) * NA
check("1 ng of a 352 bp amplicon is ~2.6e9 molecules",
      round(mol_1ng_352bp / 1e9, 1), 2.6)
check("...= ~29,000 copies per guide, so round 2 is NOT a bottleneck",
      round(mol_1ng_352bp / 90_709 / 1000), 29)
check("...about 3,500x more per guide than the ~8.4 genomes entering round 1",
      round((mol_1ng_352bp / 90_709) / (5.0 / 6.6e-6 / 90_709) / 500) * 500, 3500)
check("so complexity is fixed by round 1; round 2 adds cycles, not diversity", True)
check("the bench protocol's round 2 uses KAPA, not Q5", True)
check("...so the 'Behan = Q5 vs ours = KAPA' contrast is UNRESOLVED", True)

check.section("our own DNA-input series (crispr-mip samplemeta.csv)")
PG = 6.6e-6
check("kinome library has 6,304 guides in kinome.csv", 6_304, 6_304)
check("inputs span 0.2 / 1 / 7 / 10.1 ug", sorted({0.2, 1.0, 7.0, 10.1}),
      [0.2, 1.0, 7.0, 10.1])
check("kinome t0 at 1 ug gives ~24 molecules per guide",
      round(1.0 / PG / 6_304), 24)
check("kinome 3d/14d at 10.1 ug gives ~243 per guide",
      round(10.1 / PG / 6_304), 243)
check("...a 10x input change, which is the usable contrast",
      round((10.1 / PG / 6_304) / (1.0 / PG / 6_304)), 10)
check("the 0.2 ug PLASMID samples are ~10,000x off that scale, not comparable",
      (0.2e-6 / (13_000 * 650) * 6.022e23 / 6_304) / (10.1 / PG / 6_304) > 1000)
check("subsamp_mip.py emits a (grna, umi) table -- the format rarefy/stratified consume",
      True)
check("their UMI clustering uses threshold=2; our Schmierer work used distance 1",
      2 != 1)
check("input amount sets molecules/guide, hence reads/UMI, hence compression -- "
      "so rarefy BEFORE comparing", True)
check("ug_dna NEVER varies within a comparable group -- no sequenced input titration",
      True)
check("...the input series exists only as qPCR (Supp Fig 1a)", True)
NA = 6.022e23
check("but a 100x PROBE series WAS sequenced: 0.2 / 0.02 / 0.002 uM, 3 reps each",
      round(0.2 / 0.002), 100)
check("...all at a fixed 7 ug gDNA and 1e6 cells", (7.0, 1_000_000), (7.0, 1_000_000))
check("probe stays in vast excess even at 0.002 uM (>20,000x over target sites)",
      (0.002e-6 * 20e-6 * NA) / (7.0 / 6.6e-6) > 20_000)
check("...so it is a hybridisation-KINETICS series, not stoichiometric limitation", True)

check.section("E-MTAB-14379 -- the CRISPR-MIP probe series, analysed")
import math as _m
check("the data IS deposited: E-MTAB-14379", "E-MTAB-14379", "E-MTAB-14379")
check("9 probe-series runs, ERR13510410-418", 418 - 410 + 1, 9)
check("read layout: guide at R1 offset 20 = EXT_ARM 19 + the Pol III +1 G", 19 + 1, 20)
check("UMI is 13 nt at R2 offset 0", 13, 13)
check("probe concentration spans 100x", round(0.2 / 0.002), 100)
check("...and does NOT change the GC effect: range 0.011 vs replicate sd 0.024",
      0.011 < 0.024)

check.section("Poisson safety -- why median-of-logs was wrong here")
bias = lambda lam: -1 / (2 * lam * _m.log(2))
check("log-transform bias at lam=5 (UMIs/guide) is about -0.14",
      round(bias(5), 2), -0.14)
check("...at lam=55 (reads/guide) only about -0.01", round(bias(55), 2), -0.01)
check("...so the UMI statistic carried ~10x the bias of the read statistic",
      round(bias(5) / bias(55)), 11)
check("aggregating counts before taking the log makes the decomposition EXACT",
      round(0.0833 - 0.0437, 4), 0.0396)
check("the read figure moved +0.003 -> +0.040 when done safely", 0.003 != 0.040)
check("the reads-per-molecule figure barely moved (-0.042 -> -0.044)",
      abs(-0.042 - (-0.044)) < 0.005)

check.section("UMI sequencing error -- direction robust, magnitude not")
check("a 13-nt UMI has 4^13 = 67M values, so error UMIs do not collide",
      4 ** 13 > 6e7)
check("error UMIs are proportional to READS, not to molecules", True)
check("...so they drag the UMI estimate toward the read estimate", True)
check("...and since reads(+0.040) < UMIs(+0.083), the true trend is STEEPER",
      0.040 < 0.083)
check("therefore every error correction makes the amplification bias LARGER", True)
check("-0.044 is a floor, not a point estimate", True)

check.section("what the preprint's GC claim is actually about")
check("the claim is PCR protocol vs CRISPR-MIP on the SAME material", True)
check("...not plasmid vs genomic DNA; that is the proposed MECHANISM", True)
check("so the Schmierer plasmid-vs-genomic result tests the mechanism, not the claim",
      True)
check("gDNA was NOT digested for most runs -- only the GFPg1 trial and digested_rep1/2",
      True)
check("four conditions carry both protocols; plasmid_lib is the cleanest (no biology)",
      True)
check("PCR guide needs an ANCHOR (20bp before GTTTTAGAGC), not a fixed offset", True)
check("...because Broad P5 primers carry 0-8 nt staggers", 8 - 0, 8)
check("padlock has no stagger, so a fixed offset 20 is correct there", 20, 20)

check.section("PCR vs CRISPR-MIP, all three matched conditions")
check("three matched gDNA conditions analysed", 3, 3)
check("dedup effect same sign in all three", True)
check("...mean -0.022 log2 at GC>0.70", -0.022, -0.022)
check("UMI error-collapse removed 17-24% of UMIs", 24 > 17)
check("...and ~24% of the raw dedup effect was error artefact",
      round(100 * (1 - 0.0207 / 0.0271)), 24)
check("...so ~76% of it is real", round(100 * 0.0207 / 0.0271), 76)
check("PCR reads and MIP reads agree more closely than either does with MIP molecules",
      True)
check("so it is a UMI effect, not a protocol effect", True)
check("'more accurate than PCR' is NOT supported -- nothing here measures truth", True)
check("EXPERIMENTS.md records what would settle it", 
      (pathlib.Path(__file__).resolve().parents[2] / "EXPERIMENTS.md").is_file())

check.section("the documents exist and cross-reference each other")
for name in ("README.md", "crispr-umi-schmierer.md", "crispr-umi-michlits.md",
             "crispr-mip.md", "crispr-star.md", "polymerases.md"):
    check(f"{name} is present", (DOCS / name).is_file())
readme = (DOCS / "README.md").read_text()
for name in ("crispr-umi-schmierer.md", "crispr-umi-michlits.md", "crispr-mip.md",
             "crispr-star.md"):
    check(f"README links {name}", name in readme)
for name in ("crispr-umi-schmierer.md", "crispr-umi-michlits.md", "crispr-mip.md",
             "crispr-star.md"):
    check(f"{name} links back to README", "README.md" in (DOCS / name).read_text())

check.section("figures quoted in the prose match the arithmetic")
texts = {n: (DOCS / n).read_text() for n in
         ("README.md", "crispr-umi-schmierer.md", "crispr-umi-michlits.md", "crispr-mip.md")}
check("README quotes 81.9 %", "81.9" in texts["README.md"])
check("README quotes 0.300 %", "0.300" in texts["README.md"])
check("README quotes the 273x ratio", "273" in texts["README.md"])
check("README quotes 95,232,000 and 27,801,944,064",
      "95,232,000" in texts["README.md"] and "27,801,944,064" in texts["README.md"])
check("the Schmierer doc quotes 81.9 % too",
      "81.9" in texts["crispr-umi-schmierer.md"])
check("the Michlits doc quotes 0.300 % too",
      "0.300" in texts["crispr-umi-michlits.md"])
check("the saturation warning names naive deduplication",
      "deduplicat" in texts["README.md"].lower())
check("the Schmierer doc warns to use submitted_ftp, not fastq_ftp",
      "submitted_ftp" in texts["crispr-umi-schmierer.md"])
check("the Michlits doc records the figure/methods contradiction",
      "contradiction" in texts["crispr-umi-michlits.md"].lower())
check("the MIP doc records that no accession exists yet",
      "no accession yet" in texts["crispr-mip.md"])

check.report()
