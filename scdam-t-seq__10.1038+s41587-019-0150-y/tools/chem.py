"""Source-exact scDam&T-seq DamID and CEL-seq2 molecular states."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, MolecularState, Row, Scene, Segment, Workflow,
                      feature, revcomp)
from multibranch import illumina_branch
from pcr import PCRGeometry, PCRPrimer

TITLE = "scDam&T-seq — protein–DNA contacts and transcriptome"
NOTES = "01_scdam-t-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41587-019-0150-y">'
          'Rooijers et al., <i>Nature Biotechnology</i> (2019)</a>; detailed protocol '
          '<a href="https://doi.org/10.1038/s41596-020-0314-8">'
          '<i>Nature Protocols</i> (2020)</a>.')
SUMMARY = ("One well-specific CEL-seq2 primer labels transcript cDNA while a different "
           "well-specific forked adapter labels DpnI-cut Dam-methylated DNA. Both carry "
           "T7 promoters, so the two branches are pooled and linearly amplified together.")

T7 = "GCCGGTAATACGACTCACTATAGGG"
RA5 = "AGTTCTACAGTCCGACGATC"
UMI1, UMI2 = "NNN", "NNN"
RNA_BC1, RNA_BC2 = "GATG", "TCAT"       # protocol's printed example
DNA_BC1, DNA_BC2 = "TGCA", "TATG"       # protocol's printed example
POLYT = "T" * 25
CEL_PRIMER = T7 + RA5 + UMI1 + RNA_BC1 + UMI2 + RNA_BC2 + POLYT + "V"

DAM_FORK_TOP = "GGTGAT"
DAM_FORK_BOTTOM = "GAGCTT"
DAM_CORE_BODY = "CCGGTAATACGACTCACTATAGGGGTTCAGAGTTCTACAGTCCGACGATC"
DAM_TOP = DAM_FORK_TOP + DAM_CORE_BODY + UMI1 + DNA_BC1 + UMI2 + DNA_BC2 + "GA"
DAM_BOTTOM = ("TC" + revcomp(DNA_BC2) + UMI2 + revcomp(DNA_BC1) + UMI1
              + revcomp(DAM_CORE_BODY) + DAM_FORK_BOTTOM)

RANDOMHEX_TAIL = "GCCTTGGCACCCGAGAATTCCA"
RANDOMHEX_RT = RANDOMHEX_TAIL + "N" * 6
RNA_PCR1 = "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"
RNA_PCR_INDEX = ("CAAGCAGAAGACGGCATACGAGAT" + "N" * 6
                 + "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA")
RNA_PCR_INDEX_ANNEAL = "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"


def _umi(name):
    return Segment(name, "N" * 3, "umi",
                   feature=feature("umi", "umi", "random", group="umi",
                                   part=name))


def _barcode(name, bases, identity):
    return Segment(name, bases, "cbc",
                   feature=feature(identity, "cell_barcode", "fixed",
                                   group=identity, part=name))


def cel_primer_segments():
    return [Segment("CEL T7 promoter", T7), Segment("CEL RA5", RA5, "r1"),
            _umi("CEL UMI part 1"), _barcode("CEL barcode part 1", RNA_BC1, "cel_cell"),
            _umi("CEL UMI part 2"), _barcode("CEL barcode part 2", RNA_BC2, "cel_cell"),
            Segment("CEL poly(dT)25", POLYT),
            Segment("CEL anchoring V", "V", placeholder=True)]


def cel_rt_scene():
    sc = Scene()
    sc.strand("mRNA", [Segment("mRNA body", "X" * 40, placeholder=True),
                        Segment("poly(A)-proximal B", "B", placeholder=True),
                        Segment("poly(A) tail", "A" * 25)], label="polyadenylated mRNA")
    sc.anneal("CEL-seq2 primer", cel_primer_segments(), to="mRNA",
              pair=("CEL anchoring V", "poly(A)-proximal B"),
              label="well-specific CEL-seq2 primer")
    sc.arrow("CEL-seq2 primer", "reverse transcriptase copies the transcript")
    sc.labels("CEL-seq2 primer")
    return sc


def damid_adapter_scene():
    top = [Segment("top six-base fork", DAM_FORK_TOP),
           Segment("DamID constant core", DAM_CORE_BODY), _umi("DamID UMI part 1"),
           _barcode("DamID barcode part 1", DNA_BC1, "damid_cell"),
           _umi("DamID UMI part 2"),
           _barcode("DamID barcode part 2", DNA_BC2, "damid_cell"),
           Segment("terminal GA", "GA")]
    bottom = [Segment("terminal TC", "TC"),
              _barcode("DamID barcode part 2 bottom", revcomp(DNA_BC2), "damid_cell"),
              _umi("DamID UMI part 2 bottom"),
              _barcode("DamID barcode part 1 bottom", revcomp(DNA_BC1), "damid_cell"),
              _umi("DamID UMI part 1 bottom"),
              Segment("DamID constant-core bottom", revcomp(DAM_CORE_BODY)),
              Segment("bottom six-base fork", DAM_FORK_BOTTOM)]
    sc = Scene(); sc.strand("DamID top", top, label="DamID top oligo")
    sc.anneal("DamID bottom", bottom, to="DamID top",
              pair=("terminal TC", "terminal GA"), label="DamID bottom oligo",
              mod5="/5Phos/", unpaired=("bottom six-base fork",))
    sc.mark("DamID top", "top six-base fork", "non-complementary six-base fork")
    sc.mark("DamID bottom", "bottom six-base fork", "prevents adapter concatemers")
    sc.ligation_end("DamID top", "terminal GA", kind="blunt", owner="adapter",
                    text="blunt terminal GA/TC half-GATC end")
    sc.labels("DamID top")
    return sc


def dpni_fragment_scene():
    sc = Scene.duplex([Segment("left DpnI half-site", "TC"),
                       Segment("methylated-GATC genomic fragment", "X" * 36,
                               placeholder=True),
                       Segment("right DpnI half-site", "GA")],
                      label="DpnI fragment with blunt half-GATC ends")
    sc.ligation_end("top", "left DpnI half-site", kind="blunt", owner="substrate",
                    text="blunt DpnI substrate end")
    sc.mark("top", "left DpnI half-site", "DpnI cut G^mATC leaves TC on this end")
    sc.labels("top")
    return sc


def damid_geometry():
    return AdapterLigationGeometry.known(
        "blunt", damid_adapter_scene(), dpni_fragment_scene(),
        statement="The detailed protocol prints both adapter strands and states that DpnI leaves blunt ends; terminal adapter GA plus genomic TC restores GATC.")


def ligated_damid_scene():
    top = [Segment("top six-base fork", DAM_FORK_TOP),
           Segment("DamID constant core", DAM_CORE_BODY), _umi("DamID UMI part 1"),
           _barcode("DamID barcode part 1", DNA_BC1, "damid_cell"),
           _umi("DamID UMI part 2"),
           _barcode("DamID barcode part 2", DNA_BC2, "damid_cell"),
           Segment("adapter half-site GA", "GA"), Segment("genomic half-site TC", "TC"),
           Segment("DpnI-selected genomic DNA", "X" * 40, placeholder=True)]
    bottom = [Segment("genomic bottom", "X" * 40, placeholder=True),
              Segment("genomic half-site bottom", "GA"),
              Segment("adapter half-site bottom", "TC"),
              _barcode("DamID barcode part 2 bottom", revcomp(DNA_BC2), "damid_cell"),
              _umi("DamID UMI part 2 bottom"),
              _barcode("DamID barcode part 1 bottom", revcomp(DNA_BC1), "damid_cell"),
              _umi("DamID UMI part 1 bottom"),
              Segment("DamID constant-core bottom", revcomp(DAM_CORE_BODY)),
              Segment("bottom six-base fork", DAM_FORK_BOTTOM)]
    sc = Scene(); sc.strand("top", top, label="one adapter-ligated DpnI end")
    sc.anneal("bottom", bottom, to="top",
              pair=("genomic bottom", "DpnI-selected genomic DNA"),
              label="one adapter-ligated DpnI end",
              unpaired=("bottom six-base fork",))
    sc.junction("top", "adapter half-site GA", "genomic half-site TC",
                "blunt ligation restores GATC")
    sc.mark("top", "DamID barcode part 1", "DamID-specific well barcode")
    sc.labels("top")
    return sc


def arna_scene():
    sc = Scene(); sc.strand("DamID aRNA", [Segment("T7-derived 5-prime G", "G"),
        Segment("RA5 RNA", RA5.replace("T", "U"), "r1"), _umi("aRNA UMI part 1"),
        _barcode("aRNA barcode part 1", DNA_BC1, "damid_cell"),
        _umi("aRNA UMI part 2"), _barcode("aRNA barcode part 2", DNA_BC2, "damid_cell"),
        Segment("copied GATC-selected genomic RNA", "X" * 42, placeholder=True)],
        label="linearly amplified DamID RNA")
    sc.arrow("DamID aRNA", "T7 transcription")
    sc.labels("DamID aRNA")
    return sc


def randomhex_scene():
    sc = Scene(); sc.strand("aRNA fragment", [Segment("aRNA RA5 and identifiers", "X" * 32,
                    placeholder=True), Segment("aRNA payload", "X" * 42, placeholder=True),
                    Segment("random priming site", "N" * 6)], label="fragmented aRNA")
    sc.anneal("randomhexRT", [Segment("randomhexRT RA3 tail", RANDOMHEX_TAIL, "r2"),
                               Segment("random hexamer", "N" * 6)],
              to="aRNA fragment", pair=("random hexamer", "random priming site"),
              label="randomhexRT")
    sc.arrow("randomhexRT", "reverse transcriptase")
    sc.labels("randomhexRT")
    return sc


def pre_pcr_template():
    return Construct([Segment("RNA PCR1 annealing region", RNA_PCR1[-21:], "r1"),
                      Segment("copied barcode/UMI and RNA or DamID insert", "X" * 66,
                              placeholder=True),
                      Segment("randomhexRT-tail reverse complement",
                              revcomp(RNA_PCR_INDEX_ANNEAL), "r2")],
                     name="scDam&T randomhexRT cDNA")


def indexed_product_scene():
    parts = [Segment("P5 and RNA-PCR1", RNA_PCR1, "p5"),
             Segment("barcode/UMI and RNA or DamID insert", "X" * 66,
                     placeholder=True),
             Segment("RNA-PCR-index reverse complement", revcomp(RNA_PCR_INDEX), "p7")]
    sc = Scene.duplex(parts, label="indexed scDam&T-seq library")
    sc.mark("top", "barcode/UMI and RNA or DamID insert",
            "molecular branch is identified by its non-overlapping barcode whitelist")
    sc.labels("top")
    return sc


RNA, RP = illumina_branch("scDam&T-seq RNA", "CEL-seq2 transcript insert", umi=6)
DNA, DP = illumina_branch("scDam&T-seq DamID", "methylated-GATC genomic insert")
FINAL_LIBRARIES = (
    ("Transcriptome library", RNA, RP, "CEL-seq2 transcript-derived endpoint.",
     "Illumina primers on the transcript branch."),
    ("DamID library", DNA, DP, "DpnI-selected genomic endpoint.",
     "Illumina primers on the DamID branch."),
)


def oligos():
    return [f"<b>CEL-seq2 example primer</b> — 5′-{CEL_PRIMER}-3′",
            f"<b>DamID example top</b> — 5′-{DAM_TOP}-3′",
            f"<b>DamID example bottom</b> — 5′-/5Phos/{DAM_BOTTOM}-3′",
            f"<b>randomhexRT</b> — 5′-{RANDOMHEX_RT}-3′",
            f"<b>RNA PCR primer 1</b> — 5′-{RNA_PCR1}-3′",
            f"<b>RNA PCR index primer</b> — 5′-{RNA_PCR_INDEX}-3′"]


def workflow():
    w = Workflow(MolecularState("Dam-marked cell", (Row(chunks=[
        ("poly(A) RNA + genomic DNA containing protein-proximal G^mATC", None, False)]),)))
    w.react("Reverse-transcribe poly(A) RNA with a well-specific CEL-seq2 primer",
            cel_rt_scene().rows(), name="barcoded first-strand transcript cDNA",
            note="The printed primer contains T7, RA5, split 6-nt UMI, split 8-nt cell barcode, dT25 and anchoring V.")
    w.react("Make transcript cDNA double-stranded and cut methylated GATC with DpnI",
            [*cel_rt_scene().rows(), Row(), *dpni_fragment_scene().rows()],
            name="double-stranded transcript cDNA plus blunt DpnI fragments",
            note="DpnI cuts fully methylated GATC sites into blunt GA/TC half-sites; transcript cDNA remains in the same well.")
    w.ligate_adapter("Ligate the well-specific forked DamID adapter",
                     damid_geometry(), ligated_damid_scene().rows(),
                     name="barcoded DamID DNA and CEL-seq2 cDNA",
                     note="Adapter GA and genomic TC restore GATC at the exact covalent junction; RNA and DNA use intentionally non-overlapping barcode sets.")
    w.react("Pool wells and linearly amplify both branches by T7 IVT",
            arna_scene().rows(), name="mixed RNA- and DamID-derived aRNA",
            note="Both the CEL-seq2 primer and DamID adapter carry a T7 promoter, permitting one pooled IVT.")
    w.react("Fragment aRNA and reverse-transcribe with randomhexRT",
            randomhex_scene().rows(), name="randomhexRT cDNA",
            note="The exact 22-nt RA3-derived tail is a non-annealing 5′ primer tail; its N6 end chooses a random aRNA site.")
    geometry = PCRGeometry.exact(
        pre_pcr_template(), PCRPrimer("RNA PCR primer 1", RNA_PCR1),
        PCRPrimer("RNA PCR index primer", RNA_PCR_INDEX),
        statement="The detailed protocol prints both tailed PCR primers; their annealing suffixes face inward around either branch's randomhexRT cDNA.")
    w.pcr("Add P5, six-base sample index and P7 by terminal PCR", geometry,
          indexed_product_scene().rows(), name="indexed mixed scDam&T-seq library",
          note="Branch-specific cell barcode sets distinguish CEL-seq2 and DamID reads after sequencing.")
    return w
