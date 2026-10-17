"""DOGMA-seq molecular model.

Authority: Mimitou et al., Nat. Biotechnol. 39:1246 (2021), PMC8763625 Methods
"Cell permeabilization for DOGMA-seq" and "Transposition and barcoding for DOGMA-seq";
publisher supplementary protocol MOESM1 (DOGMA-seq protocol, Jan 2021); and the cited
10x CG000338 Rev A Multiome guide (Introduction pp17-18, Appendix p79, sequencing p68).
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

ARTICLE = '<a href="https://doi.org/10.1038/s41587-021-00927-2">Mimitou et al., <i>Nature Biotechnology</i> (2021)</a>'
SOURCE = "Defining source: " + ARTICLE + "."
GUIDE = ('<a href="https://assets.ctfassets.net/an68im79xiti/5tnYerYC4FdDwYK0hlcfWj/cc46120ef02327249b5c42f24f44d699/'
         'CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevA.pdf">10x CG000338 Rev A, Appendix p. 79</a>')

CELL = feature("dogma_cell", "cell_barcode", "whitelist", whitelist="10x Multiome gel-bead barcode")
RNA_UMI = feature("dogma_rna_umi", "umi", "random")
ADT_UMI = feature("dogma_adt_umi", "umi", "random")
ANTIBODY = feature("dogma_antibody", "feature_barcode", "whitelist",
                   whitelist="TotalSeq-A antibody barcodes (MOESM3 Supplementary Table 1)")
RNA_I5 = feature("dogma_rna_i5", "sample_index", "unknown")
RNA_I7 = feature("dogma_rna_i7", "sample_index", "unknown")
ATAC_I7 = feature("dogma_atac_i7", "sample_index", "unknown")
ADT_I7 = feature("dogma_adt_i7", "sample_index", "unknown")

# CG000338 Rev A Appendix p79.
PARTIAL_R1 = il.TRUSEQ_READ1[11:]                 # CTACACGACGCTCTTCCGATCT
ATAC_SPACER = "CGCGTCTG"
# Article Methods / MOESM1 Oligos.
ADT_ADDITIVE = "CCTTGGCACCCGAGAATTCC"             # printed CCTTGGCACCCGAGAATT*C*C
SI_PCR = il.P5 + il.TRUSEQ_READ1[4:24]            # AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC
RPI_ARM = "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
TSA_HANDLE = RPI_ARM[12:]                         # small-RNA handle present in TSA tags (Methods)
assert PARTIAL_R1 == "CTACACGACGCTCTTCCGATCT"
assert SI_PCR == "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC"
assert TSA_HANDLE.startswith(ADT_ADDITIVE)


def _cell():
    return Segment("cell barcode", "N" * 16, "cbc", placeholder=True, feature=CELL)


def rna_bead_primer():
    return Construct([Segment("partial TruSeq Read 1", PARTIAL_R1, "r1"), _cell(),
                      Segment("UMI", "N" * 12, "umi", placeholder=True, feature=RNA_UMI),
                      Segment("poly(dT)", "T" * 30), Segment("VN anchor", "VN", placeholder=True, length_bp=2)],
                     name="Multiome gel-bead RNA primer PN-2000261")


def atac_bead_primer():
    return Construct([Segment("P5", il.P5, "p5"), _cell(), Segment("10x spacer", ATAC_SPACER)],
                     name="Multiome gel-bead ATAC primer")


def rna_library():
    """CG000338 Rev A p79 Gene Expression Library, transcribed as printed (25-nt P5 prefix)."""
    return Construct([Segment("P5 flank", il.P5[:-4], "p5"),
                      Segment("i5 sample index", "N" * 10, "cbc", placeholder=True, feature=RNA_I5),
                      Segment("TruSeq Read 1", il.TRUSEQ_READ1, "r1"), *rna_bead_primer().segments[1:],
                      Segment("cDNA insert", "X" * 40, placeholder=True),
                      Segment("dA junction", "A"), Segment("TruSeq Read 2 arm", il.INDEX1_PRIMER, "r2"),
                      Segment("i7 reverse complement", "N" * 10, "cbc", placeholder=True, feature=RNA_I7),
                      Segment("P7 reverse complement", il.P7_RC, "p7")], name="DOGMA-seq gene-expression library")


def atac_precursor():
    return Construct([Segment("P5", il.P5, "p5"), _cell(), Segment("10x spacer", ATAC_SPACER),
                      Segment("Read 1N", nx.READ1_PRIMER, "r1"),
                      Segment("accessible genomic DNA", "X" * 42, placeholder=True),
                      Segment("Read 2N reverse complement", nx.INDEX1_PRIMER, "r2")],
                     name="DOGMA-seq cell-barcoded ATAC fragment")


def atac_library():
    """CG000338 Rev A p79 Multiome ATAC Library."""
    return Construct([*atac_precursor().segments,
                      Segment("i7 reverse complement", "N" * 8, "cbc", placeholder=True, feature=ATAC_I7),
                      Segment("P7 reverse complement", il.P7_RC, "p7")], name="DOGMA-seq ATAC library")


def _adt_core():
    return [Segment("partial TruSeq Read 1", PARTIAL_R1, "r1"), _cell(),
            Segment("UMI", "N" * 12, "umi", placeholder=True, feature=ADT_UMI),
            Segment("poly(dT)", "T" * 30), Segment("VN anchor", "VN", placeholder=True, length_bp=2),
            Segment("antibody barcode reverse complement", "N" * 15, "cbc", placeholder=True,
                    feature=ANTIBODY, inferred=True)]


def adt_precursor():
    return Construct([*_adt_core(), Segment("TSA small-RNA handle reverse complement", revcomp(TSA_HANDLE), "r2")],
                     name="DOGMA-seq preamplified TotalSeq-A cDNA")


def adt_library():
    return Construct([Segment("P5 flank", il.P5[:-4], "p5"), Segment("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
                      *_adt_core()[1:], Segment("RPI-x arm reverse complement", revcomp(RPI_ARM), "r2"),
                      Segment("i7 reverse complement", "N" * 8, "cbc", placeholder=True, feature=ADT_I7),
                      Segment("P7 reverse complement", il.P7_RC, "p7")], name="DOGMA-seq TotalSeq-A antibody-tag library")


RNA = rna_library()
ATAC = atac_library()
ADT = adt_library()
RNA_P = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])
ATAC_P = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["I2"], sp.NEXTERA["R2"])
ADT_P = (sp.TRUSEQ["R1"],)
for _library, _primers in ((RNA, RNA_P), (ATAC, ATAC_P), (ADT, ADT_P)):
    _problems = sp.verify(_library, _primers, required_roles=tuple(p.role for p in _primers))
    if _problems:
        raise ValueError("; ".join(_problems))

TITLE = "DOGMA-seq — RNA, ATAC and protein from one cell"
NOTES = "01_dogma-seq.html"
SUMMARY = ("TotalSeq-A-stained cells are permeabilized (LLL with light fixation, or digitonin) and run through "
           "10x Multiome ATAC + Gene Expression; an ADT additive primer in preamplification and a separate "
           "SI-PCR/RPI-x reaction turn the poly(A)-captured antibody tags into a third library.")
CAVEAT = ("The cited " + GUIDE + " prints both gel-bead primers and the final ATAC and gene-expression libraries; "
          "the article and supplementary protocol print the ADT additive, SI-PCR and RPI-x primers. "
          "<strong>INFERRED:</strong> the TotalSeq-A conjugate sequence is not printed in the inspected sources; the "
          "antibody-tag library places the 15-nt MOESM3 barcode directly after the bead poly(dT)VN and before the "
          "small-RNA handle named in Methods, omitting the conjugate's unresolved poly(A) length and any linker. "
          "The acquired read cycles are not stated in the inspected article sections (the guide's p68 cycle "
          "recommendations are not assumed to be the run), and the tag-library Read 2/Index 1 primers are not printed; "
          "no cycle counts are assigned. The gene-expression P5 prefix is transcribed as the 25 printed bases; another "
          "10x guide was reported in this repository's TEA-seq audit to print 29 bases (not re-checked here). "
          "Hashtag (HTO additive/D7xx) libraries and the optional post-cDNA-amplification tag route are not drawn.")
ADT_LIMITATION = ("Read 2 and Index 1 primers for the small-RNA arm of the tag library are not printed in the "
                  "inspected sources, so their binding sites and antibody-barcode cycles are not asserted.")
FINAL_LIBRARIES = (
    ("RNA library", RNA, RNA_P, "CG000338 Rev A gene-expression library (dual 10-nt sample indexes).",
     "TruSeq Read 1 acquires cell barcode and UMI; TruSeq Read 2 reads the cDNA insert."),
    ("ATAC library", ATAC, ATAC_P, "CG000338 Rev A ATAC library: cell barcode and 8-nt spacer lie between P5 and Read 1N.",
     "Read 1N/Read 2N read the insert; the i5 read covers the spacer and cell barcode."),
    ("Antibody-tag library", ADT, ADT_P,
     "INFERRED — TotalSeq-A tag library from SI-PCR and RPI-x; barcode adjacency to the bead poly(dT)VN is inferred.",
     "TruSeq Read 1 acquires cell barcode and UMI."),
)


def oligos():
    return [
        "<b>Multiome RNA gel-bead primer PN-2000261 (CG000338 Rev A)</b> — 5′-" + PARTIAL_R1 + "-N16-N12-" + "T" * 30 + "VN-3′ (N16: cell barcode; N12: UMI)",
        "<b>Multiome ATAC gel-bead primer (CG000338 Rev A)</b> — 5′-" + il.P5 + "-N16-" + ATAC_SPACER + "-3′",
        "<b>ADT additive, preamplification spike (Methods; MOESM1)</b> — 5′-CCTTGGCACCCGAGAATT*C*C-3′ (*: phosphorothioate)",
        "<b>SI-PCR forward (Methods; MOESM1)</b> — 5′-" + SI_PCR + "-3′",
        "<b>RPI-x, TotalSeq-A tag index primer (Methods)</b> — 5′-" + il.P7 + "xxxxxxxx" + RPI_ARM + "-3′ (x: sample index)",
    ]


def _endpoint(g, prefix, template, final, primers, reviewed, *, geometry, forward_state=None,
              reverse_state=None, limitation=""):
    forward = g.input(prefix + "_forward", prefix + " indexing forward primer", MaterialKind.OLIGO,
                      state=forward_state, inferred=forward_state is None)
    reverse = g.input(prefix + "_reverse", prefix + " indexing reverse primer", MaterialKind.OLIGO,
                      state=reverse_state, inferred=reverse_state is None)
    polymerase = g.input(prefix + "_polymerase", prefix + " indexing polymerase", MaterialKind.CATALYST)
    library = g.product(prefix + "_library", final.name, MaterialKind.LIBRARY, state=final)
    g.react(prefix + "_pcr", "Index and amplify the " + prefix + " library", StepKind.PCR,
            ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,), geometry=geometry)
    g.sequence(library, SequencingPlan(library=final, primers=tuple(primers),
                                       acquired_roles=tuple(p.role for p in primers), limitation=limitation))
    reviewed.extend((forward.id, reverse.id, polymerase.id, prefix + "_pcr", library.id))


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    cells = g.input("cells", "single-cell suspension with RNA, chromatin and surface epitopes", MaterialKind.BIOLOGICAL,
                    state=molecule("DOGMA-seq cells", (role("poly(A) RNA"), role("chromatin"), role("surface epitopes"))),
                    primary=True)
    antibodies = g.input("antibodies", "TotalSeq-A antibodies after Fc block", MaterialKind.SELECTION, inferred=True,
                         source=SOURCE + " (conjugate sequences not printed; MOESM3 lists barcodes)")
    stained = g.product("stained", "stained cells retaining antibody tags", MaterialKind.MOLECULAR,
                        state=molecule("stained cells", (role("poly(A) RNA"), role("chromatin"), role("antibody tags"))))
    unbound = g.product("unbound", "unbound antibody conjugates", MaterialKind.WASTE)
    g.react("stain", "Stain cells, wash and FACS-sort live cells", StepKind.BINDING,
            ((cells, InputRole.TEMPLATE), (antibodies, InputRole.SELECTION)), (stained, unbound))
    g.discard(unbound, "removed by washes and sorting")
    lysis = g.input("lysis", "LLL (0.1% formaldehyde 5 min, then LLL lysis) or DIG (0.01% digitonin) permeabilization",
                    MaterialKind.CONDITION)
    permeable = g.product("permeabilized", "permeabilized stained cells", MaterialKind.MOLECULAR, state=stained.state)
    g.react("permeabilize", "Permeabilize stained cells (LLL or DIG)", StepKind.PERMEABILIZATION,
            ((stained, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)), (permeable,), changes_molecule=False)
    transposome = g.input("transposome", "Multiome ATAC transposition mix", MaterialKind.ADAPTER, inferred=True,
                          source=GUIDE + " (Read 1N/Read 2N ends printed in the library; loaded-complex state undeclared)")
    tn5 = g.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = g.product("tagged", "cells with tagmented chromatin, RNA and antibody tags", MaterialKind.MOLECULAR,
                       state=stained.state)
    twaste = g.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    g.react("tagment", "Tagment accessible chromatin", StepKind.TAGMENTATION,
            ((permeable, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)),
            (tagged, twaste))
    g.discard(twaste, "does not retain an amplifiable accessibility fragment")
    droplets = g.input("droplets", "10x Multiome GEM partitioning (Chip J)", MaterialKind.CONDITION)
    occupied = g.product("occupied", "GEMs with cell, gel bead and RT master mix", MaterialKind.MOLECULAR, state=tagged.state)
    failed = g.product("failed", "empty or multiply occupied GEM products", MaterialKind.WASTE)
    g.react("partition", "Partition cells with barcoded gel beads", StepKind.PARTITION,
            ((tagged, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)), (occupied, failed), changes_molecule=False)
    g.discard(failed, "filtered computationally; not a separate library")
    rna_bead = g.input("rna_bead_primer", "Multiome RNA gel-bead primer", MaterialKind.OLIGO, state=rna_bead_primer(), source=GUIDE)
    atac_bead = g.input("atac_bead_primer", "Multiome ATAC gel-bead primer", MaterialKind.OLIGO, state=atac_bead_primer(), source=GUIDE)
    rt = g.input("rt", "RT reagents and barcoding enzyme in master mix", MaterialKind.CATALYST)
    pool = molecule("barcoded DOGMA pool", (role("cell-barcoded cDNA"), role("cell-barcoded ATAC DNA"),
                                            role("cell-barcoded TSA-tag cDNA")))
    barcoded = g.product("barcoded", "barcoded cDNA, ATAC DNA and poly(dT)-captured tag cDNA", MaterialKind.MOLECULAR, state=pool)
    g.react("gem_incubation", "GEM incubation: reverse transcription and ATAC barcoding", StepKind.REVERSE_TRANSCRIPTION,
            ((occupied, InputRole.TEMPLATE), (rna_bead, InputRole.STRUCTURAL), (atac_bead, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (barcoded,))
    silane = g.input("silane", "post-GEM silane bead cleanup", MaterialKind.SELECTION)
    recovered = g.product("recovered", "recovered barcoded products", MaterialKind.MOLECULAR, state=pool)
    spent = g.product("silane_waste", "leftover reagents and primers", MaterialKind.WASTE)
    g.react("silane_cleanup", "Break GEMs and purify with silane beads", StepKind.CLEANUP,
            ((barcoded, InputRole.TEMPLATE), (silane, InputRole.SELECTION)), (recovered, spent), changes_molecule=False)
    g.discard(spent, "removed in the silane cleanup")
    preamp_primers = g.input("preamp_primers", "10x pre-amplification primers", MaterialKind.OLIGO, inferred=True,
                             source=GUIDE + " (commercial primer sequences not printed)")
    additive = g.input("adt_additive", "1 µl 0.2 µM ADT additive CCTTGGCACCCGAGAATT*C*C", MaterialKind.OLIGO,
                       state=Construct([Segment("ADT additive (two 3′ phosphorothioates)", ADT_ADDITIVE)]),
                       source=SOURCE + " (Methods step 4.1 modification)")
    preamp_pol = g.input("preamp_polymerase", "pre-amplification mix", MaterialKind.CATALYST)
    preamplified = g.product("preamplified", "pre-amplified pool eluted in 100 µl EB", MaterialKind.MOLECULAR, state=pool)
    g.react("preamp", "Pre-amplify with spiked ADT additive, SPRI, elute 100 µl", StepKind.PCR,
            ((recovered, InputRole.TEMPLATE), (preamp_primers, InputRole.STRUCTURAL), (additive, InputRole.STRUCTURAL),
             (preamp_pol, InputRole.CATALYST)), (preamplified,),
            geometry=PCRGeometry.unavailable(pool, reason="Commercial pre-amplification primers are not printed; the pool has several precursor types."))
    atac = g.product("atac_branch", "25% for ATAC library indexing", MaterialKind.MOLECULAR, state=atac_precursor())
    cdna_in = g.product("cdna_aliquot", "35% for cDNA amplification", MaterialKind.MOLECULAR, state=pool)
    adt = g.product("adt_branch", "35% for tag indexing", MaterialKind.MOLECULAR, state=adt_precursor())
    g.react("allocate", "Allocate pre-amplified sample to three reactions", StepKind.SPLIT,
            ((preamplified, InputRole.TEMPLATE),), (atac, cdna_in, adt), changes_molecule=False)
    reviewed = ["cells", "antibodies", "stain", "lysis", "permeabilize", "transposome", "tn5", "tagment", "droplets",
                "partition", "rna_bead_primer", "atac_bead_primer", "rt", "gem_incubation", "silane", "silane_cleanup",
                "preamp_primers", "adt_additive", "preamp_polymerase", "preamp", "allocate"]
    cdna_primers = g.input("cdna_primers", "10x cDNA primer mix, forward component", MaterialKind.OLIGO, inferred=True,
                           source=GUIDE + " (commercial primer sequences not printed)")
    cdna_reverse = g.input("cdna_reverse", "10x cDNA primer mix, reverse component", MaterialKind.OLIGO, inferred=True,
                           source=GUIDE + " (commercial primer sequences not printed)")
    cdna_pol = g.input("cdna_polymerase", "cDNA amplification mix", MaterialKind.CATALYST)
    cdna = g.product("cdna", "amplified transcript cDNA", MaterialKind.MOLECULAR,
                     state=molecule("GEX cDNA", (role("cell barcode"), role("UMI"), role("transcript cDNA"))))
    g.react("amplify_cdna", "Amplify cDNA", StepKind.PCR,
            ((cdna_in, InputRole.TEMPLATE), (cdna_primers, InputRole.STRUCTURAL), (cdna_reverse, InputRole.STRUCTURAL), (cdna_pol, InputRole.CATALYST)), (cdna,),
            geometry=PCRGeometry.unavailable(pool, reason="Commercial cDNA primer sequences are not printed in the inspected sources."))
    frag = g.input("rna_fragmentation", "fragmentation, end-repair and A-tailing enzymes", MaterialKind.CATALYST)
    tailed = g.product("rna_tailed", "fragmented, end-repaired, A-tailed cDNA", MaterialKind.MOLECULAR, state=cdna.state)
    g.react("rna_fragment", "Fragment, end-repair and A-tail cDNA", StepKind.FRAGMENTATION,
            ((cdna, InputRole.TEMPLATE), (frag, InputRole.CATALYST)), (tailed,))
    adapter = g.input("rna_adapter", "10x GEX adaptor", MaterialKind.ADAPTER)
    ligase = g.input("rna_ligase", "DNA ligase", MaterialKind.CATALYST)
    adapted = g.product("rna_adapted", "adaptor-ligated GEX fragments", MaterialKind.MOLECULAR, state=cdna.state)
    g.react("rna_ligate", "Ligate the GEX adaptor", StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (adapted,),
            geometry=AdapterLigationGeometry.unavailable("The guide prints the finished library; the adaptor's ligation-ready end geometry was not inspected."))
    reviewed.extend(("cdna_primers", "cdna_reverse", "cdna_polymerase", "amplify_cdna", "rna_fragmentation", "rna_fragment",
                     "rna_adapter", "rna_ligase", "rna_ligate"))
    _endpoint(g, "rna", adapted, RNA, RNA_P, reviewed,
              geometry=PCRGeometry.unavailable(adapted.state, reason="Dual-index TT primer sequences were not inspected; the guide prints the finished library."))
    _endpoint(g, "atac", atac, ATAC, ATAC_P, reviewed,
              geometry=PCRGeometry.unavailable(atac.state, reason="Multiome ATAC sample-index primer sequences were not inspected; the guide prints the finished library."))
    fwd = PCRPrimer("SI-PCR", SI_PCR, anneal_nt=13)
    rev = PCRPrimer("RPI-x", il.P7 + "N" * 8 + RPI_ARM, anneal_nt=len(TSA_HANDLE))
    _endpoint(g, "adt", adt, ADT, ADT_P, reviewed,
              forward_state=Construct([Segment("SI-PCR", SI_PCR)]),
              reverse_state=Construct([Segment("P7", il.P7, "p7"),
                                       Segment("i7 index", "N" * 8, "cbc", placeholder=True, feature=ADT_I7),
                                       Segment("RPI-x arm", RPI_ARM, "r2")]),
              geometry=PCRGeometry.exact(adt.state, fwd, rev,
                                         statement="Methods/MOESM1 print SI-PCR and RPI-x (KAPA, 10 cycles); precursor uses the printed Multiome RNA bead primer."),
              limitation=ADT_LIMITATION)
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
