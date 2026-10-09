"""Rohland et al. 2015 partial-UDG double-stranded ancient-DNA library."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from base_conversion import PartialUDGProduct
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, Scene, feature

TITLE = "Double-stranded ancient DNA with partial UDG"
NOTES = "01_partial-udg-ds-adna.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1098/rstb.2013.0624">Rohland et al., <i>Philosophical Transactions B</i> (2015)</a>.'
SUMMARY = "USER removes internal uracils before end repair, while terminal damage is retained as an authenticity signal; two molecular barcodes are ligated before indexed PCR."
CAVEAT = ("The article establishes the two 7-nt molecular barcodes and completed dual-index Illumina geometry, but the supplementary oligo table was not available through the open full text. "
          "Adapter bases are therefore shown as an inferred standard TruSeq-compatible endpoint.")
P5_BARCODE = "ATCGATT"  # representative source-printed barcode in Table 2
P7_BARCODE = "GACTTAT"  # paired representative source-printed barcode in Table 2


def library():
    def adapter(name, bases, tag=None, **kw):
        return seg(name, bases, tag, inferred=True, **kw)
    lib = Construct([
        adapter("P5", il.P5, "p5"),
        adapter("7-nt i5", "J"*7, "cbc", placeholder=True,
                feature=feature("sample_i5", "sample_index", "unknown")),
        adapter("Read 1 arm", il.TRUSEQ_READ1, "r1"),
        seg("P5-side 7-nt molecular barcode", P5_BARCODE, "cbc",
            feature=feature("molecular_barcode_p5", "umi", "whitelist", whitelist="published molecular barcode set", group="molecular_id", part="P5 side")),
        seg("partial-UDG ancient DNA insert", "X"*36, placeholder=True),
        seg("P7-side 7-nt molecular barcode", P7_BARCODE, "cbc",
            feature=feature("molecular_barcode_p7", "umi", "whitelist", whitelist="published molecular barcode set", group="molecular_id", part="P7 side")),
        adapter("Read 2 site start", "A"),
        adapter("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        adapter("7-nt i7 reverse complement", "I"*7, "cbc", placeholder=True,
                feature=feature("sample_i7", "sample_index", "unknown")),
        adapter("P7 reverse complement", il.P7_RC, "p7"),
    ], name="partial-UDG ancient-DNA library")
    primers=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["I2"],sp.TRUSEQ["R2"])
    problems=sp.verify(lib,primers)
    if problems: raise ValueError("invalid partial-UDG library: "+"; ".join(problems))
    return lib,primers


FINAL_LIBRARY, SEQ_PRIMERS = library()
FINAL_CAPTION = "INFERRED — completed dual-index Illumina endpoint around the two source-printed 7-nt molecular barcodes."
SEQUENCING_INTRO = "Read 1 and Read 2 each begin with a molecular barcode; i5 and i7 distinguish experiments after indexing PCR."


def partial_user_rows(product):
    damaged = [
        seg("left terminal U", "U", "w1"),
        seg("left surviving interior", "ACG"),
        seg("internal U cleavage site", "U", "w1"),
        seg("right surviving interior", "TCG"),
        seg("right terminal U", "U", "w1"),
    ]
    before = Scene()
    before.strand("damaged strand", damaged)
    before.mark("damaged strand", "internal U cleavage site", "USER cleavage")

    survivor = Scene()
    survivor.strand("surviving termini", [
        seg("retained left terminal U", product.left_terminal, "w1"),
        seg("surviving interior", product.interior),
        seg("retained right terminal U", product.right_terminal, "w1"),
    ])
    survivor.note("surviving termini", "internal U-containing pieces removed")
    return [*before.rows(), *survivor.rows()]


def sections():
    survivor=PartialUDGProduct("U","ACGTCG","U")
    lig=Scene.duplex(list(FINAL_LIBRARY),label="barcoded library")
    lig.junction("top","P5-side 7-nt molecular barcode","partial-UDG ancient DNA insert","adapter ligation")
    lig.junction("top","partial-UDG ancient DNA insert","P7-side 7-nt molecular barcode","adapter ligation")
    return [
        ("Partial USER treatment", partial_user_rows(survivor),
         "UDG and Endonuclease VIII act before T4 polymerase and kinase; terminal uracils remain inefficiently removed."),
        ("End repair and barcode-adapter ligation", lig.rows(),
         "T4 polymerase/PNK blunt the surviving fragments; distinct P5- and P7-side 7-mers are ligated. ** marks both ligations."),
        ("INFERRED — fill in and indexed PCR", [Row(chunks=[("short barcoded library → fill-in → i5/i7 PCR → sequencing library",None,False)])],
         "The paper explicitly describes this order; unavailable supplementary oligo bases are not presented as source-transcribed."),
    ]
