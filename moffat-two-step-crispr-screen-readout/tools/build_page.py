#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import pooled_crispr as M
import pooled_crispr_page as P
import seqprimers as sp
from chemdraw import Construct, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

OUT = HERE.parent / "moffat-two-step-crispr-screen-readout.html"
SOURCE = 'Source: Moffat Lab <a href="https://media.addgene.org/cms/filer_public/44/9f/449f22c6-f818-4af6-8af4-2bab1d20f00a/moffat_lab_crispr_screen_protocol_v2.pdf">Pooled CRISPR Screen Protocol v2</a> (2019), LCV2::TKOv3 branch.'


# Moffat Lab protocol rev. 2019-07-23, Table 1 (PDF text lines 865-941): LCV2 PCR1 pair and
# the representative PCR2 members drawn here.  Red/blue colouring in the PDF marks the
# index and annealing parts; the printed oligos carry no modifications ("desalted").
LCV2_PCR1_F = "GAGGGCCTATTTCCCATGATTC"
LCV2_PCR1_R = "GTTGCGAAAAAGAACGTTCACGG"
D501_F = ("AATGATACGGCGACCACCGAGATCTACACTATAGCCTACACTCTTT"
          "CCCTACACGACGCTCTTCCGATCTTTGTGGAAAGGACGAAACACCG")
D701_R = ("CAAGCAGAAGACGGCATACGAGATCGAGTAATGTGACTGGAGTTCA"
          "GACGTGTGCTCTTCCGATCTACTTGCTATTTCTAGCTCTAAAAC")
PCR2_F_ANNEAL = D501_F[-22:]   # blue annealing part of every D50x-F primer
PCR2_R_ANNEAL = D701_R[-24:]   # blue annealing part of every D70x-R primer


def pcr1_product() -> Construct:
    """PCR1 product with the nested PCR2 annealing sites it must contain.

    Flanks are placeholders: the protocol gives only the ~600-bp size, not the LCV2
    sequence between the PCR1 and PCR2 sites.
    """
    return Construct([
        M.seg("PCR1 LCV2 forward site", LCV2_PCR1_F, "s5"),
        M.seg("left LCV2 flank (not printed)", "X" * 12, "w1", placeholder=True),
        M.seg("hU6 3' end / PCR2 forward site (incl. +1 G)", PCR2_F_ANNEAL, "s5"),
        M.seg("20-nt guide", M.SPACER, "cbc", placeholder=True,
              feature=feature("guide", "guide_barcode", "unknown")),
        M.seg("scaffold / PCR2 reverse site", revcomp(PCR2_R_ANNEAL), "r2"),
        M.seg("right LCV2 flank (not printed)", "X" * 12, "w1", placeholder=True),
        M.seg("PCR1 LCV2 reverse site", revcomp(LCV2_PCR1_R), "r2"),
    ], name="PCR1 enriched guide locus (~600 bp)")


def protocol_graph() -> ProtocolGraph:
    """Section 3.5.3: PCR1 enrichment, PCR2 indexing, 2% gel excision of the 200-bp band."""
    source = SOURCE
    graph = ProtocolGraph("Moffat two-step CRISPR screen readout", authority=source)
    gdna = graph.input("genomic_dna", "genomic DNA from the LCV2::TKOv3 screen (100 µg per sample)",
                       MaterialKind.BIOLOGICAL,
                       state=molecule("screen genomic DNA", (role("integrated guide cassette"),)),
                       primary=True)
    first_forward = graph.input("pcr1_forward", "PCR1 LCV2 forward primer (desalted)",
                                MaterialKind.OLIGO,
                                state=(M.seg("LCV2 forward", LCV2_PCR1_F, "s5"),))
    first_reverse = graph.input("pcr1_reverse", "PCR1 LCV2 reverse primer (desalted)",
                                MaterialKind.OLIGO,
                                state=(M.seg("LCV2 reverse", LCV2_PCR1_R, "r2"),))
    polymerase = graph.input("q5_master_mix", "2x NEBNext Ultra II Q5 Master Mix",
                             MaterialKind.CATALYST)
    intermediate = pcr1_product()
    first = graph.product("pcr1_product", "pooled ~600-bp PCR1 guide-locus product",
                          MaterialKind.MOLECULAR, state=intermediate)
    graph.react("pcr1", "PCR1: 3.5 µg gDNA per 50 µL, 25 cycles (98 °C 10 s / 66 °C 30 s / "
                "72 °C 15 s); all reactions of a sample pooled", StepKind.PCR,
                ((gdna, InputRole.TEMPLATE), (first_forward, InputRole.STRUCTURAL),
                 (first_reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (first,), geometry=PCRGeometry.exact(
                    intermediate, PCRPrimer("LCV2 forward", LCV2_PCR1_F),
                    PCRPrimer("LCV2 reverse", LCV2_PCR1_R),
                    statement="Table 1 prints both PCR1 primers; the LCV2 sequence between "
                              "them is not given, so the flanks are placeholders."))
    second_forward = graph.input("pcr2_forward", "PCR2 i5 forward primer D501-F (desalted Ultramer)",
                                 MaterialKind.OLIGO, state=(M.seg("D501-F", D501_F, "p5"),))
    second_reverse = graph.input("pcr2_reverse", "PCR2 i7 reverse primer D701-R (desalted Ultramer)",
                                 MaterialKind.OLIGO, state=(M.seg("D701-R", D701_R, "p7"),))
    final = M.two_step_library()
    amplified = graph.product("pcr2_product", "PCR2 dual-index product (~200 bp) with "
                              "side products", MaterialKind.MOLECULAR, state=final)
    graph.react("pcr2_pcr", "PCR2: 5 µL pooled PCR1, 10 cycles (98 °C 10 s / 55 °C 30 s / "
                "65 °C 15 s)", StepKind.PCR,
                ((first, InputRole.TEMPLATE), (second_forward, InputRole.STRUCTURAL),
                 (second_reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.exact(
                    final, PCRPrimer("D501-F", D501_F), PCRPrimer("D701-R", D701_R),
                    statement="Table 1 prints the full D50x-F/D70x-R Ultramers; D501/D701 "
                              "is one representative combination."))
    gel = graph.input("gel_extraction", "2% agarose gel, Dark Reader excision of the 200-bp band, "
                      "gel extraction kit", MaterialKind.SELECTION)
    library = graph.product("pcr2_library", "gel-purified dual-index sequencing library",
                            MaterialKind.LIBRARY, state=final)
    waste = graph.product("gel_waste", "remaining gel: primers and off-size products",
                          MaterialKind.WASTE)
    graph.react("gel_purify", "Run all 50 µL on a 2% gel; excise and extract the 200-bp band",
                StepKind.CLEANUP, ((amplified, InputRole.TEMPLATE), (gel, InputRole.SELECTION)),
                (library, waste), changes_molecule=False)
    graph.discard(waste, "gel outside the excised 200-bp band is discarded")
    graph.sequence(library, SequencingPlan(library=final,
                                           primers=tuple(M.two_step_sequencing_primers())))
    graph.declare_reviewed_scope(
        ("genomic_dna", "pcr1_forward", "pcr1_reverse", "q5_master_mix", "pcr1",
         "pcr1_product", "pcr2_forward", "pcr2_reverse", "pcr2_pcr", "pcr2_product",
         "gel_extraction", "gel_purify", "gel_waste", "pcr2_library"), source)
    return graph.require_publishable()


def render() -> str:
    source = SOURCE
    intermediate, final = pcr1_product(), M.two_step_library()
    body = f'''<h2>PCR1 — enrich the integrated guide locus</h2>
{P.construct_panel(intermediate, "PCR1 uses the desalted LCV2 locus primers (Table 1) and produces the ~600-bp enrichment product. No flow-cell adapters are present. Flanks are placeholders: the protocol does not print the LCV2 sequence between the PCR1 and nested PCR2 sites.")}
<h2>PCR2 — add dual-index TruSeq arms</h2>
{P.construct_panel(final, "PCR2 uses the nested hU6 and scaffold annealing parts of the D50x-F/D70x-R Ultramers (D501/D701 drawn) and yields the ~200-bp library, which is excised from a 2% gel.")}
{sp.section(final, M.two_step_sequencing_primers(), heading="Sequencing primers", intro="Section 3.5.4 recommends single-read sequencing (HiSeq 2500 SR50 or NextSeq 500 SR75 with dual index): 21 dark cycles through the constant hU6 primer region, then 26 Read 1 cycles, plus 8-cycle i7 and i5 reads. Without dark cycles, the protocol instead adds a 0–5-nt stagger to the i5 primer and lengthens Read 1; that variant is not drawn.", required_roles=("Read 1", "Index 1 (i7)", "Index 2 (i5)"))}'''
    return P.page("Moffat two-step CRISPR screen readout", "01_moffat-two-step.html",
                  source, body, graph=protocol_graph())


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
