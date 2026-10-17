"""Protocol-owned spatial assay and reaction graph."""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from protocol_families.spatial import (ident, UMI, CELL_A, SAMPLE, SECTION_HASH,
    SPACE_HASH)
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Segment, feature
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState, Scene, Workflow, duplex_rows
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
from protocol_model import chain, molecule, role, role_library

from types import SimpleNamespace

TITLE = 'sci-Space'
NOTES = '01_sci_space.html'
SOURCE = 'Defining source: <a href="https://doi.org/10.1126/science.abb9536">Srivatsan et al., Science (2021)</a>.'
SUMMARY = 'Polyadenylated DNA hashes diffuse from a registered printed grid into permeabilized tissue nuclei; sci-RNA-seq then associates each grid identity with a single-nucleus transcriptome.'
CAVEAT = ('The printed hash, section hash and transcript are separate molecules that are associated through shared nucleus indices; the final construct shows the spatial-hash branch, not a fictitious covalent fusion to a transcript. The readable sources do not establish the section-hash tail sequence, the barcoded oligo(dT) RT primer or the index members. '
          '<b>Source differences.</b> The published Supplementary Materials (Materials and Methods, pp. 3-5) differ from the linked protocols.io guides on several conditions; the model follows the study-specific Supplementary Materials and the differences remain visible here: '
          'fixation of the scraped embryo section in 4% paraformaldehyde (SM p. 3, "Spatial transfer") versus 5% PFA in the transfer guide (steps 8-10) and in the SM hash-concentration test (p. 2); '
          'indexed oligo(dT) at 25 uM (SM p. 4) versus 10 uM barcoded RT primer (two-level guide step 15); an RT ramp starting at 10 C (SM) versus 4 C (guide); '
          'index PCR 75 C 3 min then 18 cycles with 1-min extension (SM p. 5) versus 72 C 5 min then 18-22 cycles (guide step 53); '
          'pooled Zymo DNA Clean & Concentrator plus 0.8x AMPure XP (SM p. 5) versus pooled 1x AMPure (guide step 55); and only a custom n7-loaded Tn5 in the SM, whereas the guide also allows dual-loaded Illumina Tn5. '
          'These are differences between a general guide and the study methods, not established publication errors.')
FINAL_CAPTION = 'Separate transcript, printed-grid-hash and section-hash molecules share the same combinatorial nucleus indices; the representative final construct below is a spatial-hash molecule.'
SEQUENCING_INTRO = ('Separate transcript, printed-grid-hash and section-hash molecules share the same combinatorial nucleus indices; the representative final construct below is a spatial-hash molecule. '
                    'The Supplementary Materials (p. 5) state a NextSeq 500 high-output 75-cycle run with onboard primers: Read 1 18 cycles (8-nt UMI and 10-nt RT barcode), Index 1 10, Index 2 10 and Read 2 52 cycles (transcript or hash-oligo). '
                    'Hash reads are assigned after poly(A) trimming; the UMI and RT-barcode lengths are stated there but the RT-primer sequence is not printed in the inspected sources.')
READ_LENGTHS = {"Read 1": 18, "Index 1 (i7)": 10, "Index 2 (i5)": 10, "Read 2": 52}
READ_USES = {"Read 1": "extract", "Index 1 (i7)": "demultiplex",
             "Index 2 (i5)": "demultiplex", "Read 2": "include"}
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
GRID_HASH_WRITTEN = ('GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
                     '[10bp-barcode sequence]-B'
                     'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')


def oligos():
    return (
        'Spot/sector hash, 5′-' + GRID_HASH_WRITTEN + '-3′ '
        '(B = C, G or T; grid-preparation protocol step 7.1). '
        'The attached ordering workbook is listed by that protocol but is not in the fetched PDF.',
        'Indexed P5 PCR primer, 5′-AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′ '
        '(two-level sci-RNA-seq guide step 51; the 3′ TCT continues after a page break).',
        'Indexed P7 PCR primer, 5′-CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG-3′ (same step).',
    )
GUIDE = "two-level sci-RNA-seq guide, doi:10.17504/protocols.io.bquamwse"


def sci_space_library():
    """Representative copied printed-hash molecule (local; not a transcript fusion).

    The grid guide step 7.1 hash is 5'-S7-ME-[barcode]-B-A32-3'.  Copying it from an
    oligo(dT) RT primer gives [RT primer]-oligo(dT)-[barcode]'-ME'-S7', which carries the
    P7-primer site without tagmentation (derived).  The two-level guide (step 51) prints
    the P5 primer with the full TruSeq Read 1 arm; the RT primer itself is not printed in
    the read guides, so UMI, RT barcode and oligo(dT) stay role tokens.
    """
    def s(name, top, tag=None, **kw):
        return Segment(name, top, tag, **kw)
    grid = feature("sci_space.grid", "spatial_barcode", "whitelist",
                   whitelist="sci-Space printed coordinate-hash panel",
                   note="printed coordinate-hash identity registered by fluorescent fiducials")
    cell = feature("sci_space.rt", "cell_barcode", "combinatorial",
                   group="cell_barcode", part="RT round")
    umi = feature("sci_space.umi", "umi", "random")
    i5 = feature("sci_space.i5", "sample_index", "unknown")
    i7 = feature("sci_space.i7", "sample_index", "unknown")
    lib = Construct([
        s("P5", il.P5, "p5"),
        s("i5", "I" * 10, "cbc", placeholder=True, feature=i5),
        s("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
        role("UMI", tag="umi", molecular_feature=umi),
        role("RT-well cell barcode", tag="cbc", molecular_feature=cell),
        role("oligo(dT) copy of the hash poly(A)"),
        role("copied printed spatial hash and B complement", tag="cbc", molecular_feature=grid),
        s("mosaic end reverse complement", nx.ME_RC, "me"),
        s("S7 reverse complement", nx.S7_RC, "s7"),
        s("i7 reverse complement", "J" * 10, "cbc", placeholder=True, feature=i7),
        s("P7 reverse complement", il.P7_RC, "p7"),
    ], name="representative sci-Space spatial-hash library molecule")
    primers = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"],
               sp.TRUSEQ["I2"], sp.NEXTERA["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("sci-Space library: " + "; ".join(problems))
    return lib, primers


FINAL_LIBRARY, SEQ_PRIMERS = sci_space_library()
_definition = SimpleNamespace(title=TITLE, doi='10.1126/science.abb9536', citation='Srivatsan et al., Science (2021)', endpoint='Separate transcript, printed-grid-hash and section-hash molecules share the same combinatorial nucleus indices; the representative final construct below is a spatial-hash molecule.')


def _sci_space():
    return chain(molecule("fresh-frozen tissue section", (
        role("poly(A) RNA in intact nuclei"), role("anatomical position"))), (
        ("Permeabilize the section and add one polyadenylated section hash",
         molecule("section-labelled nuclei in tissue", (SECTION_HASH,
             role("poly(A) RNA in intact nuclei"))),
         "The section hash identifies the slide after sections are pooled; it is distinct from the printed coordinate hashes."),
        ("Juxtapose tissue with the agarose slide carrying a gridded hash pattern",
         molecule("spatially hashed nuclei", (SPACE_HASH, SECTION_HASH,
             role("poly(A) RNA in intact nuclei"))),
         "Single-stranded poly(A)-tailed hash oligos diffuse from each printed spot into nearby permeabilized nuclei. Fluorescent fiducials register the grid to the tissue image."),
        ("Recover nuclei, fix them and perform barcoded sci-RNA-seq reverse transcription",
         molecule("first-round indexed nuclear molecules", (
             ident("RT-well cell barcode", "sci_space.rt", "cell_barcode"),
             UMI, role("transcript cDNA or copied hash"))),
         "The same oligo(dT) reaction copies endogenous mRNA, the grid hash and the section hash, associating all three with a nucleus-specific RT-well barcode."),
        ("Pool, redistribute, second-strand synthesize, tagment and PCR-index",
         molecule("sci-Space sequencing molecules", (
             ident("RT-well cell barcode", "sci_space.rt", "cell_barcode"),
             UMI, SPACE_HASH, role("transcript or hash payload"), SAMPLE)),
         "The downstream combinatorial sci-RNA-seq indices identify one nucleus. Hash molecules and transcript molecules are separate physical reads joined computationally by that shared nucleus identity."),
    ))

workflow = _sci_space

def _sci_space_graph(d, final, primers):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    tissue = graph.input("tissue", "fresh-frozen tissue section containing intact nuclei",
                         MaterialKind.BIOLOGICAL, state=_sci_space().initial, primary=True)
    section_hash = graph.input("section_hash",
                               "slide-specific hash oligo, 5 uL of 10 uM in 500 uL permeabilization buffer with DAPI (SM p. 3; transfer guide step 1.3; sequence and tail unverified)",
                               MaterialKind.OLIGO)
    lysis = graph.input("cold_lysis_buffer", "permeabilization buffer: 10 mM Tris-HCl pH 7.4, 10 mM NaCl, 3 mM MgCl2, 1% SUPERase-In, 0.1% IGEPAL CA-630 (SM p. 3)",
                        MaterialKind.CONDITION)
    sectioned = graph.product("sectioned", "section-labelled tissue nuclei",
                              MaterialKind.MOLECULAR,
                              state=molecule("section-labelled nuclei", (SECTION_HASH,
                                  role("poly(A) RNA"))))
    graph.react("label_section", "Permeabilize the section in lysis buffer carrying the slide oligo",
                StepKind.PERMEABILIZATION,
                ((tissue, InputRole.TEMPLATE), (lysis, InputRole.CONDITION),
                 (section_hash, InputRole.STRUCTURAL)),
                (sectioned,), changes_molecule=False)
    grid_hash = graph.input("grid_hash",
                            "registered printed grid of poly(A)-tailed hash oligos, 15 uM spot + 2.5 uM sector oligo per printed well, 0.5% glycerol, SYBR Green at fiducials (SM p. 1; grid guide step 7.1)",
                            MaterialKind.OLIGO)
    spatial = graph.product("spatial_nuclei", "nuclei carrying section and grid hashes",
                            MaterialKind.MOLECULAR,
                            state=molecule("spatially hashed nuclei", (SPACE_HASH,
                                SECTION_HASH, role("poly(A) RNA"))))
    # Article/transfer guide describe diffusion of single-stranded hashes into nearby
    # permeabilized nuclei; no hybridization partner is named, so no pairing is claimed.
    graph.react("transfer_grid", "Juxtapose the grid; printed hashes diffuse into nearby nuclei",
                StepKind.GENERIC,
                ((sectioned, InputRole.TEMPLATE), (grid_hash, InputRole.STRUCTURAL)),
                (spatial,), changes_molecule=False)
    pfa = graph.input("pfa", "4% PFA (SM p. 3; the transfer guide uses 5%), sonication for slides 7-14, 15 min on ice; pooled in NSB and frozen",
                      MaterialKind.CONDITION)
    fixed = graph.product("fixed_nuclei", "scraped, fixed hashed nuclei", MaterialKind.MOLECULAR,
                          state=spatial.state)
    graph.react("fix_nuclei", "Scrape the section into PFA and fix", StepKind.CROSS_LINKING,
                ((spatial, InputRole.TEMPLATE), (pfa, InputRole.CONDITION)), (fixed,),
                changes_molecule=False)
    perm = graph.input("thaw_permeabilize",
                       "thaw, NSB + 0.25% Triton X-100 3 min, NSB resuspension, 12 s Bioruptor sonication; 5000 nuclei per RT well (SM pp. 3-4)",
                       MaterialKind.CONDITION)
    thawed = graph.product("thawed_nuclei", "permeabilized fixed hashed nuclei, 5000 per RT well",
                           MaterialKind.MOLECULAR, state=spatial.state)
    graph.react("permeabilize_nuclei", "Thaw, permeabilize, sonicate and distribute nuclei",
                StepKind.PERMEABILIZATION,
                ((fixed, InputRole.TEMPLATE), (perm, InputRole.CONDITION)), (thawed,),
                changes_molecule=False)
    rt_primer = graph.input("rt_primer",
                            "uniquely indexed oligo(dT), 1 uL of 25 uM per well, annealed 55 C 5 min (SM p. 4; guide step 15 uses 10 uM; sequence not printed)",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase",
                     "SuperScript IV, DTT, RNaseOUT; 10/20/30/40/50 C 2 min each then 55 C 15 min (SM p. 4; the guide ramp starts at 4 C)",
                     MaterialKind.CATALYST)
    first = graph.product("first_strand", "RT-well-indexed transcript and hash cDNA",
                          MaterialKind.MOLECULAR,
                          state=molecule("first-round indexed nuclear molecules", (
                              CELL_A, UMI, role("transcript cDNA or copied hash"))))
    graph.react("reverse_transcribe", "Copy endogenous RNA and both hash classes",
                StepKind.REVERSE_TRANSCRIPTION,
                ((thawed, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    # SM p. 4 names only NEB mRNA second strand synthesis buffer and enzyme;
    # no separate second-strand primer is printed, so none is declared.
    polymerase = graph.input("second_strand_enzyme",
                             "NEB mRNA second strand synthesis buffer and enzyme after pooling and FACS sorting (200/well, slides 1-6) or dilution (250/well, slides 7-14); 16 C, 150 min (SM p. 4)",
                             MaterialKind.CATALYST)
    duplex = graph.product("duplex", "double-stranded indexed transcript/hash products",
                           MaterialKind.MOLECULAR, state=first.state)
    graph.react("second_strand", "Synthesize the second strand after pool redistribution",
                StepKind.EXTENSION,
                ((first, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST)), (duplex,))
    transposome = graph.input("transposome", "custom n7-loaded Tn5, 0.01 uL per well (SM p. 4; adapter sequence not printed there)",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 in 5.74 uL 2x Nextera TD buffer; 55 C, 5 min (SM p. 4)",
                      MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented sci-Space products", MaterialKind.MOLECULAR,
                           state=duplex.state)
    waste = graph.product("tn5_waste", "candidate fragments lacking an amplifiable end pair",
                          MaterialKind.WASTE)
    graph.react("tagment", "Tagment transcript and hash products", StepKind.TAGMENTATION,
                ((duplex, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.mark_terminal(Terminal(waste.id, Disposition.UNRESOLVED,
        "No separation is described; such fragments co-purify on SPRI beads and enter PCR."))
    beads = graph.input("spri_tagment", "12 uL Zymo DNA binding buffer 5 min RT, then 36 uL AMPure XP, eluted in 17 uL EB (SM p. 4)",
                        MaterialKind.SELECTION)
    eluted = graph.product("eluted", "SPRI-purified tagmented products", MaterialKind.MOLECULAR,
                           state=duplex.state)
    sup = graph.product("spri_supernatant", "SPRI supernatant and ethanol washes", MaterialKind.WASTE)
    graph.react("tagment_cleanup", "Stop tagmentation and purify each well", StepKind.CLEANUP,
                ((tagged, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (eluted, sup))
    graph.discard(sup, "supernatant and washes removed from the beads")
    forward = graph.input("library_forward",
                          "indexed P5 primer, 5'-AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT-3' (guide step 51)",
                          MaterialKind.OLIGO)
    reverse = graph.input("library_reverse",
                          "indexed P7 primer, 5'-CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG-3' (guide step 51)",
                          MaterialKind.OLIGO)
    pcr_pol = graph.input("library_polymerase",
                          "NEBNext High-Fidelity 2x PCR Master Mix; 75 C 3 min, 98 C 30 s, 18 cycles of 98 C 10 s/66 C 30 s/72 C 1 min, 72 C 5 min (SM p. 5)",
                          MaterialKind.CATALYST)
    amplified = graph.product("amplified", "per-well indexed PCR products", MaterialKind.MOLECULAR,
                              state=final)
    graph.react("library_pcr", "Amplify each well with a unique P5/P7 index pair", StepKind.PCR,
                ((eluted, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (pcr_pol, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    final, reason="The two-level guide prints the P5/P7 architectures; the RT primer, Tn5 adapter, index members and section-hash tail are not printed in the inspected sources, so exact templates are not drawn."))
    final_beads = graph.input("ampure_final", "pooled PCR products, Zymo DNA Clean & Concentrator then 0.8x AMPure XP (SM p. 5)",
                              MaterialKind.SELECTION)
    library = graph.product("library_library", "sci-Space sequencing library",
                            MaterialKind.LIBRARY, state=final)
    final_sup = graph.product("ampure_supernatant", "AMPure supernatant and washes",
                              MaterialKind.WASTE)
    graph.react("final_cleanup", "Pool wells, concentrate and purify with 0.8x AMPure XP", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (final_beads, InputRole.SELECTION)),
                (library, final_sup))
    graph.discard(final_sup, "supernatant and washes removed from the beads")
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers),
                                           read_lengths=READ_LENGTHS, read_uses=READ_USES))
    reviewed = ["tissue", "section_hash", "cold_lysis_buffer", "label_section", "grid_hash",
                "transfer_grid", "pfa", "fix_nuclei", "thaw_permeabilize", "permeabilize_nuclei",
                "rt_primer", "reverse_transcriptase", "reverse_transcribe",
                "second_strand_enzyme", "second_strand",
                "transposome", "tn5", "tagment", "spri_tagment", "tagment_cleanup",
                "library_forward", "library_reverse", "library_polymerase", "library_pcr",
                "ampure_final", "final_cleanup", "library_library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _sci_space_graph(_definition, FINAL_LIBRARY, SEQ_PRIMERS)
