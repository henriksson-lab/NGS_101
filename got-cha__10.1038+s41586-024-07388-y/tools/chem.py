from types import SimpleNamespace
from chemdraw import Scene, Construct, Segment, feature
import illumina
import seqprimers as sp
from multimodal_protocols import Definition
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role
from rna_cell_assays import single

# Published Methods scATAC-seq with GoT-ChA: the 16-cycle i5 reports
# the 10x cell barcode, whereas the eight-base RPI-X window identifies a sample.
CELL = role("10x scATAC cell barcode", length_bp=16,
            molecular_feature=feature("cell_barcode", "cell_barcode", "unknown"))
SAMPLE = role("RPI-X sample index", length_bp=8,
              molecular_feature=feature("sample_index", "sample_index", "unknown"))
READ_LENGTHS = {"Read 1": 50, "Index 1 (i7)": 8, "Index 2 (i5)": 16, "Read 2": 50}
ACQUIRED_ROLES = tuple(READ_LENGTHS)

# Published version (Izzo, Myers et al., Nature 2024; PMC11139586) supersedes the
# bioRxiv preprint 10.1101/2022.05.11.491515 for every value below.
DEFINITION = Definition("GoT-ChA — genotyping of targeted loci with chromatin accessibility", "10.1038/s41586-024-07388-y",
        "Izzo, Myers et al., <i>Nature</i> (2024)", "Two locus-specific GoT-ChA primers (GoT-ChA_R1N with a Read 1N handle, and GoT-ChA_Rev) join the 10x scATAC in-droplet barcoding PCR, exponentially amplifying the target while ATAC fragments are linearly barcoded; a biotinylated hemi-nested PCR, streptavidin pull-down and on-bead indexing then enrich the variant without losing the ATAC cell barcode.",
        "NovaSeq 6000 paired-end: Read 1N 50 cycles, i7 index 8 cycles (sample index), i5 index 16 cycles (10x cell barcode), Read 2N 50 cycles. Primers place the variant within 50 bp of the Read 1N or Read 2 handle, linking genotype to the parent accessibility profile by cell barcode.", "01_got_cha.html")
P5 = illumina.P5  # Methods line 215 prints the 29-nt P5 primer.
RPI_X_3PRIME = "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"  # Methods line 215 / Supp. Table 1 row 33.
CAVEAT = ("INFERRED — final endpoint role order; exact kit and target-panel junctions are not reconstructed. Values follow the published article (Methods “scATAC-seq with GoT-ChA”, Supplementary Table 1); the 2022 preprint is superseded. "
          "Published discrepancies: (1) Methods (single-target) describe GoT-ChA_Rev as a plain 20–22-bp locus-specific primer and add the IS2 handle only in the multiplex section, whereas Supplementary Table 1 lists every GoT-ChA_Rev with the IS2 handle AGCAAGTGAGAAGCATCGTGTC; the drawn single-target route leaves that 5′ handle unassigned. "
          "(2) Methods print P5 as AATGATACGGCGACCACCGAGATCTACAC (29 nt); Supplementary Table 1 “Common Primers” lists p5 as AATGATACGGCGACCACCGAGA (22 nt). The model uses the Methods sequence. "
          "The polymerase for the hemi-nested and indexing PCRs is not named in that passage. The multiplex-adapted route (P5/IS2 pre-amplification, 8 cycles; hemi-nested PCR 8 cycles) and ASAP-seq integration are not drawn.")


def population_state(label, *populations):
    """Each supplied tuple is one molecule/population, not a shared covalent chain."""
    sc=Scene()
    for index,parts in enumerate(populations):
        sc.strand(f"population {index+1}",parts,label=f"{label}: population {index+1}")
    return sc


def _gotcha_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cells", "single-cell ATAC reaction containing a target genomic locus",
                        MaterialKind.BIOLOGICAL, state=_gotcha().initial, primary=True)
    bead_primer = graph.input("bead_primer", "10x gel-bead cell-barcode primer",
                              MaterialKind.OLIGO)
    locus_primer = graph.input("locus_primer",
                               "GoT-ChA_R1N primer (Read 1N handle + 20–22-bp locus sequence)",
                               MaterialKind.OLIGO)
    rev_primer = graph.input("locus_rev_primer",
                             "GoT-ChA_Rev primer (20–22-bp locus sequence; IS2 handle in Supp. Table 1)",
                             MaterialKind.OLIGO)
    polymerase = graph.input("extension_polymerase", "10x barcoding reaction polymerase",
                             MaterialKind.CATALYST)
    barcoded = graph.product("barcoded_products", "cell-barcoded ATAC and target-locus products",
                             MaterialKind.MOLECULAR,
                             state=population_state('barcoded ATAC and target molecules',(CELL,role("accessible-chromatin fragments"),),(CELL,role("target-locus amplicon"),)))
    graph.react("barcode_extension", "In-droplet barcoding PCR with 1 µL 22.5 µM GoT-ChA primer mix (exponential target, linear ATAC)",
                StepKind.EXTENSION,
                ((cells, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (locus_primer, InputRole.STRUCTURAL), (rev_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)),
                (barcoded,))
    atac = graph.product("atac_branch", "parent scATAC branch", MaterialKind.MOLECULAR,
                         state=barcoded.state)
    target = graph.product("target_branch", "target-enrichment aliquot", MaterialKind.MOLECULAR,
                           state=barcoded.state)
    graph.react("split_branches", "Post-GEM SPRIselect eluate (45.5 µL): 5 µL to GoT-ChA, 40 µL to ATAC",
                StepKind.SPLIT, ((barcoded, InputRole.TEMPLATE),), (atac, target),
                changes_molecule=False)
    reviewed = ["cells", "bead_primer", "locus_primer", "locus_rev_primer", "extension_polymerase",
                "barcode_extension", "split_branches"]
    atac_final = Construct([role("P5 / 10x scATAC end"), CELL,
                            role("accessible-chromatin insert"), SAMPLE,
                            role("P7 / 10x scATAC end")], name="GoT-ChA parent ATAC endpoint roles")
    atac_library = graph.product("atac_library", "GoT-ChA parent ATAC sequencing library",
                                 MaterialKind.LIBRARY, state=atac_final)
    graph.react("atac_library_prep", "Continue parent ATAC library construction per CG000209 Rev F",
                StepKind.MAGIC, ((atac, InputRole.TEMPLATE),), (atac_library,),
                reason="Published Methods assign the remaining 40 µL to the standard ATAC workflow; "
                       "CG000209 Rev F kit reactions and primers are not reconstructed in this graph.")
    graph.sequence(atac_library, SequencingPlan(read_lengths=READ_LENGTHS,
                   acquired_roles=ACQUIRED_ROLES, limitation="Parent scATAC branch: i5 identifies cells; "
                   "i7 identifies samples. Exact kit primer geometry remains unassigned."))
    reviewed.extend(("atac_library_prep", "atac_library"))
    nested_forward = graph.input("nested_forward", "P5 primer", MaterialKind.OLIGO,
                                 state=single("P5", (Segment("P5", P5),)))
    nested_reverse = graph.input(
        "nested_reverse",
        "GoT-ChA_nested: /5BiosG/CCTTGGCACCCGAGAATTCCA-[20–22-bp locus], partial TruSeq Small RNA Read 2 handle",
        MaterialKind.OLIGO)
    nested_polymerase = graph.input("nested_polymerase", "target PCR polymerase",
                                    MaterialKind.CATALYST)
    enriched = graph.product("enriched_target", "biotinylated variant-spanning genotype amplicon",
                             MaterialKind.MOLECULAR,
                             state=molecule("enriched genotype amplicon", (
                                 CELL, role("variant-spanning genomic insert"),
                                 role("biotinylated nested end"))))
    graph.react("nested_pcr", "Biotinylated hemi-nested PCR, P5 + GoT-ChA_nested, 15 cycles", StepKind.PCR,
                ((target, InputRole.TEMPLATE), (nested_forward, InputRole.STRUCTURAL),
                 (nested_reverse, InputRole.STRUCTURAL),
                 (nested_polymerase, InputRole.CATALYST)), (enriched,),
                geometry=PCRGeometry.unavailable(
                    target.state, reason="The locus-specific panel determines the target-primer geometry."))
    spri = graph.input("nested_spri", "1.2X SPRIselect", MaterialKind.CONDITION)
    purified = graph.product("purified_target", "SPRI-purified biotinylated amplicon",
                             MaterialKind.MOLECULAR, state=enriched.state)
    spri_waste = graph.product("nested_spri_waste", "material not retained by 1.2X SPRIselect",
                               MaterialKind.WASTE)
    graph.react("nested_spri_cleanup", "1.2X SPRIselect clean-up", StepKind.CLEANUP,
                ((enriched, InputRole.TEMPLATE), (spri, InputRole.CONDITION)),
                (purified, spri_waste), changes_molecule=False)
    graph.discard(spri_waste, "removed by 1.2X SPRIselect clean-up")
    beads = graph.input("streptavidin", "Dynabeads M-280 Streptavidin (SSPE washes, 15 min RT binding)",
                        MaterialKind.SELECTION)
    captured = graph.product("captured_target", "streptavidin-captured genotype amplicon",
                             MaterialKind.MOLECULAR, state=enriched.state)
    uncaptured = graph.product("uncaptured_target", "unbound material removed by SSPE and Tris washes",
                               MaterialKind.WASTE)
    graph.react("capture", "Bind biotinylated amplicons to M-280 streptavidin; wash 2× SSPE, 1× 10 mM Tris",
                StepKind.AFFINITY_CAPTURE,
                ((purified, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, uncaptured))
    graph.discard(uncaptured, "does not carry the selected biotinylated target amplicon")
    reviewed.extend(("nested_forward", "nested_reverse", "nested_polymerase",
                     "nested_pcr", "nested_spri", "nested_spri_cleanup", "purified_target",
                     "nested_spri_waste", "streptavidin", "capture"))
    index_forward = graph.input("genotype_forward", "P5 primer", MaterialKind.OLIGO,
                                state=single("P5", (Segment("P5", P5),)))
    index_reverse = graph.input(
        "genotype_reverse", "RPI-X indexing primer", MaterialKind.OLIGO,
        state=single("RPI-X", (
            Segment("P7", illumina.P7),
            Segment("sample index", "XXXXXXXX", "idx", placeholder=True,
                    feature=feature("sample_index", "sample_index", "unknown")),
            Segment("RPI-X 3′ handle", RPI_X_3PRIME))))
    index_polymerase = graph.input("genotype_polymerase", "indexing PCR polymerase",
                                   MaterialKind.CATALYST)
    library = graph.product("genotype_library", "GoT-ChA genotype sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("genotype_pcr", "On-bead indexing PCR, P5 + RPI-X, 6–10 cycles", StepKind.PCR,
                ((captured, InputRole.TEMPLATE), (index_forward, InputRole.STRUCTURAL),
                 (index_reverse, InputRole.STRUCTURAL), (index_polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    captured.state, reason="Locus-specific primer bases and the user-defined sample index vary with the selected target."))
    graph.sequence(library, SequencingPlan(read_lengths=READ_LENGTHS,
                   acquired_roles=ACQUIRED_ROLES, limitation=d.endpoint))
    reviewed.extend(("genotype_forward", "genotype_reverse", "genotype_polymerase",
                     "genotype_pcr", "genotype_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'got-cha'
    return _gotcha_graph(d, final)

def _gotcha():
    return chain(population_state('single-cell ATAC reaction',(role("tagmented accessible DNA"),),(role("target genomic locus"),)), (
        ("Add GoT-ChA_R1N and GoT-ChA_Rev to the in-droplet barcoding PCR",
         population_state('barcoded ATAC and target molecules',(CELL,role("ordinary accessible-chromatin fragments"),),(CELL,role("target-locus amplicon"),)),
         "GoT-ChA_R1N carries the Read 1N handle, so the gel-bead barcode primer can extend the exponentially amplified locus product with the same cell barcode as the linearly barcoded ATAC fragments."),
        ("Take 5 µL of the post-GEM eluate; hemi-nested PCR with P5 and biotinylated GoT-ChA_nested",
         molecule("enriched genotype amplicon", (CELL, role("variant-spanning genomic insert"),
                                                  role("biotinylated nested end"))),
         "15 cycles; GoT-ChA_nested binds at least 100 bp downstream of GoT-ChA_R1N in reverse orientation and adds a partial TruSeq Small RNA Read 2 handle. A 1.2X SPRIselect clean-up follows."),
        ("M-280 streptavidin pull-down, then on-bead P5/RPI-X indexing PCR (6–10 cycles)",
         molecule("GoT-ChA library", (role("P5 / Read 1 end"), CELL,
                                       role("variant-spanning insert"), SAMPLE,
                                       role("Read 2 / P7 end"))),
         "RPI-X adds the user-defined 8-nt sample index and P7; the genotype library is joined to the parent scATAC library computationally by cell barcode."),
    ))

d = DEFINITION
FINAL_LIBRARY = Construct([role("P5 / Read 1N kit end", inferred=True), CELL,
                           role("variant-spanning genomic insert"), SAMPLE,
                           role("TruSeq Small RNA Read 2 / P7 end", inferred=True)],
                          name="GoT-ChA genotype endpoint (role order inferred)")
SPEC = SimpleNamespace(
    TITLE=d.title,
    NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a> (published version of bioRxiv <a href="https://doi.org/10.1101/2022.05.11.491515">10.1101/2022.05.11.491515</a>).',
    SUMMARY=d.summary,
    CAVEAT=CAVEAT,
    FINAL_LIBRARY=FINAL_LIBRARY,
    SEQ_PRIMERS=(),
    FINAL_CAPTION=d.endpoint,
    SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=READ_LENGTHS,
    READ_SPECS=(
        sp.read_spec("Read 1", 50, "Read 1N target-locus sequence", "extract"),
        sp.read_spec("Index 1 (i7)", 8, "RPI-X sample index", "demultiplex"),
        sp.read_spec("Index 2 (i5)", 16, "10x cell barcode", "extract"),
        sp.read_spec("Read 2", 50, "variant-spanning sequence from Read 2N", "extract"),
    ),
    SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE="Target panels, pool barcodes and kit revisions vary; the defining source establishes the read roles, while this page does not invent a universal primer sequence.",
    workflow=_gotcha,
    oligos=lambda: (
        "<b>GoT-ChA_R1N</b> (Methods) — 5′-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG-[20–22-bp locus specific]-3′",
        "<b>GoT-ChA_Rev</b> (Supplementary Table 1) — 5′-AGCAAGTGAGAAGCATCGTGTC-[20–22-bp locus specific]-3′; Methods single-target text gives only the locus portion",
        "<b>GoT-ChA_nested</b> (Methods) — 5′-/5BiosG/CCTTGGCACCCGAGAATTCCA-[20–22-bp locus specific]-3′",
        "<b>P5</b> (Methods) — 5′-AATGATACGGCGACCACCGAGATCTACAC-3′",
        "<b>RPI-X</b> (Methods) — 5′-CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3′; X is the user-defined sample index",
        "<b>JAK2-V617 example</b> (Supplementary Table 1, IDT format) — GoT-ChA_R1 TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTGATGGCAGTTGCAGGTCC*A*T*A; GoT-ChA_Rev AGCAAGTGAGAAGCATCGTGTCACACTGACACCTAGCTGTGA*T*C*C; GoT-ChA_nested /5BiosG/CCTTGGCACCCGAGAATTCCAAAAGCCTGTAGTTTTACTTACTCTCG (* phosphorothioate)",
    ),
)
SPEC.protocol_graph = protocol_graph
globals().update(vars(SPEC))
