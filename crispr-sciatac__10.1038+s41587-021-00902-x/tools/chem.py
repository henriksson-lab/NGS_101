"""CRISPR-sciATAC guide and accessibility readouts."""
from chemdraw import Row
from multibranch import illumina_branch
from multimodal_spatial import finish_pcr_branch
from protocol_extensions import indexed_molecule
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
TITLE="CRISPR-sciATAC — pooled perturbations with accessibility readout"
NOTES="01_crispr-sciatac.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41587-021-00902-x">Liscovitch-Brauer et al., <i>Nature Biotechnology</i> (2021)</a>.'
SUMMARY="A lentiviral guide cassette and accessible chromatin are recovered from the same combinatorially indexed nucleus, linking each perturbation to its ATAC profile."
CELL_PARTS=(("transposition index",8),("PCR-well index",12))
ATAC,AP=illumina_branch("CRISPR-sciATAC chromatin","accessible genomic DNA",architecture="nextera",cell_parts=CELL_PARTS)
GUIDE,GP=illumina_branch("CRISPR-sciATAC guide","integrated guide cassette",cell_parts=CELL_PARTS,guide=20)
FINAL_LIBRARIES=(("Accessibility library",ATAC,AP,"Indexed ATAC fragments report open chromatin.","Nextera primers on the accessibility branch."),("Guide-identity library",GUIDE,GP,"Targeted amplification recovers the perturbation guide with the same cell-index combination.","Illumina primers on the guide branch."))
def sections():
    return [("Introduce the lentiviral guide cassette",[Row(chunks=[("provirus — U6 promoter — guide spacer — scaffold", "cbc", False)])],"The integrated cassette is the recoverable perturbation identity."),("Index accessible chromatin in successive pools",indexed_molecule(payload="tagmented accessible DNA",barcode_parts=(("transposition index",8),("PCR-well index",12),)).rows(),"The barcode combination identifies a nucleus."),("Recover guide identity with the same index combination",indexed_molecule(payload="guide-cassette amplicon",barcode_parts=(("transposition index",8),("PCR-well index",12),),target_barcode=20).rows(),"Targeted semi-suppressive PCR links the integrated guide to the cell identity."),]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    cells = g.input(
        "cells", "cells before lentiviral perturbation", MaterialKind.BIOLOGICAL,
        state=sections()[0][1], primary=True)
    guide_vector = g.input(
        "guide_vector", "lentiviral U6-guide cassette", MaterialKind.MOLECULAR)
    integration = g.input(
        "integration_machinery", "lentiviral integration machinery",
        MaterialKind.CATALYST)
    perturbed = g.product(
        "perturbed_cells", "cells carrying an integrated guide cassette",
        MaterialKind.BIOLOGICAL, state=sections()[0][1])
    g.react(
        "integrate_guide", "Introduce and integrate the lentiviral guide cassette",
        StepKind.INTEGRATION,
        ((cells, InputRole.TEMPLATE), (guide_vector, InputRole.STRUCTURAL),
         (integration, InputRole.CATALYST)), (perturbed,))

    transposome = g.input(
        "indexed_transposome", "well-indexed ATAC transposome",
        MaterialKind.ADAPTER)
    tn5 = g.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = g.product(
        "tagged_nuclei", "nuclei retaining indexed accessible DNA and guide identity",
        MaterialKind.MOLECULAR, state=sections()[1][1])
    tagment_waste = g.product(
        "tagment_waste", "non-library tagmentation products", MaterialKind.WASTE)
    g.react(
        "tagment", "Index accessible chromatin by transposition",
        StepKind.TAGMENTATION,
        ((perturbed, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
         (tn5, InputRole.CATALYST)), (tagged, tagment_waste))
    g.discard(tagment_waste,
              "does not retain an amplifiable indexed accessible-DNA fragment")

    redistribution = g.input(
        "redistribution", "pooling and PCR-well redistribution condition",
        MaterialKind.CONDITION)
    redistributed = g.product(
        "redistributed", "nuclei assigned a second well coordinate",
        MaterialKind.MOLECULAR, state=tagged.state)
    excluded = g.product(
        "excluded", "empty or colliding redistribution material", MaterialKind.WASTE)
    g.react(
        "redistribute", "Pool and redistribute nuclei for the second index",
        StepKind.PARTITION,
        ((tagged, InputRole.TEMPLATE), (redistribution, InputRole.CONDITION)),
        (redistributed, excluded), changes_molecule=False)
    g.discard(excluded, "does not preserve an interpretable two-coordinate cell identity")

    atac_template = g.product(
        "atac_template", "indexed accessibility-library template",
        MaterialKind.MOLECULAR, state=ATAC)
    guide_template = g.product(
        "guide_template", "integrated-guide targeted-amplification template",
        MaterialKind.MOLECULAR, state=GUIDE)
    g.react(
        "split", "Split each indexed nucleus into accessibility and guide readouts",
        StepKind.SPLIT, ((redistributed, InputRole.TEMPLATE),),
        (atac_template, guide_template), changes_molecule=False)

    reviewed = [
        "cells", "guide_vector", "integration_machinery", "integrate_guide",
        "indexed_transposome", "tn5", "tagment", "redistribution",
        "redistribute", "split",
    ]
    finish_pcr_branch(
        g, "atac", atac_template, ATAC, AP, reviewed,
        label="CRISPR-sciATAC accessibility",
        inferred_primers=True,
        limitation="The defining article describes indexed ATAC completion but the page model does not yet transcribe a universal exact PCR-primer pair.")
    finish_pcr_branch(
        g, "guide", guide_template, GUIDE, GP, reviewed,
        label="CRISPR-sciATAC guide identity",
        inferred_primers=True,
        limitation="The targeted semi-suppressive PCR is source-described, but exact landing geometry is not yet transcribed in this page model.")
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
