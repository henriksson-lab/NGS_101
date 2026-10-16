"""OAK overloading-and-unpacking combinatorial indexing."""
from batch_ngs import seg, truseq_library
from protocol_extensions import analyte_split, overload_then_unpack
from chemdraw import Construct, Row, feature
from multimodal_spatial import (finish_pcr_branch, overload_carriers,
                                recover_and_redistribute_carriers)
import illumina as il
import nextera as nx
import seqprimers as sp
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
TITLE="OAK — overloading and unpacking droplet barcodes"
NOTES="01_oak.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-024-53227-z">Wu et al., <i>Nature Communications</i> (2024)</a>.'
SUMMARY="Fixed cells or nuclei are overloaded into Chromium droplets for in situ first-round barcoding, recovered intact and redistributed into aliquots for a second index; RNA, ATAC and antibody-tag branches can share the resulting cell identity."
CAVEAT=("The paper prints example RNA aliquot and library PCR primers. Original 10x v3 "
        "and Multiome guides print the generic 16-base bead barcodes, 12-base RNA UMIs "
        "and final-library shells. The exact guide revisions and kit-specific oligos used "
        "in OAK have not yet been established; the graph marks their transfer chemistry "
        "and full primer families as unresolved.")
RNA_TSO_PCR="AAGCAGTGGTATCAACGCAGAGT"
RNA_ALIQUOT_I5="AATGATACGGCGACCACCGAGATCTACACAACGTGATACACTCTTTCCCTACACGACGCTCTTCCGATCT"
RNA_PARTIAL_P5="AATGATACGGCGACCACCGAGA"
RNA_LIBRARY_I7="CAAGCAGAAGACGGCATACGAGATCGCATGTTACGTGACTGGAGTTCAGACGTGT"
def oligos():
    return [f"<b>{name}</b> (published RNA route, 5′→3′): <code>{sequence}</code>"
            for name, sequence in (("TSO recognition PCR primer", RNA_TSO_PCR),
                                   ("example aliquot i5 primer", RNA_ALIQUOT_I5),
                                   ("partial P5 library primer", RNA_PARTIAL_P5),
                                   ("example library i7 primer", RNA_LIBRARY_I7))]
CELL_PARTS=(("droplet barcode",16),("aliquot index",8))
DROPLET=feature("droplet_barcode","cell_barcode","whitelist",whitelist="10x Chromium")
ALIQUOT=feature("aliquot_index","cell_barcode","combinatorial",group="cell_barcode",part="aliquot index")
UMI=feature("umi","umi","random")
RNA,RP=truseq_library([
    seg("droplet barcode","B"*16,"cbc",placeholder=True,feature=DROPLET),
    seg("UMI","U"*12,"umi",placeholder=True,feature=UMI),
    seg("poly(dT)VN","T"*30+"VN",placeholder=True),
    seg("3-prime cDNA","X"*34,placeholder=True)],
    "OAK RNA",i5_nt=8,i7_nt=10,i5_feature=ALIQUOT,i5_name="aliquot index (i5)")
MULTI_RNA=Construct(RNA.segments,name="OAK Multiome RNA")
ATAC=Construct([
    seg("P5",il.P5,"p5"),
    seg("droplet barcode","B"*16,"cbc",placeholder=True,feature=DROPLET),
    seg("10x spacer","CGCGTCTG"),seg("S5",nx.S5,"s5"),
    seg("left mosaic end",nx.ME,"me"),
    seg("accessible genomic DNA","X"*34,placeholder=True),
    seg("right mosaic end reverse complement",nx.ME_RC,"me"),
    seg("S7 reverse complement",nx.S7_RC,"s7"),
    seg("aliquot index (i7)","I"*8,"cbc",placeholder=True,feature=ALIQUOT),
    seg("P7 reverse complement",il.P7_RC,"p7")],name="OAK Multiome ATAC")
AP=tuple(sp.NEXTERA[k] for k in ("R1","I1","I2","R2"))
_atac_sites=sp.verify(ATAC,AP,required_roles=tuple(p.role for p in AP))
if _atac_sites:
    raise ValueError("OAK Multiome ATAC: "+"; ".join(_atac_sites))
FINAL_LIBRARIES=(("RNA-only library",RNA,RP,"The aliquot index is read in i5; droplet barcode and 12-base UMI follow Read 1.","Published OAK RNA run: 28/10/8/90 cycles (R1/i7/i5/R2)."),
                 ("Multiome RNA library",MULTI_RNA,RP,"The same RNA read architecture follows prior bulk transposition of intact nuclei.","RNA and ATAC identities are carried by the same nucleus."),
                 ("Multiome ATAC library",ATAC,AP,"The 16-base droplet barcode and 8-base aliquot index occupy separate ends.","Published OAK ATAC run: 50/8/24/49 cycles (R1/i7/i5/R2)."))
def sections():
    return [("Fix cells or nuclei so they survive droplet recovery",[Row(chunks=[("fixed intact cell/nucleus carrying RNA and chromatin", None, False)])],"The carrier remains intact across indexing rounds."),
            ("Transpose chromatin before GEM formation for Multiome",[Row(chunks=[("fixed nucleus → bulk Tn5 transposition → tagged chromatin retained in nucleus", "me", False)])],"OAK Multiome Methods transpose nuclei at 37 °C before loading the transposed nuclei into Chromium."),
            ("Overload Chromium droplets and perform first-round barcoding",overload_then_unpack(droplet_barcode=16,well_barcode=8)[:1],"RNA-only cells enter directly; Multiome nuclei enter after transposition."),
            ("Break droplets, mix carriers and redistribute into aliquots",overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:],"Aliquot-specific PCR supplies the orthogonal secondary index."),
            ("Separate configured modalities",analyte_split("RNA / antibody-tag library","ATAC library when Multiome is used"),"The demonstration supports RNA alone and joint RNA–ATAC; the cell identity is shared."),]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    sample=g.input("sample","cells or nuclei for OAK",MaterialKind.BIOLOGICAL,primary=True)
    fixative=g.input("fixative","carrier-preserving fixation condition",MaterialKind.CONDITION)
    fixed=g.product("fixed","fixed cells or nuclei",MaterialKind.MOLECULAR,state=sections()[0][1])
    g.react("fix","Fix carriers for the configured OAK branch",StepKind.MAGIC,
            ((sample,InputRole.TEMPLATE),(fixative,InputRole.CONDITION)),(fixed,),
            reason="RNA-only methanol fixation and Multiome formaldehyde fixation have different chemistry; this shared carrier state does not assert a single molecular transform.")
    fixed_rna=g.product("fixed_rna","fixed RNA-only cells",MaterialKind.MOLECULAR,state=fixed.state)
    fixed_multiome=g.product("fixed_multiome","fixed Multiome nuclei",MaterialKind.MOLECULAR,state=fixed.state)
    g.react("choose_assay","Route fixed carriers to RNA-only or paired Multiome workflow",StepKind.SPLIT,
            ((fixed,InputRole.TEMPLATE),),(fixed_rna,fixed_multiome),changes_molecule=False)

    # RNA-only cells enter GEMs directly. Its source-known 3-prime v3 bead architecture
    # is distinct from the still-unreconstructed kit-specific transfer geometry.
    rna_droplets,rna_bead,rna_overload=overload_carriers(
        g,fixed_rna,prefix="rna_droplet",
        state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[:1],
        bead_label="10x 3-prime v3 oligo(dT) barcode bead",
        source=SOURCE+"; 10x CG000183 Rev A Appendix")
    rt=g.input("reverse_transcriptase","in-droplet reverse transcriptase",MaterialKind.CATALYST)
    rna_copied=g.product("rna_copied","fixed RNA-only cells retaining barcoded cDNA",MaterialKind.MOLECULAR,
                         state=rna_droplets.state)
    g.react("rna_reverse_transcribe","Copy RNA and add the GEM barcode",StepKind.REVERSE_TRANSCRIPTION,
            ((rna_droplets,InputRole.TEMPLATE),(rna_bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(rna_copied,))
    rna_aliquoted,rna_recovery=recover_and_redistribute_carriers(
        g,rna_copied,prefix="rna_carrier",
        recovered_state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:],
        aliquoted_state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:])
    rna=g.product("rna_template","RNA-only cDNA released from indexed aliquots",MaterialKind.MOLECULAR,
                  state=RNA,inferred=True,source=SOURCE+" (finished duplex is a library view, not the released first strand)")
    rna_heat=g.input("rna_heat","80 °C aliquot heating condition",MaterialKind.CONDITION)
    g.react("rna_release","Heat aliquots to release first-strand cDNA",StepKind.THERMAL_TREATMENT,
            ((rna_aliquoted,InputRole.TEMPLATE),(rna_heat,InputRole.CONDITION)),(rna,))

    # OAK Multiome explicitly transposes fixed nuclei in bulk before GEM loading.
    transposome=g.input("transposome","TDE1/ATAC Buffer B transposition reagents",MaterialKind.ADAPTER,
                        inferred=True,source=SOURCE+" (loaded Tn5 adapter composition depends on the kit revision)")
    tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=g.product("tagged","fixed nuclei retaining tagged chromatin and RNA",MaterialKind.MOLECULAR,
                     state=sections()[1][1],inferred=True,source=SOURCE+" (adapter loading not assembled base for base)")
    tagment_waste=g.product("tagment_waste","non-library transposition products",MaterialKind.WASTE)
    g.react("tagment","Transpose fixed nuclei before loading them into GEMs",StepKind.TAGMENTATION,
            ((fixed_multiome,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),
            (tagged,tagment_waste),source=SOURCE+" Methods, OAK paired snRNA-Seq and snATAC-Seq")
    g.discard(tagment_waste,"does not retain a library-eligible tagged fragment")
    multi_droplets,multi_bead,multi_overload=overload_carriers(
        g,tagged,prefix="multiome_droplet",
        state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[:1],
        bead_label="10x Multiome RNA and ATAC barcode beads",
        source=SOURCE+"; 10x CG000338 Rev F Appendix")
    multi_barcoded=g.product("multi_barcoded","GEM-barcoded RNA and tagged chromatin in intact nuclei",
                             MaterialKind.MOLECULAR,state=multi_droplets.state,inferred=True,
                             source=SOURCE+" (kit-specific barcode transfer is not reconstructed)")
    g.react("multiome_barcode","Barcode RNA and tagged chromatin inside Multiome GEMs",StepKind.MAGIC,
            ((multi_droplets,InputRole.TEMPLATE),(multi_bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),
            (multi_barcoded,),reason="The defining paper delegates the Multiome GEM barcoding chemistry to the 10x kit; its exact OAK-era reagent set is under source review.")
    multi_aliquoted,multi_recovery=recover_and_redistribute_carriers(
        g,multi_barcoded,prefix="multiome_carrier",
        recovered_state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:],
        aliquoted_state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:])
    multi_released=g.product("multi_released","cDNA and ATAC fragments released from indexed nuclei",
                             MaterialKind.MOLECULAR,state=multi_aliquoted.state)
    multi_heat=g.input("multiome_heat","80 °C aliquot heating condition",MaterialKind.CONDITION)
    g.react("multiome_release","Heat aliquots to release first strands and ATAC fragments",
            StepKind.THERMAL_TREATMENT,
            ((multi_aliquoted,InputRole.TEMPLATE),(multi_heat,InputRole.CONDITION)),(multi_released,))
    multi_rna=g.product("multiome_rna_template","paired Multiome RNA template",MaterialKind.MOLECULAR,
                        state=MULTI_RNA,inferred=True,source=SOURCE+" (finished duplex is a library view)")
    atac=g.product("atac_template","paired Multiome ATAC template",MaterialKind.MOLECULAR,
                   state=ATAC,inferred=True,source=SOURCE+" (finished duplex is a library view)")
    g.react("separate_modalities","Separate released Multiome RNA and ATAC products",StepKind.SPLIT,
            ((multi_released,InputRole.TEMPLATE),),(multi_rna,atac),changes_molecule=False)

    reviewed=["sample","fixative","fix","choose_assay",*rna_overload,"reverse_transcriptase",
              "rna_reverse_transcribe",*rna_recovery,"rna_heat","rna_release","transposome","tn5","tagment",
              *multi_overload,"multiome_barcode",*multi_recovery,"multiome_heat","multiome_release","separate_modalities"]
    finish_pcr_branch(g,"rna",rna,RNA,RP,reviewed,label="OAK RNA-only",limitation="The paper prints RNA primer examples but the graph does not yet assemble every intermediate PCR product.")
    finish_pcr_branch(g,"multiome_rna",multi_rna,MULTI_RNA,RP,reviewed,label="OAK Multiome RNA",limitation="The paper delegates the complete Multiome kit primer set and its historical revision to 10x.")
    finish_pcr_branch(g,"atac",atac,ATAC,AP,reviewed,label="OAK Multiome ATAC",limitation="The generic vendor ATAC shell is known, while OAK's exact kit revision and indexed-PCR primer family need separate reconstruction.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
