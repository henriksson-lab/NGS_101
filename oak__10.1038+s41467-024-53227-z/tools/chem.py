"""OAK overloading-and-unpacking combinatorial indexing."""
from batch_ngs import seg, truseq_library
from protocol_extensions import analyte_split, overload_then_unpack
from chemdraw import Construct, Row, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from multimodal_spatial import (overload_carriers,
                                recover_and_redistribute_carriers)
import illumina as il
import nextera as nx
import seqprimers as sp
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
# Methods: aliquots are "immediately stored at −80 °C" and "a desired number" are
# processed into sub-libraries (RNA, PDF lines 594–604; Multiome, lines 680–685).
def _stored_aliquots(where):
    return ("other indexed aliquots stored at −80 °C (not processed in this representation)",
            Disposition.STORED,
            f"Methods ({where}): the aliquots are immediately stored at −80 °C, and only a "
            "desired number are thawed and converted into sub-libraries.")
TITLE="OAK — overloading and unpacking droplet barcodes"
NOTES="01_oak.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-024-53227-z">Wu et al., <i>Nature Communications</i> (2024)</a>.'
SUMMARY="Fixed cells or nuclei are overloaded into Chromium droplets for in situ first-round barcoding, recovered intact and redistributed into aliquots for a second index; RNA, ATAC and antibody-tag branches can share the resulting cell identity."
CAVEAT=("The paper prints example RNA aliquot (i5) and library i7 PCR primers; the underlined "
        "example i7 index is the 10-base CGCATGTTAC, matching the 10 acquired i7 cycles. The "
        "aliquot-primer landing on the bead's partial Read 1 uses the original 10x CG000183 Rev A "
        "v3 bead; the Multiome RNA branch is drawn with the same cDNA view because OAK states it "
        "is processed as in OAK scRNA-Seq, but the OAK-era Multiome bead revision is unverified. "
        "Standard Chromium library construction, the 10x pre-amp primer mix (PN 20002714) and "
        "ATAC Sample Index N primers are delegated to 10x and remain unresolved boundaries.")
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
RNA_READS=(("Read 1",28),("Index 1 (i7)",10),("Index 2 (i5)",8),("Read 2",90))
ATAC_READS=(("Read 1",50),("Index 1 (i7)",8),("Index 2 (i5)",24),("Read 2",49))
FINAL_LIBRARIES=(("RNA-only library",RNA,RP,"The aliquot index is read in i5; droplet barcode and 12-base UMI follow Read 1.","Published OAK RNA run: 28/10/8/90 cycles (R1/i7/i5/R2)."),
                 ("Multiome RNA library",MULTI_RNA,RP,"The same RNA read architecture follows prior bulk transposition of intact nuclei.","RNA and ATAC identities are carried by the same nucleus."),
                 ("Multiome ATAC library",ATAC,AP,"The 16-base droplet barcode and 8-base aliquot index occupy separate ends.","Published OAK ATAC run: 50/8/24/49 cycles (R1/i7/i5/R2)."))
def sections():
    return [("Fix cells or nuclei so they survive droplet recovery",[Row(chunks=[("fixed intact cell/nucleus carrying RNA and chromatin", None, False)])],"The carrier remains intact across indexing rounds."),
            ("Transpose chromatin before GEM formation for Multiome",[Row(chunks=[("fixed nucleus → bulk Tn5 transposition → tagged chromatin retained in nucleus", "me", False)])],"OAK Multiome Methods transpose nuclei at 37 °C before loading the transposed nuclei into Chromium."),
            ("Overload Chromium droplets and perform first-round barcoding",overload_then_unpack(droplet_barcode=16,well_barcode=8)[:1],"RNA-only cells enter directly; Multiome nuclei enter after transposition."),
            ("Break droplets, mix carriers and redistribute into aliquots",overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:],"Aliquot-specific PCR supplies the orthogonal secondary index."),
            ("Separate configured modalities",analyte_split("RNA / antibody-tag library","ATAC library when Multiome is used"),"The demonstration supports RNA alone and joint RNA–ATAC; the cell identity is shared."),]


# The OAK aliquot-index primer's 3' end primes on the bead-derived partial Read 1
# (vendor CG000183 Rev A bead; source transcription checked in selftest) and the TSO
# recognition primer on the template-switch end. Intermediate cDNA view only.
BEAD_PARTIAL_R1="CTACACGACGCTCTTCCGATCT"
CDNA=Construct([
    seg("partial Read 1 (bead)",BEAD_PARTIAL_R1,"r1"),
    seg("droplet barcode","B"*16,"cbc",placeholder=True,feature=DROPLET),
    seg("UMI","U"*12,"umi",placeholder=True,feature=UMI),
    seg("poly(dT)VN","T"*30+"VN",placeholder=True),
    seg("cDNA","X"*34,placeholder=True),
    seg("TSO complement",revcomp(RNA_TSO_PCR))],name="OAK first-strand cDNA (duplex view)")
CDNA_PCR=PCRGeometry.exact(CDNA,PCRPrimer("aliquot i5 primer",RNA_ALIQUOT_I5),
                           PCRPrimer("TSO recognition primer",RNA_TSO_PCR),
                           statement="OAK scRNA-Seq Methods: TSO recognition primer plus aliquot-index primer amplify first-strand cDNA; the primer's P5/i5 tail is added as a 5' overhang.")


def _rna_library(g,prefix,first_strand,final,label,reviewed):
    """cDNA PCR with the published aliquot primer, vendor library steps, then i7 PCR."""
    tso=g.input(prefix+"_tso_primer","TSO recognition primer",MaterialKind.OLIGO,state=RNA_TSO_PCR,
                source=SOURCE+" Methods, OAK scRNA-Seq")
    i5=g.input(prefix+"_aliquot_primer","aliquot-specific secondary-index (i5) primer; example printed",
               MaterialKind.OLIGO,state=RNA_ALIQUOT_I5,source=SOURCE+" Methods, OAK scRNA-Seq (example index)")
    amp_mix=g.input(prefix+"_cdna_polymerase","Chromium 3′ cDNA amplification PCR reagents",MaterialKind.CATALYST)
    cdna=g.product(prefix+"_cdna","aliquot-indexed amplified cDNA",MaterialKind.MOLECULAR)
    g.react(prefix+"_cdna_pcr","Amplify cDNA with TSO and aliquot-index primers",StepKind.PCR,
            ((first_strand,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),(tso,InputRole.STRUCTURAL),
             (amp_mix,InputRole.CATALYST)),(cdna,),geometry=CDNA_PCR)
    kit=g.input(prefix+"_library_kit","Chromium 3′ v3 library kit (PN-1000095) fragmentation, end-repair/A-tailing and adapter ligation",
                MaterialKind.CONDITION)
    ligated=g.product(prefix+"_ligated","adapter-ligated sub-library before index PCR",MaterialKind.MOLECULAR)
    g.react(prefix+"_construct","Construct sub-library by the standard Chromium protocol",StepKind.MAGIC,
            ((cdna,InputRole.TEMPLATE),(kit,InputRole.CONDITION)),(ligated,),
            reason="OAK delegates fragmentation, end repair/A-tailing and adapter ligation to the standard 10x Chromium 3′ protocol; that vendor revision is not reconstructed here.")
    p5=g.input(prefix+"_p5_primer","partial P5 primer",MaterialKind.OLIGO,state=RNA_PARTIAL_P5,
               source=SOURCE+" Methods, OAK scRNA-Seq")
    i7=g.input(prefix+"_i7_primer","i7 index primer; example printed",MaterialKind.OLIGO,state=RNA_LIBRARY_I7,
               source=SOURCE+" Methods, OAK scRNA-Seq (example index)")
    pol=g.input(prefix+"_library_polymerase","Chromium library PCR reagents",MaterialKind.CATALYST)
    library=g.product(prefix+"_library",label+" sequencing library",MaterialKind.LIBRARY,state=final)
    g.react(prefix+"_pcr","Library PCR with partial P5 and i7 index primers",StepKind.PCR,
            ((ligated,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),(i7,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),
            (library,),geometry=PCRGeometry.unavailable(final,reason="The adapter-ligated intermediate follows the vendor library protocol and is not assembled base for base."))
    g.sequence(library,SequencingPlan(library=final,primers=tuple(RP),acquired_roles=tuple(p.role for p in RP),
               read_lengths=dict(RNA_READS)))
    reviewed.extend((tso.id,i5.id,amp_mix.id,prefix+"_cdna_pcr",kit.id,prefix+"_construct",
                     p5.id,i7.id,pol.id,prefix+"_pcr",library.id))
    return library


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
    g.react("rna_reverse_transcribe","Reverse transcribe in droplets (53 °C, 45 min), adding the GEM barcode",StepKind.REVERSE_TRANSCRIPTION,
            ((rna_droplets,InputRole.TEMPLATE),(rna_bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(rna_copied,))
    rna_aliquoted,rna_recovery=recover_and_redistribute_carriers(
        g,rna_copied,prefix="rna_carrier",
        recovered_state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:],
        aliquoted_state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:],
        aliquot_remainder=_stored_aliquots("OAK scRNA-Seq"))
    rna=g.product("rna_template","RNA-only first-strand cDNA released from indexed aliquots",MaterialKind.MOLECULAR)
    rna_heat=g.input("rna_heat","80 °C aliquot heating condition",MaterialKind.CONDITION)
    g.react("rna_release","Heat aliquots to 80 °C for 5 min to release first-strand cDNA",StepKind.THERMAL_TREATMENT,
            ((rna_aliquoted,InputRole.TEMPLATE),(rna_heat,InputRole.CONDITION)),(rna,))
    rna_silane=g.input("rna_silane","Dynabeads Silane Viral NA purification",MaterialKind.SELECTION)
    rna_purified=g.product("rna_purified","purified first-strand cDNA (35 µL elution)",MaterialKind.MOLECULAR)
    rna_silane_waste=g.product("rna_silane_waste","silane-bead unbound fraction and washes",MaterialKind.WASTE)
    g.react("rna_silane_cleanup","Purify first-strand cDNA on silane beads",StepKind.CLEANUP,
            ((rna,InputRole.TEMPLATE),(rna_silane,InputRole.SELECTION)),(rna_purified,rna_silane_waste))
    g.discard(rna_silane_waste,"removed by the silane cleanup")

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
        aliquoted_state=overload_then_unpack(droplet_barcode=16,well_barcode=8)[1:],
        aliquot_remainder=_stored_aliquots("OAK paired snRNA-Seq and snATAC-Seq"))
    multi_released=g.product("multi_released","cDNA and ATAC fragments released from indexed nuclei",
                             MaterialKind.MOLECULAR,state=multi_aliquoted.state)
    multi_heat=g.input("multiome_heat","80 °C aliquot heating condition",MaterialKind.CONDITION)
    g.react("multiome_release","Heat aliquots to 80 °C for 5 min to release first strands and ATAC fragments",
            StepKind.THERMAL_TREATMENT,
            ((multi_aliquoted,InputRole.TEMPLATE),(multi_heat,InputRole.CONDITION)),(multi_released,))
    silane=g.input("multiome_silane","Dynabeads Silane Viral NA purification",MaterialKind.SELECTION)
    multi_purified=g.product("multiome_purified","purified first-strand cDNA and ATAC fragments",MaterialKind.MOLECULAR)
    silane_waste=g.product("multiome_silane_waste","silane-bead unbound fraction and washes",MaterialKind.WASTE)
    g.react("multiome_silane_cleanup","Purify released products on silane beads",StepKind.CLEANUP,
            ((multi_released,InputRole.TEMPLATE),(silane,InputRole.SELECTION)),(multi_purified,silane_waste))
    g.discard(silane_waste,"removed by the silane cleanup")
    preamp_mix=g.input("preamp_reagents","4 µL 10x pre-amp primers (PN 20002714) with NEBNext High-Fidelity 2X PCR Master Mix",
                       MaterialKind.CONDITION)
    preamped=g.product("preamplified","10-cycle pre-amplified cDNA and ATAC products",MaterialKind.MOLECULAR)
    g.react("multiome_preamp","Pre-amplify both modalities for 10 cycles",StepKind.MAGIC,
            ((multi_purified,InputRole.TEMPLATE),(preamp_mix,InputRole.CONDITION)),(preamped,),
            reason="OAK names the 10x pre-amp primer mix by part number only; its primer sequences and landing sites are not printed in the defining paper.")
    spri=g.input("preamp_spri","1.6X SPRI cleanup, 40 µL EB elution",MaterialKind.SELECTION)
    preamp_clean=g.product("preamp_clean","SPRI-recovered pre-amplified products",MaterialKind.MOLECULAR)
    spri_waste=g.product("preamp_spri_waste","SPRI supernatant and washes",MaterialKind.WASTE)
    g.react("preamp_cleanup","Clean pre-amplification product with 1.6X SPRI",StepKind.CLEANUP,
            ((preamped,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(preamp_clean,spri_waste))
    g.discard(spri_waste,"removed by SPRI")
    multi_rna=g.product("multiome_rna_template","pre-amplified product for the Multiome RNA library",MaterialKind.MOLECULAR)
    atac=g.product("atac_template","10 µL pre-amplified product for the ATAC library",MaterialKind.MOLECULAR)
    g.react("separate_modalities","Split pre-amplified product into ATAC (10 µL) and RNA libraries",StepKind.SPLIT,
            ((preamp_clean,InputRole.TEMPLATE),),(multi_rna,atac),changes_molecule=False)
    reviewed=["sample","fixative","fix","choose_assay",*rna_overload,"reverse_transcriptase",
              "rna_reverse_transcribe",*rna_recovery,"rna_heat","rna_release","rna_silane","rna_silane_cleanup",
              "transposome","tn5","tagment",
              *multi_overload,"multiome_barcode",*multi_recovery,"multiome_heat","multiome_release",
              "multiome_silane","multiome_silane_cleanup","preamp_reagents","multiome_preamp",
              "preamp_spri","preamp_cleanup","separate_modalities"]
    _rna_library(g,"rna",rna_purified,RNA,"OAK RNA-only",reviewed)
    _rna_library(g,"multiome_rna",multi_rna,MULTI_RNA,"OAK Multiome RNA",reviewed)
    atac_p5=g.input("atac_p5_primer","partial P5 primer (0.6 µL 100 µM)",MaterialKind.OLIGO,state=RNA_PARTIAL_P5,
                    source=SOURCE+" Methods, OAK paired snRNA-Seq and snATAC-Seq")
    atac_i7=g.input("atac_sample_index","10x Sample Index N primer (PN 1000212); aliquot-specific i7",MaterialKind.OLIGO,
                    source=SOURCE+" (sequence delegated to the 10x kit; not printed in OAK)")
    atac_pol=g.input("atac_polymerase","NEBNext High-Fidelity 2X PCR Master Mix",MaterialKind.CATALYST)
    atac_amp=g.product("atac_amplified","aliquot-indexed ATAC library before size selection",MaterialKind.MOLECULAR,state=ATAC,
                       inferred=True,source=SOURCE+" (generic vendor ATAC shell; kit-specific primer family not reconstructed)")
    g.react("atac_pcr","Index ATAC library: 98 °C 30 s, n×(98 °C 10 s, 67 °C 30 s, 72 °C 20 s), 72 °C 2 min",StepKind.PCR,
            ((atac,InputRole.TEMPLATE),(atac_p5,InputRole.STRUCTURAL),(atac_i7,InputRole.STRUCTURAL),(atac_pol,InputRole.CATALYST)),
            (atac_amp,),geometry=PCRGeometry.unavailable(ATAC,reason="Sample Index N sequences and the pre-amplified ATAC intermediate are delegated to the 10x kit."))
    atac_size=g.input("atac_size_selection","double-sided SPRI size selection (standard Chromium protocol)",MaterialKind.SELECTION)
    atac_lib=g.product("atac_library","OAK Multiome ATAC sequencing library",MaterialKind.LIBRARY,state=ATAC,inferred=True,
                       source=SOURCE+" (generic vendor ATAC shell)")
    atac_off=g.product("atac_size_excluded","fragments outside the double-sided size window",MaterialKind.WASTE)
    g.react("atac_size_select","Double-sided size selection",StepKind.SELECTION,
            ((atac_amp,InputRole.TEMPLATE),(atac_size,InputRole.SELECTION)),(atac_lib,atac_off))
    g.discard(atac_off,"removed by double-sided SPRI")
    g.sequence(atac_lib,SequencingPlan(library=ATAC,primers=tuple(AP),acquired_roles=tuple(p.role for p in AP),
               read_lengths=dict(ATAC_READS)))
    reviewed.extend(("atac_p5_primer","atac_sample_index","atac_polymerase","atac_pcr",
                     "atac_size_selection","atac_size_select","atac_library"))
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
