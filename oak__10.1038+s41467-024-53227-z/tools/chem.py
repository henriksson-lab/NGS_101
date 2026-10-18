"""OAK source-owned fixed-carrier recovery and mixed-input aliquot graph."""
from chemdraw import Construct,Scene,Segment,feature,revcomp
from protocol_model import molecule,role
from pcr import PCRGeometry
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind,Terminal,Disposition,Highlight
import illumina as il
import seqprimers as sp
TITLE="OAK — overloading and unpacking droplet barcodes"
NOTES="01_oak.html"
SOURCE='Wu et al., NatureCommunications(2024), <a href="https://doi.org/10.1038/s41467-024-53227-z">doi:10.1038/s41467-024-53227-z</a>; published Methods OAKscRNA-Seq and paired snRNA/snATAC.'
SUMMARY="Fixed RNA-only cells or bulk-transposed Multiome nuclei survive overloaded GEM recovery and redistribution. Aliquot-index PCR distinguishes cells sharing a droplet identity."
CAVEAT="RNA-only methanol and paired formaldehyde fixation are distinct source recipes. The source prints representative secondary-index/P5/i7 primers and positive run cycles; it does not establish full historical commercial end assembly in this model. Prior CG000183RevA3′v3 bead architecture is a separately checked reference, not proof of the OAK-era Multiome bead revision. Native first-cycle and complete ATAC/fragmentation/ligation geometry remain bounded. Mixed preamplified cDNA/ATAC enters both physical aliquots; no modality purification is inferred."
RNA_TSO_PCR="AAGCAGTGGTATCAACGCAGAGT"
RNA_ALIQUOT_I5="AATGATACGGCGACCACCGAGATCTACACAACGTGATACACTCTTTCCCTACACGACGCTCTTCCGATCT"
RNA_PARTIAL_P5="AATGATACGGCGACCACCGAGA"
RNA_LIBRARY_I7="CAAGCAGAAGACGGCATACGAGATCGCATGTTACGTGACTGGAGTTCAGACGTGT"

DROPLET=feature("droplet_barcode","cell_barcode","whitelist",whitelist="source Chromium bead; selected historical member/width unassigned")
ALIQUOT=feature("aliquot_index","cell_barcode","combinatorial",group="cell_barcode",part="aliquot index")
UMI=feature("umi","umi","random")
LIBRARY_I7=feature("oak_library_i7","sample_index","whitelist",whitelist="published examplei7")
I5_VALUE=RNA_ALIQUOT_I5[len(il.P5):len(il.P5)+8]
I7_VALUE=RNA_LIBRARY_I7[len(il.P7):len(il.P7)+10]
R1_ARM=RNA_ALIQUOT_I5[len(il.P5)+8:]
if R1_ARM!=il.TRUSEQ_READ1:raise ValueError("Published aliquot primer differs from declared R1 arm")

def identities(*,v3=False):
    if v3:return [Segment("droplet barcode","N"*16,placeholder=True,feature=DROPLET),Segment("UMI","N"*12,placeholder=True,feature=UMI)]
    return [role("droplet cell barcode; historical native width unassigned",molecular_feature=DROPLET),role("molecular UMI; historical native width unassigned",molecular_feature=UMI)]

def endpoint(label,*,v3=False):
    return Construct([Segment("P5",il.P5),Segment("aliquot i5",I5_VALUE,feature=ALIQUOT),Segment("source primer R1 arm",R1_ARM),*identities(v3=v3),role("cDNA/library interval; full historical adapter assembly unassigned",inferred=True),Segment("source i7 primer binding complement",revcomp(RNA_LIBRARY_I7[len(il.P7)+10:])),Segment("sample i7 complement",revcomp(I7_VALUE),feature=LIBRARY_I7),Segment("P7 complement",il.P7_RC)],name=label)
RNA=endpoint("OAK RNA-only source-primer-derived reference",v3=True)
MULTI_RNA=endpoint("OAK MultiomeRNA source-primer-derived reference")
ATAC=Construct([role("droplet barcode; native ATAC strand layout unassigned",molecular_feature=DROPLET,inferred=True),role("ATAC fragment and native commercial adapter assembly",inferred=True),role("aliquot sample-indexN; physical width/member unassigned",molecular_feature=ALIQUOT,inferred=True)],name="OAK MultiomeATAC bounded endpoint")
RP=(sp.TRUSEQ["R1"],);AP=()
FINAL_LIBRARIES=(("RNA-only library",RNA,RP,"Source-primer-derived ends; complete native geometry unassigned.","Published28/10/8/90 run."),("Multiome RNA library",MULTI_RNA,RP,"Same source RNA primer strategy applied to mixed Multiome aliquot; native bead revision unassigned.","RNA sequencing performed as OAKscRNA-Seq."),("Multiome ATAC library",ATAC,AP,"Native commercial ATAC geometry remains unassigned; identifiers explicit.","Published50/8/24/49 run."))
SHOW_READ_CONTENT=False
SEQUENCING_UNAVAILABLE="R1 sites derive from the published RNA aliquot primer. Full historical RNA mate/index and ATAC binding geometry are unassigned; no generic vendor endpoint supplies them."
RNA_READS=(("Read 1",28),("Index 1 (i7)",10),("Index 2 (i5)",8),("Read 2",90))
ATAC_READS=(("Read 1",50),("Index 1 (i7)",8),("Index 2 (i5)",24),("Read 2",49))
READ_SPECS={name:tuple(sp.read_spec(read,cycles,content,use) for (read,cycles),(content,use) in zip(RNA_READS if name!="Multiome ATAC library" else ATAC_READS,(("RNA droplet barcode and UMI; native widths unassigned","extract"),("library i7 sample index","demultiplex"),("aliquot i5 index","extract_and_demultiplex"),("cDNA insert","include")) if name!="Multiome ATAC library" else (("genomic insert","include"),("aliquot-specific sample indexN","extract_and_demultiplex"),("droplet identity read; native barcode/spacer geometry unassigned","extract"),("genomic insert","include")))) for name,*_ in FINAL_LIBRARIES}
READ_SPECS["RNA-only library"]=(sp.read_spec("Read 1",28,"source3′v3 N16 droplet barcode and N12 UMI; CG000183RevA component scope","extract"),*READ_SPECS["RNA-only library"][1:])
READ_SPECS_NOTE="Published acquired cycle counts are positive source facts; physical barcode widths and unplaced primer geometry are separate model boundaries."

def oligos():return tuple(name+" —5′-"+value+"-3′" for name,value in (("TSO recognition primer",RNA_TSO_PCR),("representative aliquot i5 primer",RNA_ALIQUOT_I5),("partialP5 primer",RNA_PARTIAL_P5),("representative library i7 primer",RNA_LIBRARY_I7)))

def carrier(label,*,copied=False,atac=False,duplex=False):
    sc=Scene()
    sc.strand("chromatin",[role("independent "+("tagged/ATAC-derived" if atac else "genomic")+" DNA population; native ends unassigned")],label=label+" DNA population")
    parts=[role("cellular RNA")]
    if copied:
        if not atac:
            parts=[Segment("sourcev3 partialR1","CTACACGACGCTCTTCCGATCT"),*identities(v3=True),Segment("poly(dT)30","T"*30),role("first-strand antisense cDNA"),role("copiedTSO end; native full extent unassigned")]
        else:
            parts=[role("source historical Multiome bead-derived RT end; full internal barcode/UMI order unassigned",inferred=True),*identities(),role("RNA-derived cDNA and TSOcopy; full native assembly unassigned",inferred=True)]
    sc.strand("transcripts",parts,label=label+(" amplified cDNA strand-class; complete duplex geometry unassigned" if duplex else " single-strand cDNA population" if copied else " independent RNA population"))
    if duplex:sc.strand("cDNA opposite strands",[role("complementary cDNA strand classes; complete native extent unassigned")],label="opposite amplified cDNA strands")
    return sc

def clean(g,key,template,recipe,library=False):
    selection=g.input(key+"_selection",recipe,MaterialKind.SELECTION)
    kept=g.product(key+"_clean","recovered "+template.label,MaterialKind.LIBRARY if library else MaterialKind.MOLECULAR,state=template.state)
    waste=g.product(key+"_waste","source removed liquid/wash fraction",MaterialKind.WASTE)
    g.react(key,"Recover source retained fraction",StepKind.CLEANUP,((template,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),(kept,waste),changes_molecule=False);g.discard(waste,"Source physical emulsion/bead/centrifuge/wash transfer; not computational QC")
    return kept

def recovery(g,key,template):
    recovered=clean(g,key+"_recovery",template,"125µl recoveryagent; aqueous transfer;800µl3XSSC,650g5min;1ml3XSSC wash;215µl resuspension")
    rule=g.input(key+"_aliquot_rule","10µl counting sample; remaining evenly divided into ~4000carrier aliquots",MaterialKind.CONDITION)
    aliquot=g.product(key+"_aliquot","chosen indexed carrier aliquot",MaterialKind.MOLECULAR,state=recovered.state)
    stored=g.product(key+"_stored","other aliquots immediately stored−80°C",MaterialKind.MOLECULAR,state=recovered.state)
    counting=g.product(key+"_counting","10µl carrier counting sample",MaterialKind.MOLECULAR,state=recovered.state)
    g.react(key+"_redistribute","Count carriers and divide aliquots",StepKind.SPLIT,((recovered,InputRole.TEMPLATE),(rule,InputRole.CONDITION)),(aliquot,stored,counting),changes_molecule=False)
    g.mark_terminal(Terminal(stored.id,Disposition.STORED,"Published Methods immediately store aliquots−80°C until chosen for library preparation."));g.mark_terminal(Terminal(counting.id,Disposition.UNRESOLVED,"Source uses10µl for counting; subsequent sample handling unassigned."))
    heat=g.input(key+"_heat","80°C5min release",MaterialKind.CONDITION)
    released=g.product(key+"_released","released mixed cDNA/DNA populations",MaterialKind.MOLECULAR,state=aliquot.state)
    g.react(key+"_release","Heat chosen aliquot80°C5min",StepKind.THERMAL_TREATMENT,((aliquot,InputRole.TEMPLATE),(heat,InputRole.CONDITION)),(released,))
    return clean(g,key+"_silane",released,"DynabeadsSilaneViralNA purification;35µl elution forRNA")

def aliquot_indexed_cdna(*,multiome=False):
    sc=Scene()
    parts=[Segment("installed P5",il.P5),Segment("installed aliquot i5",I5_VALUE,feature=ALIQUOT),Segment("installed source R1 arm",R1_ARM),*identities(v3=not multiome)]
    parts += [role("native cDNA interval and copied source-end extent unassigned",inferred=True),Segment("source TSO-primer complement",revcomp(RNA_TSO_PCR))]
    sc.strand("cDNA",parts,label="source-primer-derived indexed cDNA strand; full native assembly bounded")
    sc.strand("opposite cDNA",[role("complementary amplified cDNA strand; complete native extent unassigned")],label="opposite amplified strand population")
    sc.strand("carried DNA",[role("independent ATAC/chromatin background; amplification/recovery unassigned" if multiome else "independent genomic/background DNA; recovery unassigned")],label="carried independent DNA population")
    return sc

def rna_finish(g,key,template,final):
    i5=g.input(key+"_i5","source example aliquot i5 primer",MaterialKind.OLIGO,state=Construct([Segment("P5",il.P5),Segment("aliquot i5",I5_VALUE,feature=ALIQUOT),Segment("R1",R1_ARM)]))
    tso=g.input(key+"_tso","source TSO recognition primer",MaterialKind.OLIGO,state=RNA_TSO_PCR)
    pol=g.input(key+"_cdna_pol","source standard Chromium cDNA PCR reagents; exact historical formulation/cycles unassigned",MaterialKind.CATALYST)
    amplified=g.product(key+"_cdna","aliquot-indexed cDNA target plus carried mixed background",MaterialKind.MOLECULAR,state=aliquot_indexed_cdna(multiome=key=="multi_rna"))
    g.react(key+"_cdna_pcr","Amplify cDNA with source i5/TSO primers",StepKind.PCR,((template,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),(tso,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(template.state,reason="Actual aliquot input contains ss cDNA (RNA-only) or mixed preamplified cDNA/ATAC (Multiome). Expected endpoint is not the PCR template; full historical native ends unassigned."))
    kit=g.input(key+"_kit","Chromium3′v3 librarykit PN1000095; optional2–4aliquot pooling; source fragmentation/endrepair/ligation recipe delegated",MaterialKind.CONDITION)
    ligated=g.product(key+"_ligated","actual vendor-prepared pre-index cDNA fragments",MaterialKind.MOLECULAR,state=molecule("pre-index cDNA fragments",(role("fragmented/adapter-ligated cDNA; complete native ends unassigned"),)))
    g.react(key+"_construct","Prepare source cDNA sublibrary before indexing",StepKind.MAGIC,((amplified,InputRole.TEMPLATE),(kit,InputRole.CONDITION)),(ligated,),reason="Historical source-named fragmentation/endrepair/adapter-ligation and recovery are bounded; no generic completed library is substituted for precursor.")
    p5=g.input(key+"_p5","source partialP5 primer",MaterialKind.OLIGO,state=RNA_PARTIAL_P5)
    i7=g.input(key+"_i7","source representative i7 libraryprimer",MaterialKind.OLIGO,state=Construct([Segment("P7",il.P7),Segment("sample i7",I7_VALUE,feature=LIBRARY_I7),Segment("bindingarm",RNA_LIBRARY_I7[len(il.P7)+10:])] ))
    lp=g.input(key+"_lib_pol","source Chromiumlibrary PCR reagents; formulation/cycle count unassigned",MaterialKind.CATALYST)
    library=g.product(key+"_library",final.name,MaterialKind.LIBRARY,state=final)
    g.react(key+"_library_pcr","Complete source RNA-library PCR",StepKind.PCR,((ligated,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),(i7,InputRole.STRUCTURAL),(lp,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(ligated.state,reason="Actual pre-index adapter-ligated fragments retained; source primer-derived endpoint does not establish native first-cycle assembly."))
    name="RNA-only library" if key=="rna" else "Multiome RNA library"
    g.sequence(library,SequencingPlan(library=final,primers=RP,acquired_roles=("Read 1",),read_lengths={"Read 1":28},unplaced_reads=tuple(x for x in READ_SPECS[name] if x.read!="Read 1"),limitation=SEQUENCING_UNAVAILABLE))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    starting=g.input("starting","source sample configured as cells(RNA-only) or nuclei(Multiome)",MaterialKind.BIOLOGICAL,state=carrier("input"),primary=True)
    rna_start=g.product("rna_input","sample assigned toRNA-only recipe",MaterialKind.MOLECULAR,state=starting.state)
    multi_start=g.product("multi_input","sample assigned topaired Multiome recipe",MaterialKind.MOLECULAR,state=starting.state)
    g.react("assay_routes","Represent configured source assay routes; fixation recipes remain separate",StepKind.SPLIT,((starting,InputRole.TEMPLATE),),(rna_start,multi_start),changes_molecule=False)
    for key,atac in (("rna",False),("multi",True)):
        sample=rna_start if not atac else multi_start
        fix=g.input(key+"_fix_condition","methanol−20°C30min" if not atac else "0.3%formaldehyde onice10min; sourcewash",MaterialKind.CONDITION)
        fixed=g.product(key+"_fixed","source-fixed carriers",MaterialKind.MOLECULAR,state=sample.state)
        g.react(key+"_fix","Fix configured source carriers",StepKind.CROSS_LINKING if atac else StepKind.CHEMICAL_CONVERSION,((sample,InputRole.TEMPLATE),(fix,InputRole.CONDITION)),(fixed,))
        current=fixed
        if atac:
            adapter=g.input("multi_transposome","TDE1 kit20034197 loadedadapter; native composition unassigned",MaterialKind.ADAPTER)
            tn5=g.input("multi_tn5","3µlTDE1+7µlATACbufferB/12000nuclei;37°C1h",MaterialKind.CATALYST)
            tagged=g.product("multi_tagged","fixed nuclei retaining tagged chromatin and RNA",MaterialKind.MOLECULAR,state=carrier("transposed",atac=True))
            other=g.product("multi_tn5_other","other transposition end classes",MaterialKind.MOLECULAR)
            g.react("multi_tagment","Bulk transpose before GEM formation",StepKind.TAGMENTATION,((fixed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,other))
            g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,"Amplification compatibility does not establish physical removal; exact recovery through nuclei pelleting unassigned."))
            current=clean(g,"multi_tagment_recovery",tagged,"source500g5min remove mostsupernatant; retain15µl nuclei")
        condition=g.input(key+"_overload_rule","source overloaded Chromium GEMgeneration; multiplefixedcarriers may share dropletbarcode",MaterialKind.CONDITION)
        bead=g.input(key+"_bead","ChromiumMultiomeRNA/ATAC bead; exactOAK-era revision unassigned" if atac else "source Chromium3′v3 barcodebead; complete native geometry unassigned",MaterialKind.OLIGO)
        gem=g.product(key+"_gem","overloaded GEM independent carrier/molecule populations",MaterialKind.MOLECULAR,state=carrier("GEM",atac=atac))
        other=g.product(key+"_gem_other","other GEM occupancy outcomes",MaterialKind.MOLECULAR)
        g.react(key+"_partition","Overload fixed carriers into GEMs",StepKind.PARTITION,((current,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(condition,InputRole.CONDITION)),(gem,other),changes_molecule=False)
        g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,"Occupancy/cell-QC classification is not physical discard."))
        rt=g.input(key+"_gem_reagents","source standard MultiomeGEMbarcoding/quench recipe" if atac else "source in-GEMRT53°C45min",MaterialKind.CONDITION if atac else MaterialKind.CATALYST)
        copied=g.product(key+"_copied","barcoded cDNA beside independent chromatin populations",MaterialKind.MOLECULAR,state=carrier("GEMbarcoded",copied=True,atac=atac))
        g.react(key+"_gem_barcode","Barcode configured GEM analytes",StepKind.MAGIC if atac else StepKind.REVERSE_TRANSCRIPTION,((gem,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CONDITION if atac else InputRole.CATALYST)),(copied,),**({"reason":"Historical Multiome RNA/ATAC transfer/quench geometry is source-delegated; no later kit revision substituted."} if atac else {}))
        recovered=recovery(g,key,copied)
        if not atac:rna_finish(g,"rna",recovered,RNA)
        else:
            pre=g.input("multi_pre_primers","4µlPN20002714 preampprimer mix; species unassigned",MaterialKind.OLIGO)
            pre_other=g.input("multi_pre_other","same PN20002714 commercial mix opposite binding role; species unassigned",MaterialKind.OLIGO)
            pol=g.input("multi_pre_pol","NEBNextHighFidelity2×50µl/100µl reaction;10cycles",MaterialKind.CATALYST)
            amp=g.product("multi_preamp","pre-amplified mixed cDNA/ATAC populations",MaterialKind.MOLECULAR,state=carrier("preamplified mixed products; native duplex ends unassigned",copied=True,atac=True,duplex=True))
            g.react("multi_preamp_pcr","Pre-amplify both modalities10cycles",StepKind.PCR,((recovered,InputRole.TEMPLATE),(pre,InputRole.STRUCTURAL),(pre_other,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(recovered.state,reason="Named commercial primer mix on actual mixed ss-cDNA/ATAC inputs; species/first-cycle geometry unassigned. One mix does not establish one molecular primer species."))
            amp=clean(g,"multi_preamp_cleanup",amp,"1.6×SPRI;40µlEB")
            ra=g.product("multi_rna_aliquot","RNA preparation aliquot; mixed cDNA/ATAC retained",MaterialKind.MOLECULAR,state=amp.state)
            aa=g.product("multi_atac_aliquot","10µlATAC preparation aliquot; mixed cDNA/ATAC retained",MaterialKind.MOLECULAR,state=amp.state)
            g.react("multi_split","Split preamplified mixed material for library preparation",StepKind.SPLIT,((amp,InputRole.TEMPLATE),),(ra,aa),changes_molecule=False)
            rna_finish(g,"multi_rna",ra,MULTI_RNA)
            p5=g.input("atac_p5","source partialP5 primer",MaterialKind.OLIGO,state=RNA_PARTIAL_P5)
            ip=g.input("atac_i7","sampleIndexN PN1000212; species/physical width unassigned",MaterialKind.OLIGO)
            pol2=g.input("atac_pol","NEBNextHighFidelity;98/67/72°C, nsource-recommended cycles",MaterialKind.CATALYST)
            ap=g.product("atac_amplified","ATAC indexed candidate plus carried background",MaterialKind.MOLECULAR,state=ATAC,inferred=True)
            g.react("atac_pcr","Index ATAC preparation aliquot",StepKind.PCR,((aa,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),(ip,InputRole.STRUCTURAL),(pol2,InputRole.CATALYST)),(ap,),geometry=PCRGeometry.unavailable(aa.state,reason="Actual mixed preamplified aliquot enters ATAC PCR; complete SampleIndexN/native ATAC assembly remains unassigned."))
            ap=clean(g,"atac_size",ap,"source double-sidedSPRI selection; precise kitrevision/cutoffs unassigned",library=True)
            # The physically selected material is the library endpoint.
            g.sequence(ap,SequencingPlan(library=ATAC,acquired_roles=tuple(k for k,_ in ATAC_READS),read_lengths=dict(ATAC_READS),limitation=SEQUENCING_UNAVAILABLE))
    g.annotate(Highlight("rna_gem_barcode","RNA-template fate unassigned","RT output tracks cDNA beside DNA; source does not quantitatively assign nativeRNA recovery/degradation in this model."))
    g.declare_reviewed_scope(REVIEWED_TARGETS,SOURCE)
    return g.require_publishable()

REVIEWED_TARGETS=("starting","assay_routes","rna_fix_condition","rna_fix","rna_overload_rule","rna_bead","rna_partition","rna_gem_reagents","rna_gem_barcode","rna_recovery_selection","rna_recovery","rna_aliquot_rule","rna_redistribute","rna_heat","rna_release","rna_silane_selection","rna_silane","rna_i5","rna_tso","rna_cdna_pol","rna_cdna_pcr","rna_kit","rna_construct","rna_p5","rna_i7","rna_lib_pol","rna_library_pcr","rna_library","multi_input","multi_fix_condition","multi_fix","multi_transposome","multi_tn5","multi_tagment","multi_tagment_recovery_selection","multi_tagment_recovery","multi_overload_rule","multi_bead","multi_partition","multi_gem_reagents","multi_gem_barcode","multi_recovery_selection","multi_recovery","multi_aliquot_rule","multi_redistribute","multi_heat","multi_release","multi_silane_selection","multi_silane","multi_pre_primers","multi_pre_other","multi_pre_pol","multi_preamp_pcr","multi_preamp_cleanup_selection","multi_preamp_cleanup","multi_split","multi_rna_i5","multi_rna_tso","multi_rna_cdna_pol","multi_rna_cdna_pcr","multi_rna_kit","multi_rna_construct","multi_rna_p5","multi_rna_i7","multi_rna_lib_pol","multi_rna_library_pcr","multi_rna_library","atac_p5","atac_i7","atac_pol","atac_pcr","atac_size_selection","atac_size","atac_size_clean")
