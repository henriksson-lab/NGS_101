"""Private source-qualified partial preparation; no complete molecular/absence claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
class B:
    def __init__(self,title,authority,label,kind=K.MOLECULAR):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',label,kind,primary=True)
    def stock(self,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cat(self,key,label):return self.stock(key,label,K.CATALYST,R.CATALYST)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind=S.MAGIC,inputs=(),product='Actual source preparation; whole occupied substrates/products/side products and conversion unassembled',**kw):
        y=self.g.product(key+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Positive source operation; whole occupied molecules, end distributions, detailed molecular products and conversion unassembled, not a source-absence claim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,inputs=()):
        y=self.g.product(key+'_kept','Actual retained physical source fraction, carryover/recovery unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual excluded workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(inputs),(y,w));self.g.discard(w,'Actual source workup fraction excluded; no quantitative purification claim.');self.x=y;return y
    def aliquot(self,key,label,kept='Selected source aliquot',rest='Untransferred stock, amount/fate unassigned',qc=False):
        y=self.g.product(key+'_selected',kept,K.MOLECULAR);r=self.g.product(key+'_remainder',rest,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r))
        if qc:self.g.mark_terminal(Terminal(y.id,Disposition.QC,kept));self.x=r;return r
        self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,rest));self.x=y;return y




from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from illumina import TRUSEQ_READ1,TRUSEQ_READ2
AUTHORITY='PublishedPMC10864177 Strategy/Methods186–246/ExtendedData1selectedRTaseandbiotinNTcontext pluspublishedSI1MOESM1actualTables1–3p3. Standard40kHepG2PTBP1branch; noTn5/CLIP/lowinput/tissue/controlrecipeinheritance. Twooriginals; partialnotabsenceaudit.'
RT_ARM=TRUSEQ_READ2[-21:];ADAPTER_FIXED=revcomp(TRUSEQ_READ1[-18:])+'GT'
SOURCE_ORDERS={'rt':'AGACGTGTGCTCTTCCGATCTNNNNNNNNNN','adapter':'5′Phos-NNNNNNNNAGATCGGAAGAGCGTCGTGT-3′SpC3'}
NT_COMPOSITION=(('biotin_dutp','.05mMbiotin-16-dUTP; DNAdeoxyuridinenucleotide, notRNAUTP'),('biotin_dctp','.05mMbiotin-16-dCTP'),('dttp','.05mMdTTP'),('dctp','.05mMdCTP'),('datp','.1mMdATP'),('dgtp','.1mMdGTP'))
def oligo_scene(name):
    s=Scene()
    if name=='rt':s.strand('oligo',[Segment('source_RT_handle',RT_ARM),Segment('random_priming_N10','N'*10,note='Source10ntrandompriming; notassignedUMI/cellbarcode')])
    else:
        s.strand('oligo',[Segment('adapter_UMI','N'*8,feature=feature('artr.molecule','umi','random',note='Mainprimarydata-processing8ntUMIdeduplication; notRTN10randomprimer')),Segment('source_adapter_fixed',ADAPTER_FIXED)],mod5='P',mod3='SpC3');s.ligation_end('oligo','adapter_UMI',kind='single-strand',polarity='5-prime',owner='adapter',text='KnownsourcecDNAadapter5P/N8donor; SpC3isnonnucleotide3block. WholeoccupiedcDNA/sealedjunctionunassembled.')
    return s

def antibody_stock(b,key,label):
    return b.stock(key,label+' inactualblockingbuffer1mgmlBSA/.2UulRNaseOUT; volumeunassigned',chemical_state={'actual_molecular_components':['Sourceantibodyprotein','SourceBSAcarrier1mgml','SourceRNaseOUTprotein.2Uul'],'packaging':'Actualsourcepreparedantibodydilution; notfreshfreeBSAaddedtwice'})

def protocol_graph():
    b=B('ARTR-seq selectedstandardHepG2/PTBP1 antibody-directedin-situRT preparation',AUTHORITY,'Actualsource40kHepG2cells/standardPTBP1preparation; upstreamcultureunassembled',K.BIOLOGICAL);g=b.g
    g.add_source_item(SourceItem('proximity_scope','MainStrategy/ResultsARTRusesantibodylocalizedpAGRTtoidentifydirectorpotentiallynearbyindirectRNA.','NotallreadsprovedirectcovalentRBPcontacts. NoUV/IP/RTandTagoligodT/Tn5constructororinputcontrolperfectspecificityinferred.'))
    g.add_source_item(SourceItem('identifier_scope','MainstandardRTprimerhasN10randompriming; cDNAadapterN8isUMIperprimarydataprocessing.','Optionaltissue/lowinputRTprimeradds8ntsamplebarcode; notinheritedasbaselinecellbarcode/UMI. SI1qPCRprimers/siRNAsnotlibraryPCRprimers.'))
    b.op('fix','Source1.5percentPFA10minRT',inputs=(b.stock('pfa','Actual1.5percentPFAfixationchemicalcontribution'),b.cond('fix_condition','Actualsource10minRT')),product='Sourcefixedcell/RNPclasses; wholecrosslinks/occupancyunassembled')
    b.op('glycine','Source125mMglycine5minRT',S.QUENCHING,(b.stock('glycine_reagent','Actualglycinequenchchemicalcontribution'),b.cond('glycine_conditions','Actualsource125mM/5minRT')));b.clean('fix_wash','Source2xDPBS3minRT')
    b.op('permeabilize','Source.5percentTritonDPBSice10min',S.PERMEABILIZATION,(b.cond('triton','Actualsourcepermeabilizationbuffer'),));b.clean('permeabilization_wash','Source2xDPBS3minRT')
    b.op('block','Sourceblockingbuffer30minRT',inputs=(b.stock('blocking_buffer','Actual1mgmlBSA/.2UulRNaseOUTinDPBSblockingbuffer',chemical_state={'source_molecular_carrier':'BSAprotein1mgml','source_inhibitor':'RNaseOUTprotein.2Uul'}),))
    b.op('primary','SourcePTBP1primaryantibody1hRT; IFmanufacturer-guideddilutionor1to200fallback',inputs=(antibody_stock(b,'primary_ab','ActualPTBP1primaryantibodydilution; executedratioaslocallyguidedorfallback, notassumedalways1to200'),));b.clean('primary_wash','Source3xDPBS3minRTafterprimarystain')
    b.op('secondary','Sourcefluorophoresecondary1to500/30minRT',inputs=(antibody_stock(b,'secondary_ab','Actualfluorophore-labeledsecondaryantibody1to500'),));b.clean('secondary_wash','Source3xDPBS3minRTaftersecondarystain')
    b.op('localize_rt','Source10nMpAGRTaseinblockingbuffer30minRT',inputs=(b.stock('fusion_stock','Actualprepared10nMpAG-engineeredMMLV25to497/30aalinkerprotein in1mgmlBSA/.2UulRNaseOUTblockingbuffer',chemical_state={'enzyme_provenance':'Purifiedsourcefusionstock, geneticconstruction/purificationoutsidegraph','source_molecular_carrier':'BSA1mgmlwithRNaseOUTproteinintheactualdiluent'}),),product='SourceantibodylocalizedpAGRTase/fixedcellcomplexes; unboundandtargetoccupancyunassembled');b.clean('fusion_wash','Source3xDPBS3minRTafterpAGRTasestain')
    template=g.product('localized_cell_template','Templatecomponentofactualwashedcells/RNA/antibodycomplexes',K.MOLECULAR);retained=g.product('retained_pag_rtase','SameactualpAGRTasealreadylocalized/washed, amountretainedunassigned',K.CATALYST,chemical_state={'provenance':'Producedcomponentofwashedlocalize_rtpreparation, notfreshRTenzymeaddition','occupied_antibody_linkage':'Unassembledsourceproximitytargetingcomplex'})
    g.react('localized_components','Same-cellcomponentbookkeeping, nophysicalRTaseisolationordissociation',S.GENERIC,((b.x,R.TEMPLATE),),(template,retained),changes_molecule=False);b.x=template
    primer=g.input('rt_primer','Actual2uMadapterN10RTprimer',K.OLIGO,state=oligo_scene('rt'),chemical_state={'source_ordering':SOURCE_ORDERS['rt'],'N10meaning':'Randomprimingsequence, notUMI'})
    nt_inputs=[]
    for key,label in NT_COMPOSITION:nt_inputs.append(b.stock('rt_'+key,'Actualsource '+label))
    rnaseout=g.input('rt_inhibitor','ActualRTmixtureRNaseOUT1Uul',K.MOLECULAR);master=g.product('prepared_rt_mix','ActualRTmixturepreparedin50ulbuffer; finalvolumeconvention/stockvolumesunassigned',K.MOLECULAR)
    g.react('prepare_rt_mix','SourceRTmixturepreparedin50ulbufferbeforeadditiontocells',S.GENERIC,((primer,R.TEMPLATE),)+tuple(nt_inputs)+((rnaseout,R.STRUCTURAL),b.cond('rt_mix_buffer','Source50ulDPBSbuffer+3mMMg; exactfinalvolumeconventionunassigned')),(master,),changes_molecule=False)
    rp=g.product('retained_rt_primer','SameactualN10primercomponentofpreparedRTmix',K.OLIGO,state=oligo_scene('rt'),chemical_state={'source_ordering':SOURCE_ORDERS['rt'],'N10meaning':'RandomprimingnotUMI'});nt=g.product('rt_mix_nt','SameactualsixnucleotideclassesofRTmix',K.MOLECULAR,chemical_state={'source_nucleotides':dict(NT_COMPOSITION),'scope':'Biotin-dU/dCchemicalclasses, notallRNAorallDNAresiduesmodified'});other=g.product('rt_mix_other','OtheractualRTmixRNaseOUT/DPBS/Mgcomponents',K.MOLECULAR)
    g.react('rt_mix_components','Sameactualpreparedmixcomponents, nofreshduplicateprimer/NTadditions',S.GENERIC,((master,R.TEMPLATE),),(rp,nt,other),changes_molecule=False)
    rt=g.product('rt_product','Sourcebiotin-bearingcDNA/RNAhybridclassesinsidefixedcells',K.MOLECULAR,chemical_state={'polymer_class':'SourceDNAcontainsbiotin-16-deoxyuridineandbiotin-16-dCsubsetalongwithregulardeoxynucleotides','chemical_vs_readout':'DeoxyuridineisnotRNAuridine; exactpositionsandpostPCRbase/strandcopiesunassembled','primer':'ActualexternalN10randompriming, sourceN10notUMI'})
    g.react('rt','Sourcein-situRT37C30min',S.REVERSE_TRANSCRIPTION,((b.x,R.TEMPLATE),(retained,R.CATALYST),(rp,R.STRUCTURAL),(nt,R.STRUCTURAL),(other,R.STRUCTURAL)),(rt,));b.x=rt
    b.op('rt_stop','Source20mMEDTA/10mMEGTA3minRT',S.QUENCHING,(b.cond('rt_stop_conditions','Actualsourcechelationstop'),))
    b.op('imaging','Sourceanti-biotinAlexa4881to2001h/Hoechst1ugml15min/imaging',inputs=(b.stock('imaging_antibody','Actualanti-biotinBK1/39Alexa488conjugateprotein'),b.stock('hoechst','ActualHoechst33342moleculardyecontribution'),b.cond('imaging_conditions','Sourcewhole-sampleconfocalimaging, notadestructiveQCaliquot')),product='Actualimagedcells/cDNA/RNA/fluorescentantibody/dyepreparation; carryoverunassigned')
    b.op('pk','Source1mgmlPK37C2h',S.DEPROTEINIZATION,(b.cat('proteinase_k','ActualsourcePK1mgml'),b.cond('pk_conditions','Actualsource37C2h')))
    b.clean('phenol','Sourcephenol-chloroformpH8nucleicacidrecovery');b.clean('nucleic_precip','Sourceethanolconcentration; formulation/volumesunassigned')
    b.op('rna_digest','SourceRNaseH/RNaseAT1in50ul37C1h',S.DIGESTION,(b.cat('rnase_h','ActualRNaseHsource.2Uulwording'),b.cat('rnase_at1','Actual1to20RNaseA/T1'),b.cond('rna_digest_buffer','Source50mMTris7.5/75mMKCl/10mMMg/10mMDTT')),product='SourcecDNApreparationafterRNAdigestion; completeRNAremoval/biotinoccupancyunassembled')
    raw=g.input('c1_stock','ActualC1beadstockforpreblocking; preparationvolumeunassigned',K.MOLECULAR);blocked=g.product('preblocked_c1','ActualpreblockedC1preparationincludingmolecularblockers; carryover/occupancyunassigned',K.MOLECULAR)
    g.react('preblock_c1','SourceC1preblocking30minRT',S.MAGIC,((raw,R.TEMPLATE),b.stock('bead_bsa','Actual1ugulUltraPureBSAproteinblocker'),b.stock('bead_salmon_dna','Actual1ugulUltraPureSalmonSpermDNAblocker'),b.stock('bead_yeast_trna','Actual1ugulYeasttRNAblocker')),(blocked,),reason='Actualbeadblockingmolecules; fulladsorption/stockvolumes/unboundblockersunassembled. No source-unstatedpreblockwashinvented.')
    portion=g.product('c1_ten_ul','Actual10ulpreblockedC1portionforenrichment',K.SELECTION);rest=g.product('c1_untransferred','AnyuntransferredpreblockedC1preparation; amountmayzero/unassigned',K.MOLECULAR);g.react('c1_allocate','Actual10ulpreblockedbeadallocation',S.SPLIT,((blocked,R.TEMPLATE),),(portion,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Preblockedstockremainderamount/fateunassigned'))
    bound=g.product('c1_bound','ActualselectedC1associatedbiotin-cDNAclasses',K.MOLECULAR);sup=g.product('c1_supernatant','Actualexcludedcaptureunboundfraction',K.WASTE);g.react('capture','Sourcebiotin-cDNAenrichment20minRT',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),(portion,R.SELECTION)),(bound,sup));g.discard(sup,'Actualunboundcapturefractionexcluded; nocompleteRNAremoval/captureyieldclaim');b.x=bound
    adapter=g.input('cdna_adapter','Actualsource5uM5P/N8/SpC3cDNAadapter',K.OLIGO,state=oligo_scene('adapter'),chemical_state={'source_ordering':SOURCE_ORDERS['adapter'],'five_prime':'P','three_prime':'SpC3nonnucleotideblocker','N8meaning':'SourceUMIperdata-processingMethods'})
    b.op('adapter_ligate','Sourcebead-associatedcDNA3adapterligation25C16h',S.LIGATION,((adapter,R.STRUCTURAL),b.cat('rnl1','ActualT4RNAligase1at1Uul'),b.stock('ligation_atp','Actual1mMATP'),b.cond('ligation_conditions','Source50mMTris7.5/10mMMg/25percentPEG8000/5uMadapter')),geometry=AdapterLigationGeometry.known('single-strand',oligo_scene('adapter'),statement='Knownsourceadapter5P/N8/3SpC3localdonorend; fullcDNA3acceptor/occupiedbeadcomplexandsealedproductunassembled.'))
    free=g.product('released_cdna','Actualformamide/heat-releasedcDNApreparation',K.MOLECULAR);res=g.product('released_bead_residual','Actualresidualbeads/retainedmaterialafterrelease',K.MOLECULAR);g.react('release','Source95percentformamide/10mMEDTA8/95C10min',S.MAGIC,((b.x,R.TEMPLATE),b.cond('release_conditions','Actualsourceformamide/EDTAheat')),(free,res),reason='PositivephysicalcDNArelease; completeelutionyield/occupiedbiotinclassesunassembled.');g.mark_terminal(Terminal(res.id,Disposition.UNRESOLVED,'Residualbead/retainedmaterialfateunassigned'));b.x=free
    b.clean('cdna_precip','Sourceethanolprecipitation/dissolve10ulwater')
    pcrmaster=g.input('pcr_master','ActualNEBNextUltraIIQ5MasterMixcontribution, final1xin40ul; stockvolumeunassigned',K.MOLECULAR);enz=g.product('pcr_enzyme','SameactualQ5masterDNApolymerasecomponent',K.CATALYST);pn=g.product('pcr_nt','SameactualmasterdNTPcomponent',K.MOLECULAR);pb=g.product('pcr_buffer','Sameactualmasterbuffercomponent',K.CONDITION);g.react('pcr_components','Namedactualmastercomponents, nofreshseparateenzyme/NTadditions',S.GENERIC,((pcrmaster,R.TEMPLATE),),(enz,pn,pb),changes_molecule=False)
    primer_inventory=g.input('pcr_primer_inventory','Actual.5uMIlluminasequencingprimeraggregate supplies; sourceE7335Sisexample, stockcount/packaging/wholeorders/indexvaluesunassembled',K.MOLECULAR)
    pf=g.product('pcr_forward','ActualPCRforward-roleprimerclass; wholeorder/unassignedspeciescount',K.OLIGO);pr=g.product('pcr_reverse','ActualPCRreverse-roleprimerclass; wholeorder/unassignedspeciescount',K.OLIGO);g.react('pcr_primer_components','Aggregateactualprimersuppliesdecomposedintoforward/reversefunctionalclasses; noprovenpremixornewextraaddition',S.GENERIC,((primer_inventory,R.TEMPLATE),),(pf,pr),changes_molecule=False)
    pcr_scene=Scene();pcr_scene.strand('source_cDNA',[Segment('sourceadapterligatedcDNA','[wholeadapter-ligatedcDNAunassembled]',placeholder=True)])
    b.op('pcr','Source40ul/all10ulcDNA/13PCRcycles',S.PCR,inputs=((pf,R.STRUCTURAL),(pr,R.STRUCTURAL),(enz,R.CATALYST),(pn,R.STRUCTURAL),(pb,R.CONDITION),b.cond('pcr_conditions','Source98C30s;13x98C10s/60C30s/65C45s;65C5min')),geometry=PCRGeometry.unavailable(pcr_scene,reason='SourceactualPCRwithtwofunctionalprimerclasses; wholehistoricalprimerorders/landing/occupiedtemplateandampliconunassembled. E7335Sisexample, noconstructororabsenceinherited.'),product='Sourceamplifiedlibraryclasses; originalbiotin-strands/copiedDNA/readlayoutunassembled')
    b.op('library_gel','Source6percentNovexTBEgel180to400bpselection',inputs=(b.cond('gel_conditions','Actualsource180to400bpbandwindow'),))
    b.aliquot('library_band','Actual180to400bpbandexcision','Actualselected180to400bplibrarygelpiece','Unselectedgelbands/material; fateunassigned')
    b.op('gel_recovery','Sourcegelpurifiedsequencinglibrary; internalextractionconditionsunassembled',inputs=(b.stock('gel_supplies','Actualgel-purificationreagents/physicalworkupinventory; identities/amountsunassembled'),),product='Sourcegel-purifiedlibrarypreparation; fullmolecules/recoveryunassembled')
    b.aliquot('loading','Actualsequencingloadingallocation; amountunassigned')
    lib=g.product('library','SourceARTRsequencinglibraryclass; exactfullconstruct/readlandingunassembled',K.LIBRARY);g.react('load','SourceNextSeq550/2000orNovaSeq6000',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceplatformsNextSeq550/2000/NovaSeq6000; exactstandardrunlengths/fullPCRprimer/strand/index/readlandingunassembled. Optionalbarcodedlowinput/tissueRead2processingnotautomaticallyassignedtobaselinerun.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('SourcewholeRTprimerretained',''.join(x.top for x in oligo_scene('rt').strands['oligo'].segs),SOURCE_ORDERS['rt'])
    c('AdapterSpC3notextraneousSbase',''.join(x.top for x in oligo_scene('adapter').strands['oligo'].segs),'N'*8+'AGATCGGAAGAGCGTCGTGT')
    c('N10randomprimerisnotUMI',oligo_scene('rt').strands['oligo'].segs[1].feature,None)
    c('N8adapterhasUMImeaning',oligo_scene('adapter').strands['oligo'].segs[0].feature.role,'umi')
    c('RTconsumesproducedretainedfusion',g.steps['rt'].inputs[1].material,'retained_pag_rtase')
    c('Fusioncomesfromwashedtargetingcomplex',g.steps['localized_components'].inputs[0].material,'fusion_wash_kept')
    c('ActualsixNTclassesaddedoncebeforeRTmix',len(g.steps['prepare_rt_mix'].inputs[1:7]),6)
    c('RTusespreparedmixprimercomponent',g.steps['rt'].inputs[2].material,'retained_rt_primer')
    c('BeadpreblockcontainsphysicalDNAandtRNA',tuple(x.material for x in g.steps['preblock_c1'].inputs[2:]),('bead_salmon_dna','bead_yeast_trna'))
    c('Captureusesproducedpreblocked10ulportion',g.steps['capture'].inputs[1].material,'c1_ten_ul')
    c('PKprecedesRNAdigestion',g.steps['pk'].kind,S.DEPROTEINIZATION)
    c('Releasehasactualresidualbeadfraction','released_bead_residual'in g.terminals,True)
    c('PCRconsumesbothproducedfunctionalprimerclasses',tuple(x.material for x in g.steps['pcr'].inputs[1:3]),('pcr_forward','pcr_reverse'))
    c('PCRtypedwithunassembledlandinggeometry',g.steps['pcr'].kind,S.PCR)
    c('PCRuseswhole10ulrecoveredsourcecDNA',g.steps['pcr'].inputs[0].material,'cdna_precip_kept')
    c('Loadingusesselectedaliquotonly',g.steps['load'].inputs[0].material,'loading_selected')
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41592-023-02146-w_PMC10864177.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41592-023-02146-w')):
        t=p.read_text();flat=''.join(t.split())
        c('WholeoriginalRTprimer',SOURCE_ORDERS['rt'] in t,True);c('Wholeoriginalmodifiedadapter',SOURCE_ORDERS['adapter'] in t,True)
        c('ActualRNAdigestionbeforebiotincDNAenrichment',t.index('RNA was digested')<t.index('followed by biotinylated cDNA enrichment'),True)
        c('ActualbiotinDUTPsource','.05mMbiotin-16-dUTP'in flat or '0.05mMbiotin-16-dUTP'in flat,True)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
