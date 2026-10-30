from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture including incomplete/unreacted/off-target classes; occupancy/recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; complete occupied molecular inputs/products and historical formulation unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,extra=()):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(extra),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder));self.x=y;return y

ORDERS={'CD151': 'GCTGATGTAGTCACTCTTGA',
 'p5_index': 'AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATC*T',
 'p7_index': '/5’Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNCGATGTATCTCGTATGCCGTCTTCTGCTTG/3’AmMO/',
 'rc_p5': 'GATCGGAAGAGCGTCGTG',
 'rc_p7': 'GCTCTTCCGATC*T'}
from illumina import STEM_COMPLEMENT,P5,TRUSEQ_READ1
AUTHORITY='PublishedPMC8654851 mainMethods104–110/Fig1actualPDFp2; publishedMOESM2SupplementaryProtocol0–H actualadaptertablep1; originalTableS1MOESM3ZIPXML. Threeoriginalsclosed, no preprintorendpointrecipeinheritance.'
def bases(name):
    import re
    return re.sub('/[^/]+/','',ORDERS[name]).replace('*','')
def single(name):
    s=Scene()
    if name=='CD151':
        s.strand('source',[Segment('Source20nttargetingregion',ORDERS[name].replace('T','U')),Segment('WholefunctionalRNAguidebodyunassembled','[unassembledguidescaffold]',placeholder=True)]);s.note('source','SourceTableS1RNAguideentry20bases is NOTwhole100bsgRNA; no canonicalscaffold imported.');return s
    seq=bases(name);segs=[]
    if 'N'in seq:
        n=6 if name=='p7_index'else 8;a,c=seq.split('N'*n);segs=[Segment('sourceprefix',a),Segment('Sourceindexwindow','N'*n,placeholder=True,feature=feature('rgen.'+name,'sample_index','unknown',note='Sourceindexadapterwindow, actualstudyindexvalues/readorientationunassigned; NOTUMI')),Segment('sourcesuffix',c)]
    else:segs=[Segment('sourcewholeorder',seq)]
    s.strand(name,segs,mod5='P'if name=='p7_index'else('OH — sourceFig1dephosphorylatedp5L'if name=='rc_p5'else''),mod3='AmMO — source3primeaminomodifier'if name=='p7_index'else'')
    if '*'in ORDERS[name]:s.note(name,'Sourcephosphorothioatebondbeforeterminal3primeT; rawwholeorder '+ORDERS[name])
    return s

def adapter(which):
    s=Scene()
    if which=='p7':
        seq=bases('p7_index');stem=revcomp(STEM_COMPLEMENT);assert seq.startswith(stem);tail=seq[len(stem):];left,right=tail.split('N'*6)
        s.strand('p7_index',[Segment('junctionstem',stem),Segment('sourcearm',left),Segment('p7index','N'*6,placeholder=True,feature=feature('rgen.p7index','sample_index','unknown',note='SourceN6indexwindow; notUMI')),Segment('sourcedistalarm',right)],mod5='P',mod3='AmMO')
        s.anneal('rc_p7',[Segment('stemcomplement',STEM_COMPLEMENT),Segment('3primeT','T')],to='p7_index',pair=('stemcomplement','junctionstem'),unpaired=('3primeT',))
        s.ligation_end('rc_p7','3primeT',kind='cohesive',polarity='3-prime',owner='adapter',text='Sourceprotected3T/PSbond, localadapterendonly')
        s.ligation_end('p7_index','junctionstem',kind='cohesive',polarity='5-prime',owner='adapter',text='Source5P atligatingadapterend; distal3AmMOretained')
        s.note('rc_p7','TerminalTjoinedthroughsourcePSbond; doesnotassertwholegenomicjunctionoccupancy.')
    elif which=='p5':
        seq=P5+'N'*8+TRUSEQ_READ1;assert seq==bases('p5_index');stem=revcomp(bases('rc_p5'));assert seq.endswith(stem+'T');prefix=seq[:-len(stem)-1];a,c=prefix.split('N'*8)
        s.strand('p5_index',[Segment('P5arm',a),Segment('p5index','N'*8,placeholder=True,feature=feature('rgen.p5index','sample_index','unknown',note='SourceN8indexwindow; actualstudyindexuse/readorientationunassigned')),Segment('sourcearm',c),Segment('junctionstem',stem),Segment('3primeT','T')])
        s.anneal('rc_p5',[Segment('stemcomplement',bases('rc_p5'))],to='p5_index',pair=('stemcomplement','junctionstem'),mod5='OH — sourceFig1dephosphorylatedp5L')
        s.ligation_end('p5_index','3primeT',kind='cohesive',polarity='3-prime',owner='adapter',text='SourcePSprotected3T; oppositeadapter5OHnotphosphorylated')
        s.note('p5_index','SourcephosphorothioatebondbeforeterminalT; partialadapter-readyinterface, notsealedgenomicproduct.')
    else:raise ValueError('Sourcep7orp5adapteronly')
    return s

def allocated(g,source,key,label,remainder):
    y=g.product(key+'_selected',label,source.kind,state=source.state,chemical_state=source.chemical_state);r=g.product(key+'_remainder',remainder,source.kind,state=source.state,chemical_state=source.chemical_state);g.react(key,label,S.SPLIT,((source,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder));return y

def adapter_stock(g,which):
    names=('p7_index','rc_p7')if which=='p7'else('p5_index','rc_p5');xs=[g.input(n,ORDERS[n],K.OLIGO,state=single(n),chemical_state={'whole_source_ordering':ORDERS[n],'phosphorothioate':'Source*beforeterminal3T whereprinted','source_end5':'p7_index5P;rc_p5dephosphorylated5OH;otherendsnotinferred'})for n in names];out=g.product(which+'_annealed','Actualsourceannealed'+which+'Ladapterstock',K.ADAPTER,state=adapter(which),chemical_state={'source_preparation':'1ul100uMeacholigo;95C3min/slowcool','stock_volume_ul':66.5 if which=='p7'else 6.6,'duplex_buffer_ul':64.5 if which=='p7'else 4.6,'source_concentration_uM':1.5 if which=='p7'else 15})
    buf=g.input(which+'_duplex_buffer','ActualsourceIDTduplexbufferallocation',K.CONDITION);g.react(which+'_anneal','Sourceadapterstockannealing95C3min/slowcool',S.HYBRIDIZATION,((xs[0],R.TEMPLATE),(xs[1],R.STRUCTURAL),(buf,R.CONDITION)),(out,));return allocated(g,out,which+'_adapter_aliquot','Actual2.5ulsourceadapterstockallocation','Actualuntransferredannealedadapterstock')
def clean(b,key,label):b.clean(key,label,extra=(b.stock(key+'_beads','ActualsourceSPRIselectallocation '+label,K.SELECTION,R.SELECTION),))
def qc(b,key,label):
    g=b.g;y=g.product(key+'_analytical',label,K.MOLECULAR);r=g.product(key+'_stock','ActualsourcepreparationremainingafterQC; amountsunassigned',K.MOLECULAR);g.react(key,label,S.SPLIT,((b.x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(y.id,Disposition.QC,label));b.x=r
def protocol_graph():
    b=B('RGEN-seq selectedSupplementaryProtocol enzymaticFSreaction',AUTHORITY,'Sourceoptionalblockingroute: selected3uggDNAin38ulTE; postcleanupmass/recoveryunassigned');g=b.g
    b.op('block','SourceoptionalCIP/TdT/ddATPblocking37C60min/75C20min',S.MAGIC,(b.stock('cip','Actual2ulQuickCIP'),b.stock('tdt','Actual4.5ulTerminalTransferase'),b.stock('ddatp','Actual.5ul10mMddATPterminatingnucleotides',K.MOLECULAR,R.STRUCTURAL),b.cond('block_buffer','Actual5ul10xCutSmart,total50ul')),reason='Sourceblocks pre-existingrandomends; fullphosphate/enddistribution/ddAaddition and residualenzymeclassesunassembled. Not claimingperfectcompleteblocking.')
    clean(b,'block_clean','Add50ulTE/80ulSPRIselect/350ul85percentEtOH/elute39ulTE');b.split('block_transfer','Retain37ulof39uleluate','Actual37ulblockedgDNAstock','Untransferred2ulblockeluate')
    guide=g.input('guide','ActualsourcehumanCD151full100bsgRNApreparation; TableS1onlytarget20ntspecified',K.OLIGO,state=single('CD151'),chemical_state={'source_targeting_region_as_printed':ORDERS['CD151'],'polymer':'RNA; sourceTnotationprojectedtoU forknownspacer','whole_sgRNA':'Source100bfullsizeguide, scaffold/endmods/wholeorderunassembled; NOT20ntwholeguide'})
    cas=g.input('cas9','ActualS.pyogenesCas9allocation',K.CATALYST);rnp=g.product('rnp','ActualsourcepreincubatedRNPpreparation; bindingstoichiometry/occupancyunassembled',K.MOLECULAR);buf=g.input('rnp_buffer','Actual10ul1xCutSmart;sgRNA:Cas9 1.5:1,RT10min',K.CONDITION);g.react('rnp_assemble','SourceRNPpreincubation',S.MAGIC,((guide,R.TEMPLATE),(cas,R.CATALYST),(buf,R.CONDITION)),(rnp,),reason='Positive100bguide/Cas9RNPassembly; fullguidescaffold/occupiedcomplexunassembled. CD151targetingsequenceisnotfullfunctionalguide.')
    b.op('rgen_cut','AddpreincubatedRNPto40ulgDNA/CutSmart37C60min',S.MAGIC,((rnp,R.STRUCTURAL),b.cond('cut_conditions','SourceA6specifies1–4uggDNAin40ul1xCutSmart, whereasoptionalblockretains37ulTEeluate. RecoveredDNAmass/aliquot-volume/preparationinterface/finalbufferconventionunassembled; NOTcalculated3ulbufferadditionornewgDNAstock. SourceRGEN16–1024nMstudyrange; selecteddoseunassigned. No guessedexacttargetends.')),product='Sourcecut/uncut/blockedpre-existingDNAclasses; fullgenomiccutsiteendsoffproductsunassembled')
    b.op('stop','Source65C10mintermination',S.QUENCHING,(b.cond('stop_conditions','Actual65C10min'),))
    b.op('rnase','Add.5ulRNasecocktail37C30min',S.DIGESTION,(b.stock('rnase_enzyme','Actual.5ulRNasecocktail'),))
    b.op('proteinase','Add2ulProteinaseK50C30min',S.DEPROTEINIZATION,(b.stock('proteinase_enzyme','Actual2ulProteinaseK'),))
    clean(b,'cut_clean','Add47.5ulTEto52.5ul/70ulSPRIselect/85percentEtOH/elute33ulTE');b.split('cut_transfer','Retain31ulof33uleluate','Actual31ulcutgDNAstock','Untransferred2ulcuteluate');b.split('tail_template','Select30.25ulof31ulforKlenow','Actual30.25ulKlenowtemplate','Untransferred.75ulcutgDNAstock')
    b.op('tail1','Klenowexo-minus/dNTPprocessing37C30min/75C20min',S.MAGIC,(b.stock('klenow','Actual3ulKlenow3prime5primeexo-minusM0212'),b.stock('dntp','Actual.5ul10mMdNTPmix',K.MOLECULAR,R.STRUCTURAL),b.cond('tail1_buffer','Actual3.75ul10xCutSmart;total37.5ul')),reason='Sourceendprocessing/A-tailingoperationwithallfournucleotidepool, not onlydATP. Source-knownadapter3Tends separatelydeclared; wholeheterogeneousgenomicproductsunassembled.')
    p7=adapter_stock(g,'p7')
    b.op('ligate_p7','SourcefirstligationRT30min/total50ul',S.MAGIC,((p7,R.STRUCTURAL),b.stock('ligase1','Actual5ulNEBNextQuickT4DNAligase'),b.stock('atp1','Actual3.5ul10mMATP',K.MOLECULAR,R.STRUCTURAL),b.cond('ligase1_buffer','Actual1.5ul10xT4DNAligationbuffer;37.5uldA-tailedDNA+2.5uladapter')),geometry=AdapterLigationGeometry.known('cohesive',adapter('p7'),statement='Sourceknownlocaladapter3T/5P/PS/distal3AmMOonly; fullgenomicsubstrateoccupancy/sealedjunctionsunassembled.'))
    clean(b,'ligation1_clean','Add50ulTE/70ulSPRIselect/85percentEtOH/elute32ulTE');b.split('ligation1_transfer','Retain30ulof32uleluate','Actual30ulfirst-ligatedgDNAstock','Untransferred2ulfirstligationeluate');b.split('fs_template','Select26ul/150–500ngforFS','Actual26ulfragmentationtemplate; source150–500ngrange, exactmassunassigned','Untransferred4ulfirst-ligatedgDNAstock')
    mix=g.input('fs_mix','Actual2ulNEBNextUltraIIFSEnzymeMix',K.MOLECULAR);enz=g.product('fs_enzymes','ActualsameFS-mixfragmentation/endrepair/A-tailingcatalystcomponents; identitiesunassembled',K.CATALYST);other=g.product('fs_other','ActualsameFS-mixothercomponents; exactnucleotidecofactors/formulationunassembled',K.MOLECULAR);g.react('fs_mix_parts','Exposecopresentcomponentsofactual2ulFSmix, no freshseparateenzymes',S.GENERIC,((mix,R.TEMPLATE),),(enz,other),changes_molecule=False)
    b.op('fs','Concurrentfragmentation/endrepair/Atailing37C12min/65C30min',S.MAGIC,((enz,R.CATALYST),(other,R.STRUCTURAL),b.cond('fs_buffer','Actual7ulFSReactionBuffer/26ultemplate,total35ul')),reason='SourceconcurrentFSModuleoperation; wholefragmentends/lengths/retainedp7occupancyunassembled. Covarisalternativeroute excluded.')
    p5=adapter_stock(g,'p5');lm=g.input('ligation_master','Actual30ulNEBNextUltraIILigationMasterMixallocation',K.MOLECULAR);lig=g.product('ligation_master_enzyme','Ligasecomponentofsameactualmastermix',K.CATALYST);ls=g.product('ligation_master_supplies','Othercopresentactualmastermixbuffer/cofactor/nucleotidecomponents; exactformulationunassembled',K.MOLECULAR);g.react('ligation_master_parts','Same30ulphysicalmastermixcomponents, notextraexternaladditions',S.GENERIC,((lm,R.TEMPLATE),),(lig,ls),changes_molecule=False)
    b.op('ligate_p5','SecondligationRT20min/total68.5ul',S.MAGIC,((p5,R.STRUCTURAL),(lig,R.CATALYST),(ls,R.STRUCTURAL),b.stock('enhancer','Actual1ulNEBNextLigationEnhancer',K.MOLECULAR,R.STRUCTURAL)),geometry=AdapterLigationGeometry.known('cohesive',adapter('p5'),statement='Sourceknownlocalprotected3T/dephosphorylated5OHadapteronly. DoNOTinvent5P orclaimfullsealedgenomicjunctions.'))
    clean(b,'final_clean','Add31.5ulTE/70ulSPRIselect/350ul85percentEtOH/elute25ulTE');b.split('final_transfer','Retain23ulof25uleluate','Actual23ulfinishedlibrarypreparation; productive/unproductiveclassescoexist','Untransferred2ulfinaleluate')
    qc(b,'size_qc','ActualBioanalyzerDNA1000analyticalallocation; amountunassigned');qc(b,'quant_qc','ActualqPCRlibraryquantificationanalyticalallocation; PCRisQCnotlibraryamplification')
    b.split('sample_loading','Sourcecompletedlibraryallocationforsequencing; amountunassigned','Selectedlibraryloadingallocation','UntransferredfinishedRGENlibrarystock')
    # A single protocol reaction is represented. Actual experimental replicates/controls
    # and study pool composition are not invented from the three-replicate overview.
    productive=g.product('productive_library','SourceP5/P7cluster-producinglibraryclass; fullmoleculesunassembled',K.LIBRARY);nonproductive=g.product('nonproductive_loaded','LoadednonclusteredP5/P5/otherclass; completecomposition/recoveryunassembled',K.MOLECULAR);surface=g.input('flowcell','ActualIlluminaflowcellsurface; oligoidentity/occupancyunassembled',K.SELECTION)
    g.react('clustering','Sourceflowcellhybridization/clusteringboundary, notpreparationcleanup',S.MAGIC,((b.x,R.TEMPLATE),(surface,R.SELECTION)),(productive,nonproductive),reason='SourcepaperdescribesP5/P7cluster-producingversusP5/P5nonproductiveclassesandnonclusteredwashaway. No completeendcombination/sequenceproof, perfectspecificity or quantitativeyieldclaim.');g.discard(nonproductive,'Actualsourceflowcellnonclusteredwashfraction, NOTremovedbyearlierSPRIcleanup.')
    g.sequence(productive,SequencingPlan(read_lengths={'R1':26,'R2':151},limitation='Mainreported26/151 selected; publishedSIH25R1/100–300R2 retainedasseparatesourcevalues. Wholegenomicconstruct/indexassignment/primerlandingunassembled. R2mapping iscomputational.'))
    g.add_source_item(SourceItem('block_to_digest_interface','SI0retains37ulTEeluate;SIA6specifies1–4uggDNAin40ul1xCutSmart.','ActualrecoveredDNAamount/selection/bufferpreparationandfinalbufferconventionunassembled; no inferred3ulbufferaddition ornewDNA. Knownupstream37ulrecoveryanddownstream40ulrecipearenotacompletevolumetricbridge.'))
    g.add_source_item(SourceItem('readcount_scope','MainMethods110:26R1/151R2;SupplementaryProtocolH:25R1/100–300R2.','Differentstatedsettingsretained; no silent25→26repairorclaimidenticalrunidentity.'))
    g.add_source_item(SourceItem('physical_scope',AUTHORITY,'Protocolreactionpartialgraph; actualexperimentalreplicate/controlpoolidentityunassembled, no copiedwholegenomicproducts or sourceabsence.'));return g
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import xml.etree.ElementTree as E,zipfile
    c=Check();g=protocol_graph()
    c('PreparedguidekeepsunknownwholeRNAregion',g.materials['guide'].state.strands['source'].segs[1].placeholder,True)
    c('OptionalblockactualrecoveredeluatefeedsRGEN',g.steps['rgen_cut'].inputs[0].material,'block_transfer_selected')
    c('Klenowusesselected30.25ulof31ul',g.steps['tail1'].inputs[0].material,'tail_template_selected')
    c('FSuses26ulselectionof30ul',g.steps['fs'].inputs[0].material,'fs_template_selected')
    c('P7adapterligationconsumesactual2.5ulstockaliquot',any(p.material=='p7_adapter_aliquot_selected'for p in g.steps['ligate_p7'].inputs),True)
    c('P5adapterligationconsumesdifferentactualpreparedstock',any(p.material=='p5_adapter_aliquot_selected'for p in g.steps['ligate_p5'].inputs),True)
    c('P7distal3aminoblockerretained',adapter('p7').strands['p7_index'].mod3,'AmMO')
    c('P5oppositeendnotassumed5phosphate','OH' in adapter('p5').strands['rc_p5'].mod5,True)
    c('Bothadaptersretainknownlocal3Tends',all(any(e.polarity=='3-prime'for e in adapter(x).ligation_ends)for x in ('p7','p5')),True)
    c('SameFSstocksuppliesactualcatalystclass',g.producer_of('fs_enzymes').id,'fs_mix_parts')
    c('Same30ulligationmasterstocksuppliesligase',g.producer_of('ligation_master_enzyme').id,'ligation_master_parts')
    c('Finalretained23ulprecedesQC',g.steps['size_qc'].inputs[0].material,'final_transfer_selected')
    c('QuantificationPCRisQCfraction',g.terminals['quant_qc_analytical'].disposition,Disposition.QC)
    c('LibraryconstructionhasnoPCRstep',any(s.kind==S.PCR for s in g.steps.values()),False)
    c('Nonproductiveclassremovedatflowcellnotearlier',g.producer_of('nonproductive_loaded').id,'clustering')
    c('MainreportedrunselectedwithoutchangingSI',g.sequencing['productive_library'].read_lengths,{'R1':26,'R2':151})
    path=sources_root()/Path(__file__).resolve().parents[1].name/'s41598-021-03160-8_supp_41598_2021_3160_MOESM3_ESM.xlsx'
    if have(c,Source(path,'https://doi.org/10.1038/s41598-021-03160-8')):
        ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
        with zipfile.ZipFile(path) as z:
            ss=[''.join(x.itertext())for x in E.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)];rows={}
            for r in E.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//s:row',ns):
                vals=[]
                for cell in r.findall('s:c',ns):
                    v=cell.find('s:v',ns);s=v.text if v is not None else'';vals.append(ss[int(s)]if cell.get('t')=='s'else s)
                if len(vals)>1:rows[vals[0]]=vals[1]
            for name,value in ORDERS.items():c('OriginalTableS1wholeordering '+name,rows.get(name),value)
    c('Wholegenomiclibrarynotfabricated',g.materials['productive_library'].state,None);c('Privatepartialnotpublicationready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
