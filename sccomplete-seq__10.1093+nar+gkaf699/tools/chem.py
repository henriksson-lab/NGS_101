"""Private source-qualified partial preparation; no complete molecular/absence claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2
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
    def clean(self,key,label):
        y=self.g.product(key+'_kept','Actual retained physical source fraction, carryover/recovery unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual excluded workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label)),(y,w));self.g.discard(w,'Actual source workup fraction excluded; no quantitative purification claim.');self.x=y;return y
    def aliquot(self,key,label,kept='Selected source aliquot',rest='Untransferred stock, amount/fate unassigned',qc=False):
        y=self.g.product(key+'_selected',kept,K.MOLECULAR);r=self.g.product(key+'_remainder',rest,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r))
        if qc:self.g.mark_terminal(Terminal(y.id,Disposition.QC,kept));self.x=r;return r
        self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,rest));self.x=y;return y
    def pcr(self,key,label,forward,reverse,**kw):
        # Source one bulk preparation: produced components are copresent, not physical fractions.
        chemicals=kw.pop('primer_chemicals',{})
        stock=self.g.input(key+'_stock','Actual source prepared PCR reagent stock; quantitative internal recipe specified in step, components same tube',K.MOLECULAR)
        specs=[('f',K.OLIGO,forward,R.STRUCTURAL),('r',K.OLIGO,reverse,R.STRUCTURAL),('enzyme',K.CATALYST,None,R.CATALYST),('nt',K.MOLECULAR,None,R.STRUCTURAL),('buffer',K.CONDITION,None,R.CONDITION)]
        out=[];ins=[]
        for suffix,kind,state,role in specs:
            z=self.g.product(key+'_'+suffix,'Same prepared PCR stock '+suffix,kind,state=state,chemical_state=chemicals.get(suffix));out.append(z);ins.append((z,role))
        self.g.react(key+'_components','Copresent components of ONE source PCR stock, not extra additions',S.GENERIC,((stock,R.TEMPLATE),),tuple(out),changes_molecule=False)
        return self.op(key,label,S.PCR,ins,geometry=PCRGeometry.unavailable(unknown_scene(),reason='Whole occupied primer landing and heterogeneous amplicon products unassembled; source exact primers retained where assigned.'),**kw)
    def finish(self,key,label,plan):
        y=self.g.product(key,label,K.LIBRARY);self.g.react(key+'_load',label,S.GENERIC,((self.x,R.TEMPLATE),),(y,),changes_molecule=False);self.g.sequence(y,SequencingPlan(limitation=plan));return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('whole source preparation unassembled','[unassembled]',placeholder=True)]);return s

def oligo(seq,name='whole source oligo',mod5=None,mod3=None):
    import re
    parts=[Segment(name+' '+str(i),part,placeholder=bool(re.search('[NRYKMSWBDHV]',part))) for i,part in enumerate(re.split('([NRYKMSWBDHV]+)',seq)) if part]
    s=Scene();s.strand('oligo',parts,mod5=mod5,mod3=mod3);return s

LNA_PREFIX='UGACUGGAGUUCCUUGGCACCCGAGAAUUCCA'
LNA_ORDER='/5Biosg/'+LNA_PREFIX+'rGrG+G'
PRE_F='C'+TRUSEQ_READ1[-21:]
PRE_R='GTGACTGGAGTTCCTTGGCA'
HTO_ADD=TRUSEQ_READ2
LIB_F=P5+TRUSEQ_READ1[4:24]
SHORT_R_TAIL='GTGACTGGAGTTCCT'
LONG_R_TAIL=TRUSEQ_READ2[:21]
SOURCE_ORDERINGS={'LNA':LNA_ORDER,'preF':PRE_F,'preR':PRE_R,'HTOadd':HTO_ADD[:-2]+'*'+HTO_ADD[-2]+'*'+HTO_ADD[-1],'libF':LIB_F[:-2]+'*'+LIB_F[-2]+'*'+LIB_F[-1],'shortR':P7+'[8bpi7index]'+SHORT_R_TAIL,'HTOr':P7+'[8bpi7index]'+HTO_ADD[:-2]+'*'+HTO_ADD[-2]+'*'+HTO_ADD[-1],'longR':P7+'[8bpi7index]'+LONG_R_TAIL}
AUTHORITY='Published PMC12276002Methods62–89 plus actualSI TableS2p2/FigS1p4. PreparedPBMCbranch withRNA/ADT/HTO, source four physical25ulaliquots. No historical10xwholeconstructor inherited.'
def tso_scene():
    s=Scene();s.strand('tso',[Segment('literalsourceU-bearingbody',LNA_PREFIX),Segment('source terminalrG2','GG',note='Two explicitlyprintedrG residues'),Segment('source terminalLNA','G',note='Source +G LNA, not ordinaryDNA residue')],mod5='biotin');return s

def indexed_scene(tail):
    s=Scene();s.strand('reverse',[Segment('sourceP7',P7),Segment('source8nti7','N'*8,placeholder=True,feature=feature('sccomplete.i7','sample_index','unknown')),Segment('source right body',tail)]);return s

def protocol_graph():
    b=B('scComplete-seq selectedPBMC totalRNA/ADT/HTO preparation',AUTHORITY,'Actual preparedPBMCconditionpool aftersourceFcblock/0.5ugpermillionhashtag/0.125ugpermillionADTcoc ktail and0.02%FBSwashes; wholeconjugatesunassembled',K.BIOLOGICAL);g=b.g
    b.aliquot('cells_input','Actualsource21.5ulcellpoolselectedfor75ulGEMmix','Actual21.5ulpreparedPBMCmixture','UntransferredpreparedPBMCpool, fateunassigned')
    tso=g.input('lna_tso','Source3ul100uMLNA-TSO',K.OLIGO,state=tso_scene(),chemical_state={'source_ordering':LNA_ORDER,'five_prime':'biotin','literal_U_retained':True,'body_sugar_partition':'unassigned from literal U spelling alone','terminal_LNA':True})
    b.op('gem','Source75ulv3.1mixtureencapsulationonChipG',inputs=((tso,R.STRUCTURAL),b.stock('rt_reagent','Actual18.8ulRTReagentB molecularstock, internalNT/wholecaptureformulationunassembled'),b.cat('rt_enzyme','Source8.7ulRTEnzymeC'),b.cat('pap','Source3ul50UulPAP'),b.stock('atp','Source18ul10mMATP covalentsubstrate'),b.stock('capture_beads','Actualhistoricalv3.1gelbead/capturearmstock, wholeoligosunassembled',chemical_state={'cell_barcode':feature('sccomplete.cell','cell_barcode','unknown'),'umi':feature('sccomplete.umi','umi','unknown'),'lengths':'unassembled'}),b.cond('reducing','Source2ulReducingAgentB')))
    b.op('pap37','Source37C25min permitsRNApolyadenylation, allGEMstocksretained',inputs=(b.cond('pap_program','Source37C25min'),),product='SourceRNA/capture/protein-tagmixture withPAPtailingclasses; occupiedends/lengths/recoveryunassembled')
    b.op('rt53','Source53C45minRT/TSthen85C5min',inputs=(b.cond('rt_program','Source53C45min/85C5min, originalRTenzyme/PAP/TSO/beads/NTstocksretained; nofreshprimers'),),product='SourcecDNAofRNA plusADT/HTOreporterDNAs, free/reactionresidualclasses unassembled')
    b.clean('recover','SourceRecoveryAgentdemulsification/Dynabeads/silanecleanup followinghistoricalv3.1; eluteDNAstock')
    b.aliquot('udg_input','Source30uleluteselected','Actual30ulelutedDNA','Otherrecoveredelutedstock, amount/fateunassembled')
    b.op('udg','Source7ulUDGmix/37C30min85C5min',S.DIGESTION,(b.cat('udg_enzyme','Source3.5ulUDG-DNAglycosylaseNEBM0280S'),b.cond('udg_buffer','Source3.5ulUDGbuffer')),product='Source37ulDNA/cleavagepreparation; unreactedTSOremovalintended, cleavagepositions/completeproductsunassembled')
    # Same source100ulpreparedpreampmix contributes three actual whole oligos and molecularNT/enzyme components.
    stock=g.input('preamp_stock','Actualsourcepreampreagentstock:50ul2xAmpMix, final1.5uMforward/reverse/.005uMHTOadd,100ultotalwith37ulDNA',K.MOLECULAR)
    f=g.product('preamp_f','Samepreampforward',K.OLIGO,state=oligo(PRE_F));r=g.product('preamp_r','Samepreampreverse',K.OLIGO,state=oligo(PRE_R));h=g.product('preamp_h','SamepreampHTOadditive',K.OLIGO,state=oligo(HTO_ADD),chemical_state={'source_ordering':SOURCE_ORDERINGS['HTOadd'],'PS_last_two_bonds':True});e=g.product('preamp_enzyme','SameAmpMixenzyme',K.CATALYST);nt=g.product('preamp_nt','SameAmpMixnucleotidecomponent',K.MOLECULAR);buf=g.product('preamp_buffer','Samepreparedpreampconditions',K.CONDITION)
    g.react('preamp_components','Copresentcomponentsofonepreparedpreampstock, notphysicalfractions',S.GENERIC,((stock,R.TEMPLATE),),(f,r,h,e,nt,buf),changes_molecule=False)
    b.op('preamp','Source100ul12cyclepreamp98C15s63C20s72C1min;98C3mininitial/final72C1min',S.PCR,((f,R.STRUCTURAL),(r,R.STRUCTURAL),(h,R.STRUCTURAL),(e,R.CATALYST),(nt,R.STRUCTURAL),(buf,R.CONDITION)),geometry=PCRGeometry.unavailable(unknown_scene(),reason='WholeoccupiedRNA/ADT/HTOcDNAprimersites/molecularproductsunassembled; exactthreepreampoligosretained.'))
    sm=g.product('small_aliquot','Actual25ulsmallcDNAbranch',K.MOLECULAR);lg=g.product('large_aliquot','Actual25ullargeRNA/ADT/HTObranch',K.MOLECULAR);st1=g.product('stored1','Actual25ulstoredaliquot1',K.MOLECULAR);st2=g.product('stored2','Actual25ulstoredaliquot2',K.MOLECULAR)
    g.react('four_aliquots','Source100ulamplifiedpreparation dividedintofourequal25ulaliquots',S.SPLIT,((b.x,R.TEMPLATE),),(sm,lg,st1,st2));g.mark_terminal(Terminal(st1.id,Disposition.UNRESOLVED,'Source−20Clongtermstorage'));g.mark_terminal(Terminal(st2.id,Disposition.UNRESOLVED,'Source−20Clongtermstorage'))
    b.x=lg
    large=g.product('large_fraction','Source.6xpelletlargecDNA',K.MOLECULAR);sn=g.product('protein_supernatant','Actual.6xsupernatantcontinuesto2xselection',K.MOLECULAR)
    g.react('large_select','Source.6xSPRIphysicalpelletandsupernatant',S.SPLIT,((b.x,R.TEMPLATE),b.cond('large_beads','Source.6xSPRI')),(large,sn));b.x=sn
    b.clean('protein_select','Actual.6xsupernatant→2xSPRIprotein/HTOfraction; elutionvolumeunassembled')
    hto=g.product('hto_input','Actual5ulHTOtemplate',K.MOLECULAR);adt=g.product('adt_input','Actual5ulADTtemplate',K.MOLECULAR);rest=g.product('protein_remainder','Untransferredprotein/HTOfraction, amount/fateunassigned',K.MOLECULAR)
    g.react('protein_allocations','Source5ulHTO plusadditional5ulADT andactualstockremainder',S.SPLIT,((b.x,R.TEMPLATE),),(hto,adt,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferredprotein/HTOstock'))
    finished=[]
    for name,start,tail,cycles,program in [('small',sm,SHORT_R_TAIL,11,'65C30s72C1.5min final72C3min'),('adt',adt,SHORT_R_TAIL,15,'60C30s72C20s final72C5min'),('hto',hto,HTO_ADD,15,'64C30s72C20s final72C5min')]:
        b.x=start
        if name=='small':b.clean('small_select','Source1xSPRIretainedsmallcDNA, notlargefraction fragmentation')
        chem={'f':{'source_ordering':SOURCE_ORDERINGS['libF'],'PS_last_two_bonds':True},'r':{'source_ordering':SOURCE_ORDERINGS['HTOr'if name=='hto'else'shortR']}}
        if name=='hto':chem['r']['PS_last_two_bonds']=True
        b.pcr(name+'_pcr',f'Source1xKAPAHiFi/1uMeachwholeprimer;95C3min/{cycles}cycles98C20s/'+program,oligo(LIB_F),indexed_scene(tail),primer_chemicals=chem)
        finished.append((name,b.x))
    b.x=large
    large_stock=g.input('large_construct_stock','Actualhistoricalv3.1largeconstructionpreparedstocks; molecularprecursors/formulationsunassembled, knownsourcePCRoligosretained',K.MOLECULAR)
    lf=g.product('large_construct_f','SourceTableS2libraryforward fromsameconstructioninventory',K.OLIGO,state=oligo(LIB_F),chemical_state={'source_ordering':SOURCE_ORDERINGS['libF'],'PS_last_two_bonds':True})
    lr=g.product('large_construct_r','SourceTableS2fragmented-libraryindexedreverse fromsameinventory',K.OLIGO,state=indexed_scene(LONG_R_TAIL),chemical_state={'source_ordering':SOURCE_ORDERINGS['longR']})
    other=g.product('large_construct_other','Otheractualhistoricalconstructioninventory, enzyme/NT/adapterintermediatesunassembled',K.MOLECULAR)
    g.react('large_stock_components','CopresentknownsourcePCRoligos andotherhistoricalconstructioninventory, notextraadditions',S.GENERIC,((large_stock,R.TEMPLATE),),(lf,lr,other),changes_molecule=False)
    b.op('large_construct','Sourcehistoricalv3.1largecDNAfragmentation/endrepair/dAtail/adapterligation/indexing',inputs=((lf,R.STRUCTURAL),(lr,R.STRUCTURAL),(other,R.STRUCTURAL)),product='Sourcelargefractionindexedlibrarypreparation; wholeends/adapters/productsunassembled')
    finished.append(('large',b.x));selected=[]
    for name,prepared in finished:
        b.x=prepared
        if name in ('small','large'):
            b.aliquot(name+'_deplete_input','Source50ngselectedRNAcDNAlibraryfor twoSEQuoiaRibodepletionreactions','Actualselected50ngRNAcDNAlibrary','UntransferredRNAcDNAlibrarystock')
            b.op(name+'_deplete','SourceSEQuoiaRibodepletiontwo-reactionworkup, wholevendorprocedureunassembled',inputs=(b.stock(name+'_deplete_stock','ActualSEQuoiamolecularreactionstock; wholevendorrevision/formulationunassembled'),),product='SourcepostdepletionRNAcDNAlibrarypreparation; rRNA/otherclasses/yieldunassigned')
        b.aliquot(name+'_qc','SourceBioanalyzeranalyticalsize/concentrationQC amountunassigned',qc=True)
        selected.append(b.aliquot(name+'_pool_input','Sourceactualindexedlibrarypoolallocation, volumeunassigned','Selected'+name+'librarypoolaliquot','Untransferred'+name+'librarystock'))
    pool=g.product('library_pool','ActualsourceindexedlargeRNA/smallRNA/ADT/HTOsequencingpool',K.MOLECULAR);g.react('pool','Sourcephysicalindexedfourlibrarypool, exactamounts/indexchoicesunassigned',S.POOL,tuple((x,R.TEMPLATE)for x in selected),(pool,));b.x=pool
    b.finish('library','scComplete-seqselectedPBMCfourlibrarypool','PublishedNovaSeq6000PE150; FigureS1assignsR1cellbarcode/UMI andR2transcript, exactwholebeadarms/identifierlengths/runmask/primersitesunassembled. No genomicDNA or completeallRNAclassrecovery claim.')
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('LiteralUretained',g.materials['lna_tso'].chemical_state['literal_U_retained'],True)
    c('OnlyoneTSOaddition',sum(any(p.material=='lna_tso'for p in step.inputs)for step in g.steps.values()),1)
    c('PreamplifiedDNAfourphysicalaliquots',len(g.steps['four_aliquots'].outputs),4)
    c('Secondstoredstocknotsequenced',g.terminals['stored2'].disposition,Disposition.UNRESOLVED)
    c('Proteinselectionusesactual.6supernatant',g.steps['protein_select'].inputs[0].material,'protein_supernatant')
    c('HTOandADTeachactual5ul',len(g.steps['protein_allocations'].outputs),3)
    c('NoRNAdepletionimposedonADT','adt_deplete' in g.steps,False)
    c('Fouractualpoolcontributions',len(g.steps['pool'].inputs),4)
    c('ShortRNAdepletionafterindexPCR',g.steps['small_deplete_input'].inputs[0].material,'small_pcr_product')
    p=sources_root()/Path(__file__).resolve().parents[1].name/'nar_gkaf699_supp_gkaf699_supplemental_files_Supplementary_scCompleteseq.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1093/nar/gkaf699')):
        text=''.join(p.read_text().split())
        for name,seq in SOURCE_ORDERINGS.items():c('Wholemodifiedsourceordering '+name,seq in text,True)
    c('Partialnotpublishedready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
