"""Published selected preparation; private partial molecular model, not a complete audit."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from circular import CircularizationGeometry

TITLE='MethylSaferSeqS selected paired LA/LM WGS preparation'
AUTHORITY='PublishedPMC10104560 main244/250–274/398–410;publishedSIwholeTableS1 andoriginalp4inspected;2originals. SaferSeqSbasekitinternalsnotinherited.'
INITIAL='Source cfDNA molecular preparation; extraction/sourcecohort identityunassigned; originalcytosinemethylationretained'
READOUT='SourceNovaSeq6000PE100 WGS. LMandLA arephysicallydifferent libraries fromsharedoriginalUIDs; nativeLMbisulfiteconvertschemicalCtoU thenPCRreportsT, whilecopiedLAhasunmodifiedCs. No quantitativeperfectduplexyield/strandpartition/conversionguarantee.'
DEFERRED='HistoricalSaferSeqSpreligation/fullligation/primerlandinggeometry andvendororiginals;sourcePanel48/hemi-nestedvalidationbranchesnotassembled. SelectedWGSbothbranchesretained,notmistakenforoneconcatenatedmolecule.'

def inventory_state(label):return None
def inventory_chemistry(label):return None

from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2,NEBNEXT_ARM_READ2
from chemdraw import feature
R2_ARM=TRUSEQ_READ2[2:]
def ordered_sequence(name):
    return {'Source3primeadapter':'ATAAAACGACGGC'+'N'*14+'A'+NEBNEXT_ARM_READ2,'Source5primeadapter':TRUSEQ_READ1,'Sourcebiotin-dUcopyprimer':'T'*10+'U'+R2_ARM,'SourceLAprimer1':TRUSEQ_READ1,'SourceLAprimer2':R2_ARM,'SourceLMprimer1':TRUSEQ_READ1[:22],'SourceLMprimer2':R2_ARM[:21],'SourceWGSprimer1':P5+'N'*8+TRUSEQ_READ1,'SourceWGSprimer2':P7+'N'*8+TRUSEQ_READ2}.get(name)
def inventory_chemistry(name):
    seq=ordered_sequence(name)
    if not seq:return None
    d={'polymer':'DNA','sequence':seq}
    if name in {'Source3primeadapter','Source5primeadapter'}:d['five_methyl_c_positions']=tuple(i for i,b in enumerate(seq) if b=='C')
    if name=='Source3primeadapter':d.update(five_prime='P',three_prime='SpC3',phosphorothioate_after=(len(seq)-3,len(seq)-2),uid_start=13,uid_length=14)
    if name=='Sourcebiotin-dUcopyprimer':d.update(five_prime='dual biotin /52-Bio/',spacer='iSp18',deoxy_u_positions=(10,))
    if name in {'SourceLAprimer2','SourceLMprimer2'}:d['five_prime']='dual biotin /52-Bio/'
    if name.startswith('SourceWGS'):d['phosphorothioate_after']=(len(seq)-3,len(seq)-2)
    return d
def inventory_state(name):
    seq=ordered_sequence(name)
    if not seq:return None
    sc=Scene();d=inventory_chemistry(name);parts=[]
    if name=='Source3primeadapter':parts=[Segment('sourceprefix',seq[:13]),Segment('source14UID',seq[13:27],placeholder=True,feature=feature('methsafer_uid','umi','random')),Segment('source3adapterarm',seq[27:])]
    elif name.startswith('SourceWGS'):
        n=len(P5) if name.endswith('1')else len(P7);parts=[Segment('sourceflowcellarm',seq[:n]),Segment('sourceWGSwellindex',seq[n:n+8],placeholder=True,feature=feature('methsafer_wgs_index','sample_index','unknown')),Segment('sourceprimerarm',seq[n+8:])]
    else:parts=[Segment(name,seq,note='Source ordering inventory; modifications retained in structured chemical_state, fullreactiongeometryunassembled')]
    end=d.get('five_prime','');end+=('/iSp18/'if 'spacer'in d else'');sc.strand(name,parts,mod5=end,mod3=d.get('three_prime',''));return sc

SOURCE_ORDERINGS={
'Source3primeadapter':{'source':'pnas.2220704120_supp_pnas.2220704120.sapp.pdf p4; .txt124–127','raw':'/5Phos/ATAAAA/iMe-dC/GA/iMe-dC/GG/iMe-dC/NNNNNNNNNNNNNNAGAT/iMe-dC/GGAAGAG/iMe-dC/A/iMe-dC/A/iMe-dC/GT/iMe-dC/TGAA/iMe-dC/T/iMe-dC//iMe-dC/AG*T*/iMe-dC//3SpC3/'},
'SourceWGSprimer1':{'source':'pnas.2220704120_supp_pnas.2220704120.sapp.pdf p4; .txt142–144','raw':'AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGAT*C*T'},
'SourceWGSprimer2':{'source':'pnas.2220704120_supp_pnas.2220704120.sapp.pdf p4; .txt145–147','raw':'CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT*C*T'}}
SOURCE_OLIGO_EXCEPTIONS={'base_ligate_oligo_0':'PublishedSIwholeTableS1originalp4 visually verified; mixed/iMe-dC/and/iMedC/andwrappedPS/blocker notation fragmentstextmatching. Local normalized-original wholeordering check retains all modifications separately.','lm_wgs_pcr_oligo_0':'PublishedSIoriginalTableS1p4 WGS1 fullrowverified; extraction interleaves purpose-column prose before wrapped C*T tail. Localwholeprimer fixture andstructuredPS/indexdata checked.','lm_wgs_pcr_oligo_1':'PublishedSIoriginalTableS1p4 WGS2fullrowverified; purpose-column prose interrupts wrapped terminal T. Localwholeprimer fixture/PS/indexdata checked.','la_wgs_pcr_oligo_0':'Same sourceTableS1p4 WGS1ordering asLM; separatephysicalPCR addition, locallychecked.','la_wgs_pcr_oligo_1':'Same sourceTableS1p4 WGS2ordering asLM; separatephysicalPCR addition, locallychecked.'}
def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('geometry','[source geometry unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('input',INITIAL,K.MOLECULAR,primary=True)
    def op(key,label,kind,reagent,product,oligos=()):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.CIRCULARIZATION,S.LIGATION,S.PCR,S.DEPROTEINIZATION,S.TAGMENTATION,S.TRANSCRIPTION,S.REVERSE_TRANSCRIPTION,S.CHEMICAL_CONVERSION}
        c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION)
        inputs=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]
        for i,a in enumerate(oligos):
            o=g.input(key+'_oligo_'+str(i),a,K.OLIGO,state=inventory_state(a),chemical_state=inventory_chemistry(a));inputs.append((o,R.STRUCTURAL))
        kw={'changes_molecule':False} if kind==S.GENERIC else {}
        if kind==S.MAGIC:kw['reason']='Actual published component; full source-specific molecular transition unassembled.'
        if kind==S.CIRCULARIZATION:kw['geometry']=CircularizationGeometry.unavailable('Source intramolecular closure known; full substrate and closure geometry unassembled, not a claim of source absence.')
        if kind==S.LIGATION:kw['geometry']=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Source ligation operation known; full ends and sealed junction unassembled.')
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown_scene(),reason='Source primer inventories/whole amplicon geometry not assembled in this selected model.')
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,inputs,(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Physical source recovered preparation; carryover and recovery unquantified',K.MOLECULAR);w=g.product(key+'_waste','Actual source physically separated fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Physical workup, not a perfect molecular specificity assertion.');x=y
    def aliquot(key,label):
        nonlocal x
        y=g.product(key+'_used',label,K.MOLECULAR);r=g.product(key+'_rest','Physicalsource remainder; mass/fateunassigned',K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Not implicitly discarded or counted as sequenced.'));x=y
    def qc(key,label):
        nonlocal x
        y=g.product(key+'_retained','Production preparation retained after independent analytical allocation',K.MOLECULAR);q=g.product(key+'_assay',label,K.QC);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,q));g.mark_terminal(Terminal(q.id,Disposition.QC,'Analytical aliquot; not returned to production.'));x=y
    op('base_preligation','SourceSaferSeqSpreligationcomponent; no endrepair/dAtail reaction assigned',S.MAGIC,'SourcecitedSaferSeqSpreparation; exacthistoricalinternalsunreviewed','SourcecfDNApreligationmixture; actualpreparativeoperationsandendsunassembled')
    op('base_ligate','SourceYadapterligationwithmethylatedcustomadaptors',S.LIGATION,'SourcebaseSaferSeqSligationcomponentunreviewed','SourceadapterligatedoriginalswithmethylatedadapterCsoutsideUID;partial/unligatedclassesremain',('Source3primeadapter','Source5primeadapter'))
    op('copy','Source3denaturation/copyingrounds;singleprimer5uM,notexponentialPCR',S.EXTENSION,'KAPAHiFiHotStartUracilPlus7959052001','Sourceoriginal/copiedhybridsplusadditionalunpairedcopiedstrands/freeprimer/reagents; nativeDNAmarksnotcopied',('Sourcebiotin-dUcopyprimer',))
    beads=g.input('capture_beads','SourceStreptavidinMyOneDynabeadsT1',K.SELECTION);bound=g.product('bound','Sourcebead-associateddualbiotincopiesincludingoriginal/copiedhybrids',K.MOLECULAR);unbound=g.product('capture_unbound','Physicalunretainedsourcefraction; exactcontent/fateunassigned',K.MOLECULAR);g.react('capture','Sourcecopiedstrandcapture',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,unbound));g.mark_terminal(Terminal(unbound.id,Disposition.UNRESOLVED,'No perfectcaptureorunboundpurityassumed.'))
    h=g.input('heat_water','Source22.5ulwater75C2mintwice; samebeadpreparationsequentially',K.CONDITION);original1=g.product('original1','Sourcefirst22.5ulsupernatantoriginal-strandenrichedfraction;copycarryoverunquantified',K.MOLECULAR);remain1=g.product('remain1','Sourcebead-retainedcopiedstrandsplusremainingoriginals',K.MOLECULAR);g.react('heat1','Sourcefirstheat/magneticseparation',S.STRAND_SEPARATION,((bound,R.TEMPLATE),(h,R.CONDITION)),(original1,remain1));original2=g.product('original2','Sourcesecond22.5ulsupernatant',K.MOLECULAR);remain2=g.product('remain2','Sourcebead-retainedcopiesaftertwoheatreleases',K.MOLECULAR);g.react('heat2','Sourcesecondheat/magneticseparation',S.STRAND_SEPARATION,((remain1,R.TEMPLATE),(h,R.CONDITION)),(original2,remain2));native=g.product('native_pool','Sourcecombined45ulsupernatantforLM',K.MOLECULAR);g.react('pool_originals','Actualtwoeluatescombined',S.POOL,((original1,R.TEMPLATE),(original2,R.TEMPLATE)),(native,))
    x=native
    op('bisulfite','SourceD5004:afterdenaturationadd100ulCTreagent;dark50C3h/4C10min',S.CHEMICAL_CONVERSION,'SourceEZDNAMethylationKitCTreagent','Sourcebisulfite-treatednativeDNAwithU/conversionintermediatesandprotectedmethylatedCs; exactconversionfractionunassigned')
    clean('bind_wash','SourceSiliconAbindingplate/600ulMbinding/2800g5min/500ulMwash;physicalwaste')
    op('desulfonate','Source200ulMDesulphonationRT15min',S.CHEMICAL_CONVERSION,'Sourcekitdesulphonationbuffer','SourcechemicalU-containingnative-derivedDNA andunprocessedclasses')
    clean('lm_wash','Source2800g5min/600ulthen400ulMwash;two21ulelutionscombined~42ul')
    op('lm_pcr','Source100ulUracilPlus/primers5uM;98C45s,8x98C15s60C30s72C30s',S.PCR,'KAPAHiFiHotStartUracilPlus7959052001','SourceLMlibraryPCRcopies;UreportedT,originalDNAmethylmarksnotchemicallycopied',('SourceLMprimer1','SourceLMprimer2'))
    clean('lm_clean','Source180ulSPRISelect/elute100ulEB')
    aliquot('lm_wgs_input','Source100ngLMforWGS;remainingLMretainedoutsidechosenreaction')
    op('lm_wgs_pcr','Source50ulUltraQ5/primers2uM;98C30s,7x98C10s65C75s',S.PCR,'SourceUltraQ5','SourceLMWGSindexedpreparation',('SourceWGSprimer1','SourceWGSprimer2'))
    lm=x
    x=remain2
    op('user','Source3ulUSER/37ulEB;40ul37C60min',S.DIGESTION,'USERM5505S','SourcecleavedcopiedDNAreleasedfromdualbiotin-dUprimerplusbead-boundprimerremnants/uncutclasses')
    clean('la_release','Sourceactualmagneticseparation;40ulreleasedsupernatantretained')
    op('la_pcr','Source100ulHiFiHotStartK2601/primers5uM;98C45s,8x98C15s60C30s72C30s',S.PCR,'KAPAHiFiHotStartK2601 (notLMUracilPlus)','SourceLAgeneticlibraryPCRcopies andoriginalreactioncompanions',('SourceLAprimer1','SourceLAprimer2'))
    clean('la_clean','Source180ulSPRISelect/elute100ulEB')
    aliquot('la_wgs_input','Source100ngLAforWGS;remainingLAretainedoutsidechosenreaction')
    op('la_wgs_pcr','Source50ulUltraQ5/primers2uM;98C30s,7x98C10s65C75s',S.PCR,'SourceUltraQ5','SourceLAWGSindexedpreparation',('SourceWGSprimer1','SourceWGSprimer2'))
    la=x
    for key,material in [('lm',lm),('la',la)]:
        lib=g.product(key+'_library','Source'+key.upper()+'WGSsequencingpreparation',K.LIBRARY);g.react(key+'_interface','Physical'+key+'libraryidentity',S.GENERIC,((material,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation=READOUT))

    g.add_source_item(SourceItem('pending',DEFERRED,'Selected original sections read; linked-source/whole molecular geometry review incomplete. No accepted absence claim.'))
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Copying uses singleprimerextension',g.steps['copy'].kind,S.EXTENSION)
    c('Secondheatretainsfirstbeadfraction',g.steps['heat2'].inputs[0].material,'remain1')
    c('LMreceivestwooriginaleluates',g.steps['bisulfite'].inputs[0].material,'native_pool')
    c('LAcopyreleaseusesremainingbeads',g.steps['user'].inputs[0].material,'remain2')
    c('LMhasseparateuracilPCR',g.steps['lm_pcr'].inputs[0].material,'lm_wash_kept')
    c('LAphysicallyreleasedbeforePCR',g.steps['la_pcr'].inputs[0].material,'la_release_kept')
    c('Twoactualsequencinglibraries',set(g.sequencing),{'lm_library','la_library'})
    c('SourceWGS1wholeordering',ordered_sequence('SourceWGSprimer1'),'AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT')
    c('SourceWGS2wholeordering',ordered_sequence('SourceWGSprimer2'),'CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')
    for name,entry in SOURCE_ORDERINGS.items():
        raw=entry['raw'].replace('/5Phos/','').replace('/3SpC3/','').replace('/iMe-dC/','C').replace('*','');c('Wholemodifiedorderingfixture '+name,ordered_sequence(name),raw)
    c('Sourcecopyprimer retains dualbiotin' ,inventory_chemistry('Sourcebiotin-dUcopyprimer')['five_prime'],'dual biotin /52-Bio/')
    c('SourceUIDdoesnotbecomeWGSindex',inventory_chemistry('Source3primeadapter')['uid_length'],14)
    from checks import Source,have
    from storage import sources_root
    import re
    src=Source(sources_root()/'methylsaferseqs__10.1073+pnas.2220704120'/'pnas.2220704120_supp_pnas.2220704120.sapp.pdf.txt','PublishedSITableS1originalp4')
    if have(c,src,label='Published3adapter'):
        raw=src.path.read_text().split('/5Phos/',1)[1].split('/3SpC3/',1)[0];raw=re.sub(r'\s+','',raw);raw=raw.replace('/iMe-dC/','C').replace('/iMedC/','C').replace('*','');c('Whole3adapter normalizedsource ordering',ordered_sequence('Source3primeadapter'),raw)
    c('Private partial model cannot be publicly promoted',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()
