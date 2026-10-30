"""Private partial sc-rDSeq preparation; inherited pre-IVT chemistry unassembled."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from chemdraw import Scene,Segment,feature
from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2,INDEX2_PRIMER_RC
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
from pcr import PCRGeometry
AUTHORITY='PublishedPMC13076218main70–76/98–142; SIactualFigS1p2/TableS1p20; nestedrDS_220.xlsxfirstsheetrows1–6. Closed3originals.'
ANCHOR='CGATGACGTAATACGACTCACTATAGGGATACCACCATGGCTCTTTCCCTACACGACGCTCTTC'
SP1='AGCAATCACTC'
SP2='TAGTAGACGAC'
RDS1='TGCGCA'
RUN_LENGTHS={'Read 1':50,'Read 2':70,'Index 1 (i7)':8,'Index 2 (i5)':8}
SOURCE_OLIGO_EXCEPTIONS={
    'round2_primer':'PublishedSIactualTableS1p20 prints whole NNNTAGTAGACGACNNN[8bp_Barcode2]AGCAATCACTC. Generic scraper mis-expands this named barcode macro; local source fixture checks the whole printed ordering and model retains exactly8 unknown cell-barcode positions. No source/chemistry waiver.'}
def scene(parts,**mods):
    sc=Scene();sc.strand('oligo',parts,**mods);return sc
def anchor_scene():return scene([Segment('source anchor',ANCHOR)],mod5='acrydite/iSpPC attachment/photocleavable spacer')
def umi(part):return Segment('UMI '+part,'NNN',placeholder=True,feature=feature('scrdseq.umi.'+part,'umi','random',group='scrdseq.umi',part=part,note='Complete combined identity/overlap/read spans unassembled'))
def bc1_scene():return scene([Segment('Spacer1 star',SP1),Segment('BC1','[8–10nt BC1]',placeholder=True,feature=feature('scrdseq.bc1','cell_barcode','combinatorial',group='scrdseq.cell',part='round1')),Segment('PE1 star',INDEX2_PRIMER_RC[:30])])
def bc2_scene():return scene([umi('round2_left'),Segment('Spacer2 star',SP2),umi('round2_right'),Segment('BC2','NNNNNNNN',placeholder=True,feature=feature('scrdseq.bc2','cell_barcode','combinatorial',group='scrdseq.cell',part='round2')),Segment('Spacer1 star',SP1)])
def capture_scene(which='rds',member=None):
    if which=='rds':
        if member is not None and(len(member)!=6 or set(member)-set('ACGT')):raise ValueError('rDS member must have six DNA bases')
        parts=[Segment('rDS star',member or '[6nt rDS member]',placeholder=member is None)]
    elif which=='polyt':parts=[Segment('captureN8','NNNNNNNN'),Segment('C5','CCCCC'),Segment('A12','A'*12)]
    else:raise ValueError('Unknown source capture branch')
    return scene(parts+[umi('round3'),Segment('Spacer2 star',SP2)])
def random_rt_scene():return scene([Segment('PE2 arm',TRUSEQ_READ2[-21:]),Segment('random RT hexamer','NNNNNN',note='Priming, not UMI')])
def index_scene(which):
    if which not in ('pe1','pe2'):raise ValueError('Select PE1/PE2')
    return scene([Segment('flowcell arm',P5 if which=='pe1'else P7),Segment('sample index','NNNNNNNN',placeholder=True,feature=feature('scrdseq.sample.'+which,'sample_index','unknown')),Segment('read arm',TRUSEQ_READ1 if which=='pe1'else TRUSEQ_READ2)])
def unknown_scene():return scene([Segment('whole products','[unassembled products]',placeholder=True)])
def seq(sc):return ''.join(s.top for st in sc.strands.values()for s in st.segs)

def protocol_graph():
    g=ProtocolGraph('sc-rDSeq three-round factory/selectedsinglecellRNA library',authority=AUTHORITY,representation='published_instance')
    cells=g.input('cells','Sourcecells300000perml/PBS10percentOptiPrep; exactidentity/transcriptsunassigned',K.BIOLOGICAL,primary=True)
    x=g.input('beads','Sourcepolymerizedacrydite-primerhydrogels; exactanchorcomponent only, wholehydrogel/inheritedpolymerizationopaque',K.MOLECULAR,state=anchor_scene(),chemical_state={'modifications':('5acrydite','iSpPC non-nucleotide spacer')})
    def inp(k,label,kind=K.CONDITION,role=R.CONDITION,**kw):return g.input(k,label,kind,**kw),role
    def op(k,label,kind=S.MAGIC,ports=(),product='Sourceheterogeneousmixture; fullsubstrates/products/off-products/occupancyunassembled',**kw):
        nonlocal x
        y=g.product(k+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Positive source operation; fullmolecularproducts/endclasses/recoveryunassembled, no sourceabsenceclaim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(ports),(y,),**kw);x=y
    def clean(k,label):
        nonlocal x
        y=g.product(k+'_kept','Actual retainedfraction; imperfectrecovery/carryoverunassigned',K.MOLECULAR);w=g.product(k+'_waste','Actual removedworkupfractions',K.WASTE)
        g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),inp(k+'_conditions',label)),(y,w));g.discard(w,'Physicalworkup');x=y
    def allocate(k,label,qc=False):
        nonlocal x
        y=g.product(k+'_selected',label,K.MOLECULAR);r=g.product(k+'_remainder','Actualretained/untransferredstock; fateunassigned',K.MOLECULAR)
        g.react(k,label,S.SPLIT,((x,R.TEMPLATE),),(y,r))
        g.mark_terminal(Terminal(y.id if qc else r.id,Disposition.QC if qc else Disposition.UNRESOLVED,label if qc else 'Untransferredstocknotautomaticallydiscarded'));x=r if qc else y
    def prepared(k,label,spec):
        st=g.input(k,label,K.MOLECULAR);out=[];ports=[]
        for suffix,kind,role,text in spec:
            m=g.product(k+'_'+suffix,text,kind);out.append(m);ports.append((m,role))
        g.react(k+'_components','Representcopresentcomponentsofsourcepreparedmix, noadditionalfreshsupplies',S.GENERIC,((st,R.TEMPLATE),),tuple(out),changes_molecule=False);return tuple(ports)
    for i in (1,2,3):
        k='round'+str(i)
        if i<3:
            ports=[inp(k+'_primer','Source9ul50uMBCprimer; actualwellbarcodeunassigned',K.OLIGO,R.STRUCTURAL,state=bc1_scene()if i==1 else bc2_scene())]
        else:
            ports=[inp('rds_pool','Source220memberrDSprimertemplatepool; fullmembersunassembled',K.OLIGO,R.STRUCTURAL,state=capture_scene()),inp('polyt_template','SourcepolyT5C8Nstar template',K.OLIGO,R.STRUCTURAL,state=capture_scene('polyt'))]
        op(k+'_load','Source6ulbeadmixture/2.5xbuffer/.85mMdNTP;85C2min60C20min;384wellBCrounds or commonthird4:1pool',ports=tuple(ports)+(inp(k+'_loading_nts','Nucleotidecomponentcopresentin6ulbeadmixture/.85mM, notanadditionalfreshdose',K.MOLECULAR,R.STRUCTURAL),inp(k+'_conditions','Sourcepriming/beadmixtureconditions; exactrouting/stockpackagingunassembled'),))
        ext=prepared(k+'_mix','Actual5ulextensionmix1.8Upolymerase/.3mMdNTP/1xbuffer',[
            ('enzyme',K.CATALYST,R.CATALYST,'SourceNEBM0537component'),('nts',K.MOLECULAR,R.STRUCTURAL,'Nucleotidecomponentofsamemix'),('buffer',K.CONDITION,R.CONDITION,'Buffercopresentinsame5ulmix')])
        op(k+'_extension','Source60C60minextension; fullbeadproductsunassembled',ports=ext)
        op(k+'_stop','Source20ulstop100mMKCl/10mMTris8/50mMEDTA/.1percentTween/RT10min',S.GENERIC,(inp(k+'_stop_buffer','Actualsource20ulstopbuffer'),))
        allocate(k+'_qc','Source.5ulbeadsin9.5ulwaterbeforewashing; laterUV10min/2percentagaroseQC',qc=True)
        clean(k+'_wash','Sourcebeadcollection/single-strandwash/split-poolrouting; inheritedwashdetailsunassembled')
    op('exoi','Sourcepostprimingwash twiceT10E.1T.1/ExoI.25Uul37C60mindark',ports=(inp('exoi_enzyme','SourceExoI',K.CATALYST,R.CATALYST),inp('exoi_conditions','Source1xbuffer/temperature; failedprimerremovalnotperfectspecificityclaim')))
    op('stop_exoi','Source10mlstop/ice10min',S.GENERIC,(inp('exoi_stop','Actualsourcestopbuffer'),))
    clean('factory_final','Sourcewash/two70umstrainers; completeandfailedprimerclassesunquantified')
    allocate('bead_feed','Actual100–150ulbeadstockinletallocation;1500rcf2min')
    clean('bead_pellet','Sourceaspirateresidualbufferfrompackedbeads')
    feed=g.input('rt_feed_stock','Actual150ulfirstRTlysisstock:30buffer/15enzyme/10inhibitor/15Tris/10DTT/9IGEPAL/6dNTP/55water',K.MOLECULAR)
    feed_used=g.product('rt_feed_used','ActualRTlysisfeedallocation, volumeunassigned',K.MOLECULAR);feed_rest=g.product('rt_feed_rest','Untransferred150ulRTlysisstockremainder; amountunassigned',K.MOLECULAR)
    g.react('rt_feed_split','Allocateactualinletstock, not150uldropletvolume',S.SPLIT,((feed,R.TEMPLATE),),(feed_used,feed_rest));g.mark_terminal(Terminal(feed_rest.id,Disposition.UNRESOLVED,'Untransferredfeedstock'))
    fe=g.product('first_rt_enzyme','EnzymecomponentofallocatedRTlysisstock',K.CATALYST);fn=g.product('first_rt_nts','NucleotidecomponentofallocatedRTlysisstock',K.MOLECULAR);fb=g.product('first_rt_buffer','Buffer/Tris/DTT/inhibitor/IGEPAL/watercomponents',K.CONDITION)
    g.react('rt_feed_components','Copresentcomponentsoftheactualpreparedfeedallocation, notnewadditions',S.GENERIC,((feed_used,R.TEMPLATE),),(fe,fn,fb),changes_molecule=False)
    op('droplets','Sourcecells/beads/RTlysis/oil→1–3nldroplets;selectedsc-rDwell50uloil+20ulmineralaftercollection',ports=((cells,R.TEMPLATE),(fe,R.CATALYST),(fn,R.STRUCTURAL),(fb,R.CONDITION),inp('oil','Source3mlHFE7500/2percentsurfactantstock; exactallocation/remainderunassigned',K.MOLECULAR,R.STRUCTURAL)))
    op('uv','Source365nm6.5Jcm-2/ice10minprimerrelease',ports=(inp('uv_conditions','Sourcephotocleavageconditions; remainingbead/linkerfatesunassembled'),))
    op('first_rt','Source20rampedRTcycles8C12s15C45s20C45s30C30s50C15s; enzyme/NT/releasedprimerscopresent',ports=(inp('first_rt_conditions','Sourcethermalprogramme, nofreshprimers'),))
    firstwell=g.product('first_well','One actual selected collected well; identity/volume unassigned',K.MOLECULAR)
    otherwells=g.product('other_wells','Other actual1–2 selected collected wells; exact count/identities/volumes unassigned',K.MOLECULAR)
    unselected=g.product('unselected_wells','Unallocated well-array material, if any; existence/amount/fate unassigned',K.MOLECULAR)
    g.react('well_allocation','Allocate actual2–3 collected-well contributions from firstRT well array; no fixed3-well or volume assignment',S.SPLIT,((x,R.TEMPLATE),),(firstwell,otherwells,unselected))
    g.mark_terminal(Terminal(unselected.id,Disposition.UNRESOLVED,'Any unselected array material is not implicitly pooled; no positive remainder quantity asserted'))
    pooled=g.product('well_pool_product','Actual pooled2–3well sample/1000–3000cells; identities/individual volumes unassigned',K.MOLECULAR)
    g.react('well_pool','Physically combine one selected well plus actual remaining1–2 selected wells',S.POOL,((firstwell,R.TEMPLATE),(otherwells,R.TEMPLATE)),(pooled,));x=pooled
    op('through_ivt','SourceinDropsref18routeafterfirstRTthroughIVT, full inheritedrecipeUNASSEMBLED',ports=(inp('inherited_supplies','Unassembledinheritedreagentinventory/packaging, no otherprotocolrecipeimported',K.MOLECULAR,R.STRUCTURAL),),product='SourceaRNA/amplificationmixture; completeDNAsecondstrand/dropletbreakup/IVTproductsunassembled')
    clean('arna_clean','Source1.3xAMPureafterIVT/elute20ul10mMTris')
    saved=g.product('arna_backup','Actual9ulaRNAbackup-80C',K.MOLECULAR);used=g.product('arna_used','Actual9ulaRNAproduction',K.MOLECULAR);rest=g.product('arna_rest','Nominal2ulaRNAunallocated',K.MOLECULAR)
    g.react('arna_split','Actual20uleluate→9backup/9production/nominal2unassigned',S.SPLIT,((x,R.TEMPLATE),),(saved,used,rest));g.mark_terminal(Terminal(saved.id,Disposition.STORED,'SourceaRNAbackup'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'UnallocatedaRNA'));x=used
    op('fragment_arna','Source9RNA+1AM8740/70C1minice;40ulstop(5stop+35Tris)',ports=(inp('fragment_reagent','Actual1ulfragmentationreagent',K.MOLECULAR,R.STRUCTURAL),inp('fragment_conditions','Sourcefragmentation/stopconditions')))
    clean('fragment_clean','Source1.3xAMPure/elute10ulTris')
    op('second_rt_input','Source10uleluate→specified12ulRNAinput, unresolvedstockorigin/dilution/allocation',reason='Sourceallocationtraceunreconciled; no2ultopup, invented physicalsplit/remainder ornormalization.',product='Source-specifiedsecondRT RNAinput; actual12ulstockorigin/reconciliationwith10uleluateunassembled')
    op('second_rt','Source20ulsecondRT:12RNA/2primer/1dNTP/2buffer/1DTT/1inhibitor/1enzyme;25C10min50C60min80C10min',ports=(inp('pe2_random','Source2ul100uMPE2-N6-v2',K.OLIGO,R.STRUCTURAL,state=random_rt_scene()),inp('rt2_enzyme','Source1ulStellarScript200Uul',K.CATALYST,R.CATALYST),inp('rt2_nts','Actual1ul10mMdNTP',K.MOLECULAR,R.STRUCTURAL),inp('rt2_conditions','Sourcebuffer/DTT/inhibitor; distinctfrom150ulfirstRTstock')))
    op('dilute','Source+20ulwaterto40ul',S.GENERIC,(inp('dilution_water','Actual20ulwater'),))
    keep=g.product('cdna_production','Actual20ulproductionhalf',K.MOLECULAR);save=g.product('cdna_backup','Actual20ulbackup-80C',K.MOLECULAR)
    g.react('cdna_split','Actual40ul→20backup/20production',S.SPLIT,((x,R.TEMPLATE),),(keep,save));g.mark_terminal(Terminal(save.id,Disposition.STORED,'SourcecDNAbackup'));x=keep
    clean('cdna_clean','Source1.2xAMPure/elute11.5ulTris')
    primers=g.input('index_primer_mix','ActualONE1ul25uM PE1/PE2 index-primer mix',K.MOLECULAR)
    f=g.product('pe1_index','PE1 oligo copresent in same1ul primer mix',K.OLIGO,state=index_scene('pe1'));r=g.product('pe2_index','PE2 oligo copresent in same1ul primer mix',K.OLIGO,state=index_scene('pe2'))
    g.react('index_primer_components','Expose two copresent oligos of ONE source1ul mix; no fresh second addition',S.GENERIC,((primers,R.TEMPLATE),),(f,r),changes_molecule=False)
    pm=prepared('kapa_mix','Actual12.5ul2xKAPAHiFiHotStartReadyMix',[('enzyme',K.CATALYST,R.CATALYST,'SameReadyMixpolymerase'),('nts',K.MOLECULAR,R.STRUCTURAL,'SameReadyMixnucleotides'),('buffer',K.CONDITION,R.CONDITION,'SameReadyMixconditions')])
    op('pcr','Source11.5cDNA+12.5ReadyMix+1primer-mix;2lowanneal+8highannealcycles',S.PCR,((f,R.STRUCTURAL),(r,R.STRUCTURAL))+pm,geometry=PCRGeometry.unavailable(unknown_scene(),reason='Sourceprintedprimerarms/indexwidthsknown; actualindexvalues/occupiedtemplates/fullproducts/readlandingandfeaturecyclesunassembled.'))
    clean('library_clean','Source.7xAMPure/elute20ulwater')
    allocate('library_qc','SourceTapeStationD1000QCaliquot, volumeunassigned',qc=True)
    allocate('loading','Actualloadingallocation; untransferredlibrarystockunassigned')
    lib=g.product('library','Selectedsinglecellsc-rDSeqNGSlibrary, wholemoleculesunassembled',K.LIBRARY)
    g.react('load','SourceNovaSeq6000S1kitlabel100cycles; actualR1/R2/i1/i2=50/70/8/8',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(read_lengths=RUN_LENGTHS,limitation='Sourceactualsinglecell50/70/8/8; kitlabel100cycles retainedwithoutkiterrorclaim. No55/65bulkmerge. Exactprimerlanding/featurecycles/barcodesunassembled; rRNAalignmentfilter/dedup/trimmingcomputational.'))
    g.add_source_item(SourceItem('volume_trace','MainfragmentaRNAeluate10ul; nextsecondRTcalls12ulRNAin20ulreaction.','Retainbothsourcevalues; stockorigin/dilution/allocationunassembled, no2ultopuporpublicationerrorclaim.150ulfirstRTinletstockisdifferentstage.'))
    g.add_source_item(SourceItem('inheritance','Main70–76/98–118 defines inDrops-based three-round rDS/polyT5C8N capture and delegates post-firstRT through IVT to ref18, identified at main416–420 as Zilionis2017,10.1038/nprot.2016.154.','Noautomaticfamilyrecipeinheritance orlaterBiomethodsrevisionmerge; completeinheritedoperationsunassembled.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();p=sources_root()/HERE.parent.name/'nar_gkag312_supp_gkag312_supplemental_files_supplementary_data.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1093/nar/gkag312')):
        t=''.join(p.read_text().split()).upper()
        for name,value in [('modifiedanchor','/5ACryd/iSpPC/'+ANCHOR),('BC2template','NNN'+SP2+'NNN[8bp_Barcode2]'+SP1),('randomRT',seq(random_rt_scene())),('PE1index',P5+'[8bp_Index1]'+TRUSEQ_READ1),('PE2index',P7+'[8bp_Index2]'+TRUSEQ_READ2)]:c('Wholepublishedordering '+name,value.upper()in t,True)
    z=sources_root()/HERE.parent.name/'nar_gkag312_supp_gkag312_supplemental_files.zip'
    if have(c,Source(z,'https://doi.org/10.1093/nar/gkag312')):
        import zipfile,io,xml.etree.ElementTree as E
        with zipfile.ZipFile(z)as outer:data=outer.read('Supplementary File 2.sc-rDSeq Barcode primers.zip')
        with zipfile.ZipFile(io.BytesIO(data))as inner:data=inner.read('rDS_220.xlsx')
        with zipfile.ZipFile(io.BytesIO(data))as book:
            ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
            strings=[''.join(n.itertext())for n in E.fromstring(book.read('xl/sharedStrings.xml')).findall('s:si',ns)]
            row=E.fromstring(book.read('xl/worksheets/sheet1.xml')).findall('s:sheetData/s:row',ns)[1]
            cells=[]
            for cell in row:
                v=cell.find('s:v',ns);cells.append(strings[int(v.text)]if cell.attrib.get('t')=='s'else v.text)
            c('PublishednestedworkbookrDS1wholeordering',''.join(cells[1].split()),RDS1+'NNN')
    c('AcryditePCaddsnoDNAbase',seq(anchor_scene()),ANCHOR)
    c('BC1lengthunassigned',bc1_scene().strands['oligo'].segs[1].length_bp,None)
    c('Cellbarcodestructured',bc2_scene().strands['oligo'].segs[3].feature.role,'cell_barcode')
    c('UMIstructured',bc2_scene().strands['oligo'].segs[0].feature.role,'umi')
    c('SecondRTN6notUMI',random_rt_scene().strands['oligo'].segs[1].feature,None)
    try:capture_scene('rds','AAA');bad=False
    except ValueError:bad=True
    c('RejectmalformedrDSmember',bad,True)
    c('FactoryQCnotproduction',g.terminals['round1_qc_selected'].disposition,Disposition.QC)
    c('aRNAbackupstored',g.terminals['arna_backup'].disposition,Disposition.STORED)
    c('aRNArestnotwaste',g.terminals['arna_rest'].disposition,Disposition.UNRESOLVED)
    c('FirstRTafterUV',g.steps['first_rt'].inputs[0].material,'uv_product')
    c('SecondRTsourceallocationtrace',g.steps['second_rt'].inputs[0].material,'second_rt_input_product')
    c('cDNAbackupstored',g.terminals['cdna_backup'].disposition,Disposition.STORED)
    c('Cleanupusesproductionhalf',g.steps['cdna_clean'].inputs[0].material,'cdna_production')
    c('Oneactualindexmix',g.steps['index_primer_components'].inputs[0].material,'index_primer_mix')
    c('PCRusescopresentindexcomponents',set(p.material for p in g.steps['pcr'].inputs if p.material in ('pe1_index','pe2_index')),{'pe1_index','pe2_index'})
    c('Physicalwellpool',g.steps['well_pool'].kind,S.POOL)
    c('Pooledwellcontributions',set(p.material for p in g.steps['well_pool'].inputs),{'first_well','other_wells'})
    c('WholePCRproductsunassembled',g.materials['pcr_product'].state,None)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
