from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
import nextera
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

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s
AUTHORITY='Published2023 Chem-map PMC10497411 Methods131–159/actualFig1PDFp2; explicitly delegated ref24 publishedPMC8654944 Methods87 +actualSI1p6–7 lines143–151 for loadedpA-Tn5 component ONLY. Threeoriginals; no G4 cellular recipe inherited.'
SOURCE_ME={'rev':'CTGTCTCTTATACACATCT','A':'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG','B':'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'}
# The three complete orders are verified in delegated ref24 SI1 p6–7;
# reuse only these matching canonical strings, not a Nextera reaction recipe.
ME={'rev':nextera.ME_RC,'A':nextera.ADAPTOR_S5,'B':nextera.ADAPTOR_S7}

def adapter(name):
    sc=Scene();sc.strand(name,[Segment(name,ME[name])],mod5='phosphate'if name=='rev'else'');return sc

def duplex(name):
    sc=Scene();s=ME[name];core=revcomp(ME['rev']);sc.strand(name,[Segment(name+' tail',s[:-len(core)]),Segment('ME core',core)]);sc.anneal('rev',[Segment('ME rev',ME['rev'])],to=name,pair=('ME rev','ME core'),mod5='phosphate');return sc

def allocate(g,key,stock,n,label,remaining):
    xs=[g.product(key+'_'+str(i),label+' '+str(i+1),K.MOLECULAR,state=stock.state,chemical_state=stock.chemical_state)for i in range(n)];r=g.product(key+'_remainder',remaining,K.MOLECULAR,state=stock.state,chemical_state=stock.chemical_state);g.react(key,'Actual selected aliquots; no duplicated whole stock',S.SPLIT,((stock,R.TEMPLATE),),tuple(xs)+(r,));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remaining));return xs

def components(g,key,label,classes):
    st=g.input(key+'_stock',label,K.MOLECULAR);xs=[g.product(key+'_'+name,caption,kind)for name,caption,kind,role in classes];g.react(key+'_parts','Copresent components of SAME actual physical stock; no extra fresh doses',S.GENERIC,((st,R.TEMPLATE),),tuple(xs),changes_molecule=False);return tuple((x,row[3])for x,row in zip(xs,classes))

def qc(b,key):
    g=b.g;q=g.product(key+'_qc','Actual TapeStation HSD1000 analytical allocation; amountunassigned',K.MOLECULAR);r=g.product(key+'_stock','Stock remaining after actual QC allocation',K.MOLECULAR);g.react(key,'Actual library size/concentration QC allocation',S.SPLIT,((b.x,R.TEMPLATE),),(q,r));g.mark_terminal(Terminal(q.id,Disposition.QC,'Analytical QC allocation'));b.x=r

def prepare_loaded(b):
    g=b.g;dups=[]
    for name in ('A','B'):
        rev=g.input('me_'+name.lower()+'_rev','Actual200ul200uM Tn5MErev /5Phos/ order',K.OLIGO,state=adapter('rev'));long=g.input('me_'+name.lower()+'_long','Actual200ul200uM Tn5ME-'+name,K.OLIGO,state=adapter(name));d=g.product('me_'+name.lower()+'_duplex','Actual400ul preannealed100uM duplex inventory',K.MOLECULAR,state=duplex(name),chemical_state={'whole_long':ME[name],'whole_rev':ME['rev'],'rev_5prime':'phosphate','whole_occupied_transposome':'unassembled'})
        g.react('me_'+name.lower()+'_anneal','Ref24 explicit MEDS oligo annealing in10mMTris8/50mMNaCl/1mMEDTA',S.HYBRIDIZATION,((rev,R.TEMPLATE),(long,R.STRUCTURAL),b.cond('me_'+name.lower()+'_buffer','Source annealing buffer; thermal programme unassigned')),(d,));dups.append(allocate(g,'me_'+name.lower()+'_take',d,1,'Actual8ul100uM duplex allocation','Untransferred annealed duplex, nominal392ul')[0])
    prot=g.input('pa_tn5','Actual100ul purified FLAG-removed pA-Tn5 fusionprotein from explicitly citedref24; upstreamexpression/purification outside selected component scope',K.CATALYST)
    loaded=g.product('loaded_transposome','Actual116ul transposome loading preparation; complete occupancy/end-complex stoichiometry unassembled',K.MOLECULAR,chemical_state={'known_adapters':ME,'rev_5prime':'phosphate','protein':'FLAG-removed pA-Tn5','occupied_molecular_geometry':'unassembled'})
    g.react('load_transposome','Ref24 actual8ul ME-A/rev +8ul ME-B/rev +100ul pA-Tn5,1hRT rotator/store-20C',S.MAGIC,tuple((x,R.TEMPLATE)for x in dups)+((prot,R.CATALYST),),(loaded,),reason='Specific loadedcomponent inheritance explicitly established byChem-mapref24; both duplex components declared, occupiedcomplex/stoichiometry unassembled. No G4 cellular assay/whole library recipe inherited.')
    return allocate(g,'loaded_take',loaded,1,'Actual5ul pA-Tn5 transposome allocation forfiveChem-map samples; main specifies protein2uM','Untransferred loadedtransposome stock; historical concentration-adjustment interface unassembled')[0]

def protocol_graph():
    b=B('Chem-map JQ1-biotin recruitment: five K562 technical samples',AUTHORITY,'ActualK562 sourcecell suspension forfive100ul samples at1500cells/ul; chosen JQ1 sodium-buffer route');g=b.g
    b.op('fix','0.1%formaldehydePBS2minRT',S.MAGIC,(b.stock('formaldehyde','Actual0.1%formaldehydePBS allocation',K.MOLECULAR,R.STRUCTURAL),))
    b.op('quench','Glycine75mMfinal',S.MAGIC,(b.stock('glycine','Actual glycine75mMfinal',K.MOLECULAR,R.STRUCTURAL),))
    b.clean('cells_wash','600g4min; washbuffer20mMHEPES7.5/150mMNaCl/.5mMspermidine/PIC')
    cells=allocate(g,'cell_allocations',b.x,5,'Actual100ul fixedcell sample','Untransferred sourcecell suspension; amountunassigned')
    # Prepared primary complex is a real source batch, with actual BSA-containing buffer.
    probe=g.input('jq1_biotin','Actual20ul10uM JQ1-btn small-molecule probe in antibodybuffer; NOT nucleotide',K.MOLECULAR);ab=g.input('antibiotin','Actual16.7ul~10uM D5A7rabbit anti-biotin5597',K.MOLECULAR)
    primary=g.product('primary_complex','Actualfive-sample probe-primaryAb preparation; occupiedcomplex stoichiometry unassembled',K.MOLECULAR)
    g.react('primary_assemble','20ulprobe+16.7ulAb+200ulantibodybuffer,ice1h; +300ulantibodybuffer',S.MAGIC,((probe,R.TEMPLATE),(ab,R.STRUCTURAL),b.stock('primary_buffer','Actual500ul antibodybuffer total, including2mMEDTA/.1%BSA/.05%digitonin/washbuffer',K.MOLECULAR,R.STRUCTURAL)),(primary,),reason='Sourceprecomplex smallmolecule/antibiotin; source nominal0.4uM alongside536.7ul summedingredients retained, not silently corrected or DNAencoded.')
    primaryparts=allocate(g,'primary_allocations',primary,5,'Actual100ul probe-primarycomplex sample allocation','Untransferred preparedprimarycomplex; nominal36.7ul if source additive volumes exhaustive')
    loaded=prepare_loaded(b)
    second=g.product('secondary_complex','Actualfive-sample secondaryAb-pA-Tn5 complexpreparation; wholeoccupiedgeometry unassembled',K.MOLECULAR,chemical_state={'known_adapter_orders':ME,'rev_5prime':'phosphate','components':'ActualsecondaryAb +loadedFLAG-removedpA-Tn5; carries both source duplex classes'})
    g.react('secondary_assemble','2.5ul~8uM secondaryABIN101961+5ul2uMpA-Tn5+200ulDig300,ice1h; +300ulantibodybuffer',S.MAGIC,((loaded,R.TEMPLATE),b.stock('secondary_ab','Actual2.5ul secondaryAbABIN101961',K.MOLECULAR,R.STRUCTURAL),b.stock('secondary_buffers','Actual200ulDig300 +300ulantibodybuffer includingBSA; source2:1Ab:pA-Tn5 nominalratio',K.MOLECULAR,R.STRUCTURAL)),(second,),reason='Known sourceprecomplex and loadedadapterinventory, not independentfreeTn5 addition to eachcell. Actual507.5ul batch, occupancy unassembled.')
    seconds=allocate(g,'secondary_allocations',second,5,'Actual100ul secondaryAb-transposome aliquot','Untransferred preparedsecondarycomplex; nominal7.5ul if additivevolumes exhaustive')
    bead=g.input('con_a','Actual50ul ConAbeads forfiveJQ1samples',K.MOLECULAR);washed=g.product('con_a_washed','Actual washedConA stock resuspended75ul bindingbuffer',K.MOLECULAR)
    g.react('bead_prepare','Two1ml bindingbuffer washes; resuspend75ul',S.MAGIC,((bead,R.TEMPLATE),b.cond('bead_buffer','20mMHEPES7.5/10mMKCl/1mMCaCl2/1mMMnCl2')),(washed,),reason='Sourcebead wash/resuspension; beadrecovery/washcarryoverunassigned.')
    beads=allocate(g,'bead_allocations',washed,5,'Actual10ul prewashedConA allocation','Untransferred washedConAstock; nominal25ul')
    shares=[]
    for i,cell in enumerate(cells):
        k='t'+str(i+1);b.x=cell
        b.op(k+'_con_a','100ulcells+10ulConA,25C10min600rpm',S.MAGIC,((beads[i],R.STRUCTURAL),))
        b.clean(k+'_cell_wash','Twowashbuffer washes; retainbead-boundcells')
        b.op(k+'_probe_bind','100ul probe-primaryprecomplex,4Covernight600rpm',S.MAGIC,((primaryparts[i],R.STRUCTURAL),),reason='ActualJQ1-biotin drug recruitment throughanti-biotin; not directanti-BRD4 antibody binding or nucleotidebarcode.')
        b.clean(k+'_digwash','Three500ulDigwash washes,.05%digitonin')
        b.op(k+'_secondary_bind','100ul secondary-Tn5precomplex,25C1h600rpm',S.MAGIC,((seconds[i],R.STRUCTURAL),),reason='ActualsecondaryAb recruits loadedpA-Tn5; molecule-specificoccupancy/adapter-genomegeometryunassembled.')
        b.clean(k+'_dig300','Three500ulDig300 washes')
        b.op(k+'_tagment','300ultagmentationbuffer10mMMgCl2inDig300,37C1h600rpm',S.MAGIC,(b.cond(k+'_mg','ActualMgactivationprogramme; existingboundTn5 catalyst carried in priorcomplex, not newexternalenzyme'),),reason='Source adapterintegration near drug-boundchromatin; knownadapterduplexes/5P separate fromunassembledoccupiedgenomicstrand-transferproducts/offproducts. No defaultTn5junction orgenomicsequence borrowed.')
        b.clean(k+'_taps','Two500ulTAPS/.2mMEDTA washes')
        pk=components(g,k+'_extraction','Actual150ul extractionbuffer .5mg/mlProteinaseK/.5%SDS/10mMTris8',(('pk','ProteinaseK ofsame150ul extractionbuffer',K.CATALYST,R.CATALYST),('other','ActualSDS/Tris components ofsamebuffer',K.MOLECULAR,R.STRUCTURAL)))
        b.op(k+'_extract','55C1h800rpm extraction',S.DEPROTEINIZATION,pk)
        b.clean(k+'_pci','Actual150ulPCI/phase-lock16000g3min; RETAINtopaqueous layer')
        b.clean(k+'_chloroform','Actual150ulchloroform/invert10x/spin; RETAINtopaqueous layer')
        b.op(k+'_precipitate','Actual6ul5MNaCl+375ulcoldEtOH,-20Covernight',S.MAGIC,(b.stock(k+'_precip_supply','Actual6ul5MNaCl and375ulEtOH supplies; separatepackaging notassigned',K.MOLECULAR,R.STRUCTURAL),))
        b.clean(k+'_pellet','21130g4C30min, retainDNApellet;1mlcold100%EtOHwash/spin2min/airdry')
        rn=components(g,k+'_resuspension','Actual25ul ElutionBuffer1 allocation:Tris8/1mMEDTA/1:400RNaseA',(('rnase','RNaseA component ofsame25ulbuffer',K.CATALYST,R.CATALYST),('buffer','Same25ulTris/EDTA/water inventory',K.MOLECULAR,R.STRUCTURAL)))
        b.op(k+'_rnase','25ulresuspension,37C10min800rpm',S.DIGESTION,rn)
        b.split(k+'_pcr_template','Actual21ulof25ulextractedDNA','Actual21ulPCRtemplate','Untransferred4ulextractedDNA')
        q5=components(g,k+'_q5','Actual25ul NEBNextUltraIIQ5 2xmastermixM0544',(('pol','Q5DNApolymerase ofsame25ulmastermix',K.CATALYST,R.CATALYST),('nt','Nucleotide components ofsame25ulmastermix',K.MOLECULAR,R.STRUCTURAL),('buffer','Buffercomponents ofsame25ulmastermix',K.CONDITION,R.CONDITION)))
        ps=[]
        for end in ('1','2'):
            sc=Scene();sc.strand('v2Ad'+end,[Segment('Wholev2Ad'+end,'[whole_v2Ad'+end+'_unassembled]',placeholder=True,feature=feature('chemmap.sample.'+end,'sample_index','unknown',note='Sourceuniquelybarcoded v2Ad1/Ad2; wholebodies/indexwidth/selectedindex values unassembled, not UMI'))]);ps.append(b.stock(k+'_ad'+end,'Actual2ul10uMv2Ad'+end+'.x primer, explicitlycitedref40wholeorder source ordinaryunread',K.OLIGO,R.STRUCTURAL,state=sc))
        b.op(k+'_pcr','50ulPCR:72C5min/98C30s;10x98C10s/63C10s;72C1min',S.PCR,q5+tuple(ps),geometry=PCRGeometry.unavailable(unknown_scene(),reason='KnownsourceactualPCR inputs/cycles; exactref40primers/selectedbarcodes andoccupiedadapter-genomegeometry unassembled. No canonicalNextera primer borrowed.'))
        b.clean(k+'_pcr_clean','65ul1.3xAMPure10min/two80%EtOH/elute25ul10mMTris5min')
        qc(b,k+'_qc');b.split(k+'_pool_select','Balanced selectedlibrary contribution','Actual balancedpool contribution','Untransferred QCedlibrary stock');shares.append(b.x)
    pool=g.product('balanced_pool','Actual balancedfive-library samplepool',K.MOLECULAR);g.react('balanced_pool_step','Pool actualfive selectedsamplelibrary contributions',S.POOL,tuple((x,R.TEMPLATE)for x in shares),(pool,));b.x=pool
    # Explicit first supernatant and second bead-bound physical fates.
    y=g.product('first_supernatant','Actual RETAINED0.4xsupernatant',K.MOLECULAR);w=g.product('first_beads','Actual excludedfirstbead-boundfraction',K.WASTE);g.react('first_size','0.4xAMPure15minRT; RETAINsupernatant',S.CLEANUP,((b.x,R.TEMPLATE),b.stock('first_ampure','Actual0.4xAMPure allocation',K.SELECTION,R.SELECTION)),(y,w));g.discard(w,'Actualfirstbeadboundfraction excluded');b.x=y
    b.clean('second_size','Add1.3xAMPuretosupernatant15min,RETAINbeads/two80%EtOH/elute40ul10mMTris; ratio-referencevolumeconvention unassigned')
    b.split('pool_load','Selected pooledlibrary loading allocation','Actual finalpooledloadingaliquot','Untransferred40ul pooledlibrarystock afterselectedloadingallocation')
    lib=g.product('library','ActualselectedJQ1Chem-map pooledloadingmixture; complete genomicproducts/readlandingunassembled',K.LIBRARY);g.react('load','NextSeq500HighOutput paired-end36',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(read_lengths={'R1':36,'R2':36},limitation='ActualChem-map PE36 run. Source-knownloadedadapterduplexes do not resolvefullPCRprimers/selectedindexes/occupiedgenomejunctions/readlanding.'))
    g.add_source_item(SourceItem('scope',AUTHORITY,'Privatepositivepartial molecular graph. Only source-qualifiedloadedtransposomecomponent inherited, no G4buffer/antibody/cell/10xrecipe.'))
    g.add_source_item(SourceItem('volume_scope','Primaryrecipe536.7ul summedingredients plusstated0.4uM; secondaryrecipe507.5ul; bothforfive100ul samples.','Actual sourcequantities retained, nominalexcessesnotroundedaway; volumeadditivity/nominalconcentration applicabilityunassigned.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Canonical component reuse matches independent whole source-order fixture',ME,SOURCE_ME)
    for name in ('A','B'):
        sc=duplex(name);c('Whole '+name+' transferstrand retained',''.join(x.top for x in sc.strands[name].segs),ME[name]);c('Complementary19bp core '+name,revcomp(ME['rev']),ME[name][-19:]);c('Duplex aliquot retains molecular state '+name,g.materials['me_'+name.lower()+'_take_0'].state is g.materials['me_'+name.lower()+'_duplex'].state,True)
    c('Known phosphorylation survives selected load',g.materials['loaded_take_0'].chemical_state['rev_5prime'],'phosphate')
    c('Five actual PCR categories',len([s for s in g.steps.values()if s.kind==S.PCR]),5)
    c('Five allocated samples use each actual batch once',len(g.steps['primary_allocations'].outputs),6)
    c('Balanced pool uses selected fractions',[x.material for x in g.steps['balanced_pool_step'].inputs],['t'+str(i)+'_pool_select_selected'for i in range(1,6)])
    c('First size fraction keeps supernatant',g.steps['second_size'].inputs[0].material,'first_supernatant')
    c('Whole final sequencingpool not implicitly loaded',g.steps['load'].inputs[0].material,'pool_load_selected')
    c('Actual PE36',g.sequencing['library'].read_lengths,{'R1':36,'R2':36})
    for i in range(1,6):
        k='t'+str(i);c('PCR uses21ulselectedDNA '+k,g.steps[k+'_pcr'].inputs[0].material,k+'_pcr_template_selected');c('Q5catalyst from samephysicalstock '+k,g.producer_of(k+'_q5_pol').id,k+'_q5_parts');c('ProteinaseK from same150ulextractionbuffer '+k,g.producer_of(k+'_extraction_pk').id,k+'_extraction_parts')
    root=sources_root()/Path(__file__).resolve().parents[1].name;p=root/'s41598-021-02943-3_supp_41598_2021_2943_MOESM1_ESM.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41598-021-02943-3')):
        t=''.join(p.read_text().split())
        for name,s in ME.items():c('Whole actual ref24 source ordering '+name,s in t,True)
        c('Original explicit phosphorylation','[phos]'+ME['rev']in t,True)
    p=root/'s41587-022-01636-0_PMC10497411.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41587-022-01636-0')):
        t=p.read_text();c('Source explicit component delegation','pA–Tn5 assembled with DNA adapters were prepared as described in ref. 24' in t,True);c('Actual recipe Q5 master identity','NEBNext Ultra II Q5 2x PCR master mix' in t,True)
    c.report()
if __name__=='__main__':selftest()
