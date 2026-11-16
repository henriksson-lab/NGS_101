from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature,complement_segments,MolecularState,Workflow,workflow_panel
from dataclasses import replace
from protocol_model import role
from illumina import P5,P7
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SourceCoverage,SequencingPlan,Terminal,Disposition
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
AUTHORITY='Published2023 Chem-map PMC10497411 Methods131–159/actualFig1PDFp2; explicitly delegated ref24 publishedPMC8654944 Methods87 +actualSI1p6–7 lines143–151 for loadedpA-Tn5 component ONLY. Publishedref40Nature2015Table1MOESM36 verifiesindexedprimerfamilies in a boundedfollowup. No G4/scATACcellularrecipe inherited; actualindexmembers/occupiedgenomicproductsunassembled.'
SOURCE_ME={'rev':'CTGTCTCTTATACACATCT','A':'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG','B':'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'}
# The three complete orders are verified in delegated ref24 SI1 p6–7;
# reuse only these matching canonical strings, not a Nextera reaction recipe.
ME={'rev':nextera.ME_RC,'A':nextera.ADAPTOR_S5,'B':nextera.ADAPTOR_S7}

# Chem-map Methods explicitly cites Buenrostro2015 ref40 for these v2 families.
# Published Nature Table1 MOESM36 original XLSX establishes the constant bodies
# and eight-base ordered index windows. Actual Chem-map .x members remain unknown.
PRIMER_FOOT={'1':nextera.S5+nextera.ME[:10],'2':nextera.S7+nextera.ME[:7]}
REF40_URL='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature14590/MediaObjects/41586_2015_BFnature14590_MOESM36_ESM.xlsx'


def indexed_primer_scene(end,index=None):
    if end not in {'1','2'}:raise ValueError('Source v2 family is Ad1 or Ad2')
    if index is not None and (len(index)!=8 or any(base not in 'ACGT'for base in index)):
        raise ValueError('A selected source-family ordered index must have eight DNA bases')
    seq=index if index is not None else 'N'*8
    sc=Scene();sc.strand('v2Ad'+end,[Segment('P5'if end=='1'else'P7',P5 if end=='1'else P7),
        Segment('ordered '+('i5'if end=='1'else'i7')+' window',seq,'cbc',placeholder=index is None,
            feature=feature('chemmap.sample.'+end,'sample_index','fixed'if index is not None else'unknown',
                note='Ref40 Table1 eight-base ordered primer index; Chem-map selected member and acquired index-read convention unassigned')),
        Segment('source v2Ad'+end+' foot',PRIMER_FOOT[end])],label='Source v2Ad'+end+'.x family')
    return sc


def adapter(name):
    sc=Scene();sc.strand(name,[Segment(name,ME[name])],mod5='phosphate'if name=='rev'else'');return sc

def duplex(name):
    sc=Scene();s=ME[name];core=revcomp(ME['rev']);sc.strand(name,[Segment(name+' tail',s[:-len(core)]),Segment('ME core',core)]);sc.anneal('rev',[Segment('ME rev',ME['rev'])],to=name,pair=('ME rev','ME core'),mod5='phosphate');return sc

def allocate(g,key,stock,n,label,remaining):
    xs=[g.product(key+'_'+str(i),label+' '+str(i+1),K.MOLECULAR,state=stock.state,chemical_state=stock.chemical_state)for i in range(n)];r=g.product(key+'_remainder',remaining,K.MOLECULAR,state=stock.state,chemical_state=stock.chemical_state);g.react(key,'Actual selected aliquots; no duplicated whole stock',S.SPLIT,((stock,R.TEMPLATE),),tuple(xs)+(r,));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remaining));return xs

def components(g,key,label,classes):
    st=g.input(key+'_stock',label,K.MOLECULAR);xs=[g.product(key+'_'+name,caption,kind)for name,caption,kind,role in classes];g.react(key+'_parts','Copresent components of SAME actual physical stock; no extra fresh doses',S.GENERIC,((st,R.TEMPLATE),),tuple(xs),changes_molecule=False);return tuple((x,row[3])for x,row in zip(xs,classes))

def qc(b,key):
    g=b.g;q=g.product(key+'_qc','Actual TapeStation HSD1000 analytical allocation; amountunassigned',K.QC);r=g.product(key+'_stock','Stock remaining after actual QC allocation',K.MOLECULAR);g.react(key,'Actual library size/concentration QC allocation',S.SPLIT,((b.x,R.TEMPLATE),),(q,r));g.mark_terminal(Terminal(q.id,Disposition.QC,'Analytical QC allocation'));b.x=r

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
            sc=indexed_primer_scene(end);ps.append(b.stock(k+'_ad'+end,'Actual2ul10uMv2Ad'+end+'.x source-defined family; selected8ntmember unassigned',K.OLIGO,R.STRUCTURAL,state=sc))
        b.op(k+'_pcr','50ulPCR:72C5min/98C30s;10x98C10s/63C10s;72C1min',S.PCR,q5+tuple(ps),geometry=PCRGeometry.unavailable(unknown_scene(),reason='KnownsourceactualPCR inputs/cycles; Ref40 v2 constant bodies and8ntorderedindex windows verified; selectedmembers/occupiedadapter-genomegeometry andwholeamplicons unassembled.'))
        b.clean(k+'_pcr_clean','65ul1.3xAMPure10min/two80%EtOH/elute25ul10mMTris5min')
        qc(b,k+'_qc');b.split(k+'_pool_select','Balanced selectedlibrary contribution','Actual balancedpool contribution','Untransferred QCedlibrary stock');shares.append(b.x)
    pool=g.product('balanced_pool','Actual balancedfive-library samplepool',K.MOLECULAR);g.react('balanced_pool_step','Pool actualfive selectedsamplelibrary contributions',S.POOL,tuple((x,R.TEMPLATE)for x in shares),(pool,));b.x=pool
    # Explicit first supernatant and second bead-bound physical fates.
    y=g.product('first_supernatant','Actual RETAINED0.4xsupernatant',K.MOLECULAR);w=g.product('first_beads','Actual excludedfirstbead-boundfraction',K.WASTE);g.react('first_size','0.4xAMPure15minRT; RETAINsupernatant',S.CLEANUP,((b.x,R.TEMPLATE),b.stock('first_ampure','Actual0.4xAMPure allocation',K.SELECTION,R.SELECTION)),(y,w));g.discard(w,'Actualfirstbeadboundfraction excluded');b.x=y
    b.clean('second_size','Add1.3xAMPuretosupernatant15min,RETAINbeads/two80%EtOH/elute40ul10mMTris; ratio-referencevolumeconvention unassigned')
    b.split('pool_load','Selected pooledlibrary loading allocation','Actual finalpooledloadingaliquot','Untransferred40ul pooledlibrarystock afterselectedloadingallocation')
    lib=g.product('library','ActualselectedJQ1Chem-map pooledloadingmixture; complete genomicproducts/readlandingunassembled',K.LIBRARY);g.react('load','NextSeq500HighOutput paired-end36',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(read_lengths={'R1':36,'R2':36},limitation='ActualChem-map PE36 run. Source-knownloadedadapterduplexes do not resolveselectedPCRindexmembers/occupiedgenomejunctions/readlanding.'))
    g.add_source_item(SourceItem('ref40_indexed_primers',
        'Chem-map published Methods cites Buenrostro2015(ref40), published Nature Table1 original MOESM36 XLSX.',
        'v2Ad1/v2Ad2 constant bodies and8ntorderedindex windows verified; actualChem-map.x member assignments/read-index masks remain unassembled. Ad2 row names can use the reverse complement of the ordered index; no names are substituted for molecular bases.'))
    g.cover_source_item(SourceCoverage('ref40_indexed_primers','modelled','t1_pcr'))
    declare_page_scope(g)
    g.add_source_item(SourceItem('scope',AUTHORITY,'Source-qualified partial molecular graph. Only source-qualifiedloadedtransposomecomponent inherited, no G4buffer/antibody/cell/10xrecipe.'))
    g.add_source_item(SourceItem('volume_scope','Primaryrecipe536.7ul summedingredients plusstated0.4uM; secondaryrecipe507.5ul; bothforfive100ul samples.','Actual sourcequantities retained, nominalexcessesnotroundedaway; volumeadditivity/nominalconcentration applicabilityunassigned.'))
    g.cover_source_item(SourceCoverage('scope','modelled','source'))
    g.cover_source_item(SourceCoverage('volume_scope','modelled','primary_assemble'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Typed publication gate',g.publication_problems(),[])
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
    import zipfile,xml.etree.ElementTree as E
    name='nature14590_41586_2015_BFnature14590_MOESM36_ESM.xlsx'
    p=root/name
    if not p.exists():p=sources_root()/'doi_10.1038+nature14590'/name
    if have(c,Source(p,REF40_URL)):
        manifest=p.parent/'MANIFEST.tsv'
        identity=manifest.is_file()and any(row.split('\t')[:2]==[p.name,REF40_URL]for row in manifest.read_text().splitlines())
        c('Ref40 cached file belongs to the exact published source URL',identity,True)
        if identity:
            with zipfile.ZipFile(p)as z:
                ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
                ss=[''.join(x.itertext())for x in E.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)]
                original={}
                for row in E.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//s:row',ns):
                    cells={}
                    for cell in row.findall('s:c',ns):
                        v=cell.find('s:v',ns)
                        if v is not None and cell.get('t')=='s':cells[cell.get('r')[0]]=ss[int(v.text)]
                    if cells.get('A','').startswith(('v2_Ad1.','v2_Ad2.')):original[cells['A']]=cells['B']
            for name,end,index in (('v2_Ad1.1_TAGATCGC','1','TAGATCGC'),('v2_Ad2.1_TAAGGCGA','2','TCGCCTTA')):
                c('Original example row name and ordered index '+name,original[name],''.join(x.top for x in indexed_primer_scene(end,index).strands['v2Ad'+end].segs))
            for end,prefix in (('1',P5),('2',P7)):
                members=[seq for name,seq in original.items()if name.startswith('v2_Ad'+end+'.')]
                c('Original Table1 contains indexed v2Ad'+end+' members',bool(members),True)
                c('All original Table1 v2Ad'+end+' whole orders match declared body/window',
                  all(seq==''.join(x.top for x in indexed_primer_scene(end,seq[len(prefix):len(prefix)+8]).strands['v2Ad'+end].segs)for seq in members),True)
    c.report()

# Curator-authored source inventory. The five values below are the source's
# five-sample preparation, not a traversal of whatever the graph contains.
def declare_page_scope(g):
    g.declare_reviewed_scope((
        'source','formaldehyde','fix','glycine','quench','cells_wash_condition',
        'cells_wash','cell_allocations','jq1_biotin','antibiotin','primary_buffer',
        'primary_assemble','primary_allocations','me_a_rev','me_a_long','me_a_buffer',
        'me_a_anneal','me_a_take','me_b_rev','me_b_long','me_b_buffer','me_b_anneal',
        'me_b_take','pa_tn5','load_transposome','loaded_take','secondary_ab',
        'secondary_buffers','secondary_assemble','secondary_allocations',
        'con_a','bead_buffer','bead_prepare','bead_allocations','balanced_pool_step',
        'first_ampure','first_size','second_size_condition','second_size',
        'pool_load','load','library',
    ),AUTHORITY)
    sample_targets=(
        'con_a','cell_wash_condition','cell_wash','probe_bind','digwash_condition',
        'digwash','secondary_bind','dig300_condition','dig300','mg','tagment',
        'taps_condition','taps','extraction_stock','extraction_parts','extract',
        'pci_condition','pci','chloroform_condition','chloroform','precip_supply',
        'precipitate','pellet_condition','pellet','resuspension_stock',
        'resuspension_parts','rnase','pcr_template','q5_stock','q5_parts',
        'ad1','ad2','pcr','pcr_clean_condition','pcr_clean','qc','pool_select',
    )
    for sample in (1,2,3,4,5):
        g.declare_reviewed_scope(tuple(f't{sample}_{target}'for target in sample_targets),
                                'Published Chem-map Methods131–159: each of five actual technical samples; '+AUTHORITY)

PAGE_FILE = 'chem-map.html'

def molecular_class(kind):
    sc=Scene()
    if kind=='library':
        first=[replace(x,inferred=True)for x in indexed_primer_scene('1').strands['v2Ad1'].segs]
        second=[replace(x,inferred=True)for x in complement_segments(indexed_primer_scene('2').strands['v2Ad2'].segs)]
        parts=first+[role('remaining ME / genomic interval; unassembled',inferred=True)]+second
        sc.strand('reference',parts,label='INFERRED A/B PCR reference')
        sc.note('reference','Conditional intended primer-derived boundaries only. Complete occupied tagmentation/repair products, selected index members and read landing remain unassembled; this is not every reaction product or a source-validated whole genomic library.')
        return sc.rows()
    body=role('source genomic DNA interval',inferred=True)
    parts=[body]
    if kind=='tagged':parts=[role('transferred adapter end; A/B choice and geometry unassembled',inferred=True),body,role('other transferred end; unresolved',inferred=True)]
    if kind=='tagged':
        sc.strand('tagged',parts,label='INFERRED tagged DNA')
        sc.note('tagged','Reference component view only. An occupied complementary strand is not assembled here: exact transferred-strand junctions, stagger, repair extent and fragment ends remain unresolved. Known reagent duplexes do not establish a gap-filled genomic product.')
    else:
        bottom=complement_segments(parts);sc.strand('dna',parts,label='INFERRED genomic DNA class')
        sc.anneal('other',bottom,to='dna',pair=(bottom[-1].name,parts[0].name),unpaired=tuple(x.name for x in bottom),label='DNA complement class')
    return sc.rows()


def molecular_walkthroughs():
    reagent=Workflow(MolecularState('Known source ME adapter orders; separate annealing preparations',tuple(adapter('A').rows()+adapter('B').rows()+adapter('rev').rows())),track_graph=False)
    reagent.react('1R · Prepare the two source ME duplexes independently',duplex('A').rows()+duplex('B').rows(),
        note='Explicitly delegated ref24 component only: anneal 200 µl 200 µM MErev with the same amount of ME-A, and separately with ME-B. The 19-bp cores are computed from whole source oligos; MErev has a 5′ phosphate. The two preparations are not one covalent adapter or a G4 cellular recipe.')
    reagent.react('2R · Load adapter duplex aliquots onto pA–Tn5',duplex('A').rows()+duplex('B').rows(),
        note='Use 8 µl of each prepared duplex with 100 µl pA–Tn5 for 1 h at room temperature; retain the source stock boundary. These rows show loaded-adapter inventory, not physically free adapters or a measured dimer occupancy. Exact occupied transposome/genomic geometry remains unassembled. Only the selected 5 µl transposome contribution enters the five-sample Chem-map precomplex.')
    sample=Workflow(MolecularState('K562 cells before selected JQ1-biotin Chem-map preparation',tuple(molecular_class('genomic'))),track_graph=False)
    def step(title,kind,note):sample.react(title,molecular_class(kind),name=title,note=note)
    step('1 · Fix K562 cells and capture five separate aliquots on ConA beads','genomic',
         '0.1% formaldehyde for 2 min at room temperature, final75 mM glycine; collect and wash. Selected JQ1 route uses sodium buffers and 1,500 cells/µl. Five100-µl samples receive10 µl prepared ConA beads each,25°C10 min, then wash. G4 potassium-buffer and doxorubicin alternatives remain separate; cell capture is not DNA adapter ligation.')
    step('2 · Bind the preassembled JQ1-biotin / primary anti-biotin complex','genomic',
         'Five-sample batch:20 µl10 µM probe +16.7 µl approximately10 µM rabbit anti-biotin +200 µl antibody buffer,1 h on ice; add300 µl antibody buffer. Incubate each cell sample in100 µl precomplex overnight at4°C. Source nominal0.4 µM and summed536.7-µl recipe are both retained; their concentration/volume conventions are not silently normalized. The small-molecule drug is not an oligo or sequence barcode.')
    step('3 · Recruit the separate secondary-antibody / loaded pA–Tn5 complex','genomic',
         'Five-sample batch:2.5 µl secondary antibody +5 µl loaded transposome +200 µl Dig-300,1 h on ice; add300 µl antibody buffer. Source total507.5 µl is retained. After drug-complex washes, incubate cells in100 µl secondary complex at25°C1 h and wash Dig-300 three times. Recruitment is through the drug/antibody system, not direct anti-BRD4 CUT&Tag or fresh independent Tn5 additions to every cell.')
    step('4 · INFERRED — Mg²⁺ activates nearby DNA tagmentation','tagged',
         '300 µl tagmentation buffer,10 mM MgCl₂ in Dig-300,37°C1 h. Source-loaded adapter classes can enter nearby genomic DNA. The class view does not force one A/B orientation, complete occupancy, a measured cut position or filled-in ends. Adapter bodies are known reagents; the exact genomic strand-transfer/repair products remain unassembled.')
    step('5 · Recover each sample independently and select 21 µl DNA for PCR','tagged',
         'TAPS washes, then150 µl proteinaseK/SDS extraction buffer,55°C1 h. Phenol/chloroform phase extraction, chloroform wash and source ethanol precipitation recover DNA. Resuspend in25 µl elution buffer with RNaseA; transfer21 µl to PCR and keep4 µl stock explicit. The five samples are not physically pooled before independent indexing.')
    step('6 · INFERRED — source v2Ad1.x / v2Ad2.x PCR adds the indexed boundaries','library',
         '50 µl:21 µl DNA +25 µl NEBNext UltraII Q5 2×mix +2 µl each10 µM source-indexed v2 primer. Programme72°C5 min /98°C30 s;10 cycles98°C10 s /63°C10 s;72°C1 min. Direct ref40 original Table1 now establishes both constant bodies and eight-base ordered index windows; actual study members and full genomic/repair/read geometry remain unassembled. The PCR-reference diagram is a conditional boundary class, not a complete source-validated product.')
    step('7 · Clean and QC each library; pool selected balanced contributions','library',
         '1.3×AMPure (65 µl),10 min, two80% ethanol washes, elute25 µl Tris. QC consumes separate analytical allocations. Select balanced library contributions from the five actual samples and retain untransferred stocks; the final sample pool is a physical combination, unlike computational duplicate/read filtering.')
    step('8 · Two-stage size selection and NextSeq500 PE36 loading','library',
         'First0.4×AMPure retains supernatant and excludes its bead-bound fraction. The following1.3×addition captures DNA on beads; its ratio-reference convention remains unassigned. Elute40 µl, select the actual loading aliquot and sequence NextSeq500 High Output paired-end36. Selected index values, acquired index-read orientation/masks and full primer landing remain unresolved.')
    return reagent,sample


def render_page():
    from page import head,info,caveat,table
    from chemdraw import panel
    from protocol_views import assemble_page,WorkflowSection,PlanLimitationSection
    from chemistry_validation import check_graph,Status
    g=protocol_graph();g.require_publishable()
    overview=(
        '<div class="wrap" data-protocol-complete="true"><h1>Chem-map</h1>'
        '<p>Map where a small molecule binds chromatin by recruiting adapter-loaded pA–Tn5 through antibodies. '
        'This page follows JQ1-biotin in K562 cells, including the five-sample source preparation.</p>'
        +info('Defining source: <a href="https://doi.org/10.1038/s41587-022-01636-0">Yu et al., Nature Biotechnology (2023)</a>. '
              'The publication explicitly delegates the loaded pA–Tn5 component to '
              '<a href="https://doi.org/10.1038/s41598-021-02943-3">Hui et al. (2021), Supplementary Information</a>. '
              'Only that component is inherited; the earlier G4 cellular assay is not substituted for Chem-map.')
        +caveat('Whole ME-A, ME-B and phosphorylated MErev oligos are known. Occupied transposome/genome junctions, '
                'the selected v2Ad1.x/Ad2.x index members and precise read landing remain unassembled. Their constant bodies and eight-base ordered index windows are verified from the directly cited ref40 Table1. '
                'A passing workflow gate accounts for source operations and material fates; molecular validation remains scoped.')
        +'<h2>How recruitment becomes a DNA library</h2>'
        +table(('Stage','Source-defined operation','Molecular meaning'),(
            ('Prepare cells','Fix K562 cells and allocate five 100 µl samples; bind washed ConA beads.','The selected JQ1 route uses sodium buffers. G4 and doxorubicin branches are separate.'),
            ('Prepare probe complex','Preassemble JQ1-biotin with anti-biotin antibody.','The biotinylated drug is a small molecule, not an oligo or sequence barcode.'),
            ('Recruit Tn5','Bind the probe complex, then a separately prepared secondary-antibody/pA–Tn5 complex.','Drug–target binding recruits the transposome. This is not direct anti-BRD4 CUT&Tag.'),
            ('Tagment','Activate the bound transposome with MgCl₂.','Known loaded adapter classes enter nearby genomic DNA; full occupied products remain unassembled.'),
            ('Recover and PCR','Extract each sample independently, select 21 µl of its 25 µl eluate, and amplify for ten cycles.','Q5 master-mix components come from the same physical 25 µl allocation.'),
            ('Pool and sequence','QC each library, select balanced contributions and perform two-stage bead selection.','Only a loading aliquot is sequenced; unused sample and pool stocks remain explicit.'),
        ))
        +'<h2>Five-sample stock accounting</h2>'
        +caveat('The primary recipe lists 20 µl probe + 16.7 µl antibody + 200 µl buffer, then 300 µl buffer: '
                '536.7 µl in additive volumes, alongside a stated nominal 0.4 µM probe concentration. The secondary complex lists '
                '2.5 µl antibody + 5 µl transposome + 200 µl buffer, then 300 µl buffer: 507.5 µl. '
                'Five 100 µl allocations and untransferred stock remainders are retained; volume-additivity and nominal-concentration '
                'conventions are unassigned rather than silently rounded to 500 µl.')
        +info('The 0.4× AMPure step retains supernatant. The following 1.3× addition captures DNA on beads; its reference-volume '
              'convention remains unassigned. The source elutes the final pool in 40 µl before a selected loading aliquot.')
        +'<h2>Read setup</h2>'
        +table(('Item','Declaration'),(
            ('Platform','NextSeq500 High Output; paired-end 36 + 36 bases.'),
            ('PCR','21 µl DNA + 25 µl NEBNext Ultra II Q5 2× mix + 2 µl of each 10 µM indexed primer; ten cycles.'),
            ('Indices / primer landing','Ref40 v2Ad1/Ad2 constant bodies and 8-nt ordered index windows are known. Selected Chem-map members, acquired index-read orientation/masks and genomic landing geometry remain unassigned.'),
        ))
        +'<h2>Known loaded-adapter component</h2>'
        +'<p>The directly delegated component source specifies 200 µl of 200 µM MErev annealed with the same amount of ME-A or ME-B. '
        'Eight µl of each prepared duplex is mixed with 100 µl pA–Tn5 for one hour. MErev has a source-specified 5′ phosphate. '
        'The figures below show the actual whole orders and computed 19-base-pair cores, not a completed genomic library.</p>'
    )
    primer_views='<h2>Source v2 indexing-primer families</h2><p>Chem-map directly cites <a href="https://doi.org/10.1038/nature14590">Buenrostro et al. (2015), original Supplementary Table1</a>. The published workbook defines v2Ad1 as P5 / an eight-base i5 window / its fixed priming foot, and v2Ad2 as P7 / an eight-base i7 window / its different fixed foot. N8 below is an unassigned selected fixed index, not random synthesis or a UMI.</p>'
    primer_views+=''.join(panel(indexed_primer_scene(end).rows(),caption='Verified ref40 v2Ad'+end+'.x family; actual Chem-map index member unassigned')for end in ('1','2'))
    primer_views+='<p>For example, Table1 names v2_Ad1.1_TAGATCGC and orders TAGATCGC in its index window. It names v2_Ad2.1_TAAGGCGA but orders TCGCCTTA, the reverse complement, in that primer. These are ordered molecular sequences, not an automatic samplesheet convention; neither example is assigned to the five Chem-map samples.</p>'
    primer_views+=''.join(panel(indexed_primer_scene(end,index).rows(),caption='Ref40 family example only; not a selected Chem-map study index')for end,index in (('1','TAGATCGC'),('2','TCGCCTTA')))
    molecules=''.join('<h3>'+name+'</h3>'+panel(adapter(name).rows(),caption='Whole source oligo; other occupied complex/end states remain unassembled.')for name in ('rev','A','B'))
    molecules+=''.join(panel(duplex(name).rows(),caption='Computed ME-'+name+'/MErev duplex from the source orders; no occupied transposome is implied.')for name in ('A','B'))
    claims=check_graph(g)
    audit=('<h2>Molecular validation scope</h2>'+table(('Claim status','Computed count'),tuple((s.value,str(sum(c.status==s for c in claims)))for s in Status))
           +'<p>Local adapter-pairing claims do not establish all enzymatic transitions. The complete five-sample workflow below includes validation warnings, side fractions and untransferred stock fates.</p>')
    return assemble_page(g,[head('Chem-map — small-molecule recruitment and DNA library preparation'),overview,
        '<h2>Step by step</h2><p>INFERRED — the source-qualified classes below distinguish independent adapter loading and drug/antibody recruitment from tagmentation and later indexed PCR. Whole occupied transposome/genomic/repair products and actual index assignments remain unassembled.</p>',
        *(workflow_panel(w,caption=caption)for w,caption in zip(molecular_walkthroughs(),('Prepare the explicitly delegated loaded-adapter component','Five actual JQ1-biotin Chem-map samples through library pooling'))),
        primer_views,molecules,audit,'<details><summary>Complete preparation, fractions and model warnings</summary>',WorkflowSection('Complete preparation and material fates'),'</details>',PlanLimitationSection('library'),'<p><a href="01_chem-map.html">Research note and source boundaries</a></p></div>'])

if __name__=='__main__':selftest()
