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
from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2
from protocol_model import role
import re
SOURCE_ORDERS={'ad5_A':'TTCCAACTGGATACTGrGrA','ad5_G':'TTCCAACTGGATACTGrGrG','ad5_U':'TTCCAACTGGATACTGrGrU','ad5_C':'TTCCAACTGGATACTGrGrC','ad3':"/5'Phos/GTATCCAGTTGGAATTCTCGGGTGCCAAGG/3ddC/",'rt':"/5'phos/NNGATCGTCGGACTGTAGAACTCTGAACGTGT/iSp18/CCTTGGCACCCGAGAATT",'pcr_f':'AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA','pcr_r':'CAAGCAGAAGACGGCATACGAGAT[BARCODE]GTGACTGGAGTTCCTTGGCACCCGAGA'}
AD5_STEM='TTCCAACTGGATAC'
AD3_STEM=revcomp(AD5_STEM)
AD3_BODY=re.sub(r'/[^/]+/','',SOURCE_ORDERS['ad3'])
RT_LEFT=SOURCE_ORDERS['rt'].split('/iSp18/')[0].rsplit('/',1)[-1]
RT_BIND=SOURCE_ORDERS['rt'].split('/iSp18/')[1]
AUTHORITY='PublishedPMC7428014 main59/189/193 andFig1p4; publishedData1MOESM5 actualmodifiedlibraryorders andSI1. UntreatedHEK293 branch; noYAMAT/DM/mim constructor inheritance.'
SOURCE_OLIGO_EXCEPTIONS={'ad3':'Published Data1 whole order ends /3ddC/, a physical terminal dideoxy-C nucleotide. Model retains that additional C; generic source parser can remove the modification without retaining its nucleotide. Narrow raw-order extraction exception backed by whole original workbook fixture, not a chemistry or source-absence exemption.'}
def canon(text):
    # Reuse exact canonical sequences only when present in the source body.
    for x in sorted((P5,P7,TRUSEQ_READ1,TRUSEQ_READ2),key=len,reverse=True):
        if x in text:return [Segment('source fixed prefix',text[:text.index(x)])] * bool(text[:text.index(x)])+[Segment('source canonical body',x)]+([Segment('source fixed suffix',text[text.index(x)+len(x):])]if text[text.index(x)+len(x):]else[])
    return [Segment('source fixed body',text)]
def oligo_scene(name):
    s=Scene()
    if name.startswith('ad5_'):
        base=name[-1];parts=[Segment('ad5 stem',AD5_STEM),Segment('ad5 DNA TG','TG'),Segment('ad5 RNA rG','G'),Segment('ad5 terminal ribobase',base)]
        s.strand('oligo',parts);s.ligation_end('oligo','ad5 terminal ribobase',kind='splinted',polarity='3-prime',owner='adapter',text='Source3OH terminal ribonucleotide ligation leg; occupied tRNA junction unassembled')
    elif name=='ad3':
        s.strand('oligo',[Segment('ad3 stem',AD3_STEM),Segment('ad3 fixed tail',AD3_BODY[len(AD3_STEM):]),Segment('physical terminal ddC','C',note='Additional terminal dideoxycytidine nucleotide')],mod5='P',mod3='ddC')
        s.ligation_end('oligo','ad3 stem',kind='splinted',polarity='5-prime',owner='adapter',text='Source5P ligation leg, full tRNA junction unassembled')
    elif name=='rt':
        s.strand('oligo',[Segment('source RT NN','NN',placeholder=True,note='Purpose unassigned, not automatically UMI'),Segment('source RT left DNA',RT_LEFT[2:]),role('nonnucleotide iSp18, not18bases'),Segment('source RT binding DNA',RT_BIND)],mod5='P')
    elif name=='pcr_r':
        left,right=SOURCE_ORDERS[name].split('[BARCODE]');s.strand('oligo',canon(left)+[Segment('sample index','[BARCODE]',placeholder=True,feature=feature('quantm.sample','sample_index','unknown',note='Source template24compatibleindexes; actualvalue/lengthunassigned'))]+canon(right))
    else:s.strand('oligo',canon(SOURCE_ORDERS[name]))
    return s

def adapter_scene():
    s=Scene();s.strand('ad5',[Segment('ad5 stem',AD5_STEM),Segment('CCA/discriminator overhang','TGGA')]);s.anneal('ad3',[Segment('ad3 stem',AD3_STEM),Segment('ad3 unpaired tail',AD3_BODY[len(AD3_STEM):]),Segment('ad3 ddC','C')],to='ad5',pair=('ad3 stem','ad5 stem'),mod5='P',mod3='ddC',unpaired=('ad3 unpaired tail','ad3 ddC'))
    s.ligation_end('ad5','CCA/discriminator overhang',kind='splinted',polarity='3-prime',owner='adapter',text='Representative source rG/rA3OH leg; maturetRNA partner unassembled')
    s.ligation_end('ad3','ad3 stem',kind='splinted',polarity='5-prime',owner='adapter',text='Source5P leg; actual heterogeneous tRNA ends/product unassembled');return s

def source_oligo(b,key,name):
    chem={'source_ordering':SOURCE_ORDERS[name]}
    if name.startswith('ad5_'):chem['sugar_blocks']=[{'start':0,'end':16,'sugar':'DNA'},{'start':16,'end':18,'sugar':'RNA','raw':'rGr'+name[-1]}]
    if name=='ad3':chem.update({'terminal_nucleotide':'ddC additional C','five_prime':'P','three_prime':'dideoxy'})
    if name=='rt':chem.update({'five_prime':'P','non_nucleotide_parts':[{'identity':'iSp18','after_base':len(RT_LEFT),'bases':False}],'NN_purpose':'unassigned, noUMIclaim'})
    return b.g.input(key.lower(),'Source whole modified '+name,K.OLIGO,state=oligo_scene(name),chemical_state=chem)

def protocol_graph():
    b=B('QuantM-tRNAseq selected matureCCA capture and truncatedcDNA recovery',AUTHORITY,'SourceHEK293totalRNA, selecteduntreatedbranch');g=b.g
    g.add_source_item(SourceItem('rna_variant_notation','Main189 writes5prime-adapter variant rGrT; actualData1 named5primeAD ordering has rGrU.','Model uses explicit Data1 ordering and preserves both source notations; actual ordered/used variant unverified, no silentT/U repair.'))
    g.add_source_item(SourceItem('boundary','Main59/Fig1: BOTHtRNAends captured using matureCCA-compatible adapter, PAGE retains truncated+fullcDNAs before circularization.','Combination of capture and postRT product fates, not circularization alone or a tissue panel. Mapping10nt threshold is computational, not physicalRNA sizecut.'))
    b.aliquot('input_qc','ActualNanoDrop/integritygel analytical allocations; amounts/scheduleunassembled',qc=True)
    b.op('deacylate','Source20mMTris9/37C45min at1ugul',inputs=(b.cond('deacyl_conditions','Actualsourcebasicdeacylation'),),product='Source deacylated totalRNA mixture; matureCCA5P/3OH subset and other RNA classes, conversion unassembled')
    b.aliquot('one_ug','Actual1ug deacylatedRNA selected','Actual1ug selected RNA','Untransferred deacylated RNA stock')
    ad5=[source_oligo(b,'ad5_'+x,'ad5_'+x)for x in ('A','G','U','C')];ad3=source_oligo(b,'ad3','ad3')
    duplex=g.product('annealed_adapters','ActualRNA+10pmol3adapter/10pmoltotal4member5adapter mixture; occupiedRNA/complexstatesunassembled',K.MOLECULAR)
    g.react('adapter_anneal','SourceRNA+adapters95C2min/37C15min5mMTris8/.5mMEDTA/10mMMgCl2',S.HYBRIDIZATION,((b.x,R.TEMPLATE),)+tuple((x,R.STRUCTURAL)for x in ad5)+((ad3,R.STRUCTURAL),b.cond('adapter_hyb_conditions','Actual5mMTris8/.5mMEDTA/10mMMgCl2')),(duplex,));b.x=duplex
    rna_component=g.product('heated_rna_component','RNA component of same heated RNA/adapter mixture; occupiedcomplexstateunassembled',K.MOLECULAR)
    adapter_component=g.product('heated_adapter_component','Copresent adapter components of same heated mixture; no newfreeadapteraddition, occupiedstatesunassembled',K.ADAPTER)
    g.react('heated_components','Same-tube heated mixture component inventory, not physicalfractions or fresh additions',S.GENERIC,((duplex,R.TEMPLATE),),(rna_component,adapter_component),changes_molecule=False);b.x=rna_component
    b.op('capture_ligation','Source5UulRnl2/37C60min then4C60min, BOTHtRNAends',S.LIGATION,((adapter_component,R.STRUCTURAL),b.cat('rnl2','ActualsourceRnl2'),b.cond('rnl2_conditions','SourceRNA-ligasebuffer; fullmanufacturerrecipeunassembled')),geometry=AdapterLigationGeometry.known('splinted',adapter_scene(),statement='Known source14bp duplex and two ligation legs; representative rGrA/actual4variantmix, whole tRNA/occupied/sealed products unassembled.'),product='Source tRNA/adapter ligation preparation plusunligated/offtargetclasses; fullcovalentproductsunassembled')
    b.clean('lig_precip','Sourceethanol/two75percentethanol/10ulwater',inputs=(b.stock('glycoblue_lig','ActualGlycoBluecarrier forligationprecipitation'),))
    rt=source_oligo(b,'rt_primer','rt')
    b.op('rt_anneal','Source10pmolRTprimer/.5pmolul70C2min/.1Cs→37C',S.HYBRIDIZATION,((rt,R.STRUCTURAL),b.cond('rt_anneal_conditions','Actualsourceprogramme')))
    rtstock=g.input('rt_inventory','Aggregate sourceSSIV reactioninventory; actualpackaging/NTdoseunassigned',K.MOLECULAR);enzyme=g.product('rt_ssiv','SourceSSIV component',K.CATALYST);nt=g.product('rt_nt','SourceRTnucleotideinventory, formulationunassembled',K.MOLECULAR)
    g.react('rt_inventory_components','Inventory grouping, notclaimedsourcepremix or extra enzymeaddition',S.GENERIC,((rtstock,R.TEMPLATE),),(enzyme,nt),changes_molecule=False)
    b.op('rt','SourceSSIV55C60min, carriedactualRTprimer',inputs=((enzyme,R.CATALYST),(nt,R.STRUCTURAL)),product='SourceRNA/cDNAhybrids includingfull-length and modification/structure-truncatedcDNAs; wholeproductsunassembled')
    b.op('rna_hydrolyze','Source.1NNaOH98C20min RNA hydrolysis',inputs=(b.cond('naoh','ActualNaOHprogramme'),),product='SourcecDNA preparation afterRNAhydrolysis; freeRNAhydrolysisproducts/carryoverunassembled')
    b.clean('post_rt_precip','Sourceethanol/two75percentwashes/12ulwater',inputs=(b.stock('glycoblue_rt','ActualGlycoBluecarrier'),))
    b.op('page','Source7Murea6percentPAGE/SYBRGold15min/UVexcision of alltRNAderivedcDNA bands',inputs=(b.stock('sybr_gold','ActualSYBRGold stain'),b.cond('page_conditions','Actual gel/excision conditions')),product='SourceexcisedfullandtruncatedcDNA region; no forcedfull-lengthselection')
    b.clean('extract_cdna','Sourcecrush400ul300mMNaCl/10mMTris8/1mMEDTA/dryice30min/RTnight/isopropanol/two75percentwashes/12ulwater',inputs=(b.stock('glycoblue_extract','ActualGlycoBluecarrier'),))
    b.op('circle','SourceCircLigase.5Uul60C1h/80C20min',inputs=(b.cat('circligase','ActualCircLigase'),b.cond('circle_conditions','Sourcemanufacturerreactionconditionsunassembled')),product='SourcecDNAcircularizationpreparation; RT5P/iSp18 retained in intendedfull/truncatedclasses, sealedwholeproductsunassembled')
    b.clean('circle_precip','Sourceethanol/two75percentwashes/12.5ulwater',inputs=(b.stock('glycoblue_circle','ActualGlycoBluecarrier'),))
    pf=source_oligo(b,'pcr_forward','pcr_f');pr=source_oligo(b,'pcr_reverse','pcr_r')
    master=g.input('q5_master','ActualsourceNEBNextUltraIIQ5commercialmaster; detailedformulationunassembled',K.MOLECULAR)
    pol=g.product('q5_enzyme','SameactualcommercialmasterQ5polymerase',K.CATALYST);pnt=g.product('q5_nt','SameactualcommercialmasterdNTP inventory; doseunassigned',K.MOLECULAR);pbuf=g.product('q5_buffer','Sameactualcommercialmasterconditions',K.CONDITION)
    g.react('q5_components','Copresent actual commercial-master components, not extra independent additions',S.GENERIC,((master,R.TEMPLATE),),(pol,pnt,pbuf),changes_molecule=False)
    b.op('pcr','SourceNEBNextUltraIIQ5master/7cycles selectedHEK293',inputs=((pf,R.STRUCTURAL),(pr,R.STRUCTURAL),(pol,R.CATALYST),(pnt,R.STRUCTURAL),(pbuf,R.CONDITION)),product='Source indexed amplified full/truncatedcDNA librarypreparation, wholePCRgeometryunassembled')
    b.op('agarose','Source2percentagarose/.05mgmlEtBr/100–250bp excision',inputs=(b.stock('etbr','ActualEtBrstain'),b.cond('agarose_condition','Actualsourcegelcut')),product='Sourceexcised100–250bpPCRlibraryfraction; fullinternalproductsunassembled')
    b.clean('final_extract','SourceQiaquickgelcleanupatRT/alloptionalsteps/ethanol/two75percentwashes/10ulwater',inputs=(b.stock('glycoblue_final','ActualGlycoBluecarrier'),))
    b.aliquot('final_qc','ActualQubit/DNAHSBioanalyzer/qPCRanalyticalallocations, scheduleunassembled',qc=True)
    b.aliquot('pool_aliquot','ActualqPCR-directedmultiplexingselectedcontribution, amountunassigned')
    other=g.input('other_library_aliquots','Actualotherqualifiedlibrarycontributions; stocks/normalizationvolumesunassembled',K.MOLECULAR)
    pool=g.product('library_pool','Actualmultiplexedlibrarypool, proportionsunassigned',K.MOLECULAR);g.react('pool','SourceqPCR-directedmultiplexing, notassumed equimolar',S.POOL,((b.x,R.TEMPLATE),(other,R.TEMPLATE)),(pool,));b.x=pool
    b.aliquot('loading','ActualNextSeqloadingallocation, amountunassigned')
    lib=g.product('library','SelectedQuantM sourcepoolNGSlibrary; wholemolecules/readlandingunassembled',K.LIBRARY);g.react('load','SourceNextSeq550v2.5 single-end110cycles',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(read_lengths={'R1':110},limitation='SourceSE110; actualsampleindexvalues/lengthcycles/primerlandingunassembled. NoRNAphysical10ntselectioninferredfrommapping.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import zipfile,xml.etree.ElementTree as E
    c=Check();g=protocol_graph();adapter_scene()
    c('PhysicalddCnucleotideretained',len(AD3_BODY)+1,31)
    c('RTNNnotidentifier',g.materials['rt_primer'].chemical_state['NN_purpose'],'unassigned, noUMIclaim')
    c('RNAactuallyheatedwithadapters',g.steps['adapter_anneal'].inputs[0].material,'one_ug_selected')
    c('LigationusesheatedRNAmixtureonce',g.steps['heated_components'].inputs[0].material,'annealed_adapters')
    c('Ligationadapterisretainedheatedcomponent',g.steps['capture_ligation'].inputs[1].material,'heated_adapter_component')
    c('ActualRTprimeronlyoneexternalinput',sum(x.id=='rt_primer' for x in g.materials.values()),1)
    c('Circligaseaftertruncatedproductselection',g.steps['circle'].inputs[0].material,'extract_cdna_kept')
    c('Multiplexinghasotherlibrarycontributor',g.steps['pool'].inputs[1].material,'other_library_aliquots')
    c('Loadingstockremainderpreserved','loading_remainder'in g.terminals,True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41467-020-17879-x_supp_41467_2020_17879_MOESM5_ESM.xlsx'
    if have(c,Source(p,'https://doi.org/10.1038/s41467-020-17879-x')):
        z=zipfile.ZipFile(p);ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};root=E.fromstring(z.read('xl/sharedStrings.xml'));strings=[''.join(x.itertext())for x in root.findall('s:si',ns)];actual={''.join(x.split())for x in strings}
        for k,v in SOURCE_ORDERS.items():c('Wholemodifiedsourceorder '+k,''.join(v.split())in actual,True)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
