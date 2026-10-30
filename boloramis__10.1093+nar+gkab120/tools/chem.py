"""Source-owned selected BOLORAMIS component within actual multiplex preparation."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from dataclasses import dataclass
from chemdraw import Scene,Segment,Construct,revcomp,feature
from circular import circularize_ssdna,rolling_circle,CircularizationGeometry
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition,MolecularRelation as MR
AUTHORITY='Published2021PMC8191787 main181/227/237/239, originalSupplementaryTablesFinal workbookS1/S3 and SI FigS8actualPDFp12. SelectedACTN1A1 component withinactual96mRNA/77miRNA multiplexexperiment; otherpoolmembers copresent/unexpanded, not a fabricated singleplex.'
SOURCE_PROBE='/5Phos/AGGTCTGGACAATCTTGGcttcATTAAATCGGTGGCGCGGGCTCTgcaacttcACGTGGT'
ANCHOR='ATTAAATCGGTGGCGCGGGCTCT'
RCA_PRIMER='AGAGCCCGCGCCACCGATTTAAT'
MIRNA_RCA='TCAAGGCGCGGGCAG'
MIRNA_ANCHOR='ATAGCGATCTGCCCGCGCCTTGA'
LONG='AGGTCTGGACAATCTTGG';SHORT='ACGTGGT';CODE='GCAACTTC'
def probe():
    return Construct([Segment('5prime longRNA-bindingarm',LONG),Segment('4base copiedtargetcode',CODE[-4:],feature=feature('boloramis.actn1.copy','feature_barcode','fixed',group='boloramis.actn1',part='copied_5_to_8',note='Duplicate last4of8 target-code; not additional molecule ID')),Segment('universal23anchor',ANCHOR),Segment('8base targetcode',CODE,feature=feature('boloramis.actn1.code','feature_barcode','fixed',group='boloramis.actn1',part='whole_targetcode',note='ACTN1targetidentity sharedbyprobes, notUMI/cellbarcode')),Segment('3prime shortRNA-bindingarm',SHORT)],'SourceACTN1_1932-1957_A1 whole60nt5P DNApadlock')
def padlock_scene():
    sc=Scene();sc.strand('padlock',list(probe()),mod5='P',mod3='OH');sc.ligation_end('padlock','5prime longRNA-bindingarm',kind='single-strand',polarity='5-prime',owner='substrate');sc.ligation_end('padlock','3prime shortRNA-bindingarm',kind='single-strand',polarity='3-prime',owner='substrate');return sc

def strand(seq,name,mod5='',mod3=''):
    sc=Scene();sc.strand(name,[Segment(name,seq)],mod5=mod5,mod3=mod3);return sc
@dataclass(frozen=True)
class ConditionalRNAMatch:
    inferred:bool=True
    @property
    def sequence(self):return revcomp(SHORT+LONG).replace('T','U')
    def scene(self):
        # Two displayed DNA end tracts belong to ONE covalent padlock with off-panel backbone.
        sc=Scene();sc.strand('RNA',[Segment('RNA opposite longarm',revcomp(LONG).replace('T','U'),inferred=self.inferred),Segment('RNA opposite shortarm',revcomp(SHORT).replace('T','U'),inferred=self.inferred)])
        sc.anneal('probe5arm',[Segment('longarm',LONG)],to='RNA',pair=('longarm','RNA opposite longarm'),mod5='P',label='ONE probe5arm; backbone continues off-panel')
        sc.anneal('probe3arm',[Segment('shortarm',SHORT)],to='RNA',pair=('shortarm','RNA opposite shortarm'),mod3='OH',label='SAME probe3arm; closure nick next to5P')
        sc.ligation_end('probe5arm','longarm',kind='splinted',polarity='5-prime',owner='substrate');sc.ligation_end('probe3arm','shortarm',kind='splinted',polarity='3-prime',owner='substrate');return sc

def circle():return circularize_ssdna(probe(),five_prime_phosphate=True,three_prime_oh=True)
def rca():return rolling_circle(circle(),RCA_PRIMER,copies=3)
def rca_scene():
    product=rca();unit=product.unit;needle=revcomp(CODE);assert unit.count(needle)==1
    start=unit.index(needle);copy_start=len(RCA_PRIMER);copy_end=copy_start+len(CODE[-4:])
    assert unit[copy_start:copy_end]==revcomp(CODE[-4:])
    pieces=[]
    for i in range(product.copies):
        pieces.extend([Segment('repeat'+str(i)+' primer',unit[:copy_start]),Segment('repeat'+str(i)+' copied4code',unit[copy_start:copy_end],feature=feature('boloramis.actn1.copy','feature_barcode','fixed',group='boloramis.actn1',part='copied_5_to_8',note='Copiedcomplementofsource4ntredundanttargetcode, notUMI')),Segment('repeat'+str(i)+' arms',unit[copy_end:start]),Segment('repeat'+str(i)+' targetbarcode',needle,feature=feature('boloramis.actn1.code','feature_barcode','fixed',group='boloramis.actn1',part='whole_targetcode',note='Complementarycopy ofsourceACTN1code, notnewUMI')),Segment('repeat'+str(i)+' right',unit[start+len(needle):])])
    sc=Scene();sc.strand('RCA',[x for x in pieces if x.top],label='DerivedDNA backbone repeat example; aminoallyl-dU positions/occupancyunassembled');return sc

def sequencing_anchor_scene(minus=False):
    product=rca();u=product.unit;site=revcomp(ANCHOR);assert u.count(site)==1
    i=u.index(site);sc=Scene();sc.strand('RCA',[Segment('beforeanchor',u[:i]),Segment('anchorbinding',site),Segment('afteranchor',u[i+len(site):])]);sc.anneal('anchor',[Segment('source23ntanchor',ANCHOR)],to='RCA',pair=('source23ntanchor','anchorbinding'),mod5='P'if minus else'')
    sc.ligation_end('anchor','source23ntanchor',kind='splinted',polarity='5-prime'if minus else'3-prime',owner='substrate');return sc

def nonamer_scene(minus=False):
    sc=Scene();sc.strand('nonamer',[Segment('source9nt sequencing reagent','[whole9ntorderingunassembled]',placeholder=True,length_bp=9)],mod5='sourcefluorophore'if minus else'P',mod3='OH'if minus else'sourcefluorophore');sc.ligation_end('nonamer','source9nt sequencing reagent',kind='splinted',polarity='3-prime'if minus else'5-prime',owner='adapter');return sc

class B:
    def __init__(self):
        self.g=ProtocolGraph('BOLORAMIS actualmultiplex preparation, selectedACTN1component',authority=AUTHORITY,representation='published_instance');self.x=self.g.input('source','SourceiNGN/HMC3co-culture withtargetRNA; cells/otherRNAs retained, notpureACTN1',K.BIOLOGICAL,primary=True)
    def stock(self,k,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):return(self.g.input(k,label,kind,**kw),role)
    def cond(self,k,label):return self.stock(k,label,K.CONDITION,R.CONDITION)
    def op(self,k,label,kind=S.MAGIC,ins=(),product='Actualin-situ sample; occupiedcomplexes/productstoichiometry/conversionunassembled',**kw):
        y=self.g.product(k+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Actualsourceoperation; completeoccupiedmolecularstates/recoveryunassembled, notsourceabsence.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        self.g.react(k,label,kind,((self.x,R.TEMPLATE),)+tuple(ins),(y,),**kw);self.x=y;return y
    def wash(self,k,label,state=None,chemical=None,extra=()):
        y=self.g.product(k+'_kept','Actualretainedin-situ preparation; quantitativecarryoverunassigned',K.MOLECULAR,state=state,chemical_state=chemical);w=self.g.product(k+'_removed','Actualremovedwashfraction',K.WASTE)
        self.g.react(k,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(k+'_conditions',label))+tuple(extra),(y,w));self.g.discard(w,'Sourcewashfractionremoved; no perfect molecular selectivity assertion.');self.x=y

def protocol_graph():
    b=B();g=b.g
    b.op('fix','4percentPFA15min',S.CROSS_LINKING,(b.stock('pfa','Actualparaformaldehydecrosslinker'),b.cond('fix_conditions','Source4percentPFA15min')))
    b.wash('pbs','SourcePBSrinse')
    b.op('permeabilize','Source.25percentTritonX10015min',S.PERMEABILIZATION,(b.cond('permeabilization','ActualsourceTritonprogramme'),))
    exact=g.input('exact_padlock','ActualACTN1A1 component ofsourceprobe pool,5phosphorylated',K.OLIGO,state=padlock_scene(),chemical_state=probe())
    primer=g.input('rca_primer','ActualmRNARCAprimerT2S3*; sequence23nt, asteriskpartofname notsourcePSordering',K.OLIGO,state=strand(RCA_PRIMER,'RCAprimer'))
    mirca=g.input('mirna_rca_primer','Actualsourceco-pooledmiRNARCAprimerR23.3short,15nt; bodyknown, physicalsupplypackagingunassigned',K.OLIGO,state=strand(MIRNA_RCA,'miRNA RCAprimer'))
    other=g.input('other_probes','Otheractualmultiplexpoolmembers:956mRNA+77miRNAprobes, sourcewholeordersavailablebutoutsidechosencomponentcuration; notabsent',K.OLIGO)
    b.op('hybridize','Actualmultiplex10uMpadlockpool+RCAprimers1:1moleratio,1h37C',S.HYBRIDIZATION,((exact,R.STRUCTURAL),(primer,R.STRUCTURAL),(mirca,R.STRUCTURAL),(other,R.STRUCTURAL),b.cond('hyb_buffer','10percentformamide/6xSSC')))
    b.wash('hyb_wash','Three5minhybridizationbuffer washes then1xSplintRbuffer; excessunhybridizedprobe/primerremoved')
    ctx=g.product('rna_context','Samehybridizedcell/RNAcontext; selectedperfectRNA-matchclassinferred, completeactualtarget/occupiedmoleculesunassembled',K.MOLECULAR,state=ConditionalRNAMatch().scene())
    retained=g.product('retained_padlock','SameACTN1A1boundlinearprobe component afterwashing',K.OLIGO,state=padlock_scene(),chemical_state=probe())
    prim=g.product('retained_primer','SamepriorhybridizationRCAprimer component, notfresholigoaddition',K.OLIGO,state=strand(RCA_PRIMER,'RCAprimer'))
    proteins=g.product('retained_proteins','Samein-situ surroundingproteincomponent retainedfromoriginalcells for sourceAA-dU/BSPEG9anchoring, notafreshproteinaddition',K.MOLECULAR)
    others=g.product('retained_other_probes','Sameotherretainedpanelprobeclasses pluscopresentmiRNARCAprimer; unexpandedmolecularassemblies',K.MOLECULAR)
    g.react('bound_inventory','Copresentin-situ components, notphysicalfractionation',S.GENERIC,((b.x,R.TEMPLATE),),(ctx,retained,prim,others,proteins),changes_molecule=False)
    circ=g.product('selected_circle','SelectedsuccessfulACTN1padlockcircleclass; notallinputprobesconverted',K.MOLECULAR,chemical_state=circle())
    off=g.product('ligation_offproducts','Remaining/unligated/off-targetsourceclasses; amountunassigned',K.MOLECULAR)
    g.react('circle','SplintR210nM1:50in1xbuffer,1hRT; directRNA-splintedDNAclosure(noRT)',S.CIRCULARIZATION,((ctx,R.TEMPLATE),(retained,R.TEMPLATE),b.stock('splintr','ActualSplintRligase210nM',K.CATALYST,R.CATALYST),b.cond('splintr_buffer','Actual1xSplintRbuffer')),(circ,off),geometry=CircularizationGeometry.known('splinted single-strand ligation',statement='Source5P60ntpadlock closesonRNA. Localperfectmatch25ntRNAisderived/inferred; enzymeyield/specificity andactualoccupiedRNAunassembled.'),relations=(MR.CIRCLE_CLOSURE,))
    g.mark_terminal(Terminal(off.id,Disposition.UNRESOLVED,'Actualunligated/offtargetclasses outside selectedsuccessfulcomponent; subsequentwashremovalnotassumedquantitative.'));b.x=circ
    b.wash('postlig_wash','Threehybridizationbuffer rinses removeunligatedproductsandenzyme',chemical=circle())
    mix=g.input('rca_mix','ONEactualsourceRCAmix in1xPhi29buffer',K.MOLECULAR)
    specs=[('phi29',K.CATALYST,'Source1UulNxGenPhi29',R.CATALYST),('dntp',K.MOLECULAR,'Source250uMdNTPcomponents',R.STRUCTURAL),('aadutp',K.MOLECULAR,'Source40uMaminoallyldUTPchemicalcontributor',R.STRUCTURAL),('bsa',K.MOLECULAR,'Source200ugmlBSAcomponent',R.STRUCTURAL),('buffer',K.CONDITION,'Source1xPhi29buffer',R.CONDITION)]
    ms=[g.product('rca_'+a,l,k)for a,k,l,r in specs];g.react('rca_mix_parts','Copresentcomponents ofoneactualRCAmix, notfreshextras',S.GENERIC,((mix,R.TEMPLATE),),tuple(ms),changes_molecule=False)
    product=g.product('rca_product','SelectedRCAbackbonecopyexample; sourceaminoallyldUoccupancy/completecellcontextunassembled',K.MOLECULAR,state=rca_scene(),chemical_state=rca())
    g.react('rca','SourcePhi29RCA2h37C',S.EXTENSION,((b.x,R.TEMPLATE),(prim,R.STRUCTURAL),(others,R.STRUCTURAL))+tuple((m,spec[3])for m,spec in zip(ms,specs)),(product,),relations=(MR.ROLLING_CIRCLE_COPY,));b.x=product
    b.wash('rca_wash','Three5min2xSSCwashes',state=rca_scene())
    b.wash('mes','Gentlesingle10mMMESpH6.5rinse')
    b.op('anchor_fix','BeginactualBSPEG9ampliconfixation: first15minRT, aminoallyldU-to-retainedsurroundingproteins',ins=((proteins,R.STRUCTURAL),b.stock('bspeg9','Actual4ulBSPEG9ininitial172ulwater+4ul500mMMES; sourceBSPEG9stockconcentrationunassigned'),b.cond('initial_fix_mix','Initialsource180ulbeforehalfwayborateaddition, notcomplete200ulatstart')))
    b.op('borate','At15minactual20ul1Mborateaddition; continueanother15mintotal30min',ins=(b.stock('borate_stock','Actual20ul1Mborate atreactionhalfway'),))
    b.op('quench','Source1MTrisHCl8quench',S.QUENCHING,(b.cond('tris','Actual1MTris8quench'),))
    b.wash('fixed_pbs','Three1xPBSrinses')
    # Eight actual chemistry/imaging rounds, not eight consecutive extensions of one anchor.
    for j in range(8):
        k='cycle'+str(j+1)
        anchors=[b.stock(k+suffix,'Actualsource23ntsequencinganchor '+suffix,K.OLIGO,R.STRUCTURAL,state=strand(ANCHOR,'anchor',mod5='P'if suffix=='minus'else''))for suffix in ('plus','minus')]
        anchors += [b.stock(k+'_mirna_'+suffix,'Actualco-pooledmiRNAsequencinganchor '+suffix,K.OLIGO,R.STRUCTURAL,state=strand(MIRNA_ANCHOR,'miRNAanchor',mod5='P'if suffix=='minus'else''))for suffix in ('plus','minus')]
        b.op(k+'_anchor','Sequencinganchormix1:100from100uMstocks,15minRT',S.HYBRIDIZATION,tuple(anchors)+(b.cond(k+'_hyb','6xSSC/10percentformamide'),))
        b.wash(k+'_wash','Three1min6xSSCwashes,one1xT4buffer rinse')
        supplies=g.input(k+'_ligation_supplies','Actual200ulnonamermix:168ulwater/10ulT4ligase/2ul100uMnonamerstock/20ul10xbuffer',K.MOLECULAR)
        lig=g.product(k+'_ligase','T4DNAligaseofsame200ulmix',K.CATALYST);non=g.product(k+'_nonamers','Sameactualplus/minusfour-dye9ntreagentinventory; wholeorders/phaseunassembled, sourceS8knownendsretained',K.OLIGO,state=nonamer_scene(),chemical_state={'plus':'5P/3dye','minus':'5dye/3OH','dyes':['Cy5','TexasRed','Cy3','FAM'],'whole_nine_nucleotide_patterns':'unassembled; source21ordinaryunread'})
        buffer=g.product(k+'_buffer','SameT4ligasebuffer/water',K.CONDITION);g.react(k+'_mix_parts','Sameactualnonamermixcomponents, noextraadditions',S.GENERIC,((supplies,R.TEMPLATE),),(lig,non,buffer),changes_molecule=False)
        b.op(k+'_ligate','SourceT4nonamerSBL45–60minRT; sourceknownendpolarity/anchorbinding, wholeoccupiedjunctionunassembled',S.LIGATION,ins=((lig,R.CATALYST),(non,R.STRUCTURAL),(buffer,R.CONDITION)),geometry=AdapterLigationGeometry.known('splinted',sequencing_anchor_scene(),nonamer_scene(),sequencing_anchor_scene(True),nonamer_scene(True),statement='ActualFigS8 ±direction knownanchor/nonamer5P versus3OH interfaces and dyeends; nonamerwholepatterns/cycleoffset/completeoccupiedjunctionandsealedproductsunassembled.'))
        b.wash(k+'_postwash','Three1min6xSSCwashes')
        b.op(k+'_image','Actualin-situ fluorescenceimage, no physical aliquot removed',S.GENERIC,product='SAMEfixedsample afteropticalreadout; notconsumedDNA/QCstock')
        if j<7:
            b.op(k+'_strip','Source80percentformamide/.0001percentTriton80C,three3minstrips',S.THERMAL_TREATMENT,(b.cond(k+'_strip_buffer','Source800ulformamide+10ul.01percentTriton+190ulwater; washprotocol80C'),))
            b.wash(k+'_rinse','Three6xSSCrinsesbeforenextfreshanchorround')
    lib=g.product('library','ActualfixedBOLORAMISsample readbyeightin-situSBLrounds; completefluorescentligatedstatesunassembled',K.LIBRARY)
    g.react('readout','Sourceeightcycleopticalsequencingendpoint',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(read_lengths={'in_situ_barcode':8},limitation='Source8roundsnonamerSBL ±direction/figureS8. Exactcycle-anchoroffset/dye-to-template-to-reportedbase mapping/readmasknotinferredfromrawprobeDNA; codetargetidentitynotUMI. In-situphysicalsample, notIllumina loadinglibrary.'))
    g.add_source_item(SourceItem('scope',AUTHORITY,'Source-knownmoleculeinterfaces/circlebackbone/RCAcopycheckedlocally; targetRNAperfectmatchinferred, wholesample/aa-dUplacement/BSPEG9adductstoichiometry/sequencingjunctions remainunassembled.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import zipfile,xml.etree.ElementTree as E
    c=Check();g=protocol_graph()
    c('DerivedconditionalRNAdefaultsuncertain',ConditionalRNAMatch().inferred,True)
    c('SourcepadlockexactwholeDNA',probe().top(),SOURCE_PROBE.removeprefix('/5Phos/').upper())
    c('Copiedfourtargetcodebasesretained',CODE[-4:],'CTTC')
    c('RNA-splintedlocalarmsconstructwithoutpairingerror',len(ConditionalRNAMatch().scene().strands),3)
    c('ClosedcirclepreservesonepadlocknotcopiedRNA',circle().sequence,probe().top())
    c('RCAprimerlandsandrepeatlengthfollowscircle',len(rca().unit),len(circle()))
    c('SourceRCAprimerretainedfromhybridization',next(i.material for i in g.steps['rca'].inputs if i.material=='retained_primer'),'retained_primer')
    c('NoextraRCAprimerstock',sum(m.id=='rca_primer' for m in g.materials.values()),1)
    c('Otheractualpoolmemberscopresent',any(i.material=='other_probes' for i in g.steps['hybridize'].inputs),True)
    c('AA-dUTPactualmolecularcomponent',g.materials['rca_aadutp'].kind,K.MOLECULAR)
    c('SBLusesknownpartialinterfaces',g.steps['cycle1_ligate'].geometry.kind,'splinted')
    c('All8sourcechemistryrounds',sum(s.id.endswith('_ligate') for s in g.steps.values()),8)
    c('Strippingbetweenroundsnotbeforefirst',g.steps['cycle2_anchor'].inputs[0].material,'cycle1_rinse_kept')
    c('Noanalyticalaliquotinventedforopticalreadout',any(t.disposition==Disposition.QC for t in g.terminals.values()),False)
    c('Sourceplussubstratehas3OHend',sequencing_anchor_scene().ligation_ends[0].polarity,'3-prime')
    c('Sourceminusanchorhas5Pdonor',sequencing_anchor_scene(True).strands['anchor'].mod5,'P')
    p=sources_root()/Path(__file__).resolve().parents[1].name/'nar_gkab120_supp_gkab120_supplemental_files_Supplementary_Tables_Final.xlsx'
    if have(c,Source(p,'https://doi.org/10.1093/nar/gkab120')):
        ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
        with zipfile.ZipFile(p)as z:
            ss=[''.join(x.itertext())for x in E.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',ns)];rows=E.fromstring(z.read('xl/worksheets/sheet3.xml')).findall('m:sheetData/m:row',ns);cells={cell.attrib['r']:ss[int(cell.find('m:v',ns).text)]for r in rows for cell in r.findall('m:c',ns)if cell.attrib.get('t')=='s'}
            for label,got,cell in [('WholepadlockoriginalA1',SOURCE_PROBE,'C12'),('WholeRCAprimeroriginal',RCA_PRIMER,'C3'),('Wholeplusanchororiginal',ANCHOR,'C6'),('Wholeminusanchor5Poriginal','/5Phos/'+ANCHOR,'C7'),('WholemiRNARCAoriginal',MIRNA_RCA,'C4'),('WholemiRNAplusanchororiginal',MIRNA_ANCHOR,'C8'),('WholemiRNAminusanchororiginal','/5Phos/'+MIRNA_ANCHOR,'C9')]:c(label,got,cells[cell])
    c('Endpointretainsprimarybiologicalprovenance',any('lacks primary-template provenance' in x for x in g.publication_problems()),False)
    c('Privateincompletegraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
