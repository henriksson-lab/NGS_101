"""Source-owned selected BOLORAMIS component within actual multiplex preparation."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from dataclasses import dataclass
from chemdraw import Scene,Segment,Construct,revcomp,feature
from circular import circularize_ssdna,rolling_circle,CircularizationGeometry
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition,SourceCoverage,MolecularRelation as MR
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


# Curator-reviewed fixed inventories follow source main181/237/239, S3 and S8; never graph reflection.
PREPARATION_SCOPE = ('source', 'fix', 'pfa', 'fix_conditions', 'pbs', 'pbs_conditions', 'permeabilize', 'permeabilization', 'exact_padlock', 'rca_primer', 'mirna_rca_primer', 'other_probes', 'hybridize', 'hyb_buffer', 'hyb_wash', 'hyb_wash_conditions', 'bound_inventory', 'circle', 'splintr', 'splintr_buffer', 'postlig_wash', 'postlig_wash_conditions', 'rca_mix', 'rca_mix_parts', 'rca', 'rca_wash', 'rca_wash_conditions', 'mes', 'mes_conditions', 'anchor_fix', 'bspeg9', 'initial_fix_mix', 'borate', 'borate_stock', 'quench', 'tris', 'fixed_pbs', 'fixed_pbs_conditions', 'readout', 'library')
SEQUENCING_SCOPE = ('cycle1_anchor', 'cycle1plus', 'cycle1minus', 'cycle1_mirna_plus', 'cycle1_mirna_minus', 'cycle1_hyb', 'cycle1_wash', 'cycle1_wash_conditions', 'cycle1_ligation_supplies', 'cycle1_mix_parts', 'cycle1_ligate', 'cycle1_postwash', 'cycle1_postwash_conditions', 'cycle1_image', 'cycle1_strip', 'cycle1_strip_buffer', 'cycle1_rinse', 'cycle1_rinse_conditions', 'cycle2_anchor', 'cycle2plus', 'cycle2minus', 'cycle2_mirna_plus', 'cycle2_mirna_minus', 'cycle2_hyb', 'cycle2_wash', 'cycle2_wash_conditions', 'cycle2_ligation_supplies', 'cycle2_mix_parts', 'cycle2_ligate', 'cycle2_postwash', 'cycle2_postwash_conditions', 'cycle2_image', 'cycle2_strip', 'cycle2_strip_buffer', 'cycle2_rinse', 'cycle2_rinse_conditions', 'cycle3_anchor', 'cycle3plus', 'cycle3minus', 'cycle3_mirna_plus', 'cycle3_mirna_minus', 'cycle3_hyb', 'cycle3_wash', 'cycle3_wash_conditions', 'cycle3_ligation_supplies', 'cycle3_mix_parts', 'cycle3_ligate', 'cycle3_postwash', 'cycle3_postwash_conditions', 'cycle3_image', 'cycle3_strip', 'cycle3_strip_buffer', 'cycle3_rinse', 'cycle3_rinse_conditions', 'cycle4_anchor', 'cycle4plus', 'cycle4minus', 'cycle4_mirna_plus', 'cycle4_mirna_minus', 'cycle4_hyb', 'cycle4_wash', 'cycle4_wash_conditions', 'cycle4_ligation_supplies', 'cycle4_mix_parts', 'cycle4_ligate', 'cycle4_postwash', 'cycle4_postwash_conditions', 'cycle4_image', 'cycle4_strip', 'cycle4_strip_buffer', 'cycle4_rinse', 'cycle4_rinse_conditions', 'cycle5_anchor', 'cycle5plus', 'cycle5minus', 'cycle5_mirna_plus', 'cycle5_mirna_minus', 'cycle5_hyb', 'cycle5_wash', 'cycle5_wash_conditions', 'cycle5_ligation_supplies', 'cycle5_mix_parts', 'cycle5_ligate', 'cycle5_postwash', 'cycle5_postwash_conditions', 'cycle5_image', 'cycle5_strip', 'cycle5_strip_buffer', 'cycle5_rinse', 'cycle5_rinse_conditions', 'cycle6_anchor', 'cycle6plus', 'cycle6minus', 'cycle6_mirna_plus', 'cycle6_mirna_minus', 'cycle6_hyb', 'cycle6_wash', 'cycle6_wash_conditions', 'cycle6_ligation_supplies', 'cycle6_mix_parts', 'cycle6_ligate', 'cycle6_postwash', 'cycle6_postwash_conditions', 'cycle6_image', 'cycle6_strip', 'cycle6_strip_buffer', 'cycle6_rinse', 'cycle6_rinse_conditions', 'cycle7_anchor', 'cycle7plus', 'cycle7minus', 'cycle7_mirna_plus', 'cycle7_mirna_minus', 'cycle7_hyb', 'cycle7_wash', 'cycle7_wash_conditions', 'cycle7_ligation_supplies', 'cycle7_mix_parts', 'cycle7_ligate', 'cycle7_postwash', 'cycle7_postwash_conditions', 'cycle7_image', 'cycle7_strip', 'cycle7_strip_buffer', 'cycle7_rinse', 'cycle7_rinse_conditions', 'cycle8_anchor', 'cycle8plus', 'cycle8minus', 'cycle8_mirna_plus', 'cycle8_mirna_minus', 'cycle8_hyb', 'cycle8_wash', 'cycle8_wash_conditions', 'cycle8_ligation_supplies', 'cycle8_mix_parts', 'cycle8_ligate', 'cycle8_postwash', 'cycle8_postwash_conditions', 'cycle8_image')
def reviewed_inventory(g):
    g.declare_reviewed_scope(PREPARATION_SCOPE,'PublishedPMC8191787main181/227/237; originalSupplementaryTableS3C3/C4/C6–C9/C12 andactualFigS8p12. Actualmultiplexpoolretained; selectedACTN1perfectRNAalignmentderived/inferred, notgenomic-verified. Wholeotherpoolmembersavailable/unexpanded.')
    g.declare_reviewed_scope(SEQUENCING_SCOPE,'Publishedmain239 andactualS8p12: eightactualanchor/nonamerligation/imagecycles; sevenbetween-roundstrip/rinses. BothmRNA/miRNAanchorssourceS3; fullnonamerpatterns/cycle-offsetnotassembled. No implicitkitinputs.')
    g.cover_source_item(SourceCoverage('scope','modelled','library'))

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
    g.react('bound_inventory','INFERRED — RNA-match view within actual copresent in-situ components; not physical fractionation',S.GENERIC,((b.x,R.TEMPLATE),),(ctx,retained,prim,others,proteins),changes_molecule=False)
    circ=g.product('selected_circle','SelectedsuccessfulACTN1padlockcircleclass; notallinputprobesconverted',K.MOLECULAR,chemical_state=circle())
    off=g.product('ligation_offproducts','Remaining/unligated/off-targetsourceclasses; amountunassigned',K.MOLECULAR)
    g.react('circle','INFERRED — RNA-splint match example; source SplintR210nM1:50/1hRT directRNA-splintedDNAclosure(noRT)',S.CIRCULARIZATION,((ctx,R.TEMPLATE),(retained,R.TEMPLATE),b.stock('splintr','ActualSplintRligase210nM',K.CATALYST,R.CATALYST),b.cond('splintr_buffer','Actual1xSplintRbuffer')),(circ,off),geometry=CircularizationGeometry.known('splinted single-strand ligation',statement='Source5P60ntpadlock closesonRNA. Localperfectmatch25ntRNAisderived/inferred; enzymeyield/specificity andactualoccupiedRNAunassembled.'),relations=(MR.CIRCLE_CLOSURE,))
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
    reviewed_inventory(g)
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
    c('Source-mappedpublicationinventory',g.publication_problems(),[])
    g.input('unreviewed_added_stock','Regression-only unreviewed input',K.CONDITION)
    c('Unreviewedadditioncannotautopublish',any('unreviewed_added_stock' in p and 'source-inventory' in p for p in g.publication_problems()),True);c.report()

PAGE_FILE = 'boloramis.html'
def render_page():
    from page import head,info
    from chemdraw import panel
    from protocol_views import assemble_page,WorkflowSection,OligoSection,OligoEntry,PlanLimitationSection
    g=protocol_graph().require_publishable()
    entries=[OligoEntry('ACTN1 A1 padlock',tuple(probe()),mods='5′ phosphate'),OligoEntry('mRNA RCA primer',(Segment('RCA primer',RCA_PRIMER),)),OligoEntry('mRNA sequencing anchor +',(Segment('anchor',ANCHOR),)),OligoEntry('mRNA sequencing anchor −',(Segment('anchor',ANCHOR),),mods='5′ phosphate'),OligoEntry('miRNA RCA primer',(Segment('RCA primer',MIRNA_RCA),)),OligoEntry('miRNA sequencing anchor +',(Segment('anchor',MIRNA_ANCHOR),)),OligoEntry('miRNA sequencing anchor −',(Segment('anchor',MIRNA_ANCHOR),),mods='5′ phosphate')]
    return assemble_page(g,[head('BOLORAMIS — RNA-splinted padlocks'),'<div class="wrap" data-protocol-complete="true"><h1>BOLORAMIS</h1>',
      '<p>Identify RNA targets directly in fixed cells using barcoded DNA padlocks. SplintR closes probes on RNA without reverse transcription; Phi29 rolling-circle amplification produces localized copies for eight rounds of nonamer sequencing by ligation.</p>',
      info('Sources: <a href="https://doi.org/10.1093/nar/gkab120">Liu et al., NAR 2021</a>, published main, Supplementary Tables and Figure S8. The actual multiplex experiment contains 957 mRNA probes and 77 miRNA probes. One ACTN1 probe is shown in detail; the other real pool members remain copresent and unexpanded in this model.'),
      '<h2>Probe and barcode</h2><p>The selected 60-nt probe has an 18-nt arm, a copied 4-nt target code, a 23-nt anchor, the complete 8-nt target code, and a 7-nt arm. The copied code repeats the last four bases of the eight-base code. These codes identify a target gene, not a cell or a unique RNA molecule.</p>',OligoSection('Source oligos',tuple(entries)),
      '<h2>RNA-splinted closure</h2>',info('INFERRED — the illustrated 25-nt RNA is calculated from perfect complementarity to the source probe arms. Its transcript sequence and actual occupied complex have not been independently verified. Inferred RNA bases are marked below. The two displayed DNA end tracts belong to one covalent padlock whose backbone continues off-panel.'),
      panel(ConditionalRNAMatch().scene().rows(),caption='INFERRED — conditional perfect-match RNA splint; source probe 5′ phosphate and 3′ hydroxyl meet at the closure nick.'),
      '<p>The circle and RCA backbone are derived from the exact padlock and primer. This checks strand order and primer landing; it does not predict ligation specificity or conversion yield. Three displayed repeats are an illustration, not an experimental copy count. Aminoallyl-dU positions and occupancy remain unresolved.</p>',panel(rca_scene().rows(),caption='Derived RCA backbone example; target-code meanings are retained in both complementary copies.'),
      '<h2>Preparation and in-situ reads</h2><p>Fix and permeabilize cells, hybridize the multiplex padlocks with their RCA primers, and wash before SplintR ligation. The retained primer supplies Phi29 RCA with dNTPs, aminoallyl-dUTP and BSA. BS(PEG)9 then links aminoallyl groups to surrounding proteins; borate is added halfway through the 30-minute fixation. Fresh sequencing anchors and fluorescent nonamers are added for each imaging round, with stripping between rounds.</p>',
      info('Known local anchor pairing, phosphate/hydroxyl polarities and dye ends are displayed in the full graph. Whole nonamer patterns, cycle offsets, exact reported-base mapping and complete sealed junctions remain unresolved. MiRNA reads were excluded during downstream analysis; this does not remove those reagents from the physical preparation.'),PlanLimitationSection('library'),
      '<details><summary>Full preparation, all eight sequencing rounds and chemistry warnings</summary>',WorkflowSection('Source-qualified multiplex preparation'),'</details>',
      '<p><a href="01_boloramis.html">Source notes and remaining questions</a></p></div>'])

if __name__=='__main__':selftest()
