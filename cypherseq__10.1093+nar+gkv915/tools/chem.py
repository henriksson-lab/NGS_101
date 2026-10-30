from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from dataclasses import dataclass,replace
from chemdraw import Scene,Segment,Construct,feature,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from page import head,info
PAGE_FILE='chem.html'

class Operations:
    def __init__(self,g):self.g=g
    def op(self,key,x,label,reagent=None,state=None,chemical_state=None):
        g=self.g;inputs=[(x,R.TEMPLATE)]
        if reagent:
            y=g.input(key+'_condition',reagent,K.CONDITION);inputs.append((y,R.CONDITION))
        y=g.product(key+'_out',label+'; preparation includes unremoved companion species',K.MOLECULAR,state=state,chemical_state=chemical_state)
        g.react(key,label,S.GENERIC,inputs,(y,));return y
    def join(self,key,inputs,label,state=None,chemical_state=None):
        g=self.g;y=g.product(key+'_out',label+'; reaction components remain until an explicit separation',K.MOLECULAR,state=state,chemical_state=chemical_state)
        g.react(key,label,S.GENERIC,tuple((x,R.TEMPLATE if i==0 else R.STRUCTURAL) for i,x in enumerate(inputs)),(y,));return y
    def clean(self,key,x,label):
        g=self.g;z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Source physically retained fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Source physically removed fractions; exact composition/purity unassigned',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Only the source physical separation; no perfect selectivity asserted.');return y
    def split(self,key,x,label,other_label,other=Disposition.UNRESOLVED):
        g=self.g;y=g.product(key+'_production',label,K.MOLECULAR)
        z=g.product(key+'_other',other_label,K.MOLECULAR)
        g.react(key,'Source physical allocation: '+label,S.SPLIT,((x,R.TEMPLATE),),(y,z));g.terminals[z.id]=Terminal(z.id,other,other_label);return y,z
    def qc(self,key,x,label):
        g=self.g;y=g.product(key+'_production','Production preparation retained independently from analytical aliquot',K.MOLECULAR)
        z=g.product(key+'_assay',label,K.QC)
        g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,z));g.terminals[z.id]=Terminal(z.id,Disposition.QC,'Analytical aliquot does not become the sequencing substrate.');return y
    def opaque(self,key,x,label):
        g=self.g;y=g.product(key+'_out',label,K.MOLECULAR)
        g.react(key,label,S.MAGIC,((x,R.TEMPLATE),),(y,),reason='Selected originals establish this operational continuation; full molecular inventory/geometry remains unassembled.');return y

def finish(g,x,readout):
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=readout))
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unassembled source inputs, molecular geometry and directly linked targets',LIMIT));return g

def render_page():
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_chem.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(protocol_graph(),ENDPOINT))+'</main>'


from checks import Source,have
from storage import sources_root
import illumina as il
TITLE='CypherSeq — barcoded circular vector and targeted RCA route'
SOURCE='PublishedPMC4756847 Methods57–100/Fig1; publishedFile012SIworkbookTablesSI1/4/6/8,File007captionsS1–S2'
LIMIT='Private incomplete source graph. Source195ntcassetteN501/N701orN702sampleindexes differfromtwo7ntmoleculartags;pairedtagsidentifyduplexfamilies computationally. SelectedgraphmodelsCaOV3RCAenrichmentroute,notallTP53calibration/yeastinputs. DetailedMethodsalternativeT/AtailingusesKlenow whereasFig1caption saysTaq;thisdifferentvariantisnotmergedintogenomicblunt-endpreparationwithunassignedcloning/amplificationroute. Fullbackbone/genomic/covalentend/readgeometry andentrylibrarycloningbranchallocationunassembled;no perfectclonecapture/errorcorrection/yield.'
ENDPOINT='load_out'
@dataclass(frozen=True)
class TaggedCircleClass:
    left:str='N'*7
    right:str='N'*7
    insert:str|None=None
    def __post_init__(self):
        if any(len(x)!=7 or set(x)-set('ACGTN') for x in (self.left,self.right)):raise ValueError('Bothsource7nucleotidebarcodesrequired')
    @property
    def pair(self):return self.left,self.right
    def rca_copy(self):return RCAClass(self)
@dataclass(frozen=True)
class RCAClass:
    original:TaggedCircleClass
    new_strand_biotin_at_primer:bool=True
    biotin_per_repeat:bool=False
    original_circle_remains:bool=True
    def __post_init__(self):
        if self.biotin_per_repeat:raise ValueError('Biotincomesfromprimer,notNTPincorporationeachrepeat')

def cassette_inventory():
    sc=Scene();sc.strand('cassette',[Segment('unassembledlegacyNextera5arm','[legacy Nextera v1 arm]',placeholder=True),Segment('sourceleftmoleculartag','N'*7,placeholder=True,feature=feature('cypher.left','umi','random',note='Pairedwithright7nttagfor14basefamilyidentity')),Segment('sourceSmaIsite','CCCGGG'),Segment('sourcerightmoleculartag','N'*7,placeholder=True,feature=feature('cypher.right','umi','random',note='Pairedwithleft7nttag')),Segment('unassembledlegacyNextera3arm','[legacy Nextera v1 arm]',placeholder=True)],label='Sourcecassettecomponentinventory;sampleindexes/full195ntgeometryunassembled')
    return sc

def rca_primer():
    sc=Scene();sc.strand('primer',[Segment('sourcep53RCAprimer','GTGCAGCTGTGGGTTGATT')],mod5='biotin');return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);o=Operations(g)
    cassette=g.input('ultramer','Source195ntPAGEUltramerN701orN702;actualRCAindexallocationunreviewed',K.OLIGO,state=cassette_inventory())
    x=o.op('cassette_ds',cassette,'SourcePfuUltratwocyclesmakesdsDNAsequencingcassette;originalssDNA/primersremain','PfuUltra/Nexteraprimerpair')
    x=o.qc('smai_diagnostic',x,'SeparatesourceSmaIdigest/gelverificationaliquot;productioncassettecontinues')
    x=o.clean('cassette_clean',x,'ZymoDNAclean/concentrator5purification')
    x=o.op('cassette_cut',x,'SourceEcoRI-HF/BamHI-HFrestrictiondigestcreatescloning-readycassetteends','EcoRI-HF/BamHI-HF')
    vec=g.input('puc19','SourcepUC19backboneinventory;completeversion/sequenceunassembled',K.MOLECULAR)
    vec=o.op('vector_cut',vec,'SourcepUC19EcoRI/BamHIrestrictiondigest','EcoRI/BamHI')
    x=o.clean('cassette_gel',x,'Source1.5percentlowmeltagarosecassettebandextraction/Zymoclean')
    vec=o.clean('vector_gel',vec,'SourcepUC19backbonebandextraction/Zymoclean')
    x=o.join('empty_ligate',(vec,x),'SourceT4DNAligase1:6backbone/cassette16°Covernight;fullsealedgeometryunassembled')
    x=o.clean('empty_ethanol',x,'Sourceethanolprecipitationthenwaterresuspension')
    x=o.op('empty_transform',x,'SourceDH10Belectroporation,carbenicillinplates37°Covernight,scrapecolonies;newcellularcopiesnotsameinputmolecules')
    x=o.clean('empty_miniprep',x,'SourceQIAquickSpinMiniprepplasmidrecovery')
    x=o.qc('empty_measure',x,'Separatesourcespectrophotometry/diversitysequencingQC;productionemptyvectorsretained')
    x=o.op('insert_site',x,'SourceSmaIopensinsertionsite;actualRCAinputbranchblunt-versusT/Aprepallocationunassembled','SmaI')
    x=o.opaque('vector_prepare',x,'Sourcecloningpreparationallocationunassembled:blunt/phosphatase+bacterialrouteversusnondephosphorylatedT/A-PCRalternative;nohybridrecipe')
    dna=g.input('genomic','SourceCaOV3genomicDNAforRCAexperiment;independentp53WTspikelibrarypreparedseparately',K.MOLECULAR,primary=True)
    dna=o.op('shear150',dna,'CovarisS220targets150bpCaOV3fragments','CovarisS220')
    dna=o.clean('genomic_gel',dna,'SourcegelisolationofCaOV3fragmentpopulation')
    dna=o.op('blunt',dna,'SourceQuickBluntingKitendpreparation;completegenomicendgeometryunassembled','QuickBluntingKit')
    x=o.join('insert_ligate',(dna,x),'SourcecloneCaOV3genomicfragmentsintoCypherSeqvector;selectedRCAparagraphdoesnotreassignspecificcloning/amplificationbranch')
    x=o.opaque('cloned_library',x,'SourcepreparedcircularCaOV3genomiclibraryaftercloning;actualbacterialversusPCRroute/physicalworkupallocationunassembled')
    spike=g.input('wt_p53','IndependentpreviouslypreparedSKOV3WTp53exon5CypherSeqlibrary;calibrationconstruct/TOPOmanufactureoutsidethisRCAentryscope',K.MOLECULAR)
    x=o.join('spike',(x,spike),'SourceWTp53librarydilutionintoCaOV3genomiclibraryto1per10000targets;sourceinputrationotpredictedcapturepurity')
    p=g.input('rca_primer','PublishedSI6p53RCAprimer5biotin',K.OLIGO,state=rca_primer(),chemical_state='Source5primebiotin;newRCAstrandprimerincorporation,notbiotinoneveryrepeat')
    x=o.join('denature',(x,p),'Source1µMbiotinp53primer/AccuTaqbuffer,95°C4minthengradualcoolto63°C')
    x,rest=o.split('rca_allocate',x,'Source4µLdenaturedDNAinto16µLprewarmedRCAmixture','OtherdenaturedDNA/primerpreparation;laterfateunassigned')
    x=o.op('rca',x,'Bst2.0WarmStart/2mMeachdNTP/BSA63°C96hcreatesnewbiotin-primerinitiatedconcatemerstrands;originalcirclesandunreactedclassesremain','Bst2.0WarmStart120U',chemical_state=TaggedCircleClass().rca_copy())
    x=o.clean('affinity',x,'SourcekilobaseBINDERaffinitypurificationofbiotin-bearingRCApreparation;exactbinding/wash/elutiongeometryunreviewed')
    x=o.op('shear1kb',x,'SourceCovarisS220shearsRCAsamplesinto1kbfragmentspriortoquantification;separatesconcatemercopiesphysically','CovarisS220')
    x=o.qc('ddpcr',x,'Separatesourcep53/totalconstructddPCRanalyticalaliquots;productionmaterialretained')
    f=g.input('library_f','PublishedSI4LibraryPrimer1',K.OLIGO,state=Construct([Segment('P5sourceportion',il.P5[:20])]))
    r=g.input('library_r','PublishedSI4LibraryPrimer2',K.OLIGO,state=Construct([Segment('P7sourceportion',il.P7[:21])]))
    x=o.join('library_pcr',(x,f,r),'SourcesamplepreparedbylimitedPCRwithlibraryprimers;RCAparagraphdoesnotassignexactcyclecount')
    x=o.opaque('load',x,'SourceMiSeqlibraryloadingwithNexterav1.0Read1/Index/Read2primerspikedincartridge;exactpostPCRworkupandreadcyclesunassembled')
    return finish(g,x,'PublishedMiSeqwith6.8µL50µMNexterav1.0Read1/Index/Read2spikedintocartridgepositions12/13/14;fullprimerlandingandRCAreadcycleallocationunassembled. Source300-cycleNexteraXTreference-WGSrunnotassignedthisRCAexperiment.')

def selftest():
    from checks import Check
    c=Check();g=protocol_graph();original=TaggedCircleClass('ACGTACG','TGCATGC');copy=original.rca_copy()
    c('RCAcopypreservestagpair',(copy.original.left,copy.original.right),original.pair)
    c('Primerbiotinnotcopiedoneveryconcatemerrepeat',copy.biotin_per_repeat,False)
    c('Originalcirclekeptdistinctfromnewstrand',copy.original_circle_remains,True)
    c.raises('Wrongpairedtagwidthrejected',lambda:TaggedCircleClass('N'*6,'N'*7))
    c.raises('Chemicalbiotinoneveryrepeatnotallowed',lambda:RCAClass(original,biotin_per_repeat=True))
    c('DNAfragmentationoccursafteractualcapture',g.steps['shear1kb'].inputs[0].material,'affinity_kept')
    c('ProductionlibraryexcludesddPCRmeasurement',g.steps['library_pcr'].inputs[0].material,'ddpcr_production')
    c('SourceRCAreceivesactual4µLallocation',g.steps['rca'].inputs[0].material,'rca_allocate_production')
    src=Source(sources_root()/'cypherseq__10.1093+nar+gkv915'/'nar_gkv915_supp_supp_gkv915_nar-01384-met-k-2015-File012.xlsx.txt','PublishedSI4/SI6originalworkbook')
    if have(c,src,label='SourceCypherSeqprimertranscription'):
        t=src.path.read_text();rows=dict(x.split('	',1) for x in t.splitlines() if x.startswith(('Library Primer 1	','Library Primer 2	')));c('WholeSI4libraryprimerbodiesmatchsource',(il.P5[:20],il.P7[:21]),(rows['Library Primer 1'],rows['Library Primer 2']))
        c('WholeSI6biotinprimerretained','5’-biotin-GTGCAGCTGTGGGTTGATT-3’' in t,True)
    c('Fullgeometry/sourceallocationpreventfinishedclaim',g.publishable,False);c.report()
