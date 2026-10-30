"""Protocol-owned selected chemical classes; incomplete source-qualified graph."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, revcomp, feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(name,text):
    sc=Scene();sc.strand(name,[Segment(name,text,placeholder=True)]);return sc

def clean(g,key,x,label):
    condition=g.input(key+'_condition',label,K.CONDITION)
    kept=g.product(key+'_kept','Source recovered material; no perfect purity or quantitative partition assumed',K.MOLECULAR)
    waste=g.product(key+'_waste','Actually separated cleanup fractions with unquantified molecular carryover',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(kept,waste));g.discard(waste,'Source physical separation, not perfect molecular classification.');return kept

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source mixture',catalyst=False,product=None,state=None,chemical_state=None):
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else K.CONDITION)
    y=g.product(key+'_product',product or label+'; all physical reaction companions retained until separation',K.MOLECULAR,state=state,chemical_state=chemical_state)
    kwargs={}
    if kind==S.MAGIC:kwargs['reason']='Source operation established; exact internal chemistry/geometry not assembled from unreviewed historical originals.'
    if kind==S.PCR:kwargs['geometry']=PCRGeometry.unavailable(Scene(),reason='Exact source template/primer endpoint geometry unassembled, no ordinary complete library asserted.')
    if kind in (S.LIGATION,S.PROXIMITY_LIGATION):kwargs['geometry']=AdapterLigationGeometry('unavailable','unavailable',(scene('unassembled','[source-known ligation; full substrate/end geometry unassembled]'),),'Unreviewed geometry; not a claim of unavailable disclosure.')
    inputs=[(x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION)]
    if kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.PROXIMITY_LIGATION,S.REVERSE_TRANSCRIPTION,S.PCR} and not catalyst:
        enzyme=g.input(key+'_enzyme','Explicit catalytic component of source mixture: '+reagent,K.CATALYST);inputs.append((enzyme,R.CATALYST))
    if kind==S.REVERSE_TRANSCRIPTION:
        primer=g.input(key+'_primer','SourceRTprimer;exacthistoricalgeometryunassembled',K.OLIGO,state=scene(key+'_primer','[source RT primer]'));inputs.append((primer,R.STRUCTURAL))
    if kind==S.LIGATION:
        adapter=g.input(key+'_adapter','Source ligating molecular adapter; exact historical oligo/end inventory unassembled',K.ADAPTER,state=scene(key+'_adapter','[source adapter geometry unassembled]'));inputs.append((adapter,R.STRUCTURAL))
    if kind==S.PCR:
        for arm in ('f','r'):
            primer=g.input(key+'_primer_'+arm,'Source PCR primer family; exact oligo/endpoint unassembled',K.OLIGO,state=scene(key+'_primer_'+arm,'[source primer body/index geometry unassembled]'));inputs.append((primer,R.STRUCTURAL))
    g.react(key,label,kind,inputs,(y,),**kwargs);return y

def qc_split(g,key,x,label):
    assay=g.product(key+'_assay','Physically separate source analytical aliquot',K.QC)
    kept=g.product(key+'_retained','Unassayed retained preparation for downstream library work',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(assay,kept));g.mark_terminal(Terminal(assay.id,Disposition.QC,'Source analytical aliquot never feeds retained library.'));return kept
SOURCE='Published Science352aad9822 first-party finalPDF, originaltitle/DOI/datechecked, Methods1531–1668; associatedAMPMC5756622 Methods124–146 selectedcomparisononly. TableS9 originalgated.'
@dataclass(frozen=True)
class TaggedEnd:
    rna:str
    adapter:str
    def __post_init__(self):
        if not self.rna or set(self.rna)-set('ACGU') or not self.adapter or set(self.adapter)-set('ACGT'):raise ValueError('ExplicitRNAandDNAadapterrequired')
    def selected_first_copy(self,primer):
        # Caller supplies actual primer; only selected RNA extension upstream of adapter binding.
        if not primer or set(primer)-set('ACGT') or not primer.endswith(revcomp(self.adapter)):raise ValueError('ExplicitprimerfullycomplementarytoselectedDNAadapterregisterrequired')
        return primer+revcomp(self.rna.replace('U','T'))
    def selected_rna_fragment(self,start):
        if type(start)is not int or not 0<=start<len(self.rna):raise ValueError('Explicitretained3endfragmentrequired')
        return (self.rna[:start],TaggedEnd(self.rna[start:],self.adapter))
def allocate(g,key,x,label):
    used=g.product(key+'_used','Sourceallocatedfraction',K.MOLECULAR);rest=g.product(key+'_remainder','Sourceunallocatedremainder;fateunassigned',K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(used,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No source storage/discardinstructionassumed.'));return used

def graph_for_kit(kit='Ribo-ZeroMRZB12424'):
    if kit not in {'Ribo-ZeroMRZB12424','MICROBExpressAM1905'}:raise ValueError('Actualsourcekitbranchrequired')
    g=ProtocolGraph('Term-seq — final2016mainMethods/'+kit,authority=SOURCE,representation='published_instance')
    x=g.input('bacterial_pellet','Sourcefrozenbacterialpellet;organism/applicationnotindependentrecipe',K.MOLECULAR,primary=True)
    x=op(g,'extract',x,'SourceFastprep/FastRNAPROblueRNAextraction',reagent='Source116025050kit;historicalinternalsunreviewed')
    x=qc_split(g,'rna_qc',x,'SourceQubitRNA/Tapestationanalyticalaliquots')
    x=op(g,'dnase',x,'SourceTURBODNasetreatment',S.DIGESTION,reagent='TURBODNaseAM2238',catalyst=True,product='RNAplusdigestedDNAandenzymaticremnantsuntilactualcleanup')
    x=op(g,'first_ligation',x,'SourceRNA3endDNAadapterligation23C2.5h',S.LIGATION,reagent='SourceT4RNL1/ATP/DMSO/PEG8000;TableS9DNAadapterunreviewed',product='RNA-DNAadapterclassesplusunligatedRNA/freeDNAadapter/digestedDNA/reagentsuntilSPRI')
    x=clean(g,'first_spri',x,'Source2.2xSPRI,70%EtOH,elute5–10ul')
    x=op(g,'fragment',x,'SourceAmbionRNAfragmentation72C1.5min',S.FRAGMENTATION,reagent='Sourcefragmentationbuffer',product='Adapter-retaining3terminalRNAfragmentsplusotherRNAfragmentsandDNAadapter/remnants;fragmentationdoesnotdiscardunadaptedfragments')
    x=clean(g,'fragment_spri',x,'Source2.2xSPRIelute28ul')
    x=op(g,'ribodeplete',x,'Source'+kit+'rRNAdepletion',reagent='Historicalsourcekitinternals/separationfatesunassembled',product='SourceprocessedrRNA-depletedpreparationandunresolvedcompanionfractions;nomechanism/perfectdepletionorremnantdiscardinferred')
    x=op(g,'rt',x,'SourceadapterprimedAffinityScriptRT42C45min/75C15min',S.REVERSE_TRANSCRIPTION,reagent='SourceAffinityScript600559;TableS9primerunreviewed',product='cDNAplusoriginalRNA/DNAadaptertemplates,unextendedprimerandkitremnantsbeforeRNaseH')
    x=op(g,'rnase_h',x,'SourceRNaseH37C30min',S.DIGESTION,reagent='SourceRNaseHM0297',catalyst=True,product='cDNAplusRNAhydrolysisproducts,DNAadaptertemplates/freeprimersandkitremnantsuntilSPRI')
    x=clean(g,'cdna_spri',x,'Source2.2xSPRIelute5.5ul')
    x=allocate(g,'cdna_allocate',x,'Source5ulof5.5ulcDNAintosecondligation')
    x=op(g,'second_ligation',x,'SourcecDNA3endadapterligation23C4–8h',S.LIGATION,reagent='SourceT4RNL1conditionsasabove;distinctTableS9cDNAadapterunreviewed',product='Double-adaptercDNAclassesplusoriginalDNA/freeadapters/ligationremnants')
    x=clean(g,'second_spri',x,'Source1.8xSPRIelute23ul')
    x=allocate(g,'pcr_allocate',x,'Source22ulof23ulligatedcDNAforPCR')
    x=op(g,'pcr',x,'SourceKAPAHiFiKK2601PCR16–18cycles',S.PCR,reagent='SourceKAPAready-mix/actualTableS9primerfamiliesunreviewed',product='SourcePCRcopiesplusoriginalcDNA/DNAadapter/unextendedprimersandreagents')
    x=clean(g,'final_spri',x,'Source0.9xSPRI')
    x=qc_split(g,'library_qc',x,'SourceQubitdsDNA/Tapestationanalyticalaliquots')
    lib=g.product('library','SourceTerm-seqphysicalDNAsequencinglibrary;fullgeometryunassembled',K.LIBRARY);g.react('identity','Sourcefinalphysicallibraryidentity',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='PublishedmainlistsMiSeq/HiSeq1500/NextSeq500acrosslibraries;actualper-libraryreadcycles/index/primergeometryunreviewed. Mappedfirstbaseandreadthrougharecomputationalreadout,notphysicalRNAcleavageorproofofterminationorigin.'))
    g.add_source_item(SourceItem('pending','ActualTableS9/SOMandnamedhistoricalkitoriginals','TableS9gated;kitmechanism/geometryunreviewed,noacceptedabsence ornewvendorrecipe.'));return g
def protocol_graph():return graph_for_kit()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a=TaggedEnd('ACGUACGU','GACT');lost,b=a.selected_rna_fragment(3)
    c('FragmentationconservesRNAand3endDNAadapter',(lost+b.rna,b.adapter),(a.rna,a.adapter))
    c('FirstcDNAextendsprimerwithoutcopyingDNAadaptertwice',b.selected_first_copy('AGTC'),'AGTCACGTA')
    g=protocol_graph();h=graph_for_kit('MICROBExpressAM1905')
    c('RNaseHactsafterRTandbeforeactualSPRI',(g.steps['rnase_h'].inputs[0].material,g.steps['cdna_spri'].inputs[0].material),('rt_product','rnase_h_product'))
    c('LibraryPCRusesallocatedpreparation',g.steps['pcr'].inputs[0].material,'pcr_allocate_used')
    c('Unreviewedkitdoesnotinventaffinitydiscard',h.steps['ribodeplete'].kind,S.MAGIC)
    c('UnusedcDNAfateremainsunresolved',g.terminals['cdna_allocate_remainder'].disposition,Disposition.UNRESOLVED)
    c('Knownoperationsdonotpromoteincompletelibrary',bool(g.publication_problems()),True)
    src=Source(sources_root()/'term-seq__10.1126+science.aad9822'/'source_pdf_9a78b21ecf9cf7bedd9a.pdf.txt','get_sources --slugs term-seq__10.1126+science.aad9822')
    if have(c,src,label='publishedfinalMethods'):
        t=src.path.read_text();c('PublishedfinaldefinesdistinctRNAandcDNAadapterligations','cDNA-specific ligation adapter' in t and '3′ end–specific ligation' in t,True)
    c.report()
if __name__=='__main__':selftest()
