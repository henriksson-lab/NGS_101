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
    if kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.PROXIMITY_LIGATION,S.REVERSE_TRANSCRIPTION} and not catalyst:
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
SOURCE='AssociatedacceptedmanuscriptPMC3057491 Methods3P-SeqLibraries84–86 andoverview38–42. Finalfirst-partyPDFstandardretrievalcurl60;publishedfinalapplicabilityunreviewed. NoV2detailedprotocolimport.'
RNA_ARM='AGCGUGUAGGGCACCAU'
DNA_ARM='GCACATAC'
SPLINT='ATGGTGCCCTACACGCTTTTTTTTT'
RT_PRIMER='GTATGTGCATGGTGCCCTACACGCT'
from illumina import INDEX2_PRIMER_RC
SECOND=INDEX2_PRIMER_RC
@dataclass(frozen=True)
class MixedAdapter:
    rna:str=RNA_ARM
    dna:str=DNA_ARM
    def __post_init__(self):
        if not self.rna or set(self.rna)-set('ACGU') or not self.dna or set(self.dna)-set('ACGT'):raise ValueError('OrderedRNAandDNAblocksrequired')
    def scene(self):
        sc=Scene();sc.strand('adapter',[Segment('RNAarm',self.rna),Segment('DNAarm',self.dna)],mod5='phosphate',mod3='biotin');return sc
    def flatten_dna(self):raise ValueError('ActualRNAarmcannotbecomeordinaryDNA')
    def selected_rt_extension(self,polyA_length):
        if type(polyA_length)is not int or polyA_length<0:raise ValueError('ExplicitpolyAlengthrequired')
        return RT_PRIMER+'T'*polyA_length
@dataclass(frozen=True)
class SelectedRNaseHRelease:
    upstream:str
    original_tail:int
    retained_tail:int
    def __post_init__(self):
        if not self.upstream or set(self.upstream)-set('ACGU') or type(self.original_tail)is not int or type(self.retained_tail)is not int or not 0<=self.retained_tail<=self.original_tail:raise ValueError('ExplicitselectedRNA/tailcleavageboundaryrequired;notyieldprediction')
    @property
    def released_rna(self):return self.upstream+'A'*self.retained_tail
    @property
    def tail_residues(self):return 'A'*(self.original_tail-self.retained_tail)
def protocol_graph():
    g=ProtocolGraph('3P-seq — associatedoriginal2010/2011libraryMethods',authority=SOURCE,representation='published_instance')
    x=g.input('total_rna','Source30ugnematodetotalRNA;growth/extractionreference29unreviewed',K.MOLECULAR,primary=True)
    x=op(g,'poly_a',x,'SourceDynabeadsOligo(dT)25polyAenrichment',reagent='Historicalnamedkitinternals/fractiongeometryunreviewed',product='SourcepolyA-enrichedpreparation;kitremnantfatesunassembled')
    a=MixedAdapter();adapter=g.input('capture_adapter','Source5P-RNA/DNAchimera-3biotin',K.ADAPTER,state=a.scene(),chemical_state=a);splint=g.input('splint','SourceDNApolydTsplint',K.OLIGO,state=Scene());g.materials[splint.id].state.strand('splint',[Segment('sourceSplint',SPLINT)])
    enzyme=g.input('rnl2','SourceRnlIIRNAligase1U16h',K.CATALYST);mix=g.product('ligated','SourceRNA3endmixedadapterclassesplusunligatedRNA/adapters/splint/remnants',K.MOLECULAR)
    g.react('splint_ligate','SourcepolyAendsplintligation',S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(splint,R.STRUCTURAL),(enzyme,R.CATALYST)),(mix,),geometry=AdapterLigationGeometry('unavailable','unavailable',(a.scene(),),'ActualmixedadapterandDNA-splintknown;completeend/substratejunctiongeometryunassembled.'))
    x=clean(g,'ligation_workup',mix,'Sourcephenolextractionandprecipitation')
    x=op(g,'t1',x,'SourcepartialRNaseT122C20min',S.DIGESTION,reagent='SourceRNaseT1cutsafterG',catalyst=True,product='Adapter-bearing3terminalRNAfragmentsplusotherRNApieces/splint/adapter/remnants;nodigestiondiscard')
    beads=g.input('capture_beads','Source100ulM280streptavidinbeads',K.SELECTION);bound=g.product('bound','Sourcebiotin-adapterassociatedRNA/beadcomplexesplusunquantifiedcarryover',K.MOLECULAR);other=g.product('capture_other','Sourcephysicallyseparatedunboundfraction',K.WASTE)
    g.react('capture','SourceM280captureRT15min',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,other));g.discard(other,'Sourceboundfractionselected,notperfectpolyAclassification.')
    x=clean(g,'wash_bound',bound,'SourceBwashthen2Wwashes50C/equilibrateRTbuffer')
    primer=g.input('tail_rt_primer','Sourceadapter-complementaryRTprimer',K.OLIGO,state=scene('primer','source exact primer below'))
    # Actual primer is source transcription; complete strand placement unassembled.
    primer.state.strands.clear();primer.state.strand('primer',[Segment('sourceRT',RT_PRIMER)])
    dttp=g.input('dttp','SourcedTTPonlynucleotideforinitialRT',K.MOLECULAR);cat=g.input('ssiii','SourceSuperScriptIII',K.CATALYST);rt=g.product('tail_hybrid','Sourceselectedprimer-dTcopiesonRNApolyAtailplusoriginalRNA/mixedadapter/beads,unextendedprimers/remnants;upstreamRNAuncopiedinselectedtail-onlyclass',K.MOLECULAR)
    g.react('tail_rt','Source48C20mindTTP-onlyRT',S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(dttp,R.STRUCTURAL),(cat,R.CATALYST)),(rt,))
    x=op(g,'rnase_h',rt,'SourceRNaseH37C25minreleasetail-shortenedRNA',S.DIGESTION,reagent='SourceRNaseH',catalyst=True,product='ReleasedRNAclassespluscleavedRNAresidues,DNAprimer,mixedadapter/biotin/beadsandremnantsbeforephysicalcollection')
    released=g.product('released','SourcecollectedreleasedRNA;remainingterminalAsselectednotfixedcount',K.MOLECULAR);remain=g.product('bead_residue','Bead-associatedmixedadapter/primer/RNAremnantcompanions;laterfateunassigned',K.MOLECULAR);g.react('collect_release','Sourcephysicalreleasefractioncollection',S.CLEANUP,((x,R.TEMPLATE),),(released,remain));g.mark_terminal(Terminal(remain.id,Disposition.UNRESOLVED,'No laterbead/remnantfateinvented.'))
    x=clean(g,'rna_precip',released,'SourceprecipitationofreleasedRNA')
    a=g.input('second_adapter','Source5AppDNA3C3spaceradapter',K.ADAPTER,state=Scene());a.state.strand('second_adapter',[Segment('sourceSecondAdapter',SECOND)],mod5='rApp',mod3='C3spacer');cat=g.input('t4rnl1','SourceT4RNL110U22C2h',K.CATALYST);y=g.product('second_ligated','RNA/secondadapterclassesplusfreeadapter/unligatedRNA/remnants',K.MOLECULAR)
    g.react('second_ligation','SourcepreadenylatedDNAadapterligation',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(a.state,),'Sourceadaptermodsretained;fullRNAend/junctiongeometryunassembled.'))
    x=clean(g,'rna_gel',y,'Sourcegelpurify75–300ntRNAs')
    x=op(g,'library_finish',x,'Sourcecitedstrand-specificmRNAseqprotocol31finish',reagent='Actualcitedpredecessor/detailedhistoricalprotocolunreviewed;nomodernkitrecipe',product='SourcephysicalDNAsequencinglibrarypreparation;exactRT/PCR/adapterstepsunassembled')
    lib=g.product('library','Source3P-seqDNAsequencinglibrary;fullgeometryunassembled',K.LIBRARY);g.react('identity','Sourcefinallibraryidentity',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='AssociatedoriginalexplicitIllumina;source-sequencingbeginsatRNA3endadenylatecomplements,reverse-complementreadsusedastags. Exactrunlength/primergeometryunreviewed;latermiddle-imagebasecallingimprovementseparatefromoriginalruns.'))
    g.add_source_item(SourceItem('pending','Finalpublishedmainapplicability,citedprotocol31,historicaldetailed3Pprotocolrevision','AuthorfinalPDFretrievalcurl60;sourceAMexplicit. NoV2webprotocolsilentlysubstituted.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a=MixedAdapter();c('SourcechimeraactuallycontainsRNAandDNA',(a.rna[-1],a.dna),('U','GCACATAC'))
    try:a.flatten_dna();bad=False
    except ValueError:bad=True
    c('CannoteraseRNAarmtoordinaryDNA',bad,True)
    c('Tail-onlycopyretainssourceprimerandlimiteddTTPextension',a.selected_rt_extension(4),RT_PRIMER+'TTTT')
    r=SelectedRNaseHRelease('ACGU',12,2);c('Selectedreleaseconservestailresidueinventory',(r.released_rna+r.tail_residues),'ACGU'+'A'*12)
    g=protocol_graph();c('RNaseHmixtureretaineduntilphysicalcollection',g.steps['collect_release'].inputs[0].material,'rnase_h_product')
    c('Unassignedbeadresidueunresolved',g.terminals['bead_residue'].disposition,Disposition.UNRESOLVED)
    c('Unreviewedreference31notinventedPCR','pcr' not in g.steps,True)
    c('Incompleteversion/geometrypreventsfullpublication',bool(g.publication_problems()),True)
    src=Source(sources_root()/'3p-seq__10.1038+nature09616'/'nature09616_PMC3057491.xml.txt','get_sources --slugs3p-seq__10.1038+nature09616')
    if have(c,src,label='associatedAMMethods'):
        t=src.path.read_text();c('Actualmixedoligowrittenwithsugarcase','p-agcguguagggcaccauGCACATAC-Biotin' in t and SPLINT in t and RT_PRIMER in t,True)
    c.report()
if __name__=='__main__':selftest()
