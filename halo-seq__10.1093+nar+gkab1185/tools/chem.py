"""Protocol-owned source-qualified partial population graph, no finished marker."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition

def clean(g,key,x,procedure,retained='source-recovered preparation; exact molecular partition/carryover unquantified',*,molecular_reagent=False):
    condition=g.input(key+'_condition',procedure,K.MOLECULAR if molecular_reagent else K.CONDITION)
    out=g.product(key+'_kept',retained,K.MOLECULAR)
    waste=g.product(key+'_separated','Source physically separated fractions; no perfect molecular recovery/removal claimed',K.WASTE)
    g.react(key,procedure,S.CLEANUP,((x,R.TEMPLATE),(condition,R.STRUCTURAL if molecular_reagent else R.CONDITION)),(out,waste));g.discard(waste,'Actual source physical separation only.');return out

def op(g,key,x,label,kind,*,reagent=None,catalyst=False,product=None,reason='',changes_molecule=True,molecular_reagent=False,additional=()):
    inputs=[(x,R.TEMPLATE)]
    if reagent:
        stock=g.input(key+'_stock',reagent,K.CATALYST if catalyst else (K.MOLECULAR if molecular_reagent else K.CONDITION));inputs.append((stock,R.CATALYST if catalyst else (R.STRUCTURAL if molecular_reagent else R.CONDITION)))
    out=g.product(key+'_product',product or label+'; reaction products/remnants/carryover remain unassigned',K.MOLECULAR)
    g.react(key,label,kind,tuple(inputs)+tuple(additional),(out,),reason=reason,changes_molecule=changes_molecule);return out
from pcr import PCRGeometry
SOURCE='Engel et al.DOI10.1093/nar/gkab1185; published non-manuscript PMC8887463 Methods83–115,Results159–165; published SIlegends/figurelabels. Selected5minRNAseq photolabeling, not7mincomparison or later detailedprotocol.'
@dataclass(frozen=True)
class SelectedRNA:
    sequence:str
    marks:tuple[str,...]=()
    def __post_init__(self):
        if not self.marks:object.__setattr__(self,'marks',('native',)*len(self.sequence))
        if not self.sequence or set(self.sequence)-set('ACGU') or len(self.sequence)!=len(self.marks) or any(v not in {'native','alkynylated','clicked_biotin'} for v in self.marks):raise ValueError('explicit supported RNA sequence and retained per-residue chemical classes required')
    def alkynylate_selected(self,positions):
        selected=frozenset(positions)
        if any(type(i)is not int or not 0<=i<len(self.sequence) or self.marks[i]!='native' for i in selected):raise ValueError('selected supplied native residues only; this is not a prediction of photoactivity/site selectivity')
        return replace(self,marks=tuple('alkynylated' if i in selected else m for i,m in enumerate(self.marks)))
    def click_selected(self,positions):
        selected=frozenset(positions)
        if any(type(i)is not int or not 0<=i<len(self.sequence) or self.marks[i]!='alkynylated' for i in selected):raise ValueError('selected azide–alkyneclick requires actual current alkynylated sites')
        return replace(self,marks=tuple('clicked_biotin' if i in selected else m for i,m in enumerate(self.marks)))

def protocol_graph():
    g=ProtocolGraph('Halo-seq — published5minphotolabeling route',authority=SOURCE,representation='published_instance')
    current=g.input('cells','Source HeLa cells with appropriate HaloTagfusion; doxycycline1µg/ml48h expression, fusionlocalization not predicted',K.MOLECULAR,primary=True)
    dbf_environment=g.input('dbf_environment','Source1uM DBF-Halo ligand incubation15min37°C',K.CONDITION)
    current=op(g,'dbf',current,'DBFHaloligand1µM15min37°C',S.CHEMICAL_CONVERSION,reagent='Actual source DBF Halo ligand in HBSS',molecular_reagent=True,additional=((dbf_environment,R.CONDITION),),product='Selected covalent HaloTag–DBF proteinligand population plus unboundDBF/RNA/background; no RNA-biotin assigned')
    current=clean(g,'media_washes',current,'Source twicecomplete-mediawash10min37°C')
    current=op(g,'pa',current,'Propargylamine1mM5min37°C',S.MAGIC,reagent='Actual source1mM propargylamine/HBSS molecular contribution',molecular_reagent=True,reason='Source exposure, exact distribution unassembled.')
    current=op(g,'photolabel',current,'Source greenlight5min with HaloDBF/propargylamine',S.CHEMICAL_CONVERSION,reagent='Source green-light exposure5min; DBF and PA already supplied upstream, no fresh addition',product='Selected oxidized/PA-alkynylatedRNA classes plus unlabeled/background molecules and photochemistry products; adductstructure/site selectivity unassembled')
    current=clean(g,'trizol1',current,'SourceTrizolRNAextraction/20gaugeneedlehomogenization10passes',molecular_reagent=True)
    current=op(g,'dnase',current,'SourceDNaseI37°C30min',S.DIGESTION,reagent='source DNaseI',catalyst=True,product='RNA preparation plus DNAdigestionfragments, no physical discard at enzyme step')
    current=clean(g,'trizol2',current,'SourceTrizolrecovery afterDNase',molecular_reagent=True)
    supplies=g.input('click_supplies','Actual source CuAAC reagent inventory; packaging and number of stocks unassembled',K.MOLECULAR)
    biotin=g.product('click_biotin','Actual2mM biotin azide molecular contribution',K.MOLECULAR)
    ascorbate=g.product('click_ascorbate','Actual10mM freshly prepared sodium ascorbate contribution',K.MOLECULAR)
    thpta=g.product('click_thpta','Actual2mM THPTA ligand contribution',K.MOLECULAR)
    copper=g.product('click_copper','Actual100uM copper sulfate catalyst supply; active oxidation-state/speciation unassembled',K.CATALYST)
    tris=g.product('click_tris','Actual10mM Tris pH7.5 reaction-buffer component',K.CONDITION)
    g.react('click_components','Expose actual reaction-supply classes; inventory grouping, not a documented premix',S.GENERIC,((supplies,R.STRUCTURAL),),(biotin,ascorbate,thpta,copper,tris),changes_molecule=False)
    reacted=g.product('click_product','Selected clicked-biotinRNAadducts, residual alkynylated/unmodifiedRNA and freeclickreagent/products; no automatic physical cleanup or nucleotide deletion',K.MOLECULAR)
    environment=g.input('click_environment','Actual source30min dark25°C CuAAC environment',K.CONDITION)
    g.react('click','Actual source extracted-RNA CuAAC labeling',S.CHEMICAL_CONVERSION,((current,R.TEMPLATE),(biotin,R.STRUCTURAL),(ascorbate,R.STRUCTURAL),(thpta,R.STRUCTURAL),(copper,R.CATALYST),(tris,R.CONDITION),(environment,R.CONDITION)),(reacted,))
    current=reacted
    current=clean(g,'ethanol',current,'Source sodiumacetate/ethanolrecovery forlargescaleRNA; QuickRNAMinitestalternative notmerged',molecular_reagent=True)
    beads=g.input('beads','SourcePiercestreptavidinmagneticbeads2h4°C then3B&Wwashes5min',K.SELECTION)
    bound=g.product('bound','SourcebiotinRNA-associatedbeads plus possible nonspecific/carryovermolecules',K.MOLECULAR);waste=g.product('capture_waste','Physicalcaptureflowthrough/washfractions; exact RNArecovery unquantified',K.WASTE)
    g.react('capture','SourceRNAstreptavidincapture',S.AFFINITY_CAPTURE,((current,R.TEMPLATE),(beads,R.SELECTION)),(bound,waste));g.discard(waste,'Actual sourcephysicalmagneticseparation/washes.')
    current=op(g,'release',bound,'SourcePBS/Trizol37°C5–10min release',S.MAGIC,reagent='Actual source50ulPBS plus150ulTRIZOL; molecular supply grouping, packagingunassigned',molecular_reagent=True,reason='Source release known; exact biotin/adduct/bead dissociation geometry unassembled.')
    current=clean(g,'bead_separation',current,'SourceelutedRNArecovery withoutmagneticbeads')
    current=clean(g,'directzol',current,'SourceDirectZolkitRNAcleanup,10µlwaterelution; internal stepsopaque')
    library_input=g.product('rna_100ng','Actual100ng RNA allocated from source10ul eluate to library preparation',K.MOLECULAR)
    untransferred=g.product('rna_untransferred','Untransferred RNA stock from source10ul eluate; mass and later use unassembled',K.MOLECULAR)
    g.react('library_input_allocation','Actual100ng RNA library-input allocation',S.SPLIT,((current,R.TEMPLATE),),(library_input,untransferred))
    g.mark_terminal(Terminal(untransferred.id,Disposition.UNRESOLVED,'No total eluate RNA mass, volume allocation or subsequent stock fate inferred.'))
    current=library_input
    kit=g.input('kapa','Source rRNA-depletedRNAHyperPreplibrary(KAPA/Roche),100ngRNA; exact edition/rRNA-removal mechanism/internal stepsunreviewed',K.CONDITION)
    pre=g.product('pre_pcr','Sourcepre-PCRRNA-derivedpreparation; preciseDNA/adapter/strand/endgeometryunassembled and physicalremnantfatesunresolved',K.MOLECULAR)
    g.react('kit_handoff','SourcehistoricalrRNA-depletedRNAHyperPrep throughpre-PCRboundary',S.MAGIC,((current,R.TEMPLATE),(kit,R.CONDITION)),(pre,),reason='No proprietaryrRNA-removal/internalassemblymechanism or separationfates supplied.')
    f=g.input('pcr_forward','SourcekitPCRforwardmember, exactidentityunreviewed',K.OLIGO);r=g.input('pcr_reverse','SourcekitPCRreversemember, exactidentityunreviewed',K.OLIGO);enzyme=g.input('pcr_enzyme','SourcekitPCRcatalyst, exactidentityunreviewed',K.CATALYST)
    library=g.product('library','Source14cycleamplifiedlibrary plus originaltemplates/reactioncarryover; exact selectiongeometryunassembled',K.LIBRARY)
    g.react('pcr','Source14PCRcycles',S.PCR,((pre,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST)),(library,),geometry=PCRGeometry.unavailable(Scene(),reason='Historicalkitoligos/strand/endgeometryandcleanupunassembled.'))
    g.sequence(library,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':150,'Read 2':150},limitation='SourceNovaSeqPE150; exactfullRNA-derivedlibrary/runprimer/indexgeometryunassembled.'))
    g.add_source_item(SourceItem('pending_originals','Namedphotochemistrypredecessorsandhistoricalvendororiginals','MainMethods/Results159–165; KAPARNAHyperPrep,DirectZol and clickedRNAadductgeometryunreviewed'))
    return g

def selftest():
    from checks import Check
    c=Check();x=SelectedRNA('ACGU').alkynylate_selected({1,2});y=x.click_selected({2})
    c('alkynylation and click preserve physicalRNAnucleotidepositions',y.sequence,'ACGU')
    c('click converts only selected currentalkyne into biotinadduct',y.marks,('native','alkynylated','clicked_biotin','native'))
    c('freealkyne no longer exists as currentmark at clickedsite',y.marks[2]=='alkynylated',False)
    try:SelectedRNA('ACGU').click_selected({2});bad=False
    except ValueError:bad=True
    c('click cannot manufactureRNA-biotin withoutalkynylatedsubstrate',bad,True)
    c('untreatedRNAstate is not overwritten by selectedlabeling',SelectedRNA('ACGU').marks,('native',)*4)
    g=protocol_graph()
    c('DBF and PA are actual molecular contributors',(g.materials['dbf_stock'].kind,g.materials['pa_stock'].kind),(K.MOLECULAR,K.MOLECULAR))
    c('CuAAC contributions originate from actual grouped supplies',all(g.producer_of(k).id=='click_components' for k in ('click_biotin','click_ascorbate','click_thpta','click_copper','click_tris')),True)
    c('Photolabel step does not add fresh DBF or PA',tuple(p.material for p in g.steps['photolabel'].inputs),('pa_product','photolabel_stock'))
    c('Historical kit consumes actual100ng allocation',g.steps['kit_handoff'].inputs[0].material,'rna_100ng')
    c('Untransferred RNA stock has explicit unresolved fate',g.terminals['rna_untransferred'].disposition,Disposition.UNRESOLVED)
    c('historicalkit/readboundariespreventpublication',bool(protocol_graph().publication_problems()),True);c.report()
if __name__=='__main__':selftest()
