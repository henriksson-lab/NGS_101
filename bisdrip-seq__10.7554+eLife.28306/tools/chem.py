"""Source-qualified bisDRIP operations, private and incomplete."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition

@dataclass(frozen=True)
class SelectedDeamination:
    """Selected final DNA conversion witness, not reactivity prediction."""
    dna:str
    positions:tuple
    inferred:bool=True
    def __post_init__(self):
        if not self.dna or set(self.dna)-set('ACGT') or len(set(self.positions))!=len(self.positions) or any(type(i)is not int or not 0<=i<len(self.dna) or self.dna[i]!='C' for i in self.positions):raise ValueError('Concrete DNA and explicit distinct cytosine coordinates required')
    def dna_product(self):return ''.join('U' if i in self.positions else b for i,b in enumerate(self.dna))
    def selected_reported_reference(self):return self.dna_product().replace('U','T')

def protocol_graph():
    g=ProtocolGraph('bisDRIP-seq selected standard IP route',authority='Published eLife28306v2 Methods514–546,584–592; controls556–582 identified separately. Source-boundary partial graph; historical Pico-Methyl kit unreviewed.',representation='published_instance')
    x=g.input('cells','Independent source MCF7 dish; untreated/estrogen/triptolide study treatment choice unassigned here',K.BIOLOGICAL,primary=True)
    def op(key,label,kind,reagent,product):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.DEPROTEINIZATION};c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION);y=g.product(key+'_product',product,K.MOLECULAR)
        kw={'reason':'Actual named source component; complete molecular transition unassembled.'} if kind==S.MAGIC else {}
        g.react(key,label,kind,((x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)),(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Source physically recovered nucleic-acid preparation; yield/carryover unquantified',K.MOLECULAR);w=g.product(key+'_removed','Physical source separated workup fractions',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Actual physical source separation, no perfect molecular partition.');x=y
    def aliquot(key,label,selected,remaining):
        nonlocal x
        y=g.product(key+'_used',selected,K.MOLECULAR);rest=g.product(key+'_rest',remaining,K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Actual remaining source material; current selected route does not assign its ultimate fate.'));x=y
    def qc(key,label):
        nonlocal x
        y=g.product(key+'_retained','Source production material after diagnostic allocation',K.MOLECULAR);q=g.product(key+'_assay',label,K.QC);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,q));g.mark_terminal(Terminal(q.id,Disposition.QC,'Separate source analytical aliquot; no return to production.'));x=y
    clean('cell_wash','Source three PBS washes/scrape/pellet/PBS removal')
    op('bisulfite_lysis','Source bisulfite-containing lysis buffer1.5ml',S.MAGIC,'Source ammonium sulfite/bisulfite/hydroquinone/EDTA in degassedwater','Bisulfite-exposed cellular lysate with DNA/RNA/protein and possibleconversionintermediates; no perfecthybridprotection imposed')
    aliquot('lysis_tube','Selected one of two source475ul tube allocations','Selected475ul lysate','Other475ul parallel-source tube and550ul unassigned residual of1.5ml; no implicit disposal or source-wide pooling')
    op('initial_conversion','Source add25ul lithiumdodecylsulfate;37C30min1100rpm dark',S.CHEMICAL_CONVERSION,'Source25ulLDS;stockconcentrationunassigned','Selectedbisulfite-exposed preparation with uracil/sulfonatedintermediates and protected/unconvertedclasses; RNA/proteins/reagentsremain')
    op('concentrated_conversion','Source higherbisulfite37C2h1100rpm dark',S.CHEMICAL_CONVERSION,'Source1mlammoniumbisulfitesolution/26.8ul20mMhydroquinone','Further sourceconversionmixture; exactconversionfraction notpredicted')
    clean('dialysis','SourceAmiconUltra4 30kMW nominal200fold dialysis then5000fold intoPKbuffer;1mleluate')
    op('proteinase','Source proteinaseK overnight37C',S.DEPROTEINIZATION,'50ul20percentSDS/200ul20ugulProteinaseK','Proteolyzed DNA/RNA mixture with retainedproteinfragments/reagents')
    clean('organic_physical','Source two-tube phenol/chloroform/chloroform workup,2h4Cethanolprecipitation,transferandcombineprecipitatedDNA/70percentethanolwash;internalallocationsopaque')
    op('restriction','Source900ul restrictiondigest37Covernight',S.DIGESTION,'SourceHindIII/EcoRI/BsrGI/XbaI/SspI inNEB3.1','Restrictiondaughters andunprocessedDNA/RNA/reagents; full genomiccutcoordinatesunknown')
    qc('digest_gel','Source0.8percentagarosecheck;IPusedonlyifmostly1–10kb;allocationunassigned')
    qc('digest_quant','SourceQuantIT/SpectraMaxconcentrationassay;allocationunassigned')
    aliquot('ip_input','Source20ugDNAallocation tostandardIP','Source20ugdigestedDNA/RNA','RemainingdigestafterQC/IPallocation;massunassigned')
    antibody=g.input('s96','Source20ugS9.6antibody',K.SELECTION);cond=g.input('ip_buffer','Source1mlIPbuffer/overnight4Cendoverend',K.CONDITION);y=g.product('ip_complex','SourceS9.6associatedhybridsplusunboundDNA/RNA/antibodymixture',K.MOLECULAR);g.react('ip_bind','SourceS9.6binding',S.BINDING,((x,R.TEMPLATE),(antibody,R.SELECTION),(cond,R.CONDITION)),(y,));x=y
    beads=g.input('beads','Source150ulwashedDynabeadsProteinG',K.SELECTION);b=g.product('captured','Sourcebead-associatedIPfractionwithcarryover',K.MOLECULAR);w=g.product('capture_unbound','Actualsourceunboundfraction',K.WASTE);g.react('capture','Sourcebeadbinding2h4C/magneticseparation',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(b,w));g.discard(w,'Actual magnetic source separation; no perfect S9.6 specificity assumed.');x=b
    clean('capture_wash','Source3x750ulIPbufferwash10min4C/magneticseparation')
    op('complete_conversion','Source NaOHcompletion37C20min',S.CHEMICAL_CONVERSION,'150ul0.3NNaOH','Sourcecompletedbisulfite-treatedDNAclasses with RNA/beads/antibodies/reagents; C-to-UchemicalproductdistinctfromTreportedbysequencing')
    op('neutralize','Source neutralization',S.QUENCHING,'150ul0.3NHCl/17.5ul1MTrispH8','NeutralizedsourceDNA/RNA/matrixmixture')
    op('rnase_i','SourceRNaseI50C20min',S.DIGESTION,'0.5ul100UulRNaseI','SourceDNA plusoriginal/residualRNAandRNAfragments/matrix/remnants')
    op('rnase_h','SourceRNaseH37C20min',S.DIGESTION,'1ul0.5MMgCl2/3ul5UulRNaseH','SourceDNA andRNAcleavageproducts/matrix/reagents; no instantaneous physicalRNAremoval')
    op('elute_pk','SourceproteinaseKelution50C1h',S.DEPROTEINIZATION,'2ul20ugulProteinaseK','ElutedsourceDNAplusproteolysis/RNAdigestionproductsandbeadmatrix')
    clean('supernatant','Source3minmagnet/supernatanttransfer;matrixseparated')
    clean('post_ip_phenol','Sourceadd150ulwater/500ulphenolchloroform/5min18400g aqueousrecovery')
    aliquot('aqueous','Source400ulaqueoustransfer','Source400ulaqueousfraction','Untransferredaqueousfraction;totalvolume/recoveryunassigned')
    clean('post_ip_chloroform','Source400ulchloroform/2min18400g/supernatanttransfer')
    clean('precipitate','Source20ugglycogen-assistedethanolprecipitation48h')
    op('shear','SourceCovarisS2fragmentation~300bp',S.FRAGMENTATION,'150ul1xTE;4C/intensity5/duty10percent/200cyclesburst/140s','SonicatedDNAclasses containingselectedbisulfiteconversions plusothernucleicacids/remnants')
    op('kit','SourcePico-MethylSeq random-prime/index/amplificationcomponent',S.MAGIC,'HistoricalZymoPico-MethylSeqLibraryPrep;sourceexplicitlyOMITSkitbisulfite-treatmentsteps','IndexedDNA source librarypreparation; exactendgeometry/primers/copycyclesunassembled,no secondbisulfiteexposure imposed')
    qc('bioanalyzer','SourceAgilent2100BioanalyzerQC;allocationunassigned')
    qc('library_quant','SourceQuantIT/SpectraMaxQC;allocationunassigned')
    other=g.input('other_libraries','Sourceotherindependentlypreparedlibraries;5or6samplesperlane,membership/amountunassigned',K.LIBRARY,primary=True);pool=g.product('pool','PhysicalsourceIlluminalane loadingpool',K.LIBRARY);g.react('pool_load','Source5or6librariesperlane',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pool,))
    g.sequence(pool,SequencingPlan(limitation='SourceHiSeq2500PE100single-index. Exacthistoricalkitprimer/indexgeometryunassembled. BismarkreportsC-to-T/G-to-AoutcomesfromchemicalDNAuracil; computationalfiltering/deconvolutionneverphysicalDNAremoval.'))
    g.add_source_item(SourceItem('pending','HistoricalPico-MethylSeq guide/primer/index/endgeometry;source-linkedGinno2012/ Yu2003/Hayatsu2004','Source main standardIProuteinspected;controlroutesidentifiedbutnotfullymodeled;noacceptedabsence.'))
    return g

def selftest():
    from checks import Check
    c=Check();r=SelectedDeamination('ACGCC',(1,4));c('Physical product contains uracil not thymine',r.dna_product(),'AUGCU');c('Reported reference distinct from chemical product',r.selected_reported_reference(),'ATGCT');c('Original native DNA unchanged',r.dna,'ACGCC');c('Selected conversion inferred by default',r.inferred,True)
    for pos in ((0,),(5,),(1,1)):
        try:SelectedDeamination('ACGCC',pos);bad=False
        except ValueError:bad=True
        c('Invalid conversion coordinate rejected '+str(pos),bad,True)
    g=protocol_graph();c('RNAcleavage followschemicalcompletionandneutralization',g.steps['rnase_i'].inputs[0].material,'neutralize_product');c('RNaseHsubstrate retains RNaseIreaction',g.steps['rnase_h'].inputs[0].material,'rnase_i_product');c('Source20ugallocationleavesphysicalremainder',g.terminals['ip_input_rest'].disposition,Disposition.UNRESOLVED);c('ActualmatrixremovedonlyafterPKelution',g.steps['supernatant'].inputs[0].material,'elute_pk_product');c('Unreviewedkitpreventspromotion',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
