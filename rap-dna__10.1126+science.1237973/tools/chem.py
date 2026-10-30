"""Protocol-owned selected chemical classes; incomplete source-qualified graph."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment, revcomp, feature
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
    if kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.PROXIMITY_LIGATION} and not catalyst:
        enzyme=g.input(key+'_enzyme','Explicit catalytic component of source mixture: '+reagent,K.CATALYST);inputs.append((enzyme,R.CATALYST))
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
SOURCE='PublishedScience3411237973 first-partyCaltechPDF source_pdf_61caebc669c26ad8bbc5.pdf printedtitle/DOI/date verified, Methods655–699 (pp5–6)/Fig1. AssociatedAMPMC3778663 Methods138 comparedonlyselectedscope. ActualSI .doc gated.'
@dataclass(frozen=True)
class AntisenseRNAProbe:
    region:str
    biotin_u_positions:tuple
    def __post_init__(self):
        if not self.region or set(self.region)-set('ACGU'):raise ValueError('ExplicitRNAcapture-region,notDNAprobe required')
        if len(set(self.biotin_u_positions))!=len(self.biotin_u_positions) or any(type(i)!=int or not 0<=i<len(self.region) or self.region[i]!='U' for i in self.biotin_u_positions):raise ValueError('SelectedbiotinUTPpositions mustidentifyRNAU; incorporationfractionunpredicted')
    def target_reference(self):return revcomp(self.region.replace('U','T')).replace('T','U')
    def capture_association(self,dna_id,rna_id):
        if not dna_id or not rna_id:raise ValueError('IndependentcapturedDNAandendogenousRNAidentitiesrequired')
        return {'probe_region':self.region,'probe_biotin_u':self.biotin_u_positions,'captured_dna':dna_id,'endogenous_rna':rna_id,'dna_direct_biotin':()}

def protocol_graph():
    g=ProtocolGraph('RAP-DNA — published2013RNA-probe capture',authority=SOURCE,representation='published_instance')
    cells=g.input('cells','Sourcecrosslinkedsamplepreparation; XistRNA targetcontext, exactsampleunassigned',K.MOLECULAR,primary=True)
    templates=g.input('probe_dna_templates','SourcemicroarrayDNAoligopool with120ntantisensetargetregions,15nttiling; exactmembers/sequence/endgeometryunreviewed',K.MOLECULAR)
    p=op(g,'probe_pcr',templates,'SourcePCRincorporatesT7promoter intoDNAtargetpool',S.PCR,reagent='SourceprobePCRprimer/polymerasecomponentsunreviewed')
    utp=g.input('biotin_utp','Sourcebiotin-UTP structuralRNAbuildingblock; fractionandpreciseadductunreviewed',K.MOLECULAR);t7=g.input('t7','SourceT7IVTpolymerasecomponent; exactrecipeunreviewed',K.CATALYST);probe=g.product('rna_probe_mix','Biotin-UTP-containingantisenseRNAprobes plusuntranscribedDNA/nucleotide/enzymeproducts untilactualprobepreparationworkup; notssDNA',K.MOLECULAR)
    g.react('probe_ivt','SourceIVTbiotinylatedRNAprobes',S.TRANSCRIPTION,((p,R.TEMPLATE),(utp,R.STRUCTURAL),(t7,R.CATALYST)),(probe,))
    probe=op(g,'probe_preparation',probe,'SourceRNAcaptureprobestockpreparation; detailedoriginalSIunreviewed',reagent='ActualprobeDNase/purification/fragmentends remainunassembled; nopuremixtureinferred')
    x=op(g,'dsg',cells,'Source2mMDSG45min',S.CROSS_LINKING,reagent='SourceDSG')
    x=op(g,'fa',x,'Source3%formaldehyde10min',S.CROSS_LINKING,reagent='Sourceformaldehyde')
    x=op(g,'lyse',x,'Sourcecelllysis/chromatinfragmentationpreparation',S.EXTRACTION,reagent='Sourcebufferdetails inactualSIunreviewed')
    x=op(g,'sonicate',x,'SourcesonicationplusTURBODNase100–300bpchromatinfragments',S.FRAGMENTATION,reagent='Sourcesonication/TURBODNase',product='FragmentedcrosslinkedDNA/RNA/protein withDNAcleavageproducts andothercompanions; no exactcut/enddistributionasserted')
    x=op(g,'denaturing_conditions',x,'Source3Mguanidinethiocyanatehybridizationbuffer',S.GENERIC,reagent='Sourceguanidinethiocyanate/buffer')
    beads=g.input('preclear_beads','Sourcestreptavidinmagneticpreclearbeads',K.SELECTION);kept=g.product('precleared','SourcepreclearedlysateforRNAprobehybridization; noidealnonbindingpopulationasserted',K.MOLECULAR);bound=g.product('preclear_bound','Physicalpreclearbead-associatedcompanionmaterial, laterfateunreviewed',K.MOLECULAR)
    g.react('preclear','Source45C20minstreptavidinpreclear',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(kept,bound));g.mark_terminal(Terminal(bound.id,Disposition.UNRESOLVED,'Actualpreclearcompanionfraction; exactfate/unwantedcompositionunassembled.'))
    condition=g.input('hyb_conditions','Sourceheatedlysate/RNAprobes45C2h',K.CONDITION);y=g.product('hybridized','RNAprobes/endogenousRNA/chromatinDNAcomplexpopulation plusfreeprobes/background/remnants',K.MOLECULAR)
    g.react('hybridize','SourceRNAantisensecapturehybridization',S.HYBRIDIZATION,((kept,R.TEMPLATE),(probe,R.STRUCTURAL),(condition,R.CONDITION)),(y,))
    beads=g.input('capture_beads','SourcestreptavidincoatedbeadsforbiotinRNAprobecomplexes',K.SELECTION);yes=g.product('bound','Sourceprobe-boundRNA/chromatinDNA/proteinbeadpreparation',K.MOLECULAR);no=g.product('capture_unbound','Physicallyseparatedunboundcapturecompanion, sourceexactfateunassigned',K.MOLECULAR)
    g.react('capture','Sourceprobecomplexstreptavidincapture',S.AFFINITY_CAPTURE,((y,R.TEMPLATE),(beads,R.SELECTION)),(yes,no));g.mark_terminal(Terminal(no.id,Disposition.UNRESOLVED,'Noidealpurity orsourceunstateddownstreamfate.'))
    x=clean(g,'wash',yes,'Sourcesix45Cwashes retainingbeadprobe/chromatinpreparation')
    x=op(g,'proteinase_release',x,'SourceProteinaseK65Covernightrelease/de-crosslinking',S.DEPROTEINIZATION,reagent='SourceProteinaseK/65Cconditions',product='ReleasedchromatinDNAplusendogenousRNA/RNAprobes/protein/peptide/beadcompanions untilactualworkup; capturedDNAdoesnotinheritprobebiotin')
    x=op(g,'dna_library',x,'SourcestandardIlluminaDNA-librarypreparationincludingrequiredDNAsampleworkup; detailedSIunreviewed',reagent='ActualDNAisolation/endprep/adapter/PCR/indexgeometryunreviewed; noRNA-libraryorRAPRNA2014recipeborrowed')
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation='PublishedMethods25bp paired-endreads, >5millionpersample; exactinstrument/index/primer/endgeometryunassembled. DNAcontacts areRNA-associatedenrichment, notproofdirectRNA-DNAbasepairing.'))
    g.add_source_item(SourceItem('pending','PublishedScienceSIandactualassociatedSupplementary_Material.doc; authorsRAPprotocol andhistoricalIlluminakitguides','Corefinalpublishedoperationsread, completeprobe/manufacturing/end/librarygeometry andcontrolworkupsremainunassembled; actualSIroutePMC-gated, noacceptedabsence.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();p=AntisenseRNAProbe('AUGCUU',(1,4));relation=p.capture_association('dna','rna')
    c('AntisenseRNAprobecomputesRNAcomplement',p.target_reference(),'AAGCAU')
    c('CapturedoesnottransferprobeBiotinontoDNA',relation['dna_direct_biotin'],())
    c('RNAprobehasselectedmodifieduridines',relation['probe_biotin_u'],(1,4))
    try:AntisenseRNAProbe('ATGC',());bad=False
    except ValueError:bad=True
    c('OriginalRAPRNAprobecannotbesilentlyreplacedbyDNA',bad,True)
    try:AntisenseRNAProbe('ACGU',(0,));bad=False
    except ValueError:bad=True
    c('BiotinUTPmarkmustidentifyU',bad,True)
    g=protocol_graph();c('PreclearbeforeRNAhybridization',g.steps['hybridize'].inputs[0].material,'precleared')
    c('2013releaseusesPKnot2014RNaseH','proteinase_release' in g.steps and 'rnaseh_elute' not in g.steps,True)
    q=Source(sources_root()/'rap-dna__10.1126+science.1237973'/'source_pdf_61caebc669c26ad8bbc5.pdf.txt','get_sources --slugs rap-dna__10.1126+science.1237973; exactcoordinator-declaredfirstpartyPDF')
    if have(c,q,label='publishedfinal'):
        t=q.path.read_text();c('FinalpublishedRNAprobe/pairedreadsource','antisense RNA' in t and '25-bp paired-end' in t and 'proteinase K' in t,True)
    c.report()
if __name__=='__main__':selftest()
