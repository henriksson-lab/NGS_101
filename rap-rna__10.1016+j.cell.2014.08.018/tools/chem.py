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
SOURCE='Associated accepted-manuscript PMC4177037 Methods137–163,overview69–91 andFig1caption403; finalpublicationapplicability unverified, actualSIgated. DOI10.1016/j.cell.2014.08.018.'
@dataclass(frozen=True)
class CrosslinkedRNA:
    sequence:str
    mode:str='amt'
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGU') or self.mode not in {'amt','fa','fa_dsg','minus_amt','noncrosslinked'}:raise ValueError('Explicit selected RNA sequence and sourcebranch required')
    def selected_rt_copy(self,copied_residues):
        if type(copied_residues)!=int or not 0<=copied_residues<=len(self.sequence):raise ValueError('Caller-selected copied extent required, no inferredRTkinetics')
        return RetainedRNAWithCopy(self,revcomp(self.sequence[len(self.sequence)-copied_residues:].replace('U','T')))
@dataclass(frozen=True)
class RetainedRNAWithCopy:
    template:CrosslinkedRNA
    dna_extension:str
    def __post_init__(self):
        if set(self.dna_extension)-set('ACGT'):raise ValueError('DNAcopy required')
    def scene(self):
        sc=Scene();sc.strand('original_rna',[Segment('originalRNA',self.template.sequence,note='Original RNA retained; crosslinksite/partnerunassigned')]);sc.strand('cdna_extension',[Segment('DNAextension',self.dna_extension,note='Conditional copiedextension only, primerbody/endgeometryunassembled')]);return sc

def graph_for_branch(branch='amt',short_target=True):
    if branch not in {'amt','fa','fa_dsg','minus_amt','noncrosslinked'}:raise ValueError('Actual sourcecrosslink/control branch required')
    g=ProtocolGraph('RAP-RNA — '+branch,authority=SOURCE,representation='published_instance')
    x=g.input('cells','Associated-manuscript sourceEScells; selectedRNA target unassigned',K.MOLECULAR,primary=True)
    if branch in {'amt','minus_amt'}:
        x=op(g,'photo',x,'Source0.5mg/mlAMT(orPBSminusAMTcontrol),longwaveUV7min',S.CROSS_LINKING,reagent='AMT+UV' if branch=='amt' else 'SourcePBS+UVcontrol; noAMTintroduced',product='SourceAMT RNAcrosslinkpopulation withunreactedRNA andcellcompanions retained' if branch=='amt' else 'SourceUVcontrolcellpopulation; noinstalledAMTadductasserted')
        x=op(g,'nuclear_enrich',x,'SourceAMT nuclearenrichment; exactphysicalpartitiondetails in gatedSI',reagent='SourceunreviewedAMTnuclearenrichmentprocedure')
        x=op(g,'trizol',x,'SourceAMTTRIzol isolation removesprotein/DNAassociation beforeRNAcapture',S.EXTRACTION,reagent='SourceTRIzolworkup')
        x=clean(g,'trizol_recover',x,'SourceAMTnucleicacidextraction/recovery')
        x=op(g,'fragment',x,'SourceAMTRNAfragmentation tosource100ntmedian',S.FRAGMENTATION,reagent='Sourcefragmentationconditions unreviewedSI')
        x=op(g,'dnase',x,'SourceAMTdigestresidualDNA afterfragmentation',S.DIGESTION,reagent='SourceDNase, completeends unassembled',product='SourcecrosslinkedRNA plusdigestedDNAproducts andotherreactioncompanions retained')
        x=op(g,'capture_input',x,'Source2ug resultingpurifiedRNA for15pmolprobecapture; cleanup detailsinSI',reagent='ExactDNasecleanup/purifiedinputpartitionunreviewed')
    else:
        if branch=='fa_dsg':x=op(g,'dsg',x,'Source2mMDSG before3%FA',S.CROSS_LINKING,reagent='SourceDSG',product='Sourceprotein-proteincrosslinkedcellclasses plusunreactedcompanions')
        if branch!='noncrosslinked':x=op(g,'fa',x,'Source2%FA (FAbranch) or3%FA (FA-DSGbranch)',S.CROSS_LINKING,reagent='Sourceformaldehyde',product='SourceRNA-protein/proteincomplexes plusunreactedcellmaterial; no directRNA-RNApsoralenbridgeasserted')
        x=op(g,'solubilize',x,'SourceFAstrongsonication orFA-DSGlightsonication',S.FRAGMENTATION,reagent='Sourcechosenbranchsonication; captureinputworkupgeometryunassembled')
        if branch=='fa_dsg':x=op(g,'dsg_dnase',x,'SourceFA-DSGDNaseI afterlights sonication',S.DIGESTION,reagent='SourceDNaseI',product='SourceRNA/proteincomplexes plusDNAcleavageproducts untilactualworkup')
    if short_target:
        probe=g.input('antisense_probe','Source3–4 purchased50nt5primebiotinDNAprobes toshortRNA; exacttargets/sequencesunreviewedTableS4',K.OLIGO,state=scene('biotin_probe','[5primebiotin antisenseDNAprobe; exactstudysequenceunassembled]'))
    else:
        precursor=g.input('probe_templates','Sourcebarcodedtiled120ntcomplementarityprobeDNAtemplates; exactTableS4/S5unreviewed',K.MOLECULAR)
        p=op(g,'probe_amplify',precursor,'SourcePCRamplifiedprobesubsets',S.PCR,reagent='SourceprobePCRcomponentsunassembled')
        p=op(g,'probe_transcribe',p,'Sourceprobeproduction IVT',reagent='SourceIVTreagents/polymeraseunreviewed')
        p=op(g,'probe_reverse',p,'Source5primebiotinDNAprimerRT generatesantisensessDNAprobes',reagent='SourceprobeRTreagents/biotinprimerunreviewed')
        probe=p
    hyb=g.product('hybridized','RNAcomplexes withbiotinDNAhybridprobe plusunboundRNA/DNAprobe/reactioncompanions',K.MOLECULAR)
    condition=g.input('hyb_conditions','SourceAMTalternate denaturingbuffers (SIunread)' if branch in {'amt','minus_amt'} else 'Source37C2h denaturingbuffer;5millioncells/50pmolprobe',K.CONDITION)
    g.react('hybridize','SourcespecificRNAantisensehybridization',S.HYBRIDIZATION,((x,R.TEMPLATE),(probe,R.STRUCTURAL),(condition,R.CONDITION)),(hyb,))
    beads=g.input('beads','SourceFig1/method159streptavidin captureofbiotinDNAprobecomplexes',K.SELECTION);bound=g.product('bound','Probe/bead-associatedRNAcomplexpopulation; no transferofprobeBiotin ontoRNA',K.MOLECULAR);rest=g.product('unbound','Physically separatedcapturefraction; specificcomposition/laterfateunresolved',K.MOLECULAR)
    g.react('capture','Sourcecaptureandwashinterface',S.AFFINITY_CAPTURE,((hyb,R.TEMPLATE),(beads,R.SELECTION)),(bound,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Sourceworkupcompanionfraction; no perfectspecificity ordownstreamuse inferred.'))
    x=clean(g,'wash',bound,'SourceFAbranches six45C washes;AMTalternatewashbufferdetailsunread')
    if branch in {'amt','minus_amt'}:
        x=op(g,'amt_release',x,'SourceAMTcapturedRNArecoveryinterface; releasechemistryunreviewedSI',reagent='ExactAMTrelease/endgeometrynotread; noRNaseHelutionborrowedfromFA')
        x=op(g,'amt_library_entry',x,'SourceAMTRNAtoRTsubstrate preparation; completeprimer/adapterorder inSIunreviewed',reagent='SourceRNAadapter/primerinventoryunreviewed')
        x=op(g,'rt',x,'SourceFig1AMTRTwithoutreversalofcrosslinks; recoveryofstoppedcDNA',reagent='SourceRT/primeridentity andends unreviewed',product='Selectedconditionalstopped/runoffcDNA andoriginalRNA/crosslinks/primers/reactioncompanions retained; stopdoesnotdeleteRNA')
        x=op(g,'cdna_adapter',x,'SourceFig1secondadapterligation tocDNA3primeend',S.LIGATION,reagent='SourcecDNA3primeadapter/ligase inventoryunreviewed',product='Sourceend-adaptedcDNApopulation withoriginalRNA/unligatedadapters/reactionremnants untilsourceworkup')
    else:
        x=op(g,'rnaseh_elute',x,'SourceFAbranchesRNaseHelution oftargetRNA/co-purifyingcomplexes',S.DIGESTION,reagent='SourceRNaseH',product='Releasedcomplexes withRNA-DNAhybridcleavageproducts/probe/beads andenzymes untilactualseparation')
        x=op(g,'proteinase',x,'SourceFAbranchesProteinaseKdigestion',S.DEPROTEINIZATION,reagent='SourceProteinaseK',product='RNA/probe/DNA/proteinremnantmixture; exactsourcecrosslinkreversalunassembled')
        x=clean(g,'recover_rna',x,'Sourceisolationofnucleicacids afterRNaseH/ProteinaseK')
    x=op(g,'library_finish',x,'SourceRNAsequencinglibrary completion; detailedSIandpublishedfinalgeometryunreviewed',reagent='ActualsourceRNAadapter/RT/PCR/workup detailsunreviewed, not inferredfromRAP-DNA')
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation='AssociatedMSMethods163 reports10–50millionreadpairs; sourceinstrument/cycles/primers unreviewed. AMTstoppedcDNA identifiesRNAinteractionproximity, notphysicalRNAdeletion. Finalpublishedapplicability unverified.'))
    g.add_source_item(SourceItem('pending','Publishedfinalarticle; actualSupplementalExperimentalProceduresDOCX, TablesS4/S5; authorsRAPprotocol andEngreitz2013 predecessor','ActualSIlinksPMC-gated; finalsource/version andfullprobe/end/library/acquisitiongeometry unresolved.'))
    return g

def protocol_graph():return graph_for_branch()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();rna=CrosslinkedRNA('ACGUUGCA');copy=rna.selected_rt_copy(3)
    c('SelectedRTstopretainsfullphysicalRNA',copy.template.sequence,'ACGUUGCA')
    c('ConditionalRTcopiesthreeRNApositionsintoDNA',copy.dna_extension,'TGC')
    c('NoRTbasescopiedproducesemptyextensionnotdeletedtemplate',rna.selected_rt_copy(0).dna_extension,'')
    a=graph_for_branch();f=graph_for_branch('fa');d=graph_for_branch('fa_dsg')
    c('AMTproteinremovalprecedesRNAcapture',tuple(a.steps).index('trizol_recover')<tuple(a.steps).index('capture'),True)
    c('FARNaseHelutiondoesnotgetborrowedintoAMT','rnaseh_elute' in f.steps and 'rnaseh_elute' not in a.steps,True)
    c('DSGbranchalonehasprotein-proteinprestep','dsg' in d.steps and 'dsg' not in f.steps,True)
    c('UVcontrolcontainsnoAMTinput','noAMTintroduced' in graph_for_branch('minus_amt').materials['photo_reagent'].label,True)
    p=Source(sources_root()/'rap-rna__10.1016+j.cell.2014.08.018'/'j.cell.2014.08.018_PMC4177037.xml.txt','get_sources --slugs rap-rna__10.1016+j.cell.2014.08.018')
    if have(c,p,label='associatedMSsource'):
        t=p.path.read_text();c('SourcecDNAadapterrecoverywithoutRNAcrosslinkreversal','without reversal of crosslinks' in t and '3’ end of the cDNA' in t,True)
    c.report()
if __name__=='__main__':selftest()
