"""Published platinum Pt-seq; not phosphorothioate PT-seq."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,revcomp
from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2,STEM,STEM_COMPLEMENT,TRIM_SEEN_IN_READ1,INDEX2_PRIMER_RC
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

@dataclass(frozen=True)
class OrderedDNA:
    sequence:str
    ps_after:tuple=()
    five_prime:str='unassigned'
    three_prime:str='unassigned'
    index_start:int|None=None
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGTN') or any(type(i)is not int or not 0<=i<len(self.sequence)-1 for i in self.ps_after) or len(set(self.ps_after))!=len(self.ps_after):raise ValueError('Actual DNA and internal distinct PS bond coordinates required')
        if self.five_prime not in {'P','unassigned'} or self.three_prime not in {'NH2 C7','unassigned'}:raise ValueError('Source-qualified terminal chemistry required')
        if self.index_start is not None and (type(self.index_start)is not int or self.sequence[self.index_start:self.index_start+6]!='NNNNNN'):raise ValueError('Declared source6nt index window required')
    def starred(self):return ''.join(b+('*' if i in self.ps_after else '') for i,b in enumerate(self.sequence))
    def scene(self,name):
        sc=Scene();parts=[Segment(name,self.sequence,note='Source PS bonds after 0-based positions '+str(self.ps_after))]
        if self.index_start is not None:
            i=self.index_start;parts=[Segment('P7 arm',self.sequence[:i]),Segment('sample index',self.sequence[i:i+6],placeholder=True,feature=feature('pt_sample','sample_index','unknown')),Segment('source primer arm',self.sequence[i+6:])]
        sc.strand(name,parts,mod5='P' if self.five_prime=='P' else '',mod3='NH2 C7' if self.three_prime=='NH2 C7' else '');return sc

def ordering_oligo(name):
    u='T'*9+'CCTGCAGG'+STEM_COMPLEMENT+'T';b=TRIM_SEEN_IN_READ1[1:]
    return {'first_u':OrderedDNA(u,tuple(range(6))+(len(u)-2,)), 'first_b':OrderedDNA(b,tuple(range(len(b)-7,len(b)-1)),'P','NH2 C7'), 'extension':OrderedDNA(TRUSEQ_READ2[2:]), 'second_u':OrderedDNA(TRUSEQ_READ1+'N'*6,three_prime='NH2 C7'), 'second_b':OrderedDNA(INDEX2_PRIMER_RC,five_prime='P',three_prime='NH2 C7'), 'universal':OrderedDNA(P5+TRUSEQ_READ1[4:20]), 'index':OrderedDNA(P7+'N'*6+TRUSEQ_READ2,index_start=len(P7))}[name]

def first_adapter_scene():
    """Matched inventory core only; unpaired arms stay in full oligo inventories."""
    sc=Scene();sc.strand('u',[Segment('u_stem',STEM_COMPLEMENT)],label='Source first-adapter matched core ONLY; unpaired arms and ends are in separate ordered inventories')
    sc.anneal('b',[Segment('b_stem',STEM)],to='u',pair=('b_stem','u_stem'));return sc

def unknown_scene():
    sc=Scene();sc.strand('unknown',[Segment('geometry','[source full substrate/adapter/end geometry unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph('Pt-seq — platinum-DNA adduct mapping',authority='Published2025 PMC12571083 Results111–113/Methods193–205;publishedSI FigS1 originalp2;TableS1 originalxlsx rows5–11. DistinctfromphosphorothioatePT-seq.',representation='published_instance')
    x=g.input('genomic','Source extracted500–1000ng genomic DNA from drug/control cells; selected drug/cell/run mapping unassigned',K.MOLECULAR,primary=True)
    def op(key,label,kind,reagent,product,structural=()):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.THERMAL_TREATMENT,S.PCR,S.LIGATION,S.DEPROTEINIZATION} and kind!=S.THERMAL_TREATMENT
        c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION);inputs=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]+[(a,R.STRUCTURAL) for a in structural]
        kw={}
        if kind==S.MAGIC:kw['reason']='Actual source operation/component; complete molecular products and geometry unassembled.'
        if kind==S.LIGATION:kw['geometry']=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Source exact oligo/end chemistry inventory known; full ligation-ready substrate and sealed junction unassembled.')
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown_scene(),reason='Whole sourcetemplate/primer endpoint geometry unassembled; no predicted universal amplicon.')
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,inputs,(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Source retained physical preparation; recovery/carryover unquantified',K.MOLECULAR);w=g.product(key+'_removed','Actual physical source separated fractions',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Source physical cleanup; no perfect product classification.');x=y
    def reagent(name):
        a=ordering_oligo(name);return g.input('oligo_'+name,'PublishedTableS1 '+name+' DNA ordering inventory',K.OLIGO,state=a.scene(name),chemical_state=a)
    def assemble(name,a,b):
        c=g.input(name+'_conditions','Source gradual-cooling annealing; temperature programme unassigned',K.CONDITION);p=g.product(name+'_adapter','Source two-oligo adapter inventory; fullendgeometry unassembled',K.ADAPTER,state=first_adapter_scene() if name=='first' else None);g.react(name+'_anneal','Source adapter preparation',S.HYBRIDIZATION,((a,R.TEMPLATE),(b,R.STRUCTURAL),(c,R.CONDITION)),(p,));return p
    def qc(key,label):
        nonlocal x
        y=g.product(key+'_retained','Source production material',K.MOLECULAR);q=g.product(key+'_aliquot',label,K.QC);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,q));g.mark_terminal(Terminal(q.id,Disposition.QC,'Independent analytical aliquot; no return to library.'));x=y
    op('fragment','Source CovarisS220fragmentation median300bp',S.FRAGMENTATION,'Sourceultrasonication','Genomic fragments plusunfragmentedDNA; adductlocation/state unassigned')
    model=g.input('model','Source1ng modelsequencepreparation;exactspecies/platinatedfractionunassigned;independent preparation boundary',K.MOLECULAR,primary=True);pool=g.product('model_pool','Physical genomicDNA/modelsequencepool',K.MOLECULAR);g.react('model_spike','Source1ngmodeladdition beforelibrarypreparation',S.POOL,((x,R.TEMPLATE),(model,R.TEMPLATE)),(pool,));x=pool
    op('end_prep','Source NEBE7546 endrepair/dAtailingcomponent',S.END_REPAIR,'HistoricalE7546 module;conditions/reactioninternalsunreviewed','Source end-conditioned DNA plusunprocessedsubstrates/reagents; fulladduct/endsgeometryunassembled')
    first=assemble('first',reagent('first_u'),reagent('first_b'));op('first_ligation','Source NEBE7595 restriction-siteadapterligation',S.LIGATION,'SourceE7595ligationmodule', 'First-adapterligatedclasses plusunligatedDNA/freeadapters/reagents; PSbondand3amineinventoryretained,ligationsuccessunpredicted', (first,))
    op('denature_first','Source denaturation/ice beforeIP',S.THERMAL_TREATMENT,'Sourcetemperatures/durationunassigned','Source single-strandclasses plusresidualduplexes/freeadapters/remnants; no newchemicalcleavage')
    used=g.product('ip_used','Source DNAallocationtoantiPtIP;amountunassigned',K.MOLECULAR);ctrl=g.product('input_control','Source nonimmunoprecipitatedinputcontrol',K.MOLECULAR);g.react('input_split','Source saveinput beforeimmunoprecipitation',S.SPLIT,((x,R.TEMPLATE),),(used,ctrl));g.mark_terminal(Terminal(ctrl.id,Disposition.STORED,'Source savedinput at minus20C; not pooledbackintoIP.'));x=used
    ab=g.input('anti_pt','Source1ugantiPtDNAantibodyAbcamab103261',K.SELECTION);c=g.input('ip_conditions','SourceIPbuffer/overnight4Cgentlerotation',K.CONDITION);y=g.product('ip_complex','Sourceantibody/DNAmixture; nonPtbackgroundandunboundDNAremain',K.MOLECULAR);g.react('antibody','Source antibody binding',S.BINDING,((x,R.TEMPLATE),(ab,R.SELECTION),(c,R.CONDITION)),(y,));x=y
    beads=g.input('protein_g','Sourceprewashed0.1percentBSA/PBSblockedDynabeadsProteinG',K.SELECTION);b=g.product('bead_bound','SourcebeadassociatedDNAmaterial',K.MOLECULAR);w=g.product('ip_unbound','Sourcephysicallyseparatedunboundfraction',K.WASTE);g.react('beads','Source2h4Cbeadcapture',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(b,w));g.discard(w,'Source physical IP separation; specificity/recovery not idealized.');x=b
    clean('ip_wash','Source3xIPbufferwashes')
    op('pk_elute','SourceproteinaseKelution56C2h',S.DEPROTEINIZATION,'ProteinaseK/freshelutionbuffer','ElutedsourceDNAplusmatrix/proteolysisproducts/reagents')
    clean('elution_purification','Sourcepostproteinaseelutionpurification;fullinternalsunreviewed')
    op('exonucleases','SourceRecJf/exonucleaseI treatment',S.DIGESTION,'SourceRecJfM0264/exonucleaseIM0293','Sourceprotected/unprocessedDNAclasses plusdegradedDNAfragmentsandnucleotides/reagents; no physicalinstantremovalorperfectPSprotection')
    primer=reagent('extension');op('primer_extension','SourceQ5primerextension',S.EXTENSION,'SourceQ5M0544/dNTPconditionsunassigned','SelectedPt-stallednewDNAandnonstalledcopies plusoriginalDNA/unextendedprimer/reagents; no predictedstopposition',(primer,))
    op('primer_exo','Source exonucleaseI excessprimerdigestion',S.DIGESTION,'SourceexonucleaseIM0293','ExtendedandoriginalDNAplusprimer-digestproducts/reagents; no interveningcleanup invented')
    op('denature_second','Source second denaturation/ice',S.THERMAL_TREATMENT,'Sourcetemperatures/durationunassigned','Source separatedDNAclasses andreactionremnants')
    second=assemble('second',reagent('second_u'),reagent('second_b'));op('second_ligation','Source secondadapterligation16Covernight',S.LIGATION,'SourceT4DNAligaseM0202/PEG8000','Second-adapterligatedDNAplusunligatedDNA/freeadapters/3amineblockedoligos/reactionremnants; exactNN6annealinggeometryunassembled',(second,))
    f=reagent('universal');r=reagent('index');op('pcr1','SourceQ5 initial2PCRcycles',S.PCR,'SourceQ5M0544','Source amplified/newDNAplusoriginalDNA/primers/digestfragments/reagents',(f,r))
    clean('pcr1_cleanup','Source1.2xAMPureXPcleanup')
    op('sbf1','SourceFIRSTSbfIdigest37C1h',S.DIGESTION,'SourceSbfIR3642','Sourcecut/uncutDNAandallfragments/reagents; background recognitionrequiresduplexsite, no perfectPtclassificationorautomaticphysicalfragmentloss')
    clean('sbf1_cleanup','Source1.2xAMPureXPcleanup')
    op('pcr2','SourceQ5 second6–8PCRcycles',S.PCR,'SourceQ5M0544','SourceamplifiedDNAplusoriginalsubstrates/unamplifiedDNA/primers/remnants',(f,r))
    op('sbf2','Source SECONDSbfIdigest',S.DIGESTION,'SourceSbfI;secondincubationconditionsnotrepeatedinMethods','Sourcecut/uncutDNAanddigestfragments; no postdigestPCRinvented')
    clean('final_cleanup','Sourcefinal1.2xAMPureXPcleanup')
    qc('qubit','SourceQubit4/EqualbitdsDNAHSquantitation;allocationunassigned');qc('tapestation','SourceAgilent4150TapeStationQC;allocationunassigned')
    lib=g.product('library','SourcePtseqDNAsequencinglibrary; fullendpointsunassembled',K.LIBRARY);g.react('identity','Sourcephysicalsequencinginterface',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='PublishedNovaSeq6000/MGISEQ2000/MGISEQT7 PE150; specimen-to-runassignmentandhistoricalprimers/indexlandingunassembled. SbfI5prime-readfilteringoccursafteracquisition, not molecular removal. NN6secondadapterwindowfunctionnotassignedUMI;separateHGPRTvalidationUMIpipelineisnotthisgenomiclibrary.'))
    g.add_source_item(SourceItem('pending','HistoricalE7546/E7595modules/completeadapterligationandprimerlanding/readgeometry;sourcecitedIP/dedicationpredecessors','Main/SI/TableS1selectedscope read;no fullsource/absenceaudit/publicpromotion.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a=ordering_oligo('first_u');b=ordering_oligo('first_b');c('TerminalPSbondkeepslastnucleotide',(a.sequence[-1],a.ps_after[-1]),('T',len(a.sequence)-2));c('Bterminalamineandphosphatepersist',(b.five_prime,b.three_prime),('P','NH2 C7'));c('SecondNNwindowhasunassignedmeaning',ordering_oligo('second_u').index_start,None);c('FirstadapterstemannealsbyScene',bool(first_adapter_scene().rows()),True)
    for bonds in ((-1,),(99,),(0,0)):
        try:OrderedDNA('AC',bonds);bad=False
        except ValueError:bad=True
        c('InvalidPSbondrejected '+str(bonds),bad,True)
    g=protocol_graph();c('Inputcontrolseparatestored',g.terminals['input_control'].disposition,Disposition.STORED);c('Primerdigestremainsuntilsecondligation',g.steps['denature_second'].inputs[0].material,'primer_exo_product');c('SecondPCRfollowedbysecondSbfIbeforecleanup',g.steps['sbf2'].inputs[0].material,'pcr2_product');c('NoPCRafterfinaldigest',g.steps['final_cleanup'].inputs[0].material,'sbf2_product');c('Unknownwholegeometrypreventspromotion',bool(g.publication_problems()),True)
    source=Source(sources_root()/'pt-seq-platinum__10.1126+sciadv.adx6809'/'sciadv.adx6809_supp_sciadv.adx6809_tables_s1_to_s5_adx6809_tables_s1_to_s5.xlsx.txt','get_sources --slugs pt-seq-platinum__10.1126+sciadv.adx6809;TableS1rows5–11')
    if have(c,source,label='PublishedTableS1'):
        text=source.path.read_text();c('WholemodifiedfirstUsourceinventory',a.starred() in text,True);c('WholemodifiedfirstBsourceinventory',b.starred() in text,True);c('Allsevenorderedoligosreadfromsource',all(ordering_oligo(n).starred() in text for n in ('first_u','first_b','extension','second_u','second_b','universal','index')),True)
    c.report()
if __name__=='__main__':selftest()
