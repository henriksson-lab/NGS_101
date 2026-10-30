"""Protocol-owned source-qualified incomplete graph; no publication marker."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

def scene(label):
    sc=Scene();sc.strand('preparation',[Segment('population','['+label+']',placeholder=True)]);return sc

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; geometry unassembled',product=None,state=None,catalyst=False):
    stock=g.input(key+'_stock',reagent,K.CATALYST if catalyst else K.CONDITION)
    out=g.product(key+'_product',product or label+'; molecular mixture/remnants unassembled',K.MOLECULAR,state=state)
    kw={'reason':'Source operation established; exact molecular ends, companion species and geometry unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION)),(out,),**kw);return out

def cleanup(g,key,x,label):
    stock=g.input(key+'_condition',label,K.CONDITION)
    out=g.product(key+'_kept','Source recovered preparation; recovery and carryover unquantified',K.MOLECULAR)
    waste=g.product(key+'_waste','Actual source-separated workup fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(stock,R.CONDITION)),(out,waste));g.discard(waste,'Source physical cleanup, no perfect partition implied.');return out

def uncertain_ligation():return AdapterLigationGeometry('unavailable','unavailable',(scene('unreviewed exact ends and junction geometry'),),'Known source ligation, complete molecular geometry not assembled.')



from chemdraw import revcomp
from illumina import P5,TRUSEQ_READ1
SOURCE='Published2010GenomeBiologyPMC2884539 Methods103–111/Results37–51/Fig1;publishedAdditionalfile1TableS9 originalp218rendered/source107474–107501. Historicalkit/fullgenomicjunction/RTgeometryunassembled.'
T7='cagtaatacgactcactatagggagaaggct'
PE='CAAGCAGAAGACGGCATACGACAGCAG'
AE1='A'+T7.upper()+PE
AE2='CG'+revcomp(AE1[1:])
AS1=TRUSEQ_READ1
AS2=revcomp(AS1[:-1])
PS=P5+TRUSEQ_READ1[4:]
@dataclass(frozen=True)
class OligoInventory:
    sequence:str
    five_prime:str='unassigned'
    ps_after:tuple=()
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGT') or self.five_prime not in {'P','unassigned'} or any(type(i)is not int or not 0<=i<len(self.sequence)-1 for i in self.ps_after):raise ValueError('ActualorderedDNA/endandPSbondpositionsrequired')
    def scene(self,name):
        sc=Scene();sc.strand(name,[Segment(name,self.sequence)],mod5='P' if self.five_prime=='P' else '',mod3='sourcePSbondbeforeterminalT' if self.ps_after else '');return sc

def source_oligo(name):
    return {'ae1':OligoInventory(AE1),'ae2':OligoInventory(AE2,'P'),'as1':OligoInventory(AS1,ps_after=(len(AS1)-2,)),'as2':OligoInventory(AS2,'P'),'ps':OligoInventory(PS),'pe':OligoInventory(PE),'pseq':OligoInventory(TRUSEQ_READ1)}[name]
@dataclass(frozen=True)
class SelectedIVTReference:
    template_5to3:str
    begin:int
    end:int
    inferred:bool=True
    def __post_init__(self):
        if not self.template_5to3 or set(self.template_5to3)-set('ACGT') or any(type(i)is not int for i in (self.begin,self.end)) or not 0<=self.begin<self.end<=len(self.template_5to3):raise ValueError('CallerexplicitDNATEMPLATEandtranscribedintervalrequired;noT7startdefaults')
    @property
    def new_rna(self):return revcomp(self.template_5to3[self.begin:self.end]).replace('T','U')
    def selected_first_cdna(self,primer,covered):
        if not primer or set(primer)-set('ACGT') or type(covered)is not int or not 0<covered<=len(self.new_rna) or not primer.endswith(revcomp(self.new_rna[-covered:].replace('U','T'))):raise ValueError('CallerexplicitvalidRNAprimerregisterrequired,notanoperationalRTassignment')
        return primer+revcomp(self.new_rna[:-covered].replace('U','T'))

def adapter(g,name):
    names=('ae1','ae2') if name=='ae' else ('as1','as2');xs=[]
    for n in names:
        inv=source_oligo(n);xs.append(g.input(n,'PublishedTableS9 '+n+' actualDNA/end/modification inventory',K.OLIGO,state=inv.scene(n),chemical_state=inv))
    a=g.product(name+'_adapter','Sourcecustom'+name.upper()+'two-oligoadapterinventory;annealingprogramme/fullmoleculargeometryunassembled',K.ADAPTER);cond=g.input(name+'_assembly','Sourceadapterduplexpreparation; exacttemperatureprogrammeunreviewed',K.CONDITION);g.react(name+'_assemble','Sourcecustomtwo-oligoadapterpreparation',S.HYBRIDIZATION,((xs[0],R.TEMPLATE),(xs[1],R.STRUCTURAL),(cond,R.CONDITION)),(a,));return a

def graph_for_channel(channel='hpaii'):
    if channel not in {'hpaii','mspi'}:raise ValueError('SourceindependentHpaIIorMspIchannelrequired')
    g=ProtocolGraph('HELP-tagging — '+channel+'sourcechannel',authority=SOURCE,representation='published_instance')
    x=g.input('genomic','SourceH1WA01genomicDNA;unboundednativeCpG/methylationstatesnotassigned',K.MOLECULAR,primary=True)
    x=op(g,'digest',x,'Source5ugDNA/'+('HpaII' if channel=='hpaii' else 'MspI')+'in200ul;time/bufferunreviewed',S.DIGESTION,reagent='Source'+('HpaII' if channel=='hpaii' else 'MspI'),catalyst=True,product='CutanduncutnativeDNA/fragments/nucleasemixture;nocytosineconversionorfaultlessspecificityassumed')
    x=cleanup(g,'digest_workup',x,'Sourcephenol/chloroform extraction and ethanol precipitation')
    used=g.product('ligation_input','Source'+('3ugHpaII' if channel=='hpaii' else '1ugMspI')+'digestallocationtoAEligation',K.MOLECULAR);rest=g.product('digest_remainder','Independentunuseddigest;mass/fateunassigned',K.MOLECULAR);g.react('allocate','SourceactualAEligationmassallocation',S.SPLIT,((x,R.TEMPLATE),),(used,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Noimplicitlydiscarded/storedDNA.'))
    ae=adapter(g,'ae');cat=g.input('quick_ae','SourceQuickLigase3ul/2xQuickbuffer',K.CATALYST);y=g.product('ae_ligated','AE-genomicclasseswithsourceT7/EcoP15Isequenceplusfreeadapters/otherDNA/ligationremnants;wholejunctionunassembled',K.MOLECULAR);g.react('ae_ligate','SourceAEQuickligation50ul;temperature/timeunreviewed',S.LIGATION,((used,R.TEMPLATE),(ae,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=uncertain_ligation())
    x=cleanup(g,'ae_cleanup',y,'SourceAMPurecleanup')
    x=op(g,'ecop15i',x,'SourceEcoP15Idigestion',S.DIGESTION,reagent='SourceEcoP15I',catalyst=True,product='SelectedAE-bearingflankingtagDNAplusotherEcoP15Idaughters/uncutDNAandadapters/remnants;nosilentdistaldiscard')
    x=op(g,'repair',x,'Sourcefragmentend-repair,sourceantiadapterdimerizationpurpose',S.END_REPAIR,reagent='Sourceend-repaircatalyticcomponentidentity/geometryunreviewed',catalyst=True,product='SelectedrepairedDNAendsandallDNA/offproducts/remnants;exactrepairchemistryunassembled')
    x=op(g,'a_tail',x,'Sourcesingle3dAtailing',S.TAILING,reagent='SourcedA-tailingcatalyticcomponentidentityunreviewed',catalyst=True,product='Selected3dA-endedDNAplusunprocessedDNA/remnants')
    a=adapter(g,'as');cat=g.input('quick_as','SourceQuickLigase2.5ul/2xQuickbuffer',K.CATALYST);y=g.product('as_ligated','AE-insert-AS andAS-insert-AS/otherDNAclasseswithactualsourcePSterminalT/phosphatesplusfreeadapters/reagents',K.MOLECULAR);g.react('as_ligate','SourceASQuickligation50ul',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=uncertain_ligation())
    x=cleanup(g,'minelute',y,'SourceMinElutePCRpurification')
    x=op(g,'ivt',x,'SourceAmbionMEGAshortscriptT7IVT',S.TRANSCRIPTION,reagent='SourceMEGAshortscriptT7RNApolymerase/NTPmixture;start/terminationgeometryunassembled',catalyst=True,product='SelectedAE-promoter-derivedNEWRNAplusoriginalnative/copiedDNA,untranscribedAS-ASclasses/reagents;selectionnotphysicalinstantdiscard')
    x=cleanup(g,'rneasy',x,'SourceRNeasyRNAcleanup;exactDNA/RNA/remnantpartitionunassembled')
    primer=g.input('rt_primer','SourceSuperScriptIIIRTprimer;exactassignmentunreviewed,TableS9doesnotbyitselfestablishlanding',K.OLIGO);cat=g.input('ssiii','SourceSuperScriptIIIkit',K.CATALYST);y=g.product('cdna','SourcefirststrandDNAcopiesplusoriginalIVTRNA/nativeDNAcarryover/unextendedprimer/reagents;noinventedRNAhydrolysis',K.MOLECULAR);g.react('rt','SourceSuperScriptIIIfirststrandRT',S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(cat,R.CATALYST)),(y,))
    ps=g.input('primer_ps','SourceTableS9PSPCRprimer',K.OLIGO,state=source_oligo('ps').scene('ps'));pe=g.input('primer_pe','SourceTableS9PEPCRprimer',K.OLIGO,state=source_oligo('pe').scene('pe'));cat=g.input('pcr_enzyme','SourcePCRcatalyticcomponentidentityunreviewed',K.CATALYST);z=g.product('pcr_product','NewDNAcopiesplusoriginalcDNA/RNA/DNAcarryover/freeprimers and reactionremnants',K.MOLECULAR);g.react('pcr','Source96C2min,18cycles96C15s/72C15s,72C5min',S.PCR,((y,R.TEMPLATE),(ps,R.STRUCTURAL),(pe,R.STRUCTURAL),(cat,R.CATALYST)),(z,),geometry=PCRGeometry.unavailable(scene('SourcefullAE/insert/AS/RNA/RT/PCRbackbonegeometryunassembled'),reason='ExactoligoinventorydoesnotassignwholeoperationalRT/PCRgeometry.'))
    x=cleanup(g,'qiaquick',z,'SourceQIAquickPCRcleanup')
    lib=g.product('library','SourceHELP-taggingDNAsequencinglibrary;wholegeometryunassembled',K.LIBRARY);g.react('identity','Sourcephysicalfinallibrary',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='PublishedResultsGenomeAnalyzer36bpSE. TableS9PSeqexactinventoryknown,butbinding/readgeometryunassembled. EcoP15Iflankingtags27ntsourceexpectation,notinferredfixedgenomiccoordinates;HpaII/MspIcountnormalizationandtagfilteringcomputational. NoUMI/indexassignmentorperfectIVTselectionassumed.'))
    g.add_source_item(SourceItem('pending','Historicalkit/end-repair/T7start/RT/fullPCR/readgeometry andotherSIsections','Definingpublishedmain+SIoligotable read;massivetabledata/otherprimer-validationinventoriesnotcompleteabsenceaudit. NoMSCCrecipeimport/default.'))
    return g

def protocol_graph():return graph_for_channel()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a=source_oligo('as1');c('TerminalPSbondretainsactualTnucleotide',(a.sequence[-1],a.ps_after),('T',(len(a.sequence)-2,)))
    c('AEsource5CGoverhangandcomplementbody',AE2, 'CG'+revcomp(AE1[1:]))
    r=SelectedIVTReference('ACGTTCAG',1,7);c('IVTcreatesRNAcomplementwithoutrewritingDNA',(r.new_rna,r.template_5to3),('UGAACG','ACGTTCAG'))
    c('ExplicitRTregisterextendsrealprimeronly',r.selected_first_cdna('CG',2),'CGTTCA')
    try:r.selected_first_cdna('AA',2);bad=False
    except ValueError:bad=True
    c('InvalidRTregisterrejected',bad,True)
    g=protocol_graph();h=graph_for_channel('mspi');c('Unuseddigestremainsphysical',g.terminals['digest_remainder'].disposition,Disposition.UNRESOLVED)
    c('IVTRNAiscreatedfromretainedDNAproduct',g.steps['ivt'].inputs[0].material,'minelute_kept')
    c('NoRNAhydrolysisinvented','rnase_h' not in g.steps,True)
    c('BothsourcechannelsretainIVTselection','ivt' in h.steps,True)
    c('Sourcegeometrylimitpreventsfinishedmarker',bool(g.publication_problems()),True)
    src=Source(sources_root()/'help-tagging__10.1186+gb-2010-11-4-r36'/'gb-2010-11-4-r36_supp_gb-2010-11-4-r36-S1.PDF.txt','get_sources --slugs help-tagging__10.1186+gb-2010-11-4-r36;TableS9originalp218rendered')
    if have(c,src,label='PublishedTableS9'):
        text=''.join(src.path.read_text().split());c('WholemodifiedAS1sourceverified',AS1[:-1]+'-x-T' in text,True);c('WholeAE1caseandAE2modificationpreserved',('A'+T7+PE) in text and ('p-'+AE2[:2]+AE2[2:-len(T7)]+T7.translate(str.maketrans('acgt','tgca'))[::-1]) in text,True)
    c.report()
if __name__=='__main__':selftest()
