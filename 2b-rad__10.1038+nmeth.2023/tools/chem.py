"""Protocol-owned source-qualified incomplete graph; no publication marker."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
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

def ligation_geometry():
    statement='SourceAlfIlongstrands have3primeNNcohesive ends;shortinvertedTcap andcompletegenomicjunction geometry remainunassembled. Known adapter end doesnotassertfullcovalentproduct.'
    ends=[reported_adapter_end('SourceAlfIadapter'+str(n),kind='cohesive',polarity='3-prime',terminal_sequence='NN',body_sequence=body,statement=statement) for n,body in ((1,LONG1),(2,LONG2))]
    return AdapterLigationGeometry.known('cohesive',*(sc for end in ends for sc in end.scenes),scene('Full genomic substrate and inverted-cap junction geometry unassembled'),statement=statement)

from chemdraw import Construct,feature,revcomp
SOURCE='Published2012SI DOI10.1038/nmeth.2023 MOESM591 Table6/8/10(originalp11rendered),DetailedProtocol848–1029,Appendix1031–1067. Definingmainnotobtained;privateincompleteSOLiDAlfIscope.'
LONG1='CTCTCTATGGGCAGTCGGTGAT'
LONG2='CTGCTGTACGGCCAAGGCGAA'
P3='CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT'
P4_PREFIX='CTGCCCCGGGTTCCTCATTCTCT'
P4_TAIL='CTGCTGTACGGCCAAGGCG'
P5='CCACTACGCCTCCGCTTTCCTCTCTATG'
P6='CTGCCCCGGGTTCCTCATTCT'
SOURCE_OLIGO_EXCEPTIONS={'p4':'PublishedTable10originalp11 andAppendix1031–1067 both print exact CTGCCCCGGGTTCCTCATTCTCT-[bc]-CTGCTGTACGGCCAAGGCG family;Table10 explicitly allows6–10nt generally andactualstudy6nt. Thismodelpreserveswholefixedflanks andstructuredunassigned6ntsampleindex,notaconcretemember. Rootindependentlyverifiedwholeoriginaltable/appendix2026-10-29 andapprovedthissingleinputexemption because genericbc macro width is not represented by stricttextmatcher. Fiveotherexactstrandsremainchecked;no sourceabsenceorchemistryexemption.'}
BC=feature('two_b_rad.sample','sample_index','unknown',note='Sourceactualstudy6ntsamplebarcode;memberunassigned,notUMI.')
@dataclass(frozen=True)
class InvertedTCapInventory:
    core:str
    written:str
    table10_written:str
    def __post_init__(self):
        if not self.core or set(self.core)-set('ACGT') or self.written!=self.core+'T(invT)' or '(invT)'+self.core[::-1] not in self.table10_written:raise ValueError('Explicitnormal5to3coreandpostfixannotatedterminalTnotationrequired')
    @property
    def nucleotide_count(self):return len(self.core)+1
    @property
    def capped_nucleotide(self):return 'T'
    @property
    def extendable_3oh(self):return False
    @property
    def bond_geometry(self):return 'Table10selectedterminalinverted-dTnoncanonicalpolarity;exactbondgeometryunassembled'
    @property
    def appendix_equivalence(self):return 'Unreviewed:postfixT(invT)mayannotateterminalTorindicateaddition;fullchemicalequivalence/sourceapplicabilitynotasserted'
    def selected_extension(self,bases):raise ValueError('BlockedterminalinvertedTcannotbeordinary3OHextensionprimer')
SHORT1=InvertedTCapInventory('ATCACCGACTGCCCA','ATCACCGACTGCCCAT(invT)',"3'-(invT)ACCCGTCAGCCACTA-5'")
SHORT2=InvertedTCapInventory('TTCGCCTTGGCCG','TTCGCCTTGGCCGT(invT)',"3'-(invT)GCCGGTTCCGCTT-5'")
@dataclass(frozen=True)
class TwoSidedCutReference:
    top_pieces:tuple
    bottom_pieces:tuple
    coordinates:tuple
    inferred:bool=True
    # Eachstrandis5to3;bottomcoordinatesareinbottom5to3order,notcopiedtopaxis.
    def rejoin(self):return (''.join(self.top_pieces),''.join(self.bottom_pieces))
def selected_two_sided_cut(top,bottom,top_left,top_right,bottom_left,bottom_right,*,inferred=True):
    if not top or set(top)-set('ACGT') or bottom!=revcomp(top):raise ValueError('ExplicitfullycomplementaryDNAstrandsrequired')
    for left,right in ((top_left,top_right),(bottom_left,bottom_right)):
        if type(left)is not int or type(right)is not int or not 0<=left<right<=len(top):raise ValueError('Fourorderedcallerexplicitstrandcutcoordinatesrequired')
    if max(top_left,len(top)-bottom_right)>=min(top_right,len(top)-bottom_left):raise ValueError('Central strands must share a duplex interval to form a selected tag')
    return TwoSidedCutReference((top[:top_left],top[top_left:top_right],top[top_right:]),(bottom[:bottom_left],bottom[bottom_left:bottom_right],bottom[bottom_right:]),(top_left,top_right,bottom_left,bottom_right),inferred)
def oligo(body):return Construct([Segment('source oligo',body)])
def p4_family():return Construct([Segment('source p4 prefix',P4_PREFIX),Segment('source study sample barcode','N'*6,feature=BC),Segment('source p4 tail',P4_TAIL)])
def protocol_graph():
    g=ProtocolGraph('2b-RAD — publishedSOLiDAlfIdetailedprotocol',authority=SOURCE,representation='published_instance')
    x=g.input('genomic','Source1ugintactgenomicDNAin4ul;sourcequality/dry-resuspensionworkupbeforeentry',K.MOLECULAR,primary=True)
    x=op(g,'digest',x,'Source6uldigest:bufferR/SAM/1UAlfI37C1h then65C20min,holdonice',S.DIGESTION,reagent='AlfI/bufferR/SAM;motifoffsets andcut/yielddistributionnotinferred',catalyst=True,product='TypeIIBtwo-sidedselectedtagclassesplusflankingDNA,uncutDNA,enzymesandcompanions; allphysicalDNAretained')
    assay=g.product('digest_qc','Source1ul digestanalyticalgelaliquot;intactcomparisonDNAseparate',K.QC);remain=g.product('digest_remaining','Source5ulremainingdigest;allclasses/remnantsretained',K.MOLECULAR)
    g.react('digest_qc_split','Source1ulQCversus5ulretaineddigest',S.SPLIT,((x,R.TEMPLATE),),(assay,remain));g.mark_terminal(Terminal(assay.id,Disposition.QC,'Sourceagarosecheckisnotpreparativesequencingmaterial.'))
    duplexes=[]
    for number,body,cap in ((1,LONG1,SHORT1),(2,LONG2,SHORT2)):
        a=g.input('adapter'+str(number)+'_long','SourceAlfIlongoligo;3primeNNcohesiveend,notUMI',K.OLIGO,state=oligo(body+'NN'))
        b=g.input('adapter'+str(number)+'_short','Sourceanti-sld-ada'+str(number)+'normalcore+invertedTnucleotidecap;exactnoncanonicalbondgeometryunassembled',K.OLIGO,chemical_state=cap)
        d=g.product('adapter'+str(number),'Source4uMeachpartialduplexAlfIadapterstock;terminalinvertedTretained,fullstrand/overhangpairingunassembled',K.MOLECULAR)
        g.react('prepare_adapter'+str(number),'Sourcecombinelongandcappedshortoligosseparatelyforeachadapterstock',S.HYBRIDIZATION,((a,R.TEMPLATE),(b,R.STRUCTURAL)),(d,));duplexes.append(d)
    ligase=g.input('ligase','SourceATP/T4DNAligase/buffer;20ulmastermix+5ulretaineddigest',K.CATALYST);ligated=g.product('ligated','Source25ulligationmixture:adapter-tagclassesplusflanking/uncutDNA,freeadapters,cappedshortstrandsandreagents;nodepletionassumed',K.MOLECULAR)
    g.react('ligate','SourceAlfIligation16C1h thenholdonice;neverabove65CbeforePCR',S.LIGATION,((remain,R.TEMPLATE),*( (a,R.STRUCTURAL) for a in duplexes),(ligase,R.CATALYST)),(ligated,),geometry=ligation_geometry())
    test=g.product('test_input','Source4ulligationfortestPCR',K.MOLECULAR);production=g.product('prep_input','Source20ulligationforpreparatoryPCR',K.MOLECULAR);rest=g.product('ligation_remainder','Unallocatedphysicalremainderof25ulligation;source laterfateunassigned',K.MOLECULAR)
    g.react('allocate_ligation','Sourceindependent4ultest/20ulprepallocationwithunallocatedremainder',S.SPLIT,((ligated,R.TEMPLATE),),(test,production,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Notinferredstoredordiscarded.'))
    primers=[]
    for key,body in (('p3',P3),('p5',P5),('p6',P6)):
        primers.append(g.input(key,'Sourceappendixsld-'+key[1:]+'usedasdetailedprotocolsld-'+key,K.OLIGO,state=oligo(body)))
    primers.append(g.input('p4','Sourceappendixsld-4/detailedprotocolsld-p4;6ntstudysampleindexmemberunassigned',K.OLIGO,state=p4_family()))
    for key,substrate,label in (('test_pcr',test,'Source20ultestPCR;sample5ulat6/8/10/12cycles for2%gel'),('prep_pcr',production,'Source100ulpreparatoryPCR;QCselectedminimumvisible130bpcyclecount,notuniversally12')):
        enzyme=g.input(key+'_phusion','SourcePhusion/HFbuffer/dNTPs;fourprimersatprintedconcentrations',K.CATALYST);out=g.product(key+'_mixture','SourcePCRcopiesplusoriginalDNA,unextendedprimers,cappedadapterstrands,reagents/remnants;fullcopy/endgeometryunassembled',K.QC if key=='test_pcr' else K.MOLECULAR)
        g.react(key,label,S.PCR,((substrate,R.TEMPLATE),*((a,R.STRUCTURAL) for a in primers),(enzyme,R.CATALYST)),(out,),geometry=PCRGeometry.unavailable(scene('Sourcefourprimer/cappedadapter/duplexproductgeometryunassembled'),reason='SelectedPCRknown;nativeDNA/endmarks arenotassumedfaithfullyretainedonnewcopies.'))
        if key=='test_pcr':g.mark_terminal(Terminal(out.id,Disposition.QC,'Sourceentire20ultestPCRsampled4x5ulforcycle-choicegels;analyticalmaterialnotproduction.'));test_out=out
        else:x=out
    # CyclechoiceisaQCobservation,notphysicalmixingoftestPCRintoprepPCR.
    x=cleanup(g,'gel',x,'SourceentireprepPCRon2%gel;recover~130bpband,excludehighMWand<=90bpdimers')
    x=op(g,'elute',x,'Sourcegelslice+40ulwater,spin2min,4Covernight',reagent='Sourcewaterelution;commercialgelextractionalternativeoutside selectedbranch',product='GelpluselutedDNAandcarryover;elutionnotyetphysicaltransfer')
    kept=g.product('eluate','Source~50ultransferredeluentreadyforsequencing;exactpurityandadaptergeometryunreviewed',K.LIBRARY);gel=g.product('residual_gel','Sourcephysicalgel/retainedmaterialaftereluenttransfer;laterfateunassigned',K.MOLECULAR)
    g.react('collect_eluent','Sourcenext-dayeluenttransferto newtube',S.CLEANUP,((x,R.TEMPLATE),),(kept,gel));g.mark_terminal(Terminal(gel.id,Disposition.UNRESOLVED,'Untransferredgel/retainedDNAarephysicallyseparate;laterfateunassigned.'))
    g.sequence(kept,SequencingPlan(limitation='ActualSOLiDAlfIlibraries verifiedinTable6;exactrun/cycles/indexacquisition/fulllibrarygeometryunreviewed. Computationaltag/allelecountingdoesnotcreateuniformtagsorclearnon-targetDNA. OverhangNNnotUMI;study6ntsamplebarcodeunassigned.'))
    g.add_source_item(SourceItem('pending','Definingmain,historicalSOLiDread/adaptor/capbondgeometryandfullsourcefigures','OnepublishedSIread;Table10p1/p2versusdetailed/appendixp5/p6namesandsetskeptscoped,notdeclaredapublicationerror. BsaXI/Illumina/RTRbranchesnotimported;RTRsource55.6%targetingnotperfectspecificity.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();top='ACGTGCTAGCTAGGCAATGC';bottom=revcomp(top);cut=selected_two_sided_cut(top,bottom,3,12,5,14)
    c('selectedtwosidedcutconservesbothwholestrands',cut.rejoin(),(top,bottom))
    c('fourindependentcutsretained',cut.coordinates,(3,12,5,14))
    c('tagclassdoesnoteraseflankingDNA',len(cut.top_pieces),3)
    try:selected_two_sided_cut(top,bottom,3,2,5,14);bad=False
    except ValueError:bad=True
    c('invalidcallercoordinatesrejected',bad,True)
    c.raises('nonoverlappingcentralstrandsdonotformaduplextag',lambda:selected_two_sided_cut(top,bottom,0,2,0,2))
    c('invertedTcountedasanucleotide',(SHORT1.nucleotide_count,SHORT2.nucleotide_count),(16,14))
    c('rawtable10andappendixformsbothretained',('(invT)ACCCGTCAGCCACTA' in SHORT1.table10_written and SHORT1.written.endswith('T(invT)')),True)
    c('invertedTnotordinaryextendableT',(SHORT1.capped_nucleotide,SHORT1.extendable_3oh),('T',False))
    try:SHORT1.selected_extension('A');bad=False
    except ValueError:bad=True
    c('extensioncannoteraseorignorecap',bad,True)
    g=protocol_graph();c('digestQCnotpreparative',g.terminals['digest_qc'].disposition,Disposition.QC)
    c('testreactionisanalyticalQCkind',g.materials['test_pcr_mixture'].kind,K.QC)
    c('testPCRnotpreparative',{i.material for i in g.steps['prep_pcr'].inputs if i.role==R.TEMPLATE},{'prep_input'})
    c('unallocatedligationfractionunresolved',g.terminals['ligation_remainder'].disposition,Disposition.UNRESOLVED)
    c('fullgeometrypreventspromotion',bool(g.publication_problems()),True)
    source=Source(sources_root()/'2b-rad__10.1038+nmeth.2023'/'nmeth.2023_41592_2012_BFnmeth2023_MOESM591_ESM.pdf.txt','PublishedTable10andAppendix;originaltablep11rendered')
    if have(c,source,label='publishedSI'):
        t=source.path.read_text();c('wholesourceappendixprimersretained',all(x in t for x in (P3,P4_PREFIX+'[bc]'+P4_TAIL,P5,P6)),True)
        c('bothwholewritteninvertedTnotationsretained',SHORT1.written in t and SHORT2.written in t,True)
    c.report()
if __name__=='__main__':selftest()
