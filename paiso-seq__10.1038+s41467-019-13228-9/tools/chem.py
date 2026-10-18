"""PAIso-seq2019 source-owned mixed RNA/DNA extension, ssPCR and nativeSMRTbell boundary."""
from chemdraw import Construct,Scene,Segment,MolecularState,feature,revcomp,complement_segments,circle_rows
from protocol_model import role
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan,SourceClaim,SourceDiscrepancy
import seqprimers as sp
TITLE="PAIso-seq — poly(A)-inclusive full-length isoform sequencing"
NOTES="01_paiso-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-019-13228-9">Liu et al., <i>Nature Communications</i> (2019)</a>; detailed protocol <a href="https://doi.org/10.1038/s41596-022-00704-8">Nat. Protocols (2022)</a>.'
SUMMARY="A dU-interrupted oligo templates end extension across the complete RNA poly(A) tail, USER releases the copied product, template switching captures the 5-prime end, and PacBio CCS reads the entire isoform plus native tail."
CAVEAT=('The 2019 paper (Methods, “PAIso-seq library construction”) names a PacBio '
        'SMRTbell Template Prep Kit without its version. The later <a href="https://doi.org/10.1038/s41596-022-00704-8">'
        '2022 author protocol</a> specifies kit 1.0-SPv3 (Materials and step 29), '
        'but that does not establish the 2019 experiment\'s kit. The <a href="https://www.pacb.com/wp-content/uploads/2015/09/Insert-SMRTbell-Template-Prep-Kit-1.0.pdf">'
        'PacBio 1.0 insert, p. 1</a> lists a Blunt Adapter, while the <a href="https://www.pacb.com/wp-content/uploads/SMRTbell-Express-TPK-2.0-Overview-Large-Insert-gDNA-Library-Prep-Customer-Training.pdf">'
        'Express 2.0 overview, pp. 6 and 8</a> specifies a T-overhang adapter and A/T ligation. '
        'The model preserves the closed SMRTbell topology but leaves this experiment\'s '
        'end-preparation and insert–adapter junction geometry unresolved. Native hairpin roles '
        'carry no assigned bases or molecular length. '
        'The later protocol also adds 5′ iCiGiC and 3′ amine blocking to its TSO '
        '(Table 2); this model uses the 2019 Supplementary Table 1 TSO instead.')

PCR_HANDLE="AAGCAGTGGTATCAACGCAGAGT"
GUIDE_LINKER="AC"
TSO_LINKER="ACAT"
TSO_ORDERING=PCR_HANDLE+TSO_LINKER+"rGrG+G"
RT_PRIMER=PCR_HANDLE+GUIDE_LINKER
PCR_PRIMER=PCR_HANDLE
BARCODES = (
    "TCAGACGATGCGTCAT", "CTATACATGACTCTGC", "TACTAGAGTAGCACTC",
    "TGTGTATCAGTACATG", "GATCTCTACTATATGC", "ACAGTCTATACTGCTG",
    "ATGATGTGCTACATCT", "CTGCGTGCTCTACGAC", "GCGCGATACGATGACT",
    "CGCGCTCAGCTGATCG", "GCGCACGCACTACAGA", "ACACTGACGTCGCGAC",
    "CGTCTATATACGTATA", "ATAGAGACTCAGAGCT", "TAGATGCGAGAGTAGA",
)
DU_TAIL="dUTTTTTTTTdUTTTTTTTT"
SAMPLE_BARCODE=feature("sample_barcode","sample_index","whitelist",whitelist="original2019TableS1 2dU-BC1..15")
REPRESENTATIVE_BARCODE=BARCODES[2]
def seg(name,seq,tag=None,**kw):return Segment(name,seq,tag,**kw)
def extension_oligo(sample):
 if not 1<=sample<=15:raise ValueError('sourceBCmember mustbe1..15')
 return RT_PRIMER+BARCODES[sample-1]+DU_TAIL

def extension_oligo_state(sample=3):
 return Construct([seg('RTprimerDNA',RT_PRIMER),seg('BC'+str(sample),BARCODES[sample-1],'cbc',feature=SAMPLE_BARCODE),seg('firstdU','U',note='DNAdeoxyuridine; USERsite'),seg('firstT8','T'*8),seg('seconddU','U',note='DNAdeoxyuridine; USERsite'),seg('secondT8','T'*8)],name='source2dUguideBC'+str(sample))
def tso_state():return Construct([seg('sourceTSODNAprefix',PCR_HANDLE+TSO_LINKER),seg('rG1','G',note='physicalRNA G'),seg('rG2','G',note='physicalRNA G'),seg('LNA-G','G',note='physical+G lockednucleicacid')],name='original2019TSO')

def first_parts(switched=False):
 parts=[seg('sourceRTprimer',RT_PRIMER),seg('copiedBC3',REPRESENTATIVE_BARCODE,'cbc',feature=SAMPLE_BARCODE),role('DNAcopyofnativeRNApolyAtail; nonAcontent/lengthretained'),role('antisense transcriptDNA')]
 if switched:parts += [role('native switchjunction/terminalC extentunassigned'),seg('copiedTSODNAprefix',revcomp(PCR_HANDLE+TSO_LINKER))]
 return parts
INSERT=Construct([seg('5-prime PCR handle',PCR_HANDLE),seg('TSO linker',TSO_LINKER),role('DNAcopyofnative switchjunction; Gextentunassigned'),role('full-length transcriptDNAcopy'),role('DNAcopyofnativeRNApolyAtail; variablelength/nonA'),seg('BC3 reverse complement',revcomp(REPRESENTATIVE_BARCODE),'cbc',feature=SAMPLE_BARCODE),seg('guide linker reverse complement',revcomp(GUIDE_LINKER)),seg('3-prime PCR handle reverse complement',revcomp(PCR_HANDLE))],name='PAIso2019BC3 PCRdaughter')
FINAL_LIBRARY=None
SEQ_PRIMERS=()
SHOW_READ_CONTENT=False
SEQUENCING_ENDING='Source2019SequelSMRTbell CCS; loadedprimer/hairpin/kitversion unassigned. CCS≥10passes is computationalselection, not physical purification.'
READ_SPECS=(sp.read_spec('Sequencing output','unassigned','Full cDNA and nativepolyAtail viaPacBioCCS;16ntsamplebarcode; no invented indexcycles','include'),)

def native_state(stage):
 sc=Scene()
 if stage in ('rna','extended','user'):
  parts=[role('original RNAtranscript; physicalU'),role('original RNApolyAtail; variablelength and nonAcontent')]
  if stage!='rna':parts += [seg('DNAcopiedBC3',revcomp(REPRESENTATIVE_BARCODE),'cbc',feature=SAMPLE_BARCODE),seg('DNAcopiedRTguide',revcomp(RT_PRIMER))]
  sc.strand('originalRNA',parts,label='originalRNA'+(' with dNTP-derived DNA3′extension' if stage!='rna' else ''))
  if stage=='extended':sc.strand('guide',list(extension_oligo_state()),label='independent DNAguide; bothdUphysical, pairingextentunassigned')
  if stage=='user':sc.strand('guidefragments',[role('USER-cleaved dUguideDNA fragments; chemistrynotphysicalremoval')])
 elif stage in ('first','switched'):
  sc.strand('cDNA',first_parts(stage=='switched'),label='actual '+('switched' if stage=='switched' else 'first')+' ssDNA')
  sc.strand('originalRNA',[role('independent originalRNA–DNAextendedtemplates; postRTchemicalfateunassigned')])
 elif stage in ('amplified','selected'):
  parts=list(INSERT);sc.strand('daughter',parts,label='PCRdaughter DNA; originalRNA/dU/TSOrG/LNA modifications not copied')
  sc.anneal('daughteropposite',complement_segments(parts),to='daughter',pair=(parts[0].name+"'",parts[0].name),unpaired=tuple(p.name+"'" for p in parts if p.is_role_token()))
 elif stage=='smrtbell':
  parts=[*list(INSERT),role('nativehairpinA; bases/length/ligationendsunassigned'),*complement_segments(list(INSERT)),role('nativehairpinB; bases/length/ligationendsunassigned')]
  rows=list(circle_rows(Construct(parts,name='sourceclosedSMRTbell; topologyonly'),closure_label='two nativehairpins join oppositeinsertstrands; no fixedhairpinlength'))
  residual=Scene();residual.strand('residuals',[role('possible originalDNA/RNA/guide/TSO/PCR and unclosedclasses; kit-specific recoveryunassigned')]);rows.extend(residual.rows())
  return MolecularState('nativeSMRTbell closedcircuit; hairpinendgeometry unassigned',tuple(rows))
 sc.strand('otherlineages',[role('independent otherRNA species / possible guidefragments, residualoligos and originalcDNA; classrecoveryunassigned')])
 return sc

def oligos():
 return [f'TSO—5′-{TSO_ORDERING}-3′',f'RTprimer—5′-{RT_PRIMER}-3′',f'PCRprimer—5′-{PCR_PRIMER}-3′',*[f'2dU-BC{i}—5′-{extension_oligo(i)}-3′' for i in range(1,16)]]

def protocol_graph():
 graph=ProtocolGraph(TITLE,authority=SOURCE)
 def inp(i,l,k=M.CATALYST,state=None):return graph.input(i,l,k,state=state)
 def prod(i,l,stage,kind=M.MOLECULAR):return graph.product(i,l,kind,state=native_state(stage))
 def clean(i,l,x,state=None):
  sel=inp(i+'_selection',l,M.SELECTION);y=graph.product(i+'_out',l+'retained',x.kind,state=state or x.state);w=graph.product(i+'_waste','physicalnonretainedfraction',M.WASTE)
  graph.react(i,l,K.CLEANUP,((x,I.TEMPLATE),(sel,I.SELECTION)),(y,w));graph.discard(w,'Sourcephysical nonretainedfraction; complete classpuritynotasserted.');return y
 rna=graph.input('rna','selectedbulkGV-oocyte totalRNA; BC3guide branch',M.BIOLOGICAL,state=native_state('rna'),primary=True)
 guide=inp('extension_oligo','originalTableS1BC3DNAguide with2dU/T8spacing',M.OLIGO,extension_oligo_state())
 kl=inp('extension_polymerase','Klenowexo− plusdNTP;37°C1h')
 extended=prod('extended','RNA–DNA chimeric3′extension plus originalDNAguide','extended')
 graph.react('extend_tail','Template originalRNA3′end extension with dNTP-derivedDNA',K.EXTENSION,((rna,I.TEMPLATE),(guide,I.STRUCTURAL),(kl,I.CATALYST)),(extended,))
 user=inp('user','USER37°C30min');digested=prod('digested','endextendedRNA plus cleavedguideDNA fragments','user')
 graph.react('open_du','Cleave guideDNA at two physicaldU sites',K.DIGESTION,((extended,I.TEMPLATE),(user,I.CATALYST)),(digested,))
 released=clean('clean_rna','RNA Clean &Concentrator5;>200ntselection;6–8µlwater',digested)
 rp=inp('rt_primer','TableS1RTprimer',M.OLIGO,Construct([seg('RTprimer',RT_PRIMER)]));rtase=inp('reverse_transcriptase','template-switchingRT; source sayssimilarlySmartseq2; exactpredecessorconditionsunread')
 first=prod('first','actualtail-inclusive firstssDNA plusoriginalRNA','first');graph.react('rt','Reverse-transcribe endextendedRNA withsourceRTprimer',K.REVERSE_TRANSCRIPTION,((released,I.TEMPLATE),(rp,I.STRUCTURAL),(rtase,I.CATALYST)),(first,))
 tso=inp('tso','original2019TSO physicalrGrG+G',M.OLIGO,tso_state());switched=prod('switched','actualswitchedssDNA withnativejunctionbounded','switched')
 graph.react('switch','Template-switch attranscript5′end',K.TEMPLATE_SWITCHING,((first,I.TEMPLATE),(tso,I.STRUCTURAL),(rtase,I.CATALYST)),(switched,))
 pp=inp('pcr_primer','originalTableS1PCRsingleprimer species',M.OLIGO,Construct([seg('PCRprimer',PCR_PRIMER)]));pol=inp('pcr_polymerase','KAPAHiFiHotStartReadyMix; nativecyclesunassigned')
 amp=prod('amplified','tail-inclusivePCRdaughters plusresidualtemplates','amplified')
 graph.react('pcr','Amplify actualnative switchedssDNA',K.PCR,((switched,I.TEMPLATE),(pp,I.STRUCTURAL),(pol,I.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(switched.state,reason='ActualssDNA withsourceprimerhandles, originalRNA andresidualguides entersPCR; switchjunction/firstcycletopology not inventedfromfinalamplicon.',symmetric_species=True))
 a=prod('allocated_200','portionallocated to>200bpselection','amplified');b=prod('allocated_2000','independentportionallocated to>2kbselection','amplified')
 graph.react('allocate','Allocate source two portions forindependentbeadselections',K.SPLIT,((amp,I.TEMPLATE),),(a,b),changes_molecule=False)
 a=clean('select_200','PurePB1×;>200bp fraction',a);b=clean('select_2000','PurePB0.4×;>2kb fraction',b)
 molarity=inp("equimolar_rule","sourceequalmolarcontributions; takevolumesunassigned",M.CONDITION)
 chosen=[]
 for name,x in (('shortcutoff',a),('longcutoff',b)):
  take=prod(name+'_equimolar','equimolarcontribution from'+name,'selected');left=graph.product(name+'_remainder','unallocatedfraction; physicalfateunassigned',M.MOLECULAR,state=x.state)
  graph.react(name+'_sample','Takeequimolarcontribution from'+name,K.SAMPLING,((x,I.TEMPLATE),(molarity,I.CONDITION)),(take,left),changes_molecule=False)
  from protocol_graph import Terminal,Disposition
  graph.mark_terminal(Terminal(left.id,Disposition.UNRESOLVED,'Not assignedtopooledSMRTbellinput; no physicaldisposalasserted.'));chosen.append(take)
 pool=prod('selected','equimolarpool of independentlyselectedcDNApopulations','selected')
 graph.react('pool','Recombine two selectedpopulations at equalmolarity',K.POOL,tuple((x,I.TEMPLATE) for x in chosen),(pool,),changes_molecule=False)
 kit=inp('native_kit','unversioned2019SMRTbellTemplatePrepKit; vendorversionunassigned',M.CONDITION);lib=prod('library','sourceSMRTbell closedtemplate pluspossibleoriginal/unclosedclasses','smrtbell',M.LIBRARY)
 graph.react('native_smrtbell','Native unversionedSMRTbell preparationboundary',K.MAGIC,((pool,I.TEMPLATE),(kit,I.CONDITION)),(lib,),reason='2019Methods namesSMRTbellkit andclosedCCS output, not nativeversion/repair/nuclease/hairpinendorder. Later2022/kit1.0/Express2.0 comparators remainseparate; no canonicalA/T/blunt/hairpinwidthdefaults.')
 anneal=inp('sequencing_primer','source sequencingprimer; sequence/memberunassigned',M.OLIGO);ann=prod('annealed_library','sequencingprimer-annealedSMRTbell','smrtbell',M.LIBRARY)
 graph.react('anneal_sequencing','Anneal source sequencingprimer',K.HYBRIDIZATION,((lib,I.TEMPLATE),(anneal,I.STRUCTURAL)),(ann,))
 polymerase=inp('sequencing_polymerase','sourcePacBio polymerase',M.MOLECULAR);loaded=prod('polymerase_loaded','polymerase-boundSMRTbell','smrtbell',M.LIBRARY)
 graph.react('bind_polymerase','Bind source sequencingpolymerase',K.BINDING,((ann,I.TEMPLATE),(polymerase,I.STRUCTURAL)),(loaded,))
 mag=inp('magbeads','sourceMagbeads',M.SELECTION);final=prod('sequencing_library','Magbead-boundPacBioSequeltemplate','smrtbell',M.LIBRARY)
 graph.react('load_magbeads','Bind polymerase–template toMagbeads',K.BINDING,((loaded,I.TEMPLATE),(mag,I.SELECTION)),(final,))
 graph.sequence(final,SequencingPlan(acquired_roles=('Sequencing output',),limitation=SEQUENCING_ENDING))
 graph.declare_reviewed_scope(('equimolar_rule','rna','extension_oligo','extension_polymerase','extend_tail','user','open_du','clean_rna_selection','clean_rna','rt_primer','reverse_transcriptase','rt','tso','switch','pcr_primer','pcr_polymerase','pcr','allocate','select_200_selection','select_200','select_2000_selection','select_2000','shortcutoff_sample','longcutoff_sample','pool','native_kit','native_smrtbell','sequencing_primer','anneal_sequencing','sequencing_polymerase','bind_polymerase','magbeads','load_magbeads','sequencing_library'),SOURCE)
 graph.annotate(SourceDiscrepancy("extension_oligo","Guide-oligo dU/T spacing differs between Results and Supplementary Table 1",
     (SourceClaim("Results, first PAIso-seq paragraph (PMC6876564 JATS text line 54; article PDF text line 158)","3′ end 5′-dUTTTTTTTdUTTTTTTT-3′ (seven T after each dU)"),
         SourceClaim("Supplementary Table 1, rows 2dU-BC1–2dU-BC15","3′ end dUTTTTTTTTdUTTTTTTTT (eight T after each dU)")),
        "The drawn ordering strings use the Supplementary Table 1 oligos.",
        "The table is the base-level ordering list for every barcode member; the Results value is an in-text summary.",
        "The inspected 2019 article and supplement do not reconcile the two values; the later 2022 protocol was not used to settle it."))
 graph.annotate(SourceDiscrepancy("rt_primer","Primer handles described as 'TSO minus triple G'",
     (SourceClaim("Results, first PAIso-seq paragraph (JATS text line 54)","guide, RT and PCR primers carry the TSO sequence minus the triple G"),
         SourceClaim("Supplementary Table 1, rows TSO / RT primer / PCR primer","TSO …GAGTACATrGrG+G; RT primer ends …GAGTAC; PCR primer ends …GAGT")),
        "The exact Supplementary Table 1 primers are drawn.",
        "Removing rGrG+G from the TSO leaves …GAGTACAT, which neither printed primer matches.",
        "The Results wording is likely descriptive; no source states why the RT and PCR primers stop short of ACAT."))
 return graph.require_publishable()
