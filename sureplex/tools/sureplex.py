"""SurePlex15053626v02 and VeriSeq15052877v05 native molecular routes."""
from chemdraw import Scene, Segment, feature
from protocol_model import molecule,role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
import seqprimers as sp
TITLE='SurePlex with VeriSeq PGS library preparation'
NOTES='01_sureplex.html'
SOURCE='Illumina <a href="https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/veriseq-pgs/sureplex-protocol-reference-guide-15053626-02.pdf">SurePlex15053626v02, September2020</a>; <a href="https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/veriseq-pgs/veriseq-pgs-library-prep-reference-guide-15052877-05.pdf">VeriSeqPGS15052877v05, August2021</a>; <a href="https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/veriseq-pgs/veriseq-pgs-guidance-note-1000000138837-00.pdf">GuidanceNote1000000138837v00</a>.'
SUMMARY='Two commercial amplification stages precede native tagmentation/NT neutralization and indexPCR. Bead-normalized products are single-stranded before LNS1 transfer, pooling and HT1 preparation.'
CAVEAT='Source primer pools and native library ends are unassigned; later canonical Nextera strands and index widths are not substituted. No-cell controls are parallel reactions. Nonproductive end classes, original templates and possible artifacts are molecular populations with unresolved recovery; unmapped-read filtering is computational.'
FINAL_LIBRARY=None
SEQ_PRIMERS=()
SEQUENCING_ENDING='GuidanceNote describes the36-cycle short-read assay with single- or dual-index options. This page selects the namedS503/S504 plusN701–N712 dual-index route; native primer sites and acquired index cycle counts remain unassigned.'
SOURCE_LIMITATION='Final library is source single-stranded normalization product. Exact native tagged-end/PCR/sequence-primer assembly is unassigned; no completed duplex is supplied as its physical state.'
SEQUENCING_UNAVAILABLE=''  # No duplicate duplex endpoint: final ssDNA is drawn by the graph.
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec('Read 1','36-cycle assay','SourceGuidanceNote fast short-read recipe; sample-derived genomic observations','include'),sp.read_spec('Index 1 (i7)','unassigned','NamedN701–N712 family; selectedmember/nativewidth/runmask unassigned','demultiplex'),sp.read_spec('Index 2 (i5)','unassigned','NamedS503/S504 family; selectedmember/nativewidth/runmask unassigned','demultiplex'))
I5=feature('sureplex_i5','sample_index','unknown',note='S503/S504 sourcefamily; nativewidth/member unassigned')
I7=feature('sureplex_i7','sample_index','unknown',note='N701–N712 sourcefamily; nativewidth/member unassigned')

def ss(label,parts):return molecule(label,tuple(parts))
def native(stage,indexed=False,single=False):
 sc=Scene()
 parts=[role('native left boundary; exactgeometry unassigned',molecular_feature=I5 if indexed else None),role('SurePlex-derived DNA fragment: '+stage),role('native opposite boundary; exactgeometry unassigned',molecular_feature=I7 if indexed else None)]
 sc.strand('strandA',parts,label='normalized ssDNA; nativeorientationunassigned' if single else 'nativeDNA strandA')
 if not single:sc.strand('strandB',[role('opposite native fragmentstrand: '+stage)],label='separate oppositeDNA strand; pairing/gap/extent unassigned')
 return sc

def cleanup(g,key,current,label,selection=False):
 beads=g.input(key+'_reagent',label,M.SELECTION)
 kept=g.product(key+'_retained','source recovered '+current.label,M.MOLECULAR,state=current.state)
 waste=g.product(key+'_waste','physically nonretained fraction and washes; exactcomposition unassigned',M.WASTE)
 g.react(key,label,K.SELECTION if selection else K.CLEANUP,((current,I.TEMPLATE),(beads,I.SELECTION)),(kept,waste),changes_molecule=False)
 g.discard(waste,'Physically removed source fraction; no species-specific complete depletion claimed.')
 return kept

def commercial_stages(g,key,template):
 pool=g.input(key+'_preamp_pool','source SurePlex Pre-amp Buffer; oligo components and their location unassigned',M.CONDITION)
 enzyme=(g.materials['sample_preamp_enzyme'] if key=='control' else g.input(key+'_preamp_enzyme','SurePlexPreAmpEnzyme stock; shared reagent source across parallel reactions',M.CATALYST))
 pre=g.product(key+'_preamp','12-cycle '+key+' proprietarypreamp material',M.MOLECULAR,state=ss(key+' preamp',(role('possible no-cell amplificationproducts; identity unassigned' if key=='control' else 'genomicDNA-derived native12cycle preamp molecules; ends unassigned'),)))
 g.react(key+'_preamp_stage','NativeSurePlex '+key+' preamplification12cycles',K.MAGIC,((template,I.TEMPLATE),(pool,I.CONDITION),(enzyme,I.CATALYST)),(pre,),reason='Source supplies12-cycle program and named commercialbuffer/enzyme; native primerpool species, first-copy chemistry and strand geometry remain unassigned. This does not assert one symmetric primer species.')
 pool2=g.input(key+'_amp_pool','source SurePlex Amplification Buffer; oligo components and their location unassigned',M.CONDITION)
 enzyme2=(g.materials['sample_amp_enzyme'] if key=='control' else g.input(key+'_amp_enzyme','SurePlexAmplificationEnzyme stock; shared reagent source across parallel reactions',M.CATALYST))
 wga=g.product(key+'_wga','14-cycle '+key+' nativeWGA material',M.MOLECULAR,state=ss(key+' WGA',(role('possible no-cell/background artifact products; no genomicidentity assigned' if key=='control' else 'native genomicDNA-derived WGA population; exact primerjunctions unassigned'),)))
 if key=='sample':wga.state.strand('oppositeWGA',[role('opposite native WGA DNAstrand; terminal andpairing coordinatesunassigned')],label='source double-stranded WGA strandclass')
 g.react(key+'_wga_stage','NativeSurePlex '+key+' amplification14cycles',K.MAGIC,((pre,I.TEMPLATE),(pool2,I.CONDITION),(enzyme2,I.CATALYST)),(wga,),reason='Source supplies14-cycle program and commercialAmplificationBuffer/Enzyme. Native pool/strand architecture and product-class chemistry are unassigned; no inferred isothermal mechanism or symmetricPCR is supplied.')
 return wga

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cell=g.input('cell','singlecell in2.5µlPBS; sourceSurePlexv02',M.BIOLOGICAL,state=ss('cell DNA',(role('cellular genomicDNA'),)),primary=True)
 lysis=g.input('lysis','CellExtractionBuffer andExtractionEnzyme/DilutionBuffer; sourceheatprogram',M.CONDITION)
 dna=g.product('released','released source genomicDNA',M.MOLECULAR,state=ss('released DNA',(role('native extracted genomicDNA'),)))
 g.react('extract','SurePlexcell lysis/extraction',K.EXTRACTION,((cell,I.TEMPLATE),(lysis,I.CONDITION)),(dna,))
 blank=g.input('no_cell','parallel1XPBS no-cell reaction; no samplemixing',M.MOLECULAR,state=ss('parallelno-cell control',(role('PBScontrol; possible backgroundtemplate unassigned'),)))
 wga=commercial_stages(g,'sample',dna)
 artifact=commercial_stages(g,'control',blank)
 g.mark_terminal(Terminal(artifact.id,Disposition.UNRESOLVED,'GuidanceNote allows possible no-cell artifacts to generate unmapped sequencingreads which are computationally filtered. Their exact material/decision to sequence/recovery remains unassigned; this is not physical discard.'))
 amount=g.input('input_amount','5µl diluted0.2ng/µl unpurifiedSurePlex product;1nginput',M.CONDITION)
 chosen=g.product('chosen','1ngWGA aliquot',M.MOLECULAR,state=wga.state)
 rest=g.product('wga_remainder','WGA outside selected1nginput',M.MOLECULAR,state=wga.state)
 g.react('sample_input','Dilute/select1ng SurePlexinput',K.SAMPLING,((wga,I.TEMPLATE),(amount,I.CONDITION)),(chosen,rest),changes_molecule=False)
 g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Material outside selected libraryaliquot; disposition unassigned.'))
 donor=g.input('transposome','VeriSeqATM loaded adapters; nativeoligos unassigned',M.ADAPTER)
 tn5=g.input('tn5','nativeATM transposase;55C5min',M.CATALYST)
 tagged=g.product('tagged','native partially tagged SurePlexDNA',M.MOLECULAR,state=native('tagmented; nativegap andtransferredstrand geometry unassigned'))
 other=g.product('other_tagged','other nativeTn5 endclasses',M.MOLECULAR)
 g.react('tagment','Tagment SurePlexDNA withTD/ATM',K.TAGMENTATION,((chosen,I.TEMPLATE),(donor,I.STRUCTURAL),(tn5,I.CATALYST)),(tagged,other))
 g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Endeligibility is not physical removal orproof of zeroamplification; nativeclass recovery through beads unassigned.'))
 nt=g.input('nt','5µlNT;roomtemperature5min',M.CONDITION)
 neutral=g.product('neutralized','native NT-neutralized taggedDNA',M.MOLECULAR,state=tagged.state)
 g.react('neutralize','Neutralize tagmentation withNT',K.QUENCHING,((tagged,I.TEMPLATE),(nt,I.CONDITION)),(neutral,),changes_molecule=False)
 f=g.input('index5','namedS503/S504 primerfamily; nativewidth/member/strand unassigned',M.OLIGO)
 r=g.input('index7','namedN701–N712 primerfamily; nativewidth/member/strand unassigned',M.OLIGO)
 pol=g.input('npm','NexteraPCRMasterMix',M.CATALYST)
 extended=g.product('extended','native first-extension products before thermalPCR',M.MOLECULAR,state=native('72C3min nativeextension; exact gap/end derivation unassigned'))
 g.react('initial_extension','Initial72C3min extension before denaturation',K.EXTENSION,((neutral,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST)),(extended,))
 amplified=g.product('amplified','sourceindexed double-stranded library population',M.MOLECULAR,state=native('indexedPCR copiedDNA',indexed=True))
 residual=g.product('pcr_residual','co-carriedoriginals and off-targetPCR products',M.MOLECULAR,state=extended.state)
 g.react('index_pcr','95C then12cycles95/55/72C; final72C5min',K.PCR,((extended,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST)),(amplified,residual),geometry=PCRGeometry.unavailable(extended.state,reason='Actual NT-treated and initiallyextended native taggedDNA; complete kitprimer/gap/strand geometry unassigned. No final indexedduplex is used as template.'))
 g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'PCR does not erase originaltemplates or physically remove off-targets; later bead recovery unassigned.'))
 cleaned=cleanup(g,'ampure',amplified,'45µlAMPure per45µlPCR1:1; two80%EtOH washes andRSBelution')
 bound=cleanup(g,'norm_bind',cleaned,'LNA1/LNB1 normalizationbeads;30minbinding; discard unboundlibrary',True)
 washed1=cleanup(g,'norm_wash1',bound,'First45µlLNW1wash; discard washsolution')
 washed2=cleanup(g,'norm_wash2',washed1,'Second45µlLNW1wash; discard washsolution')
 alkali=g.input('naoh','30µl0.1NNaOH;5min normalizationelution',M.CONDITION)
 normalized=g.product('normalized_ss','alkaline-eluted source single-stranded library',M.MOLECULAR,state=native('normalized source library',indexed=True,single=True))
 bead_remaining=g.product('normalization_remaining','normalizationbead-associated/untransferredDNA; exactstrandclass unassigned',M.MOLECULAR)
 g.react('elute_ss','Elute source ssDNA fromnormalizationbeads withNaOH',K.STRAND_SEPARATION,((washed2,I.TEMPLATE),(alkali,I.CONDITION)),(normalized,bead_remaining))
 g.mark_terminal(Terminal(bead_remaining.id,Disposition.UNRESOLVED,'Source transferssupernatant; bead-associated counterpart/remainder molecularidentity unassigned.'))
 stop=g.input('lns1','25µlLNS1 innewplate',M.CONDITION)
 transferred=g.product('lns1_transferred','25µlnormalizationsupernatant transferredintoLNS1',M.MOLECULAR,state=normalized.state)
 untransferred=g.product('untransferred','source alkalineeluate outside25µltransfer',M.MOLECULAR,state=normalized.state)
 g.react('transfer_lns1','Transfer25µlss eluate into25µlLNS1',K.SAMPLING,((normalized,I.TEMPLATE),(stop,I.CONDITION)),(transferred,untransferred),changes_molecule=False)
 g.mark_terminal(Terminal(untransferred.id,Disposition.UNRESOLVED,'Outside source transfer; physicalfate unassigned.'))
 pool=g.input('pool_context','other independentlyprepared normalizedss libraries; equalvolumes/sample indices retained',M.MOLECULAR,state=native('other source-normalized libraries; members unassigned',indexed=True,single=True))
 pooled=g.product('pooled','equal-volume normalizedlibrary pool',M.MOLECULAR,state=transferred.state)
 g.react('pool','Pool equalvolumes of normalizedlibraries',K.POOL,((transferred,I.TEMPLATE),(pool,I.TEMPLATE)),(pooled,),changes_molecule=False)
 ht1=g.input('ht1','HT1 dilution and sourceheatdenaturation beforeMiSeq loading',M.CONDITION)
 library=g.product('library','normalized pooledsingle-stranded library forMiSeq',M.LIBRARY,state=pooled.state)
 g.react('prepare_loading','Dilute pooledlibraries inHT1; sourceheatdenaturation beforeloading',K.THERMAL_TREATMENT,((pooled,I.TEMPLATE),(ht1,I.CONDITION)),(library,),changes_molecule=False)
 g.sequence(library,SequencingPlan(limitation=SOURCE_LIMITATION))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()
REVIEWED=('cell','lysis','extract','no_cell','sample_preamp_pool','sample_preamp_enzyme','sample_preamp_stage','sample_amp_pool','sample_amp_enzyme','sample_wga_stage','control_preamp_pool','control_preamp_stage','control_amp_pool','control_wga_stage','control_wga','input_amount','sample_input','transposome','tn5','tagment','nt','neutralize','index5','index7','npm','initial_extension','index_pcr','ampure_reagent','ampure','norm_bind_reagent','norm_bind','norm_wash1_reagent','norm_wash1','norm_wash2_reagent','norm_wash2','naoh','elute_ss','lns1','transfer_lns1','pool_context','pool','ht1','prepare_loading','library')
