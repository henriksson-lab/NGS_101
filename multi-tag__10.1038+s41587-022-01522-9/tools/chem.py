"""Meers2022 single-cell MulTI-Tag: local source-owned ordered graph."""
from chemdraw import Construct,Scene,Segment,feature,revcomp
from protocol_model import role
from protocol_graph import ProtocolGraph,InputRole as I,MaterialKind as M,StepKind as K,Terminal,Disposition,SequencingPlan
from pcr import PCRGeometry
import seqprimers as sp
TITLE='MulTI-Tag — source single-cell sequential target barcoding'
NOTES='01_multi_tag.html'
SOURCE='<a href="https://doi.org/10.1038/s41587-022-01522-9">Meers et al., Nature Biotechnology2022</a>, published Methods107–121 and published SupplementaryTable1 (manifest-proven byte-identical cached workbook).'
SUMMARY='Sequential H3K4me2/H3K36me3/H3K27me3 conjugate transposomes mark target identity. Secondary-antibody i7 transposomes mark96wells, followed by pooled nuclei processing and ICELL8 indexed PCR.'
CAVEAT='Source examples are selected family members, not universal sequences. The selected single-cell target order keeps H3K27me3 last. Exact native insertion gaps/strands, ICELL8 first-cycle extension and PCR programme remain unassigned: Methods cites predecessor62, which was not inspected in this capped pass. Bulk58°C/72°C gap-fill conditions are not assigned to microwells. Chemical DNase cleavage and nonamplification are not physical disposal. Original target lineages remain independent and original/native co-carriage is unresolved.'
SEQ_PRIMERS=()
TARGET=feature('target_barcode','feature_barcode','whitelist',whitelist='published target-associated P5_i5 family')
WELL=feature('i7_transposome_well','cell_barcode','whitelist',whitelist='published96 P7_i7 wells')
I5=feature('icell8_i5_primer','cell_barcode','whitelist',whitelist='published72 i5 primers')
I7=feature('icell8_i7_primer','cell_barcode','whitelist',whitelist='published72 i7 primers')
def identifier(name,f):return Segment(name,'N'*8,'cbc',placeholder=True,feature=f)
TARGET_BC=identifier('target barcode',TARGET)
TN5_WELL=identifier('i7 well barcode',WELL)
I5_PCR=identifier('microwell i5 barcode',I5)
I7_PCR=identifier('microwell i7 barcode',I7)
INSERT=role('target-proximal genomic DNA; native extent unassigned')
FINAL_LIBRARY=Construct([role('native primer-derived P5 end'),I5_PCR,TARGET_BC,INSERT,TN5_WELL,I7_PCR,role('native primer-derived P7 end')],name='MulTI-Tag source-bounded daughter library')
FINAL_CAPTION='Source identifiers8nt each; complete native primer-to-insert geometry remains bounded.'
SEQUENCING_INTRO='Methods: HiSeq PE25×25; native read/index primer landing and acquired barcode/index cycles remain unassigned.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec('Read 1',25,'genomic fragment; barcode entry unassigned','include'),sp.read_spec('Read 2',25,'opposite genomic fragment; barcode entry unassigned','include'),sp.read_spec('Index 1 (i7)','unassigned','Known8nt identifiers do not establish acquisition allocation','undetermined'),sp.read_spec('Index 2 (i5)','unassigned','Native acquisition allocation unassigned','undetermined'))
# Supplementary Table 1 example rows. Each is one member of a printed family; the
# experiment chooses the row (target, well or nanowell), so none is universal.
TN5MEREV = "CTGTCTCTTATACACATCT"            # README: [phos]CTGTCTCTTATACACATCT
P5_I5_1 = "TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG"  # H3K27me3, TATAGCCT
P7_I7_13 = "GTCTCGTGGGCTCGGCTGTCCCTGTCCGTCTACATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG"  # first single-cell row, GTCTACAT
SCI_AD1_1 = "AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTCTCCACGC"  # i5 primer, TAGATCGC
SCI_AD2_1 = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGCTGTCCCTGTCC"   # i7 primer, TCGCCTTA

P5_I5_3='TCGTCGGCAGCGTCTCCACGCCCTATCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG'
P5_I5_6='TCGTCGGCAGCGTCTCCACGCTAATCTTAGCGATCGAGGACGGCAGATGTGTATAAGAGACAG'
TARGETS=(('H3K4me2',P5_I5_3),('H3K36me3',P5_I5_6),('H3K27me3',P5_I5_1))
def adapter_scene(seq,*,target=False,conjugated=False):
 prefix=len("TCGTCGGCAGCGTCTCCACGC") if target else 24
 # Source prefixes precede the barcode; retain every source base, including fixed bridge.
 barcode=seq[prefix:prefix+8]
 if target and barcode not in ('CCTATCCT','TAATCTTA','TATAGCCT'):raise ValueError('source target barcode boundary drift')
 if not target:
  prefix=seq.index('GTCTACAT');barcode=seq[prefix:prefix+8]
 parts=[Segment('native5′adapter prefix',seq[:prefix]),Segment('source barcode',barcode,feature=TARGET if target else WELL),Segment('fixed adapter bridge',seq[prefix+8:-19]),Segment('mosaic end',seq[-19:])]
 sc=Scene();sc.strand('transfer',parts,label='source target/well transfer strand',mod5='antibody–amino conjugate' if conjugated else 'amino' if target else '')
 sc.anneal('nontransfer',[Segment('Tn5MErev',TN5MEREV)],to='transfer',pair=('Tn5MErev','mosaic end'),mod5='p',label='independent source nontransfer strand')
 return sc

def population(stage,targets=0,well=False,indexed=False):
 sc=Scene()
 for k,(label,seq) in enumerate(TARGETS[:targets]):
  parts=[role('native target-adapter boundary; '+stage),Segment(label+' target barcode',seq[len("TCGTCGGCAGCGTCTCCACGC"):len("TCGTCGGCAGCGTCTCCACGC")+8],feature=TARGET),role(label+' genomic fragment; native extent/gap/strand geometry unassigned')]
  if well:parts.append(identifier(label+' i7 well barcode',WELL))
  if indexed:parts=[identifier(label+' installed microwell i5',I5),*parts,identifier(label+' installed microwell i7',I7)]
  sc.strand('target'+str(k),parts,label='independent '+label+' tagged/copy fragment population',mod5='primer-derived daughter end' if indexed else 'possible antibody–amine adapter conjugate; free-adapter contribution unassigned')
 if not targets:sc.strand('genome',[role('nuclear genomic DNA; independent target loci not covalently joined')],label='nuclear chromatin population')
 residual='possible original antibody–adapter-covalent DNA, free-adapter DNA, nontransfer strands and other nonproductive material; physical recovery unassigned' if targets else 'native nuclear protein/other material; physical recovery unassigned; no introduced adapter lineage yet'
 sc.strand('residual',[role(residual)],label='independent possible residual lineages')
 return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 def inp(i,l,k=M.CONDITION,state=None):return g.input(i,l,k,state=state)
 def op(i,l,k,x,extra,state):
  y=g.product(i+'_out',l,M.MOLECULAR,state=state)
  products=[y]
  if k is K.PARTITION:
   other=g.product(i+'_other','other physical partitions outside represented lineage; disposition unassigned',M.MOLECULAR,state=x.state);products.append(other)
  g.react(i,l,k,((x,I.TEMPLATE),*extra),products,changes_molecule=k is not K.GENERIC)
  if k is K.PARTITION:g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Other physical wells outside this example lineage; no disposal inferred.'))
  return y
 def wash(i,l,x):
  sel=inp(i+'_selection',l,M.SELECTION);y=g.product(i+'_out',l+' retained fraction',M.MOLECULAR,state=x.state);w=g.product(i+'_waste','physical removed fraction',M.WASTE)
  g.react(i,l,K.CLEANUP,((x,I.TEMPLATE),(sel,I.SELECTION)),(y,w));g.discard(w,'Source physical wash; class partition/purity unassigned.');return y
 x=g.input('nuclei','1.8million nuclei for source single-cell route',M.BIOLOGICAL,state=population('unprocessed'),primary=True)
 beads=inp('wga_beads','biotinWGA on streptavidinT1 beads;100µl prebound to nuclei',M.SELECTION);x=op('bind_wga','Bind nuclei to WGA beads',K.BINDING,x,((beads,I.SELECTION),),x.state)
 for n,(target,seq) in enumerate(TARGETS,1):
  prefix='target'+str(n)
  oligo=inp(prefix+'_oligo','source5′aminated '+target+' P5_i5 transfer oligo',M.OLIGO,seq)
  me=inp(prefix+'_me','source5′phosphorylated Tn5MErev',M.OLIGO,Scene())
  me.state.strand('MErev',[Segment('MErev',TN5MEREV)],mod5='p')
  annealed=g.product(prefix+'_annealed','annealed '+target+' adapter',M.MOLECULAR,state=adapter_scene(seq,target=True))
  g.react(prefix+'_anneal','Anneal equimolar transfer/MErev strands',K.HYBRIDIZATION,((oligo,I.TEMPLATE),(me,I.STRUCTURAL)),(annealed,))
  antibody=inp(prefix+'_antibody',target+' antibody;Abcam amino-conjugation kit',M.CONDITION)
  conj=op(prefix+'_conjugate','Conjugate5′amino adapter to antibody; source~10:1excess',K.CHEMICAL_CONVERSION,annealed,((antibody,I.CONDITION),),adapter_scene(seq,target=True,conjugated=True))
  free=inp(prefix+'_free','independent16pmol unconjugated annealed same adapter',M.ADAPTER,adapter_scene(seq,target=True));tn=inp(prefix+'_tn5','pA–Tn5~5µM;load30–60minRT',M.CATALYST)
  bound=op(prefix+'_bind','Assemble '+target+' conjugate/free-adapter transposome and bind nuclei',K.BINDING,x,((conj,I.STRUCTURAL),(free,I.STRUCTURAL),(tn,I.CATALYST)),x.state)
  bound=wash(prefix+'_bindwash','Clear supernatant after final conjugate incubation' if n==3 else 'Three Wash300 washes before tagmentation',bound)
  if n==3:
   final_conj,final_free,final_tn=conj,free,tn;x=bound;continue
  mg=inp(prefix+'_mg','Wash30010mMMgCl2;37°C1h')
  x=op(prefix+'_tagment','Tagment '+target+' targeted chromatin',K.TAGMENTATION,bound,((conj,I.STRUCTURAL),(free,I.STRUCTURAL),(tn,I.CATALYST),(mg,I.CONDITION)),population('tagmented',n))
  x=wash(prefix+'_edta','Source5mMEDTA post-tagment wash',x)
 sec=inp('secondary','species-matched secondary antibodies;1hRT',M.SELECTION);x=op('secondary_bind','Bind secondary antibodies',K.BINDING,x,((sec,I.SELECTION),),x.state);x=wash('secondary_wash','Two WashBuffer washes',x)
 well=inp('well_adapter','source P7_i7_13 example among96wells',M.ADAPTER,adapter_scene(P7_I7_13));tn=inp('well_tn5','pA–Tn5 loaded only i7adapter;32pmol/µl enzyme',M.CATALYST)
 wells=inp('well_partition','distribute equally across96 i7transposome wells');x=op('distribute','Distribute nuclei across i7barcode wells',K.GENERIC,x,((wells,I.CONDITION),),x.state)
 x=op('well_bind','Bind i7transposomes1hRT',K.BINDING,x,((well,I.STRUCTURAL),(tn,I.CATALYST)),x.state);x=wash('well_wash','Three Wash300 washes',x)
 mg=inp('well_mg','10mMMgCl2;37°C1h');x=op('well_tag','Final i7 CUT&Tag tagmentation',K.TAGMENTATION,x,((well,I.STRUCTURAL),(final_conj,I.STRUCTURAL),(final_free,I.STRUCTURAL),(final_tn,I.CATALYST),(tn,I.CATALYST),(mg,I.CONDITION)),population('i5/i7 tagmented',3,True))
 pooling=inp('pooling','reaggregate96well nuclei');x=op('pool','Pool nuclei after finaltag',K.GENERIC,x,((pooling,I.CONDITION),),x.state);x=wash('pre_dnase_wash','Two100µlTAPS washes',x)
 dn=inp('dnase_enzyme','RQ1DNase0°C10min thenStopBuffer',M.CATALYST);x=op('dnase','Digest accessible/extranuclear DNA; physical class recovery unassigned',K.DIGESTION,x,((dn,I.CATALYST),),x.state)
 x=wash('dnase_wash','SourceTAPS wash afterDNaseStop',x)
 strainer=inp('strainers','two20µm strainers300g10min',M.SELECTION);passed=g.product('passed','strainer flowthrough nuclei',M.MOLECULAR,state=x.state);ret=g.product('strainer_retained','physically retained clumps/material',M.WASTE)
 g.react('strain','Pass nuclei through strainers',K.SELECTION,((x,I.TEMPLATE),(strainer,I.SELECTION)),(passed,ret));g.discard(ret,'Physical retained fraction; not computational cellQC.')
 dispenser=inp('icell8','DAPI/ICELL8SecondDiluent;857nuclei/µl,640µl dispense');x=op('dispense','Dispense nuclei intoICELL8 nanowells',K.GENERIC,passed,((dispenser,I.CONDITION),),passed.state)
 sds=inp('sds','source microwellSDSdenaturation; conditions delegatedref62');x=op('release','SDS release native taggedDNA',K.DEPROTEINIZATION,x,((sds,I.CONDITION),),population('native released actualPCRprecursor',3,True))
 tr=inp('triton','source microwellTritonX100 neutralization');x=op('neutralize','Neutralize SDS',K.CHEMICAL_CONVERSION,x,((tr,I.CONDITION),),x.state)
 f=inp('i5_primer','source sci_Ad1.1_TAGATCGC example;72family',M.OLIGO,SCI_AD1_1);r=inp('i7_primer','source sci_Ad2.1_TAAGGCGA example;72family',M.OLIGO,SCI_AD2_1);pol=inp('polymerase','native ICELL8 PCR setup/cycling delegatedref62; unassigned',M.CATALYST)
 lib=g.product('amplified','native microwell-indexed daughter/original mixture',M.MOLECULAR,state=population('amplified daughters/native originals unresolved',3,True,True))
 g.react('pcr','Amplify actual released nativeTaggedDNA; first-cycle repair bounded',K.PCR,((x,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST)),(lib,),geometry=PCRGeometry.unavailable(x.state,reason='Actual native tagged strand/gap substrate; published microwell PCR delegatedref62, no bulk programme or completed endpoint substituted.'))
 collect=inp('collect_condition','reaggregate microwell contents');x=op('collect','Pool microwell products',K.GENERIC,lib,((collect,I.CONDITION),),lib.state)
 x=wash('cleanup1','First1.3×HighPrep cleanup',x);x=wash('cleanup2','Second1.3×HighPrep cleanup;20µlTris/EDTA',x)
 final=g.product('library',FINAL_LIBRARY.name,M.LIBRARY,state=x.state);g.react('retain','Retain recovered libraries',K.GENERIC,((x,I.TEMPLATE),),(final,),changes_molecule=False)
 g.sequence(final,SequencingPlan(library=FINAL_LIBRARY,acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':25,'Read 2':25},limitation=SEQUENCING_ENDING))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()

# Explicit curator-written source targets, including all3 selected target rounds.
REVIEWED=('nuclei','wga_beads','bind_wga','target1_oligo','target1_me','target1_anneal','target1_antibody','target1_conjugate','target1_free','target1_tn5','target1_bind','target1_bindwash','target1_bindwash_selection','target1_mg','target1_tagment','target1_edta','target1_edta_selection','target2_oligo','target2_me','target2_anneal','target2_antibody','target2_conjugate','target2_free','target2_tn5','target2_bind','target2_bindwash','target2_bindwash_selection','target2_mg','target2_tagment','target2_edta','target2_edta_selection','target3_oligo','target3_me','target3_anneal','target3_antibody','target3_conjugate','target3_free','target3_tn5','target3_bind','target3_bindwash','target3_bindwash_selection','secondary','secondary_bind','secondary_wash','secondary_wash_selection','well_adapter','well_tn5','well_partition','distribute','well_bind','well_wash','well_wash_selection','well_mg','well_tag','pooling','pool','pre_dnase_wash','pre_dnase_wash_selection','dnase_enzyme','dnase','dnase_wash','dnase_wash_selection','strainers','strain','icell8','dispense','sds','release','triton','neutralize','i5_primer','i7_primer','polymerase','pcr','collect_condition','collect','cleanup1','cleanup1_selection','cleanup2','cleanup2_selection','retain','library')
def oligos():
 return ('Tn5MErev —5′-/5Phos/'+TN5MEREV+'-3′',)+tuple(name+' —5′-'+seq+'-3′' for name,seq in (('P5_i5_1 amino H3K27me3',P5_I5_1),('P5_i5_3 amino H3K4me2',P5_I5_3),('P5_i5_6 amino H3K36me3',P5_I5_6),('P7_i7_13',P7_I7_13),('sciAd1.1',SCI_AD1_1),('sciAd2.1',SCI_AD2_1)))
