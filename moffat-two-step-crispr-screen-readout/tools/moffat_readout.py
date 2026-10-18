from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Scene,Segment,feature,revcomp,complement_segments
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry,PCRPrimer
import illumina as il
import seqprimers as sp
TITLE='Moffat two-step CRISPR screen readout'
NOTES='01_moffat-two-step.html'
SOURCE='MoffatLab <a href="https://media.addgene.org/cms/filer_public/44/9f/449f22c6-f818-4af6-8af4-2bab1d20f00a/moffat_lab_crispr_screen_protocol_v2.pdf">PooledCRISPRScreenProtocolv2, printedREV.20190723</a>, §3.5.3/Table1/§3.5.4, LCV2::TKOv3 branch.'
SUMMARY='Parallel genomic-locus enrichment PCRs are pooled, measured on an aliquot and sampled into a nested dual-index PCR. The source200bp gel fraction supplies the final library.'
CAVEAT='Source LCV2 branch only; the pLCKO2 Table2 primer family and optional one-PCR/stagger variant are separate. D501/D701 is a source representative combination. PCR1 native flank lengths remain unassigned. The PCR2 geometry receives the actual enriched precursor without flowcell arms; only its source-known nested sites are placed.'
SEQUENCING_ENDING='Source run:21 dark cycles then26 acquired Read1 cycles, Index1=8 andIndex2=8. The compatible primer-binding reference precedes the initial21 dark additions; physical identifiers and acquired-cycle assignments remain distinct. No Read2.'
READ_LENGTHS={'Index 1 (i7)':8,'Index 2 (i5)':8}
READ_SPECS_NOTE='Read1 acquisition is26 after21dark additions; exact acquiredstart is not represented by the ordinary primer-site geometry.'
SEQ_PRIMERS=(sp.TRUSEQ['I1'],sp.TRUSEQ['I2'])
GUIDE=feature('guide','guide_barcode','unknown')
I5=feature('sample_i5','sample_index','whitelist',whitelist='sourceTable1 D50x')
I7=feature('sample_i7','sample_index','whitelist',whitelist='sourceTable1 D70x')
LCV2_PCR1_F = "GAGGGCCTATTTCCCATGATTC"
LCV2_PCR1_R = "GTTGCGAAAAAGAACGTTCACGG"
D501_F = ("AATGATACGGCGACCACCGAGATCTACACTATAGCCTACACTCTTT"
          "CCCTACACGACGCTCTTCCGATCTTTGTGGAAAGGACGAAACACCG")
D701_R = ("CAAGCAGAAGACGGCATACGAGATCGAGTAATGTGACTGGAGTTCA"
          "GACGTGTGCTCTTCCGATCTACTTGCTATTTCTAGCTCTAAAAC")
PCR2_F_ANNEAL = D501_F[-22:]   # blue annealing part of every D50x-F primer
PCR2_R_ANNEAL = D701_R[-24:]   # blue annealing part of every D70x-R primer


def pcr1_product():
 return Construct([Segment('PCR1 F site',LCV2_PCR1_F),role('variable-length native leftLCV2 flank'),Segment('nested hU6 site',PCR2_F_ANNEAL),Segment('guide','N'*20,placeholder=True,feature=GUIDE),Segment('nested scaffold site rc',revcomp(PCR2_R_ANNEAL)),role('variable-length native rightLCV2 flank'),Segment('PCR1 R site rc',revcomp(LCV2_PCR1_R))],name='native enrichedPCR1 precursor; sourceapprox600bp')

def pcr1_mixture():
 sc=Scene.duplex(list(pcr1_product()),label='native PCR1 daughter component; variable flanks',unpaired=tuple(s.name+"'" for s in pcr1_product() if s.is_role_token()))
 sc.strand('possible_genomes',[role('possible original genomic DNA; guide identity retained, recovery unassigned',molecular_feature=GUIDE)],label='independent possible genomic originals')
 sc.strand('possible_other',[role('possible co-carried RNA, non-target DNA and primers; no physical removal before PCR2')])
 return sc

def final_library():
 lib=Construct([Segment('P5',il.P5),Segment('i5',D501_F[len(il.P5):len(il.P5)+8],feature=I5),Segment('Read1',il.TRUSEQ_READ1),Segment('nestedF',PCR2_F_ANNEAL),Segment('guide','N'*20,placeholder=True,feature=GUIDE),Segment('nestedR rc',revcomp(PCR2_R_ANNEAL)),Segment('Read2 rc',revcomp(il.TRUSEQ_READ2)),Segment('i7 rc',revcomp(D701_R[len(il.P7):len(il.P7)+8]),feature=I7),Segment('P7 rc',il.P7_RC)],name='sourceD501/D701 representativePCR2 daughter')
 if not lib.top().startswith(D501_F) or not lib.top().endswith(revcomp(D701_R)):raise ValueError('published PCR2 primer decomposition differs from canonical arms')
 return lib

def pcr2_mixture():
 sc=Scene.duplex(list(final_library()),label='representative PCR2 daughter component')
 sc.strand('possible_precursors',[role('possible original PCR1/genomic DNA carrying guide identity; recovery unassigned',molecular_feature=GUIDE)])
 sc.strand('possible_other',[role('possible RNA, other DNA and primers; source gel fraction does not establish complete purity')])
 return sc

FINAL_LIBRARY=final_library()

def oligos():
 for name in ('LCV2_PCR1_F','LCV2_PCR1_R','D501_F','D701_R'):yield name+':5′-'+globals()[name]+'-3′'

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 state=Scene();state.strand('genome',[role('native genomicDNA containing integratedLCV2::TKOv3 cassette; guide coordinates and native extent unassigned',molecular_feature=GUIDE)]);state.strand('other',[role('possible co-extracted RNA/otherDNA; purityunassigned')])
 gdna=g.input('genomic_dna','screen gDNA;source100µg/sample',M.BIOLOGICAL,state=state,primary=True)
 allocation=g.input('allocation','3.5µg/50µl reaction; parallel reactions tosourcecoverage',M.CONDITION)
 a=g.product('reaction_class','one representativePCR1 reaction aliquot',M.MOLECULAR,state=state);peers=g.product('parallel_reactions','remaining independently amplifiedPCR1 reaction classes',M.MOLECULAR,state=state)
 g.react('allocate','Allocate parallelPCR1 reaction classes',K.SPLIT,((gdna,I.TEMPLATE),(allocation,I.CONDITION)),(a,peers),changes_molecule=False)
 f=g.input('f1','sourceLCV2 forward',M.OLIGO,state=Construct([Segment('F1',LCV2_PCR1_F)]));r=g.input('r1','sourceLCV2 reverse',M.OLIGO,state=Construct([Segment('R1',LCV2_PCR1_R)]));pol=g.input('polymerase','2×NEBNextUltraIIQ5MasterMix',M.CATALYST)
 amps=[]
 for key,x in [('one',a),('peers',peers)]:
  out=g.product('pcr1_'+key+'_out','native sourceapprox600bpPCR1 daughters',M.MOLECULAR,state=pcr1_mixture())
  g.react('pcr1_'+key,'SourcePCR1:98°C30s;25cycles98°C10s/66°C30s/72°C15s;72°C2min',K.PCR,((x,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST)),(out,),geometry=PCRGeometry.unavailable(x.state,reason='Native genomicDNA is theactual precursor; no finishedPCR1/2 product is substituted for it. Actual genomicflanks remain unassigned.'))
  amps.append(out)
 qc_rule=g.input('qc_rule','source2µlPCR1sample to1%analyticalgel',M.CONDITION)
 qc=g.product('qc_aliquot','2µlPCR1QC aliquot',M.MOLECULAR,state=amps[0].state);rem=g.product('qc_remainder','remainingPCR1reaction',M.MOLECULAR,state=amps[0].state)
 g.react('qc_sample','Take sourcePCR1QC aliquot',K.SAMPLING,((amps[0],I.TEMPLATE),(qc_rule,I.CONDITION)),(qc,rem),changes_molecule=False)
 g.mark_terminal(Terminal(qc.id,Disposition.UNRESOLVED,'Analyticalgel product is not recombined; later physicalfateunassigned.'))
 pool=g.product('pcr1_pool','pooled nativePCR1 daughters andpossibleco-carriedoriginals',M.MOLECULAR,state=pcr1_mixture())
 g.react('pool_pcr1','Pool parallelPCR1 reactions forone sample',K.POOL,((rem,I.TEMPLATE),(amps[1],I.TEMPLATE)),(pool,),changes_molecule=False)
 sample_rule=g.input('pcr2_allocation','source5µlpooledPCR1into50µlPCR2',M.CONDITION)
 sample=g.product('pcr2_input','actual5µlPCR1precursor aliquot',M.MOLECULAR,state=pool.state);rest=g.product('pcr1_rest','remainingPCR1pool',M.MOLECULAR,state=pool.state)
 g.react('sample_pcr2','Allocate5µlPCR1pool toPCR2',K.SAMPLING,((pool,I.TEMPLATE),(sample_rule,I.CONDITION)),(sample,rest),changes_molecule=False);g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'UnsampledPCR1remainder; notsource-prescribeddiscard.'))
 f2=g.input('f2','sourceD501-F',M.OLIGO,state=Construct([Segment('P5',D501_F[:len(il.P5)]),Segment('i5',D501_F[len(il.P5):len(il.P5)+8],feature=I5),Segment('forward arm and site',D501_F[len(il.P5)+8:])]));r2=g.input('r2','sourceD701-R',M.OLIGO,state=Construct([Segment('P7',D701_R[:len(il.P7)]),Segment('i7',D701_R[len(il.P7):len(il.P7)+8],feature=I7),Segment('reverse arm and site',D701_R[len(il.P7)+8:])]))
 amp=g.product('pcr2_daughters','dual-indexPCR2daughters',M.MOLECULAR,state=pcr2_mixture())
 geom=PCRGeometry.unavailable(sample.state,reason='Actual pooled PCR1 mixture includes enriched daughters and possible co-carried genomes, RNA and primers; no physical purification precedes PCR2. The separate target-component panel locates only nested sites on the represented PCR1 daughter, never on a completed PCR2 product.')
 g.react('pcr2','SourcePCR2:98°C30s;10cycles98°C10s/55°C30s/65°C15s;65°C5min',K.PCR,((sample,I.TEMPLATE),(f2,I.STRUCTURAL),(r2,I.STRUCTURAL),(pol,I.CATALYST)),(amp,),geometry=geom)
 gel=g.input('gel','source2%agarose;DarkReader;gelextractionkit',M.SELECTION)
 library=g.product('library','source200bp excisedPCR2library fraction',M.LIBRARY,state=pcr2_mixture());w=g.product('gel_waste','unexcisedgel fractions',M.WASTE)
 g.react('gel_select','Run source50µlPCR2 andexcise200bp band',K.CLEANUP,((amp,I.TEMPLATE),(gel,I.SELECTION)),(library,w),changes_molecule=False);g.discard(w,'Physicalgel outsideexcisedband; doesnotprove speciespuritywithinband.')
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,acquired_roles=('Index 1 (i7)','Index 2 (i5)'),read_lengths=READ_LENGTHS,unplaced_reads=(sp.read_spec('Read 1',26,'Source21dark additions precede26 acquiredcycles; the model doesnotconvertthat offset into an ordinary primer-start cyclemap.','include'),),limitation='SourceRead1 acquisition26 follows21 initialdark additions; the acquiredstart is leftunplaced bythis primer-start geometry. D501/D701 target component of the source gel fraction; possible co-migrating classes remain unresolved. Source index reads8/8; no Read2.'))
 g.declare_reviewed_scope(('genomic_dna','allocation','allocate','f1','r1','polymerase','pcr1_one','pcr1_peers','qc_rule','qc_sample','pool_pcr1','pcr2_allocation','sample_pcr2','f2','r2','pcr2','gel','gel_select','library'),SOURCE)
 return g.require_publishable()

FINAL_CAPTION='Source Table1 D501/D701 representative nested-PCR daughter; prior variable-lengthPCR1 flanks are outside the selected molecule.'

SEQUENCING_INTRO='Index sites locate on the final molecule. Read1 acquisition after21darkcycles is declared separately from ordinary primer-start geometry.'
TOPOLOGY_PANELS=(('Compatible Read1 binding before the source dark additions',sp.section(FINAL_LIBRARY,(sp.TRUSEQ['R1'],),required_roles=(),show_read_content=False,intro='This binding reference starts before21dark additions; it doesnotassign acquiredRead1featurecycles.')), )

from chemdraw import panel
TOPOLOGY_PANELS += (('Nested primer sites on the PCR1 target component',panel(PCRGeometry.exact(pcr1_product(),PCRPrimer('D501-F',D501_F),PCRPrimer('D701-R',D701_R),statement='Target daughter component only; variable native flanks and co-carried molecular classes are not assigned exact PCR2 mixture geometry.').rows(),caption='Computed sites on the actual enriched PCR1 component; its pool still contains possible original genomes, RNA and primers.')), )
