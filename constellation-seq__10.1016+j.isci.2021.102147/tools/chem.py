from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition

def unknown_scene(label='Unassembled whole molecular product'):
    sc=Scene();sc.strand('unassembled',[Segment('unassembled','[unassembled]',placeholder=True)]);return sc

def add(g,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):
    return (g.input(key,label,kind,**kw),role)
def condition(g,key,label):return add(g,key,label,K.CONDITION,R.CONDITION)
def opaque(g,key,label,x,ins=(),kind=S.MAGIC):
    y=g.product(key+'_product',label+'; whole molecular identities/junctions unassembled',K.MOLECULAR)
    kw={'reason':'Source-positive operation; whole kit composition, substrate context or product geometry unassembled, not publication absence.'} if kind==S.MAGIC else {}
    g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(ins),(y,),**kw);return y
def cleanup(g,key,label,x):
    y=g.product(key+'_kept','Source retained preparation; recovery/carryover unassigned',K.MOLECULAR);w=g.product(key+'_removed','Source workup fraction',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),condition(g,key+'_condition',label)),(y,w));g.discard(w,'Physical workup, not perfect molecular selectivity.');return y
def aliquot(g,key,label,x,remainder_label='Untransferred source stock; later fate unassigned'):
    y=g.product(key+'_chosen',label,K.MOLECULAR);r=g.product(key+'_remainder',remainder_label,K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder_label));return y
def qc(g,key,label,x):
    a=g.product(key+'_analytical',label,K.MOLECULAR);k=g.product(key+'_kept','Retained production stock',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(a,k));g.mark_terminal(Terminal(a.id,Disposition.QC,'Analytical aliquot not returned; amount unassigned.'));return k
TITLE='Constellation-Seq: selected published Drop-seq library branch'
AUTHORITY='Published mainPMC7900351 Results91/97/111; transparentMethods mmc1 lines761–813; published TableS2 mmc3 CDS52 whole orderings. Three originals; no predecessor recipe inheritance.'
HANDLE='TATCAACGCAGAGT'
SOURCE_PRIMERS=[('LDHB', 'TATCAACGCAGAGTGGTGTGAATGTGGCAGGTGTT'), ('IL7R', 'TATCAACGCAGAGTTGTCTGAGCCTCTGTCTGTTACTG'), ('S100A9', 'TATCAACGCAGAGTGCAACATAGAGACCATCATCAACA'), ('LYZ', 'TATCAACGCAGAGTGATGTCCGTCAGTATGTTCAAGGT'), ('CD79A', 'TATCAACGCAGAGTGAAACGATGGCAGAACGA'), ('MS4A1', 'TATCAACGCAGAGTTGAACGGTGAGGTGTCAGGTC'), ('CCL5', 'TATCAACGCAGAGTGAGAAGAAATGGGTTCGGGAGTA'), ('GZMK', 'TATCAACGCAGAGTAGTTACTACAACGGCGACCCTTT'), ('FCGR3A', 'TATCAACGCAGAGTAGCCGTAGTGGAATTAACAGGAA'), ('LST1', 'TATCAACGCAGAGTTCCCTGATCCCTGACCTAATGTT'), ('GZMB', 'TATCAACGCAGAGTGAGCCGACCCAGCAGTTTATC'), ('PRF1', 'TATCAACGCAGAGTCCAACGCAAATTCGCAAACTTTCT'), ('FCER1A', 'TATCAACGCAGAGTAGTGAACCTGTGTACCTGGAAGTC'), ('CLEC10A', 'TATCAACGCAGAGTTCACGACCTCCTCGCAAGAC'), ('CD207', 'TATCAACGCAGAGTCATCAGAACCGTGACAGGACAG'), ('CD1A', 'TATCAACGCAGAGTGCCTCCTCGTCACCCTTCTC'), ('HLA-DRB1', 'TATCAACGCAGAGTTGGAATGGAGAGCACGGTCT'), ('CD74', 'TATCAACGCAGAGTGGCTCTTTCGTCACCCTTGGA'), ('CD86', 'TATCAACGCAGAGTGGATGAATGGAAGGAGGCTTAGGA'), ('CD274', 'TATCAACGCAGAGTGTCATCACTACACAGCCCTCCTAA'), ('PSME1', 'TATCAACGCAGAGTAGCAGGAGAAGGAAGACAAGGA'), ('STAT3', 'TATCAACGCAGAGTACAGGTGAGGCAGAACAGC'), ('SQSTM1', 'TATCAACGCAGAGTCCCTGTCCCTCCTAACAAGTGTA'), ('NFKB1', 'TATCAACGCAGAGTGCTCCGAGACAGTGACAGTG'), ('SOCS5', 'TATCAACGCAGAGTGCAAAATGGGGATTGAAGGGACTT'), ('IRF4', 'TATCAACGCAGAGTAGTGAAATGAGCAGCCCTTACAG'), ('SOCS1', 'TATCAACGCAGAGTACGCACTTCCGCACATTCC'), ('SOCS3', 'TATCAACGCAGAGTGGAGGAGACGGGACATCTTTCA'), ('SOCS6', 'TATCAACGCAGAGTTCACACTACCAGCAATAGACAGTT'), ('IDO1', 'TATCAACGCAGAGTGCCAGCAGCCAAAGGAGAA'), ('MYC', 'TATCAACGCAGAGTTCTCAGACGACGGTAGGAATCAG'), ('CCL22', 'TATCAACGCAGAGTAGGAGGCAAAGAGTAGGGTGTAAT'), ('IL10RA', 'TATCAACGCAGAGTTGCCAAGACAGTATCAGACACAGC'), ('CD80', 'TATCAACGCAGAGTAAAGCAGTAGTCAAGGAGGAAAGG'), ('CD70', 'TATCAACGCAGAGTGCAAGTAGACGCCCAGGAG'), ('IRF1', 'TATCAACGCAGAGTATCACGAGGTCAGGAGATCAAGA'), ('CYBB', 'TATCAACGCAGAGTACTTAGGTCCAGCCTGTTCACA'), ('CD207', 'TATCAACGCAGAGTCATCAGAACCGTGACAGGACAG'), ('CD1A', 'TATCAACGCAGAGTGCCTCCTCGTCACCCTTCTC'), ('HLA-DRB1', 'TATCAACGCAGAGTTGGAATGGAGAGCACGGTCT'), ('CD74', 'TATCAACGCAGAGTGGCTCTTTCGTCACCCTTGGA'), ('CD86', 'TATCAACGCAGAGTGGATGAATGGAAGGAGGCTTAGGA'), ('CD274', 'TATCAACGCAGAGTGTCATCACTACACAGCCCTCCTAA'), ('LYN', 'TATCAACGCAGAGTTGATTGCCTGAACACCTGAACATC'), ('MALT1', 'TATCAACGCAGAGTAAGCAGGAAGGCAGAGCATCA'), ('BIRC2', 'TATCAACGCAGAGTCCTTGTGGTCATCTGGTAGTATGC'), ('HLA-DQA1', 'TATCAACGCAGAGTATCGCCATCTACAGGAGCAGAAG'), ('HLA-DPA1', 'TATCAACGCAGAGTATCGCCATCTACAGGAGCAGAAG'), ('NFKBIZ', 'TATCAACGCAGAGTTGATCCTCCTACCTGAGCCTTCT'), ('MLANA', 'TATCAACGCAGAGTCGAACCTTGACCGACATGAACT'), ('TYRP1', 'TATCAACGCAGAGTAGTCTCAGGTTCAAGGCTACAACA'), ('IGFBP7', 'TATCAACGCAGAGTGCTGTGAGGTCATCGGAATCC')]
def primer_inventory():
    sc=Scene()
    for i,(name,seq) in enumerate(SOURCE_PRIMERS):sc.strand('primer_'+str(i),[Segment('handle_'+str(i),HANDLE),Segment('target_'+str(i),seq[len(HANDLE):])])
    return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    stock=g.input('bead_stock','Source extracted/purified Drop-seq beads after no-TSO reverse transcription; prior bead oligo/RNA/RT product geometry unassembled',K.BIOLOGICAL,primary=True)
    x=aliquot(g,'bead_aliquot','Source aliquot of 2000 beads, approximately100 STAMPs',stock)
    pool=add(g,'target_primers','Source52-target TableS2 hybrid-primer pool;10uM pool wording retained without stock/final concentration inference',K.OLIGO,R.STRUCTURAL,state=primer_inventory())
    mix=g.input('linear_ready_mix','Source25ul2x KAPA HiFi HotStart ReadyMix in50ul amplification preparation',K.MOLECULAR)
    enzyme=g.product('retained_kapa','Same-addition KAPA enzyme component',K.CATALYST);nt=g.product('retained_nt','Same-addition ReadyMix nucleotide component; dose unassigned',K.MOLECULAR);buf=g.product('retained_buffer','Same-addition ReadyMix buffer component',K.CONDITION)
    g.react('linear_mix_components','Declared co-present components of one ReadyMix addition, not physical fractionation',S.GENERIC,((mix,R.STRUCTURAL),),(enzyme,nt,buf),changes_molecule=False)
    x=opaque(g,'linear_extension','Source20 rounds of single-primer linear amplification at60C; no SMART reverse-primer PCR at this stage',x,(pool,(enzyme,R.CATALYST),(nt,R.STRUCTURAL),(buf,R.CONDITION),condition(g,'linear_water','Source24.6ulwater/50ul final; primer dosing convention unassigned')),S.EXTENSION)
    x=opaque(g,'standard_pcr','Source then continues standard Drop-seq PCR/library workflow; historical PCR stock/components and primer landings unassembled',x,(add(g,'pcr_preparation','Source PCR preparation with enzyme/nucleotide/primer components unassembled'),))
    x=cleanup(g,'ampure1','Source AMPureXP1:0.6 first purification',x);x=cleanup(g,'ampure2','Source AMPureXP1:0.6 second purification',x);x=qc(g,'bioanalyzer','Source Bioanalyser aliquot before tagmentation',x)
    x=opaque(g,'tagmentation','Source tagmentation/NGS library preparation; kit-specific products unassembled',x,(add(g,'tagment_preparation','Source Nextera preparation and subsequent library amplification stocks; complete composition unassembled'),))
    lib=g.product('library','Source selected Constellation Drop-seq sequencing library',K.LIBRARY);g.react('loading','Source NextSeq loading preparation',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='Selected published Drop-seq branch NextSeq500 paired20x50bp. Whole bead/cell barcode/UMI/cDNA/PCR/tagmentation geometry unassembled; no whole-library or read-landing claim.'))
    g.add_source_item(SourceItem('unassembled','Historical upstream Drop-seq source, full library geometry and production PCR orderings ordinary unreviewed. Separate10X branch not imported.','Capped positive preparation scope; not absence or full chemistry audit.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('Targeted linear extension precedes opaque PCR',list(g.steps).index('linear_extension')<list(g.steps).index('standard_pcr'),True)
    c('Linear amplification is not declared PCR',g.steps['linear_extension'].kind,S.EXTENSION)
    c('No TSO addition borrowed',any('tso' in x.id.lower()for x in g.materials.values()),False)
    c('Actual bead stock remainder retained','bead_aliquot_remainder'in g.terminals,True)
    c('Same ReadyMix supplies structural nucleotides',any(p.material=='retained_nt'and p.role==R.STRUCTURAL for p in g.steps['linear_extension'].inputs),True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.isci.2021.102147_supp_mmc3.xlsx.txt'
    if have(c,Source(p,'https://doi.org/10.1016/j.isci.2021.102147')):
        source=p.read_text()
        for i,(name,seq)in enumerate(SOURCE_PRIMERS):c('Whole original TableS2 row '+str(i+4)+' '+name,seq in source,True)
    c('Private partial graph remains incomplete',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
