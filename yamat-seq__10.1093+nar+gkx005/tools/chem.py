"""Published YAMAT native chemistry; unknown kit members remain local roles."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import Disposition, InputRole as R, MaterialKind as M, ProtocolGraph, SequencingPlan, StepKind as K, Terminal
from protocol_model import molecule, role
import seqprimers as sp

TITLE='YAMAT-seq — Y-shaped adapter mature-tRNA sequencing'
NOTES='01_yamat_seq.html'
SOURCE='Shigematsu et al., Nucleic Acids Research (2017); <a href="https://doi.org/10.1093/nar/gkx005">doi:10.1093/nar/gkx005</a>, Methods and Results Y-adapter description.'
SUMMARY='Deacylate mature tRNAs, anneal discriminator-matched Y adapters, seal both nicks with Rnl2, reverse-transcribe, amplify and select the mature-tRNA library by native PAGE.'
CAVEAT=('Published discrepancy: Methods, “Y-shaped adapter design”, specifies Y-3′-AD '
        '5′-/5Phos/GTATCCAGTTGGAATTCTCGGGTGCCAAGG/3ddC/-3′; Results, Y-AD description, '
        'instead says its 3′ end is OH. The model uses the explicit Methods ordering specification; '
        'the physical reagent remains uncertain. The annealing view checks the source adapter stem '
        'and a representative A-discriminator terminal pair, not a reconstructed complete tRNA fold. '
        'INFERRED — kit primer boundaries and full cDNA extent remain unassigned.')
FINAL_LIBRARY=None
SEQ_PRIMERS=()
SEQUENCING_ENDING=('The defining experiment used 100-nt single-read HiSeq2500 sequencing after PhiX clustering. '
                    'The selected kit members, native index width and library sequencing-primer landing remain unassigned in this source pass.')
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec('Read 1',100,'Source 100-nt single-read acquisition','include'),
 sp.read_spec('Read 2','not acquired','Source single-read sequencing','not_applicable'),
 sp.read_spec('Index 1 (i7)','unassigned','Native index acquisition allocation remains unassigned','undetermined'),
 sp.read_spec('Index 2 (i5)','unassigned','Native index acquisition allocation remains unassigned','undetermined'))
ACQUIRED_ROLES=('Read 1',)
READ_LENGTHS={'Read 1':100}
Y3_SEQUENCE='GTATCCAGTTGGAATTCTCGGGTGCCAAGG'
Y5_DNA='GTTCAGAGTTCTACAGTCCGACGATCACTGGATACTG'
Y3=Segment('Y-3-prime-AD',Y3_SEQUENCE,note='DNA; source ordered 5P')
DDC=Segment('terminal ddC','C',note='physical additional 3′ dideoxycytidine from /3ddC/')
# Distinguish physical RNA U from DNA T; the last two residues are RNA in all four oligos.
Y5_A=Segment('Y-5-prime-AD-A',Y5_DNA+'GA',note='last two residues RNA')
Y5_G=Segment('Y-5-prime-AD-G',Y5_DNA+'GG',note='last two residues RNA')
Y5_C=Segment('Y-5-prime-AD-C',Y5_DNA+'GC',note='last two residues RNA')
Y5_U=Segment('Y-5-prime-AD-U',Y5_DNA+'GU',note='last two residues RNA; physical terminal U')

def oligos():
    return ('Y-3′-AD — 5′-/5Phos/'+Y3_SEQUENCE+'/3ddC/-3′',)+tuple(
        'Y-5′-AD-'+base+' — 5′-'+Y5_DNA+'g'+base.lower()+'-3′ (last two residues RNA)'
        for base in 'AGCU')


def adapter_set_scene():
    sc=Scene()
    for base in 'AGCU':
        sc.strand('y5_'+base,[Segment('Y5 '+base+' DNA',Y5_DNA),
                                Segment('Y5 '+base+' RNA','G'+base,note='RNA')],label='Y-5′-AD-'+base)
    sc.strand('y3',[Y3,DDC],label='Y-3′-AD',mod5='/5Phos/',mod3='/3ddC/')
    return sc


def trna_scene(acylated=False):
    sc=Scene();sc.strand('trna',[role('mature tRNA body, RNA'),role('discriminator N'),
                                Segment('CCA','CCA',note='RNA')],label='mature tRNA',mod5='p',
                         mod3='aminoacyl ester (population includes uncharged RNA)' if acylated else 'OH')
    return sc


def annealed_scene():
    # An explicit representative for the source's A-discriminator/Y5-U class.
    full=Y5_U.top
    candidates=[n for n in range(1,len(Y3_SEQUENCE)+1) if full[:-4].endswith(revcomp(Y3_SEQUENCE[:n]))]
    n=max(candidates)
    sc=Scene()
    sc.strand('y5',[Segment('Y5 DNA arm',full[:-(n+4)]),Segment('Y5 stem',full[-(n+4):-4]),
                     Segment('Y5 overhang DNA','TG'),Segment('Y5 overhang RNA','GU',note='RNA')],label='Y5-U representative')
    sc.anneal('y3',[Segment('Y3 stem',Y3_SEQUENCE[:n]),Segment('Y3 tail',Y3_SEQUENCE[n:]),DDC],to='y5',
               pair=('Y3 stem','Y5 stem'),unpaired=('Y3 tail',DDC.name),mod5='p',mod3='/3ddC/',label='Y3')
    sc.anneal('trna',[role('representative tRNA RNA body; full acceptor fold not laid out'),
                       Segment('A-discriminator CCA','ACCA',note='RNA')],to='y5',
               pair=('A-discriminator CCA','Y5 overhang DNA'),
               unpaired=('representative tRNA RNA body; full acceptor fold not laid out',),mod5='p',mod3='OH')
    sc.note('trna','Terminal pairing view only; the tRNA 5′ acceptor-stem nick is not positioned in this linear view.')
    return sc


def ligated_scene():
    sc=Scene();parts=[Segment('Y5 DNA',Y5_DNA),role('matched Y5 terminal RNA dinucleotide',length_bp=2),
                       role('mature tRNA RNA body'),role('discriminator RNA',length_bp=1),Segment('RNA CCA','CCA',note='RNA'),Y3,DDC]
    sc.strand('ligated',parts,label='both nicks sealed; tRNA fold not laid out',mod3='/3ddC/')
    sc.junction('ligated','matched Y5 terminal RNA dinucleotide','mature tRNA RNA body','Rnl2 ligation')
    sc.junction('ligated','RNA CCA',Y3.name,'Rnl2 ligation')
    return sc


def cdna_scene():
    sc=Scene()
    sc.strand('cdna',[role('kit RTP-derived DNA, extent unassigned',inferred=True),
                       role('Y3-derived DNA copy, RTP boundary unassigned',inferred=True),
                       role('tRNA DNA copy, RT extent unassigned',inferred=True),
                       role('Y5-derived DNA copy, extent unassigned',inferred=True)],label='INFERRED — first-strand cDNA')
    # RT does not consume or hydrolyze the RNA template.
    sc.strand('rna',[role('original Y-adapted tRNA population; hybrid association not positioned')],label='retained RNA/DNA template with physical terminal ddC; extent unassigned')
    return sc


def library_scene():
    parts=[role('kit P5/read-entry DNA',inferred=True),role('Y5/tRNA/Y3-derived cDNA',inferred=True),
           role('sample index, width unassigned',molecular_feature=feature('yamat_index','sample_index','unknown'),inferred=True),
           role('kit P7 DNA',inferred=True)]
    return Scene.duplex(parts,label='INFERRED — selected native kit PCR library',
                         unpaired=tuple(s.name+"'" for s in parts))


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    t=g.input('trna','total RNA including aminoacylated mature tRNA',M.BIOLOGICAL,state=trna_scene(True),primary=True)
    alk=g.input('alkaline','Tris pH9;37°C40min',M.CONDITION)
    free=g.product('deacylated','deacylated RNA',M.MOLECULAR,state=trna_scene())
    aa=g.product('released_amino_acids','hydrolytically released amino acids',M.MOLECULAR)
    g.react('deacylate','Remove aminoacyl groups',K.CHEMICAL_CONVERSION,((t,R.TEMPLATE),(alk,R.CONDITION)),(free,aa))
    etoh=g.input('ethanol','ethanol precipitation',M.SELECTION)
    recovered=g.product('recovered','ethanol-recovered deacylated RNA',M.MOLECULAR,state=free.state)
    sol=g.product('precipitation_supernatant','discarded precipitation supernatant; released amino acids/reagents',M.WASTE)
    g.react('precipitate','Ethanol precipitation after deacylation',K.CLEANUP,
             ((free,R.TEMPLATE),(aa,R.TEMPLATE),(etoh,R.SELECTION)),(recovered,sol))
    g.discard(sol,'supernatant physically removed during RNA recovery; component-specific recovery not quantified')
    adapter=g.input('adapter','four Y5 DNA/RNA variants plus 5P/3ddC Y3',M.ADAPTER,state=adapter_set_scene())
    anneal=g.input('annealing','90°C2min then37°C15min with annealing buffer',M.CONDITION)
    complex_=g.product('complex','noncovalent Y-adapter/tRNA complex',M.MOLECULAR,state=annealed_scene())
    g.react('anneal','Anneal discriminator-matched Y adapters',K.HYBRIDIZATION,
             ((recovered,R.TEMPLATE),(adapter,R.STRUCTURAL),(anneal,R.CONDITION)),(complex_,))
    ligase=g.input('rnl2','T4 RNA ligase2;37°C1h then4°Covernight',M.CATALYST)
    sealed=g.product('sealed','Y5–tRNA–Y3 covalent strand; native modifications retained',M.MOLECULAR,state=ligated_scene())
    g.react('seal','Seal both RNA/RNA and RNA/DNA nicks',K.NICK_SEALING,
             ((complex_,R.TEMPLATE),(ligase,R.CATALYST)),(sealed,),
             geometry=AdapterLigationGeometry('unavailable','unavailable',(annealed_scene(),),
              'Source-known Y-adapter stem/terminal NCCA pairing is drawn; the complete acceptor fold and both spatial nick geometries are not reconstructed. This is a model boundary, not source absence.'))
    rtp=g.input('rtp','TruSeq Small RNA kit RTP; native boundary unassigned',M.OLIGO,state=molecule('RTP',(role('native RTP primer'),)))
    rt=g.input('ssiii','SuperScriptIII;55°C60min after70°C2min primer annealing',M.CATALYST)
    first=g.product('first_strand','INFERRED — cDNA with retained RNA template',M.MOLECULAR,state=cdna_scene(),inferred=True)
    g.react('rt','Reverse-transcribe ligated tRNA',K.REVERSE_TRANSCRIPTION,
             ((sealed,R.TEMPLATE),(rtp,R.STRUCTURAL),(rt,R.CATALYST)),(first,))
    f=g.input('forward','TruSeq Small RNA forward library primer; kit member unassigned',M.OLIGO)
    r=g.input('reverse','TruSeq Small RNA indexed reverse library primer; kit member unassigned',M.OLIGO)
    pol=g.input('pml','TruSeq Small RNA PML PCR enzyme mix;11cycles',M.CATALYST)
    amplified=g.product('amplified','INFERRED — native indexed PCR product',M.MOLECULAR,state=library_scene(),inferred=True)
    residual=g.product('residual_templates','original RNA and first-strand templates after PCR; physical recovery unassigned',M.MOLECULAR,state=first.state)
    g.react('pcr','INFERRED — eleven-cycle library PCR on actual first-strand input',K.PCR,
             ((first,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(pol,R.CATALYST)),(amplified,residual),
             geometry=PCRGeometry.unavailable(first.state,reason='The inspected defining Methods names the library kit primers without assigning native member sequences or first-strand boundaries; assay-control RT-PCR primers are not substituted.'))
    gel=g.input('gel','8% native PAGE; excise mature-tRNA amplified region',M.SELECTION)
    library=g.product('library','INFERRED — native-PAGE-purified mature-tRNA library',M.LIBRARY,state=amplified.state,inferred=True)
    outside=g.product('gel_remainder','material outside excised native-PAGE region',M.WASTE)
    residual_after=g.product('residual_after_selection','residual-template recovery through selected slice unresolved',M.MOLECULAR,state=residual.state)
    g.react('page','Gel-purify mature-tRNA amplified region',K.SELECTION,
             ((amplified,R.TEMPLATE),(residual,R.TEMPLATE),(gel,R.SELECTION)),(library,outside,residual_after))
    g.discard(outside,'outside physically selected gel region; no assignment of a specific RNA class solely from gel selection')
    g.mark_terminal(Terminal(residual_after.id,Disposition.UNRESOLVED,'Source does not assign component-specific recovery of original RNA/cDNA templates.'))
    g.sequence(library,SequencingPlan(acquired_roles=ACQUIRED_ROLES,read_lengths=READ_LENGTHS,limitation=SEQUENCING_ENDING))
    g.declare_reviewed_scope(('trna','alkaline','deacylate','ethanol','precipitate','adapter','annealing','anneal',
                              'rnl2','seal','rtp','ssiii','rt','forward','reverse','pml','pcr','gel','page','library'),SOURCE)
    return g.require_publishable()
