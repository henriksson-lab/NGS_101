"""Protocol-owned representative two-step 4C-seq declaration."""
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Segment, Scene, circle_rows, feature, oligo, revcomp
from circular import CircularizationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import Drawing, molecule, role

TITLE='4C-seq — two-step inverse-PCR contact libraries'
NOTES='01_4c_seq.html'
SOURCE='<a href="https://doi.org/10.1016/j.ymeth.2019.07.014">Krijger et al., Methods (2020)</a>, original author-hosted protocol, Sections 3–6 and Table 1.'
SUMMARY='Primary restriction/proximity ligation creates 3C contacts. A second digest and dilute ligation make circles; outward viewpoint PCR followed by indexed universal PCR prepares sequencing libraries.'
CAVEAT=('INFERRED — the diagram represents the published DpnII/Csp6I example combination and the Table 1 index 1–12 primer architecture. These are explicit representative choices, not a universal enzyme pair or an assigned experiment. Viewpoint genomic arms, optional spacer and selected six-base index are unassigned. '
        'Other printed reverse primers 13–27 contain two extra physical bases that are not index cycles; those variants are not silently collapsed into this endpoint. '
        'Circle/contact assemblies are role-level representations; no genomic restriction-site positions or primer-specific circle geometry are claimed. Self-ligation and undigested products can amplify; their analytical exclusion is not physical removal.')
READING_TAIL='TACACGACGCTCTTCCGATCT'
NONREADING_TAIL='ACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
FORWARD2=il.P5+il.TRUSEQ_READ1[4:]
REVERSE2_ARM=il.TRUSEQ_READ2[:23]
INDEX=Segment('six-base sample i7 index','N'*6,'cbc',placeholder=True,feature=feature('sample_i7','sample_index','unknown'),inferred=True)
FORWARD1=Construct((Segment('published reading tail',READING_TAIL,'r1'),role('optional 0–5-base spacer; actual sequence/length unassigned'),role('viewpoint reading genomic arm')),name='reading inverse-PCR oligo')
REVERSE1=Construct((Segment('published non-reading tail',NONREADING_TAIL,'r2'),role('viewpoint non-reading genomic arm')),name='non-reading inverse-PCR oligo')
FIRST_AMPLICON=Construct((Segment('reading tail',READING_TAIL,'r1'),role('optional spacer unassigned'),role('viewpoint reading genomic arm'),role('captured contact DNA'),role('opposite viewpoint genomic region'),Segment('non-reading tail reverse complement',revcomp(NONREADING_TAIL),'r2')),name='INFERRED — locus-unassigned first inverse-PCR product')
FORWARD2_OLIGO=Construct((Segment('source second-round forward',FORWARD2,'p5'),),name='source universal PCR2 forward')
REVERSE2=Construct((Segment('P7',il.P7,'p7'),INDEX,Segment('source reverse primer arm',REVERSE2_ARM,'r2')),name='Table 1 index 1–12 reverse family; selected member unassigned')
FINAL_LIBRARY=Construct((Segment('source forward completion',FORWARD2,'p5',inferred=True),*FIRST_AMPLICON.segments[1:-1],Segment('complete read2-side reverse complement',revcomp(il.TRUSEQ_READ2),'r2',inferred=True),INDEX,Segment('P7 reverse complement',il.P7_RC,'p7',inferred=True)),name='INFERRED — representative single-index 4C library')
FINAL_CAPTION='INFERRED — two-stage completion with Table 1 index 1–12 family. Read 1 enters optional spacer, viewpoint genomic primer and captured reporter. No particular viewpoint/index/run is assigned.'
SEQ_PRIMERS=()
SEQUENCING_ENDING='Section 6 describes single-read 50-base HiSeq or 75-base MiniSeq/NextSeq runs as generally sufficient, depending on viewpoint primer placement. The index is six bases in Table 1.'
SEQUENCING_UNAVAILABLE='Exact loaded sequencing oligos and a particular acquired run are unassigned; source PCR-site architecture does not itself establish the loaded primer.'
READ_SPECS=(sp.read_spec('Read 1','50 or 75 source examples','optional spacer, viewpoint genomic arm, then captured reporter','include'),sp.read_spec('Index 1 (i7)','6','sample identity; selected index unassigned','demultiplex'))

def oligos():
    return ('Reading inverse-PCR tail: 5′-'+READING_TAIL+'-[optional 0–5-base spacer]-[viewpoint arm]-3′',
            'Non-reading inverse-PCR tail: 5′-'+NONREADING_TAIL+'-[viewpoint arm]-3′',
            'Table 1 second-round universal forward: 5′-'+FORWARD2+'-3′',
            'Table 1 index 1–12 reverse family: 5′-'+il.P7+'-[6-base index]-'+REVERSE2_ARM+'-3′')

def separate_classes(label, names):
    sc=Scene()
    for i,name in enumerate(names): sc.strand(str(i),[role(name)],label=name)
    return sc

def contact_scene():
    return molecule('INFERRED — representative 3C contact',(role('viewpoint fragment'),role('reporter fragment')),junctions=(('viewpoint fragment','reporter fragment','primary proximity-ligation junction'),))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority='Krijger original Sections 3–6/Table 1; representative DpnII/Csp6I and index 1–12 branch')
    fixed=g.input('fixed_chromatin','formaldehyde-fixed chromatin',MaterialKind.BIOLOGICAL,state=separate_classes('fixed chromatin',('viewpoint region','proximally crosslinked chromatin')),primary=True)
    re1=g.input('primary_restriction','DpnII, representative source example',MaterialKind.CATALYST)
    dig=g.product('digested','primary restriction fragments',MaterialKind.MOLECULAR,state=separate_classes('digested chromatin',('viewpoint restriction fragment','other restriction fragments')))
    g.react('primary_digest','Digest fixed chromatin with chosen primary enzyme',StepKind.DIGESTION,((fixed,InputRole.TEMPLATE),(re1,InputRole.CATALYST)),(dig,))
    lig=g.input('ligase','source T4 DNA ligase',MaterialKind.CATALYST)
    contact=g.product('contact','representative proximity-ligated 3C DNA',MaterialKind.MOLECULAR,state=contact_scene(),inferred=True)
    g.react('proximity_ligation','INFERRED — representative contact after source proximity ligation',StepKind.PROXIMITY_LIGATION,((dig,InputRole.TEMPLATE),(lig,InputRole.CATALYST)),(contact,))
    reversal=g.input('crosslink_reversal','source Proteinase K crosslink-reversal condition',MaterialKind.CONDITION)
    free=g.product('deproteinized','de-crosslinked 3C DNA',MaterialKind.MOLECULAR,state=contact.state)
    g.react('reverse_crosslinks','Reverse crosslinks',StepKind.DEPROTEINIZATION,((contact,InputRole.TEMPLATE),(reversal,InputRole.CONDITION)),(free,))
    beads=g.input('template_cleanup','source P-bead template purification',MaterialKind.SELECTION)
    clean=g.product('clean_3c','purified 3C DNA',MaterialKind.MOLECULAR,state=free.state)
    wash=g.product('template_wash','template purification wash material',MaterialKind.WASTE)
    g.react('purify_3c','Purify de-crosslinked 3C DNA',StepKind.CLEANUP,((free,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(clean,wash));g.discard(wash,'Source P-bead procedure discards washes before DNA elution.')
    re2=g.input('secondary_restriction','Csp6I, representative source example',MaterialKind.CATALYST)
    trimmed=g.product('secondary_fragments','secondary-digested 3C DNA',MaterialKind.MOLECULAR,state=contact.state)
    g.react('secondary_digest','Trim 3C DNA with chosen second enzyme',StepKind.DIGESTION,((clean,InputRole.TEMPLATE),(re2,InputRole.CATALYST)),(trimmed,))
    circle_parts=Construct((role('viewpoint-side region'),role('primary contact junction'),role('reporter region'),role('secondary-ligation junction')),name='INFERRED — representative 4C circle')
    circ=g.product('circle','INFERRED — representative recircularized 4C template',MaterialKind.MOLECULAR,state=Drawing(tuple(circle_rows(circle_parts))),inferred=True)
    g.react('recircularize','INFERRED — representative circle after dilute ligation',StepKind.CIRCULARIZATION,((trimmed,InputRole.TEMPLATE),(lig,InputRole.CATALYST)),(circ,),geometry=CircularizationGeometry.known('double-strand intramolecular ligation',statement='Source Section 4.7 dilute-ligates secondary fragments into 4C circles; genomic strand-end bases unassigned.'))
    ready=g.product('clean_circle','purified 4C circle template',MaterialKind.MOLECULAR,state=circ.state,inferred=True)
    wash2=g.product('circle_wash','circle purification wash material',MaterialKind.WASTE)
    g.react('purify_circle','Purify 4C template',StepKind.CLEANUP,((circ,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(ready,wash2));g.discard(wash2,'Source Section 4.8 repeats template purification.')
    f1=g.input('reading_primer','source reading inverse-PCR primer; genomic arm unassigned',MaterialKind.OLIGO,state=FORWARD1)
    r1=g.input('nonreading_primer','source non-reading inverse-PCR primer; genomic arm unassigned',MaterialKind.OLIGO,state=REVERSE1)
    pol=g.input('polymerase','Expand Long Template Polymerase mix',MaterialKind.CATALYST)
    amp=g.product('first_amplicon',FIRST_AMPLICON.name,MaterialKind.MOLECULAR,state=FIRST_AMPLICON,inferred=True)
    off=g.product('amplified_byproducts','source observed self-ligation/undigested amplicon classes',MaterialKind.MOLECULAR,state=separate_classes('byproduct amplicon classes',('viewpoint self-ligation amplicon','undigested-neighbor amplicon')))
    g.react('inverse_pcr','INFERRED — outward viewpoint PCR; byproducts can also amplify',StepKind.PCR,((ready,InputRole.TEMPLATE),(f1,InputRole.STRUCTURAL),(r1,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amp,off),geometry=PCRGeometry.unavailable(circle_parts,reason='Actual genomic viewpoint arms/sites are unassigned; no exact outward circle binding geometry declared.'))
    g.mark_terminal(Terminal(off.id,Disposition.UNRESOLVED,'Original Figure 3/note 12 and Section 7 describe amplified self-ligation/undigested products. Their analytical exclusion from metrics is not a physical discard; this representative branch does not trace their complete recovery.'))
    ampure=g.input('ampure','source 0.8x AMPure first-PCR cleanup',MaterialKind.SELECTION)
    ampclean=g.product('first_clean','purified first-PCR amplicon',MaterialKind.MOLECULAR,state=FIRST_AMPLICON,inferred=True)
    pw=g.product('first_pcr_wash','first-PCR cleanup wash material',MaterialKind.WASTE)
    g.react('cleanup_first_pcr','Clean first-PCR products',StepKind.CLEANUP,((amp,InputRole.TEMPLATE),(ampure,InputRole.SELECTION)),(ampclean,pw));g.discard(pw,'Source Section 5.2 performs AMPure cleanup; no selective removal of every byproduct asserted.')
    f2=g.input('forward2','source universal second-PCR forward',MaterialKind.OLIGO,state=FORWARD2_OLIGO,inferred=True)
    r2=g.input('reverse2','source index 1–12 second-PCR reverse family',MaterialKind.OLIGO,state=REVERSE2,inferred=True)
    lib=g.product('library',FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True)
    geom=PCRGeometry.exact(FIRST_AMPLICON,PCRPrimer('second forward',FORWARD2,anneal_nt=len(READING_TAIL)),PCRPrimer('second reverse family',REVERSE2.top(),anneal_nt=20),statement='Source universal tails/Table 1; exact inward binding recomputed on locus-unassigned first amplicon, selected index unknown.')
    g.react('second_pcr','INFERRED — source second PCR completes representative index 1–12 library',StepKind.PCR,((ampclean,InputRole.TEMPLATE),(f2,InputRole.STRUCTURAL),(r2,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(lib,),geometry=geom)
    column=g.input('final_cleanup','source QIAquick PCR purification',MaterialKind.SELECTION)
    out=g.product('purified_library',FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True)
    cw=g.product('final_wash','column flow-through/washes',MaterialKind.WASTE)
    g.react('purify_library','Purify completed library',StepKind.CLEANUP,((lib,InputRole.TEMPLATE),(column,InputRole.SELECTION)),(out,cw));g.discard(cw,'Source second-PCR purification discards column flow-through/washes.')
    g.sequence(out,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(('fixed_chromatin','primary_restriction','primary_digest','ligase','proximity_ligation','crosslink_reversal','reverse_crosslinks','template_cleanup','purify_3c','secondary_restriction','secondary_digest','recircularize','purify_circle','reading_primer','nonreading_primer','polymerase','inverse_pcr','ampure','cleanup_first_pcr','forward2','reverse2','second_pcr','library','final_cleanup','purify_library','purified_library'),g.authority)
    return g.require_publishable()
