"""Published QuRIE-seq branches; source statements do not imply exact products."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
from chemdraw import Construct, Scene, Segment, feature, protocol_graph_panels
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S,
    SequencingPlan, Terminal, Disposition, SourceDiscrepancy, SourceClaim, AnnotationSlot)
from protocol_views import render_annotations
from page import head, info

PAGE_FILE='qurie.html'
SOURCE='https://doi.org/10.1016/j.crmeth.2021.100070; published STAR Methods and Supplementary Table S3'
PE2_N6_AS_WRITTEN=il.PE_READ2_ARM+'NNNNN'
PE2=il.P5+il.PE_READ2_PRIMER[1:23]
PE1_TAIL=il.TRUSEQ_READ1[3:19]
CUSTOM_READ1=il.PE_READ2_ARM[2:]
CUSTOM_INDEX=il.INDEX2_PRIMER_RC[:30]
CUSTOM_READ2=il.TRUSEQ_READ1[3:]
RNA_IVT_C=37
PROTEIN_IVT_C=50

def pe1():
    return Construct([Segment('P7',il.P7),
        Segment('six-base library index','N'*6,placeholder=True,
                feature=feature('qurie.sample','sample_index','unknown')),
        Segment('PE1 annealing arm',PE1_TAIL)],name='Source PE1 index-primer pattern')

def random_primer():
    return Construct([Segment('PE2-N6 arm',il.PE_READ2_ARM),
        Segment('five Ns printed in Table S3','N'*5,placeholder=True)],
        name='PE2-N6 as printed; source name/length discrepancy unresolved')

def _one(g,key,label,kind,template,reagent,role):
    reagent_kind={R.CATALYST:K.CATALYST,R.CONDITION:K.CONDITION,R.STRUCTURAL:K.OLIGO}[role]
    r=g.input(key+'_reagent',reagent,reagent_kind)
    out=g.product(key+'_product',label,K.MOLECULAR)
    g.react(key,label,kind,((template,R.TEMPLATE),(r,role)),(out,))
    return out

def _cleanup(g,key,template,reagent):
    selection=g.input(key+'_beads',reagent,K.SELECTION)
    kept=g.product(key+'_kept',key.replace('_',' ')+' retained material',K.MOLECULAR)
    waste=g.product(key+'_waste',key.replace('_',' ')+' unretained material/reagents',K.WASTE)
    g.react(key,key.replace('_',' '),S.CLEANUP,
            ((template,R.TEMPLATE),(selection,R.SELECTION)),(kept,waste))
    g.discard(waste,'Removed during this source-declared purification')
    return kept

def _branch_finish(g,prefix,template,p1,p2):
    # Source explicitly consumes a qPCR aliquot to determine the final PCR cycles.
    reserve=g.product(prefix+'_pcr_input','Remainder for '+prefix+' library PCR',K.MOLECULAR)
    qc=g.product(prefix+'_qpcr','qPCR calibration aliquot for '+prefix+' PCR cycle selection',K.QC)
    g.react(prefix+'_qpcr_split','Reserve aliquot for qPCR cycle determination',S.SPLIT,
            ((template,R.TEMPLATE),),(reserve,qc))
    g.mark_terminal(Terminal(qc.id,Disposition.QC))
    enzyme=g.input(prefix+'_kapa','KAPA HiFi HotStart PCR mix; cycles selected by qPCR',K.CATALYST)
    amplified=g.product(prefix+'_amplified',prefix+' PE1/PE2-amplified library material',K.MOLECULAR)
    g.react(prefix+'_pcr','Amplify '+prefix+' branch with PE1/PE2',S.PCR,
        ((reserve,R.TEMPLATE),(p1,R.STRUCTURAL),(p2,R.STRUCTURAL),(enzyme,R.CATALYST)),(amplified,),
        geometry=PCRGeometry.unavailable(Scene(),reason='Full branch templates/bead-promoter and RT-primer geometry remain unassigned in this capped reading; individual source oligos do not determine the final molecule.'))
    ratio='0.7× AMPure XP' if prefix=='rna' else '1.2× AMPure XP'
    selection=g.input(prefix+'_final_beads',ratio,K.SELECTION)
    lib=g.product(prefix+'_library',prefix+' QuRIE-seq library',K.LIBRARY)
    waste=g.product(prefix+'_final_waste','Unretained '+prefix+' PCR/cleanup material',K.WASTE)
    g.react(prefix+'_final_cleanup','Purify '+prefix+' library',S.CLEANUP,
            ((amplified,R.TEMPLATE),(selection,R.SELECTION)),(lib,waste))
    g.discard(waste,'Removed with final library bead cleanup')
    g.sequence(lib,SequencingPlan(acquisition='unspecified',limitation='Source states NextSeq500 and lists custom sequencing primers in Table S3; exact run lengths, library strand and primer/cycle geometry are not reconstructed in this pass.'))

def protocol_graph():
    g=ProtocolGraph('QuRIE-seq — BJAB RNA and intracellular-protein branches',authority=SOURCE)
    cells=g.input('cells','BJAB cells after the source experimental stimulation',K.BIOLOGICAL,primary=True)
    fixed=_one(g,'fix','DSP/SPDP reversible fixation',S.CROSS_LINKING,cells,'DSP/SPDP in PBS',R.CONDITION)
    permeable=_one(g,'permeabilize','Tris quench and Triton permeabilization (combined source treatment)',S.PERMEABILIZATION,
                    fixed,'100mM Tris-HCl,150mM NaCl,0.1% TritonX100',R.CONDITION)
    conjugate=g.input('antibody_tags','DNA-conjugated antibodies prepared with DBCO-S-S-NHS and 5′ azide oligos',K.OLIGO)
    stained=g.product('stained','Cells immunostained with intracellular/surface antibody-DNA panel',K.MOLECULAR)
    g.react('stain','Block and stain cells',S.BINDING,((permeable,R.TEMPLATE),(conjugate,R.STRUCTURAL)),(stained,))
    washed=_cleanup(g,'cell_wash',stained,'Two blocking-buffer washes after staining; recover singlets')
    reduced=_one(g,'reduce','DTT reduces reversible cell cross-links in droplets',S.CHEMICAL_CONVERSION,
                 washed,'DTT in droplet RT mix',R.CONDITION)
    bead=g.input('bead_primers','Source inDrop-based primer-loaded gel beads; exact revision/design unresolved',K.OLIGO)
    uv=g.input('uv','8min UV exposure at10mW/cm²',K.CONDITION)
    primers=g.product('released_primers','Released poly(dT) cell-barcoding primers',K.OLIGO)
    g.react('release_primers','Photorelease bead primers',S.CHEMICAL_CONVERSION,
            ((bead,R.TEMPLATE),(uv,R.CONDITION)),(primers,))
    rt=g.input('superscript','SuperScript III droplet RT mixture',K.CATALYST)
    copied=g.product('barcoded_mix','Barcoded RNA-derived cDNA and antibody-derived DNA products',K.MOLECULAR)
    g.react('barcode','Prime poly(A)/poly(dA) and copy in droplets:50°C2h,70°C15min',S.REVERSE_TRANSCRIPTION,
            ((reduced,R.TEMPLATE),(primers,R.STRUCTURAL),(rt,R.CATALYST)),(copied,))
    recovered=_cleanup(g,'recover',copied,'PFO emulsion break; pellet debris; separate aqueous phase and Spin-X purification')
    exo=_one(g,'exo','Exonuclease I treatment:37°C30min',S.DIGESTION,recovered,'Exonuclease I in FastDigest buffer',R.CATALYST)
    fractionator=g.input('fractionator','0.6× AMPure XP size separation',K.SELECTION)
    rna=g.product('rna_fraction','Long RNA-derived material retained on beads',K.MOLECULAR)
    protein=g.product('protein_fraction','Short antibody-derived material retained in supernatant',K.MOLECULAR)
    g.react('fractionate','Retain both RNA and protein size fractions',S.SPLIT,
            ((exo,R.TEMPLATE),(fractionator,R.SELECTION)),(rna,protein))
    p1=g.input('pe1','PE1 with six-base sample/library index',K.OLIGO,state=pe1())
    p2=g.input('pe2','PE2 primer',K.OLIGO,state=Construct([Segment('PE2',PE2)]))
    ds=_one(g,'rna_sss','RNA-branch second-strand synthesis:16°C1.5h,65°C20min',S.EXTENSION,rna,'NEBNext second-strand synthesis module',R.CATALYST)
    rna_ivt=_one(g,'rna_ivt',f'RNA-branch T7 IVT:{RNA_IVT_C}°C15h',S.TRANSCRIPTION,ds,'HiScribe T7 mixture with ATP,CTP,GTP,UTP',R.CATALYST)
    rna_clean=_cleanup(g,'rna_ivt_cleanup',rna_ivt,'1.3× AMPure XP')
    fragmented=_one(g,'rna_fragment','Fragment amplified RNA:70°C3min',S.FRAGMENTATION,rna_clean,'RNA fragmentation reagents; STOP mix on ice',R.CONDITION)
    rna_fragment_clean=_cleanup(g,'rna_fragment_cleanup',fragmented,'Source AMPure/STOP mixture')
    random=g.input('pe2_n6','PE2-N6 as printed; effective N count unresolved',K.OLIGO,state=random_primer())
    prime=g.input('primescript','PrimeScript RT mixture',K.CATALYST)
    rna_cdna=g.product('rna_cdna','Random-primed cDNA from amplified RNA',K.MOLECULAR)
    g.react('rna_rt','Random-primer RT of RNA branch',S.REVERSE_TRANSCRIPTION,
        ((rna_fragment_clean,R.TEMPLATE),(random,R.STRUCTURAL),(prime,R.CATALYST)),(rna_cdna,))
    rna_rt_clean=_cleanup(g,'rna_rt_cleanup',rna_cdna,'1.2× AMPure XP')
    _branch_finish(g,'rna',rna_rt_clean,p1,p2)
    protein_ivt=_one(g,'protein_ivt',f'Protein-branch IVT:{PROTEIN_IVT_C}°C15h as published',S.TRANSCRIPTION,
                     protein,'HiScribe T7 mixture with ATP,CTP,GTP,UTP',R.CATALYST)
    protein_clean=_cleanup(g,'protein_ivt_cleanup',protein_ivt,'1.5× AMPure XP; elute20µL')
    backup=g.product('protein_backup','10µL post-IVT protein backup at−80°C',K.MOLECULAR)
    protein_rt_input=g.product('protein_rt_input','10µL post-IVT protein material for RT',K.MOLECULAR)
    g.react('protein_backup_split','Retain post-IVT protein backup',S.SPLIT,
            ((protein_clean,R.TEMPLATE),),(backup,protein_rt_input))
    g.mark_terminal(Terminal(backup.id,Disposition.STORED,'Source explicitly stores half of the20µL eluate at−80°C'))
    specific=g.input('protein_rt_primers','PE2-NNNN-Next1 and PE2-NNNN-BioHash2; full designs unassigned',K.OLIGO)
    protein_cdna=g.product('protein_cdna','Sequence-primed cDNA of amplified protein reporters',K.MOLECULAR)
    g.react('protein_rt','Protein-branch sequence-specific RT',S.REVERSE_TRANSCRIPTION,
            ((protein_rt_input,R.TEMPLATE),(specific,R.STRUCTURAL),(prime,R.CATALYST)),(protein_cdna,))
    protein_rt_clean=_cleanup(g,'protein_rt_cleanup',protein_cdna,'1.5× AMPure XP')
    _branch_finish(g,'protein',protein_rt_clean,p1,p2)
    g.declare_reviewed_scope((
        'cells','fix_reagent','fix_product','fix','permeabilize_reagent','permeabilize_product','permeabilize',
        'antibody_tags','stained','stain','cell_wash_beads','cell_wash_kept','cell_wash_waste','cell_wash',
        'reduce_reagent','reduce_product','reduce','bead_primers','uv','released_primers','release_primers',
        'superscript','barcoded_mix','barcode','recover_beads','recover_kept','recover_waste','recover',
        'exo_reagent','exo_product','exo','fractionator','rna_fraction','protein_fraction','fractionate','pe1','pe2',
        'rna_sss_reagent','rna_sss_product','rna_sss','rna_ivt_reagent','rna_ivt_product','rna_ivt',
        'rna_ivt_cleanup_beads','rna_ivt_cleanup_kept','rna_ivt_cleanup_waste','rna_ivt_cleanup',
        'rna_fragment_reagent','rna_fragment_product','rna_fragment','rna_fragment_cleanup_beads',
        'rna_fragment_cleanup_kept','rna_fragment_cleanup_waste','rna_fragment_cleanup','pe2_n6','primescript',
        'rna_cdna','rna_rt','rna_rt_cleanup_beads','rna_rt_cleanup_kept','rna_rt_cleanup_waste','rna_rt_cleanup',
        'rna_pcr_input','rna_qpcr','rna_qpcr_split','rna_kapa','rna_amplified','rna_pcr','rna_final_beads',
        'rna_library','rna_final_waste','rna_final_cleanup','protein_ivt_reagent','protein_ivt_product','protein_ivt',
        'protein_ivt_cleanup_beads','protein_ivt_cleanup_kept','protein_ivt_cleanup_waste','protein_ivt_cleanup',
        'protein_backup','protein_rt_input','protein_backup_split','protein_rt_primers','protein_cdna','protein_rt',
        'protein_rt_cleanup_beads','protein_rt_cleanup_kept','protein_rt_cleanup_waste','protein_rt_cleanup',
        'protein_pcr_input','protein_qpcr','protein_qpcr_split','protein_kapa','protein_amplified','protein_pcr',
        'protein_final_beads','protein_library','protein_final_waste','protein_final_cleanup'),SOURCE)
    g.annotate(SourceDiscrepancy(target='pe2_n6',title='PE2-N6 name/length and printed sequence differ',
        claims=(SourceClaim('Main RNA-library Methods and TableS3 name/length','Random hexamer/PE2-N6; Table length38bp.'),
                SourceClaim('Published rendered TableS3,PDF14,printed13','32-base arm followed by five printed Ns:37 characters.')),
        model_choice='Preserve five Ns as written; do not claim effective priming geometry.',
        evidence='Original rendered table inspected; this is not a text-extraction-only discrepancy.',
        uncertainty='Actual manufactured N count/order is unresolved; source follow-up remains open.',slot=AnnotationSlot.BEFORE_WORKFLOW))
    return g

def render_page():
    g=protocol_graph()
    return '\n'.join([head('QuRIE-seq — bounded two-branch model'),'<div class="wrap"><h1>QuRIE-seq</h1>',
        '<p class="research-notes"><a href="01_qurie_seq.html">Research notes</a></p>',
        info('Published BJAB branch. Complete bead/tag/promoter strands and read layout remain unassigned. Protein IVT50°C and RNA IVT37°C are separate published conditions; no biochemical efficiency is predicted.'),
        render_annotations(g,AnnotationSlot.BEFORE_WORKFLOW),protocol_graph_panels(g,cls='long'),'</div>'])

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Both size fractions proceed to library branches',len(g.sequencing),2)
    c('Post-IVT protein backup retained',g.terminals['protein_backup'].disposition,Disposition.STORED)
    c('Both qPCR aliquots accounted',sum(t.disposition is Disposition.QC for t in g.terminals.values()),2)
    c('Source branch-specific IVT temperatures retained',(RNA_IVT_C,PROTEIN_IVT_C),(37,50))
    from storage import sources_root
    f=sources_root()/Path(__file__).resolve().parents[1].name/'j.crmeth.2021.100070_supp_mmc1.pdf.txt'
    if f.exists():
        text=''.join(f.read_text().split())
        for seq in (PE2,CUSTOM_READ1,CUSTOM_INDEX,CUSTOM_READ2):c('Published TableS3 transcription '+seq,seq in text)
        row=f.read_text().split('PE2-N6',1)[1].split('PE1',1)[0]
        c('Source PE2-N6 arm and printed five-N line',il.PE_READ2_ARM in row and 'NNNNN' in row.splitlines())
    else:c.skip('Published primer transcription','External supplement unavailable')
    c.report(exit_on_fail=True)
