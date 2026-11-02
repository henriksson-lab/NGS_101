"""Optional authoring checks: python3 -m animation.selftest [--native]."""
from dataclasses import replace
import json
import subprocess
import tempfile
from pathlib import Path
import sys
import unittest
from .model import Movie, Event
from .examples import reference, smartseq


class SubstrateSimulationTests(unittest.TestCase):
    def test_pcr_enzymes_arrive_when_extension_starts(self):
        from .smartseq2 import movie
        from .evaluate import evaluate
        m=movie()
        for key in ['index','index_copy','index_duplex']:
            shot=next(s for s in m.shots if s.id==key)
            for phase in [0,.3,.54]:self.assertFalse(evaluate(shot,shot.establish+shot.action*phase)[1])
            self.assertEqual(len(evaluate(shot,shot.establish+shot.action*.7)[1]),len([e for e in shot.events if e.kind=='extend']))
        switch=next(s for s in m.shots if s.id=='switch')
        self.assertEqual(len(evaluate(switch,0)[1]),1)

    def test_all_ends_extend_without_primer_role_or_fixed_actor_count(self):
        from .simulator import Simulator,Polymerase
        s=Simulator();s.supply_nt('DNA','ACGT',source='test');roots=[];old=[]
        for label in ['one','two']:
            a=s.add(label,'AAAACCGG',source='test')
            b=s.add(label+'_other','TTTCC',source='test')
            s.anneal(label+'_other',label,binding=2,source='test')
            roots.append(a[0]);old.append((a,b))
        enzyme=Polymerase('test polymerase',('DNA',),'DNA')
        sites=set(s.extension_sites(enzyme,roots));self.assertEqual(len(sites),4)
        operations=s.extend_available(enzyme,roots,source='test')
        self.assertEqual({op.details['starting_3prime']for op in operations},sites)
        self.assertFalse(s.extension_sites(enzyme,roots))
        for a,b in old:
            self.assertEqual(s.chain(a[0])[:len(a)],a)
            self.assertEqual(s.chain(b[0])[:len(b)],b)
            self.assertGreater(len(s.chain(a[0])),len(a))
            self.assertGreater(len(s.chain(b[0])),len(b))

    def test_vessel_scope_and_blocked_end_change_discovered_activity(self):
        from .simulator import Simulator,Polymerase
        s=Simulator();s.supply_nt('DNA','ACGT',source='test')
        for label in ['chosen','other']:
            s.add(label,'AAAACCGG',three='P'if label=='chosen'else'OH',source='test')
            s.add(label+'_oligo','TTTCC',source='test');s.anneal(label+'_oligo',label,binding=2,source='test')
        enzyme=Polymerase('test',('DNA',),'DNA');before=s.snapshot()
        operations=s.extend_available(enzyme,['chosen'],source='test')
        self.assertEqual(len(operations),1)
        self.assertEqual(s.sequence('chosen'),'AAAACCGG')
        for id in s.chain('other')+s.chain('other_oligo'):
            self.assertEqual(s.atoms[id],before.atoms[id])
            self.assertEqual(s.next.get(id),before.next.get(id))
        self.assertEqual(len(s.extension_sites(enzyme,['other'])),2)

    def test_failed_multi_end_reaction_rolls_back_all_extensions(self):
        from .simulator import Simulator,Polymerase
        s=Simulator();s.supply_nt('DNA','GT',source='test')
        s.add('substrate','AAAACCGG',source='test');s.add('oligo','TTTCC',source='test')
        s.anneal('oligo','substrate',binding=2,source='test')
        before=s.snapshot();history=len(s.history)
        with self.assertRaises(ValueError):
            s.extend_available(Polymerase('test',('DNA',),'DNA'),['substrate'],source='test')
        self.assertEqual(s.snapshot(),before)
        self.assertEqual(len(s.history),history)

    def test_rna_primer_requires_explicit_capability_and_retains_sugars(self):
        from .simulator import Simulator,Polymerase
        s=Simulator();s.supply_nt('DNA','ACGT',source='test')
        s.add('template','AAAACCGG',source='test');rna=s.add('primer','UUUCC',sugar='RNA',source='test')
        s.anneal('primer','template',binding=2,source='test')
        self.assertFalse(s.extension_sites(Polymerase('DNA primer only',('DNA',),'DNA'),['primer']))
        enzyme=Polymerase('explicit RNA primer profile',('DNA',),'DNA',primer_sugars=('DNA','RNA'))
        operations=s.extend_available(enzyme,['primer'],source='test');self.assertEqual(len(operations),1)
        self.assertEqual(s.chain('primer')[:len(rna)],rna)
        self.assertTrue(all(s.atoms[id].sugar=='RNA'for id in rna))
        self.assertTrue(all(s.atoms[id].sugar=='DNA'for id in operations[0].added))

    def test_generic_projection_keeps_rna_primer_and_new_dna_distinct(self):
        from .simulator import Simulator,Polymerase
        from .generated_views import strand
        from .simulation_projection import extensions
        from .simulation_adapters import generated
        from .model import Shot
        s=Simulator();s.supply_nt('DNA','ACGT',source='test')
        t=s.add('template','GGGCCAA',source='test');p=s.add('primer','GGGUU',sugar='RNA',source='test')
        s.anneal('primer','template',binding=2,source='test');before=s.snapshot()
        ops=s.extend_available(Polymerase('explicit RNA-primer profile',('DNA',),'DNA',primer_sugars=('DNA','RNA')),['primer'],source='test')
        shot=Shot('mixed','RNA primer extended with DNA',[strand('template',before,t,500,440,direction=-1),strand('primer',before,p,410,464)],[],'Fixture.',source='synthetic fixture')
        extensions(shot,before,s.snapshot(),ops);generated(shot,before);shot.validate()
        product=next(t for t in shot.strands if t.id=='primer')
        self.assertEqual(product.polymer,'mixed')
        self.assertEqual(product.sugars,[(0,len(p),'RNA'),(len(p),len(product.sequence),'DNA')])

    def test_smartseq2_projects_every_recessed_end_in_index_pcr(self):
        from .smartseq2 import movie
        m=movie()
        for key,count in [('index',4),('index_copy',4),('index_duplex',2)]:
            shot=next(s for s in m.shots if s.id==key)
            extensions=[e for e in shot.events if e.kind=='extend']
            self.assertEqual(len(extensions),count)
            self.assertEqual(len({e.strand for e in extensions}),count)
            self.assertTrue(all(e.molecular_operation['details']['starting_3prime']for e in extensions))

    def test_extension_follows_cut_template_without_a_product_recipe(self):
        from .simulator import Simulator,Polymerase
        from .model import complement
        for boundary in [0,3,5]:
            s=Simulator();s.supply_nt('DNA','ACGT',source='test')
            t=s.add('template','AAAACCGG',source='test');p=s.add('primer','CC',source='test')
            s.anneal('primer','template',binding=2,source='test')
            if boundary:s.cleave(t[boundary-1],source='test')
            original=s.snapshot();op=s.extend('primer',Polymerase('test DNA polymerase',('DNA',),'DNA'),source='test')
            self.assertEqual(s.sequence('primer'),'CC'+complement('AAAACCGG'[boundary:-2][::-1]))
            self.assertEqual(op.details['stop'],'physical_template_end')
            self.assertEqual(s.chain('primer')[:2],p)
            for id in op.added:
                self.assertIn(s.atoms[id].template,original.atoms)

    def test_failed_synthesis_does_not_leave_a_partial_product(self):
        from .simulator import Simulator,Polymerase
        s=Simulator();s.supply_nt('DNA','G',source='test')
        s.add('template','AAAACCGG',source='test');s.add('primer','CC',source='test')
        s.anneal('primer','template',binding=2,source='test');before=s.snapshot()
        with self.assertRaises(ValueError):s.extend('primer',Polymerase('test',('DNA',),'DNA'),source='test')
        self.assertEqual(s.snapshot(),before)

    def test_smartseq2_changes_propagate_from_input_through_pcr(self):
        from .smartseq2_simulation import simulate
        from .protocols import load
        m=load('smart-seq2__10.1038+nmeth.2639','smartseq2.py')
        body='ACGUGCUAGCUACGACUGGUACGUACCGUAUGCUGACGUAGCUACGUA'
        r=simulate(m,body)
        variant=simulate(m,body[:20]+('A'if body[20]!='A'else'C')+body[21:])
        self.assertNotEqual(r['final_plus'],variant['final_plus'])
        self.assertEqual(len(r['final_plus']),len(variant['final_plus']))
        shifted=simulate(m,body,cut_offsets=(8,34))
        self.assertEqual(len(shifted['final_plus']),len(r['final_plus'])+2)
        for result in [r,variant,shifted]:
            for name in ['top','bottom']:
                op=result[name+'_gap']
                self.assertEqual(len(op.added),9)
                self.assertEqual(op.details['stop'],'occupied_template')
            for name in ['index_f1_op','index_r1_op']:
                self.assertEqual(result[name].details['stop'],'physical_template_end')
            state=result['final_state']
            for a,b in zip(result['final_plus_ids'],reversed(result['final_minus_ids'])):
                self.assertEqual(state.paired[a],b)
            for key in ['transfer_left','transfer_right']:
                op=result[key]
                self.assertEqual(op.before.atoms,op.after.atoms)
                self.assertEqual(op.before.paired,op.after.paired)


class RNaseHTests(unittest.TestCase):
    def hybrid(self, sugar='RNA', partner_sugar='DNA', **kwargs):
        from .simulator import Simulator
        from .model import complement
        s=Simulator();seq='AUGCAGUACCGU' if sugar=='RNA' else 'ATGCAGTACCGT'
        for lo,hi,kind in kwargs.get('sugars',[]):
            if kind=='RNA':seq=seq[:lo]+seq[lo:hi].replace('T','U')+seq[hi:]
        r=s.add('substrate',seq,sugar=sugar,source='synthetic test',**kwargs)
        opposite=complement(seq[::-1])
        if partner_sugar=='RNA':opposite=opposite.replace('T','U')
        d=s.add('partner',opposite,sugar=partner_sugar,source='synthetic test')
        s.duplex('substrate','partner',source='synthetic test')
        return s,r,d

    def test_cuts_only_rna_and_retains_hybrid_and_atoms(self):
        s,r,d=self.hybrid();before=s.snapshot()
        op=s.rnase_h1([r[5]],source='selected H1 illustration')
        self.assertEqual(s.atoms,before.atoms)
        self.assertEqual(s.paired,before.paired)
        self.assertEqual(set(before.next.items())-set(s.next.items()),{(r[5],r[6])})
        self.assertEqual(s.three[r[5]],'OH');self.assertEqual(s.five[r[6]],'P')
        self.assertEqual(s.chain('partner'),d)
        self.assertFalse(s.released);self.assertFalse(op.added)

    def test_requires_rna_dna_hybrid_and_four_rna_footprint(self):
        for sugar,partner in [('DNA','DNA'),('RNA','RNA')]:
            s,r,d=self.hybrid(sugar,partner);before=s.snapshot()
            with self.assertRaises(ValueError):s.rnase_h1([r[5]],source='test')
            self.assertEqual(s.snapshot(),before)
        for span in [(5,6,'RNA'),(4,7,'RNA')]:
            s,r,d=self.hybrid('DNA',sugars=[span]);before=s.snapshot()
            with self.assertRaisesRegex(ValueError,'two consecutive RNA'):
                s.rnase_h1([r[5]],source='test')
            self.assertEqual(s.snapshot(),before)
        s,r,d=self.hybrid();s.unpair(r,source='test');before=s.snapshot()
        with self.assertRaisesRegex(ValueError,'hybrid'):s.rnase_h1([r[5]],source='test')
        self.assertEqual(s.snapshot(),before)

    def test_rechecks_connectivity_and_rolls_back_failed_batch(self):
        s,r,d=self.hybrid();before=s.snapshot()
        with self.assertRaises(ValueError):s.rnase_h1([r[5],r[6]],source='test')
        self.assertEqual(s.snapshot(),before)
        s.cleave(d[5],source='test');before=s.snapshot()
        with self.assertRaisesRegex(ValueError,'continuous antiparallel'):
            s.rnase_h1([r[5]],source='test')
        self.assertEqual(s.snapshot(),before)

    def test_movie_replays_real_bond_changes_and_roundtrips(self):
        from .examples import rnaseh
        from .molecular_frame import molecular_frame
        movie=rnaseh();shot=movie.shots[0]
        clone=Movie.from_dict(json.loads(json.dumps(movie.to_dict())))
        initial=molecular_frame(shot,0)
        for phase,removed in [(0,0),(.45,1),(1,2)]:
            t=shot.establish+shot.action*phase
            current=molecular_frame(shot,t)
            self.assertEqual(len(initial.backbone)-len(current.backbone),removed)
            self.assertEqual(current.atoms,initial.atoms)
            self.assertEqual(current.duplex,initial.duplex)
            restored=molecular_frame(clone.shots[0],t)
            self.assertEqual(json.dumps(current.atoms,sort_keys=True),json.dumps(restored.atoms,sort_keys=True))
            self.assertEqual((current.backbone,current.duplex,current.released),(restored.backbone,restored.duplex,restored.released))


class MolecularTests(unittest.TestCase):
    def test_copy_rejects_wrong_base(self):
        m=reference();m.shots[0].strands[1].sequence='A'*len(m.shots[0].strands[1].sequence)
        with self.assertRaisesRegex(ValueError,'complement'):m.validate()

    def test_copy_rejects_blocked_end(self):
        m=reference();m.shots[0].strands[1].three_prime='ddC'
        with self.assertRaisesRegex(ValueError,'extendable'):m.validate()

    def test_mapping_rejects_wrong_coordinates(self):
        m=reference();m.shots[0].events[0].template_start=999
        with self.assertRaisesRegex(ValueError,'leaves template'):m.validate()

    def test_unknown_is_not_exact(self):
        m=reference();m.shots[0].strands[0].schematic=True
        with self.assertRaisesRegex(ValueError,'unknown classes'):m.validate()

    def test_speech_extends_hold_not_action(self):
        shot=reference().shots[0];short=shot.timing();long=replace(shot,narration_seconds=50).timing()
        self.assertEqual(short[:2],long[:2]);self.assertEqual(long[2],51)

    def test_duration_does_not_truncate(self):
        shot=replace(reference().shots[0],narration_seconds=50,duration=10)
        with self.assertRaisesRegex(ValueError,'truncate'):shot.timing()

    def test_roundtrip_and_seek(self):
        m=smartseq();r=Movie.from_dict(json.loads(json.dumps(m.to_dict())))
        for t in [0,5,25,10,0,10000]:
            self.assertEqual(m.at(t)[0].id,r.at(t)[0].id)
        self.assertIn('not SMART-seq2',r.provenance['original_protocol'])
        self.assertTrue(all(s.schematic for shot in r.shots for s in shot.strands))

    def test_invalid_speed_and_timing(self):
        for speed in [0,-1,float('nan')]:
            with self.assertRaises(ValueError):reference().duration(speed=speed)
        with self.assertRaises(ValueError):reference().duration(wpm=0)

    def test_speech_estimate_uses_chosen_rate_at_export(self):
        m=reference();shot=m.shots[0];shot.narration=' '.join(['word']*50);shot.duration=12
        m.validate()
        self.assertEqual(m.duration(wpm=300,speed=1),12)
        with self.assertRaisesRegex(ValueError,'truncate'):m.duration(wpm=100)

    def test_version_is_checked(self):
        d=reference().to_dict();d['schema']='future/99'
        with self.assertRaisesRegex(ValueError,'schema'):Movie.from_dict(d)

    def test_ligation_not_exact_without_ends(self):
        m=reference();m.shots[0].strands[1].nicks=[8];m.shots[0].events=[Event('seal','daughter',bond=8,exact=True)]
        with self.assertRaisesRegex(ValueError,'two-end'):m.validate()


class ReactionTests(unittest.TestCase):
    def test_closure_checks_actual_end_modifications(self):
        from .protocols import crispr_mip
        m=crispr_mip();m.shots[-1].strands[0].five_prime='OH'
        with self.assertRaisesRegex(ValueError,'terminal modifications'):m.validate()

    def test_blocked_junction_cannot_ligate(self):
        from .protocols import smallseq
        m=smallseq();next(e for e in m.shots[1].events if e.kind=='seal').junction.acceptor_3prime='ddC'
        with self.assertRaisesRegex(ValueError,'3primeOH'):m.validate()

    def test_long_probe_backbone_is_retained_not_forced_to_pair(self):
        from .protocols import crispr_mip
        m=crispr_mip();shot=m.shots[1];event=shot.events[0]
        self.assertTrue(event.exact)
        self.assertEqual(event.primer_binding_length,19)
        self.assertGreater(event.primer_length,event.primer_binding_length)
        sequence=shot.strands[0].sequence
        shot.strands[0].sequence=sequence[:event.primer_length]+'A'+sequence[event.primer_length+1:]
        with self.assertRaisesRegex(ValueError,'complement'):m.validate()

    def test_fill_must_precede_closure(self):
        from .protocols import crispr_mip,closure_event
        m=crispr_mip();m.shots[1].events.append(closure_event('probe','Regression end fixture',start=0))
        with self.assertRaisesRegex(ValueError,'before extension'):m.validate()

    def test_bending_does_not_create_a_covalent_circle(self):
        from .evaluate import evaluate
        from .protocols import crispr_mip
        s=crispr_mip().shots[0]
        states,_=evaluate(s,s.establish+s.action)
        self.assertFalse(states['probe'].closed)

    def test_sealing_preserves_polymer_and_closes_exactly_one_nick(self):
        from .evaluate import evaluate
        from .protocols import crispr_mip
        s=crispr_mip().shots[-1];before,_=evaluate(s,0);after,_=evaluate(s,s.establish+s.action)
        self.assertEqual(before['probe'].strand.sequence,after['probe'].strand.sequence)
        self.assertFalse(before['probe'].closed);self.assertTrue(after['probe'].closed)
        self.assertEqual(after['probe'].nicks[0],0)

    def test_circle_requires_ordered_closed_route(self):
        from .protocols import boloramis
        m=boloramis();shot=next(s for s in m.shots if s.id=='rca');shot.strands[0].path[-1]=(len(shot.strands[0].sequence),1,1)
        with self.assertRaisesRegex(ValueError,'geometrically closed'):m.validate()

    def test_rca_cannot_copy_linear_template_as_a_circle(self):
        from .protocols import boloramis
        m=boloramis();next(s for s in m.shots if s.id=='rca').strands[0].circular=False
        with self.assertRaisesRegex(ValueError,'covalently closed'):m.validate()

    def test_rca_copied_repeat_and_origin_are_checked(self):
        from .protocols import boloramis
        m=boloramis();s=next(shot for shot in m.shots if shot.id=='rca').strands[2]
        s.sequence=s.sequence[:-1]+('A'if s.sequence[-1]!='A'else'C')
        with self.assertRaisesRegex(ValueError,'complement'):m.validate()

    def test_rca_threeprime_front_stays_near_template_polymerase(self):
        from .protocols import boloramis
        from .evaluate import evaluate
        import math
        shot=next(s for s in boloramis().shots if s.id=='rca')
        for f in [0,.2,.6,1]:
            states,actors=evaluate(shot,shot.establish+shot.action*f)
            child=states['dna'];front=child.point(child.visible[0][1])
            self.assertLess(math.dist(front,actors[0].position),20)

    def test_retained_identity_cannot_silently_change_bases(self):
        from .protocols import smallseq
        m=smallseq();shot=m.shots[1]
        shot.inventory['ra3_input']='A'*len(shot.inventory['ra3_input'])
        with self.assertRaisesRegex(ValueError,'original molecule'):m.validate()

    def test_pairing_rejects_parallel_or_noncomplementary_arms(self):
        from .protocols import crispr_mip
        m=crispr_mip();m.shots[0].pairings[0].b_step=1
        with self.assertRaisesRegex(ValueError,'antiparallel'):m.validate()
        m=crispr_mip();m.shots[0].pairings[0].b_start+=1
        with self.assertRaisesRegex(ValueError,'complementary'):m.validate()

    def test_conversion_retains_protected_c_and_physical_dna_u(self):
        from .protocols import scrrbs
        from .evaluate import evaluate
        m=scrrbs();shot=m.shots[1];states,actors=evaluate(shot,1000)
        self.assertEqual(set(states['original'].converted.values()),{'U'})
        self.assertNotIn(3,states['original'].converted)
        self.assertEqual(states['original'].strand.polymer,'DNA')
        self.assertEqual(actors,[])  # Bisulfite is not a protein enzyme.

    def test_source_missing_slots_never_become_exact_claims(self):
        from .examples import EXAMPLES
        m=EXAMPLES['smartseq']()
        self.assertEqual(m.provenance['original_protocol'],'Ramskold2012; not SMART-seq2')
        self.assertIn('index',m.provenance['steps_used'])
        self.assertFalse(any(e.exact for s in m.shots for e in s.events))

    def test_ambiguous_base_family_does_not_establish_exact_copy(self):
        m=reference();m.shots[0].strands[0].sequence='N'*len(m.shots[0].strands[0].sequence)
        with self.assertRaisesRegex(ValueError,'unknown classes'):m.validate()

    def test_ordinary_model_use_does_not_require_animation_libraries(self):
        code="import sys;from animation.examples import EXAMPLES;EXAMPLES['smartseq']().validate();assert 'cairo'not in sys.modules"
        subprocess.run([sys.executable,'-S','-c',code],check=True)

    def test_all_protocols_have_explicit_source_bindings_and_roundtrip(self):
        from .examples import EXAMPLES
        for name,fn in EXAMPLES.items():
            with self.subTest(example=name):
                m=fn();r=Movie.from_dict(json.loads(json.dumps(m.to_dict())))
                self.assertEqual(json.dumps(m.to_dict(),sort_keys=True),json.dumps(r.to_dict(),sort_keys=True))
                if m.provenance.get('representation')!='illustrative_example':
                    self.assertEqual(set(s.id for s in m.shots),set(m.provenance['steps_used']))
                    self.assertTrue(m.provenance['module_sha256'])


class PresentationTests(unittest.TestCase):
    def test_default_speed_and_seek_share_multiplier(self):
        m=reference()
        self.assertAlmostEqual(m.duration()*1.5,m.duration(speed=1))
        for t in [0,1,3,5]:self.assertEqual(m.at(t),m.at(t*1.5,speed=1))

    def test_poly_a_and_primer_handle_are_physically_separate(self):
        m=smartseq();shot=m.shots[0]
        rna,primer=shot.strands;pair=shot.pairings[0]
        # Check the displayed representative member, independently of schematic body.
        from .model import complement
        for j in range(pair.length):
            self.assertEqual(primer.sequence[pair.a_start+j],complement(rna.sequence[pair.b_start-j]))
        self.assertTrue(primer.regions[0].end<=pair.a_start)
        self.assertEqual(primer.sequence[pair.a_start:pair.a_start+30],'T'*30)

    def test_rna_and_cdna_persist_through_template_switch(self):
        shots={s.id:s for s in smartseq().shots}
        for id in ['rt','terminal','switch']:
            rna=next(s for s in shots[id].strands if s.id=='rna')
            self.assertEqual(rna.sequence,shots['anchor'].strands[0].sequence)
            self.assertEqual(rna.start,shots['anchor'].strands[0].start)
        first=next(s for s in shots['firstcopy'].strands if s.id=='reverse')
        self.assertLess(first.initially_visible,len(first.sequence))
        self.assertTrue(any(e.kind=='extend'and e.strand=='reverse'for e in shots['firstcopy'].events))

    def test_transposition_has_dimers_and_retains_cleaved_pieces(self):
        from .protocols import tagmentation_shot
        from .evaluate import evaluate
        shot=tagmentation_shot();shot.source='Schematic test';shot.validate()
        state,actors=evaluate(shot,shot.establish+shot.action)
        self.assertEqual(sum(a.complex=='homodimer'for a in actors),2)
        top=state['top']
        self.assertEqual(top.visible,[(0.,140.)])
        self.assertLess(top.segment_point(30,29)[0],top.segment_point(30,30)[0])
        self.assertGreater(top.segment_point(105,106)[0],top.segment_point(105,104)[0])

    def test_index_overhang_is_not_forced_to_pair(self):
        shot=next(s for s in smartseq().shots if s.id=='index')
        event=next(e for e in shot.events if e.kind=='extend')
        self.assertLess(event.primer_binding_length,event.primer_length)
        self.assertTrue(any(r.label=='Sample index'for r in shot.strands[1].regions))
        self.assertEqual(len([e for e in shot.events if e.kind=='extend']),2)

    def test_functional_regions_and_initial_pieces_roundtrip(self):
        from .protocols import smallseq
        from .evaluate import evaluate
        m=smallseq();r=Movie.from_dict(json.loads(json.dumps(m.to_dict())))
        shot=next(s for s in r.shots if s.id=='lig3')
        start,_=evaluate(shot,0);end,_=evaluate(shot,shot.timing()[2])
        junction=shot.strands[0].nicks[0]
        self.assertGreater(start['rna_adapter'].segment_point(junction,junction)[0],start['rna_adapter'].segment_point(junction,junction-1)[0])
        self.assertEqual(end['rna_adapter'].segment_point(junction,junction),end['rna_adapter'].segment_point(junction,junction-1))
        self.assertEqual(end['rna_adapter'].nicks[junction],0)


class SequenceRevisionTests(unittest.TestCase):
    def test_pairing_state_is_independent_of_distance(self):
        from .evaluate import evaluate_pairings
        from .model import Pairing
        shot=reference().shots[0]
        shot.pairings=[Pairing('template',0,'daughter',len(shot.strands[1].sequence)-1,5,id='duplex')]
        before=evaluate_pairings(shot,0)[0].strength
        shot.strands[1].start=(10000,10000);shot.strands[1].end=(12000,10000)
        self.assertEqual(evaluate_pairings(shot,0)[0].strength,before)
        shot.events.append(Event('unpair','daughter',pairing='duplex'))
        self.assertEqual(evaluate_pairings(shot,shot.timing()[2])[0].strength,0)

    def test_incoming_primers_start_outside_viewport(self):
        from .smartseq2 import movie
        from .evaluate import evaluate
        for shot_id,strand_ids in [('anchor',['cdna']),('switch',['tso']),('firstcopy',['reverse']),('index',['index_plus','index_minus'])]:
            shot=next(s for s in movie().shots if s.id==shot_id)
            states,_=evaluate(shot,0)
            for id in strand_ids:
                st=states[id];end=st.visible[0][1]
                xs=[st.point(i)[0]for i in [0,end]]
                self.assertTrue(max(xs)<40 or min(xs)>1560,(shot_id,id,xs))

    def test_paired_bases_keep_tight_spacing_through_transitions(self):
        import math
        from .smartseq2 import movie
        from .evaluate import evaluate,evaluate_pairings
        for shot in movie().shots:
            for phase in [0,.25,.5,.75,1]:
                t=shot.establish+shot.action*phase;states,_=evaluate(shot,t)
                for bound in evaluate_pairings(shot,t):
                    if bound.strength<.999:continue
                    p=bound.pairing;a,b=states[p.a],states[p.b]
                    for j in range(p.length):
                        ia=p.a_start+j;ib=p.b_start-j
                        if a.fraction(ia)<.99 or b.fraction(ib)<.99:continue
                        distance=math.dist(a.segment_point(ia+.5,ia),b.segment_point(ib+.5,ib))
                        self.assertAlmostEqual(distance,24,places=4,msg=f'{shot.id}:{phase}:{p.key}:{j}')

    def test_only_nine_base_bridges_unpair_when_tn5_flanks_leave(self):
        from .smartseq2 import movie
        from .molecular_frame import molecular_frame
        shot=next(s for s in movie().shots if s.id=='tagment')
        initial=molecular_frame(shot,0);final=molecular_frame(shot,shot.timing()[2])
        lost=initial.duplex-final.duplex
        self.assertEqual(len(lost),18)
        self.assertFalse(final.duplex-initial.duplex)
        for event in [e for e in shot.events if e.kind=='unpair']:
            removed={tuple(sorted(p))for p in event.molecular_operation['pairing_removed']}
            self.assertEqual(len(removed),9)
            self.assertTrue(removed<=lost)

    def test_transfer_junctions_and_target_lineage(self):
        from .smartseq2 import movie
        from .evaluate import evaluate
        m=movie();tag=next(s for s in m.shots if s.id=='tagment');fill=next(s for s in m.shots if s.id=='fill')
        joins=[e for e in tag.events if e.kind=='strand_transfer']
        self.assertEqual(len(joins),4)
        self.assertTrue(all(e.note and e.other in {'target_top','target_bottom'}for e in joins))
        start,_=evaluate(fill,0);end,_=evaluate(fill,fill.timing()[2])
        for id in ['selected_top','selected_bottom']:
            self.assertLess(start[id].visible[0][1],len(start[id].strand.sequence))
            self.assertEqual(end[id].visible[0][1],len(end[id].strand.sequence))
            self.assertTrue(any(origin.startswith('target_')for a,b,origin,j in start[id].strand.lineage))
        self.assertEqual(start['selected_top'].strand.sequence,complement_for_test(start['selected_bottom'].strand.sequence[::-1]))
        self.assertFalse(any(e.kind=='seal'for e in fill.events))

    def test_flanking_dna_exits_while_central_dna_stays(self):
        from .smartseq2 import movie
        from .evaluate import evaluate
        shot=next(s for s in movie().shots if s.id=='tagment')
        states,_=evaluate(shot,shot.timing()[2])
        top=states['target_top'];moves=[e for e in shot.events if e.kind=='displace'and e.strand=='target_top']
        for e in moves:
            xs=[top.segment_point(i,(e.span[0]+e.span[1])/2)[0]for i in e.span]
            self.assertTrue(max(xs)<40 or min(xs)>1560,xs)
        lo=moves[0].span[1];hi=moves[1].span[0]
        self.assertTrue(40<top.segment_point((lo+hi)/2,(lo+hi)/2)[0]<1560)

    def test_index_first_extension_stops_at_actual_template_end(self):
        from .smartseq2 import movie
        m=movie();first=next(s for s in m.shots if s.id=='index');second=next(s for s in m.shots if s.id=='index_copy')
        by_id={s.id:s for s in first.strands};next_by_id={s.id:s for s in second.strands}
        for e in [e for e in first.events if e.kind=='extend']:
            self.assertEqual(e.copy_end,len(by_id[e.strand].sequence))
            self.assertEqual(e.molecular_operation['details']['stop'],'physical_template_end')
            self.assertEqual(e.copy_end-e.primer_length,e.template_start+1)
            original=by_id[e.strand]
            parent=next(s for s in second.strands if s.nucleotide_ids[0]==original.nucleotide_ids[0])
            self.assertEqual(parent.sequence[:e.copy_end],original.sequence)
            self.assertEqual(parent.nucleotide_ids[:e.copy_end],original.nucleotide_ids)
        for e in [e for e in second.events if e.kind=='extend']:
            out=next_by_id[e.strand];template=next_by_id[e.template]
            copied=complement_for_test(template.sequence[:e.template_start+1][::-1])
            self.assertEqual(out.sequence[e.primer_length:],copied)

    def test_final_view_shows_every_defined_component_and_nucleotide(self):
        from .smartseq2 import movie
        m=movie();final=m.shots[-1];components=m.provenance['final_components']
        self.assertEqual(len(components),9)
        upper=[s for s in final.strands if s.id.startswith('library_')and not s.id.startswith('library_bottom')]
        self.assertEqual(''.join(s.sequence for s in upper),''.join(c['sequence']for c in components))
        self.assertEqual(sum(len(s.regions)for s in upper),9)
        self.assertTrue(upper[0].continued_3prime and upper[1].continued_5prime)
        self.assertFalse(upper[0].continued_5prime or upper[1].continued_3prime)


def complement_for_test(seq):
    from .model import complement
    return complement(seq)


class NativeTests(unittest.TestCase):
    def test_unborn_transcripts_have_no_captions_or_end_markers(self):
        from .synthesis_examples import celseq
        from .render import Renderer,Style
        class Labels(Renderer):
            def text(self,c,text,x,y,size=24,color='text',align='left',max_width=None):
                if self.in_molecule_layer and (text=='RNA transcript' or text in {'5′','3′'} and color=='rna'):self.ghosts.append(text)
                return super().text(c,text,x,y,size,color,align,max_width)
        m=celseq();shot=m.shots[1];offset=m.shots[0].timing()[2];r=Labels(Style(960,540));r.ghosts=[]
        r.image(m,(offset+shot.transition/2)/1.5)
        self.assertEqual(r.ghosts,[])
        r.ghosts=[];r.image(m,(offset+shot.establish+shot.action*.5)/1.5)
        self.assertIn('RNA transcript',r.ghosts)

    def test_padlock_end_labels_do_not_overlap_molecule_captions(self):
        from .boloramis import movie
        from .render import Renderer,Style
        from dataclasses import replace
        class Ink(Renderer):
            def text(self,c,text,x,y,size=24,color='text',align='left',max_width=None):
                if self.in_molecule_layer and text:
                    c.save();c.select_font_face(self.style.font);c.set_font_size(size);e=c.text_extents(text);c.restore()
                    left=x-e[2]/2 if align=='center'else x-e[2]if align=='right'else x
                    box=(left,y+e[1],left+e[2],y+e[1]+e[3]);clip=c.clip_extents()
                    if box[0]>=clip[0]and box[1]>=clip[1]and box[2]<=clip[2]and box[3]<=clip[3]:
                        if text in {'3′','5′'}:
                            for b in self.ink:
                                if box[0]<b[2]and box[2]>b[0]and box[1]<b[3]and box[3]>b[1]:self.overlaps.append((text,box,b))
                        self.ink.append(box)
                return super().text(c,text,x,y,size,color,align,max_width)
        m=movie();r=Ink(Style(960,540))
        for shot in m.shots:
            for phase in [.4,.9]:
                r.ink=[];r.overlaps=[];r.image(replace(m,shots=[shot]),(shot.establish+shot.action*phase)/1.5)
                self.assertFalse(r.overlaps,(shot.id,phase,r.overlaps))
    def test_timeline_is_opt_in(self):
        from .render import Renderer,Style
        self.assertFalse(Style().timeline)
        plain=Renderer(Style(320,180)).image(reference(),0)
        detailed=Renderer(Style(320,180,timeline=True)).image(reference(),0)
        # The footer is clear by default; the counter/progress ink is optional.
        start=plain.get_stride()*164
        self.assertNotEqual(bytes(plain.get_data())[start:],bytes(detailed.get_data())[start:])
        pixels=bytes(plain.get_data())[start:]
        self.assertTrue(all(pixels[i:i+4]==pixels[:4]for i in range(0,len(pixels),4)))
    def test_index_polymerase_labels_remain_inside_viewport(self):
        from .smartseq2 import movie
        from .render import Renderer,Style
        from dataclasses import replace
        class Labels(Renderer):
            def text(self,c,text,x,y,size=24,color='text',align='left',max_width=None):
                if text=='DNA polymerase':
                    c.save();c.select_font_face(self.style.font);c.set_font_size(size)
                    box=c.text_extents(text);c.restore();clip=c.clip_extents()
                    left=x-box[2]/2 if align=='center' else x
                    self.polymerase_labels.append((left>=clip[0],left+box[2]<=clip[2],y+box[1]>=clip[1],y+box[1]+box[3]<=clip[3]))
                return super().text(c,text,x,y,size,color,align,max_width)
        m=movie();r=Labels(Style(960,540))
        for shot in m.shots:
            if shot.id not in {'index','index_copy','index_duplex'}:continue
            for phase in [.75,1]:
                r.polymerase_labels=[]
                r.image(replace(m,shots=[shot]),(shot.establish+shot.action*phase)/1.5)
                self.assertEqual(len(r.polymerase_labels),len([e for e in shot.events if e.kind=='extend']),(shot.id,phase))
                self.assertTrue(all(all(box)for box in r.polymerase_labels),(shot.id,phase,r.polymerase_labels))

    def test_deterministic_pixels_and_stills(self):
        from .render import Renderer,Style
        r=Renderer(Style(640,360));m=reference()
        a=bytes(r.image(m,4).get_data());r.image(m,9);b=bytes(r.image(m,4).get_data())
        self.assertEqual(a,b)
        with tempfile.TemporaryDirectory()as tmp:
            png=Path(tmp)/'frame.png';svg=Path(tmp)/'frame.svg'
            r.png(png,m,4);r.svg(svg,m,4)
            self.assertTrue(png.read_bytes().startswith(b'\x89PNG'))
            self.assertIn('<svg',svg.read_text())

    def test_native_video_and_failure(self):
        from .render import Renderer,Style
        from .export import export_mp4
        r=Renderer(Style(320,180));m=reference()
        # Shorten presentation for smoke test without altering the molecular event.
        s=m.shots[0];s.narration='Copy.';s.establish=.1;s.action=.3;s.hold=.1;s.narration_seconds=0
        with tempfile.TemporaryDirectory()as tmp:
            path=Path(tmp)/'clip.mp4';meta=export_mp4(m,path,r,fps=5)
            probe=json.loads(subprocess.check_output(['ffprobe','-v','error','-show_streams','-of','json',str(path)]))['streams'][0]
            self.assertEqual((probe['codec_name'],probe['pix_fmt']),('h264','yuv420p'))
            self.assertEqual(int(probe['nb_frames']),meta['frames'])
            self.assertTrue(path.with_suffix('.manifest.json').exists())
            with self.assertRaisesRegex(RuntimeError,'FFmpeg export failed'):
                export_mp4(m,Path(tmp)/'bad.mp4',r,fps=5,ffmpeg='false')
            self.assertFalse((Path(tmp)/'bad.mp4').exists())


if __name__=='__main__':
    from .selftest_padlock import PadlockTests
    suite=unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromTestCase(MolecularTests))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(ReactionTests))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(PresentationTests))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(SequenceRevisionTests))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(RNaseHTests))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(SubstrateSimulationTests))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(PadlockTests))
    from .selftest_capture_ligation import CaptureLigationTests
    from .selftest_synthesis_examples import SynthesisTests
    from .selftest_assembly_transposition import AssemblyTranspositionTests
    from .selftest_golden_gate import GoldenGateTests
    from .selftest_collection import CollectionTests
    for cls in [CaptureLigationTests,SynthesisTests,AssemblyTranspositionTests,GoldenGateTests,CollectionTests]:
        suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(cls))
    if '--native'in sys.argv:suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(NativeTests))
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(not result.wasSuccessful())
