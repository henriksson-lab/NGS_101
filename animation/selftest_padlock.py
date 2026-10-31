"""Causal chemistry and layout checks for the generated BOLORAMIS example."""
import collections
import math
import unittest
from .boloramis import simulate,movie,SLUG
from .protocols import load
from .molecular_frame import molecular_frame
from .evaluate import evaluate,evaluate_pairings


class PadlockTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.owner=load(SLUG,'chem.py');cls.reaction=simulate(cls.owner);cls.movie=movie()

    def test_reframing_starts_at_previous_coordinates_without_changing_chemistry(self):
        from .evaluate import presentation
        for previous,shot in zip(self.movie.shots,self.movie.shots[1:]):
            old,_=evaluate(previous,previous.timing()[2])
            incoming,_=evaluate(shot,0)
            actual,_,camera=presentation(self.movie,shot,0)
            self.assertEqual(camera,previous.camera)
            locations={id:(st,j)for st in old.values()for j,id in enumerate(st.strand.nucleotide_ids)if st.fraction(j)}
            for key,st in actual.items():
                self.assertEqual(st.molecular,incoming[key].molecular)
                for j,id in enumerate(st.strand.nucleotide_ids):
                    if st.fraction(j):
                        before,at=locations[id]
                        self.assertLess(math.dist(st.point(j+.5),before.segment_point(at+.5,at)),1e-8)
            end,_,camera=presentation(self.movie,shot,shot.transition)
            self.assertEqual(camera,shot.camera)
            expected,_=evaluate(shot,shot.transition)
            for key,st in end.items():
                self.assertEqual(st.molecular,expected[key].molecular)
                for j in range(len(st.strand.sequence)):
                    self.assertEqual(st.point(j+.5),expected[key].point(j+.5))

    def test_circle_keeps_its_size_while_reorienting(self):
        from .evaluate import presentation
        shot=self.movie.shots[-1]
        reference=evaluate(shot,0)[0]['probe']
        mid=presentation(self.movie,shot,shot.transition/2)[0]['probe']
        for a,b in [(0,18),(18,35),(35,53),(0,35)]:
            self.assertAlmostEqual(math.dist(mid.point(a),mid.point(b)),math.dist(reference.point(a),reference.point(b)),places=6)

    def test_reframing_cannot_hide_a_chemical_change(self):
        from dataclasses import replace
        from .evaluate import presentation
        m=replace(self.movie,shots=[self.movie.shots[0],self.movie.shots[2]])
        with self.assertRaisesRegex(ValueError,'unchanged molecular state'):
            presentation(m,m.shots[1],0)

    def test_ligase_departs_and_sealed_junction_remains_marked(self):
        close=next(s for s in self.movie.shots if s.id=='close')
        states,actors=evaluate(close,close.establish+close.action)
        self.assertTrue(states['probe'].closed)
        self.assertEqual(states['probe'].strand.sealed_marks,[0])
        self.assertEqual(actors[0].label,'SplintR ligase')
        self.assertEqual(actors[0].progress,1)
        for shot in self.movie.shots[3:]:
            self.assertEqual(next(s for s in shot.strands if s.id=='probe').sealed_marks,[0])
        product=self.movie.shots[-1]
        self.assertEqual(evaluate(product,0)[1][0].label,'Phi29 polymerase')

    def test_ligation_preserves_atoms_and_adds_only_the_actual_closure_bond(self):
        r=self.reaction;op=r['seal'];p=r['probe']
        self.assertEqual(op.before.atoms,op.after.atoms)
        self.assertEqual(op.before.paired,op.after.paired)
        self.assertEqual(set(op.after.next.items())-set(op.before.next.items()),{(p[-1],p[0])})
        self.assertEqual(op.before.five[p[0]],'P');self.assertEqual(op.before.three[p[-1]],'OH')
        self.assertNotIn(p[0],op.after.five);self.assertNotIn(p[-1],op.after.three)
        self.assertTrue(all(id in op.before.atoms for id in r['primer']))
        self.assertTrue(all(op.after.paired[id]==op.before.paired[id]for id in r['primer']))

    def test_rca_nucleotides_are_supplied_after_ligation(self):
        r=self.reaction
        self.assertFalse(r['initial'].ntp);self.assertFalse(r['seal'].after.ntp)
        self.assertEqual(r['closed'].ntp['DNA'],set('ACGT'))

    def test_missing_phosphate_and_unsupported_pairing_are_rejected(self):
        with self.assertRaisesRegex(ValueError,'3OH and5P'):simulate(self.owner,phosphate='OH')
        seq=self.owner.ConditionalRNAMatch().sequence
        with self.assertRaises(ValueError):simulate(self.owner,rna_sequence=('A'if seq[0]!='A'else'C')+seq[1:])

    def test_rna_displacement_capability_changes_the_product_instead_of_erasing_rna(self):
        r=simulate(self.owner,rna_displacement=False)
        self.assertEqual(r['rca'].details['stop'],'occupied_template')
        self.assertLess(len(r['rca'].added),len(self.reaction['rca'].added))
        self.assertTrue(all(id in r['after'].paired for id in r['rna']))
        for reaction in [r,self.reaction]:
            self.assertFalse(reaction['after'].released)
            for a,b in zip(reaction['rna'],reaction['rna'][1:]):self.assertEqual(reaction['after'].next[a],b)
            self.assertTrue(all(reaction['initial'].atoms[id]==reaction['after'].atoms[id]for id in reaction['rna']))
        self.assertTrue(all(id not in self.reaction['after'].paired for id in self.reaction['rna']))

    def test_rca_walks_the_actual_circle_and_keeps_the_stock_primer(self):
        r=self.reaction;op=r['rca'];counts=collections.Counter(op.after.atoms[id].template for id in op.added)
        self.assertEqual(counts,{id:3 for id in r['probe']})
        self.assertEqual(tuple(r['product'][:len(r['primer'])]),tuple(r['primer']))
        self.assertEqual(len(r['product']),len(r['primer'])+len(op.added))
        for id in op.added:
            atom=op.after.atoms[id]
            if atom.base=='T':self.assertIn('occupancy unassigned',atom.modification)

    def test_cut_circle_or_unpaired_primer_blocks_the_actual_rca_operation(self):
        from .simulator import Polymerase
        for perturbation in ['cut','unpair']:
            r=simulate(self.owner,budget=17);s=r['simulator']
            if perturbation=='cut':s.cleave(r['probe'][-1],source='synthetic circle disruption')
            else:s.unpair(s.chain('primer'),source='synthetic primer dissociation')
            before=s.snapshot();enzyme=Polymerase(**r['rca'].details['enzyme'])
            with self.assertRaisesRegex(ValueError,'covalently closed'if perturbation=='cut'else'annealed primer'):
                s.rolling_circle('primer',enzyme,max_bases=60,source='causal test')
            self.assertEqual(s.snapshot(),before)

    def test_input_perturbations_propagate_through_the_same_reaction(self):
        seq=self.owner.probe().top();changed=seq[:47]+('A'if seq[47]!='A'else'C')+seq[48:]
        alternate=simulate(self.owner,probe_sequence=changed)
        def product(r):return ''.join(r['after'].atoms[id].base for id in r['product'])
        self.assertNotEqual(product(alternate),product(self.reaction))
        shorter=simulate(self.owner,primer_sequence=self.owner.RCA_PRIMER[1:])
        self.assertEqual(len(shorter['product']),len(self.reaction['product'])-1)
        small=simulate(self.owner,budget=17)
        self.assertEqual(len(small['product']),len(small['primer'])+17)

    def test_folding_is_not_ligation_and_bound_primers_follow_the_backbone(self):
        shot=next(s for s in self.movie.shots if s.id=='fold');initial=molecular_frame(shot,0)
        for phase in [0,.2,.5,.74,1]:
            t=shot.establish+shot.action*phase;frame=molecular_frame(shot,t);states,_=evaluate(shot,t)
            self.assertEqual(frame.backbone,initial.backbone)
            self.assertFalse(states['probe'].closed)
            for bound in evaluate_pairings(shot,t,states=states):
                p=bound.pairing;a,b=states[p.a],states[p.b]
                distance=math.dist(a.point(p.a_start+.5),b.point(p.b_start+.5))
                self.assertAlmostEqual(distance,24,places=5)

    def test_unbound_short_arm_does_not_pass_through_the_long_arm(self):
        shot=next(s for s in self.movie.shots if s.id=='fold')
        n=len(self.reaction['probe']);long=len(self.owner.LONG);short=len(self.owner.SHORT)
        for step in range(21):
            state=evaluate(shot,shot.establish+shot.action*step/20)[0]['probe']
            for i in range(long):
                a=state.point(i+.5)
                for j in range(n-short,n):
                    b=state.point(j+.5)
                    if abs(a[0]-b[0])<8:self.assertGreaterEqual(abs(a[1]-b[1]),24)

    def test_rca_duplex_and_displaced_prefix_come_from_actual_bonds(self):
        shot=next(s for s in self.movie.shots if s.id=='rca')
        for phase in [0,.1,.3,.5,.75,1]:
            t=shot.establish+shot.action*phase;states,_=evaluate(shot,t)
            visible=int(states['dna'].visible[0][1]);dna_pairs=0
            for bound in evaluate_pairings(shot,t,states=states):
                p=bound.pairing
                if 'dna' in {p.a,p.b}:dna_pairs+=1
                a,b=states[p.a],states[p.b]
                self.assertAlmostEqual(math.dist(a.point(p.a_start+.5),b.point(p.b_start+.5)),24,places=5)
            self.assertEqual(dna_pairs,min(visible,len(self.reaction['probe'])))
            self.assertEqual(len(states['rna'].visible),1)
            self.assertEqual(states['rna'].visible[0][1],len(self.reaction['rna']))

    def test_final_layout_changes_no_molecules_or_bonds(self):
        copy=next(s for s in self.movie.shots if s.id=='rca');final=self.movie.shots[-1]
        a=molecular_frame(copy,copy.establish+copy.action);b=molecular_frame(final,0)
        self.assertEqual(a.backbone,b.backbone);self.assertEqual(a.duplex,b.duplex)
        self.assertEqual(a.atoms,b.atoms);self.assertEqual(a.released,b.released)
        self.assertEqual(next(s for s in final.strands if s.id=='dna').nucleotide_ids,list(self.reaction['product']))


if __name__=='__main__':unittest.main(verbosity=2)
