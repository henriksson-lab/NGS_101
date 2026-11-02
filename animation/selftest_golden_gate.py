"""Causal checks for Type IIS digestion and four-bond cohesive assembly."""
import unittest,math
from dataclasses import replace
from .golden_gate import simulate,movie
from .evaluate import evaluate,presentation
from .molecular_frame import molecular_frame

class GoldenGateTests(unittest.TestCase):
    def test_cuts_are_offset_and_both_strands_conserved(self):
        r=simulate();self.assertEqual(len(r['cutops']),4)
        for op in r['cutops']:
            self.assertEqual((op.details['top_cut'],op.details['bottom_cut']),(7,11))
            self.assertEqual(len(op.to_dict()['backbone_removed']),2)
            self.assertEqual(op.before.atoms,op.after.atoms)
        self.assertEqual(r['initial'].atoms,r['final'].atoms)
        self.assertFalse(r['final'].released)
    def test_actual_four_phosphodiester_bonds_close_both_strands(self):
        r=simulate();s=r['simulator']
        for ids in [r['top'],r['bottom']]:
            self.assertEqual(s.chain(ids[0]),ids)
            self.assertEqual(r['final'].next[ids[-1]],ids[0])
            self.assertTrue(all(r['final'].next[a]==b for a,b in zip(ids,ids[1:])))
        self.assertEqual(sum(len(op.to_dict()['backbone_added'])for op in r['seals']),4)
        seq=s.sequence(r['top'][0]);site=r['module'].BsaI.site
        self.assertNotIn(site,seq+seq);self.assertNotIn(r['module'].revcomp(site),seq+seq)
    def test_missing_phosphate_cofactor_or_mismatched_fusion_blocks_assembly(self):
        for kwargs in [dict(phosphate=False),dict(cofactors=()),dict(fusion_mutation=True)]:
            with self.assertRaises(ValueError):simulate(**kwargs)
    def test_every_atom_is_drawn_once_in_each_view(self):
        for shot in movie().shots:
            ids=[id for s in shot.strands for id in s.nucleotide_ids]
            self.assertEqual(len(ids),len(set(ids)))
            self.assertEqual(set(ids),set(shot.initial_nucleotides))
    def test_transitions_keep_molecular_state_and_start_positions(self):
        m=movie()
        for prev,shot in zip(m.shots,m.shots[1:]):
            old=evaluate(prev,prev.timing()[2])[0]
            st,_,camera=presentation(m,shot,0)
            self.assertEqual(camera,prev.camera)
            self.assertEqual(next(iter(old.values())).molecular,next(iter(st.values())).molecular)
            locations={id:(s,j)for s in old.values()for j,id in enumerate(s.strand.nucleotide_ids)}
            for s in st.values():
                for j,id in enumerate(s.strand.nucleotide_ids):
                    before,k=locations[id]
                    self.assertLess(math.dist(s.point(j+.5),before.segment_point(k+.5,k)),1e-7)
    def test_intact_reagent_enters_without_fake_synthesis(self):
        from .simulator import Simulator
        from .generated_views import strand
        from .simulation_adapters import generated
        from .model import Shot,Event
        s=Simulator();s.add('old','ACG',source='test');before=s.snapshot();ids=s.add('stock','TAGCTAGC',source='test');op=s.history[-1]
        shot=Shot('s','Stock reagent',[strand('stock',s.snapshot(),ids,-400,400,initial=0)],
            [Event('introduce','stock',fraction_start=.2,molecular_operation=op.to_dict())],'An intact stock oligo enters.',source='test')
        generated(shot,before);shot.validate()
        self.assertNotIn(ids[0],molecular_frame(shot,0).atoms)
        frame=molecular_frame(shot,shot.establish+shot.action*.2)
        self.assertTrue(set(ids)<=set(frame.atoms))
        self.assertTrue(all((a,b)in frame.backbone for a,b in zip(ids,ids[1:])))

if __name__=='__main__':unittest.main()
