"""Consistency checks for the curved-path teaching model; not validation."""
import unittest
import numpy as np
from xfem import Model, quadrature, enrichment, resistance, strike, trial_endpoint, WIDTH, HEIGHT, length


def crack(points):
    return {'points': points, 'stopped': False}


def inputs(mesh='coarse', **overrides):
    return dict(mesh=mesh, young=70, poisson=.22, toughness=10, mass=.5, height=.12, aim=.5, variation=.45, cracks=[], **overrides)


class NumericalChecks(unittest.TestCase):
    def test_kinked_subcell_area_is_conserved(self):
        model = Model('coarse')
        cracks = [crack([[.3,.073],[.27,.09],[.24,.084],[.18,.12]]), crack([[.3,.15],[.25,.14],[.2,.17]])]
        area = sum(quadrature(v, cracks, model.h)[1].sum() for v, _, _ in model.geometry)
        self.assertAlmostEqual(area, WIDTH*HEIGHT, places=12)

    def test_rotated_tip_gradient_matches_finite_difference(self):
        c = crack([[.3,.1],[.26,.08],[.22,.1]])
        p = np.array([[.20,.112],[.244,.083],[.281,.106]])
        _, gradient = enrichment(p, c, 'tip', .025)
        for axis in range(2):
            delta = np.zeros(2); delta[axis] = 1e-7
            numerical = (enrichment(p+delta, c, 'tip', .025)[0]-enrichment(p-delta, c, 'tip', .025)[0])/(2e-7)
            np.testing.assert_allclose(gradient[:,axis], numerical, rtol=2e-6, atol=1e-6)

    def test_jump_follows_polyline_not_extrapolated_tip_line(self):
        c = crack([[.3,.1],[.26,.08],[.22,.1]])
        # Actual crack at x=.28 is y=.09; tip-tangent extension is y=.07.
        face_values = enrichment([[.28,.090001],[.28,.089999]], c, 'tip', .025)[0]
        self.assertLess(face_values.prod(), 0)
        ghost_values = enrichment([[.28,.070001],[.28,.069999]], c, 'tip', .025)[0]
        self.assertGreater(ghost_values.prod(), 0)
        self.assertLess(abs(ghost_values[0]-ghost_values[1]), 1e-3)

    def test_linear_scaling_and_equilibrium_with_curved_crack(self):
        model = Model('coarse')
        cracks = [crack([[.3,.103],[.26,.09],[.22,.115],[.17,.12]])]
        a, b = model.solve(cracks,20e9,.103), model.solve(cracks,40e9,.103)
        self.assertGreater(a['enriched_dofs'], 0)
        self.assertLess(a['residual'], 1e-7)
        self.assertAlmostEqual(a['energy'], 2*b['energy'], delta=a['energy']*1e-7)
        np.testing.assert_allclose(a['values'],b['values'],rtol=1e-6)

    def test_intact_reflection_symmetry(self):
        model = Model('medium')
        a,b = model.solve([],20e9,.06),model.solve([],20e9,.14)
        self.assertLess(abs(a['energy']/b['energy']-1),.02)

    def test_uniform_and_fixed_varied_toughness(self):
        q = np.array([[.2,.1],[.1,.1],[.15,.08]])
        np.testing.assert_array_equal(resistance(q,0),np.ones(3))
        np.testing.assert_array_equal(resistance(q,.45),resistance(q,.45))
        self.assertGreater(np.ptp(resistance(q,.45)),.1)
        self.assertTrue((resistance(q,.45)>0).all())

    def test_growth_turns_and_multiple_cracks_on_all_meshes(self):
        for mesh in ['coarse','medium','fine']:
            with self.subTest(mesh=mesh):
                p=inputs(mesh); p.update(mass=.25,height=.01)
                self.assertFalse(strike(p)['accepted'])
                p.update(mass=.5,height=.12)
                r=strike(p)
                self.assertTrue(r['accepted']); self.assertGreater(len(r['events']),1)
                first=np.asarray(r['cracks'][0]['points'])
                self.assertGreater(np.ptp(first[:,1]),.001)
                self.assertTrue((np.diff(first[:,0])<0).all())
                p['cracks'],p['aim']=r['cracks'],.85
                r2=strike(p)
                self.assertEqual(len(r2['cracks']),2)
                self.assertLess(r2['residual'],1e-7)

    def test_shifted_impact_changes_loading_without_snapping(self):
        p=inputs()
        initial=strike(p)['cracks']
        p['cracks']=initial
        center=strike(p)
        p['aim']=.6
        shifted=strike(p)
        # A small load shift need not change the winning discrete direction.
        # It must change the calculated driving energy and stress field.
        self.assertNotAlmostEqual(center['ratio'],shifted['ratio'],places=4)
        self.assertFalse(np.allclose(center['field']['values'],shifted['field']['values'],rtol=1e-3))

    def test_boundary_and_intersection_cutoffs(self):
        a=np.array([.1,.19]); b=np.array([.08,.21])
        endpoint,terminal=trial_endpoint(a,b,[],0,.002)
        self.assertTrue(terminal); self.assertAlmostEqual(endpoint[1],HEIGHT-.002)
        other=crack([[.3,.1],[.1,.1]])
        endpoint,terminal=trial_endpoint(np.array([.2,.08]),np.array([.18,.12]),[other],1,.002)
        self.assertTrue(terminal); self.assertLess(endpoint[1],.1)
        p=inputs(); r=strike(p)
        for _ in range(12):
            if r['fractured']:break
            p['cracks']=r['cracks']; r=strike(p)
        self.assertTrue(r['fractured'])
        self.assertLess(r['cracks'][0]['points'][-1][1],.004)

    def test_reported_length_is_arc_length(self):
        c=crack([[.3,.1],[.26,.13],[.22,.1]])
        self.assertAlmostEqual(length(c),.1)

    def test_invalid_inputs_fail_cleanly(self):
        p=inputs()
        for key,value in [('young',0),('poisson',.6),('toughness',float('nan')),('aim',1),('variation',2),('cracks',[crack([[.3,.1],[.2,.1],[.25,.12]])])]:
            with self.subTest(key=key),self.assertRaises(ValueError):
                strike({**p,key:value})


if __name__=='__main__':
    unittest.main()
