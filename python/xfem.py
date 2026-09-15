"""Educational antiplane (Mode III) XFEM on a fixed triangular mesh.

SI units. Scalar displacement w(x,y), shear modulus mu=E/[2(1+nu)].
Shifted Heaviside and sqrt(r)*sin(theta/2) tip enrichments, subdivided
integration, fixed left boundary, equivalent right-edge tearing loads.
See README.md for the deliberately restricted fracture/impact model.
"""
import json
import math
import time
import numpy as np

WIDTH, HEIGHT, THICKNESS = .300, .200, .005
MESHES = {'coarse': (12, 8), 'medium': (18, 12), 'fine': (24, 16)}


def cross(a, b):
    return a[0]*b[1]-a[1]*b[0]


def clip(poly, axis, value, positive):
    result = []
    for a, b in zip(poly, poly[1:]+poly[:1]):
        da, db = (a[axis]-value)*positive, (b[axis]-value)*positive
        if da >= 0:
            result.append(a)
        if (da > 0 and db < 0) or (da < 0 and db > 0):
            result.append(a+(b-a)*da/(da-db))
    return result


def quadrature(vertices, cracks, h):
    polys = [list(vertices)]
    for c in cracks:
        # Split along each discontinuity, and through tips. This is integration
        # subdivision only: no change to the physical finite element mesh.
        for axis, value in [(1, c['y']), (0, WIDTH-c['length'])]:
            new = []
            for poly in polys:
                vals = [p[axis]-value for p in poly]
                if min(vals) < -1e-13 and max(vals) > 1e-13:
                    new.extend([clip(poly, axis, value, 1), clip(poly, axis, value, -1)])
                else:
                    new.append(poly)
            polys = [p for p in new if len(p) >= 3]
    triangles = []
    for p in polys:
        triangles.extend(np.array([p[0], p[i], p[i+1]]) for i in range(1, len(p)-1))
    near = any(np.linalg.norm(vertices.mean(axis=0)-[WIDTH-c['length'], c['y']]) < 2*h for c in cracks)
    if near:
        refined = []
        for a, b, c in triangles:
            ab, bc, ca = (a+b)/2, (b+c)/2, (c+a)/2
            refined.extend([np.array(t) for t in [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]])
        triangles = refined
    points, weights = [], []
    for tri in triangles:
        area = abs(cross(tri[1]-tri[0], tri[2]-tri[0]))/2
        if area < 1e-18:
            continue
        for bary in [[2/3, 1/6, 1/6], [1/6, 2/3, 1/6], [1/6, 1/6, 2/3]]:
            points.append(np.array(bary) @ tri)
            weights.append(area/3)
    return np.array(points), np.array(weights)


def enrichment(points, crack, kind, h):
    points = np.atleast_2d(points)
    if kind == 'H':
        return np.where(points[:, 1] >= crack['y'], 1., -1.), np.zeros_like(points)
    X = WIDTH-crack['length']-points[:, 0]
    Y = crack['y']-points[:, 1]
    r = np.maximum(np.hypot(X, Y), 1e-14)
    theta = np.arctan2(Y, X)
    value = np.sqrt(r/h)*np.sin(theta/2)
    gradient = np.column_stack((np.sin(theta/2), -np.cos(theta/2)))/(2*np.sqrt(r*h))[:, None]
    return value, gradient


class Model:
    def __init__(self, mesh='medium'):
        nx, ny = MESHES[mesh]
        self.mesh = mesh
        self.h = WIDTH/nx
        self.nodes = np.array([(x, y) for y in np.linspace(0, HEIGHT, ny+1) for x in np.linspace(0, WIDTH, nx+1)])
        self.triangles = []
        for j in range(ny):
            for i in range(nx):
                a = j*(nx+1)+i
                self.triangles.extend([[a, a+1, a+nx+2], [a, a+nx+2, a+nx+1]])
        self.triangles = np.array(self.triangles)
        self.support = [set() for _ in self.nodes]
        self.geometry = []
        for ids in self.triangles:
            v = self.nodes[ids]
            inv = np.linalg.inv(np.column_stack((np.ones(3), v)))
            self.geometry.append((v, inv, inv[1:, :].T))
            for n in ids:
                self.support[n].update(ids)
        self.bounds = [(self.nodes[list(s)].min(axis=0), self.nodes[list(s)].max(axis=0)) for s in self.support]

    def solve(self, cracks, mu, impact_y):
        n = len(self.nodes)
        enriched = {}
        ndof = n
        for ci, c in enumerate(cracks):
            tip = np.array([WIDTH-c['length'], c['y']])
            for i, p in enumerate(self.nodes):
                if p[0] < 1e-12:
                    continue
                lo, hi = self.bounds[i]
                kind = None
                if np.linalg.norm(p-tip) <= 2*self.h:
                    kind = 'tip'
                elif lo[0] > tip[0]+1e-12 and lo[1] < c['y'] < hi[1]:
                    kind = 'H'
                if kind:
                    enriched[(i, ci)] = (ndof, kind)
                    ndof += 1
        K = np.zeros((ndof, ndof))
        f = np.zeros(ndof)
        records = []
        for ei, ids in enumerate(self.triangles):
            vertices, inv, grads = self.geometry[ei]
            entries = [(local, ci, *enriched[(int(node), ci)]) for local, node in enumerate(ids) for ci in range(len(cracks)) if (int(node), ci) in enriched]
            q, weights = quadrature(vertices, cracks if entries else [], self.h)
            shape = np.column_stack((np.ones(len(q)), q)) @ inv
            B = np.repeat(grads[None, :, :], len(q), axis=0)
            dofs = list(ids)
            for local, ci, dof, kind in entries:
                value, grad = enrichment(q, cracks[ci], kind, self.h)
                nodal_value = enrichment(vertices[local], cracks[ci], kind, self.h)[0][0]
                ext = (value-nodal_value)[:, None]*grads[local]+shape[:, local, None]*grad
                B = np.concatenate((B, ext[:, None, :]), axis=1)
                dofs.append(dof)
            local_k = mu*THICKNESS*np.einsum('qid,qjd,q->ij', B, B, weights)
            K[np.ix_(dofs, dofs)] += local_k
            records.append((dofs, B, weights))
        # Two smooth opposite tractions represent an idealized tearing impulse.
        right = np.flatnonzero(self.nodes[:, 0] > WIDTH-1e-12)
        offset = self.nodes[right, 1]-impact_y
        width = .035
        positive = np.maximum(offset, 0)*np.exp(-(offset/width)**2)
        negative = np.maximum(-offset, 0)*np.exp(-(offset/width)**2)
        f[right] = positive/positive.sum()-negative/negative.sum()
        fixed = np.flatnonzero(self.nodes[:, 0] < 1e-12)
        free = np.setdiff1d(np.arange(ndof), fixed)
        A = K[np.ix_(free, free)]
        scale = np.sqrt(np.maximum(np.diag(A), 1e-30))
        A = A / scale[:, None] / scale[None, :]
        # Negligible diagonal stabilization for almost-vanishing cut supports.
        A.flat[::len(A)+1] += 1e-11
        rhs = f[free]/scale
        u = np.zeros(ndof)
        u[free] = np.linalg.solve(A, rhs)/scale
        residual = np.linalg.norm((K@u-f)[free])/np.linalg.norm(f[free])
        if not np.isfinite(u).all() or residual > 1e-4:
            raise ValueError('The enriched system could not be resolved. Reset or choose a coarser mesh.')
        energy = float(f@u/2)
        values = []
        for dofs, B, weights in records:
            gradient = np.einsum('qid,i->qd', B, u[dofs])
            values.append(float(np.sum(np.linalg.norm(gradient, axis=1)*weights)/sum(weights)*mu))
        return {'energy': energy, 'values': values, 'enriched_dofs': ndof-n, 'residual': float(residual)}


_models = {}


def strike(payload):
    started = time.perf_counter()
    mesh = payload.get('mesh', 'medium')
    if mesh not in MESHES:
        raise ValueError('Choose coarse, medium, or fine.')
    young, nu, gc = float(payload['young']), float(payload['poisson']), float(payload['toughness'])
    mass, height, aim = float(payload['mass']), float(payload['height']), float(payload['aim'])
    if not all(math.isfinite(v) for v in [young, nu, gc, mass, height, aim]):
        raise ValueError('All material and hammer inputs must be finite numbers.')
    if not (.01 <= young <= 1000 and -.9 <= nu <= .49 and .1 <= gc <= 100000 and .01 <= mass <= 10 and .001 <= height <= 1 and .08 <= aim <= .92):
        raise ValueError('One or more inputs are outside the supported range.')
    model = _models.setdefault(mesh, None)
    if model is None:
        model = _models[mesh] = Model(mesh)
    cracks = [{'y': float(c['y']), 'length': float(c['length'])} for c in payload.get('cracks', [])]
    if len(cracks) > 12 or any(not (0 < c['y'] < HEIGHT and 0 < c['length'] < WIDTH) for c in cracks):
        raise ValueError('Invalid crack geometry; reset the specimen.')
    mu = young*1e9/(2*(1+nu))
    impact_y = aim*HEIGHT
    energy = mass*9.81*height
    force = energy/.0005
    match = next((i for i, c in enumerate(cracks) if abs(c['y']-impact_y) < 1.5*model.h), None)
    # Avoid an exact nodal crack plane to reduce degenerate integration support.
    if match is None:
        if len(cracks) >= 12:
            raise ValueError('The specimen already contains the maximum number of separate cracks.')
        y = impact_y
        if min(abs(model.nodes[:, 1]-y)) < 1e-9:
            y += .017*model.h
    else:
        y = cracks[match]['y']
    base = model.solve(cracks, mu, y)
    trial = [dict(c) for c in cracks]
    increment = .65*model.h
    if match is None:
        increment = 1.3*model.h
        trial.append({'y': y, 'length': increment})
    else:
        increment = min(increment, WIDTH-trial[match]['length']-.2*model.h)
        trial[match]['length'] += increment
    candidate = model.solve(trial, mu, y)
    released = max(0., (candidate['energy']-base['energy'])*force**2)
    driving = released/(THICKNESS*increment)
    accepted = driving >= gc and energy >= gc*THICKNESS*increment
    fractured = False
    if accepted:
        cracks = trial
        field = candidate
        active = len(cracks)-1 if match is None else match
        if WIDTH-cracks[active]['length'] <= model.h:
            cracks[active]['length'] = WIDTH
            fractured = True
        message = 'Specimen separated. Reset to start a new experiment.' if fractured else ('New crack formed. Strike again to extend it, or aim elsewhere.' if match is None else 'Crack extended. Strike again to continue propagation.')
    else:
        field = base
        message = 'No growth: increase the drop height or use a heavier hammer.'
    return {'cracks': cracks, 'fractured': fractured, 'ratio': driving/gc, 'message': message, 'accepted': accepted, 'enriched_dofs': field['enriched_dofs'], 'residual': field['residual'], 'seconds': time.perf_counter()-started, 'field': {'nodes': model.nodes.tolist(), 'triangles': model.triangles.tolist(), 'values': [v*force/1e6 for v in field['values']]}}


def run_json(payload):
    return json.dumps(strike(json.loads(payload)), allow_nan=False)
