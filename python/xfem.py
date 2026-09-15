"""Educational Mode III XFEM with energy-selected, piecewise-linear crack paths.

SI units; scalar antiplane elasticity, not a dynamic contact calculation.
Cracks are nonintersecting, left-monotone polylines. Shifted Heaviside
and tip enrichments enter the elastic stiffness matrix. Candidate directions
are compared by virtual-extension energy release divided by local toughness.
"""
import json
import math
import time
import numpy as np

WIDTH, HEIGHT, THICKNESS = .300, .200, .005
MESHES = {'coarse': (12, 8), 'medium': (18, 12), 'fine': (24, 16)}
TURNS = (0, -25, 25, -50, 50)
MAX_ADVANCES = 3


def cross(a, b):
    return a[0]*b[1]-a[1]*b[0]


def path(crack):
    return np.asarray(crack['points'], dtype=float)


def length(crack):
    return float(np.linalg.norm(np.diff(path(crack), axis=0), axis=1).sum())


def side(points, crack):
    """Continuous graph separates faces; continuation ahead uses tip tangent."""
    points = np.atleast_2d(points)
    p = path(crack)
    yc = np.interp(points[:, 0], p[::-1, 0], p[::-1, 1])
    slope = (p[-1, 1]-p[-2, 1])/(p[-1, 0]-p[-2, 0])
    ahead = points[:, 0] < p[-1, 0]
    yc[ahead] = p[-1, 1]+slope*(points[ahead, 0]-p[-1, 0])
    # The left-pointing tangent has a downward-facing normal.
    return np.where(points[:, 1] <= yc, 1., -1.)


def clip_line(poly, normal, offset, positive):
    result = []
    for a, b in zip(poly, poly[1:]+poly[:1]):
        da, db = (a@normal-offset)*positive, (b@normal-offset)*positive
        if da >= 0:
            result.append(a)
        if da*db < 0:
            result.append(a+(b-a)*da/(da-db))
    return result


def quadrature(vertices, cracks, h):
    polys = [list(vertices)]
    lo, hi = vertices.min(axis=0), vertices.max(axis=0)
    lines = []
    for c in cracks:
        p = path(c)
        for a, b in zip(p[:-1], p[1:]):
            if np.any(np.maximum(a, b) < lo-1e-12) or np.any(np.minimum(a, b) > hi+1e-12):
                continue
            tangent = b-a
            normal = np.array([-tangent[1], tangent[0]])/np.linalg.norm(tangent)
            lines.append((normal, float(a@normal)))
            # Partition changes of the graph segment at kinks as well.
            if lo[0] < b[0] < hi[0]:
                lines.append((np.array([1., 0.]), b[0]))
        tip = p[-1]
        if np.linalg.norm(vertices.mean(axis=0)-tip) < 2*h:
            tangent = (p[-1]-p[-2])/np.linalg.norm(p[-1]-p[-2])
            lines.append((tangent, float(tip@tangent)))
    for normal, offset in lines:
        new = []
        for poly in polys:
            vals = np.array(poly)@normal-offset
            if vals.min() < -1e-13 and vals.max() > 1e-13:
                new.extend([clip_line(poly, normal, offset, 1), clip_line(poly, normal, offset, -1)])
            else:
                new.append(poly)
        polys = [p for p in new if len(p) >= 3]
    triangles = []
    for p in polys:
        triangles.extend(np.array([p[0], p[i], p[i+1]]) for i in range(1, len(p)-1))
    near = any(np.linalg.norm(vertices.mean(axis=0)-path(c)[-1]) < 2*h for c in cracks)
    if near:
        refined = []
        for a, b, c in triangles:
            ab, bc, ca = (a+b)/2, (b+c)/2, (c+a)/2
            refined.extend(np.array(t) for t in [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)])
        triangles = refined
    points, weights = [], []
    for tri in triangles:
        area = abs(cross(tri[1]-tri[0], tri[2]-tri[0]))/2
        if area < 1e-18:
            continue
        for bary in [[2/3, 1/6, 1/6], [1/6, 2/3, 1/6], [1/6, 1/6, 2/3]]:
            points.append(np.array(bary)@tri)
            weights.append(area/3)
    return np.array(points), np.array(weights)


def enrichment(points, crack, kind, h):
    points = np.atleast_2d(points)
    if kind == 'H':
        return side(points, crack), np.zeros_like(points)
    p = path(crack)
    tangent = (p[-1]-p[-2])/np.linalg.norm(p[-1]-p[-2])
    normal = np.array([-tangent[1], tangent[0]])
    relative = points-p[-1]
    X, Y = relative@tangent, relative@normal
    r = np.maximum(np.hypot(X, Y), 1e-14)
    theta = np.arctan2(Y, X)
    # Route the jump along the actual polyline, not the last segment's
    # infinitely extended line. Ahead of the tip use its straight tangent.
    face = side(points, crack)
    face[points[:, 0] < p[-1, 0]] = np.where(Y[points[:, 0] < p[-1, 0]] >= 0, 1., -1.)
    adjustment = face*np.where(Y >= 0, 1., -1.)
    value = adjustment*np.sqrt(r/h)*np.sin(theta/2)
    gradient = adjustment[:, None]*(-np.sin(theta/2)[:, None]*tangent+np.cos(theta/2)[:, None]*normal)/(2*np.sqrt(r*h))[:, None]
    return value, gradient


def crosses_box(a, b, lo, hi):
    """Liang–Barsky clipping for enrichment-support intersection."""
    low, high = 0., 1.
    for axis in range(2):
        delta = b[axis]-a[axis]
        if abs(delta) < 1e-14:
            if a[axis] < lo[axis] or a[axis] > hi[axis]:
                return False
        else:
            v0, v1 = (lo[axis]-a[axis])/delta, (hi[axis]-a[axis])/delta
            low, high = max(low, min(v0, v1)), min(high, max(v0, v1))
    return low <= high


def resistance(points, variation):
    """Fixed smooth illustrative toughness field; never generated per strike."""
    q = np.atleast_2d(points)
    x, y = q[:, 0], q[:, 1]
    structure = .6*np.sin(80*x+23*y)+.4*np.sin(44*x-95*y+1.7)
    return np.exp(variation*structure)


class Model:
    def __init__(self, mesh='medium'):
        nx, ny = MESHES[mesh]
        self.mesh, self.h = mesh, WIDTH/nx
        self.nodes = np.array([(x, y) for y in np.linspace(0, HEIGHT, ny+1) for x in np.linspace(0, WIDTH, nx+1)])
        triangles = []
        for j in range(ny):
            for i in range(nx):
                a = j*(nx+1)+i
                triangles.extend([[a, a+1, a+nx+2], [a, a+nx+2, a+nx+1]])
        self.triangles = np.array(triangles)
        self.support = [set() for _ in self.nodes]
        self.geometry = []
        for ids in self.triangles:
            v = self.nodes[ids]
            inv = np.linalg.inv(np.column_stack((np.ones(3), v)))
            self.geometry.append((v, inv, inv[1:, :].T))
            for n in ids:
                self.support[n].update(ids)
        self.bounds = [(self.nodes[list(s)].min(axis=0), self.nodes[list(s)].max(axis=0)) for s in self.support]

    def solve(self, cracks, mu, impact_y, include_field=True):
        n, enriched, ndof = len(self.nodes), {}, len(self.nodes)
        for ci, c in enumerate(cracks):
            p = path(c)
            for i, node in enumerate(self.nodes):
                if node[0] < 1e-12:
                    continue
                lo, hi = self.bounds[i]
                kind = None
                if np.linalg.norm(node-p[-1]) <= 2*self.h:
                    kind = 'tip'
                elif any(crosses_box(a, b, lo, hi) for a, b in zip(p[:-1], p[1:])):
                    corners = np.array([[lo[0], lo[1]], [lo[0], hi[1]], [hi[0], lo[1]], [hi[0], hi[1]]])
                    signs = side(corners, c)
                    if signs.min() != signs.max():
                        kind = 'H'
                if kind:
                    enriched[(i, ci)] = (ndof, kind)
                    ndof += 1
        K, f, records = np.zeros((ndof, ndof)), np.zeros(ndof), []
        for ei, ids in enumerate(self.triangles):
            vertices, inv, grads = self.geometry[ei]
            entries = [(local, ci, *enriched[(int(node), ci)]) for local, node in enumerate(ids) for ci in range(len(cracks)) if (int(node), ci) in enriched]
            relevant = sorted(set(e[1] for e in entries))
            q, weights = quadrature(vertices, [cracks[ci] for ci in relevant], self.h)
            shape = np.column_stack((np.ones(len(q)), q))@inv
            B = np.repeat(grads[None, :, :], len(q), axis=0)
            dofs = list(ids)
            for local, ci, dof, kind in entries:
                value, grad = enrichment(q, cracks[ci], kind, self.h)
                nodal_value = enrichment(vertices[local], cracks[ci], kind, self.h)[0][0]
                ext = (value-nodal_value)[:, None]*grads[local]+shape[:, local, None]*grad
                B = np.concatenate((B, ext[:, None, :]), axis=1)
                dofs.append(dof)
            K[np.ix_(dofs, dofs)] += mu*THICKNESS*np.einsum('qid,qjd,q->ij', B, B, weights)
            if include_field:
                records.append((dofs, B, weights))
        right = np.flatnonzero(self.nodes[:, 0] > WIDTH-1e-12)
        offset = self.nodes[right, 1]-impact_y
        width = .035
        positive = np.maximum(offset, 0)*np.exp(-(offset/width)**2)
        negative = np.maximum(-offset, 0)*np.exp(-(offset/width)**2)
        f[right] = positive/positive.sum()-negative/negative.sum()
        free = np.setdiff1d(np.arange(ndof), np.flatnonzero(self.nodes[:, 0] < 1e-12))
        A = K[np.ix_(free, free)]
        scale = np.sqrt(np.maximum(np.diag(A), 1e-30))
        A = A/scale[:, None]/scale[None, :]
        A.flat[::len(A)+1] += 1e-11
        u = np.zeros(ndof)
        u[free] = np.linalg.solve(A, f[free]/scale)/scale
        residual = np.linalg.norm((K@u-f)[free])/np.linalg.norm(f[free])
        if not np.isfinite(u).all() or residual > 1e-4:
            raise ValueError('The enriched system could not be resolved. Reset or choose a coarser mesh.')
        values = []
        for dofs, B, weights in records:
            gradient = np.einsum('qid,i->qd', B, u[dofs])
            values.append(float(np.sum(np.linalg.norm(gradient, axis=1)*weights)/sum(weights)*mu))
        return {'energy': float(f@u/2), 'values': values, 'enriched_dofs': ndof-n, 'residual': float(residual)}


def segment_hit(a, b, c, d):
    denominator = cross(b-a, d-c)
    if abs(denominator) < 1e-13:
        return None
    t, u = cross(c-a, d-c)/denominator, cross(c-a, b-a)/denominator
    return t if 1e-8 < t <= 1 and -1e-8 <= u <= 1+1e-8 else None


def trial_endpoint(a, b, cracks, ci, margin):
    """Keep numerical tips away from outer edges/other cracks; mark approach."""
    fraction, terminal = 1., False
    for axis, boundary in [(0, margin), (1, margin), (1, HEIGHT-margin)]:
        delta = b[axis]-a[axis]
        if abs(delta) < 1e-14:
            continue
        t = (boundary-a[axis])/delta
        crossing = b[axis] < boundary if (axis == 0 or boundary == margin) else b[axis] > boundary
        if crossing and 0 < t <= fraction:
            fraction, terminal = t, True
    for j, c in enumerate(cracks):
        if j == ci:
            continue
        p = path(c)
        for c0, c1 in zip(p[:-1], p[1:]):
            t = segment_hit(a, b, c0, c1)
            if t is not None and t <= fraction:
                fraction = max(0., t-margin/np.linalg.norm(b-a))
                terminal = True
    return a+(b-a)*fraction, terminal


_models = {}


def strike(payload):
    started = time.perf_counter()
    mesh = payload.get('mesh', 'medium')
    if mesh not in MESHES:
        raise ValueError('Choose coarse, medium, or fine.')
    young, nu, gc = float(payload['young']), float(payload['poisson']), float(payload['toughness'])
    mass, height, aim = float(payload['mass']), float(payload['height']), float(payload['aim'])
    variation = float(payload.get('variation', .45))
    if not all(math.isfinite(v) for v in [young, nu, gc, mass, height, aim, variation]):
        raise ValueError('All inputs must be finite numbers.')
    if not (.01 <= young <= 1000 and -.9 <= nu <= .49 and .1 <= gc <= 100000 and .01 <= mass <= 10 and .001 <= height <= 1 and .08 <= aim <= .92 and 0 <= variation <= .8):
        raise ValueError('One or more inputs are outside the supported range.')
    if mesh not in _models:
        _models[mesh] = Model(mesh)
    model = _models[mesh]
    cracks = []
    for c in payload.get('cracks', []):
        p = np.asarray(c.get('points', []), dtype=float)
        if p.ndim != 2 or p.shape[1] != 2 or not 2 <= len(p) <= 200 or not np.isfinite(p).all():
            raise ValueError('Invalid crack geometry; reset the specimen.')
        if abs(p[0, 0]-WIDTH) > 1e-9 or np.any(np.diff(p[:, 0]) >= -1e-10) or np.any(p[:, 0] < 0) or np.any(p[:, 0] > WIDTH) or np.any(p[:, 1] <= 0) or np.any(p[:, 1] >= HEIGHT):
            raise ValueError('Crack paths must progress into the specimen without looping.')
        cracks.append({'points': p.tolist(), 'stopped': bool(c.get('stopped', False))})
    if len(cracks) > 8:
        raise ValueError('A maximum of eight separate cracks is supported.')
    mu, impact_y = young*1e9/(2*(1+nu)), aim*HEIGHT
    impact_energy = mass*9.81*height
    force, remaining = impact_energy/.0005, impact_energy
    events, best_ratio, tested, failed = [], 0., 0, 0
    separated = any(c['stopped'] for c in cracks)
    base = model.solve(cracks, mu, impact_y, False)
    for step in range(MAX_ADVANCES):
        options = []
        targets = [i for i, c in enumerate(cracks) if not c['stopped']]
        allow_seed = step == 0 and len(cracks) < 8 and all(abs(c['points'][0][1]-impact_y) >= 1.5*model.h for c in cracks)
        if allow_seed:
            targets.append(len(cracks))
        for ci in targets:
            is_new = ci == len(cracks)
            if is_new:
                y = impact_y
                if min(abs(model.nodes[:, 1]-y)) < 1e-9:
                    y += .017*model.h
                a, heading, increment = np.array([WIDTH, y]), math.pi, 1.3*model.h
            else:
                p = path(cracks[ci]); a = p[-1]
                heading = math.atan2(*(p[-1]-p[-2])[::-1])
                increment = .65*model.h
            for turn in TURNS:
                angle = heading+math.radians(turn)
                direction = np.array([math.cos(angle), math.sin(angle)])
                if direction[0] > -math.cos(math.radians(65)):
                    continue
                b, terminal = trial_endpoint(a, a+increment*direction, cracks, ci, .15*model.h)
                da = float(np.linalg.norm(b-a))
                if da < .15*model.h:
                    continue
                candidate_cracks = [{'points': [p[:] for p in c['points']], 'stopped': c['stopped']} for c in cracks]
                if is_new:
                    candidate_cracks.append({'points': [a.tolist(), b.tolist()], 'stopped': terminal})
                else:
                    candidate_cracks[ci]['points'].append(b.tolist())
                    candidate_cracks[ci]['stopped'] = terminal
                tested += 1
                try:
                    candidate = model.solve(candidate_cracks, mu, impact_y, False)
                except (ValueError, np.linalg.LinAlgError):
                    failed += 1
                    continue
                local_gc = gc*float(np.mean(resistance(np.array([a*.75+b*.25, (a+b)/2, a*.25+b*.75]), variation)))
                release = max(0., (candidate['energy']-base['energy'])*force**2)
                ratio, cost = release/(THICKNESS*da*local_gc), THICKNESS*da*local_gc
                options.append((ratio, -abs(turn), ci, turn, cost, candidate_cracks, candidate, is_new, terminal))
        if not options:
            break
        best = max(options, key=lambda v: (v[0], v[1]))
        ratio, _, ci, turn, cost, trial, candidate, is_new, terminal = best
        best_ratio = ratio
        if ratio < 1 or cost > remaining:
            break
        cracks, base, remaining = trial, candidate, remaining-cost
        events.append({'crack': ci+1, 'turn': turn, 'new': is_new, 'ratio': ratio, 'from': cracks[ci]['points'][-2], 'to': cracks[ci]['points'][-1]})
        if terminal:
            separated = True
            break
    field = model.solve(cracks, mu, impact_y)
    for c in cracks:
        c['length'] = length(c)
        c['y'] = c['points'][0][1]
    if separated:
        message = 'A crack reached the boundary/intersection cutoff. Fracture limit reached; reset to continue.'
    elif events:
        formed = sum(e['new'] for e in events)
        ids = ', '.join('C'+str(i) for i in sorted(set(e['crack'] for e in events)))
        message = f'{len(events)} growth steps in {ids}. ' + ('New crack formed. ' if formed else '') + 'Move the hammer to change the loading or start another crack.'
    else:
        message = 'Cracks arrested: increase the drop height or move the impact to change the loading.'
    if failed:
        message += f' {failed} unresolved trial directions were excluded.'
    return {'cracks': cracks, 'fractured': separated, 'ratio': best_ratio, 'message': message, 'accepted': bool(events), 'events': events, 'tested': tested, 'failed_trials': failed, 'enriched_dofs': field['enriched_dofs'], 'residual': field['residual'], 'seconds': time.perf_counter()-started, 'field': {'nodes': model.nodes.tolist(), 'triangles': model.triangles.tolist(), 'values': [v*force/1e6 for v in field['values']]}}


def run_json(payload):
    return json.dumps(strike(json.loads(payload)), allow_nan=False)
