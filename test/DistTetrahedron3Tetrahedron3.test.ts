import { describe, it, expect } from 'vitest';
import { inContainerTetrahedron3 } from '../src/ContTetrahedron3.js';
import { DistPoint3Tetrahedron3 } from '../src/DistPoint3Tetrahedron3.js';
import { DistTetrahedron3Tetrahedron3 } from '../src/DistTetrahedron3Tetrahedron3.js';
import { Tetrahedron3 } from '../src/Tetrahedron3.js';
import { DistTriangle3Triangle3 } from '../src/DistTriangle3Triangle3.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, latticeVector, rotationFrame,
    seededRandom, wellScaledVector
} from './helpers/arbitraries.js';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// The reference tetrahedron uses the vertex ordering that Tetrahedron3 and
// inContainerTetrahedron3 expect.
function unitTetra(offset: Vector = v(0, 0, 0), s = 1): Tetrahedron3 {
    return Tetrahedron3.fromVertices(
        add(offset, v(0, 0, 0)), add(offset, v(s, 0, 0)),
        add(offset, v(0, s, 0)), add(offset, v(0, 0, s)));
}

// Brute-force minimum from points sampled on the faces of tetra0 to the solid
// tetra1, using the independently ported point-tetrahedron query.
function bruteForce(tetra0: Tetrahedron3, tetra1: Tetrahedron3,
    n: number): number {
    const pointQuery = new DistPoint3Tetrahedron3();
    let best = Infinity;
    for (let face = 0; face < 4; ++face) {
        const indices = Tetrahedron3.getFaceIndices(face);
        for (let i = 0; i <= n; ++i) {
            for (let j = 0; i + j <= n; ++j) {
                const b0 = i / n, b1 = j / n, b2 = 1 - b0 - b1;
                const p = add(add(mul(b0, tetra0.v[indices[0]]),
                    mul(b1, tetra0.v[indices[1]])),
                    mul(b2, tetra0.v[indices[2]]));
                const d = pointQuery.compute(p, tetra1).distance;
                if (d < best) {
                    best = d;
                }
            }
        }
    }
    return best;
}

describe('DistTetrahedron3Tetrahedron3', () => {
    it('reports the gap between two separated tetrahedra', () => {
        // The two tetrahedra are translates along x. The right face of the
        // first is the plane x = 0 shifted, so the closest points lie on the
        // vertex (1,0,0) of the first and (4,0,0) of the second.
        const tetra0 = unitTetra();
        const tetra1 = unitTetra(v(4, 0, 0));
        const result = new DistTetrahedron3Tetrahedron3()
            .compute(tetra0, tetra1);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.sqrDistance).toBeCloseTo(9, 12);
        expect(result.closest[0].values[0]).toBeCloseTo(1, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(4, 10);
    });

    it('reports zero for tetrahedra that share a face', () => {
        const tetra0 = unitTetra();
        const tetra1 = Tetrahedron3.fromVertices(
            v(0, 0, 0), v(0, 0, 1), v(0, 1, 0), v(-1, 0, 0));
        const result = new DistTetrahedron3Tetrahedron3()
            .compute(tetra0, tetra1);
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
        expect(length(sub(result.closest[0], result.closest[1])))
            .toBeCloseTo(0, 12);
    });

    it('reports zero for overlapping tetrahedra', () => {
        const tetra0 = unitTetra(v(0, 0, 0), 2);
        const tetra1 = unitTetra(v(0.2, 0.2, 0.2), 2);
        const result = new DistTetrahedron3Tetrahedron3()
            .compute(tetra0, tetra1);
        expect(result.distance).toBe(0);
    });

    it('reports zero when one tetrahedron is nested inside the other', () => {
        // The small tetrahedron is strictly inside the large one, so no pair
        // of faces intersects and the nesting is detected by the centroid
        // containment test.
        const tetra0 = unitTetra(v(0, 0, 0), 8);
        const tetra1 = unitTetra(v(1, 1, 1), 1);
        const query = new DistTetrahedron3Tetrahedron3();
        expect(inContainerTetrahedron3(tetra1.computeCentroid(), tetra0))
            .toBe(true);

        const result = query.compute(tetra0, tetra1);
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
        expect(length(sub(result.closest[0], result.closest[1])))
            .toBeCloseTo(0, 12);
        // The reported point is the centroid of the nested tetrahedron.
        expect(length(sub(result.closest[0], tetra1.computeCentroid())))
            .toBeCloseTo(0, 12);

        // The reversed query also reports zero.
        const reversed = query.compute(tetra1, tetra0);
        expect(reversed.distance).toBe(0);
    });

    it('reports barycentric coordinates that reproduce the closest points',
        () => {
            const tetra0 = unitTetra();
            const tetra1 = unitTetra(v(3, 1, 2), 2);
            const result = new DistTetrahedron3Tetrahedron3()
                .compute(tetra0, tetra1);

            const from = (t: Tetrahedron3, b: readonly number[]) => {
                let p = new Vector(3);
                for (let i = 0; i < 4; ++i) {
                    p = add(p, mul(b[i], t.v[i]));
                }
                return p;
            };
            expect(result.barycentric0[0] + result.barycentric0[1] +
                result.barycentric0[2] + result.barycentric0[3])
                .toBeCloseTo(1, 10);
            expect(length(sub(from(tetra0, result.barycentric0),
                result.closest[0]))).toBeCloseTo(0, 8);
            expect(length(sub(from(tetra1, result.barycentric1),
                result.closest[1]))).toBeCloseTo(0, 8);
        });

    it('is symmetric in its arguments', () => {
        const tetra0 = Tetrahedron3.fromVertices(
            v(0, 0, 0), v(2, 0.5, 0), v(0.5, 2, 0.25), v(0.25, 0.5, 3));
        const tetra1 = Tetrahedron3.fromVertices(
            v(5, 5, 5), v(7, 5.5, 5), v(5.5, 7, 5.25), v(5.25, 5.5, 8));
        const query = new DistTetrahedron3Tetrahedron3();
        const forward = query.compute(tetra0, tetra1);
        const backward = query.compute(tetra1, tetra0);
        expect(backward.distance).toBeCloseTo(forward.distance, 10);
        expect(length(sub(backward.closest[1], forward.closest[0])))
            .toBeCloseTo(0, 8);
        expect(length(sub(backward.closest[0], forward.closest[1])))
            .toBeCloseTo(0, 8);
    });

    it('agrees with brute-force sampling and reports consistent points', () => {
        let seed = 192837465;
        const rnd = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };
        const query = new DistTetrahedron3Tetrahedron3();
        const pointQuery = new DistPoint3Tetrahedron3();
        for (let trial = 0; trial < 25; ++trial) {
            // Randomly perturbed copies of the reference tetrahedron keep the
            // vertex ordering that the containment test assumes.
            const mk = (shift: Vector, scale: number) => {
                const jitter = () =>
                    v(rnd() * 0.4 - 0.2, rnd() * 0.4 - 0.2, rnd() * 0.4 - 0.2);
                const base = unitTetra(shift, scale);
                return Tetrahedron3.fromArray(
                    base.v.map(p => add(p, jitter())));
            };
            const tetra0 = mk(v(0, 0, 0), 1 + rnd());
            const tetra1 = mk(
                v(rnd() * 5 - 1, rnd() * 5 - 1, rnd() * 5 - 1), 1 + rnd());

            const result = query.compute(tetra0, tetra1);
            const brute = Math.min(bruteForce(tetra0, tetra1, 12),
                bruteForce(tetra1, tetra0, 12));
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-9);
            expect(result.distance).toBeGreaterThan(brute - 0.3);
            expect(result.sqrDistance)
                .toBeCloseTo(result.distance * result.distance, 9);
            expect(length(sub(result.closest[0], result.closest[1])))
                .toBeCloseTo(result.distance, 8);
            expect(pointQuery.compute(result.closest[0], tetra0).distance)
                .toBeCloseTo(0, 8);
            expect(pointQuery.compute(result.closest[1], tetra1).distance)
                .toBeCloseTo(0, 8);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against upstream DistTetrahedron3Tetrahedron3.h.
// ---------------------------------------------------------------------------

describe('DistTetrahedron3Tetrahedron3 verification', () => {
    const query = new DistTetrahedron3Tetrahedron3();
    const ttQuery = new DistTriangle3Triangle3();

    const tetraArb = fc.tuple(latticeVector(3, -5, 5), latticeVector(3, -5, 5),
        latticeVector(3, -5, 5), latticeVector(3, -5, 5))
        .filter(([a, b, c, d]) =>
            Math.abs(dot(sub(b, a), cross(sub(c, a), sub(d, a)))) > 4)
        .map(([a, b, c, d]) => Tetrahedron3.fromVertices(a, b, c, d));

    function faces(t: Tetrahedron3): Triangle[] {
        const out: Triangle[] = [];
        for (let f = 0; f < 4; ++f) {
            const idx = Tetrahedron3.getFaceIndices(f);
            out.push(Triangle.fromVertices(t.v[idx[0]], t.v[idx[1]],
                t.v[idx[2]]));
        }
        return out;
    }

    it('reports consistent distances and barycentric closest points', () => {
        check(fc.tuple(tetraArb, tetraArb), ([t0, t1]) => {
            const r = query.compute(t0, t1);
            // The absolute tolerance is 1e-6: these queries accumulate the
            // squared distance while clamping to faces and edges, so a
            // near-touching configuration loses about half the mantissa and
            // the distance carries an absolute error of order sqrt(eps)
            // times the coordinate scale. A translation or frame error
            // would show up as an O(1) discrepancy.
            expectClose(r.sqrDistance, r.distance * r.distance, 1e-12, 1e-12);
            const d = sub(r.closest[0], r.closest[1]);
            expectClose(Math.sqrt(dot(d, d)), r.distance, 1e-6, 1e-8);

            const verify = (b: [number, number, number, number],
                t: Tetrahedron3, c: Vector) => {
                expectClose(b[0] + b[1] + b[2] + b[3], 1, 1e-8, 1e-8);
                let rebuilt = new Vector(3);
                for (let i = 0; i < 4; ++i) {
                    expect(b[i]).toBeGreaterThanOrEqual(-1e-7);
                    rebuilt = add(rebuilt, mul(b[i], t.v[i]));
                }
                expectVectorClose(rebuilt, c, 1e-6, 1e-6);
            };
            verify(r.barycentric0, t0, r.closest[0]);
            verify(r.barycentric1, t1, r.closest[1]);
        });
    });

    it('is symmetric under argument swap', () => {
        check(fc.tuple(tetraArb, tetraArb), ([t0, t1]) => {
            const a = query.compute(t0, t1);
            const b = query.compute(t1, t0);
            expectClose(a.distance, b.distance, 1e-8, 1e-8);
        });
    });

    // When neither solid contains the other, the minimum distance is attained
    // on the boundaries, which are the four triangular faces of each solid.
    it('equals the minimum over the sixteen face pairs when not nested',
        () => {
            check(fc.tuple(tetraArb, tetraArb), ([t0, t1]) => {
                const c0 = t0.computeCentroid(), c1 = t1.computeCentroid();
                if (inContainerTetrahedron3(c0, t1)
                    || inContainerTetrahedron3(c1, t0)) {
                    return;   // nested; the face minimum is not the answer
                }
                let best = Number.MAX_VALUE;
                for (const f0 of faces(t0)) {
                    for (const f1 of faces(t1)) {
                        best = Math.min(best,
                            ttQuery.compute(f0, f1).sqrDistance);
                    }
                }
                const r = query.compute(t0, t1);
                expectClose(r.sqrDistance, best, 1e-8, 1e-8);
            });
        });

    it('reports zero when one tetrahedron is nested inside the other', () => {
        const outer = unitTetra(v(0, 0, 0), 10);
        const inner = unitTetra(v(1, 1, 1), 1);
        const r = query.compute(inner, outer);
        expect(r.distance).toBe(0);
        expect(r.sqrDistance).toBe(0);
        expectVectorClose(r.closest[0], r.closest[1], 1e-12, 1e-12);
        const rs = query.compute(outer, inner);
        expect(rs.distance).toBe(0);
    });

    it('never exceeds a face sampling of the first tetrahedron', () => {
        const rng = seededRandom(0x0ddba11);
        for (let k = 0; k < 15; ++k) {
            const p = () => v(8 * rng() - 4, 8 * rng() - 4, 8 * rng() - 4);
            const t0 = Tetrahedron3.fromVertices(p(), p(), p(), p());
            const t1 = Tetrahedron3.fromVertices(p(), p(), p(), p());
            const r = query.compute(t0, t1);
            expect(r.distance)
                .toBeLessThanOrEqual(bruteForce(t0, t1, 12) + 1e-8);
        }
    }, 30000);

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(tetraArb, tetraArb, rotationFrame(3),
            wellScaledVector(3, -4, 4)), ([t0, t1, R, tr]) => {
                const xf = (q: Vector) => add(add(add(
                    mul(q.values[0], R[0]), mul(q.values[1], R[1])),
                    mul(q.values[2], R[2])), tr);
                const xt = (t: Tetrahedron3) => Tetrahedron3.fromVertices(
                    xf(t.v[0]), xf(t.v[1]), xf(t.v[2]), xf(t.v[3]));
                const a = query.compute(t0, t1);
                const b = query.compute(xt(t0), xt(t1));
                expectClose(a.distance, b.distance, 1e-7, 1e-7);
            });
    });
});
