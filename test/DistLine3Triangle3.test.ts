import { describe, expect, it } from 'vitest';
import { DistLine3Triangle3 } from '../src/DistLine3Triangle3';
import { Line } from '../src/Line';
import { Triangle } from '../src/Triangle';
import { Vector, add, dot, mul, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistLine3Triangle3', () => {
    const query = new DistLine3Triangle3();
    // The triangle in the z = 0 plane with vertices (0,0,0), (1,0,0),
    // (0,1,0).
    const tri = Triangle.fromVertices(v(0, 0, 0), v(1, 0, 0), v(0, 1, 0));

    it('reports zero distance for a line piercing the triangle', () => {
        const result = query.compute(line([0.25, 0.25, -3], [0, 0, 1]), tri);
        expect(result.distance).toBe(0);
        expect(result.parameter).toBeCloseTo(3, 12);
        expect(result.barycentric[0]).toBeCloseTo(0.5, 12);
        expect(result.barycentric[1]).toBeCloseTo(0.25, 12);
        expect(result.barycentric[2]).toBeCloseTo(0.25, 12);
        expect(result.closest[0].values[2]).toBeCloseTo(0, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(0, 12);
    });

    it('measures a line parallel to the plane of the triangle', () => {
        const result = query.compute(line([0, 0, 4], [1, 0, 0]), tri);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(0, 12);
    });

    it('measures a line whose plane intersection is outside the triangle',
        () => {
            // The line pierces the plane at (2,2,0), outside the triangle.
            const result = query.compute(line([2, 2, -1], [0, 0, 1]), tri);
            // The closest triangle point is the midpoint region of the
            // hypotenuse, at (0.5,0.5,0); the distance is sqrt(2*1.5^2).
            expect(result.distance).toBeCloseTo(Math.sqrt(2) * 1.5, 10);
            expect(result.closest[1].values[0]).toBeCloseTo(0.5, 10);
            expect(result.closest[1].values[1]).toBeCloseTo(0.5, 10);
        });

    it('reports the barycentric coordinates of an edge closest point', () => {
        const result = query.compute(line([2, 2, -1], [0, 0, 1]), tri);
        const sum = result.barycentric[0] + result.barycentric[1]
            + result.barycentric[2];
        expect(sum).toBeCloseTo(1, 10);
        for (const b of result.barycentric) {
            expect(b).toBeGreaterThanOrEqual(-1e-12);
        }
        // The closest point reconstructed from the barycentric coordinates.
        const q = add(mul(result.barycentric[0], tri.v[0]),
            add(mul(result.barycentric[1], tri.v[1]),
                mul(result.barycentric[2], tri.v[2])));
        for (let i = 0; i < 3; ++i) {
            expect(q.values[i]).toBeCloseTo(result.closest[1].values[i], 9);
        }
    });

    it('reports zero distance for a line lying in an edge of the triangle',
        () => {
            const result = query.compute(line([0, 0, 0], [1, 0, 0]), tri);
            expect(result.distance).toBeCloseTo(0, 12);
        });

    it('handles a degenerate triangle collapsed to a segment', () => {
        const degenerate = Triangle.fromVertices(v(0, 0, 0), v(1, 0, 0),
            v(2, 0, 0));
        const result = query.compute(line([0, 3, 0], [1, 0, 0]), degenerate);
        expect(result.distance).toBeCloseTo(3, 10);
    });

    it('agrees with a dense sampling of the triangle', () => {
        const rnd = makeRandom(24689);
        const t = Triangle.fromVertices(v(0.5, -1, 0.25), v(2, 0.5, -0.5),
            v(-1, 1.5, 1));

        for (let trial = 0; trial < 30; ++trial) {
            const origin = v(6 * rnd() - 3, 6 * rnd() - 3, 6 * rnd() - 3);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const ln = Line.fromOriginDirection(origin, dir);
            const result = query.compute(ln, t);

            // The barycentric coordinates are valid and reproduce closest[1].
            const b = result.barycentric;
            expect(b[0] + b[1] + b[2]).toBeCloseTo(1, 8);
            for (const bi of b) {
                expect(bi).toBeGreaterThanOrEqual(-1e-8);
            }
            const q = add(mul(b[0], t.v[0]),
                add(mul(b[1], t.v[1]), mul(b[2], t.v[2])));
            for (let i = 0; i < 3; ++i) {
                expect(q.values[i]).toBeCloseTo(result.closest[1].values[i],
                    7);
            }

            // The reported line point matches the reported parameter.
            const onLine = add(ln.origin, mul(result.parameter, ln.direction));
            for (let i = 0; i < 3; ++i) {
                expect(onLine.values[i]).toBeCloseTo(
                    result.closest[0].values[i], 7);
            }

            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 7);

            // No sampled triangle point is closer to the line.
            const n = 60;
            const dd = dot(ln.direction, ln.direction);
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= n; ++i) {
                for (let j = 0; i + j <= n; ++j) {
                    const b1 = i / n, b2 = j / n, b0 = 1 - b1 - b2;
                    const p = add(mul(b0, t.v[0]),
                        add(mul(b1, t.v[1]), mul(b2, t.v[2])));
                    const w = sub(p, ln.origin);
                    const s = dot(w, ln.direction) / dd;
                    const f = sub(w, mul(s, ln.direction));
                    best = Math.min(best, dot(f, f));
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});
