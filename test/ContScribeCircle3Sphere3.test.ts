import { describe, it, expect } from 'vitest';
import {
    circumscribeCircle3,
    circumscribeSphere3,
    inscribeCircle3,
    inscribeSphere3
} from '../src/ContScribeCircle3Sphere3.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// Distance from C to the line through p0 and p1.
function lineDistance3(c: Vector, p0: Vector, p1: Vector): number {
    const e = sub(p1, p0);
    return length(cross(e, sub(c, p0))) / length(e);
}

// Distance from C to the plane through p0, p1, p2.
function planeDistance(c: Vector, p0: Vector, p1: Vector, p2: Vector): number {
    const n = cross(sub(p1, p0), sub(p2, p0));
    return Math.abs(dot(n, sub(c, p0))) / length(n);
}

let seed = 13572468;
function rand(): number {
    seed = (seed * 1103515245 + 12345) % 2147483648;
    return seed / 2147483648;
}
function rv(scale = 10): Vector {
    return v3(scale * (2 * rand() - 1), scale * (2 * rand() - 1),
        scale * (2 * rand() - 1));
}

describe('ContScribeCircle3Sphere3', () => {
    describe('circumscribeCircle3', () => {
        it('computes the circumcircle of a right triangle in the z = 2 plane',
            () => {
                const circle = circumscribeCircle3(v3(0, 0, 2), v3(3, 0, 2),
                    v3(0, 4, 2));
                expect(circle).not.toBeNull();
                expect(circle!.center.values[0]).toBeCloseTo(1.5, 12);
                expect(circle!.center.values[1]).toBeCloseTo(2, 12);
                expect(circle!.center.values[2]).toBeCloseTo(2, 12);
                expect(circle!.radius).toBeCloseTo(2.5, 12);
                expect(Math.abs(circle!.normal.values[2])).toBeCloseTo(1, 12);
                expect(length(circle!.normal)).toBeCloseTo(1, 12);
            });

        it('returns null for collinear points', () => {
            expect(circumscribeCircle3(v3(0, 0, 0), v3(1, 1, 1), v3(2, 2, 2)))
                .toBeNull();
        });

        it('rejects non-3D points', () => {
            expect(() => circumscribeCircle3(Vector.fromArray([0, 0]),
                v3(1, 0, 0), v3(0, 1, 0))).toThrow();
        });
    });

    describe('circumscribeSphere3', () => {
        it('computes the circumsphere of the canonical tetrahedron', () => {
            // (0,0,0), (1,0,0), (0,1,0), (0,0,1): center (0.5,0.5,0.5),
            // radius sqrt(3)/2.
            const sphere = circumscribeSphere3(v3(0, 0, 0), v3(1, 0, 0),
                v3(0, 1, 0), v3(0, 0, 1));
            expect(sphere).not.toBeNull();
            for (let i = 0; i < 3; ++i) {
                expect(sphere!.center.values[i]).toBeCloseTo(0.5, 12);
            }
            expect(sphere!.radius).toBeCloseTo(Math.sqrt(3) / 2, 12);
        });

        it('returns null for coplanar points', () => {
            expect(circumscribeSphere3(v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0),
                v3(1, 1, 0))).toBeNull();
        });
    });

    describe('inscribeCircle3', () => {
        it('computes the incircle of a 3-4-5 triangle in the z = -1 plane',
            () => {
                const circle = inscribeCircle3(v3(0, 0, -1), v3(3, 0, -1),
                    v3(0, 4, -1));
                expect(circle).not.toBeNull();
                expect(circle!.radius).toBeCloseTo(1, 12);
                expect(circle!.center.values[0]).toBeCloseTo(1, 12);
                expect(circle!.center.values[1]).toBeCloseTo(1, 12);
                expect(circle!.center.values[2]).toBeCloseTo(-1, 12);
                expect(length(circle!.normal)).toBeCloseTo(1, 12);
            });

        it('returns null for collinear points', () => {
            expect(inscribeCircle3(v3(0, 0, 0), v3(1, 0, 0), v3(2, 0, 0)))
                .toBeNull();
        });
    });

    describe('inscribeSphere3', () => {
        it('computes the insphere of the canonical tetrahedron', () => {
            // Inradius = 3V/A with V = 1/6 and A = 3/2 + sqrt(3)/2, so
            // r = 1/(3 + sqrt(3)); the incenter has all coordinates equal to r.
            const sphere = inscribeSphere3(v3(0, 0, 0), v3(1, 0, 0),
                v3(0, 1, 0), v3(0, 0, 1));
            expect(sphere).not.toBeNull();
            const r = 1 / (3 + Math.sqrt(3));
            expect(sphere!.radius).toBeCloseTo(r, 12);
            for (let i = 0; i < 3; ++i) {
                expect(sphere!.center.values[i]).toBeCloseTo(r, 12);
            }
        });

        it('returns null for degenerate tetrahedra', () => {
            expect(inscribeSphere3(v3(0, 0, 0), v3(1, 0, 0), v3(2, 0, 0),
                v3(3, 0, 0))).toBeNull();
            expect(inscribeSphere3(v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0),
                v3(1, 1, 0))).toBeNull();
        });
    });

    it('circumscribed circles/spheres pass through their defining points and '
        + 'inscribed ones are tangent to every edge/face (randomized)', () => {
        let triangles = 0;
        let tetrahedra = 0;

        for (let trial = 0; trial < 300; ++trial) {
            const a = rv(), b = rv(), c = rv();
            const areaTwice = length(cross(sub(b, a), sub(c, a)));
            if (areaTwice < 5) {
                continue;
            }
            ++triangles;

            const circum = circumscribeCircle3(a, b, c);
            expect(circum).not.toBeNull();
            for (const p of [a, b, c]) {
                expect(length(sub(p, circum!.center)))
                    .toBeCloseTo(circum!.radius, 7);
            }
            // The circumcenter lies in the plane of the triangle.
            expect(planeDistance(circum!.center, a, b, c)).toBeCloseTo(0, 7);
            expect(length(circum!.normal)).toBeCloseTo(1, 10);

            const inc = inscribeCircle3(a, b, c);
            expect(inc).not.toBeNull();
            expect(lineDistance3(inc!.center, a, b)).toBeCloseTo(inc!.radius, 7);
            expect(lineDistance3(inc!.center, b, c)).toBeCloseTo(inc!.radius, 7);
            expect(lineDistance3(inc!.center, c, a)).toBeCloseTo(inc!.radius, 7);
            expect(planeDistance(inc!.center, a, b, c)).toBeCloseTo(0, 7);
            expect(length(inc!.normal)).toBeCloseTo(1, 10);
        }

        for (let trial = 0; trial < 300; ++trial) {
            const a = rv(), b = rv(), c = rv(), d = rv();
            const volTimes6 = Math.abs(dot(sub(b, a),
                cross(sub(c, a), sub(d, a))));
            if (volTimes6 < 50) {
                continue;
            }
            ++tetrahedra;

            const circum = circumscribeSphere3(a, b, c, d);
            expect(circum).not.toBeNull();
            for (const p of [a, b, c, d]) {
                expect(length(sub(p, circum!.center)))
                    .toBeCloseTo(circum!.radius, 6);
            }

            const ins = inscribeSphere3(a, b, c, d);
            expect(ins).not.toBeNull();
            const faces: Vector[][] = [[b, c, d], [a, c, d], [a, b, d],
                [a, b, c]];
            for (const f of faces) {
                expect(planeDistance(ins!.center, f[0], f[1], f[2]))
                    .toBeCloseTo(ins!.radius, 6);
            }
            // The insphere center is inside the tetrahedron: it is a convex
            // combination of the vertices with nonnegative weights, which is
            // equivalent to lying on the interior side of every face.
            const centroid = mul(0.25,
                add(add(a, b), add(c, d)));
            for (const f of faces) {
                const n = cross(sub(f[1], f[0]), sub(f[2], f[0]));
                const sCenter = dot(n, sub(ins!.center, f[0]));
                const sCentroid = dot(n, sub(centroid, f[0]));
                expect(Math.sign(sCenter)).toBe(Math.sign(sCentroid));
            }
        }

        expect(triangles).toBeGreaterThan(100);
        expect(tetrahedra).toBeGreaterThan(100);
    });
});
