import { describe, expect, it } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3.js';
import { DistPoint3Cylinder3 } from '../src/DistPoint3Cylinder3.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function cylinder(origin: number[], direction: number[], radius: number,
    height: number): Cylinder3 {
    const d = v(...direction);
    normalize(d);
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(v(...origin), d), radius, height);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistPoint3Cylinder3', () => {
    const query = new DistPoint3Cylinder3();

    it('reports zero distance for a point inside a finite cylinder', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 2, 4);
        const result = query.compute(v(0.5, 0.5, 1), c);
        expect(result.distance).toBe(0);
        expect(result.closest[1].values[0]).toBeCloseTo(0.5, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(0.5, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
    });

    it('reports zero distance for a point on the cylinder axis', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 2, 4);
        const result = query.compute(v(0, 0, 1), c);
        expect(result.distance).toBe(0);
    });

    it('measures a point outside the wall of an infinite cylinder', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 1);
        c.makeInfiniteCylinder();
        const result = query.compute(v(5, 0, 100), c);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(100, 12);
    });

    it('reports zero distance inside an infinite cylinder at any height',
        () => {
            const c = cylinder([0, 0, 0], [0, 0, 1], 1, 1);
            c.makeInfiniteCylinder();
            const result = query.compute(v(0.5, 0, -1000), c);
            expect(result.distance).toBe(0);
        });

    it('clamps to the cap of a finite cylinder', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 2, 4);
        // Directly above the axis, beyond the +z cap at z = 2.
        const result = query.compute(v(0, 0, 6), c);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(2, 12);
    });

    it('measures the rim of a finite cylinder', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 2);
        // Radially 3 units outside the wall and 4 units above the +z cap.
        const result = query.compute(v(4, 0, 5), c);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
    });

    it('measures the -z cap of a finite cylinder', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 2, 4);
        const result = query.compute(v(0.5, 0.5, -7), c);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(-2, 12);
    });

    it('rejects a nonpositive radius', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 0, 2);
        expect(() => query.compute(v(1, 1, 1), c)).toThrow(
            /positive radius/);
    });

    it('rejects a zero height for a finite cylinder', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 0);
        expect(() => query.compute(v(1, 1, 1), c)).toThrow(
            /positive height/);
    });

    it('agrees with a dense sampling of a tilted finite cylinder', () => {
        const rnd = makeRandom(2718);
        const origin = v(0.5, -1, 0.25);
        const dir = v(1, 2, 3);
        normalize(dir);
        const radius = 1.25;
        const height = 3;
        const c = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(origin, dir), radius, height);

        // An orthonormal basis of the plane perpendicular to the axis.
        const U = v(-dir.values[1], dir.values[0], 0);
        normalize(U);
        const W = cross(dir, U);

        for (let trial = 0; trial < 30; ++trial) {
            const p = v(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4);
            const result = query.compute(p, c);

            // The reported closest point is in the solid cylinder.
            const delta = sub(result.closest[1], origin);
            const h = dot(delta, dir);
            const radial = sub(delta, mul(h, dir));
            expect(Math.abs(h)).toBeLessThanOrEqual(0.5 * height + 1e-9);
            expect(Math.sqrt(dot(radial, radial))).toBeLessThanOrEqual(
                radius + 1e-9);

            // The reported closest point realizes the reported distance.
            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 9);

            // No sampled solid-cylinder point is closer.
            const nR = 10, nA = 40, nH = 12;
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= nR; ++i) {
                const r = radius * i / nR;
                for (let j = 0; j < nA; ++j) {
                    const a = 2 * Math.PI * j / nA;
                    for (let k = 0; k <= nH; ++k) {
                        const hh = height * (k / nH - 0.5);
                        const q = add(origin, add(mul(hh, dir),
                            add(mul(r * Math.cos(a), U),
                                mul(r * Math.sin(a), W))));
                        const f = sub(p, q);
                        best = Math.min(best, dot(f, f));
                    }
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});
