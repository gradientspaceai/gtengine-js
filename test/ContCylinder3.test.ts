import { describe, it, expect } from 'vitest';
import {
    getContainerCylinder3,
    inContainerCylinder3
} from '../src/ContCylinder3.js';
import { Cylinder3 } from '../src/Cylinder3.js';
import { Line } from '../src/Line.js';
import { Vector, dot, length, mul, sub } from '../src/Vector.js';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// Independent containment check with a tolerance, since points on the
// cylinder boundary satisfy the upstream strict tests only up to rounding.
function containedWithin(point: Vector, cylinder: Cylinder3,
    tol: number): boolean {
    const diff = sub(point, cylinder.axis.origin);
    const zProj = dot(diff, cylinder.axis.direction);
    if (Math.abs(zProj) * 2 > cylinder.height + tol) {
        return false;
    }
    const xyProj = sub(diff, mul(zProj, cylinder.axis.direction));
    return length(xyProj) <= cylinder.radius + tol;
}

function makeRandom(seed: number): () => number {
    let s = seed;
    return () => {
        s = (s * 1103515245 + 12345) % 2147483648;
        return s / 2147483648;
    };
}

describe('getContainerCylinder3', () => {
    it('recovers a cylinder from points on its lateral surface', () => {
        const points: Vector[] = [];
        for (const z of [-2, 2]) {
            for (let k = 0; k < 16; ++k) {
                const a = (2 * Math.PI * k) / 16;
                points.push(v(2 * Math.cos(a), 2 * Math.sin(a), z));
            }
        }

        const cylinder = getContainerCylinder3(points);
        // The dominant variance is along z, so the fitted axis is +/- e2.
        expect(Math.abs(cylinder.axis.direction.values[2])).toBeCloseTo(1, 10);
        expect(cylinder.axis.direction.values[0]).toBeCloseTo(0, 10);
        expect(cylinder.axis.direction.values[1]).toBeCloseTo(0, 10);
        expect(cylinder.axis.origin.values[0]).toBeCloseTo(0, 10);
        expect(cylinder.axis.origin.values[1]).toBeCloseTo(0, 10);
        expect(cylinder.axis.origin.values[2]).toBeCloseTo(0, 10);
        expect(cylinder.radius).toBeCloseTo(2, 10);
        expect(cylinder.height).toBeCloseTo(4, 10);
    });

    it('contains every input point (random clouds)', () => {
        const rand = makeRandom(9871);
        for (let trial = 0; trial < 20; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 40; ++i) {
                // An elongated cloud so the least-squares axis is stable.
                points.push(v(
                    10 * (rand() - 0.5),
                    2 * (rand() - 0.5),
                    2 * (rand() - 0.5)));
            }
            const cylinder = getContainerCylinder3(points);
            expect(cylinder.radius).toBeGreaterThan(0);
            expect(cylinder.height).toBeGreaterThan(0);
            for (const p of points) {
                expect(containedWithin(p, cylinder, 1e-9)).toBe(true);
            }
        }
    });

    it('handles collinear points (zero radius)', () => {
        const points = [v(0, 0, 0), v(1, 1, 1), v(2, 2, 2), v(-1, -1, -1)];
        const cylinder = getContainerCylinder3(points);
        expect(cylinder.radius).toBeCloseTo(0, 10);
        expect(cylinder.height).toBeCloseTo(3 * Math.sqrt(3), 10);
        expect(cylinder.axis.origin.values[0]).toBeCloseTo(0.5, 10);
        expect(cylinder.axis.origin.values[1]).toBeCloseTo(0.5, 10);
        expect(cylinder.axis.origin.values[2]).toBeCloseTo(0.5, 10);
    });

    it('handles coincident points (degenerate cylinder)', () => {
        const points = [v(3, -1, 2), v(3, -1, 2), v(3, -1, 2)];
        const cylinder = getContainerCylinder3(points);
        expect(cylinder.radius).toBeCloseTo(0, 12);
        expect(cylinder.height).toBeCloseTo(0, 12);
        expect(cylinder.axis.origin.values[0]).toBeCloseTo(3, 12);
        expect(cylinder.axis.origin.values[1]).toBeCloseTo(-1, 12);
        expect(cylinder.axis.origin.values[2]).toBeCloseTo(2, 12);
    });

    it('throws for an empty point set', () => {
        expect(() => getContainerCylinder3([])).toThrow();
    });

    it('throws for non-3D points', () => {
        expect(() => getContainerCylinder3([Vector.fromArray([1, 2])]))
            .toThrow();
    });
});

describe('inContainerCylinder3', () => {
    const cylinder = Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(v(0, 0, 0), v(0, 0, 1)), 2, 4);

    it('accepts interior points', () => {
        expect(inContainerCylinder3(v(0, 0, 0), cylinder)).toBe(true);
        expect(inContainerCylinder3(v(1, 1, 1), cylinder)).toBe(true);
    });

    it('accepts boundary points (the boundary is part of the cylinder)', () => {
        expect(inContainerCylinder3(v(2, 0, 0), cylinder)).toBe(true);
        expect(inContainerCylinder3(v(0, 0, 2), cylinder)).toBe(true);
        expect(inContainerCylinder3(v(0, 2, -2), cylinder)).toBe(true);
    });

    it('rejects points outside the radius or the height', () => {
        expect(inContainerCylinder3(v(2.0001, 0, 0), cylinder)).toBe(false);
        expect(inContainerCylinder3(v(0, 0, 2.0001), cylinder)).toBe(false);
        expect(inContainerCylinder3(v(0, 0, -2.0001), cylinder)).toBe(false);
    });

    it('works for an oblique axis', () => {
        const d = v(1 / Math.sqrt(3), 1 / Math.sqrt(3), 1 / Math.sqrt(3));
        const oblique = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(v(1, 2, 3), d), 1, 2);
        // A point one unit along the axis from the origin, on the axis.
        expect(inContainerCylinder3(v(1 + d.values[0], 2 + d.values[1],
            3 + d.values[2]), oblique)).toBe(true);
        // Two units along the axis is outside the half-height of 1.
        expect(inContainerCylinder3(v(1 + 2 * d.values[0], 2 + 2 * d.values[1],
            3 + 2 * d.values[2]), oblique)).toBe(false);
    });

    it('throws for non-3D points', () => {
        expect(() => inContainerCylinder3(Vector.fromArray([1, 2]), cylinder))
            .toThrow();
    });
});
