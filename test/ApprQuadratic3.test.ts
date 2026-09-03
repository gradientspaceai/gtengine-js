import { describe, expect, it } from 'vitest';
import { ApprQuadratic3, ApprQuadraticSphere3 } from '../src/ApprQuadratic3.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector } from '../src/Vector.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// V = (1, x, y, z, x^2, x*y, x*z, y^2, y*z, z^2)
function basis(x: number, y: number, z: number): number[] {
    return [1, x, y, z, x * x, x * y, x * z, y * y, y * z, z * z];
}

// Normalize to unit length with a deterministic sign so eigenvectors can be
// compared.
function canonical(c: readonly number[]): number[] {
    let len = 0, imax = 0;
    for (let i = 0; i < c.length; ++i) {
        len += c[i] * c[i];
        if (Math.abs(c[i]) > Math.abs(c[imax])) {
            imax = i;
        }
    }
    len = Math.sqrt(len);
    const s = (c[imax] < 0 ? -1 : 1) / len;
    return c.map((v) => v * s);
}

// Points spread over an ellipsoid with semi-axes (a,b,c).
function ellipsoidPoints(a: number, b: number, c: number,
    numTheta: number, numPhi: number): Vector[] {
    const points: Vector[] = [];
    for (let i = 1; i < numPhi; ++i) {
        const phi = (Math.PI * i) / numPhi;
        for (let j = 0; j < numTheta; ++j) {
            const theta = (2 * Math.PI * j) / numTheta;
            points.push(v3(a * Math.sin(phi) * Math.cos(theta),
                b * Math.sin(phi) * Math.sin(theta),
                c * Math.cos(phi)));
        }
    }
    points.push(v3(0, 0, c), v3(0, 0, -c));
    return points;
}

describe('ApprQuadratic3', () => {
    it('recovers an ellipsoid from points sampled on it', () => {
        // x^2/4 + y^2/9 + z^2 = 1, i.e. -1 + x^2/4 + y^2/9 + z^2 = 0.
        const points = ellipsoidPoints(2, 3, 1, 8, 6);

        const result = new ApprQuadratic3().compute(points);
        expect(result.coefficients.length).toBe(10);
        expect(result.minEigenvalue).toBeLessThan(1e-14);

        let sqrLen = 0;
        for (const c of result.coefficients) {
            sqrLen += c * c;
        }
        expect(sqrLen).toBeCloseTo(1, 10);

        const expected = canonical([-1, 0, 0, 0, 1 / 4, 0, 0, 1 / 9, 0, 1]);
        const actual = canonical(result.coefficients);
        for (let i = 0; i < 10; ++i) {
            expect(actual[i]).toBeCloseTo(expected[i], 6);
        }

        for (const p of points) {
            const b = basis(p.get(0), p.get(1), p.get(2));
            let sum = 0;
            for (let i = 0; i < 10; ++i) {
                sum += result.coefficients[i] * b[i];
            }
            expect(Math.abs(sum)).toBeLessThan(1e-8);
        }
    });

    it('recovers a sheared and translated quadric', () => {
        // A unit sphere mapped by (u,v,w) -> (u + 0.3*v, v + 0.2*w, w) and
        // translated; still an exact quadric surface.
        const base = ellipsoidPoints(1, 1, 1, 10, 7);
        const points = base.map((p) => {
            const u = p.get(0), v = p.get(1), w = p.get(2);
            return v3(2 + u + 0.3 * v, -1 + v + 0.2 * w, 5 + w);
        });

        const result = new ApprQuadratic3().compute(points);
        expect(result.minEigenvalue).toBeLessThan(1e-10);
        for (const p of points) {
            const b = basis(p.get(0), p.get(1), p.get(2));
            let sum = 0;
            for (let i = 0; i < 10; ++i) {
                sum += result.coefficients[i] * b[i];
            }
            expect(Math.abs(sum)).toBeLessThan(1e-6);
        }
    });

    it('recovers a circular cylinder', () => {
        // x^2 + y^2 = 4, sampled at several heights.
        const points: Vector[] = [];
        for (let k = 0; k < 4; ++k) {
            const z = -3 + 2 * k;
            for (let j = 0; j < 8; ++j) {
                const t = (2 * Math.PI * j) / 8;
                points.push(v3(2 * Math.cos(t), 2 * Math.sin(t), z));
            }
        }

        const result = new ApprQuadratic3().compute(points);
        expect(result.minEigenvalue).toBeLessThan(1e-14);
        const expected = canonical([-4, 0, 0, 0, 1, 0, 0, 1, 0, 0]);
        const actual = canonical(result.coefficients);
        for (let i = 0; i < 10; ++i) {
            expect(actual[i]).toBeCloseTo(expected[i], 8);
        }
    });

    it('returns a positive fit measure for points not on a quadric', () => {
        const points: Vector[] = [];
        for (let i = 1; i < 8; ++i) {
            const phi = (Math.PI * i) / 8;
            for (let j = 0; j < 12; ++j) {
                const theta = (2 * Math.PI * j) / 12;
                // A radius that varies with 3*theta: not a quadric.
                const r = 1 + 0.4 * Math.cos(3 * theta) * Math.sin(phi);
                points.push(v3(r * Math.sin(phi) * Math.cos(theta),
                    r * Math.sin(phi) * Math.sin(theta),
                    r * Math.cos(phi)));
            }
        }

        const result = new ApprQuadratic3().compute(points);
        expect(result.minEigenvalue).toBeGreaterThan(1e-8);
    });
});

describe('ApprQuadraticSphere3', () => {
    it('recovers a sphere from points sampled on it', () => {
        const cx = 1, cy = 2, cz = -3, radius = 4;
        const points = ellipsoidPoints(radius, radius, radius, 8, 6).map(
            (p) => v3(cx + p.get(0), cy + p.get(1), cz + p.get(2)));

        const sphere = new Hypersphere(3);
        const measure = new ApprQuadraticSphere3().compute(points, sphere);
        expect(measure).toBeLessThan(1e-13);
        expect(sphere.center.get(0)).toBeCloseTo(cx, 8);
        expect(sphere.center.get(1)).toBeCloseTo(cy, 8);
        expect(sphere.center.get(2)).toBeCloseTo(cz, 8);
        expect(sphere.radius).toBeCloseTo(radius, 8);
    });

    it('recovers the unit sphere at the origin', () => {
        const points = ellipsoidPoints(1, 1, 1, 10, 8);
        const sphere = new Hypersphere(3);
        const measure = new ApprQuadraticSphere3().compute(points, sphere);
        expect(measure).toBeLessThan(1e-15);
        expect(sphere.center.get(0)).toBeCloseTo(0, 10);
        expect(sphere.center.get(1)).toBeCloseTo(0, 10);
        expect(sphere.center.get(2)).toBeCloseTo(0, 10);
        expect(sphere.radius).toBeCloseTo(1, 10);
    });

    it('averages symmetric radial perturbations of a sphere', () => {
        // Every direction contributes a sample at radius r-d and one at
        // r+d. The algebraic fit is constrained by Length(C') = 1, a norm
        // that is not invariant under translation, so the fitted sphere is
        // only near the true one; it must stay within the noise band.
        const cx = -2, cy = 0.5, cz = 3, radius = 2, d = 0.05;
        const base = ellipsoidPoints(1, 1, 1, 12, 8);
        const points: Vector[] = [];
        for (const p of base) {
            for (const r of [radius - d, radius + d]) {
                points.push(v3(cx + r * p.get(0), cy + r * p.get(1),
                    cz + r * p.get(2)));
            }
        }

        const sphere = new Hypersphere(3);
        const measure = new ApprQuadraticSphere3().compute(points, sphere);
        expect(measure).toBeGreaterThan(0);
        expect(Math.abs(sphere.center.get(0) - cx)).toBeLessThan(d);
        expect(Math.abs(sphere.center.get(1) - cy)).toBeLessThan(d);
        expect(Math.abs(sphere.center.get(2) - cz)).toBeLessThan(d);
        expect(sphere.radius).toBeGreaterThan(radius - d);
        expect(sphere.radius).toBeLessThan(radius + d);
    });
});
