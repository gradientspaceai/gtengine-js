import { describe, it, expect } from 'vitest';
import { slerp, slerpUsingCosAngle, slerpUsingMidpoint } from '../src/Slerp';

// Closed-form geodesic reference:
//   slerp(t,q0,q1) = [sin((1-t)*A)*q0 + sin(t*A)*q1]/sin(A), cos(A)=dot(q0,q1)
function referenceSlerp(t: number, q0: number[], q1: number[]): number[] {
    const cosA = dot(q0, q1);
    const angle = Math.acos(cosA);
    const sinA = Math.sin(angle);
    if (sinA === 0) {
        return q0.map((v, i) => (1 - t) * v + t * q1[i]);
    }
    const f0 = Math.sin((1 - t) * angle) / sinA;
    const f1 = Math.sin(t * angle) / sinA;
    return q0.map((v, i) => f0 * v + f1 * q1[i]);
}

function dot(a: number[], b: number[]): number {
    let s = 0;
    for (let i = 0; i < a.length; ++i) {
        s += a[i] * b[i];
    }
    return s;
}

function normalize(v: number[]): number[] {
    const len = Math.sqrt(dot(v, v));
    return v.map((x) => x / len);
}

function expectClose(actual: number[], expected: number[], tol: number): void {
    expect(actual.length).toBe(expected.length);
    for (let i = 0; i < expected.length; ++i) {
        expect(Math.abs(actual[i] - expected[i])).toBeLessThanOrEqual(tol);
    }
}

// Deterministic pseudo-random generator so failures are reproducible.
function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (s * 1664525 + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

// A random unit 4-vector.
function randomQuat(rand: () => number): number[] {
    return normalize([
        2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1
    ]);
}

describe('slerp', () => {
    it('traverses a quarter circle at constant speed', () => {
        // q0 and q1 are orthogonal, so A = pi/2 and the arc is the unit
        // circle in the (x,y) plane traversed with angle t*pi/2.
        const q0 = [1, 0, 0, 0];
        const q1 = [0, 1, 0, 0];
        for (let i = 0; i <= 16; ++i) {
            const t = i / 16;
            const angle = t * Math.PI / 2;
            expectClose(slerp(t, q0, q1),
                [Math.cos(angle), Math.sin(angle), 0, 0], 1e-14);
        }
    });

    it('reproduces the endpoints', () => {
        const q0 = normalize([1, 2, 3, 4]);
        const q1 = normalize([4, -1, 2, 1]);
        expectClose(slerp(0, q0, q1), q0, 1e-14);
        expectClose(slerp(1, q0, q1), q1, 1e-14);
    });

    it('produces the arc midpoint (q0+q1)/|q0+q1| at t = 1/2', () => {
        const q0 = normalize([1, 0.2, -0.3, 0.5]);
        const q1 = normalize([0.4, 1, 0.1, -0.2]);
        const sum = q0.map((v, i) => v + q1[i]);
        expectClose(slerp(0.5, q0, q1), normalize(sum), 1e-14);
    });

    it('produces unit-length results for unit-length inputs', () => {
        const rand = makeRandom(12345);
        for (let k = 0; k < 200; ++k) {
            const q0 = randomQuat(rand);
            let q1 = randomQuat(rand);
            if (dot(q0, q1) < 0) {
                q1 = q1.map((v) => -v);
            }
            if (dot(q0, q1) >= 1) {
                continue;
            }
            for (const t of [0.1, 0.37, 0.5, 0.82]) {
                const r = slerp(t, q0, q1);
                expect(Math.abs(Math.sqrt(dot(r, r)) - 1)).toBeLessThan(1e-12);
            }
        }
    });

    it('matches the closed-form geodesic for random pairs', () => {
        const rand = makeRandom(987);
        for (let k = 0; k < 300; ++k) {
            const q0 = randomQuat(rand);
            let q1 = randomQuat(rand);
            if (dot(q0, q1) < 0) {
                q1 = q1.map((v) => -v);
            }
            if (dot(q0, q1) >= 1) {
                continue;
            }
            for (const t of [0, 0.13, 0.5, 0.77, 1]) {
                expectClose(slerp(t, q0, q1), referenceSlerp(t, q0, q1), 1e-12);
            }
        }
    });

    it('travels the arc with constant angular speed', () => {
        const q0 = normalize([1, 0.3, 0.1, -0.2]);
        const q1 = normalize([-0.2, 0.9, 0.4, 0.3]);
        const total = Math.acos(dot(q0, q1));
        for (let i = 0; i <= 10; ++i) {
            const t = i / 10;
            const r = slerp(t, q0, q1);
            const angle = Math.acos(Math.min(1, Math.max(-1, dot(q0, r))));
            expect(Math.abs(angle - t * total)).toBeLessThan(1e-12);
        }
    });

    it('degenerates to linear interpolation for identical inputs', () => {
        // cosA = 1 exactly, so the ratios are the removable-singularity
        // limits {1-t, t}.
        const q0 = [0, 0, 0, 1];
        for (const t of [0, 0.25, 0.5, 1]) {
            expectClose(slerp(t, q0, q0), q0, 0);
        }
    });

    it('is continuous in the nearly-parallel limit', () => {
        // As the angle shrinks the result must approach the lerp of the
        // endpoints (renormalized), with no blowup from sin(A) -> 0.
        const q0 = [1, 0, 0, 0];
        for (const eps of [1e-4, 1e-6, 1e-8, 1e-10]) {
            const q1 = normalize([Math.cos(eps), Math.sin(eps), 0, 0]);
            const r = slerp(0.5, q0, q1);
            expectClose(r, normalize([1 + q1[0], q1[1], 0, 0]), 1e-12);
            expect(Math.abs(Math.sqrt(dot(r, r)) - 1)).toBeLessThan(1e-12);
        }
    });

    it('extrapolates outside [0,1] along the same great arc', () => {
        const q0 = [1, 0, 0, 0];
        const q1 = [0, 1, 0, 0];
        for (const t of [-0.5, 1.5]) {
            const angle = t * Math.PI / 2;
            expectClose(slerp(t, q0, q1),
                [Math.cos(angle), Math.sin(angle), 0, 0], 1e-14);
        }
    });

    it('works for dimensions other than 4', () => {
        const q0 = [1, 0];
        const q1 = [0, 1];
        expectClose(slerp(0.5, q0, q1),
            [Math.SQRT1_2, Math.SQRT1_2], 1e-15);
        const p0 = normalize([1, 1, 0]);
        const p1 = normalize([0, 1, 1]);
        expectClose(slerp(0.25, p0, p1), referenceSlerp(0.25, p0, p1), 1e-14);
    });

    it('rejects bad dimensions and the angle pi', () => {
        expect(() => slerp(0.5, [1], [1])).toThrow(/Invalid dimension/);
        expect(() => slerp(0.5, [1, 0], [1, 0, 0])).toThrow(/Mismatched dimensions/);
        expect(() => slerp(0.5, [1, 0], [-1, 0])).toThrow(/Invalid angle/);
    });
});

describe('slerpUsingCosAngle', () => {
    it('matches slerp when given the exact dot product', () => {
        const rand = makeRandom(555);
        for (let k = 0; k < 100; ++k) {
            const q0 = randomQuat(rand);
            let q1 = randomQuat(rand);
            if (dot(q0, q1) < 0) {
                q1 = q1.map((v) => -v);
            }
            const cosA = dot(q0, q1);
            if (cosA >= 1) {
                continue;
            }
            for (const t of [0, 0.3, 0.5, 0.9, 1]) {
                expect(slerpUsingCosAngle(t, q0, q1, cosA))
                    .toEqual(slerp(t, q0, q1));
            }
        }
    });

    it('matches the closed-form geodesic', () => {
        const q0 = [1, 0, 0, 0];
        const q1 = normalize([1, 1, 0, 0]);
        const cosA = dot(q0, q1);
        for (let i = 0; i <= 8; ++i) {
            const t = i / 8;
            expectClose(slerpUsingCosAngle(t, q0, q1, cosA),
                referenceSlerp(t, q0, q1), 1e-14);
        }
    });

    it('validates its inputs', () => {
        expect(() => slerpUsingCosAngle(0.5, [1], [1], 0.5))
            .toThrow(/Invalid dimension/);
        expect(() => slerpUsingCosAngle(0.5, [1, 0], [0, 1], -1))
            .toThrow(/Invalid angle/);
    });
});

describe('slerpUsingMidpoint', () => {
    function midpointData(q0: number[], q1: number[]): { qh: number[]; cosAH: number } {
        const cosA = dot(q0, q1);
        const cosAH = Math.sqrt((1 + cosA) / 2);
        const qh = q0.map((v, i) => (v + q1[i]) / (2 * cosAH));
        return { qh, cosAH };
    }

    it('agrees with slerp for angles in [0,pi/2]', () => {
        const q0 = [1, 0, 0, 0];
        const q1 = [0, 1, 0, 0];
        const { qh, cosAH } = midpointData(q0, q1);
        for (let i = 0; i <= 20; ++i) {
            const t = i / 20;
            expectClose(slerpUsingMidpoint(t, q0, q1, qh, cosAH),
                referenceSlerp(t, q0, q1), 1e-13);
        }
    });

    it('handles angles approaching pi, where the plain form is ill-conditioned', () => {
        // A = 3*pi/4: the half-angle form uses cosAH = cos(3*pi/8) > 0.
        const angle = 3 * Math.PI / 4;
        const q0 = [1, 0, 0, 0];
        const q1 = [Math.cos(angle), Math.sin(angle), 0, 0];
        const { qh, cosAH } = midpointData(q0, q1);
        for (let i = 0; i <= 20; ++i) {
            const t = i / 20;
            const r = slerpUsingMidpoint(t, q0, q1, qh, cosAH);
            const expected = [Math.cos(t * angle), Math.sin(t * angle), 0, 0];
            expectClose(r, expected, 1e-13);
        }
    });

    it('reproduces the endpoints and the midpoint', () => {
        const q0 = normalize([1, 0.4, -0.2, 0.7]);
        const q1 = normalize([-0.3, 1, 0.5, 0.1]);
        const { qh, cosAH } = midpointData(q0, q1);
        expectClose(slerpUsingMidpoint(0, q0, q1, qh, cosAH), q0, 1e-13);
        expectClose(slerpUsingMidpoint(0.5, q0, q1, qh, cosAH), qh, 1e-13);
        expectClose(slerpUsingMidpoint(1, q0, q1, qh, cosAH), q1, 1e-13);
    });

    it('is continuous across the t = 1/2 branch switch', () => {
        const q0 = normalize([1, 0.4, -0.2, 0.7]);
        const q1 = normalize([-0.3, 1, 0.5, 0.1]);
        const { qh, cosAH } = midpointData(q0, q1);
        const below = slerpUsingMidpoint(0.5 - 1e-9, q0, q1, qh, cosAH);
        const above = slerpUsingMidpoint(0.5 + 1e-9, q0, q1, qh, cosAH);
        expectClose(below, above, 1e-8);
    });

    it('matches the closed-form geodesic for random pairs', () => {
        const rand = makeRandom(4242);
        for (let k = 0; k < 200; ++k) {
            const q0 = randomQuat(rand);
            let q1 = randomQuat(rand);
            if (dot(q0, q1) < 0) {
                q1 = q1.map((v) => -v);
            }
            if (dot(q0, q1) >= 1) {
                continue;
            }
            const { qh, cosAH } = midpointData(q0, q1);
            for (const t of [0.2, 0.5, 0.6, 0.95]) {
                expectClose(slerpUsingMidpoint(t, q0, q1, qh, cosAH),
                    referenceSlerp(t, q0, q1), 1e-12);
            }
        }
    });

    it('validates its inputs', () => {
        expect(() => slerpUsingMidpoint(0.5, [1], [1], [1], 0.9))
            .toThrow(/Invalid dimension/);
        expect(() => slerpUsingMidpoint(0.5, [1, 0], [0, 1], [1, 0, 0], 0.9))
            .toThrow(/Mismatched dimensions/);
    });
});
