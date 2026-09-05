import { describe, it, expect } from 'vitest';
import { slerp, slerpUsingCosAngle, slerpUsingMidpoint } from '../src/Slerp.js';
import { chebyshevRatiosUsingCosAngle } from '../src/ChebyshevRatio.js';
import { check, fc, scaled, unitVector } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream Slerp.h.
// ---------------------------------------------------------------------------

describe('Slerp verification', () => {
    // Unit vectors whose components are well scaled, so no coordinate is
    // subnormal and the dot product carries full precision.
    const unitArray = (n: number): fc.Arbitrary<number[]> =>
        unitVector(n).map(v => [...v.values]);
    const quatArb = unitArray(4);
    // Pairs separated by an angle bounded away from 0 and pi, where the
    // 1/sin(A) in the weights is well conditioned.
    const pairArb = fc.tuple(quatArb, quatArb).filter(([a, b]) => {
        const c = dot(a, b);
        return c > -0.99 && c < 0.99;
    });
    const tArb = scaled(0, 1);

    it('reproduces the endpoints exactly, not just closely', () => {
        // f = {sin(A)/sin(A), 0} at t = 0 and {0, sin(A)/sin(A)} at t = 1, so
        // the same double is divided by itself and the weights are exactly
        // 1 and 0.
        check(pairArb, ([q0, q1]) => {
            const at0 = slerp(0, q0, q1);
            const at1 = slerp(1, q0, q1);
            for (let i = 0; i < 4; ++i) {
                if (at0[i] + 0 !== q0[i] + 0) { return false; }
                if (at1[i] + 0 !== q1[i] + 0) { return false; }
            }
            return true;
        });
    });

    it('stays on the unit hypersphere', () => {
        check(fc.tuple(pairArb, tArb), ([[q0, q1], t]) => {
            const r = slerp(t, q0, q1);
            return Math.abs(Math.sqrt(dot(r, r)) - 1) <= 1e-12;
        });
    });

    it('moves along the arc at constant angular speed', () => {
        // The defining property of slerp: the angle from q0 grows linearly.
        // The angle is read off with atan2 in the plane's own orthonormal
        // basis rather than with acos(dot(q0,r)); acos has an infinite
        // derivative at 1, so it turns the 1e-16 of the dot product into an
        // error of 1e-16/(t*A) radians and says nothing near t = 0.
        check(fc.tuple(pairArb, tArb), ([[q0, q1], t]) => {
            const c = dot(q0, q1);
            const angle = Math.atan2(Math.sqrt(Math.max(0, 1 - c * c)), c);
            const e1 = q1.map((v, i) => v - c * q0[i]);
            const u1 = e1.map(v => v / Math.sqrt(dot(e1, e1)));
            const r = slerp(t, q0, q1);
            const measured = Math.atan2(dot(r, u1), dot(r, q0));
            return Math.abs(measured - t * angle) <= 1e-12;
        });
    });

    it('stays in the plane spanned by the two inputs', () => {
        // Any component orthogonal to span{q0,q1} would mean the weights were
        // applied to the wrong vectors.
        check(fc.tuple(pairArb, tArb), ([[q0, q1], t]) => {
            const r = slerp(t, q0, q1);
            // Gram-Schmidt the plane, then project the residual out.
            const c = dot(q0, q1);
            const e1 = q1.map((v, i) => v - c * q0[i]);
            const n1 = Math.sqrt(dot(e1, e1));
            const u1 = e1.map(v => v / n1);
            const residual = r.map((v, i) => v - dot(r, q0) * q0[i]
                - dot(r, u1) * u1[i]);
            return Math.sqrt(dot(residual, residual)) <= 1e-12;
        });
    });

    it('is symmetric under reversing the arc', () => {
        // slerp(t, q0, q1) traverses the same arc as slerp(1-t, q1, q0).
        check(fc.tuple(pairArb, tArb), ([[q0, q1], t]) => {
            const a = slerp(t, q0, q1);
            const b = slerp(1 - t, q1, q0);
            for (let i = 0; i < 4; ++i) {
                if (Math.abs(a[i] - b[i]) > 1e-12) { return false; }
            }
            return true;
        });
    });

    it('is exactly the precomputed-cosine overload for the exact dot', () => {
        // The two overloads differ only in whether they compute the dot
        // product; feeding it back must reproduce the result bit for bit.
        check(fc.tuple(pairArb, tArb), ([[q0, q1], t]) => {
            const a = slerp(t, q0, q1);
            const b = slerpUsingCosAngle(t, q0, q1, dot(q0, q1));
            for (let i = 0; i < 4; ++i) {
                if (a[i] + 0 !== b[i] + 0) { return false; }
            }
            return true;
        });
    });

    it('agrees with the midpoint overload over the whole arc', () => {
        // The midpoint form halves the angle each side of t = 1/2, which is
        // what lets it handle angles up to pi; on the shared domain the two
        // must agree.
        check(fc.tuple(pairArb, tArb), ([[q0, q1], t]) => {
            const cosA = dot(q0, q1);
            const cosAH = Math.sqrt((1 + cosA) / 2);
            const qh = q0.map((v, i) => (v + q1[i]) / (2 * cosAH));
            const a = slerp(t, q0, q1);
            const b = slerpUsingMidpoint(t, q0, q1, qh, cosAH);
            for (let i = 0; i < 4; ++i) {
                if (Math.abs(a[i] - b[i]) > 1e-10) { return false; }
            }
            return true;
        });
    });

    it('reduces to the Chebyshev weights of the two inputs', () => {
        // The result is exactly f[0]*q0 + f[1]*q1 with f from ChebyshevRatio;
        // this pins the port's loop against the file it delegates to.
        check(fc.tuple(pairArb, tArb), ([[q0, q1], t]) => {
            const f = chebyshevRatiosUsingCosAngle(t, dot(q0, q1));
            const r = slerp(t, q0, q1);
            for (let i = 0; i < 4; ++i) {
                if (r[i] + 0 !== f[0] * q0[i] + f[1] * q1[i] + 0) {
                    return false;
                }
            }
            return true;
        });
    });

    it('works in every dimension from 2 up', () => {
        for (const n of [2, 3, 5, 7]) {
            check(fc.tuple(unitArray(n), unitArray(n), tArb),
                ([q0, q1, t]) => {
                    if (!(dot(q0, q1) > -0.99 && dot(q0, q1) < 0.99)) {
                        return true;
                    }
                    const r = slerp(t, q0, q1);
                    return r.length === n
                        && Math.abs(Math.sqrt(dot(r, r)) - 1) <= 1e-12;
                });
        }
    });

    it('rejects mismatched and degenerate dimensions', () => {
        expect(() => slerp(0.5, [1], [1])).toThrow('Invalid dimension.');
        expect(() => slerp(0.5, [1, 0], [1, 0, 0]))
            .toThrow('Mismatched dimensions.');
        expect(() => slerpUsingCosAngle(0.5, [1, 0], [1, 0, 0], 1))
            .toThrow('Mismatched dimensions.');
        expect(() => slerpUsingMidpoint(0.5, [1, 0], [0, 1], [1, 0, 0], 1))
            .toThrow('Mismatched dimensions.');
    });
});
