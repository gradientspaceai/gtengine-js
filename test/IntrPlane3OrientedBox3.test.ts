import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane.js';
import {
    IntrPlane3OrientedBox3TI,
    defaultIntrPlane3OrientedBox3TIResult
} from '../src/IntrPlane3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, mul, normalize } from '../src/Vector.js';
import {
    check, fc, orientedBox as arbOrientedBox, plane as arbPlane,
    rotationFrame, wellScaledVector
} from './helpers/arbitraries.js';

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

// An oriented box rotated by 'angle' about the z-axis.
function box(center: number[], angle: number, extent: number[]): OrientedBox {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(Vector.fromArray(center),
        [Vector.fromArray([c, s, 0]), Vector.fromArray([-s, c, 0]),
            Vector.fromArray([0, 0, 1])],
        Vector.fromArray(extent));
}

// An independent test: the box vertices straddle (or touch) the plane.
function bruteForceIntersect(P: Hyperplane, B: OrientedBox): boolean {
    let numNeg = 0;
    let numPos = 0;
    for (const V of B.getVertices()) {
        const sd = dot(P.normal, V) - P.constant;
        if (sd < 0) {
            ++numNeg;
        }
        else if (sd > 0) {
            ++numPos;
        }
        else {
            return true;
        }
    }
    return numNeg > 0 && numPos > 0;
}

const ti = new IntrPlane3OrientedBox3TI();

describe('IntrPlane3OrientedBox3', () => {
    it('defaults to no intersection', () => {
        expect(defaultIntrPlane3OrientedBox3TIResult().intersect).toBe(false);
    });

    it('detects a plane cutting through an axis-aligned box', () => {
        const B = box([0, 0, 0], 0, [1, 1, 1]);
        expect(ti.test(plane([0, 0, 1], [0, 0, 0]), B).intersect).toBe(true);
        expect(ti.test(plane([0, 0, 1], [0, 0, 0.5]), B).intersect).toBe(true);
    });

    it('detects the tangent plane at a box face', () => {
        const B = box([0, 0, 0], 0, [1, 1, 1]);
        expect(ti.test(plane([0, 0, 1], [0, 0, 1]), B).intersect).toBe(true);
        expect(ti.test(plane([0, 0, 1], [0, 0, 1.0001]), B).intersect)
            .toBe(false);
    });

    it('accounts for the box orientation', () => {
        // A square rotated 45 degrees about z has half-width sqrt(2) in x.
        const B = box([0, 0, 0], Math.PI / 4, [1, 1, 1]);
        expect(ti.test(plane([1, 0, 0], [1.4, 0, 0]), B).intersect).toBe(true);
        expect(ti.test(plane([1, 0, 0], [1.5, 0, 0]), B).intersect).toBe(false);
        // The unrotated box would not reach x = 1.4.
        const A = box([0, 0, 0], 0, [1, 1, 1]);
        expect(ti.test(plane([1, 0, 0], [1.4, 0, 0]), A).intersect).toBe(false);
    });

    it('agrees with a vertex-sign test on random inputs', () => {
        let state = 13579;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numHits = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const P = plane([rand(), rand(), rand() + 0.001],
                [rand() * 2, rand() * 2, rand() * 2]);
            const B = box([rand() * 3, rand() * 3, rand() * 3],
                rand() * Math.PI,
                [0.2 + Math.abs(rand()), 0.2 + Math.abs(rand()),
                    0.2 + Math.abs(rand())]);
            const actual = ti.test(P, B).intersect;
            expect(actual).toBe(bruteForceIntersect(P, B));
            if (actual) {
                ++numHits;
            }
        }
        expect(numHits).toBeGreaterThan(50);
        expect(numHits).toBeLessThan(350);
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrPlane3OrientedBox3.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrPlane3OrientedBox3 verification', () => {
    const q = new IntrPlane3OrientedBox3TI();

    function corners(box: OrientedBox): Vector[] {
        const out: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            let p = box.center.clone();
            for (let d = 0; d < 3; ++d) {
                const s = ((i >> d) & 1) === 0 ? -1 : +1;
                p = add(p, mul(s * box.extent.values[d], box.axis[d]));
            }
            out.push(p);
        }
        return out;
    }

    it('agrees with the signed distances of the eight corners', () => {
        // A box is the convex hull of its corners, so the extreme values of
        // the linear function Dot(N,X) - c over the box are attained at
        // corners. This is an independent derivation of the radius formula.
        check(fc.tuple(arbPlane(3), arbOrientedBox(3)), ([P, B]) => {
            let lo = Number.POSITIVE_INFINITY;
            let hi = Number.NEGATIVE_INFINITY;
            for (const X of corners(B)) {
                const sd = dot(P.normal, X) - P.constant;
                lo = Math.min(lo, sd);
                hi = Math.max(hi, sd);
            }
            const got = q.test(P, B).intersect;
            const scale = 1 + Math.abs(lo) + Math.abs(hi);
            // Skip configurations that are within rounding of tangency; the
            // two computations of the extreme distance differ there in the
            // last bits.
            if (Math.abs(lo) < 1e-12 * scale || Math.abs(hi) < 1e-12 * scale) {
                return;
            }
            expect(got).toBe(lo <= 0 && hi >= 0);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(arbPlane(3), arbOrientedBox(3), rotationFrame(3),
            wellScaledVector(3)),
            ([P, B, R, t]) => {
                const rot = (v: Vector) => Vector.fromArray([
                    dot(R[0], v), dot(R[1], v), dot(R[2], v)]);
                const map = (v: Vector) => add(rot(v), t);
                const P2 = Hyperplane.fromNormalOrigin(rot(P.normal),
                    map(P.origin));
                const B2 = OrientedBox.fromCenterAxisExtent(map(B.center),
                    B.axis.map(rot), B.extent);
                const radius =
                    Math.abs(B.extent.values[0] * dot(P.normal, B.axis[0]))
                    + Math.abs(B.extent.values[1] * dot(P.normal, B.axis[1]))
                    + Math.abs(B.extent.values[2] * dot(P.normal, B.axis[2]));
                const d = Math.abs(dot(P.normal, B.center) - P.constant);
                if (Math.abs(d - radius) < 1e-9 * (1 + radius)) { return; }
                expect(q.test(P2, B2).intersect).toBe(q.test(P, B).intersect);
            });
    });

    it('reports a plane touching a single box corner', () => {
        const B = OrientedBox.fromCenterAxisExtent(Vector.zero(3),
            [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
                Vector.fromArray([0, 0, 1])],
            Vector.fromArray([1, 1, 1]));
        const n = Vector.fromArray([1, 1, 1]);
        normalize(n);
        const touching = Hyperplane.fromNormalOrigin(n,
            Vector.fromArray([1, 1, 1]));
        expect(q.test(touching, B).intersect).toBe(true);
        const outside = Hyperplane.fromNormalOrigin(n,
            Vector.fromArray([1.001, 1.001, 1.001]));
        expect(q.test(outside, B).intersect).toBe(false);
    });

    it('reports a degenerate (zero-extent) box as a point test', () => {
        check(fc.tuple(arbPlane(3), wellScaledVector(3), rotationFrame(3)),
            ([P, c, axis]) => {
                const B = OrientedBox.fromCenterAxisExtent(c, axis,
                    Vector.zero(3));
                const d = dot(P.normal, c) - P.constant;
                expect(q.test(P, B).intersect).toBe(d === 0);
            });
    });

    it('reports a box flattened into the plane as an intersection', () => {
        const B = OrientedBox.fromCenterAxisExtent(Vector.zero(3),
            [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
                Vector.fromArray([0, 0, 1])],
            Vector.fromArray([2, 3, 0]));
        const P = Hyperplane.fromNormalOrigin(Vector.fromArray([0, 0, 1]),
            Vector.zero(3));
        expect(q.test(P, B).intersect).toBe(true);
    });
});
