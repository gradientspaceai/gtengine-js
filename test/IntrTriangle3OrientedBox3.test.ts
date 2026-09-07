import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import {
    IntrTriangle3AlignedBox3TI,
    IntrTriangle3AlignedBox3FI
} from '../src/IntrTriangle3AlignedBox3.js';
import {
    IntrTriangle3OrientedBox3TI,
    IntrTriangle3OrientedBox3FI,
    defaultIntrTriangle3OrientedBox3FIResult,
    defaultIntrTriangle3OrientedBox3TIResult
} from '../src/IntrTriangle3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Triangle } from '../src/Triangle.js';
import { DistTriangle3OrientedBox3 } from '../src/DistTriangle3OrientedBox3.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc,
    rotationFrame as randomRotationFrame, wellScaledVector
} from './helpers/arbitraries.js';
import { cross } from '../src/Vector3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function tri(v0: number[], v1: number[], v2: number[]): Triangle {
    return Triangle.fromVertices(Vector.fromArray(v0), Vector.fromArray(v1),
        Vector.fromArray(v2));
}

// A right-handed orthonormal frame from a rotation about the given axis.
function rotationFrame(axis: Vector, angle: number): Vector[] {
    const u = axis.clone();
    normalize(u);
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return [Vector.unit(3, 0), Vector.unit(3, 1), Vector.unit(3, 2)].map(
        (e) => add(add(mul(c, e), mul(s, cross(u, e))),
            mul((1 - c) * dot(u, e), u)));
}

// True when the point is in the oriented box (with tolerance).
function inOBB(B: OrientedBox, X: Vector, tol: number): boolean {
    const diff = sub(X, B.center);
    for (let i = 0; i < 3; ++i) {
        if (Math.abs(dot(B.axis[i], diff)) > B.extent.values[i] + tol) {
            return false;
        }
    }
    return true;
}

function planeDistance(T: Triangle, X: Vector): number {
    const n = cross(sub(T.v[1], T.v[0]), sub(T.v[2], T.v[0]));
    return dot(n, sub(X, T.v[0])) / Math.sqrt(dot(n, n));
}

function inTriangle(T: Triangle, X: Vector, tol: number): boolean {
    const n = cross(sub(T.v[1], T.v[0]), sub(T.v[2], T.v[0]));
    const area2 = dot(n, n);
    for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
        const e = sub(T.v[i1], T.v[i0]);
        const d = sub(X, T.v[i0]);
        if (dot(cross(e, d), n) < -tol * area2) {
            return false;
        }
    }
    return true;
}

const ti = new IntrTriangle3OrientedBox3TI();
const fi = new IntrTriangle3OrientedBox3FI();
const abTI = new IntrTriangle3AlignedBox3TI();
const abFI = new IntrTriangle3AlignedBox3FI();

const identityAxes = [Vector.unit(3, 0), Vector.unit(3, 1), Vector.unit(3, 2)];
const unitOBB = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0), identityAxes,
    vec(1, 1, 1));

describe('IntrTriangle3OrientedBox3', () => {
    it('default-constructs the results as empty', () => {
        expect(defaultIntrTriangle3OrientedBox3TIResult())
            .toEqual({ intersect: false });
        const r = defaultIntrTriangle3OrientedBox3FIResult();
        expect(r.insidePolygon).toEqual([]);
        expect(r.outsidePolygons).toEqual([]);
    });

    it('keeps a triangle fully inside the box', () => {
        const T = tri([-0.5, -0.5, 0], [0.5, -0.5, 0], [0, 0.5, 0]);
        expect(ti.test(T, unitOBB).intersect).toBe(true);
        const result = fi.find(T, unitOBB);
        expect(result.insidePolygon.length).toBe(3);
        expect(result.outsidePolygons.length).toBe(0);
    });

    it('rejects a triangle entirely outside the box', () => {
        const T = tri([5, 5, 5], [6, 5, 5], [5, 6, 5]);
        expect(ti.test(T, unitOBB).intersect).toBe(false);
        expect(fi.find(T, unitOBB).insidePolygon.length).toBe(0);
    });

    it('clips a large triangle to a rotated box cross section', () => {
        // The box is rotated 45 degrees about z, so the slice at z = 0 is the
        // square of "radius" 1 rotated by 45 degrees, that is, the diamond
        // with vertices at distance sqrt(2) along the axes.
        const box = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0),
            rotationFrame(vec(0, 0, 1), Math.PI / 4), vec(1, 1, 1));
        const T = tri([-10, -10, 0], [10, -10, 0], [0, 10, 0]);
        const result = fi.find(T, box);
        expect(result.insidePolygon.length).toBe(4);
        for (const p of result.insidePolygon) {
            expect(inOBB(box, p, 1e-10)).toBe(true);
            expect(p.values[2]).toBeCloseTo(0, 12);
            expect(Math.abs(p.values[0]) + Math.abs(p.values[1]))
                .toBeCloseTo(Math.SQRT2, 10);
        }
    });

    it('matches the aligned-box query when the axes are the standard basis', () => {
        const B = AlignedBox.fromMinMax(vec(-1, -1, -1), vec(1, 1, 1));
        let seed = 33445566;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (a: number, b: number): number => a + (b - a) * rand();

        for (let trial = 0; trial < 600; ++trial) {
            const T = tri(
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)],
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)],
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)]);
            expect(ti.test(T, unitOBB).intersect).toBe(abTI.test(T, B).intersect);
            const obbResult = fi.find(T, unitOBB);
            const abResult = abFI.find(T, B);
            expect(obbResult.insidePolygon.length)
                .toBe(abResult.insidePolygon.length);
            for (let i = 0; i < obbResult.insidePolygon.length; ++i) {
                for (let j = 0; j < 3; ++j) {
                    expect(obbResult.insidePolygon[i].values[j])
                        .toBeCloseTo(abResult.insidePolygon[i].values[j], 10);
                }
            }
        }
    });

    it('is invariant under a rigid motion of triangle and box', () => {
        const axes = rotationFrame(vec(2, -1, 0.5), 0.9);
        const translation = vec(3, -2, 4);
        const mapPoint = (p: Vector): Vector => {
            let q = translation.clone();
            for (let i = 0; i < 3; ++i) {
                q = add(q, mul(p.values[i], axes[i]));
            }
            return q;
        };
        const mapDirection = (p: Vector): Vector => {
            let q = Vector.zero(3);
            for (let i = 0; i < 3; ++i) {
                q = add(q, mul(p.values[i], axes[i]));
            }
            return q;
        };

        let seed = 90210777;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (a: number, b: number): number => a + (b - a) * rand();

        let numInside = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const T = tri(
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)],
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)],
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)]);
            const before = fi.find(T, unitOBB);
            if (before.insidePolygon.length > 0) {
                ++numInside;
            }
            const movedT = Triangle.fromVertices(mapPoint(T.v[0]),
                mapPoint(T.v[1]), mapPoint(T.v[2]));
            const movedBox = OrientedBox.fromCenterAxisExtent(
                mapPoint(unitOBB.center), unitOBB.axis.map(mapDirection),
                unitOBB.extent);
            const after = fi.find(movedT, movedBox);
            expect(after.insidePolygon.length)
                .toBe(before.insidePolygon.length);
            for (let i = 0; i < before.insidePolygon.length; ++i) {
                const expected = mapPoint(before.insidePolygon[i]);
                for (let j = 0; j < 3; ++j) {
                    expect(after.insidePolygon[i].values[j])
                        .toBeCloseTo(expected.values[j], 8);
                }
            }
            expect(ti.test(movedT, movedBox).intersect)
                .toBe(ti.test(T, unitOBB).intersect);
        }
        expect(numInside).toBeGreaterThan(50);
    });

    it('agrees with TI and produces valid geometry (randomized)', () => {
        let seed = 11223344;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (a: number, b: number): number => a + (b - a) * rand();

        const box = OrientedBox.fromCenterAxisExtent(vec(0.25, -0.5, 0.75),
            rotationFrame(vec(1, 2, -1), 0.8), vec(1, 0.5, 1.5));

        let numInside = 0;
        let numClipped = 0;
        for (let trial = 0; trial < 1200; ++trial) {
            const T = tri(
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)],
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)],
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)]);
            const tiResult = ti.test(T, box).intersect;
            const fiResult = fi.find(T, box);
            if (fiResult.insidePolygon.length > 0) {
                ++numInside;
                expect(tiResult).toBe(true);
                if (fiResult.insidePolygon.length > 3) {
                    ++numClipped;
                }
                for (const p of fiResult.insidePolygon) {
                    expect(inOBB(box, p, 1e-9)).toBe(true);
                    expect(Math.abs(planeDistance(T, p))).toBeLessThan(1e-7);
                    expect(inTriangle(T, p, 1e-7)).toBe(true);
                }
            }
            for (const poly of fiResult.outsidePolygons) {
                for (const p of poly) {
                    expect(Math.abs(planeDistance(T, p))).toBeLessThan(1e-7);
                }
            }
        }
        expect(numInside).toBeGreaterThan(80);
        expect(numClipped).toBeGreaterThan(30);
    });
});

// ---------------------------------------------------------------------------
// Verification (V35): cross-checks against the triangle-box distance query,
// against the aligned-box query for an identity frame, and under rigid
// motions of the whole configuration.
// ---------------------------------------------------------------------------

describe('IntrTriangle3OrientedBox3 verification', () => {
    const tiv = new IntrTriangle3OrientedBox3TI();
    const fiv = new IntrTriangle3OrientedBox3FI();
    const alignedTI = new IntrTriangle3AlignedBox3TI();
    const alignedFI = new IntrTriangle3AlignedBox3FI();
    const distQuery = new DistTriangle3OrientedBox3();

    function twiceArea3(poly: readonly Vector[]): number {
        if (poly.length < 3) {
            return 0;
        }
        const acc = new Vector(3);
        for (let i = 1; i + 1 < poly.length; ++i) {
            const c = cross(sub(poly[i], poly[0]), sub(poly[i + 1], poly[0]));
            for (let k = 0; k < 3; ++k) {
                acc.values[k] += c.values[k];
            }
        }
        return length(acc);
    }

    const boxArb = fc.tuple(wellScaledVector(3, -4, 4),
        randomRotationFrame(3),
        fc.array(fc.double({ min: 0.2, max: 4, noNaN: true }),
            { minLength: 3, maxLength: 3 }))
        .map(([c, axis, e]) => OrientedBox.fromCenterAxisExtent(c, axis,
            Vector.fromArray(e)));

    const triArb = fc.array(wellScaledVector(3, -8, 8),
        { minLength: 3, maxLength: 3 })
        .map(vs => Triangle.fromVertices(vs[0], vs[1], vs[2]))
        .filter(t => length(cross(sub(t.v[1], t.v[0]), sub(t.v[2], t.v[0])))
            > 1);

    // A classification is stable when growing and shrinking the box by delta
    // leaves it unchanged. Grazing configurations are not stable, and for
    // those the answer is decided by rounding, so properties comparing two
    // computations of the same predicate are asserted only on stable input.
    function grow(b: OrientedBox, delta: number): OrientedBox {
        return OrientedBox.fromCenterAxisExtent(b.center, b.axis,
            Vector.fromArray([b.extent.values[0] + delta,
                b.extent.values[1] + delta, b.extent.values[2] + delta]));
    }

    function stable(t: Triangle, b: OrientedBox, delta = 1e-6): boolean {
        return tiv.test(t, grow(b, delta)).intersect
            === tiv.test(t, grow(b, -delta)).intersect;
    }

    it('TI is true exactly when the triangle-box distance is zero', () => {
        check(fc.tuple(triArb, boxArb), ([t, b]) => {
            const d = distQuery.compute(t, b).distance;
            if (d === 0) {
                if (stable(t, b)) {
                    expect(tiv.test(t, b).intersect).toBe(true);
                }
            }
            else if (d > 1e-9) {
                expect(tiv.test(t, b).intersect).toBe(false);
            }
        });
    });

    it('equals the aligned-box query for a box with the standard frame',
        () => {
            check(fc.tuple(triArb, wellScaledVector(3, -4, 4),
                fc.array(fc.double({ min: 0.2, max: 4, noNaN: true }),
                    { minLength: 3, maxLength: 3 })), ([t, c, e]) => {
                    const abb = AlignedBox.fromMinMax(
                        sub(c, Vector.fromArray(e)),
                        add(c, Vector.fromArray(e)));
                    // Derive the oriented box from the aligned box exactly as
                    // the aligned-box query does, so that both queries see
                    // bit-identical centres and extents; recomputing them from
                    // c and e would differ in the last bit and turn grazing
                    // configurations into spurious disagreements.
                    const center = mul(0.5, add(abb.max, abb.min));
                    const extent = mul(0.5, sub(abb.max, abb.min));
                    const obb = OrientedBox.fromCenterAxisExtent(center,
                        [Vector.unit(3, 0), Vector.unit(3, 1),
                            Vector.unit(3, 2)], extent);
                    if (!stable(t, obb)) {
                        return;
                    }
                    expect(tiv.test(t, obb).intersect)
                        .toBe(alignedTI.test(t, abb).intersect);
                    const ro = fiv.find(t, obb);
                    const ra = alignedFI.find(t, abb);
                    expect(ro.insidePolygon.length > 0)
                        .toBe(ra.insidePolygon.length > 0);
                    expectClose(twiceArea3(ro.insidePolygon),
                        twiceArea3(ra.insidePolygon), 1e-7, 1e-7);
                });
        });

    it('FI conserves the triangle area and stays inside the box', () => {
        check(fc.tuple(triArb, boxArb), ([t, b]) => {
            const r = fiv.find(t, b);
            let total = twiceArea3(r.insidePolygon);
            for (const poly of r.outsidePolygons) {
                total += twiceArea3(poly);
            }
            expectClose(total,
                length(cross(sub(t.v[1], t.v[0]), sub(t.v[2], t.v[0]))),
                1e-7, 1e-7);
            for (const x of r.insidePolygon) {
                const diff = sub(x, b.center);
                for (let k = 0; k < 3; ++k) {
                    expect(Math.abs(dot(b.axis[k], diff)))
                        .toBeLessThan(b.extent.values[k] + 1e-8);
                }
            }
        });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(triArb, boxArb, randomRotationFrame(3),
            wellScaledVector(3, -5, 5)), ([t, b, frame, tr]) => {
                if (!stable(t, b, 1e-5)) {
                    return;
                }
                const xfDir = (v: Vector): Vector => {
                    const w = new Vector(3);
                    for (let i = 0; i < 3; ++i) {
                        w.values[i] = frame[0].values[i] * v.values[0]
                            + frame[1].values[i] * v.values[1]
                            + frame[2].values[i] * v.values[2];
                    }
                    return w;
                };
                const xf = (v: Vector): Vector => add(xfDir(v), tr);
                const t2 = Triangle.fromVertices(xf(t.v[0]), xf(t.v[1]),
                    xf(t.v[2]));
                const b2 = OrientedBox.fromCenterAxisExtent(xf(b.center),
                    [xfDir(b.axis[0]), xfDir(b.axis[1]), xfDir(b.axis[2])],
                    b.extent);
                expect(tiv.test(t2, b2).intersect)
                    .toBe(tiv.test(t, b).intersect);
                const r = fiv.find(t, b);
                const r2 = fiv.find(t2, b2);
                // A vertex lying on a face plane can survive the motion or
                // not, so compare the clipped area rather than the vertices.
                expect(r2.insidePolygon.length > 0)
                    .toBe(r.insidePolygon.length > 0);
                expectClose(twiceArea3(r2.insidePolygon),
                    twiceArea3(r.insidePolygon), 1e-7, 1e-7);
            });
    });
});
