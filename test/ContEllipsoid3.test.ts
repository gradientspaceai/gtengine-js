import { describe, it, expect } from 'vitest';
import {
    getContainerEllipsoid3,
    inContainerEllipsoid3,
    mergeContainersEllipsoid3
} from '../src/ContEllipsoid3.js';
import { Hyperellipsoid, type Ellipsoid3 } from '../src/Hyperellipsoid.js';
import { Line } from '../src/Line.js';
import { projectEllipsoid3 } from '../src/Projection.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, fc, rotationFrame, wellScaledVector
} from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// A right-handed orthonormal frame: the standard basis rotated about 'axis'.
function frame(axis: Vector, angle: number): Vector[] {
    const w = axis.clone();
    normalize(w);
    const rot = (u: Vector) => {
        const c = Math.cos(angle), s = Math.sin(angle);
        return add(add(mul(c, u), mul(s, cross(w, u))),
            mul((1 - c) * dot(w, u), w));
    };
    return [rot(v3(1, 0, 0)), rot(v3(0, 1, 0)), rot(v3(0, 0, 1))];
}

function ellipsoid3(center: Vector, axis: Vector[], e: Vector): Ellipsoid3 {
    return Hyperellipsoid.fromCenterAxisExtent(center, axis, e);
}

// Independent evaluation of the ellipsoid quadratic form: <= 1 inside.
function quadratic(p: Vector, e: Ellipsoid3): number {
    const d = sub(p, e.center);
    let sum = 0;
    for (let j = 0; j < 3; ++j) {
        const t = dot(d, e.axis[j]) / e.extent.values[j];
        sum += t * t;
    }
    return sum;
}

let seed = 777333111;
function rand(): number {
    seed = (seed * 1103515245 + 12345) % 2147483648;
    return seed / 2147483648;
}

describe('ContEllipsoid3', () => {
    describe('inContainerEllipsoid3', () => {
        it('matches the unit sphere', () => {
            const e = new Hyperellipsoid(3);
            expect(inContainerEllipsoid3(v3(0, 0, 0), e)).toBe(true);
            expect(inContainerEllipsoid3(v3(1, 0, 0), e)).toBe(true);
            expect(inContainerEllipsoid3(v3(0, 0, -1), e)).toBe(true);
            expect(inContainerEllipsoid3(v3(1.0001, 0, 0), e)).toBe(false);
            // sqrt(3) * 0.577 < 1 < sqrt(3) * 0.578.
            expect(inContainerEllipsoid3(v3(0.577, 0.577, 0.577), e))
                .toBe(true);
            expect(inContainerEllipsoid3(v3(0.578, 0.578, 0.578), e))
                .toBe(false);
        });

        it('matches a rotated, translated, anisotropic ellipsoid', () => {
            const e = ellipsoid3(v3(2, -3, 1), frame(v3(1, 2, -1), 0.8),
                v3(3, 1, 0.5));
            for (let j = 0; j < 3; ++j) {
                // Just inside the axis endpoint (exactly on the boundary the
                // comparison is at the mercy of round-off).
                const tip = add(e.center,
                    mul(e.extent.values[j] * 0.999, e.axis[j]));
                expect(inContainerEllipsoid3(tip, e)).toBe(true);
                const beyond = add(e.center,
                    mul(e.extent.values[j] * 1.001, e.axis[j]));
                expect(inContainerEllipsoid3(beyond, e)).toBe(false);
            }
        });

        it('rejects mismatched dimensions', () => {
            expect(() => inContainerEllipsoid3(Vector.fromArray([0, 0]),
                new Hyperellipsoid(3))).toThrow();
        });
    });

    describe('getContainerEllipsoid3', () => {
        it('contains every input point and is tight on at least one', () => {
            const points: Vector[] = [];
            for (let i = 0; i < 8; ++i) {
                points.push(v3((i & 1) ? 4 : -4, (i & 2) ? 2 : -2,
                    (i & 4) ? 1 : -1));
            }
            points.push(v3(0, 0, 0), v3(1, 0.5, 0.25));
            const e = getContainerEllipsoid3(points);
            expect(e).not.toBeNull();
            for (const p of points) {
                expect(quadratic(p, e!)).toBeLessThanOrEqual(1 + 1e-9);
            }
            expect(Math.max(...points.map(p => quadratic(p, e!))))
                .toBeCloseTo(1, 9);
        });

        it('produces a sphere for a spherically symmetric set', () => {
            const points: Vector[] = [];
            const r = 3;
            for (const s of [-1, 1]) {
                points.push(v3(s * r, 0, 0), v3(0, s * r, 0), v3(0, 0, s * r));
            }
            const e = getContainerEllipsoid3(points);
            expect(e).not.toBeNull();
            for (let j = 0; j < 3; ++j) {
                expect(e!.center.values[j]).toBeCloseTo(0, 10);
                expect(e!.extent.values[j]).toBeCloseTo(r, 8);
            }
        });

        it('yields an infinite extent for coplanar points (ill-conditioned)',
            () => {
                const points = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0),
                    v3(1, 1, 0), v3(0.3, 0.7, 0)];
                const e = getContainerEllipsoid3(points);
                expect(e).not.toBeNull();
                expect([...e!.extent.values].some(x => !Number.isFinite(x)))
                    .toBe(true);
            });

        it('throws on an empty point set', () => {
            expect(() => getContainerEllipsoid3([])).toThrow();
        });

        it('rejects non-3D points', () => {
            expect(() => getContainerEllipsoid3([Vector.fromArray([0, 0])]))
                .toThrow();
        });

        it('contains every input point for random clouds', () => {
            for (let trial = 0; trial < 30; ++trial) {
                const axis = frame(v3(2 * rand() - 1, 2 * rand() - 1,
                    2 * rand() - 1 + 0.1), 2 * Math.PI * rand());
                const c = v3(8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4);
                const ext = [0.5 + 3 * rand(), 0.5 + 3 * rand(),
                    0.5 + 3 * rand()];
                const points: Vector[] = [];
                for (let i = 0; i < 60; ++i) {
                    let p = c.clone();
                    for (let j = 0; j < 3; ++j) {
                        p = add(p, mul(ext[j] * (2 * rand() - 1), axis[j]));
                    }
                    points.push(p);
                }
                const e = getContainerEllipsoid3(points);
                expect(e).not.toBeNull();
                for (const p of points) {
                    expect(quadratic(p, e!)).toBeLessThanOrEqual(1 + 1e-9);
                }
                for (let j = 0; j < 3; ++j) {
                    expect(length(e!.axis[j])).toBeCloseTo(1, 10);
                }
                expect(dot(e!.axis[0], e!.axis[1])).toBeCloseTo(0, 10);
                expect(dot(e!.axis[0], e!.axis[2])).toBeCloseTo(0, 10);
                expect(dot(e!.axis[1], e!.axis[2])).toBeCloseTo(0, 10);
            }
        });
    });

    describe('mergeContainersEllipsoid3', () => {
        it('reproduces the input when both ellipsoids are the same', () => {
            const e = ellipsoid3(v3(1, 2, -1), frame(v3(0, 0, 1), 0.4),
                v3(3, 2, 1));
            const m = mergeContainersEllipsoid3(e, e);
            for (let j = 0; j < 3; ++j) {
                expect(m.center.values[j]).toBeCloseTo(e.center.values[j], 8);
            }
            const sorted = [...m.extent.values].sort((a, b) => a - b);
            expect(sorted[0]).toBeCloseTo(1, 8);
            expect(sorted[1]).toBeCloseTo(2, 8);
            expect(sorted[2]).toBeCloseTo(3, 8);
        });

        it('spans the projection intervals of both inputs (randomized)', () => {
            // The upstream algorithm builds the bounding box of the two
            // projected intervals along the averaged axes and uses its
            // half-widths as the merged extents.
            for (let trial = 0; trial < 40; ++trial) {
                const e0 = ellipsoid3(
                    v3(6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3),
                    frame(v3(2 * rand() - 1, 2 * rand() - 1,
                        2 * rand() - 1 + 0.1), 2 * Math.PI * rand()),
                    v3(0.5 + 2 * rand(), 0.5 + 2 * rand(), 0.5 + 2 * rand()));
                const e1 = ellipsoid3(
                    v3(6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3),
                    frame(v3(2 * rand() - 1, 2 * rand() - 1,
                        2 * rand() - 1 + 0.1), 2 * Math.PI * rand()),
                    v3(0.5 + 2 * rand(), 0.5 + 2 * rand(), 0.5 + 2 * rand()));
                const m = mergeContainersEllipsoid3(e0, e1);

                for (let j = 0; j < 3; ++j) {
                    expect(length(m.axis[j])).toBeCloseTo(1, 9);
                }
                expect(dot(m.axis[0], m.axis[1])).toBeCloseTo(0, 9);
                expect(dot(m.axis[0], m.axis[2])).toBeCloseTo(0, 9);
                expect(dot(m.axis[1], m.axis[2])).toBeCloseTo(0, 9);

                for (let j = 0; j < 3; ++j) {
                    const line = Line.fromOriginDirection(m.center, m.axis[j]);
                    const p0 = projectEllipsoid3(e0, line);
                    const p1 = projectEllipsoid3(e1, line);
                    const lo = Math.min(p0.smin, p1.smin);
                    const hi = Math.max(p0.smax, p1.smax);
                    expect(0.5 * (lo + hi)).toBeCloseTo(0, 8);
                    expect(m.extent.values[j]).toBeCloseTo(0.5 * (hi - lo), 8);
                }
            }
        });

        it('rejects mismatched dimensions', () => {
            expect(() => mergeContainersEllipsoid3(new Hyperellipsoid(3),
                new Hyperellipsoid(2))).toThrow();
        });
    });
});

// ---------------------------------------------------------------------------
// Verification pass (VERIFYING.md): property-based cross-checks of the port
// against the upstream ContEllipsoid3.h semantics.
// ---------------------------------------------------------------------------

describe('ContEllipsoid3 verification', () => {
    // Anisotropic lattice cloud so the fitted frame is unambiguous (the three
    // principal variances are about 4, 0.7 and 0.07).
    const baseGrid: Vector[] = [];
    for (let i = -3; i <= 3; ++i) {
        for (let j = -1; j <= 1; ++j) {
            for (let k = -1; k <= 1; ++k) {
                baseGrid.push(v3(i, j, 0.3 * k));
            }
        }
    }

    const quadratic = (ellipsoid: Ellipsoid3, p: Vector): number => {
        const diff = sub(p, ellipsoid.center);
        let sum = 0;
        for (let j = 0; j < 3; ++j) {
            const s = dot(diff, ellipsoid.axis[j]) / ellipsoid.extent.values[j];
            sum += s * s;
        }
        return sum;
    };

    it('contains every input point and is tight on one of them', () => {
        check(fc.array(wellScaledVector(3, -6, 6), { minLength: 3, maxLength: 12 }),
            (points: Vector[]) => {
                const ellipsoid = getContainerEllipsoid3(points);
                if (ellipsoid === null || !ellipsoid.extent.values.every(
                    e => Number.isFinite(e) && e > 0)) {
                    // Coplanar or collinear input: upstream's zero eigenvalue
                    // gives an infinite extent (issue #292).
                    return;
                }
                // The construction makes max_i Q(X_i) exactly 1 in exact
                // arithmetic, so the extremal point sits on the boundary and
                // the strict '<= 1' of InContainer can reject it by one ulp.
                // That is upstream behavior; the property is that Q is at
                // most 1 and that every point strictly inside is accepted.
                let maxQ = 0;
                for (const p of points) {
                    const q = quadratic(ellipsoid, p);
                    maxQ = Math.max(maxQ, q);
                    expect(q).toBeLessThanOrEqual(1 + 1e-9);
                    if (q < 1 - 1e-9) {
                        expect(inContainerEllipsoid3(p, ellipsoid)).toBe(true);
                    }
                }
                expectClose(maxQ, 1, 1e-9, 1e-9);
            });
    });

    it('is equivariant under rigid motions', () => {
        const reference = getContainerEllipsoid3(baseGrid)!;
        check(fc.tuple(rotationFrame(3), wellScaledVector(3)),
            ([frame, t]: [Vector[], Vector]) => {
                const xform = (p: Vector): Vector =>
                    add(add(add(mul(p.get(0), frame[0]), mul(p.get(1), frame[1])),
                        mul(p.get(2), frame[2])), t);
                const moved = getContainerEllipsoid3(baseGrid.map(xform))!;
                const want = [0, 1, 2].map(k => reference.extent.get(k))
                    .sort((a, b) => a - b);
                const got = [0, 1, 2].map(k => moved.extent.get(k))
                    .sort((a, b) => a - b);
                for (let k = 0; k < 3; ++k) {
                    expectClose(got[k], want[k], 1e-8, 1e-8);
                }
                expect(length(sub(moved.center, xform(reference.center))))
                    .toBeLessThanOrEqual(1e-8);
            });
    });

    // Upstream issue #292 item 1, the 3D case: coplanar input gives a zero
    // eigenvalue, an infinite extent, and an InContainer that accepts
    // everything. Preserved; pinned.
    it('coplanar input yields an infinite extent that accepts everything (#292)',
        () => {
            const points: Vector[] = [];
            for (let i = -2; i <= 2; ++i) {
                for (let j = -2; j <= 2; ++j) {
                    points.push(v3(i, j, 0));
                }
            }
            const ellipsoid = getContainerEllipsoid3(points)!;
            expect([0, 1, 2].some(k => !Number.isFinite(ellipsoid.extent.get(k))))
                .toBe(true);
            expect(inContainerEllipsoid3(v3(0, 0, 1e6), ellipsoid)).toBe(true);
        });

    // Upstream issue #292 item 2, the 3D case: the merged ellipsoid is
    // inscribed in the box of the projected intervals, so it does not contain
    // its inputs. Two unit spheres one unit either side of the origin.
    it('merge does not contain its inputs (#292)', () => {
        const unit = (cx: number): Ellipsoid3 =>
            Hyperellipsoid.fromCenterAxisExtent(v3(cx, 0, 0),
                [v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)], v3(1, 1, 1));
        const merge = mergeContainersEllipsoid3(unit(-1), unit(1));
        expectClose(merge.extent.get(0), 2, 1e-9, 1e-9);
        expectClose(merge.extent.get(1), 1, 1e-9, 1e-9);
        expectClose(merge.extent.get(2), 1, 1e-9, 1e-9);
        expect(inContainerEllipsoid3(v3(-1, 1, 0), unit(-1))).toBe(true);
        expect(inContainerEllipsoid3(v3(-1, 1, 0), merge)).toBe(false);
    });

    // What the merge does guarantee: its extents span the projected intervals
    // of both inputs along the merged axes, and the merged axes are
    // orthonormal.
    it('merge spans the projection intervals of both inputs', () => {
        const ellipsoidArb = fc.tuple(wellScaledVector(3, -5, 5),
            rotationFrame(3),
            fc.double({ min: 0.2, max: 3, noNaN: true }),
            fc.double({ min: 0.2, max: 3, noNaN: true }),
            fc.double({ min: 0.2, max: 3, noNaN: true }))
            .map(([c, frame, e0, e1, e2]) =>
                Hyperellipsoid.fromCenterAxisExtent(c, frame,
                    Vector.fromArray([e0, e1, e2])));
        check(fc.tuple(ellipsoidArb, ellipsoidArb),
            ([e0, e1]: [Ellipsoid3, Ellipsoid3]) => {
                const merge = mergeContainersEllipsoid3(e0, e1);
                for (let a = 0; a < 3; ++a) {
                    for (let b = a; b < 3; ++b) {
                        expectClose(dot(merge.axis[a], merge.axis[b]),
                            a === b ? 1 : 0, 1e-9, 1e-9);
                    }
                }
                for (let j = 0; j < 3; ++j) {
                    const line = Line.fromOriginDirection(merge.center,
                        merge.axis[j]);
                    const p0 = projectEllipsoid3(e0, line);
                    const p1 = projectEllipsoid3(e1, line);
                    const lo = Math.min(p0.smin, p1.smin);
                    const hi = Math.max(p0.smax, p1.smax);
                    expectClose(lo + hi, 0, 1e-7, 1e-7);
                    expectClose(merge.extent.get(j), 0.5 * (hi - lo),
                        1e-7, 1e-7);
                }
            });
    });
});
