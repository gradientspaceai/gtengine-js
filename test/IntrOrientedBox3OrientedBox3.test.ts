import { describe, it, expect } from 'vitest';
import {
    IntrOrientedBox3OrientedBox3TI
} from '../src/IntrOrientedBox3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { check, fc, positive, rotationFrame, wellScaled } from './helpers/arbitraries.js';

const ti = new IntrOrientedBox3OrientedBox3TI();

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The columns of Rz(c)*Ry(b)*Rx(a), an orthonormal right-handed frame.
function rotationAxes(a: number, b: number, c: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    const cc = Math.cos(c), sc = Math.sin(c);
    const m = [
        [cc * cb, cc * sb * sa - sc * ca, cc * sb * ca + sc * sa],
        [sc * cb, sc * sb * sa + cc * ca, sc * sb * ca - cc * sa],
        [-sb, cb * sa, cb * ca]
    ];
    return [
        v3(m[0][0], m[1][0], m[2][0]),
        v3(m[0][1], m[1][1], m[2][1]),
        v3(m[0][2], m[1][2], m[2][2])
    ];
}

const stdAxes = [v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)];

function makeBox(center: Vector, axes: Vector[], extent: number[]): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(center, axes,
        Vector.fromArray(extent));
}

// The projection radius of a box onto the unit-length direction 'dir'.
function radius(box: OrientedBox, dir: Vector): number {
    let r = 0;
    for (let k = 0; k < 3; ++k) {
        r += box.extent.values[k] * Math.abs(dot(box.axis[k], dir));
    }
    return r;
}

// An independent separating-axis oracle that projects the box vertices onto
// the 15 candidate axes. It returns null when a projection is within
// 'tolerance' of touching, so borderline configurations can be skipped.
function oracleSeparated(box0: OrientedBox, box1: OrientedBox,
    tolerance: number): boolean | null {
    const v0 = box0.getVertices();
    const v1 = box1.getVertices();
    const axes: Vector[] = [];
    for (let i = 0; i < 3; ++i) {
        axes.push(box0.axis[i]);
        axes.push(box1.axis[i]);
    }
    for (let i = 0; i < 3; ++i) {
        for (let j = 0; j < 3; ++j) {
            const c = cross(box0.axis[i], box1.axis[j]);
            if (Math.sqrt(dot(c, c)) > 1e-6) {
                normalize(c);
                axes.push(c);
            }
        }
    }

    let separated = false;
    let borderline = false;
    for (const axis of axes) {
        let min0 = Infinity, max0 = -Infinity;
        for (const p of v0) {
            const d = dot(p, axis);
            min0 = Math.min(min0, d);
            max0 = Math.max(max0, d);
        }
        let min1 = Infinity, max1 = -Infinity;
        for (const p of v1) {
            const d = dot(p, axis);
            min1 = Math.min(min1, d);
            max1 = Math.max(max1, d);
        }
        const gap = Math.max(min0 - max1, min1 - max0);
        if (Math.abs(gap) < tolerance) {
            borderline = true;
        }
        if (gap > 0) {
            separated = true;
        }
    }

    if (borderline) {
        return null;
    }
    return separated;
}

describe('IntrOrientedBox3OrientedBox3TI basic configurations', () => {
    const unit0 = makeBox(v3(0, 0, 0), stdAxes, [1, 1, 1]);

    it('reports coincident boxes as intersecting', () => {
        expect(ti.test(unit0, unit0).intersect).toBe(true);
    });

    it('reports a contained box as intersecting', () => {
        const inner = makeBox(v3(0.1, -0.2, 0.05), stdAxes, [0.25, 0.25, 0.25]);
        expect(ti.test(unit0, inner).intersect).toBe(true);
        expect(ti.test(inner, unit0).intersect).toBe(true);
    });

    it('treats face-touching boxes as intersecting (solid, closed boxes)', () => {
        const touching = makeBox(v3(2, 0, 0), stdAxes, [1, 1, 1]);
        expect(ti.test(unit0, touching).intersect).toBe(true);
    });

    it('reports separated boxes and the separating face normal', () => {
        const apart = makeBox(v3(2.001, 0, 0), stdAxes, [1, 1, 1]);
        const result = ti.test(unit0, apart);
        expect(result.intersect).toBe(false);
        expect(result.separating).toEqual([0, -1]);
    });

    it('handles zero-extent boxes (points)', () => {
        const insidePoint = makeBox(v3(0.5, 0.5, 0.5), stdAxes, [0, 0, 0]);
        expect(ti.test(unit0, insidePoint).intersect).toBe(true);
        const outsidePoint = makeBox(v3(1.5, 0, 0), stdAxes, [0, 0, 0]);
        const result = ti.test(unit0, outsidePoint);
        expect(result.intersect).toBe(false);
        expect(result.separating).toEqual([0, -1]);
        const boundaryPoint = makeBox(v3(1, 0, 0), stdAxes, [0, 0, 0]);
        expect(ti.test(unit0, boundaryPoint).intersect).toBe(true);
    });

    it('is symmetric in its arguments', () => {
        const rand = makeRandom(42);
        for (let trial = 0; trial < 200; ++trial) {
            const axes1 = rotationAxes(rand(), rand(), rand());
            const box1 = makeBox(
                v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2), axes1,
                [0.3 + rand(), 0.3 + rand(), 0.3 + rand()]);
            expect(ti.test(unit0, box1).intersect)
                .toBe(ti.test(box1, unit0).intersect);
        }
    });
});

describe('IntrOrientedBox3OrientedBox3TI face-normal axes', () => {
    const axes1 = rotationAxes(0.3, 0.5, 0.7);
    const extent0 = [1, 1, 1];
    const extent1 = [0.7, 1.3, 0.9];
    const box0 = makeBox(v3(0, 0, 0), stdAxes, extent0);

    for (let k = 0; k < 3; ++k) {
        it(`separates on box0.axis[${k}]`, () => {
            const dir = stdAxes[k];
            const proto = makeBox(v3(0, 0, 0), axes1, extent1);
            const threshold = extent0[k] + radius(proto, dir);
            const separated = makeBox(mul(threshold * 1.001, dir), axes1, extent1);
            const result = ti.test(box0, separated);
            expect(result.intersect).toBe(false);
            expect(result.separating).toEqual([k, -1]);
            // Well inside the threshold the boxes intersect.
            const overlapping = makeBox(mul(threshold * 0.5, dir), axes1, extent1);
            expect(ti.test(box0, overlapping).intersect).toBe(true);
        });
    }

    for (let k = 0; k < 3; ++k) {
        it(`separates on box1.axis[${k}]`, () => {
            const dir = axes1[k];
            const threshold = radius(box0, dir) + extent1[k];
            const separated = makeBox(mul(threshold * 1.001, dir), axes1, extent1);
            const result = ti.test(box0, separated);
            expect(result.intersect).toBe(false);
            expect(result.separating).toEqual([-1, k]);
            const overlapping = makeBox(mul(threshold * 0.5, dir), axes1, extent1);
            expect(ti.test(box0, overlapping).intersect).toBe(true);
        });
    }
});

describe('IntrOrientedBox3OrientedBox3TI edge-cross axes', () => {
    // Each entry gives box1 Euler angles and extents for which the axis
    // Cross(box0.axis[i], box1.axis[j]) is the first separating axis when the
    // boxes are pushed just past the critical separation distance.
    const cases: Array<{ i: number, j: number, angles: number[], extent: number[] }> = [
        { i: 0, j: 0, angles: [0.3, 0.5, 0.7], extent: [0.7, 1.3, 0.9] },
        { i: 0, j: 1, angles: [0.3, 0.5, 0.7], extent: [0.7, 1.3, 0.9] },
        { i: 0, j: 2, angles: [0.9, 0.35, 0.15], extent: [0.7, 1.3, 0.9] },
        { i: 1, j: 0, angles: [0.3, 0.5, 0.7], extent: [1.5, 0.5, 1.0] },
        { i: 1, j: 1, angles: [0.3, 0.5, 0.7], extent: [0.7, 1.3, 0.9] },
        { i: 1, j: 2, angles: [0.3, 0.5, 0.7], extent: [0.7, 1.3, 0.9] },
        { i: 2, j: 0, angles: [0.3, 0.5, 0.7], extent: [0.7, 1.3, 0.9] },
        { i: 2, j: 1, angles: [0.3, 0.5, 0.7], extent: [0.7, 1.3, 0.9] },
        { i: 2, j: 2, angles: [0.3, 0.5, 0.7], extent: [1.5, 0.5, 1.0] }
    ];

    const box0 = makeBox(v3(0, 0, 0), stdAxes, [1, 1, 1]);

    for (const c of cases) {
        it(`separates on Cross(box0.axis[${c.i}], box1.axis[${c.j}])`, () => {
            const axes1 = rotationAxes(c.angles[0], c.angles[1], c.angles[2]);
            const proto = makeBox(v3(0, 0, 0), axes1, c.extent);
            const L = cross(stdAxes[c.i], axes1[c.j]);
            normalize(L);
            const threshold = radius(box0, L) + radius(proto, L);

            const separated = makeBox(mul(threshold * 1.001, L), axes1, c.extent);
            const result = ti.test(box0, separated);
            expect(result.intersect).toBe(false);
            expect(result.separating).toEqual([c.i, c.j]);

            // Just inside the threshold the boxes intersect (no face normal
            // separates there either).
            const overlapping = makeBox(mul(threshold * 0.999, L), axes1, c.extent);
            expect(ti.test(box0, overlapping).intersect).toBe(true);
        });
    }
});

describe('IntrOrientedBox3OrientedBox3TI epsilon handling', () => {
    const box0 = makeBox(v3(0, 0, 0), stdAxes, [1, 1, 1]);
    const axes1 = rotationAxes(0.3, 0.5, 0.7);
    const extent1 = [0.7, 1.3, 0.9];
    const L = (() => {
        const c = cross(stdAxes[0], axes1[0]);
        normalize(c);
        return c;
    })();
    const proto = makeBox(v3(0, 0, 0), axes1, extent1);
    const threshold = radius(box0, L) + radius(proto, L);
    const separated = makeBox(mul(threshold * 1.001, L), axes1, extent1);

    it('reports the edge-cross separation with the default epsilon', () => {
        expect(ti.test(box0, separated).intersect).toBe(false);
    });

    it('clamps a negative epsilon to zero (same as the default)', () => {
        const withNegative = ti.test(box0, separated, -0.5);
        const withZero = ti.test(box0, separated, 0);
        expect(withNegative.intersect).toBe(withZero.intersect);
        expect(withNegative.separating).toEqual(withZero.separating);
    });

    it('skips the edge-cross axes when epsilon makes the axes "parallel"', () => {
        // With epsilon = 1 the cutoff is 0, so every axis pair looks parallel
        // and only the six face normals are tested. These boxes are separated
        // only by an edge-cross axis, so the query now reports intersection.
        expect(ti.test(box0, separated, 1).intersect).toBe(true);
    });

    it('does not affect face-normal separation', () => {
        const apart = makeBox(v3(5, 0, 0), stdAxes, [1, 1, 1]);
        expect(ti.test(box0, apart, 1).intersect).toBe(false);
    });

    it('takes the parallel-pair shortcut for axis-aligned boxes', () => {
        // Both boxes use the standard axes, so |Dot(A0[i],A1[i])| = 1 > cutoff
        // for epsilon = 0 and the edge-cross axes are skipped.
        const other = makeBox(v3(1.5, 1.5, 1.5), stdAxes, [1, 1, 1]);
        expect(ti.test(box0, other).intersect).toBe(true);
    });
});

describe('IntrOrientedBox3OrientedBox3TI randomized cross-check', () => {
    it('agrees with a vertex-projection separating-axis oracle', () => {
        const rand = makeRandom(1234567);
        let intersectCount = 0, separatedCount = 0;
        for (let trial = 0; trial < 800; ++trial) {
            const axes0 = rotationAxes(6 * rand(), 6 * rand(), 6 * rand());
            const axes1 = rotationAxes(6 * rand(), 6 * rand(), 6 * rand());
            const box0 = makeBox(
                v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1), axes0,
                [0.3 + rand(), 0.3 + rand(), 0.3 + rand()]);
            const box1 = makeBox(
                v3(3 * rand() - 1.5, 3 * rand() - 1.5, 3 * rand() - 1.5), axes1,
                [0.3 + rand(), 0.3 + rand(), 0.3 + rand()]);

            const expected = oracleSeparated(box0, box1, 1e-9);
            if (expected === null) {
                continue;
            }

            const result = ti.test(box0, box1);
            expect(result.intersect).toBe(!expected);
            if (result.intersect) {
                ++intersectCount;
            }
            else {
                ++separatedCount;
            }
        }
        expect(intersectCount).toBeGreaterThan(20);
        expect(separatedCount).toBeGreaterThan(20);
    });

    it('reports intersection whenever a common point exists', () => {
        const rand = makeRandom(2468013);
        let found = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const axes0 = rotationAxes(6 * rand(), 6 * rand(), 6 * rand());
            const axes1 = rotationAxes(6 * rand(), 6 * rand(), 6 * rand());
            const e0 = [0.3 + rand(), 0.3 + rand(), 0.3 + rand()];
            const box0 = makeBox(
                v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1), axes0, e0);
            const box1 = makeBox(
                v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1), axes1,
                [0.3 + rand(), 0.3 + rand(), 0.3 + rand()]);

            // Sample points of box0 and test them against box1.
            let common = false;
            for (let s = 0; s < 60 && !common; ++s) {
                let p = box0.center.clone();
                for (let k = 0; k < 3; ++k) {
                    p = add(p, mul(e0[k] * (2 * rand() - 1), box0.axis[k]));
                }
                const delta = sub(p, box1.center);
                let inside = true;
                for (let k = 0; k < 3; ++k) {
                    if (Math.abs(dot(delta, box1.axis[k]))
                        > box1.extent.values[k]) {
                        inside = false;
                        break;
                    }
                }
                common = inside;
            }

            if (common) {
                ++found;
                expect(ti.test(box0, box1).intersect).toBe(true);
            }
        }
        expect(found).toBeGreaterThan(10);
    });
});

describe('IntrOrientedBox3OrientedBox3 verification', () => {
    const boxArb = fc.tuple(
        fc.array(wellScaled(-4, 4), { minLength: 3, maxLength: 3 }),
        rotationFrame(3),
        fc.array(positive(3), { minLength: 3, maxLength: 3 }))
        .map(([c, axes, ext]) => makeBox(Vector.fromArray(c), axes, ext));

    it('matches the 15-axis oracle whenever the oracle is decisive', () => {
        check(fc.tuple(boxArb, boxArb), ([box0, box1]) => {
            const expected = oracleSeparated(box0, box1, 1e-9);
            if (expected === null) {
                return;
            }
            expect(ti.test(box0, box1).intersect).toBe(!expected);
        });
    });

    it('is symmetric under argument swap away from the touching boundary', () => {
        check(fc.tuple(boxArb, boxArb), ([box0, box1]) => {
            if (oracleSeparated(box0, box1, 1e-9) === null) {
                return;
            }
            expect(ti.test(box1, box0).intersect)
                .toBe(ti.test(box0, box1).intersect);
        });
    });

    it('the reported separating axis really separates the two boxes', () => {
        check(fc.tuple(boxArb, boxArb), ([box0, box1]) => {
            const r = ti.test(box0, box1);
            if (r.intersect) {
                return;
            }
            const [i0, i1] = r.separating;
            let axis: Vector;
            if (i1 < 0) {
                axis = box0.axis[i0];
            } else if (i0 < 0) {
                axis = box1.axis[i1];
            } else {
                axis = cross(box0.axis[i0], box1.axis[i1]);
                if (Math.sqrt(dot(axis, axis)) < 1e-8) {
                    // A degenerate edge-edge axis cannot report separation
                    // (both sides of the test are zero), so nothing to check.
                    return;
                }
                normalize(axis);
            }
            const c0 = dot(box0.center, axis), c1 = dot(box1.center, axis);
            const r0 = radius(box0, axis), r1 = radius(box1, axis);
            const scale = 1 + Math.abs(c0) + Math.abs(c1) + r0 + r1;
            expect(Math.abs(c1 - c0)).toBeGreaterThan(r0 + r1 - 1e-11 * scale);
        });
    });

    it('a point common to both boxes forces intersect = true', () => {
        const rnd = makeRandom(0x51ed270b);
        check(fc.tuple(boxArb, boxArb), ([box0, box1]) => {
            if (ti.test(box0, box1).intersect) {
                return;
            }
            const e = box1.extent.values;
            for (let k = 0; k < 150; ++k) {
                const s = [(2 * rnd() - 1) * e[0], (2 * rnd() - 1) * e[1],
                    (2 * rnd() - 1) * e[2]];
                let p = box1.center.clone();
                for (let d = 0; d < 3; ++d) {
                    p = add(p, mul(s[d], box1.axis[d]));
                }
                const delta = sub(p, box0.center);
                let inside = true;
                for (let d = 0; d < 3; ++d) {
                    if (Math.abs(dot(delta, box0.axis[d]))
                        > box0.extent.values[d]) {
                        inside = false;
                        break;
                    }
                }
                expect(inside).toBe(false);
            }
        }, 50);
    }, 30000);

    it('is equivariant under a rigid motion applied to both boxes', () => {
        check(fc.tuple(boxArb, boxArb, rotationFrame(3),
            fc.array(wellScaled(-3, 3), { minLength: 3, maxLength: 3 })),
            ([box0, box1, R, t]) => {
                if (oracleSeparated(box0, box1, 1e-8) === null) {
                    return;
                }
                const move = (b: OrientedBox): OrientedBox => {
                    const apply = (v: Vector): Vector => {
                        const out = [0, 0, 0];
                        for (let d = 0; d < 3; ++d) {
                            out[d] = R[0].values[d] * v.values[0]
                                + R[1].values[d] * v.values[1]
                                + R[2].values[d] * v.values[2];
                        }
                        return Vector.fromArray(out);
                    };
                    const c = apply(b.center);
                    return makeBox(
                        Vector.fromArray([c.values[0] + t[0],
                            c.values[1] + t[1], c.values[2] + t[2]]),
                        [apply(b.axis[0]), apply(b.axis[1]), apply(b.axis[2])],
                        [...b.extent.values]);
                };
                expect(ti.test(move(box0), move(box1)).intersect)
                    .toBe(ti.test(box0, box1).intersect);
            });
    });

    it('a negative epsilon is clamped to zero', () => {
        const b0 = makeBox(v3(0, 0, 0), stdAxes, [1, 1, 0.1]);
        const axes = rotationAxes(0, 0, Math.PI / 4);
        const b1 = makeBox(v3(0, 0, 0.3), axes, [1, 1, 0.1]);
        const zero = ti.test(b0, b1, 0).intersect;
        expect(ti.test(b0, b1, -1).intersect).toBe(zero);
        expect(ti.test(b0, b1, -1e-12).intersect).toBe(zero);
    });

    it('a large epsilon restricts separation to the six face normals', () => {
        const rnd = makeRandom(0xbeef1234);
        for (let trial = 0; trial < 400; ++trial) {
            const b0 = makeBox(v3(0, 0, 0), rotationAxes(2 * rnd(), 2 * rnd(),
                2 * rnd()), [1, 1, 1]);
            const b1 = makeBox(v3(6 * rnd() - 3, 6 * rnd() - 3, 6 * rnd() - 3),
                rotationAxes(2 * rnd(), 2 * rnd(), 2 * rnd()),
                [1 + rnd(), 1 + rnd(), 1 + rnd()]);
            const r = ti.test(b0, b1, 2);
            if (!r.intersect) {
                expect(r.separating[0] < 0 || r.separating[1] < 0).toBe(true);
            }
        }
    });
});
