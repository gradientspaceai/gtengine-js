import { describe, it, expect } from 'vitest';
import { OrientedBox } from '../src/OrientedBox.js';
import { Sector2 } from '../src/Sector2.js';
import { Vector, add, mul, normalize } from '../src/Vector.js';
import { IntrOrientedBox2Sector2TI } from '../src/IntrOrientedBox2Sector2.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function box(center: number[], angle: number, extent: number[]): OrientedBox {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(vec(center),
        [vec([c, s]), vec([-s, c])], vec(extent));
}

function sector(vertex: number[], radius: number, direction: number[],
    angle: number): Sector2 {
    const d = vec(direction);
    normalize(d);
    return Sector2.fromVertexRadiusDirectionAngle(vec(vertex), radius, d, angle);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrOrientedBox2Sector2', () => {
    const ti = new IntrOrientedBox2Sector2TI();

    // A wedge of half-angle pi/6 about +x, radius 3, vertex at the origin.
    const wedge = sector([0, 0], 3, [1, 0], Math.PI / 6);

    it('reports an intersection when the sector vertex is inside the box', () => {
        const b = box([0, 0], 0.3, [1, 1]);
        expect(ti.test(b, wedge).intersect).toBe(true);
    });

    it('reports an intersection for a box inside the wedge', () => {
        const b = box([2, 0], 0, [0.2, 0.2]);
        expect(ti.test(b, wedge).intersect).toBe(true);
    });

    it('rejects a box beyond the sector radius', () => {
        const b = box([6, 0], 0, [0.5, 0.5]);
        expect(ti.test(b, wedge).intersect).toBe(false);
    });

    it('rejects a box behind the sector vertex', () => {
        const b = box([-4, 0], 0, [1, 1]);
        expect(ti.test(b, wedge).intersect).toBe(false);
    });

    it('rejects a box outside the angular wedge', () => {
        // The half-angle is pi/6, so the boundary at x = 2 is at
        // y = 2*tan(pi/6) ~ 1.1547. A small box at y = 3 is well outside.
        const b = box([2, 3], 0, [0.3, 0.3]);
        expect(ti.test(b, wedge).intersect).toBe(false);
    });

    it('accepts a box that straddles a wedge boundary', () => {
        const b = box([2, 1.6], 0, [0.6, 0.6]);
        expect(ti.test(b, wedge).intersect).toBe(true);
    });

    it('works for a half-angle of exactly pi/2 (a halfplane wedge)', () => {
        const half = sector([0, 0], 3, [1, 0], Math.PI / 2);
        expect(ti.test(box([1, 2], 0.3, [0.2, 0.2]), half).intersect).toBe(true);
        expect(ti.test(box([-1, 2], 0.3, [0.2, 0.2]), half).intersect).toBe(false);
        expect(ti.test(box([1, 5], 0.3, [0.2, 0.2]), half).intersect).toBe(false);
    });

    it('is only valid for half-angles up to pi/2 (upstream limitation)', () => {
        // The Sector2 default has half-angle pi, which is a full disk, but
        // the query treats the wedge as the intersection of two halfplanes.
        // A point of the box is inside the disk, yet the query rejects it.
        const disk = new Sector2();  // vertex (0,0), radius 1, angle pi
        const inside = box([0.5, 0.5], 0.4, [0.1, 0.1]);
        expect(disk.contains(inside.center)).toBe(true);
        expect(ti.test(inside, disk).intersect).toBe(false);

        // A box outside the disk radius is still rejected, as it should be.
        const outside = box([5, 5], 0.4, [0.1, 0.1]);
        expect(ti.test(outside, disk).intersect).toBe(false);
    });

    it('keeps the polygon when the second clip needs no clipping', () => {
        // Regression test for the upstream bug fixed in the port: after the
        // first boundary clip, the remaining polygon can lie entirely inside
        // the second boundary halfplane. IntrHalfspace2Polygon2 then reports
        // intersect = true with an empty polygon, and upstream would discard
        // the polygon and report "no intersection".
        const s = Sector2.fromVertexRadiusDirectionAngle(
            vec([1.4273504056036472, 1.3781593022868037]),
            2.2784554166719317,
            vec([0.9202844644704273, -0.39124992582028567]),
            1.1654300159387743);
        const b = OrientedBox.fromCenterAxisExtent(
            vec([1.064378826878965, 0.5832185461185873]),
            [vec([-0.6157519184012874, 0.7879400833725456]),
                vec([-0.7879400833725456, -0.6157519184012874])],
            vec([0.4722247305791825, 0.3057352866046131]));
        // Part of the box really does lie in the sector.
        let found = false;
        for (let i = 0; i <= 40 && !found; ++i) {
            const a0 = b.extent.values[0] * (-1 + i / 20);
            for (let j = 0; j <= 40; ++j) {
                const a1 = b.extent.values[1] * (-1 + j / 20);
                const p = add(b.center,
                    add(mul(a0, b.axis[0]), mul(a1, b.axis[1])));
                if (s.contains(p)) {
                    found = true;
                    break;
                }
            }
        }
        expect(found).toBe(true);
        expect(ti.test(b, s).intersect).toBe(true);
    });

    it('reports an intersection whenever a sampled box point is in the sector', () => {
        const rand = makeRandom(19937);
        let hits = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const s = sector(
                [4 * rand() - 2, 4 * rand() - 2],
                0.5 + 2 * rand(),
                [2 * rand() - 1, 2 * rand() - 1],
                // The algorithm is valid only for half-angles up to pi/2.
                0.2 + rand() * (Math.PI / 2 - 0.2));
            const b = box(
                [6 * rand() - 3, 6 * rand() - 3],
                rand() * Math.PI,
                [0.2 + rand(), 0.2 + rand()]);
            const result = ti.test(b, s).intersect;

            // Sample the box on a fine grid; any sample inside the sector
            // forces the query to report an intersection.
            let found = false;
            const n = 60;
            for (let i = 0; i <= n && !found; ++i) {
                const a0 = -b.extent.values[0] + (2 * b.extent.values[0] * i) / n;
                for (let j = 0; j <= n; ++j) {
                    const a1 = -b.extent.values[1]
                        + (2 * b.extent.values[1] * j) / n;
                    const p = add(b.center,
                        add(mul(a0, b.axis[0]), mul(a1, b.axis[1])));
                    if (s.contains(p)) {
                        found = true;
                        break;
                    }
                }
            }

            if (found) {
                ++hits;
                expect(result).toBe(true);
            }
        }
        expect(hits).toBeGreaterThan(20);
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, rotationFrame, unitVector, wellScaledVector, seededRandom
} from './helpers/arbitraries.js';
import { dot, sub, length } from '../src/Vector.js';

// Upstream models the wedge as the intersection of two halfplanes, which is
// valid only for a half-angle of at most pi/2.
const sectorV31 = fc.tuple(wellScaledVector(2, -4, 4), unitVector(2),
    fc.double({ min: 0.3, max: 4, noNaN: true }),
    fc.double({ min: 0.05, max: Math.PI / 2, noNaN: true }))
    .map(([v, d, r, a]) =>
        Sector2.fromVertexRadiusDirectionAngle(v, r, d, a));

const boxV31 = fc.tuple(wellScaledVector(2, -4, 4), rotationFrame(2),
    fc.double({ min: 0.1, max: 2, noNaN: true }),
    fc.double({ min: 0.1, max: 2, noNaN: true }))
    .map(([c, R, e0, e1]) => OrientedBox.fromCenterAxisExtent(c, R,
        Vector.fromArray([e0, e1])));

// Membership in the solid sector, from the definition in Sector2.h.
function inSectorV31(s: Sector2, p: Vector): boolean {
    const q = sub(p, s.vertex);
    const len = length(q);
    if (len > s.radius) {
        return false;
    }
    if (len === 0) {
        return true;
    }
    return dot(s.direction, q) >= s.cosAngle * len;
}

describe('IntrOrientedBox2Sector2 verification', () => {
    const tiq = new IntrOrientedBox2Sector2TI();

    it('reports an intersection whenever a sampled box point is inside', () => {
        check(fc.tuple(boxV31, sectorV31,
            fc.double({ min: -1, max: 1, noNaN: true }),
            fc.double({ min: -1, max: 1, noNaN: true })),
            ([b, s, s0, s1]) => {
                const p = add(b.center,
                    add(mul(s0 * b.extent.values[0], b.axis[0]),
                        mul(s1 * b.extent.values[1], b.axis[1])));
                if (!inSectorV31(s, p)) {
                    return;
                }
                expect(tiq.test(b, s).intersect).toBe(true);
            });
    });

    it('agrees with a dense sampling of the box', () => {
        const rnd = seededRandom(0x5b71d0c6);
        for (let trial = 0; trial < 200; ++trial) {
            const ang = rnd() * 2 * Math.PI;
            const axis = [
                Vector.fromArray([Math.cos(ang), Math.sin(ang)]),
                Vector.fromArray([-Math.sin(ang), Math.cos(ang)])];
            const b = OrientedBox.fromCenterAxisExtent(
                Vector.fromArray([rnd() * 6 - 3, rnd() * 6 - 3]), axis,
                Vector.fromArray([0.2 + rnd(), 0.2 + rnd()]));
            const sa = rnd() * 2 * Math.PI;
            const s = Sector2.fromVertexRadiusDirectionAngle(
                Vector.fromArray([rnd() * 6 - 3, rnd() * 6 - 3]),
                0.5 + rnd() * 3,
                Vector.fromArray([Math.cos(sa), Math.sin(sa)]),
                0.1 + rnd() * (Math.PI / 2 - 0.1));
            const got = tiq.test(b, s).intersect;
            let anyInside = false;
            let minGap = Infinity;
            for (let i = 0; i <= 40 && !anyInside; ++i) {
                for (let j = 0; j <= 40; ++j) {
                    const p = add(b.center, add(
                        mul((2 * i / 40 - 1) * b.extent.values[0], b.axis[0]),
                        mul((2 * j / 40 - 1) * b.extent.values[1],
                            b.axis[1])));
                    if (inSectorV31(s, p)) {
                        anyInside = true;
                        break;
                    }
                    // Distance from the sampled point to the sector vertex,
                    // used only as a coarse separation witness.
                    minGap = Math.min(minGap,
                        length(sub(p, s.vertex)) - s.radius);
                }
            }
            if (anyInside) {
                expect(got).toBe(true);
            }
            if (!anyInside && minGap > 0.5) {
                // Every sampled box point is well beyond the sector radius,
                // and the box is convex, so the whole box is outside.
                expect(got).toBe(false);
            }
        }
    }, 30000);

    it('a box containing the sector vertex always intersects', () => {
        check(fc.tuple(sectorV31, rotationFrame(2),
            fc.double({ min: 0.1, max: 2, noNaN: true }),
            fc.double({ min: 0.1, max: 2, noNaN: true }),
            fc.double({ min: -0.8, max: 0.8, noNaN: true }),
            fc.double({ min: -0.8, max: 0.8, noNaN: true })),
            ([s, R, e0, e1, f0, f1]) => {
                // Place the box so that the sector vertex is strictly inside.
                const center = sub(s.vertex, add(mul(f0 * e0, R[0]),
                    mul(f1 * e1, R[1])));
                const b = OrientedBox.fromCenterAxisExtent(center, R,
                    Vector.fromArray([e0, e1]));
                expect(tiq.test(b, s).intersect).toBe(true);
            });
    });

    it('a box beyond the sector radius never intersects', () => {
        check(fc.tuple(sectorV31, rotationFrame(2), unitVector(2),
            fc.double({ min: 0.1, max: 1, noNaN: true }),
            fc.double({ min: 0.1, max: 1, noNaN: true }),
            fc.double({ min: 0.5, max: 3, noNaN: true })),
            ([s, R, d, e0, e1, extra]) => {
                const dist = s.radius + Math.hypot(e0, e1) + extra;
                const b = OrientedBox.fromCenterAxisExtent(
                    add(s.vertex, mul(dist, d)), R,
                    Vector.fromArray([e0, e1]));
                expect(tiq.test(b, s).intersect).toBe(false);
            });
    });

    it('is invariant under a common rigid motion', () => {
        check(fc.tuple(boxV31, sectorV31,
            fc.double({ min: -Math.PI, max: Math.PI, noNaN: true }),
            wellScaledVector(2, -4, 4)), ([b, s, ang, tr]) => {
            const ca = Math.cos(ang), sa = Math.sin(ang);
            const rot = (p: Vector): Vector => Vector.fromArray([
                ca * p.values[0] - sa * p.values[1],
                sa * p.values[0] + ca * p.values[1]]);
            const xf = (p: Vector): Vector => add(tr, rot(p));
            const b2 = OrientedBox.fromCenterAxisExtent(xf(b.center),
                [rot(b.axis[0]), rot(b.axis[1])], b.extent);
            const s2 = Sector2.fromVertexRadiusDirectionAngle(xf(s.vertex),
                s.radius, rot(s.direction), s.angle);
            const a = tiq.test(b, s).intersect;
            const c = tiq.test(b2, s2).intersect;
            if (a === c) {
                return;
            }
            // A disagreement is acceptable only for a configuration that is
            // (numerically) touching: some box corner sits on the sector
            // boundary. Verify that with an independent membership test.
            let minAbs = Infinity;
            for (const [i, j] of [[-1, -1], [1, -1], [1, 1], [-1, 1]]) {
                const p = add(b.center,
                    add(mul(i * b.extent.values[0], b.axis[0]),
                        mul(j * b.extent.values[1], b.axis[1])));
                const q = sub(p, s.vertex);
                const len = length(q);
                minAbs = Math.min(minAbs, Math.abs(len - s.radius),
                    Math.abs(dot(s.direction, q) - s.cosAngle * len));
            }
            expect(minAbs).toBeLessThan(1e-6);
        });
    });

    it('keeps the clipped polygon when a clip is a no-op (upstream #200)', () => {
        // Regression for the fixed upstream defect: when the first clip
        // leaves the polygon entirely inside the second halfplane, upstream
        // overwrites the polygon with the empty "no clipping necessary"
        // result and reports a false negative. Sweep a band of boxes that
        // sit inside the wedge and confirm every one is reported.
        const s = Sector2.fromVertexRadiusDirectionAngle(
            Vector.fromArray([0, 0]), 10, Vector.fromArray([1, 0]),
            Math.PI / 3);
        const axis = [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])];
        for (let k = 0; k <= 40; ++k) {
            const cx = 1 + (7 * k) / 40;
            const b = OrientedBox.fromCenterAxisExtent(
                Vector.fromArray([cx, 0]), axis, Vector.fromArray([0.3, 0.3]));
            expect(tiq.test(b, s).intersect).toBe(true);
        }
    });
});
