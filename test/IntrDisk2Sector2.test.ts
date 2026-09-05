import { describe, it, expect } from 'vitest';
import { Hypersphere } from '../src/Hypersphere.js';
import { Sector2 } from '../src/Sector2.js';
import { Vector, add, mul, sub, length, normalize } from '../src/Vector.js';
import { IntrDisk2Sector2TI } from '../src/IntrDisk2Sector2.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function disk(cx: number, cy: number, r: number): Hypersphere {
    return Hypersphere.fromCenterRadius(v2(cx, cy), r);
}

function sector(vx: number, vy: number, radius: number, dirAngle: number,
    halfAngle: number): Sector2 {
    const d = v2(Math.cos(dirAngle), Math.sin(dirAngle));
    normalize(d);
    return Sector2.fromVertexRadiusDirectionAngle(v2(vx, vy), radius, d,
        halfAngle);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrDisk2Sector2', () => {
    const ti = new IntrDisk2Sector2TI();

    // A quarter-plane style sector: vertex at the origin, axis +x, half
    // angle 45 degrees, radius 2.
    const s = sector(0, 0, 2, 0, Math.PI / 4);

    it('reports intersection when the disk contains the sector vertex', () => {
        expect(ti.test(disk(0, 0, 0.5), s).intersect).toBe(true);
        expect(ti.test(disk(-0.25, 0, 0.5), s).intersect).toBe(true);
    });

    it('reports intersection when the disk is inside the sector', () => {
        expect(ti.test(disk(1, 0, 0.2), s).intersect).toBe(true);
    });

    it('reports no intersection when the disk is behind the vertex', () => {
        expect(ti.test(disk(-2, 0, 0.5), s).intersect).toBe(false);
    });

    it('reports no intersection when the disk is beyond the sector arc',
        () => {
            expect(ti.test(disk(4, 0, 1), s).intersect).toBe(false);
            // Just touching the arc at (2,0).
            expect(ti.test(disk(3, 0, 1), s).intersect).toBe(true);
        });

    it('reports no intersection when the disk is outside the cone', () => {
        // Straight up from the vertex, well away from the 45-degree cone.
        expect(ti.test(disk(0, 3, 0.5), s).intersect).toBe(false);
        // The perpendicular distance from (0.5,1) to the left boundary ray
        // (the line y = x) is 0.3536, so radius 0.3 misses and 0.4 reaches.
        expect(ti.test(disk(0.5, 1.0, 0.3), s).intersect).toBe(false);
        expect(ti.test(disk(0.5, 1.0, 0.4), s).intersect).toBe(true);
    });

    it('agrees with a brute-force sampling of the sector', () => {
        const rnd = makeRandom(24680);
        let numHit = 0, numMiss = 0;
        for (let k = 0; k < 400; ++k) {
            const sec = sector(rnd() * 4 - 2, rnd() * 4 - 2,
                0.3 + rnd() * 2, rnd() * 2 * Math.PI,
                0.05 + rnd() * (Math.PI / 2 - 0.05));
            const d = disk(rnd() * 6 - 3, rnd() * 6 - 3, 0.1 + rnd() * 1.5);
            const intersect = ti.test(d, sec).intersect;

            // Dense sampling of the solid sector: a sample inside the disk
            // forces the query to report an intersection.
            let sampleHit = false;
            const base = Math.atan2(sec.direction.values[1],
                sec.direction.values[0]);
            for (let i = 0; i <= 60 && !sampleHit; ++i) {
                const a = base - sec.angle + (2 * sec.angle * i) / 60;
                const dir = v2(Math.cos(a), Math.sin(a));
                for (let j = 0; j <= 60 && !sampleHit; ++j) {
                    const p = add(sec.vertex,
                        mul((sec.radius * j) / 60, dir));
                    if (length(sub(p, d.center)) <= d.radius) {
                        sampleHit = true;
                    }
                }
            }

            if (sampleHit) {
                expect(intersect).toBe(true);
                ++numHit;
            } else if (!intersect) {
                ++numMiss;
            }
        }
        expect(numHit).toBeGreaterThan(0);
        expect(numMiss).toBeGreaterThan(0);
    });

    it('does not report intersections that a fine sampling refutes', () => {
        // The converse check: whenever the query says the objects intersect,
        // a fine sampling of the disk should find a point in the sector (up
        // to sampling resolution). We only assert the strong direction, that
        // a clearly-separated configuration is reported as disjoint.
        const rnd = makeRandom(13579);
        for (let k = 0; k < 300; ++k) {
            const sec = sector(0, 0, 1, 0, Math.PI / 4);
            const c = v2(rnd() * 10 - 5, rnd() * 10 - 5);
            const r = 0.1 + rnd() * 0.5;
            // The sector is contained in the disk of radius 1 at the origin,
            // so a disk whose closest point to the origin is farther than 1
            // cannot intersect.
            if (length(c) - r > 1 + 1e-9) {
                expect(ti.test(Hypersphere.fromCenterRadius(c, r), sec)
                    .intersect).toBe(false);
            }
        }
    });

    it('handles a zero-radius disk as a point containment test', () => {
        // A point strictly inside the sector.
        expect(ti.test(disk(1, 0, 0), s).intersect).toBe(true);
        // A point outside the cone.
        expect(ti.test(disk(0, 1, 0), s).intersect).toBe(false);
        // A point beyond the arc.
        expect(ti.test(disk(3, 0, 0), s).intersect).toBe(false);
    });

    it('handles a sector with half angle exactly pi/2 (a half disk)', () => {
        const half = sector(0, 0, 1, 0, Math.PI / 2);
        expect(ti.test(disk(0.5, 0, 0.1), half).intersect).toBe(true);
        expect(ti.test(disk(0, 0.5, 0.1), half).intersect).toBe(true);
        expect(ti.test(disk(-0.5, 0, 0.1), half).intersect).toBe(false);
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, unitVector, wellScaledVector, seededRandom
} from './helpers/arbitraries.js';
import { dot } from '../src/Vector.js';

// Upstream requires the sector to be convex, so the half-angle is in
// (0, pi/2].
const sectorArb = fc.tuple(wellScaledVector(2, -4, 4), unitVector(2),
    fc.double({ min: 0.2, max: 4, noNaN: true }),
    fc.double({ min: 0.05, max: Math.PI / 2, noNaN: true }))
    .map(([v, d, r, a]) =>
        Sector2.fromVertexRadiusDirectionAngle(v, r, d, a));

const diskArb = fc.tuple(wellScaledVector(2, -4, 4),
    fc.double({ min: 0.1, max: 3, noNaN: true }))
    .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));

// Membership in the solid sector, from the definition in Sector2.h.
function inSector(s: Sector2, p: Vector): boolean {
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

// The minimum distance from the disk centre to the sector, by sampling the
// sector densely. This is an upper bound on the true distance.
function sampledDistanceToSector(s: Sector2, c: Vector, nr: number,
    na: number): number {
    let best = Infinity;
    const perp = Vector.fromArray([-s.direction.values[1],
        s.direction.values[0]]);
    for (let i = 0; i <= nr; ++i) {
        const r = (s.radius * i) / nr;
        for (let k = 0; k <= na; ++k) {
            const a = -s.angle + (2 * s.angle * k) / na;
            const p = add(s.vertex, add(mul(r * Math.cos(a), s.direction),
                mul(r * Math.sin(a), perp)));
            best = Math.min(best, length(sub(p, c)));
        }
    }
    return best;
}

describe('IntrDisk2Sector2 verification', () => {
    const tiq = new IntrDisk2Sector2TI();

    it('agrees with a dense sampling of the sector', () => {
        const rnd = seededRandom(0x18ae6d54);
        for (let trial = 0; trial < 200; ++trial) {
            const s = Sector2.fromVertexRadiusDirectionAngle(
                Vector.fromArray([rnd() * 6 - 3, rnd() * 6 - 3]), 0.3 + rnd() * 3,
                (() => {
                    const a = rnd() * 2 * Math.PI;
                    return Vector.fromArray([Math.cos(a), Math.sin(a)]);
                })(), 0.05 + rnd() * (Math.PI / 2 - 0.05));
            const disk = Hypersphere.fromCenterRadius(
                Vector.fromArray([rnd() * 6 - 3, rnd() * 6 - 3]),
                0.2 + rnd() * 1.5);
            const got = tiq.test(disk, s).intersect;
            const d = sampledDistanceToSector(s, disk.center, 60, 90);
            // The sampled distance is an upper bound on the true distance.
            if (d <= disk.radius) {
                expect(got).toBe(true);
            }
            if (d > disk.radius + 0.2) {
                expect(got).toBe(false);
            }
        }
    }, 30000);

    it('a disk containing a sampled sector point always intersects', () => {
        // Strictly interior fractions: a point generated exactly on the arc
        // or on a boundary ray rounds to just outside the sector. The radial
        // fraction is also bounded away from zero, because a subnormal-scale
        // offset from the sector vertex is absorbed by the vertex
        // coordinates and lands the sample on the far side of the wedge.
        check(fc.tuple(sectorArb,
            fc.double({ min: 0.05, max: 0.95, noNaN: true }),
            fc.double({ min: -0.95, max: 0.95, noNaN: true }),
            fc.double({ min: 0.05, max: 2, noNaN: true })),
            ([s, rf, af, radius]) => {
                const perp = Vector.fromArray([-s.direction.values[1],
                    s.direction.values[0]]);
                const r = rf * s.radius;
                const a = af * s.angle;
                const p = add(s.vertex, add(mul(r * Math.cos(a), s.direction),
                    mul(r * Math.sin(a), perp)));
                expect(inSector(s, p)).toBe(true);
                const disk = Hypersphere.fromCenterRadius(p, radius);
                expect(tiq.test(disk, s).intersect).toBe(true);
            });
    });

    it('a disk containing the sector vertex always intersects', () => {
        check(fc.tuple(sectorArb, unitVector(2),
            fc.double({ min: 0, max: 0.9, noNaN: true }),
            fc.double({ min: 0.2, max: 2, noNaN: true })),
            ([s, d, frac, radius]) => {
                const disk = Hypersphere.fromCenterRadius(
                    add(s.vertex, mul(frac * radius, d)), radius);
                expect(tiq.test(disk, s).intersect).toBe(true);
            });
    });

    it('a disk beyond the sector radius never intersects', () => {
        check(fc.tuple(sectorArb, unitVector(2),
            fc.double({ min: 0.1, max: 2, noNaN: true }),
            fc.double({ min: 0.5, max: 3, noNaN: true })),
            ([s, d, radius, extra]) => {
                const far = add(s.vertex,
                    mul(s.radius + radius + extra, d));
                const disk = Hypersphere.fromCenterRadius(far, radius);
                expect(tiq.test(disk, s).intersect).toBe(false);
            });
    });

    it('is invariant under a common rigid motion', () => {
        check(fc.tuple(diskArb, sectorArb,
            fc.double({ min: -Math.PI, max: Math.PI, noNaN: true }),
            wellScaledVector(2, -4, 4)), ([disk, s, ang, tr]) => {
            const ca = Math.cos(ang), sa = Math.sin(ang);
            const rot = (p: Vector): Vector => Vector.fromArray([
                ca * p.values[0] - sa * p.values[1],
                sa * p.values[0] + ca * p.values[1]]);
            const xf = (p: Vector): Vector => add(tr, rot(p));
            const s2 = Sector2.fromVertexRadiusDirectionAngle(xf(s.vertex),
                s.radius, rot(s.direction), s.angle);
            const d2 = Hypersphere.fromCenterRadius(xf(disk.center),
                disk.radius);
            const a = tiq.test(disk, s).intersect;
            const b = tiq.test(d2, s2).intersect;
            // Skip near-touching configurations, where the decision flips.
            const dist = sampledDistanceToSector(s, disk.center, 24, 36);
            if (Math.abs(dist - disk.radius) < 1e-6) {
                return;
            }
            expect(a).toBe(b);
        });
    });

    it('a zero-radius disk is the point-in-sector test', () => {
        check(fc.tuple(sectorArb, wellScaledVector(2, -5, 5)), ([s, p]) => {
            const disk = Hypersphere.fromCenterRadius(p, 0);
            const expected = inSector(s, p);
            const got = tiq.test(disk, s).intersect;
            const q = sub(p, s.vertex);
            const len = length(q);
            // Skip the boundary band, where the two formulations round
            // differently (the query works with sines and cosines of the
            // half-angle, the reference with the cosine constraint).
            const nearArc = Math.abs(len - s.radius) < 1e-9;
            const nearRay = Math.abs(dot(s.direction, q) - s.cosAngle * len)
                < 1e-9;
            if (nearArc || nearRay) {
                return;
            }
            expect(got).toBe(expected);
        });
    });
});
