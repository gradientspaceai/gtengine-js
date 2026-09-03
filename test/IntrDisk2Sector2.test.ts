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
