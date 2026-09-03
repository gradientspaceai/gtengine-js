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
