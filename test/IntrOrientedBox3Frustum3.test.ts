import { describe, it, expect } from 'vitest';
import { Frustum3 } from '../src/Frustum3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';
import { IntrOrientedBox3Frustum3TI } from '../src/IntrOrientedBox3Frustum3.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// An axis-aligned oriented box.
function alignedBox(center: number[], extent: number[]): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(vec(center),
        [vec([1, 0, 0]), vec([0, 1, 0]), vec([0, 0, 1])], vec(extent));
}

// An oriented box whose first axis is the given direction.
function tiltedBox(center: number[], axis0: number[], extent: number[]):
    OrientedBox {
    const a0 = vec(axis0);
    normalize(a0);
    const basis = [a0, new Vector(3), new Vector(3)];
    computeOrthogonalComplement3(1, basis);
    return OrientedBox.fromCenterAxisExtent(vec(center), basis, vec(extent));
}

// Independent point-in-frustum test.
function insideFrustum(f: Frustum3, x: Vector): boolean {
    const diff = sub(x, f.origin);
    const d = dot(diff, f.dVector);
    if (d < f.dMin || d > f.dMax) {
        return false;
    }
    const scale = d / f.dMin;
    return Math.abs(dot(diff, f.rVector)) <= f.rBound * scale
        && Math.abs(dot(diff, f.uVector)) <= f.uBound * scale;
}

describe('IntrOrientedBox3Frustum3', () => {
    const ti = new IntrOrientedBox3Frustum3TI();

    // The default frustum: origin (0,0,0), view direction (0,0,1),
    // dMin = 1, dMax = 2, uBound = rBound = 1.
    const frustum = new Frustum3();

    it('accepts a small box inside the frustum', () => {
        const b = alignedBox([0, 0, 1.5], [0.1, 0.1, 0.1]);
        expect(ti.test(b, frustum).intersect).toBe(true);
    });

    it('rejects a box behind the eye point', () => {
        const b = alignedBox([0, 0, -5], [1, 1, 1]);
        expect(ti.test(b, frustum).intersect).toBe(false);
    });

    it('rejects a box beyond the far plane', () => {
        const b = alignedBox([0, 0, 10], [1, 1, 1]);
        expect(ti.test(b, frustum).intersect).toBe(false);
    });

    it('rejects a box entirely to the right of the frustum', () => {
        const b = alignedBox([10, 0, 1.5], [1, 1, 1]);
        expect(ti.test(b, frustum).intersect).toBe(false);
    });

    it('rejects a box entirely above the frustum', () => {
        const b = alignedBox([0, 10, 1.5], [1, 1, 1]);
        expect(ti.test(b, frustum).intersect).toBe(false);
    });

    it('accepts a box straddling the near plane', () => {
        const b = alignedBox([0, 0, 1], [0.5, 0.5, 0.5]);
        expect(ti.test(b, frustum).intersect).toBe(true);
    });

    it('accepts a large box that encloses the frustum', () => {
        const b = alignedBox([0, 0, 0], [20, 20, 20]);
        expect(ti.test(b, frustum).intersect).toBe(true);
    });

    it('accepts a tilted box poking into the frustum side', () => {
        const b = tiltedBox([1.6, 0, 1.5], [1, 0, 1], [0.6, 0.2, 0.2]);
        expect(ti.test(b, frustum).intersect).toBe(true);
    });

    it('handles a degenerate box with zero extents as a point', () => {
        const insidePoint = alignedBox([0, 0, 1.5], [0, 0, 0]);
        expect(ti.test(insidePoint, frustum).intersect).toBe(true);
        const outsidePoint = alignedBox([0, 0, 0.5], [0, 0, 0]);
        expect(ti.test(outsidePoint, frustum).intersect).toBe(false);
    });

    it('reports an intersection whenever a sampled point lies in both', () => {
        const rand = makeRandom(1414213);
        const f = Frustum3.fromParameters(
            vec([0.5, -1, 0.25]), vec([0, 0, 1]), vec([0, 1, 0]),
            vec([1, 0, 0]), 1, 4, 0.8, 1.2);
        let hits = 0;
        for (let trial = 0; trial < 250; ++trial) {
            const b = tiltedBox(
                [6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 1],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1],
                [0.2 + rand(), 0.2 + rand(), 0.2 + rand()]);
            const result = ti.test(b, f).intersect;

            // Sample the box on a grid; any sample inside the frustum forces
            // the query to report an intersection.
            let found = false;
            const n = 12;
            for (let i = 0; i <= n && !found; ++i) {
                const a0 = -b.extent.values[0] + (2 * b.extent.values[0] * i) / n;
                for (let j = 0; j <= n && !found; ++j) {
                    const a1 = -b.extent.values[1]
                        + (2 * b.extent.values[1] * j) / n;
                    for (let k = 0; k <= n; ++k) {
                        const a2 = -b.extent.values[2]
                            + (2 * b.extent.values[2] * k) / n;
                        const p = add(b.center, add(mul(a0, b.axis[0]),
                            add(mul(a1, b.axis[1]), mul(a2, b.axis[2]))));
                        if (insideFrustum(f, p)) {
                            found = true;
                            break;
                        }
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

    it('reports an intersection when a frustum vertex is inside the box', () => {
        const rand = makeRandom(99991);
        const vertices = frustum.computeVertices();
        let hits = 0;
        for (let trial = 0; trial < 200; ++trial) {
            const b = tiltedBox(
                [4 * rand() - 2, 4 * rand() - 2, 4 * rand()],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1],
                [0.1 + rand(), 0.1 + rand(), 0.1 + rand()]);
            let found = false;
            for (const v of vertices) {
                const diff = sub(v, b.center);
                if (Math.abs(dot(diff, b.axis[0])) <= b.extent.values[0]
                    && Math.abs(dot(diff, b.axis[1])) <= b.extent.values[1]
                    && Math.abs(dot(diff, b.axis[2])) <= b.extent.values[2]) {
                    found = true;
                    break;
                }
            }
            if (found) {
                ++hits;
                expect(ti.test(b, frustum).intersect).toBe(true);
            }
        }
        expect(hits).toBeGreaterThan(5);
    });
});
