import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import { OrientedBox } from '../src/OrientedBox';
import { Vector } from '../src/Vector';
import { IntrAlignedBox3OrientedBox3TI } from '../src/IntrAlignedBox3OrientedBox3';

function alignedBox(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(Vector.fromArray(min), Vector.fromArray(max));
}

// A rotation built from Euler angles, applied as Rz*Ry*Rx. The columns are the
// box axes, which are orthonormal by construction.
function rotationAxes(ax: number, ay: number, az: number): Vector[] {
    const cx = Math.cos(ax), sx = Math.sin(ax);
    const cy = Math.cos(ay), sy = Math.sin(ay);
    const cz = Math.cos(az), sz = Math.sin(az);
    const m = [
        [cz * cy, cz * sy * sx - sz * cx, cz * sy * cx + sz * sx],
        [sz * cy, sz * sy * sx + cz * cx, sz * sy * cx - cz * sx],
        [-sy, cy * sx, cy * cx]
    ];
    return [
        Vector.fromArray([m[0][0], m[1][0], m[2][0]]),
        Vector.fromArray([m[0][1], m[1][1], m[2][1]]),
        Vector.fromArray([m[0][2], m[1][2], m[2][2]])
    ];
}

function orientedBox(center: number[], axes: Vector[], extent: number[]): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(Vector.fromArray(center), axes,
        Vector.fromArray(extent));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Support value of the aligned box in direction d.
function supportAligned(box: AlignedBox, d: number[]): number {
    let s = 0;
    for (let i = 0; i < 3; ++i) {
        s += d[i] > 0 ? d[i] * box.max.values[i] : d[i] * box.min.values[i];
    }
    return s;
}

// Support value of the oriented box in direction d.
function supportOriented(box: OrientedBox, d: number[]): number {
    let s = d[0] * box.center.values[0] + d[1] * box.center.values[1] +
        d[2] * box.center.values[2];
    for (let i = 0; i < 3; ++i) {
        const u = box.axis[i].values;
        s += box.extent.values[i] *
            Math.abs(d[0] * u[0] + d[1] * u[1] + d[2] * u[2]);
    }
    return s;
}

describe('IntrAlignedBox3OrientedBox3', () => {
    const ti = new IntrAlignedBox3OrientedBox3TI();
    const identity = rotationAxes(0, 0, 0);

    it('matches aligned-aligned behavior when the oriented box is axis aligned', () => {
        const b0 = alignedBox([0, 0, 0], [2, 2, 2]);
        expect(ti.test(b0, orientedBox([3, 1, 1], identity, [1, 1, 1])).intersect)
            .toBe(true);   // face contact at x = 2
        const r = ti.test(b0, orientedBox([3.5, 1, 1], identity, [1, 1, 1]));
        expect(r.intersect).toBe(false);
        expect(r.separating).toEqual([0, -1]);
    });

    it('reports separation by a face normal of the oriented box', () => {
        // A thin plate rotated 45 degrees about z, offset along the diagonal.
        const b0 = alignedBox([-1, -1, -1], [1, 1, 1]);
        const axes = rotationAxes(0, 0, Math.PI / 4);
        const box1 = orientedBox([2.5, 2.5, 0], axes, [0.1, 10, 10]);
        const r = ti.test(b0, box1);
        expect(r.intersect).toBe(false);
        expect(r.separating).toEqual([-1, 0]);
        // Moving it onto the origin makes the boxes intersect.
        expect(ti.test(b0, orientedBox([0, 0, 0], axes, [0.1, 10, 10])).intersect)
            .toBe(true);
    });

    it('reports separation by an edge-edge cross product', () => {
        // The classic edge-edge case: two long thin boxes whose long axes are
        // skew. Neither face normal separates them but their cross product
        // does. The rotations avoid any parallel face-normal pair.
        const b0 = alignedBox([-4, -0.5, -0.5], [4, 0.5, 0.5]);
        const axes = rotationAxes(0.3, 0.4, Math.PI / 2 + 0.2);
        const box1 = orientedBox([0, 0, 2.2], axes, [4, 0.5, 0.5]);
        const r = ti.test(b0, box1);
        expect(r.intersect).toBe(false);
        expect(r.separating[0]).toBeGreaterThanOrEqual(0);
        expect(r.separating[1]).toBeGreaterThanOrEqual(0);
    });

    it('skips the edge-edge axes when a parallel pair exists (epsilon)', () => {
        // Boxes with a shared axis direction. The parallel-pair short circuit
        // returns intersect = true after the six face-normal tests.
        const b0 = alignedBox([-1, -1, -1], [1, 1, 1]);
        const axes = rotationAxes(0, 0, 0.3);
        expect(ti.test(b0, orientedBox([0.5, 0.5, 0], axes, [1, 1, 1])).intersect)
            .toBe(true);
        // A large epsilon makes every axis pair "nearly parallel", so the
        // edge-edge tests are always skipped; the result can only become
        // "intersect" more often, never less.
        const skew = rotationAxes(0.3, 0.4, Math.PI / 2 + 0.2);
        const box1 = orientedBox([0, 0, 2.2], skew, [4, 0.5, 0.5]);
        const thin = alignedBox([-4, -0.5, -0.5], [4, 0.5, 0.5]);
        expect(ti.test(thin, box1, 1).intersect).toBe(true);
        // A negative epsilon is clamped to zero, giving the default result.
        expect(ti.test(thin, box1, -5).intersect)
            .toBe(ti.test(thin, box1).intersect);
    });

    it('handles degenerate boxes', () => {
        const b0 = alignedBox([0, 0, 0], [2, 2, 2]);
        const axes = rotationAxes(0.5, 0.6, 0.7);
        // A point at a corner of the aligned box.
        expect(ti.test(b0, orientedBox([2, 2, 2], axes, [0, 0, 0])).intersect)
            .toBe(true);
        expect(ti.test(b0, orientedBox([2, 2, 2.1], axes, [0, 0, 0])).intersect)
            .toBe(false);
        // A zero-extent aligned box (a point) inside the oriented box.
        const point = alignedBox([1, 1, 1], [1, 1, 1]);
        expect(ti.test(point, orientedBox([1, 1, 1], axes, [1, 1, 1])).intersect)
            .toBe(true);
    });

    it('agrees with sampling and separating-direction oracles', () => {
        const rand = makeRandom(987654321);
        let numIntersect = 0, numSeparate = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const lo: number[] = [], hi: number[] = [];
            for (let d = 0; d < 3; ++d) {
                const a = 3 * rand() - 1.5;
                lo.push(a);
                hi.push(a + 0.3 + 1.5 * rand());
            }
            const b0 = alignedBox(lo, hi);
            const axes = rotationAxes(2 * Math.PI * rand(), 2 * Math.PI * rand(),
                2 * Math.PI * rand());
            const b1 = orientedBox(
                [3 * rand() - 1.5, 3 * rand() - 1.5, 3 * rand() - 1.5], axes,
                [0.2 + rand(), 0.2 + rand(), 0.2 + rand()]);

            const intersect = ti.test(b0, b1).intersect;

            // Oracle 1 (sound proof of intersection): sample a grid of points
            // in the oriented box and test containment in the aligned box, and
            // vice versa.
            let foundCommonPoint = false;
            const n = 9;
            for (let i = 0; i <= n && !foundCommonPoint; ++i) {
                for (let j = 0; j <= n && !foundCommonPoint; ++j) {
                    for (let k = 0; k <= n && !foundCommonPoint; ++k) {
                        const a = [-1 + (2 * i) / n, -1 + (2 * j) / n,
                            -1 + (2 * k) / n];
                        const p = [0, 0, 0];
                        for (let d = 0; d < 3; ++d) {
                            p[d] = b1.center.values[d];
                            for (let m = 0; m < 3; ++m) {
                                p[d] += a[m] * b1.extent.values[m] *
                                    b1.axis[m].values[d];
                            }
                        }
                        let inside = true;
                        for (let d = 0; d < 3; ++d) {
                            if (p[d] < b0.min.values[d] || p[d] > b0.max.values[d]) {
                                inside = false;
                            }
                        }
                        if (inside) {
                            foundCommonPoint = true;
                        }
                    }
                }
            }
            if (foundCommonPoint) {
                expect(intersect).toBe(true);
            }

            // Oracle 2 (sound proof of separation): search random directions
            // for one whose projection intervals are disjoint.
            let foundSeparatingDirection = false;
            for (let s = 0; s < 400 && !foundSeparatingDirection; ++s) {
                const d = [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1];
                const len = Math.hypot(d[0], d[1], d[2]);
                if (len < 1e-6) {
                    continue;
                }
                d[0] /= len; d[1] /= len; d[2] /= len;
                const neg = [-d[0], -d[1], -d[2]];
                if (supportAligned(b0, d) < -supportOriented(b1, neg) ||
                    supportOriented(b1, d) < -supportAligned(b0, neg)) {
                    foundSeparatingDirection = true;
                }
            }
            if (foundSeparatingDirection) {
                expect(intersect).toBe(false);
            }

            if (intersect) {
                ++numIntersect;
            } else {
                ++numSeparate;
            }
        }
        expect(numIntersect).toBeGreaterThan(20);
        expect(numSeparate).toBeGreaterThan(20);
    });
});
