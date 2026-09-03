import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Cone } from '../src/Cone.js';
import { IntrAlignedBox3Cone3TI } from '../src/IntrAlignedBox3Cone3.js';
import {
    IntrOrientedBox3Cone3TI,
    defaultIntrOrientedBox3Cone3TIResult
} from '../src/IntrOrientedBox3Cone3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function cone(origin: number[], direction: number[], angle: number,
    minHeight: number, maxHeight: number): Cone {
    const C = new Cone(3);
    const d = Vector.fromArray(direction);
    normalize(d);
    C.ray = Ray.fromOriginDirection(Vector.fromArray(origin), d);
    C.setAngle(angle);
    if (maxHeight < 0) {
        if (minHeight > 0) {
            C.makeInfiniteTruncatedCone(minHeight);
        }
        else {
            C.makeInfiniteCone();
        }
    }
    else if (minHeight > 0) {
        C.makeConeFrustum(minHeight, maxHeight);
    }
    else {
        C.makeFiniteCone(maxHeight);
    }
    return C;
}

// An oriented box whose axes are the standard basis (so it coincides with the
// aligned box [center-extent, center+extent]).
function axisAlignedOBB(center: number[], extent: number[]): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(Vector.fromArray(center),
        [Vector.unit(3, 0), Vector.unit(3, 1), Vector.unit(3, 2)],
        Vector.fromArray(extent));
}

// A right-handed orthonormal frame from a rotation about the given unit axis.
function rotationFrame(axis: Vector, angle: number): Vector[] {
    const u = axis.clone();
    normalize(u);
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    const basis = [Vector.unit(3, 0), Vector.unit(3, 1), Vector.unit(3, 2)];
    return basis.map((e) => {
        // Rodrigues' rotation formula.
        const term0 = mul(c, e);
        const term1 = mul(s, cross(u, e));
        const term2 = mul((1 - c) * dot(u, e), u);
        const r = add(add(term0, term1), term2);
        return r;
    });
}

// The box in world coordinates given as a rigid transform of the OBB frame.
function boxPoint(box: OrientedBox, s: number[]): Vector {
    let p = box.center.clone();
    for (let i = 0; i < 3; ++i) {
        p = add(p, mul(s[i] * box.extent.values[i], box.axis[i]));
    }
    return p;
}

const ti = new IntrOrientedBox3Cone3TI();
const abTI = new IntrAlignedBox3Cone3TI();

const quarterPi = Math.PI / 4;

describe('IntrOrientedBox3Cone3', () => {
    it('default-constructs the result as no intersection', () => {
        expect(defaultIntrOrientedBox3Cone3TIResult()).toEqual({
            intersect: false
        });
    });

    it('matches the aligned-box query when the axes are the standard basis', () => {
        const cones = [
            cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1),
            cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4),
            cone([0, 0, 0], [0, 0, 1], quarterPi, 1, 3),
            cone([1, 2, -1], [1, 1, 1], 0.5, 0.5, 6)
        ];
        const boxes: [number[], number[]][] = [
            [[0, 0, 2], [0.5, 0.5, 0.5]],
            [[5, 5, 5], [0.5, 0.5, 0.5]],
            [[0, 0, -3], [1, 1, 1]],
            [[3, 0, 2], [1, 1, 1]],
            [[0, 0, 8], [2, 2, 2]],
            [[1, 1, 1], [3, 3, 3]]
        ];
        for (const C of cones) {
            for (const [center, extent] of boxes) {
                const obb = axisAlignedOBB(center, extent);
                const ab = AlignedBox.fromMinMax(
                    sub(Vector.fromArray(center), Vector.fromArray(extent)),
                    add(Vector.fromArray(center), Vector.fromArray(extent)));
                expect(ti.test(obb, C).intersect)
                    .toBe(abTI.test(ab, C).intersect);
            }
        }
    });

    it('detects a box straddling the cone axis', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 10);
        const box = OrientedBox.fromCenterAxisExtent(vec(0, 0, 5),
            rotationFrame(vec(1, 1, 1), 0.7), vec(0.5, 0.5, 0.5));
        expect(ti.test(box, C).intersect).toBe(true);
    });

    it('rejects a box far outside the cone', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 10);
        const box = OrientedBox.fromCenterAxisExtent(vec(50, 0, 5),
            rotationFrame(vec(0, 1, 0), 0.3), vec(0.5, 0.5, 0.5));
        expect(ti.test(box, C).intersect).toBe(false);
    });

    it('rejects a box below the cone vertex that only touches it', () => {
        // Upstream reports an intersection only for positive volume, so a box
        // whose top face touches the cone vertex is not an intersection.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 10);
        const touching = axisAlignedOBB([0, 0, -1], [1, 1, 1]);
        expect(ti.test(touching, C).intersect).toBe(false);
        const overlapping = axisAlignedOBB([0, 0, -0.9], [1, 1, 1]);
        expect(ti.test(overlapping, C).intersect).toBe(true);
    });

    it('respects the cone height range for a frustum', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 4, 6);
        // A small box on the axis below minHeight.
        expect(ti.test(axisAlignedOBB([0, 0, 2], [0.5, 0.5, 0.5]), C).intersect)
            .toBe(false);
        // On the axis inside the height range.
        expect(ti.test(axisAlignedOBB([0, 0, 5], [0.5, 0.5, 0.5]), C).intersect)
            .toBe(true);
        // On the axis above maxHeight.
        expect(ti.test(axisAlignedOBB([0, 0, 8], [0.5, 0.5, 0.5]), C).intersect)
            .toBe(false);
    });

    it('is invariant under a rigid motion of both box and cone', () => {
        // Rotating and translating the world must not change the answer.
        const axes = rotationFrame(vec(1, -2, 3), 1.1);
        const translation = vec(-2, 5, 1);
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

        let seed = 4242424;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (a: number, b: number): number => a + (b - a) * rand();

        let numHits = 0;
        for (let trial = 0; trial < 200; ++trial) {
            const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 6);
            const boxAxes = rotationFrame(
                vec(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1) + 0.01), rnd(0, 3));
            const center = vec(rnd(-5, 5), rnd(-5, 5), rnd(-2, 8));
            const extent = vec(rnd(0.2, 1.5), rnd(0.2, 1.5), rnd(0.2, 1.5));
            const box = OrientedBox.fromCenterAxisExtent(center, boxAxes,
                extent);
            const before = ti.test(box, C).intersect;
            if (before) {
                ++numHits;
            }

            const movedCone = new Cone(3);
            movedCone.ray = Ray.fromOriginDirection(mapPoint(C.ray.origin),
                mapDirection(C.ray.direction));
            movedCone.setAngle(C.angle);
            movedCone.makeFiniteCone(C.getMaxHeight());
            const movedBox = OrientedBox.fromCenterAxisExtent(mapPoint(center),
                boxAxes.map(mapDirection), extent);
            expect(ti.test(movedBox, movedCone).intersect).toBe(before);
        }
        expect(numHits).toBeGreaterThan(20);
    });

    it('agrees with a dense point sampling of the box (randomized)', () => {
        // Sample points of the box; if a sampled interior point is strictly
        // inside the solid cone, the query must report an intersection.
        let seed = 777001;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (a: number, b: number): number => a + (b - a) * rand();

        const C = cone([0, 0, 0], [0, 0, 1], 0.6, 0.5, 5);
        const strictlyInside = (X: Vector): boolean => {
            const diff = sub(X, C.ray.origin);
            const h = dot(C.ray.direction, diff);
            if (h <= C.getMinHeight() + 1e-9 || h >= C.getMaxHeight() - 1e-9) {
                return false;
            }
            return h > length(diff) * C.cosAngle + 1e-9;
        };

        let numDetected = 0;
        for (let trial = 0; trial < 150; ++trial) {
            const boxAxes = rotationFrame(
                vec(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1) + 0.01), rnd(0, 3));
            const box = OrientedBox.fromCenterAxisExtent(
                vec(rnd(-4, 4), rnd(-4, 4), rnd(-1, 7)), boxAxes,
                vec(rnd(0.2, 1.2), rnd(0.2, 1.2), rnd(0.2, 1.2)));

            let sawInside = false;
            const n = 4;
            for (let i = 0; i <= n && !sawInside; ++i) {
                for (let j = 0; j <= n && !sawInside; ++j) {
                    for (let k = 0; k <= n && !sawInside; ++k) {
                        const s = [
                            -1 + (2 * i) / n, -1 + (2 * j) / n, -1 + (2 * k) / n
                        ];
                        if (strictlyInside(boxPoint(box, s))) {
                            sawInside = true;
                        }
                    }
                }
            }

            const reported = ti.test(box, C).intersect;
            if (sawInside) {
                expect(reported).toBe(true);
                ++numDetected;
            }
        }
        expect(numDetected).toBeGreaterThan(10);
    });
});
