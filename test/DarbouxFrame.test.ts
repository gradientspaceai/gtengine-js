import { describe, it, expect } from 'vitest';
import { DarbouxFrame3 } from '../src/DarbouxFrame';
import { ParametricSurface } from '../src/ParametricSurface';
import { Vector, dot, length, sub } from '../src/Vector';
import { cross } from '../src/Vector3';

// A sphere of radius r centered at the origin,
//   X(u,v) = r*(cos(v)*cos(u), cos(v)*sin(u), sin(v)).
// The outward unit normal is X/r and both principal curvatures are 1/r for
// the sign convention of the Darboux frame.
class Sphere3Surface extends ParametricSurface {
    constructor(private r: number) {
        super(3, 0, 2 * Math.PI, -Math.PI / 2, Math.PI / 2, true);
        this.mConstructed = true;
    }

    override evaluate(u: number, v: number, order: number, jet: Vector[]): void {
        const r = this.r;
        const cu = Math.cos(u), su = Math.sin(u);
        const cv = Math.cos(v), sv = Math.sin(v);
        jet[0] = Vector.fromArray([r * cv * cu, r * cv * su, r * sv]);
        if (order >= 1) {
            jet[1] = Vector.fromArray([-r * cv * su, r * cv * cu, 0]);
            jet[2] = Vector.fromArray([-r * sv * cu, -r * sv * su, r * cv]);
        }
        if (order >= 2) {
            jet[3] = Vector.fromArray([-r * cv * cu, -r * cv * su, 0]);
            jet[4] = Vector.fromArray([r * sv * su, -r * sv * cu, 0]);
            jet[5] = Vector.fromArray([-r * cv * cu, -r * cv * su, -r * sv]);
        }
    }
}

// A circular cylinder of radius r about the z-axis,
//   X(u,v) = (r*cos(u), r*sin(u), v).
// The principal curvatures are 0 (along the axis) and 1/r.
class Cylinder3Surface extends ParametricSurface {
    constructor(private r: number) {
        super(3, 0, 2 * Math.PI, -1, 1, true);
        this.mConstructed = true;
    }

    override evaluate(u: number, v: number, order: number, jet: Vector[]): void {
        const r = this.r;
        const cu = Math.cos(u), su = Math.sin(u);
        jet[0] = Vector.fromArray([r * cu, r * su, v]);
        if (order >= 1) {
            jet[1] = Vector.fromArray([-r * su, r * cu, 0]);
            jet[2] = Vector.fromArray([0, 0, 1]);
        }
        if (order >= 2) {
            jet[3] = Vector.fromArray([-r * cu, -r * su, 0]);
            jet[4] = Vector.fromArray([0, 0, 0]);
            jet[5] = Vector.fromArray([0, 0, 0]);
        }
    }
}

// A plane through the origin: X(u,v) = u*U + v*V with {U,V} orthonormal.
class Plane3Surface extends ParametricSurface {
    constructor() {
        super(3, -1, 1, -1, 1, true);
        this.mConstructed = true;
    }

    override evaluate(u: number, v: number, order: number, jet: Vector[]): void {
        jet[0] = Vector.fromArray([u, v, 0]);
        if (order >= 1) {
            jet[1] = Vector.fromArray([1, 0, 0]);
            jet[2] = Vector.fromArray([0, 1, 0]);
        }
        if (order >= 2) {
            jet[3] = Vector.fromArray([0, 0, 0]);
            jet[4] = Vector.fromArray([0, 0, 0]);
            jet[5] = Vector.fromArray([0, 0, 0]);
        }
    }
}

function expectOrthonormalRightHanded(t0: Vector, t1: Vector, n: Vector): void {
    expect(length(t0)).toBeCloseTo(1, 12);
    expect(length(t1)).toBeCloseTo(1, 12);
    expect(length(n)).toBeCloseTo(1, 12);
    expect(dot(t0, t1)).toBeCloseTo(0, 12);
    expect(dot(t0, n)).toBeCloseTo(0, 12);
    expect(dot(t1, n)).toBeCloseTo(0, 12);
    // Right-handed: Cross(t0, t1) = n.
    const c = cross(t0, t1);
    expect(length(sub(c, n))).toBeLessThan(1e-12);
}

describe('DarbouxFrame3', () => {
    it('rejects a surface whose dimension is not 3', () => {
        class Surface2 extends ParametricSurface {
            constructor() {
                super(2, 0, 1, 0, 1, true);
            }
            override evaluate(): void {
            }
        }
        expect(() => new DarbouxFrame3(new Surface2())).toThrow();
    });

    it('exposes the surface it was constructed with', () => {
        const surface = new Sphere3Surface(2);
        const frame = new DarbouxFrame3(surface);
        expect(frame.getSurface()).toBe(surface);
    });

    it('computes an orthonormal frame on a sphere with normal along X/r', () => {
        const r = 2.5;
        const frame = new DarbouxFrame3(new Sphere3Surface(r));
        const uv: Array<[number, number]> = [
            [0, 0], [0.3, 0.7], [1.9, -0.4], [4.2, 1.1], [5.9, -1.2]
        ];
        for (const [u, v] of uv) {
            const result = frame.compute(u, v);
            expectOrthonormalRightHanded(result.tangent0, result.tangent1,
                result.normal);
            // The outward unit normal of the sphere is position/r.
            for (let i = 0; i < 3; ++i) {
                expect(result.normal.get(i)).toBeCloseTo(
                    result.position.get(i) / r, 12);
            }
            expect(length(result.position)).toBeCloseTo(r, 12);
            // tangent0 is the normalized dX/du.
            expect(dot(result.tangent0,
                Vector.fromArray([-Math.sin(u), Math.cos(u), 0]))).toBeCloseTo(1, 12);
        }
    });

    it('computes the principal curvatures of a sphere as 1/r (umbilic)', () => {
        for (const r of [0.5, 1, 3.25]) {
            const frame = new DarbouxFrame3(new Sphere3Surface(r));
            for (const [u, v] of [[0.2, 0.1], [2.0, -0.9], [5.0, 1.3]]) {
                const info = frame.getPrincipalInformation(u, v);
                // The characteristic polynomial has a double root here, so
                // sqrt(max(c1^2 - 4*c0*c2, 0)) amplifies the cancellation
                // error; a tolerance near sqrt(eps) is the best available.
                expect(info.curvature0).toBeCloseTo(1 / r, 7);
                expect(info.curvature1).toBeCloseTo(1 / r, 7);
                // The directions are orthonormal and tangent to the sphere.
                expect(length(info.direction0)).toBeCloseTo(1, 12);
                expect(length(info.direction1)).toBeCloseTo(1, 12);
                expect(dot(info.direction0, info.direction1)).toBeCloseTo(0, 12);
                const position = frame.getSurface().getPosition(u, v);
                expect(dot(info.direction0, position)).toBeCloseTo(0, 10);
                expect(dot(info.direction1, position)).toBeCloseTo(0, 10);
            }
        }
    });

    it('computes the principal information of a cylinder', () => {
        const r = 1.75;
        const frame = new DarbouxFrame3(new Cylinder3Surface(r));
        for (const u of [0, 0.7, 2.4, 5.5]) {
            const info = frame.getPrincipalInformation(u, 0.3);
            // curvature0 is the smaller root; along the axis it is zero.
            expect(info.curvature0).toBeCloseTo(0, 12);
            expect(info.curvature1).toBeCloseTo(1 / r, 12);
            // The zero-curvature direction is parallel to the axis.
            expect(Math.abs(info.direction0.get(2))).toBeCloseTo(1, 12);
            expect(info.direction0.get(0)).toBeCloseTo(0, 12);
            expect(info.direction0.get(1)).toBeCloseTo(0, 12);
            // The other direction is in the circular cross-section.
            expect(info.direction1.get(2)).toBeCloseTo(0, 12);
            expect(length(info.direction1)).toBeCloseTo(1, 12);
            expect(dot(info.direction0, info.direction1)).toBeCloseTo(0, 12);
        }
    });

    it('reports zero curvatures on a plane and falls back to dX/du', () => {
        const frame = new DarbouxFrame3(new Plane3Surface());
        const info = frame.getPrincipalInformation(0.2, -0.4);
        expect(info.curvature0).toBeCloseTo(0, 14);
        expect(info.curvature1).toBeCloseTo(0, 14);
        // Both a0 and a1 vanish, so direction0 degenerates to dX/du.
        expect(info.direction0.get(0)).toBeCloseTo(1, 14);
        expect(info.direction0.get(1)).toBeCloseTo(0, 14);
        expect(info.direction1.get(1)).toBeCloseTo(-1, 14);
    });

    it('has mean and Gaussian curvature consistent with a torus', () => {
        // A torus with tube radius a about a circle of radius b. The
        // principal curvatures are 1/a and cos(v)/(b + a*cos(v)).
        const a = 0.75, b = 2;
        class Torus3Surface extends ParametricSurface {
            constructor() {
                super(3, 0, 2 * Math.PI, 0, 2 * Math.PI, true);
                this.mConstructed = true;
            }
            override evaluate(u: number, v: number, order: number,
                jet: Vector[]): void {
                const cu = Math.cos(u), su = Math.sin(u);
                const cv = Math.cos(v), sv = Math.sin(v);
                const rr = b + a * cv;
                jet[0] = Vector.fromArray([rr * cu, rr * su, a * sv]);
                if (order >= 1) {
                    jet[1] = Vector.fromArray([-rr * su, rr * cu, 0]);
                    jet[2] = Vector.fromArray([-a * sv * cu, -a * sv * su, a * cv]);
                }
                if (order >= 2) {
                    jet[3] = Vector.fromArray([-rr * cu, -rr * su, 0]);
                    jet[4] = Vector.fromArray([a * sv * su, -a * sv * cu, 0]);
                    jet[5] = Vector.fromArray([-a * cv * cu, -a * cv * su, -a * sv]);
                }
            }
        }
        const frame = new DarbouxFrame3(new Torus3Surface());
        for (const [u, v] of [[0.4, 0.6], [2.2, 3.1], [5.0, 1.0]]) {
            const info = frame.getPrincipalInformation(u, v);
            const k1 = 1 / a;
            const k2 = Math.cos(v) / (b + a * Math.cos(v));
            const expected = [Math.min(k1, k2), Math.max(k1, k2)];
            const actual = [
                Math.min(info.curvature0, info.curvature1),
                Math.max(info.curvature0, info.curvature1)
            ];
            expect(actual[0]).toBeCloseTo(expected[0], 10);
            expect(actual[1]).toBeCloseTo(expected[1], 10);
        }
    });

    it('has a frame whose normal matches the principal-direction normal', () => {
        const frame = new DarbouxFrame3(new Sphere3Surface(1.4));
        const u = 1.23, v = 0.45;
        const result = frame.compute(u, v);
        const info = frame.getPrincipalInformation(u, v);
        // direction1 = Cross(direction0, normal), so Cross(direction1,
        // direction0) recovers the surface normal.
        const n = cross(info.direction1, info.direction0);
        expect(length(sub(n, result.normal))).toBeLessThan(1e-12);
    });
});
