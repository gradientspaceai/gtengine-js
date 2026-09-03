import { describe, expect, it } from 'vitest';
import { ConformalMapGenus0 } from '../src/ConformalMapGenus0.js';
import { Vector, dot, length, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

interface Mesh {
    positions: Vector[];
    indices: number[];
}

// An octahedron subdivided 'level' times by 1-to-4 midpoint splits. When
// 'project' is true, the vertices are pushed onto the unit sphere, so the
// mesh is a genus-0 surface that converges to the sphere as the level grows.
function makeSphereMesh(level: number, project: boolean): Mesh {
    let positions: Vector[] = [
        Vector.fromArray([1, 0, 0]), Vector.fromArray([-1, 0, 0]),
        Vector.fromArray([0, 1, 0]), Vector.fromArray([0, -1, 0]),
        Vector.fromArray([0, 0, 1]), Vector.fromArray([0, 0, -1])
    ];
    let indices: number[] = [
        0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4,
        2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5
    ];
    for (let s = 0; s < level; ++s) {
        const mid = new Map<string, number>();
        const newIndices: number[] = [];
        const midpoint = (a: number, b: number): number => {
            const key = a < b ? `${a},${b}` : `${b},${a}`;
            let m = mid.get(key);
            if (m === undefined) {
                const p = new Vector(3);
                for (let k = 0; k < 3; ++k) {
                    p.values[k] = 0.5 * (positions[a].values[k] + positions[b].values[k]);
                }
                m = positions.length;
                positions.push(p);
                mid.set(key, m);
            }
            return m;
        };
        for (let t = 0; t < indices.length / 3; ++t) {
            const i0 = indices[3 * t], i1 = indices[3 * t + 1], i2 = indices[3 * t + 2];
            const m01 = midpoint(i0, i1);
            const m12 = midpoint(i1, i2);
            const m20 = midpoint(i2, i0);
            newIndices.push(i0, m01, m20, m01, i1, m12, m12, i2, m20, m01, m12, m20);
        }
        indices = newIndices;
    }
    if (project) {
        positions = positions.map(p => {
            const q = p.clone();
            normalize(q);
            return q;
        });
    }
    return { positions, indices };
}

// The solid angle subtended at the origin by the spherical triangle (a,b,c)
// of unit vectors (Van Oosterom and Strackee). The value is signed: positive
// when (a,b,c) is counterclockwise as seen from outside the sphere.
function solidAngle(a: Vector, b: Vector, c: Vector): number {
    const numer = dot(a, cross(b, c));
    const denom = 1 + dot(a, b) + dot(b, c) + dot(c, a);
    return 2 * Math.atan2(numer, denom);
}

// The cotangent-Laplacian weight sum at each vertex, applied to a scalar
// field. This is an independent reconstruction of the linear operator the
// algorithm solves with, used to verify that the plane coordinates are
// discretely harmonic away from the punctured triangle.
function applyCotangentLaplacian(mesh: Mesh, u: number[]): number[] {
    const n = mesh.positions.length;
    const result = new Array<number>(n).fill(0);
    const weights = new Map<string, number>();
    for (let t = 0; t < mesh.indices.length / 3; ++t) {
        const v = [mesh.indices[3 * t], mesh.indices[3 * t + 1], mesh.indices[3 * t + 2]];
        for (let i = 0; i < 3; ++i) {
            const v0 = v[i], v1 = v[(i + 1) % 3], v2 = v[(i + 2) % 3];
            const E0 = sub(mesh.positions[v0], mesh.positions[v2]);
            const E1 = sub(mesh.positions[v1], mesh.positions[v2]);
            const cot = dot(E0, E1) / length(cross(E0, E1));
            const key = v0 < v1 ? `${v0},${v1}` : `${v1},${v0}`;
            weights.set(key, (weights.get(key) ?? 0) + cot);
        }
    }
    for (const [key, sum] of weights) {
        const parts = key.split(',');
        const v0 = Number(parts[0]), v1 = Number(parts[1]);
        const w = -0.5 * sum;
        result[v0] += w * u[v1];
        result[v1] += w * u[v0];
        result[v0] -= w * u[v0];
        result[v1] -= w * u[v1];
    }
    return result;
}

// The median of the per-triangle relative edge-ratio distortion between the
// input mesh and its spherical image. A conformal map preserves shapes of
// small triangles, so this measure decreases as the mesh is refined.
function medianShapeDistortion(mesh: Mesh, sphere: Vector[]): number {
    const errors: number[] = [];
    for (let t = 0; t < mesh.indices.length / 3; ++t) {
        const i0 = mesh.indices[3 * t], i1 = mesh.indices[3 * t + 1];
        const i2 = mesh.indices[3 * t + 2];
        const o0 = length(sub(mesh.positions[i1], mesh.positions[i0]));
        const o1 = length(sub(mesh.positions[i2], mesh.positions[i1]));
        const o2 = length(sub(mesh.positions[i0], mesh.positions[i2]));
        const n0 = length(sub(sphere[i1], sphere[i0]));
        const n1 = length(sub(sphere[i2], sphere[i1]));
        const n2 = length(sub(sphere[i0], sphere[i2]));
        errors.push(Math.abs((o1 / o0) / (n1 / n0) - 1));
        errors.push(Math.abs((o2 / o0) / (n2 / n0) - 1));
    }
    errors.sort((x, y) => x - y);
    return errors[Math.floor(errors.length / 2)];
}

describe('ConformalMapGenus0', () => {
    it('starts empty', () => {
        const map = new ConformalMapGenus0();
        expect(map.getPlaneCoordinates()).toEqual([]);
        expect(map.getSphereCoordinates()).toEqual([]);
        expect(map.getSphereRadius()).toBe(0);
        expect(map.getMinPlaneCoordinate().values).toEqual([0, 0]);
        expect(map.getMaxPlaneCoordinate().values).toEqual([0, 0]);
    });

    it('maps every vertex of a subdivided octahedron onto the unit sphere', () => {
        const mesh = makeSphereMesh(2, false);
        const map = new ConformalMapGenus0();
        const converged = map.compute(mesh.positions, mesh.indices, 0);
        expect(converged).toBe(true);

        const sphere = map.getSphereCoordinates();
        const plane = map.getPlaneCoordinates();
        expect(sphere.length).toBe(mesh.positions.length);
        expect(plane.length).toBe(mesh.positions.length);
        for (const s of sphere) {
            expect(s.size).toBe(3);
            expect(length(s)).toBeCloseTo(1, 12);
        }
        for (const p of plane) {
            expect(p.size).toBe(2);
        }
        expect(map.getSphereRadius()).toBeGreaterThan(0);
    });

    it('preserves the mesh topology: no folded or degenerate triangles', () => {
        for (const level of [1, 2, 3]) {
            const mesh = makeSphereMesh(level, true);
            const map = new ConformalMapGenus0();
            map.compute(mesh.positions, mesh.indices, 0);
            const sphere = map.getSphereCoordinates();

            // Every image triangle keeps the input winding, so the signed
            // volume of (origin, a, b, c) is positive for all of them.
            let total = 0;
            for (let t = 0; t < mesh.indices.length / 3; ++t) {
                const a = sphere[mesh.indices[3 * t]];
                const b = sphere[mesh.indices[3 * t + 1]];
                const c = sphere[mesh.indices[3 * t + 2]];
                expect(dot(a, cross(sub(b, a), sub(c, a)))).toBeGreaterThan(0);
                total += solidAngle(a, b, c);
            }

            // The images tile the sphere exactly once, so the solid angles
            // sum to the full 4*pi. This holds only for a bijective map.
            expect(total).toBeCloseTo(4 * Math.PI, 8);
        }
    });

    it('produces distinct images for distinct vertices', () => {
        const mesh = makeSphereMesh(2, true);
        const map = new ConformalMapGenus0();
        map.compute(mesh.positions, mesh.indices, 0);
        const sphere = map.getSphereCoordinates();
        let minSeparation = Infinity;
        for (let i = 0; i < sphere.length; ++i) {
            for (let j = i + 1; j < sphere.length; ++j) {
                minSeparation = Math.min(minSeparation,
                    length(sub(sphere[i], sphere[j])));
            }
        }
        expect(minSeparation).toBeGreaterThan(0);
    });

    it('approximately preserves angles, with the error decreasing under refinement', () => {
        const distortions: number[] = [];
        for (const level of [1, 2, 3]) {
            const mesh = makeSphereMesh(level, true);
            const map = new ConformalMapGenus0();
            map.compute(mesh.positions, mesh.indices, 0);
            distortions.push(medianShapeDistortion(mesh, map.getSphereCoordinates()));
        }
        // Refinement halves the triangle size, so the first-order shape
        // distortion of the (Moebius) map decreases monotonically.
        expect(distortions[1]).toBeLessThan(distortions[0]);
        expect(distortions[2]).toBeLessThan(distortions[1]);
        expect(distortions[2]).toBeLessThan(0.1);
    });

    it('produces plane coordinates that are discretely harmonic away from the puncture', () => {
        const mesh = makeSphereMesh(2, true);
        const punctureTriangle = 3;
        const map = new ConformalMapGenus0();
        map.compute(mesh.positions, mesh.indices, punctureTriangle);
        const plane = map.getPlaneCoordinates();

        const punctured = new Set<number>([
            mesh.indices[3 * punctureTriangle],
            mesh.indices[3 * punctureTriangle + 1],
            mesh.indices[3 * punctureTriangle + 2]
        ]);
        for (let k = 0; k < 2; ++k) {
            const u = plane.map(p => p.values[k]);
            const scale = Math.max(...u.map(Math.abs));
            const Au = applyCotangentLaplacian(mesh, u);
            for (let i = 0; i < u.length; ++i) {
                if (!punctured.has(i)) {
                    expect(Math.abs(Au[i]) / scale).toBeLessThan(1e-4);
                }
            }
        }
    });

    it('centers the plane coordinates and bounds them by the reported extremes', () => {
        const mesh = makeSphereMesh(2, true);
        const map = new ConformalMapGenus0();
        map.compute(mesh.positions, mesh.indices, 0);
        const plane = map.getPlaneCoordinates();
        const pmin = map.getMinPlaneCoordinate();
        const pmax = map.getMaxPlaneCoordinate();

        let sumX = 0, sumY = 0;
        for (const p of plane) {
            expect(p.values[0]).toBeGreaterThanOrEqual(pmin.values[0]);
            expect(p.values[0]).toBeLessThanOrEqual(pmax.values[0]);
            expect(p.values[1]).toBeGreaterThanOrEqual(pmin.values[1]);
            expect(p.values[1]).toBeLessThanOrEqual(pmax.values[1]);
            sumX += p.values[0];
            sumY += p.values[1];
        }
        // The average is subtracted before the extremes are computed.
        expect(sumX / plane.length).toBeCloseTo(0, 10);
        expect(sumY / plane.length).toBeCloseTo(0, 10);
    });

    it('depends on the choice of puncture triangle only up to a Moebius transformation', () => {
        // The mapped meshes differ, but each is a valid bijection onto the
        // sphere and each preserves the mesh topology.
        const mesh = makeSphereMesh(2, true);
        const spheres: Vector[][] = [];
        for (const punctureTriangle of [0, 5, 100]) {
            const map = new ConformalMapGenus0();
            expect(map.compute(mesh.positions, mesh.indices, punctureTriangle))
                .toBe(true);
            const sphere = map.getSphereCoordinates();
            let total = 0;
            for (let t = 0; t < mesh.indices.length / 3; ++t) {
                total += solidAngle(sphere[mesh.indices[3 * t]],
                    sphere[mesh.indices[3 * t + 1]],
                    sphere[mesh.indices[3 * t + 2]]);
            }
            expect(total).toBeCloseTo(4 * Math.PI, 8);
            spheres.push(sphere.map(s => s.clone()));
        }
        // The maps are genuinely different (a different puncture triangle
        // gives a different Moebius normalization).
        let maxDifference = 0;
        for (let i = 0; i < spheres[0].length; ++i) {
            maxDifference = Math.max(maxDifference,
                length(sub(spheres[0][i], spheres[1][i])));
        }
        expect(maxDifference).toBeGreaterThan(1e-3);
    });

    it('is insensitive to a rigid motion of the input mesh', () => {
        const mesh = makeSphereMesh(2, true);
        // A rotation by 90 degrees about the z axis followed by a
        // translation. The cotangent weights and the puncture-triangle frame
        // are rigid-motion invariant, so the plane coordinates match.
        const moved: Vector[] = mesh.positions.map(p => Vector.fromArray([
            -p.values[1] + 3, p.values[0] - 7, p.values[2] + 11
        ]));

        const map0 = new ConformalMapGenus0();
        map0.compute(mesh.positions, mesh.indices, 0);
        const map1 = new ConformalMapGenus0();
        map1.compute(moved, mesh.indices, 0);

        expect(map1.getSphereRadius()).toBeCloseTo(map0.getSphereRadius(), 8);
        const s0 = map0.getSphereCoordinates();
        const s1 = map1.getSphereCoordinates();
        for (let i = 0; i < s0.length; ++i) {
            expect(length(sub(s0[i], s1[i]))).toBeLessThan(1e-8);
        }
    });

    it('is invariant, up to scale, when the input mesh is uniformly scaled', () => {
        const mesh = makeSphereMesh(2, true);
        const scaled = mesh.positions.map(p => Vector.fromArray(
            p.values.map(x => 7 * x)));

        const map0 = new ConformalMapGenus0();
        map0.compute(mesh.positions, mesh.indices, 0);
        const map1 = new ConformalMapGenus0();
        map1.compute(scaled, mesh.indices, 0);

        const s0 = map0.getSphereCoordinates();
        const s1 = map1.getSphereCoordinates();
        for (let i = 0; i < s0.length; ++i) {
            expect(length(sub(s0[i], s1[i]))).toBeLessThan(1e-8);
        }
    });

    it('rejects invalid inputs', () => {
        const mesh = makeSphereMesh(1, false);
        const map = new ConformalMapGenus0();
        expect(() => map.compute([], mesh.indices, 0))
            .toThrow('The mesh must have positions.');
        expect(() => map.compute(mesh.positions, [0, 1, 2, 3], 0))
            .toThrow('The index array must have 3 indices per triangle.');
        expect(() => map.compute(mesh.positions, [], 0))
            .toThrow('The index array must have 3 indices per triangle.');
        expect(() => map.compute(mesh.positions, mesh.indices, -1))
            .toThrow('Invalid puncture triangle.');
        expect(() => map.compute(mesh.positions, mesh.indices,
            mesh.indices.length / 3)).toThrow('Invalid puncture triangle.');
    });

    it('rejects a mesh with a boundary', () => {
        // Remove one triangle of the octahedron, leaving three boundary
        // edges. Upstream dereferences the null second triangle pointer of a
        // boundary edge; the port reports the precondition instead.
        const mesh = makeSphereMesh(0, false);
        const indices = mesh.indices.slice(3);
        const map = new ConformalMapGenus0();
        expect(() => map.compute(mesh.positions, indices, 0))
            .toThrow('The mesh must be a closed manifold surface.');
    });
});
