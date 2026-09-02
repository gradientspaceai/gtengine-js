import { describe, it, expect } from 'vitest';
import { IntrLine2SegmentMesh2FI } from '../src/IntrLine2SegmentMesh2';
import { Line } from '../src/Line';
import { SegmentMesh } from '../src/SegmentMesh';
import { Vector, add, mul, normalize, sub } from '../src/Vector';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function line(origin: number[], direction: number[]): Line {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Line.fromOriginDirection(Vector.fromArray(origin), d);
}

// The unit square as a closed polyline.
function unitSquare(): SegmentMesh {
    return SegmentMesh.fromContiguous(
        [vec(0, 0), vec(1, 0), vec(1, 1), vec(0, 1)], false);
}

describe('IntrLine2SegmentMesh2FI', () => {
    const query = new IntrLine2SegmentMesh2FI();

    it('finds the two crossings of a horizontal line through a square', () => {
        const mesh = unitSquare();
        const result = query.find(line([-5, 0.5], [1, 0]), mesh);
        expect(result.intersections.length).toBe(2);

        // Sorted by line parameter, so the left crossing comes first.
        expect(result.intersections[0].point.values[0]).toBeCloseTo(0, 12);
        expect(result.intersections[0].point.values[1]).toBeCloseTo(0.5, 12);
        expect(result.intersections[1].point.values[0]).toBeCloseTo(1, 12);
        expect(result.intersections[1].point.values[1]).toBeCloseTo(0.5, 12);
        expect(result.intersections[0].lineParameter)
            .toBeLessThan(result.intersections[1].lineParameter);

        // The left crossing is on the segment <V[3],V[0]>, the right one is
        // on <V[1],V[2]>.
        expect(result.intersections[0].indexPair).toEqual([3, 0]);
        expect(result.intersections[1].indexPair).toEqual([1, 2]);
    });

    it('reports the mesh-segment parameter of each crossing', () => {
        const mesh = SegmentMesh.fromDisjoint([vec(0, -1), vec(0, 3)]);
        const result = query.find(line([-1, 0], [1, 0]), mesh);
        expect(result.intersections.length).toBe(1);
        // The crossing is at (0,0), which is 1/4 of the way along the
        // segment from (0,-1) to (0,3).
        expect(result.intersections[0].meshSegmentParameter)
            .toBeCloseTo(0.25, 12);
        expect(result.intersections[0].lineParameter).toBeCloseTo(1, 12);
    });

    it('reports both endpoints when the line and a mesh segment are coincident', () => {
        const mesh = SegmentMesh.fromDisjoint([vec(2, 0), vec(5, 0)]);
        const result = query.find(line([0, 0], [1, 0]), mesh);
        expect(result.intersections.length).toBe(2);
        expect(result.intersections[0].point.values[0]).toBeCloseTo(2, 12);
        expect(result.intersections[1].point.values[0]).toBeCloseTo(5, 12);
        expect(result.intersections[0].lineParameter)
            .toBeLessThan(result.intersections[1].lineParameter);
    });

    it('returns an empty list when the line misses every segment', () => {
        const mesh = unitSquare();
        const result = query.find(line([-5, 5], [1, 0]), mesh);
        expect(result.intersections.length).toBe(0);
    });

    it('sorts every intersection by line parameter', () => {
        // A zig-zag polyline crossed by a diagonal line.
        const mesh = SegmentMesh.fromContiguous([
            vec(-3, -1), vec(-2, 1), vec(-1, -1), vec(0, 1), vec(1, -1),
            vec(2, 1), vec(3, -1)
        ], true);
        const result = query.find(line([-4, 0], [1, 0]), mesh);
        expect(result.intersections.length).toBe(6);
        for (let i = 1; i < result.intersections.length; ++i) {
            expect(result.intersections[i - 1].lineParameter)
                .toBeLessThanOrEqual(result.intersections[i].lineParameter);
        }
        for (const record of result.intersections) {
            expect(record.point.values[1]).toBeCloseTo(0, 12);
        }
    });

    it('reports points consistent with both parameterizations on random lines', () => {
        let seed = 13579;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const vertices: Vector[] = [];
        for (let i = 0; i < 12; ++i) {
            const theta = (2 * Math.PI * i) / 12;
            vertices.push(vec(2 * Math.cos(theta), 1.3 * Math.sin(theta)));
        }
        const mesh = SegmentMesh.fromContiguous(vertices, false);
        const meshVertices = mesh.getVertices();

        for (let trial = 0; trial < 150; ++trial) {
            const l = line([rand() * 6 - 3, rand() * 6 - 3],
                [rand() * 2 - 1, rand() * 2 - 1]);
            const result = query.find(l, mesh);
            for (const record of result.intersections) {
                const onLine = add(l.origin,
                    mul(record.lineParameter, l.direction));
                expect(sub(onLine, record.point).values[0]).toBeCloseTo(0, 8);
                expect(sub(onLine, record.point).values[1]).toBeCloseTo(0, 8);

                const p0 = meshVertices[record.indexPair[0]];
                const p1 = meshVertices[record.indexPair[1]];
                const onSegment = add(p0,
                    mul(record.meshSegmentParameter, sub(p1, p0)));
                expect(sub(onSegment, record.point).values[0])
                    .toBeCloseTo(0, 8);
                expect(sub(onSegment, record.point).values[1])
                    .toBeCloseTo(0, 8);
                expect(record.meshSegmentParameter).toBeGreaterThanOrEqual(0);
                expect(record.meshSegmentParameter).toBeLessThanOrEqual(1);
            }
            for (let i = 1; i < result.intersections.length; ++i) {
                expect(result.intersections[i - 1].lineParameter)
                    .toBeLessThanOrEqual(
                        result.intersections[i].lineParameter);
            }
        }
    });
});
