import { describe, it, expect } from 'vitest';
import { Triangle } from '../src/Triangle.js';
import { Vector, sub, dot, length } from '../src/Vector.js';

describe('Triangle construction', () => {
    it('the default constructor is the standard corner triangle', () => {
        const triangle = new Triangle(3);
        expect(triangle.dimension).toBe(3);
        expect(triangle.v.length).toBe(3);
        expect(triangle.v[0].values).toEqual([0, 0, 0]);
        expect(triangle.v[1].values).toEqual([1, 0, 0]);
        expect(triangle.v[2].values).toEqual([0, 1, 0]);
    });

    it('the default constructor works in 2D', () => {
        const triangle = new Triangle(2);
        expect(triangle.v[0].values).toEqual([0, 0]);
        expect(triangle.v[1].values).toEqual([1, 0]);
        expect(triangle.v[2].values).toEqual([0, 1]);
    });

    it('fromVertices copies the input vectors', () => {
        const v0 = Vector.fromArray([1, 1]);
        const v1 = Vector.fromArray([2, 1]);
        const v2 = Vector.fromArray([1, 3]);
        const triangle = Triangle.fromVertices(v0, v1, v2);
        v0.set(0, 99);
        v1.set(1, 99);
        v2.set(0, 99);
        expect(triangle.v[0].values).toEqual([1, 1]);
        expect(triangle.v[1].values).toEqual([2, 1]);
        expect(triangle.v[2].values).toEqual([1, 3]);
    });

    it('fromVertexArray takes the three vertices as an array', () => {
        const triangle = Triangle.fromVertexArray([
            Vector.fromArray([0, 0]),
            Vector.fromArray([1, 0]),
            Vector.fromArray([0, 1])]);
        expect(triangle.equals(new Triangle(2))).toBe(true);
        expect(() => Triangle.fromVertexArray([new Vector(2), new Vector(2)]))
            .toThrow('Triangle: invalid number of vertices.');
    });

    it('fromVertices throws on mismatched sizes', () => {
        expect(() => Triangle.fromVertices(new Vector(2), new Vector(2),
            new Vector(3))).toThrow('Triangle: mismatched sizes.');
    });

    it('clone is a deep copy', () => {
        const triangle = new Triangle(2);
        const copy = triangle.clone();
        copy.v[1].set(0, 7);
        expect(triangle.v[1].values).toEqual([1, 0]);
        expect(copy.v[1].values).toEqual([7, 0]);
    });
});

describe('Triangle geometry derived from the vertices', () => {
    it('the default 3D triangle has edge lengths 1, sqrt(2), 1', () => {
        const t = new Triangle(3);
        expect(length(sub(t.v[1], t.v[0]))).toBe(1);
        expect(length(sub(t.v[2], t.v[1]))).toBeCloseTo(Math.SQRT2, 15);
        expect(length(sub(t.v[0], t.v[2]))).toBe(1);
    });

    it('a 3-4-5 triangle has a right angle at v0', () => {
        const t = Triangle.fromVertices(Vector.fromArray([0, 0]),
            Vector.fromArray([3, 0]), Vector.fromArray([0, 4]));
        const e1 = sub(t.v[1], t.v[0]);
        const e2 = sub(t.v[2], t.v[0]);
        expect(dot(e1, e2)).toBe(0);
        expect(length(e1)).toBe(3);
        expect(length(e2)).toBe(4);
        expect(length(sub(t.v[2], t.v[1]))).toBe(5);
    });
});

describe('Triangle comparisons', () => {
    const a = Triangle.fromVertices(Vector.fromArray([0, 0]),
        Vector.fromArray([1, 0]), Vector.fromArray([0, 1]));
    const sameAsA = Triangle.fromVertices(Vector.fromArray([0, 0]),
        Vector.fromArray([1, 0]), Vector.fromArray([0, 1]));
    // Larger first vertex.
    const b = Triangle.fromVertices(Vector.fromArray([0, 1]),
        Vector.fromArray([1, 0]), Vector.fromArray([0, 1]));
    // Same first two vertices, larger third vertex.
    const c = Triangle.fromVertices(Vector.fromArray([0, 0]),
        Vector.fromArray([1, 0]), Vector.fromArray([0, 2]));

    it('equals/notEquals compare the vertex arrays', () => {
        expect(a.equals(sameAsA)).toBe(true);
        expect(a.notEquals(sameAsA)).toBe(false);
        expect(a.equals(c)).toBe(false);
        expect(a.notEquals(b)).toBe(true);
    });

    it('orders lexicographically by the vertex array', () => {
        expect(a.lessThan(b)).toBe(true);
        expect(b.lessThan(a)).toBe(false);
        expect(a.lessThan(c)).toBe(true);
        expect(c.lessThan(a)).toBe(false);
        expect(a.lessThanOrEqual(sameAsA)).toBe(true);
        expect(a.greaterThanOrEqual(sameAsA)).toBe(true);
        expect(a.greaterThan(sameAsA)).toBe(false);
        expect(c.greaterThan(a)).toBe(true);
        expect(c.greaterThanOrEqual(a)).toBe(true);
    });
});
