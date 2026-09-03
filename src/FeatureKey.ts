// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) FeatureKey.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// FeatureKey is the base class for the keys that identify the features
// (edges, triangles, tetrahedra, ...) of a mesh. The key stores the N vertex
// indices of the feature in the V array.
//
// An ordered feature key has V[0] = min(V[]) with (V[0],V[1],...,V[N-1]) a
// permutation of N inputs with an even number of transpositions.
//
// An unordered feature key has V[0] < V[1] < ... < V[N-1].
//
// Note that the word 'order' is about the geometry of the feature, not the
// comparison order for any sorting.
//
// Port notes: upstream is the class template FeatureKey<N, Ordered> whose
// derived classes (EdgeKey, TriangleKey, TetrahedronKey) fill in V[]. The
// port is a single concrete class: the compile-time N becomes the runtime
// constructor argument 'n' and the compile-time Ordered becomes the readonly
// 'ordered' property, which derived classes pass to super(). Upstream leaves
// V[] uninitialized in the abstract base; per the port conventions the
// constructor zero-fills it. The C++ comparison operators become methods
// (equals, notEqual, lessThan, lessThanOrEqual, greaterThan,
// greaterThanOrEqual). The two upstream operator() overloads form the
// std::unordered_* support: the hash function becomes hashValue() (and the
// static FeatureKey.hashValue) and the bucket equality comparison becomes the
// static FeatureKey.equal. Because JavaScript Map and Set compare keys by
// identity rather than by a hash/equality pair, mapKey() is provided for use
// as a primitive Map/Set key; it is a port addition with no upstream
// counterpart.

import { logAssert } from './Logger.js';
import { hashCombine } from './HashCombine.js';

export class FeatureKey {
    // The vertex indices of the feature. The length is the N of the upstream
    // FeatureKey<N, Ordered>.
    V: number[];

    // The port of the upstream 'Ordered' template argument. It is not used by
    // any comparison; the derived classes use it to decide how to arrange the
    // vertex indices in V.
    readonly ordered: boolean;

    // The port of the protected upstream default constructor. It is public
    // here so that FeatureKey may be used directly for feature sizes that
    // have no dedicated derived class.
    constructor(n: number, ordered: boolean) {
        logAssert(n >= 1, 'The number of vertex indices must be positive.');
        this.V = new Array<number>(n).fill(0);
        this.ordered = ordered;
    }

    // The number N of vertex indices in the key.
    get n(): number {
        return this.V.length;
    }

    // The upstream comparisons loop over the N indices rather than using the
    // std::array comparisons, which was done to improve the speed of the
    // comparisons for the C++ Standard Library that ships with Microsoft
    // Visual Studio 2019 16.7.*. The loops are preserved here.
    equals(key: FeatureKey): boolean {
        for (let i = 0; i < this.V.length; ++i) {
            if (this.V[i] !== key.V[i]) {
                return false;
            }
        }
        return true;
    }

    notEqual(key: FeatureKey): boolean {
        return !this.equals(key);
    }

    // Lexicographical ordering of the vertex indices.
    lessThan(key: FeatureKey): boolean {
        for (let i = 0; i < this.V.length; ++i) {
            if (this.V[i] < key.V[i]) {
                return true;
            }
            if (this.V[i] > key.V[i]) {
                return false;
            }
        }
        return false;
    }

    lessThanOrEqual(key: FeatureKey): boolean {
        return !key.lessThan(this);
    }

    greaterThan(key: FeatureKey): boolean {
        return key.lessThan(this);
    }

    greaterThanOrEqual(key: FeatureKey): boolean {
        return !this.lessThan(key);
    }

    // The port of the upstream hash function operator(). The hash is an
    // unsigned 32-bit value; see the port notes of HashCombine for why the
    // values are not those of any particular C++ implementation.
    hashValue(): number {
        let seed = 0;
        for (const value of this.V) {
            seed = hashCombine(seed, value);
        }
        return seed;
    }

    // The functor form of the hash function operator().
    static hashValue(key: FeatureKey): number {
        return key.hashValue();
    }

    // The port of the upstream equality comparison operator(), used by
    // std::unordered_* for elements in the same bucket.
    static equal(key0: FeatureKey, key1: FeatureKey): boolean {
        return key0.equals(key1);
    }

    // Port addition: a primitive key suitable for JavaScript Map and Set,
    // which compare keys by identity. Two FeatureKey objects have the same
    // mapKey() if and only if equals() returns true for them (keys of
    // different lengths never collide because the length is encoded by the
    // number of separators).
    mapKey(): string {
        return this.V.join(',');
    }

    // Comparison function for sorting arrays of keys in the order that
    // std::set<FeatureKey> and std::map<FeatureKey,...> iterate.
    static compare(key0: FeatureKey, key1: FeatureKey): number {
        if (key0.lessThan(key1)) {
            return -1;
        }
        if (key1.lessThan(key0)) {
            return 1;
        }
        return 0;
    }
}
