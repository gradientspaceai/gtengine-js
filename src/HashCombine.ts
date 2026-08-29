// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) HashCombine.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Support for creating hash values for a list of values.
//
// The code here comes from the book
//   Nicolai M. Josuttis, "The C++ Standard Library: A Tutorial
//   and Reference, 2nd edition", Addison-Wesley Professional,
//   March 2012, Section 7.9.2, pp. 364-365.
// Credit for the hash_combine concept is from
//   https://www.boost.org/doc/libs/1_35_0/doc/html/hash/combine.html
// The magic number and shifts are based on the paper
//   Timothy C. Hoad and Justin Zobel, "Methods for Identifying Versioned
//   and Plagiarised Documents", Journal of the American Society for
//   Information Science and Technology, vol. 54, no. 3, February 2003.
//   https://dl.acm.org/doi/10.1002/asi.10170
//
// Port notes: upstream is generic over any T with a std::hash<T>
// specialization and mutates a size_t seed by reference. The port hashes
// 'number' values (the numeric type of the library), returns the new seed
// instead of mutating an out-parameter, and uses unsigned 32-bit seeds
// (C++ size_t hash values are implementation-defined, so no C++
// implementation's exact hash values are reproduced; the combine formula is
// preserved in 32-bit arithmetic).

const scratchBuffer = new ArrayBuffer(8);
const scratchF64 = new Float64Array(scratchBuffer);
const scratchU32 = new Uint32Array(scratchBuffer);

// Hash a number by folding the 64 bits of its IEEE binary64 representation
// into 32 bits. This is the port's stand-in for std::hash<double>. All zero
// values hash alike (+0 and -0 are equal numbers upstream as well as here).
function hashNumber(value: number): number {
    scratchF64[0] = value === 0 ? 0 : value; // normalize -0 to +0
    return (scratchU32[0] ^ scratchU32[1]) >>> 0;
}

// The port of HashCombine(seed, value): returns the updated seed.
export function hashCombine(seed: number, value: number): number {
    seed = seed >>> 0;
    return (seed ^ (hashNumber(value) + 0x9e3779b9 + ((seed << 6) >>> 0) + (seed >>> 2))) >>> 0;
}

// The port of the seeded HashValue overloads: combines each value in order
// into the seed and returns the result.
export function hashValueWithSeed(seed: number, ...values: number[]): number {
    for (const value of values) {
        seed = hashCombine(seed, value);
    }
    return seed >>> 0;
}

// The port of HashValue(arguments...): creates a hash value from a list of
// values using seed 0.
export function hashValue(...values: number[]): number {
    return hashValueWithSeed(0, ...values);
}
