// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) TIQuery.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Test-intersection queries.
//
// Upstream declares 'template <typename T, typename Type0, typename Type1>
// class TIQuery {}' and each Intr* header defines a specialization with a
// nested Result struct and an operator() member. The port uses an interface
// (see DCPQuery.ts for the rationale): each Intr* file exports a concrete
// class named after its file with a 'TI' suffix (e.g. IntrRay3Sphere3TI)
// implementing TIQuery, plus an exported result type (e.g.
// IntrRay3Sphere3TIResult). The T/Real template parameter is dropped (always
// 'number' in the port) and the query result type is a type parameter.
//
// Per PORTING.md, 'TIQuery::operator()' is ported as 'test(...)'.
export interface TIQuery<Type0, Type1, Result> {
    test(primitive0: Type0, primitive1: Type1): Result;
}
