// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DCPQuery.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Distance and closest-point queries.
//
// Upstream declares 'template <typename Real, typename Type0, typename Type1>
// class DCPQuery {}' and each Dist* header defines a specialization with a
// nested Result struct and an operator() member. TypeScript has no template
// specialization, so the port uses an interface: each Dist* file exports a
// concrete class named after its file (e.g. DistPoint3Triangle3) that
// implements DCPQuery, plus an exported Result type named after the file
// (e.g. DistPoint3Triangle3Result). An interface (rather than a base class)
// is used because the upstream primary template is empty -- there is no
// shared state or behavior to inherit -- and an interface leaves the
// implementing classes free to define their own construction. The Real
// template parameter is dropped (always 'number' in the port) and the query
// result type is added as a type parameter instead.
//
// Per PORTING.md, 'DCPQuery::operator()' is ported as 'compute(...)'.
export interface DCPQuery<Type0, Type1, Result> {
    compute(primitive0: Type0, primitive1: Type1): Result;
}
