// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Logger.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Upstream Logger.h provides the macros GTE_ASSERT, GTE_ERROR, LogAssert and
// LogError, which throw std::runtime_error with a message prefixed by the
// source file, function and line of the expansion site. This port provides
// them as functions that throw Error with the caller's message text; the
// JavaScript Error stack trace supplies the source-location information that
// the C++ macros compose manually, so no file/function/line prefix is added.
//
// All ported files use these helpers wherever upstream uses the macros:
//   LogAssert(condition, message)  ->  logAssert(condition, message)
//   LogError(message)              ->  logError(message)

// Throws Error(message) when the condition is false. The 'asserts' return
// type lets TypeScript narrow types after the call, mirroring how code after
// a C++ assertion may assume the condition holds.
export function logAssert(condition: boolean, message: string): asserts condition {
    if (!condition) {
        throw new Error(message);
    }
}

// Unconditionally throws Error(message). Used for unexpected conditions.
export function logError(message: string): never {
    throw new Error(message);
}
