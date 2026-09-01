// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ArbitraryPrecision.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Upstream ArbitraryPrecision.h is an umbrella header with no code of its own.
// It includes, in this order:
//
//   UIntegerALU32.h, UIntegerAP32.h, UIntegerFP32.h,
//   BSNumber.h, BSRational.h, BSPrecision.h
//
// The three UInteger*32 headers are the C++ 32-bit-word storage backends for
// the arbitrary-precision significands. They are recorded as `omitted` in the
// port manifest: the port's BSNumber/BSRational are built on the JavaScript
// `bigint` type (precedent set by B34), so there is no UInteger template
// parameter and no consumer for those backends. Nothing they expose has an
// analogue here.
//
// This module therefore re-exports the arbitrary-precision surface that the
// umbrella header makes available in C++, so `import { ... } from
// './ArbitraryPrecision'` pulls in the same working set as
// `#include <Mathematics/ArbitraryPrecision.h>`.
//
// Note: APConversion.h and APInterval.h are separate upstream headers that
// ArbitraryPrecision.h does not include; they are not re-exported here.

export {
  BSNumber,
  BSNumberRoundingMode,
  convertBSNumber,
} from './BSNumber';

export {
  BSRational,
  convertBSRational,
  convertBSRationalToBSNumber,
  convertBSRationalToNumber,
  convertBSRationalToFloat32,
} from './BSRational';

export {
  BSPrecision,
  BSPrecisionParameters,
  BSPrecisionType,
} from './BSPrecision';
