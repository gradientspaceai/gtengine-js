// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Constants.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The literals below are copied verbatim from upstream Constants.h. Each
// parses to exactly the same IEEE binary64 value as the C++ double constant.
// Note that some differ in the last bit from the naively computed value
// (e.g. GTE_C_INV_SQRT_2 !== Math.SQRT1_2); the upstream literals win.

// Constants involving pi.
export const GTE_C_PI = 3.1415926535897931;
export const GTE_C_HALF_PI = 1.5707963267948966;
export const GTE_C_QUARTER_PI = 0.7853981633974483;
export const GTE_C_TWO_PI = 6.2831853071795862;
export const GTE_C_INV_PI = 0.3183098861837907;
export const GTE_C_INV_TWO_PI = 0.15915494309189535;
export const GTE_C_INV_HALF_PI = 0.63661977236758138;

// Conversions between degrees and radians.
export const GTE_C_DEG_TO_RAD = 0.017453292519943295;
export const GTE_C_RAD_TO_DEG = 57.295779513082321;

// Common constants.
export const GTE_C_SQRT_2 = 1.4142135623730951;
export const GTE_C_INV_SQRT_2 = 0.7071067811865475;
export const GTE_C_LN_2 = 0.6931471805599453;
export const GTE_C_INV_LN_2 = 1.4426950408889634;
export const GTE_C_LN_10 = 2.3025850929940459;
export const GTE_C_INV_LN_10 = 0.43429448190325176;
