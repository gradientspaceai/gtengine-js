// Exact (bigint) geometric predicates for verification tests.
//
// Every finite binary64 value is m * 2^e with an integer mantissa m, so a set
// of doubles can be scaled by a single positive power of two into exact
// integers. The orientation / in-circle / in-sphere determinants below are
// homogeneous in the coordinates, so a positive uniform scale never changes
// their sign: computing them on the scaled integers with bigint arithmetic
// gives the exact sign of the floating-point predicate.

const scratch = new ArrayBuffer(8);
const scratchF = new Float64Array(scratch);
const scratchU = new BigUint64Array(scratch);

/** Decompose a finite double as m * 2^e with integer m (m = 0 for zero). */
function decompose(x: number): { m: bigint; e: number } {
    if (x === 0) {
        return { m: 0n, e: 0 };
    }
    scratchF[0] = x;
    const bits = scratchU[0];
    const sign = (bits >> 63n) === 1n ? -1n : 1n;
    const biased = Number((bits >> 52n) & 0x7ffn);
    let mantissa = bits & 0xfffffffffffffn;
    let e: number;
    if (biased === 0) {
        e = -1074;                      // subnormal
    } else {
        mantissa |= 0x10000000000000n;  // implicit leading 1
        e = biased - 1075;
    }
    return { m: sign * mantissa, e };
}

/**
 * Scale a set of finite doubles by a common positive power of two so that all
 * of them become exact integers, returned as bigints. Pass every coordinate
 * that participates in one predicate evaluation in a single call.
 */
export function exactDyadic(values: readonly number[]): bigint[] {
    let eMin = 0;
    let any = false;
    const parts = values.map(v => decompose(v));
    for (const p of parts) {
        if (p.m !== 0n) {
            if (!any || p.e < eMin) {
                eMin = p.e;
                any = true;
            }
        }
    }
    return parts.map(p => (p.m === 0n ? 0n : p.m << BigInt(p.e - eMin)));
}

function sign(x: bigint): number {
    return x > 0n ? 1 : (x < 0n ? -1 : 0);
}

/**
 * Sign of the 2D orientation determinant of (a, b, c): +1 when the triangle
 * is counterclockwise, -1 when clockwise, 0 when the points are collinear.
 * The arguments are exact integer coordinates on a common scale.
 */
export function orient2(ax: bigint, ay: bigint, bx: bigint, by: bigint,
    cx: bigint, cy: bigint): number {
    return sign((bx - ax) * (cy - ay) - (by - ay) * (cx - ax));
}

/**
 * Sign of the in-circle determinant of (a, b, c; d). When (a, b, c) is
 * counterclockwise the result is +1 when d is strictly inside the
 * circumcircle, 0 when d is on it and -1 when d is strictly outside.
 */
export function inCircle2(ax: bigint, ay: bigint, bx: bigint, by: bigint,
    cx: bigint, cy: bigint, dx: bigint, dy: bigint): number {
    const a0 = ax - dx, a1 = ay - dy, a2 = a0 * a0 + a1 * a1;
    const b0 = bx - dx, b1 = by - dy, b2 = b0 * b0 + b1 * b1;
    const c0 = cx - dx, c1 = cy - dy, c2 = c0 * c0 + c1 * c1;
    return sign(a0 * (b1 * c2 - b2 * c1)
        - a1 * (b0 * c2 - b2 * c0)
        + a2 * (b0 * c1 - b1 * c0));
}

/**
 * Sign of det[b - a, c - a, d - a] (the rows are the difference vectors):
 * +1 when (a, b, c, d) is a positively oriented tetrahedron, 0 when the four
 * points are coplanar.
 */
export function orient3(a: readonly bigint[], b: readonly bigint[],
    c: readonly bigint[], d: readonly bigint[]): number {
    const u0 = b[0] - a[0], u1 = b[1] - a[1], u2 = b[2] - a[2];
    const v0 = c[0] - a[0], v1 = c[1] - a[1], v2 = c[2] - a[2];
    const w0 = d[0] - a[0], w1 = d[1] - a[1], w2 = d[2] - a[2];
    return sign(u0 * (v1 * w2 - v2 * w1)
        - u1 * (v0 * w2 - v2 * w0)
        + u2 * (v0 * w1 - v1 * w0));
}

/**
 * Sign of the in-sphere determinant of (a, b, c, d; e). When (a, b, c, d) is
 * a positively oriented tetrahedron (orient3 > 0) the result is -1 when e is
 * strictly inside the circumsphere, 0 when e is on it and +1 when e is
 * strictly outside. (The 4x4 expansion carries the opposite sign from the
 * 3x3 in-circle determinant of inCircle2, which is +1 for a point inside.)
 */
export function inSphere3(a: readonly bigint[], b: readonly bigint[],
    c: readonly bigint[], d: readonly bigint[], e: readonly bigint[]): number {
    const row = (p: readonly bigint[]): bigint[] => {
        const x = p[0] - e[0], y = p[1] - e[1], z = p[2] - e[2];
        return [x, y, z, x * x + y * y + z * z];
    };
    const m = [row(a), row(b), row(c), row(d)];

    // Expansion of the 4x4 determinant along the first row.
    const minor = (r0: number, r1: number, r2: number,
        c0: number, c1: number, c2: number): bigint =>
        m[r0][c0] * (m[r1][c1] * m[r2][c2] - m[r1][c2] * m[r2][c1])
        - m[r0][c1] * (m[r1][c0] * m[r2][c2] - m[r1][c2] * m[r2][c0])
        + m[r0][c2] * (m[r1][c0] * m[r2][c1] - m[r1][c1] * m[r2][c0]);

    const det = m[0][0] * minor(1, 2, 3, 1, 2, 3)
        - m[0][1] * minor(1, 2, 3, 0, 2, 3)
        + m[0][2] * minor(1, 2, 3, 0, 1, 3)
        - m[0][3] * minor(1, 2, 3, 0, 1, 2);
    return sign(det);
}

/** Exact twice-signed-area of a polygon given as exact integer coordinates. */
export function twiceSignedAreaExact(poly: readonly (readonly bigint[])[]): bigint {
    let area = 0n;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
        area += poly[j][0] * poly[i][1] - poly[i][0] * poly[j][1];
    }
    return area;
}

/**
 * The IEEE 754 remainder of x and y, computed exactly: x - n*y where n is the
 * integer nearest to the exact quotient x/y, ties to even. This is what C++
 * `std::remainder` returns (it is specified to compute the exact value), and
 * it is the reference for the argument reduction inside the range-reduced
 * trigonometric estimates. The quotient is formed on the common dyadic scale
 * with bigint arithmetic, so no rounding enters before the final result,
 * which is representable because |result| <= |y|/2.
 */
export function exactRemainder(x: number, y: number): number {
    if (!Number.isFinite(x) || Number.isNaN(y) || y === 0) {
        return Number.NaN;
    }
    if (!Number.isFinite(y)) {
        return x;
    }

    const [X, Y] = exactDyadic([x, y]);
    // exactDyadic scales by the smallest exponent among the nonzero inputs.
    const scale = x === 0 ? decompose(y).e
        : Math.min(decompose(x).e, decompose(y).e);
    const absX = X < 0n ? -X : X;
    const absY = Y < 0n ? -Y : Y;
    let q = absX / absY;
    const rest = absX - q * absY;
    if (rest * 2n > absY || (rest * 2n === absY && (q & 1n) === 1n)) {
        q += 1n;
    }
    const negative = (X < 0n) !== (Y < 0n);
    const R = X - (negative ? -q : q) * Y;
    const magnitude = Number(R < 0n ? -R : R) * Math.pow(2, scale);
    // std::remainder(-x, y) = -remainder(x, y), sign of zero included.
    const negateResult = R < 0n || (R === 0n && (x < 0 || Object.is(x, -0)));
    return negateResult ? -magnitude : magnitude;
}
