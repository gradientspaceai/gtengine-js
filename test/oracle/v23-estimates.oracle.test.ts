// Replays oracle/cpp/cases/v23-estimates.cpp. Keep the two files in the same
// order.
//
// Every scalar estimate and every ChebyshevRatioEstimate/SlerpEstimate case
// is { exact: true }: the polynomial evaluations are + - * and sqrt, and the
// range reductions use floor, frexp, ldexp, remainder and fabs, which IEEE
// 754 specifies exactly. ChebyshevRatio.h and Slerp.h call std::sin and
// std::acos and carry the default 1e-12 tolerance, except for their
// exactly evaluated outputs (the angle-0 branch and the +-0 components),
// which are emitted with outRealExact.
//
// Independent reference checks (ORACLE.md, "Agreement is not correctness")
// run on every replayed record: each estimate inside its documented domain
// must lie within its published GetMaxError bound of the true function
// (Math.* in double precision, accurate to about 1 ulp, far below the 1e-2
// to 3e-14 bounds); Slerp results must be unit vectors on the great arc at
// the right angle fraction; ChebyshevRatio must satisfy
// f(t, A) sin(A) = sin(t A). Branch histograms of the range reductions are
// collected too (the RR branches are recomputed with an exact BigInt
// remainder) and the last tests require every branch to be reached.
import { describe, expect, it } from 'vitest';
import { acosEstimate, getACosEstimateMaxError } from '../../src/ACosEstimate.js';
import { asinEstimate, getASinEstimateMaxError } from '../../src/ASinEstimate.js';
import { atanEstimate, atanEstimateRR, getATanEstimateMaxError } from '../../src/ATanEstimate.js';
import {
    chebyshevRatio, chebyshevRatios, chebyshevRatiosUsingCosAngle, chebyshevRatioUsingCosAngle
} from '../../src/ChebyshevRatio.js';
import {
    chebyshevRatioEstimate, chebyshevRatioEstimateR, getChebyshevRatioEstimateMaxError,
    getChebyshevRatioEstimateRMaxError
} from '../../src/ChebyshevRatioEstimate.js';
import { GTE_C_HALF_PI, GTE_C_INV_SQRT_2, GTE_C_PI, GTE_C_QUARTER_PI, GTE_C_TWO_PI } from '../../src/Constants.js';
import { cosEstimate, cosEstimateRR, getCosEstimateMaxError } from '../../src/CosEstimate.js';
import { exp2Estimate, exp2EstimateRR, getExp2EstimateMaxError } from '../../src/Exp2Estimate.js';
import { expEstimate, expEstimateRR, getExpEstimateMaxError } from '../../src/ExpEstimate.js';
import { getInvSqrtEstimateMaxError, invSqrtEstimate, invSqrtEstimateRR } from '../../src/InvSqrtEstimate.js';
import { getLog2EstimateMaxError, log2Estimate, log2EstimateRR } from '../../src/Log2Estimate.js';
import { getLogEstimateMaxError, logEstimate, logEstimateRR } from '../../src/LogEstimate.js';
import { getSinEstimateMaxError, sinEstimate, sinEstimateRR } from '../../src/SinEstimate.js';
import { slerp, slerpUsingCosAngle, slerpUsingMidpoint } from '../../src/Slerp.js';
import {
    slerpEstimate, slerpEstimateUsingCosAngle, slerpEstimateUsingMidpoint
} from '../../src/SlerpEstimate.js';
import { getSqrtEstimateMaxError, sqrtEstimate, sqrtEstimateRR } from '../../src/SqrtEstimate.js';
import { getTanEstimateMaxError, tanEstimate, tanEstimateRR } from '../../src/TanEstimate.js';
import { OracleFamily, type OracleIO } from './harness.js';

const NUM_MODES = 8;
const range = (lo: number, hi: number, step = 1): number[] => {
    const a: number[] = [];
    for (let d = lo; d <= hi; d += step) { a.push(d); }
    return a;
};
const DEG_1_7 = range(1, 7);
const DEG_1_8 = range(1, 8);
const DEG_SIN = range(3, 11, 2);
const DEG_COS = range(2, 10, 2);
const DEG_TAN = range(3, 13, 2);
const DEG_1_12 = range(1, 12);
const DEG_1_16 = range(1, 16);

// ---- reference checks and branch histograms ------------------------------

interface RefStats { checked: number; violations: string[]; worst: number }
const refStats = new Map<string, RefStats>();
const histograms = new Map<string, Map<string, number>>();

// Records |err| / bound (the worst ratio is reported) and a violation when
// the error exceeds bound + slack.
function refCheck(key: string, err: number, bound: number, slack: number, what: () => string): void {
    let s = refStats.get(key);
    if (s === undefined) { s = { checked: 0, violations: [], worst: 0 }; refStats.set(key, s); }
    ++s.checked;
    const ratio = Math.abs(err) / bound;
    if (!(ratio <= s.worst)) { s.worst = Number.isNaN(ratio) ? Infinity : ratio; }
    if (!(Math.abs(err) <= bound + slack) && s.violations.length < 5) {
        s.violations.push(`${what()}: error ${err} bound ${bound}`);
    }
}

function count(key: string, branch: string): void {
    let h = histograms.get(key);
    if (h === undefined) { h = new Map(); histograms.set(key, h); }
    h.set(branch, (h.get(branch) ?? 0) + 1);
}

// ---- exact IEEE remainder over BigInt (independent of the port) ----------

function decompose(x: number): { m: bigint; e: number } {
    const dv = new DataView(new ArrayBuffer(8));
    dv.setFloat64(0, x);
    const bits = dv.getBigUint64(0);
    const sign = (bits >> 63n) !== 0n ? -1n : 1n;
    const biased = Number((bits >> 52n) & 0x7ffn);
    const frac = bits & 0xfffffffffffffn;
    if (biased === 0) { return { m: sign * frac, e: -1074 }; }
    return { m: sign * (frac | 0x10000000000000n), e: biased - 1075 };
}

function scale2(v: number, e: number): number {
    // v * 2^e for integer-valued v whose result is representable; the steps
    // are exact because the final value is.
    while (e > 600) { v *= 2 ** 600; e -= 600; }
    while (e < -600) { v *= 2 ** -600; e += 600; }
    return v * 2 ** e;
}

// std::remainder(x, y) for finite x and finite nonzero y, plus whether the
// quotient x/y is exactly halfway between two integers.
function exactRemainder(x: number, y: number): { r: number; tie: boolean } {
    const a = decompose(x);
    const b = decompose(y);
    const e = Math.min(a.e, b.e);
    const X = a.m << BigInt(a.e - e);
    const Y = b.m < 0n ? -(b.m << BigInt(b.e - e)) : b.m << BigInt(b.e - e);
    let q = X / Y;
    let r = X - q * Y;
    const twice = 2n * (r < 0n ? -r : r);
    const tie = twice === Y;
    if (twice > Y || (tie && (q & 1n) !== 0n)) {
        if (r > 0n) { r -= Y; ++q; } else { r += Y; --q; }
    }
    const value = scale2(Number(r), e);
    // A zero remainder takes the sign of x.
    return { r: value === 0 ? (Object.is(x, -0) || x < 0 ? -0 : 0) : value, tie };
}

// ---- scalar estimate replays ---------------------------------------------

type Estimate = (x: number, degree: number) => number;
type Check = (x: number, degree: number, y: number) => void;

// A domain check: inside [lo, hi] the estimate is within its published bound
// of truth(x) (absolute error), up to an evaluation-rounding slack.
function within(key: string, lo: number, hi: number, truth: (x: number) => number,
    maxError: (d: number) => number, slack = 4e-16): Check {
    return (x, d, y) => {
        if (!(x >= lo && x <= hi)) { return; }
        refCheck(`${key} degree ${d}`, y - truth(x), maxError(d), slack, () => `${key}(${x})`);
    };
}

// A relative check for the RR variants that scale by 2^p: the relative error
// is at most the reduced-range absolute bound divided by the smallest value
// the reduced function takes (factor = 1/min). Only finite normal results.
function withinRelative(key: string, admit: (x: number) => boolean,
    truth: (x: number) => number, maxError: (d: number) => number, factor: number,
    slack: (x: number) => number = () => 1e-15): Check {
    return (x, d, y) => {
        const z = truth(x);
        if (!admit(x) || !Number.isFinite(z) || Math.abs(z) < 2.2250738585072014e-308) { return; }
        refCheck(`${key} degree ${d}`, (y - z) / z, maxError(d) * factor, slack(x),
            () => `${key}(${x})`);
    };
}

function estimateCase(family: OracleFamily, name: string, degrees: readonly number[],
    estimate: Estimate, checks: Check[] = [], branch?: (x: number) => string): void {
    family.case(name, (io) => {
        for (let m = 0; m < NUM_MODES; ++m) {
            const x = io.real();
            if (branch !== undefined) { count(name, branch(x)); }
            for (const d of degrees) {
                const y = estimate(x, d);
                io.outReal(y);
                for (const c of checks) { c(x, d, y); }
            }
        }
    }, { exact: true });
}

function maxErrorCase(family: OracleFamily, name: string, degrees: readonly number[],
    maxError: (d: number) => number): void {
    family.case(name, (io) => {
        for (const d of degrees) { io.outReal(maxError(d)); }
    }, { exact: true });
}

// Branch labels of the range reductions, recomputed independently.
function periodicBranch(x: number, period: number, threshold: number): string {
    if (!Number.isFinite(x)) { return 'non-finite'; }
    const { r, tie } = exactRemainder(x, period);
    const side = r > threshold ? 'r>t' : (r < -threshold ? 'r<-t' : '|r|<=t');
    return tie ? `${side} tie` : side;
}

function tanBranch(x: number): string {
    if (!Number.isFinite(x)) { return 'non-finite'; }
    const { r, tie } = exactRemainder(x, GTE_C_PI);
    // r > pi/2 and r < -pi/2 are unreachable (remainder(x, pi) is in
    // [-pi/2, pi/2]; UPSTREAM-FINDINGS #57), so y = r.
    const reach = r > GTE_C_HALF_PI || r < -GTE_C_HALF_PI ? ' UNREACHABLE' : '';
    const side = Math.abs(r) <= GTE_C_QUARTER_PI ? '|y|<=pi/4'
        : (r > GTE_C_QUARTER_PI ? 'y>pi/4' : 'y<-pi/4');
    return `${side}${tie ? ' tie' : ''}${reach}`;
}

function frexpBranch(x: number): string {
    if (!Number.isFinite(x)) { return 'non-finite'; }
    if (x === 0) { return Object.is(x, -0) ? '-0' : '+0'; }
    const e = decompose(x);
    const bits = e.m < 0n ? -e.m : e.m;
    const p = e.e + bits.toString(2).length - 1;  // x = f*2^p, f in [1,2)
    const sub = Math.abs(x) < 2.2250738585072014e-308 ? 'subnormal ' : '';
    const neg = x < 0 ? 'negative ' : '';
    return `${neg}${sub}p ${(p & 1) === 0 ? 'even' : 'odd'}`;
}

function exp2Branch(x: number): string {
    if (!Number.isFinite(x)) { return 'non-finite'; }
    const p = Math.floor(x);
    const y = x - p;
    const kind = p >= 1024 ? 'overflow' : (p < -1022 ? (p < -1075 ? 'underflow to 0' : 'subnormal') : 'normal');
    return `${kind}${y === 0 ? ' y=0' : ''}`;
}

function atanBranch(x: number): string {
    if (Math.abs(x) <= 1) { return '|x|<=1'; }
    if (x > 1) { return 'x>1'; }
    return Number.isNaN(x) ? 'NaN (else arm)' : 'x<-1';
}

const sinTruth = (x: number): number => Math.sin(x);
const sinRRCheck = (key: string, truth: (x: number) => number, err: (d: number) => number): Check =>
    // Upstream reduces with the double 2*pi, whose error from the true period
    // grows with the number of periods removed; beyond |x| = 1000 the
    // reduced argument is no longer the true one to the bound's precision.
    (x, d, y) => {
        if (!(Math.abs(x) <= 1000)) { return; }
        refCheck(`${key} degree ${d}`, y - truth(x), err(d), 1e-15 * Math.abs(x) + 4e-16,
            () => `${key}(${x})`);
    };

describe('oracle: v23-estimates', () => {
    const family = new OracleFamily('v23-estimates');
    const LN_2 = Math.LN2;
    const normal = (x: number): boolean => Number.isFinite(x) && x > 0;

    // Exp2Estimate.h
    estimateCase(family, 'Exp2Estimate.estimate', DEG_1_7, exp2Estimate,
        [within('exp2', 0, 1, (x) => 2 ** x, getExp2EstimateMaxError)]);
    estimateCase(family, 'Exp2Estimate.estimateRR', DEG_1_7, exp2EstimateRR,
        [withinRelative('exp2RR', Number.isFinite, (x) => 2 ** x, getExp2EstimateMaxError, 1)],
        exp2Branch);
    maxErrorCase(family, 'Exp2Estimate.getMaxError', DEG_1_7, getExp2EstimateMaxError);

    // ExpEstimate.h. x*(1/ln 2) rounds, so the relative slack grows with |x|.
    estimateCase(family, 'ExpEstimate.estimate', DEG_1_7, expEstimate,
        [within('exp', 0, LN_2, Math.exp, getExpEstimateMaxError)]);
    estimateCase(family, 'ExpEstimate.estimateRR', DEG_1_7, expEstimateRR,
        [withinRelative('expRR', Number.isFinite, Math.exp, getExpEstimateMaxError, 1,
            (x) => 1e-15 + 4e-16 * Math.abs(x))],
        (x) => exp2Branch(x * 1.4426950408889634));
    maxErrorCase(family, 'ExpEstimate.getMaxError', DEG_1_7, getExpEstimateMaxError);

    // Log2Estimate.h. poly + p rounds at ulp(p) <= 2^-42.
    estimateCase(family, 'Log2Estimate.estimate', DEG_1_8, log2Estimate,
        [within('log2', 1, 2, Math.log2, getLog2EstimateMaxError)]);
    estimateCase(family, 'Log2Estimate.estimateRR', DEG_1_8, log2EstimateRR,
        [within('log2RR', 5e-324, Number.MAX_VALUE, Math.log2, getLog2EstimateMaxError, 2.3e-13)],
        frexpBranch);
    maxErrorCase(family, 'Log2Estimate.getMaxError', DEG_1_8, getLog2EstimateMaxError);

    // LogEstimate.h. GetLogEstimateMaxError is the log2 bound (#57); the
    // extra 'logTight' check measures the error against ln 2 times it.
    estimateCase(family, 'LogEstimate.estimate', DEG_1_8, logEstimate,
        [within('log', 1, 2, Math.log, getLogEstimateMaxError),
            within('logTight', 1, 2, Math.log, (d) => LN_2 * getLog2EstimateMaxError(d))]);
    estimateCase(family, 'LogEstimate.estimateRR', DEG_1_8, logEstimateRR,
        [within('logRR', 5e-324, Number.MAX_VALUE, Math.log, getLogEstimateMaxError, 2.3e-13)],
        frexpBranch);
    maxErrorCase(family, 'LogEstimate.getMaxError', DEG_1_8, getLogEstimateMaxError);

    // SinEstimate.h
    estimateCase(family, 'SinEstimate.estimate', DEG_SIN, sinEstimate,
        [within('sin', -GTE_C_HALF_PI, GTE_C_HALF_PI, sinTruth, getSinEstimateMaxError)]);
    estimateCase(family, 'SinEstimate.estimateRR', DEG_SIN, sinEstimateRR,
        [sinRRCheck('sinRR', sinTruth, getSinEstimateMaxError)],
        (x) => periodicBranch(x, GTE_C_TWO_PI, GTE_C_HALF_PI));
    maxErrorCase(family, 'SinEstimate.getMaxError', DEG_SIN, getSinEstimateMaxError);

    // CosEstimate.h
    estimateCase(family, 'CosEstimate.estimate', DEG_COS, cosEstimate,
        [within('cos', -GTE_C_HALF_PI, GTE_C_HALF_PI, Math.cos, getCosEstimateMaxError)]);
    estimateCase(family, 'CosEstimate.estimateRR', DEG_COS, cosEstimateRR,
        [sinRRCheck('cosRR', Math.cos, getCosEstimateMaxError)],
        (x) => periodicBranch(x, GTE_C_TWO_PI, GTE_C_HALF_PI));
    maxErrorCase(family, 'CosEstimate.getMaxError', DEG_COS, getCosEstimateMaxError);

    // TanEstimate.h. The RR bound is only checked where the polynomial is
    // used directly (|y| <= pi/4); the (1 + p)/(1 - p) arms amplify the
    // polynomial's error without bound as y approaches pi/2.
    estimateCase(family, 'TanEstimate.estimate', DEG_TAN, tanEstimate,
        [within('tan', -GTE_C_QUARTER_PI, GTE_C_QUARTER_PI, Math.tan, getTanEstimateMaxError)]);
    estimateCase(family, 'TanEstimate.estimateRR', DEG_TAN, tanEstimateRR,
        [within('tanRR', -GTE_C_QUARTER_PI, GTE_C_QUARTER_PI, Math.tan, getTanEstimateMaxError)],
        tanBranch);
    maxErrorCase(family, 'TanEstimate.getMaxError', DEG_TAN, getTanEstimateMaxError);

    // ATanEstimate.h. The RR arms pi/2 - p(1/x) keep the bound (1/x rounds
    // by half an ulp, and atan' <= 1).
    estimateCase(family, 'ATanEstimate.estimate', DEG_TAN, atanEstimate,
        [within('atan', -1, 1, Math.atan, getATanEstimateMaxError)]);
    estimateCase(family, 'ATanEstimate.estimateRR', DEG_TAN, atanEstimateRR,
        [within('atanRR', -Infinity, Infinity, Math.atan, getATanEstimateMaxError, 1e-15)],
        atanBranch);
    maxErrorCase(family, 'ATanEstimate.getMaxError', DEG_TAN, getATanEstimateMaxError);

    // ACosEstimate.h and ASinEstimate.h
    estimateCase(family, 'ACosEstimate.estimate', DEG_1_8, acosEstimate,
        [within('acos', 0, 1, Math.acos, getACosEstimateMaxError)]);
    maxErrorCase(family, 'ACosEstimate.getMaxError', DEG_1_8, getACosEstimateMaxError);
    estimateCase(family, 'ASinEstimate.estimate', DEG_1_8, asinEstimate,
        [within('asin', 0, 1, Math.asin, getASinEstimateMaxError)]);
    maxErrorCase(family, 'ASinEstimate.getMaxError', DEG_1_8, getASinEstimateMaxError);

    // SqrtEstimate.h and InvSqrtEstimate.h. RR relative bounds: sqrt(y) >= 1
    // on [1, 2) (factor 1) and 1/sqrt(y) > 1/sqrt(2) (factor sqrt(2)).
    estimateCase(family, 'SqrtEstimate.estimate', DEG_1_8, sqrtEstimate,
        [within('sqrt', 1, 2, Math.sqrt, getSqrtEstimateMaxError)]);
    estimateCase(family, 'SqrtEstimate.estimateRR', DEG_1_8, sqrtEstimateRR,
        [withinRelative('sqrtRR', normal, Math.sqrt, getSqrtEstimateMaxError, 1)],
        frexpBranch);
    maxErrorCase(family, 'SqrtEstimate.getMaxError', DEG_1_8, getSqrtEstimateMaxError);
    estimateCase(family, 'InvSqrtEstimate.estimate', DEG_1_8, invSqrtEstimate,
        [within('invSqrt', 1, 2, (x) => 1 / Math.sqrt(x), getInvSqrtEstimateMaxError)]);
    estimateCase(family, 'InvSqrtEstimate.estimateRR', DEG_1_8, invSqrtEstimateRR,
        [withinRelative('invSqrtRR', normal, (x) => 1 / Math.sqrt(x), getInvSqrtEstimateMaxError,
            Math.SQRT2)],
        frexpBranch);
    maxErrorCase(family, 'InvSqrtEstimate.getMaxError', DEG_1_8, getInvSqrtEstimateMaxError);

    // Floor(x) >= 2^31: upstream's static_cast<int32_t> is undefined
    // behaviour (MSVC: INT_MIN, so the result is +0); the port keeps the
    // double exponent and returns +inf, which is 2^x rounded. New upstream
    // suspect, see oracle/reports/v23-estimates.md.
    family.case('Exp2Estimate.estimateRR.hugeArgument', (io) => {
        const x = io.real();
        for (const d of DEG_1_7) { io.outReal(exp2EstimateRR(x, d)); }
        for (const d of DEG_1_7) { io.outReal(expEstimateRR(x, d)); }
    }, { exact: true, deviation: 'v23 report: Exp2EstimateRR int32_t conversion of floor(x) >= 2^31' });

    // ---- ChebyshevRatio.h (std::sin, std::acos) ---------------------------

    // Identity check f(t, A) sin(A) = sin(t A), with Math.sin.
    const ratioRef = (key: string, t: number, angle: number, f: number): void => {
        if (!(angle > 0 && angle < GTE_C_PI)) { return; }
        const lhs = f * Math.sin(angle);
        const rhs = Math.sin(t * angle);
        refCheck(key, lhs - rhs, 1e-13 * Math.max(1, Math.abs(rhs)), 0, () => `${key}(${t}, ${angle})`);
    };

    family.case('ChebyshevRatio.ratio', (io) => {
        for (let m = 0; m < NUM_MODES; ++m) {
            const t = io.real();
            const angle = io.real();
            const f = chebyshevRatio(t, angle);
            if (angle === 0) {
                count('ChebyshevRatio.ratio', 'angle == 0 (returns t)');
                io.outRealExact(f);
            } else {
                count('ChebyshevRatio.ratio', 'angle in (0, pi)');
                io.outReal(f);
                ratioRef('ChebyshevRatio identity', t, angle, f);
            }
        }
    });

    family.case('ChebyshevRatio.ratios', (io) => {
        for (let m = 0; m < NUM_MODES; ++m) {
            const t = io.real();
            const angle = io.real();
            const f = chebyshevRatios(t, angle);
            if (angle === 0) {
                io.outRealExact(f[0]);
                io.outRealExact(f[1]);
            } else {
                io.outReal(f[0]);
                io.outReal(f[1]);
                ratioRef('ChebyshevRatios identity', 1 - t, angle, f[0]);
                ratioRef('ChebyshevRatios identity', t, angle, f[1]);
            }
        }
    });

    const cosBranch = (c: number): string => (c < 1
        ? (c > -1 ? 'cosA in (-1, 1)' : 'cosA <= -1 (throws)')
        : (Number.isNaN(c) ? 'NaN (returns t)' : 'cosA >= 1 (returns t)'));

    family.case('ChebyshevRatio.ratioUsingCosAngle', (io) => {
        for (let m = 0; m < NUM_MODES; ++m) {
            const t = io.real();
            const c = io.real();
            count('ChebyshevRatio.ratioUsingCosAngle', cosBranch(c));
            const f = chebyshevRatioUsingCosAngle(t, c);
            if (c < 1) {
                io.outReal(f);
                ratioRef('ChebyshevRatioUsingCosAngle identity', t, Math.acos(c), f);
            } else {
                io.outRealExact(f);
            }
        }
    });

    family.case('ChebyshevRatio.ratiosUsingCosAngle', (io) => {
        for (let m = 0; m < NUM_MODES; ++m) {
            const t = io.real();
            const c = io.real();
            const f = chebyshevRatiosUsingCosAngle(t, c);
            if (c < 1) {
                io.outReal(f[0]);
                io.outReal(f[1]);
            } else {
                io.outRealExact(f[0]);
                io.outRealExact(f[1]);
            }
        }
    });

    // Conditioning: pi - A = gap ~ sqrt(2 (1 + c)) (1 + c is exact), and a
    // 1-ulp(pi) = 2^-51 difference in acos moves the ratio by 2^-51/gap
    // relative. Tolerance: 4 such ulps, never below the default 1e-12.
    family.case('ChebyshevRatio.usingCosAngle.nearPi', (io) => {
        for (let m = 0; m < NUM_MODES; ++m) {
            const t = io.real();
            const c = io.real();
            const gap = Math.sqrt(2 * (1 + c));
            const tol = Math.max(1e-12, 4 * 2 ** -51 / gap);
            io.outReal(chebyshevRatioUsingCosAngle(t, c), tol);
            const f = chebyshevRatiosUsingCosAngle(t, c);
            io.outReal(f[0], tol);
            io.outReal(f[1], tol);
        }
    });

    family.case('ChebyshevRatio.throwParity', (io) => {
        const which = io.integer();
        const t = io.real();
        const arg = io.real();
        let f: [number, number] = [0, 0];
        switch (which) {
            case 0: f[0] = chebyshevRatio(t, arg); break;
            case 1: f = chebyshevRatios(t, arg); break;
            case 2: f[0] = chebyshevRatioUsingCosAngle(t, arg); break;
            default: f = chebyshevRatiosUsingCosAngle(t, arg); break;
        }
        io.outReal(f[0]);
        io.outReal(f[1]);
    });

    // ---- ChebyshevRatioEstimate.h (exact) ----------------------------------

    // Within the published bound of the true ratio sin(tA)/sin(A), A =
    // acos(x), for t in [0, 1] and x in [xmin, 1].
    const chbRef = (key: string, t: number, x: number, xmin: number, d: number,
        f: [number, number], bound: number): void => {
        if (!(t >= 0 && t <= 1 && x >= xmin && x <= 1)) { return; }
        const A = Math.acos(x);
        const truth = (s: number): number => (A === 0 ? s : Math.sin(s * A) / Math.sin(A));
        const what = (): string => `${key}(${t}, ${x}, degree ${d})`;
        refCheck(`${key} degree ${d}`, f[0] - truth(1 - t), bound, 1e-15, what);
        refCheck(`${key} degree ${d}`, f[1] - truth(t), bound, 1e-15, what);
    };

    family.case('ChebyshevRatioEstimate.estimate', (io) => {
        for (let m = 0; m < NUM_MODES; ++m) {
            const t = io.real();
            const x = io.real();
            for (const d of DEG_1_16) {
                const f = chebyshevRatioEstimate(t, x, d);
                io.outReal(f[0]);
                io.outReal(f[1]);
                chbRef('ChebyshevRatioEstimate', t, x, 0, d, f, getChebyshevRatioEstimateMaxError(d));
            }
        }
    }, { exact: true });

    family.case('ChebyshevRatioEstimate.estimateR', (io) => {
        for (let m = 0; m < NUM_MODES; ++m) {
            const t = io.real();
            const x = io.real();
            for (const d of DEG_1_12) {
                const f = chebyshevRatioEstimateR(t, x, d);
                io.outReal(f[0]);
                io.outReal(f[1]);
                chbRef('ChebyshevRatioEstimateR', t, x, GTE_C_INV_SQRT_2, d, f,
                    getChebyshevRatioEstimateRMaxError(d));
            }
        }
    }, { exact: true });

    maxErrorCase(family, 'ChebyshevRatioEstimate.getMaxError', DEG_1_16,
        getChebyshevRatioEstimateMaxError);
    maxErrorCase(family, 'ChebyshevRatioEstimate.getMaxErrorR', DEG_1_12,
        getChebyshevRatioEstimateRMaxError);

    // ---- Slerp.h (std::acos, std::sin) and SlerpEstimate.h (exact) ---------

    const DIMS = [2, 3, 4];
    const dot = (a: readonly number[], b: readonly number[]): number => {
        let d = 0;
        for (let i = 0; i < a.length; ++i) { d += a[i] * b[i]; }
        return d;
    };
    const norm = (a: readonly number[]): number => Math.sqrt(a.reduce((s, v) => s + v * v, 0));
    // The angle between unit vectors, accurate for all angles.
    const angleBetween = (a: readonly number[], b: readonly number[]): number => {
        const dif = a.map((v, i) => v - b[i]);
        const sum = a.map((v, i) => v + b[i]);
        return 2 * Math.atan2(norm(dif), norm(sum));
    };
    const trueSlerp = (t: number, a: readonly number[], b: readonly number[]): number[] => {
        const A = angleBetween(a, b);
        if (A === 0) { return a.map((v, i) => (1 - t) * v + t * b[i]); }
        const s = Math.sin(A);
        return a.map((v, i) => (Math.sin((1 - t) * A) * v + Math.sin(t * A) * b[i]) / s);
    };

    // Emit a slerp result: bit identity where the output is decided by
    // arithmetic alone (the angle-0 branch, or a slot where both inputs of
    // the combination are +-0), the libm tolerance elsewhere.
    const emitSlerp = (io: OracleIO, r: readonly number[], exactAll: boolean,
        a: readonly number[], b: readonly number[]): void => {
        for (let i = 0; i < r.length; ++i) {
            if (exactAll || (a[i] === 0 && b[i] === 0)) { io.outRealExact(r[i]); } else { io.outReal(r[i]); }
        }
    };

    // Unit length, and the angle from q0 is t times the angle q0-q1.
    const slerpRef = (key: string, t: number, q0: readonly number[], q1: readonly number[],
        r: readonly number[]): void => {
        if (!(t >= 0 && t <= 1)) { return; }
        const A = angleBetween(q0, q1);
        const what = (): string => `${key}(t = ${t}, q0 = ${q0}, q1 = ${q1})`;
        refCheck(`${key} |r| - 1`, norm(r) - 1, 1e-12, 0, what);
        refCheck(`${key} angle(q0, r) - tA`, angleBetween(q0, r) - t * A, 1e-12, 0, what);
        refCheck(`${key} angle(r, q1) - (1-t)A`, angleBetween(r, q1) - (1 - t) * A, 1e-12, 0, what);
    };

    family.case('Slerp.slerp', (io) => {
        for (const n of DIMS) {
            const q0 = io.reals(n);
            const q1 = io.reals(n);
            const t = io.real();
            const c = dot(q0, q1);
            count('Slerp.slerp', cosBranch(c));
            const r = slerp(t, q0, q1);
            emitSlerp(io, r, !(c < 1), q0, q1);
            if (c < 1) { slerpRef('Slerp', t, q0, q1, r); }
        }
    });

    family.case('Slerp.slerpCosAngle', (io) => {
        for (const n of DIMS) {
            const q0 = io.reals(n);
            const q1 = io.reals(n);
            const t = io.real();
            const c = io.real();
            count('Slerp.slerpCosAngle', cosBranch(c));
            const r = slerpUsingCosAngle(t, q0, q1, c);
            emitSlerp(io, r, !(c < 1), q0, q1);
            if (c < 1) { slerpRef('SlerpCosAngle', t, q0, q1, r); }
        }
    });

    family.case('Slerp.slerpMidpoint', (io) => {
        for (const n of DIMS) {
            const q0 = io.reals(n);
            const q1 = io.reals(n);
            const t = io.real();
            const qh = io.reals(n);
            const cosAH = io.real();
            const first = 2 * t <= 1;
            count('Slerp.slerpMidpoint', `${cosBranch(cosAH)}, ${first ? '2t <= 1' : '2t > 1'}`);
            if (dot(q0, q1) < -0.99) { count('Slerp.slerpMidpoint', 'dot(q0, q1) < -0.99'); }
            const r = slerpUsingMidpoint(t, q0, q1, qh, cosAH);
            emitSlerp(io, r, !(cosAH < 1), first ? q0 : qh, first ? qh : q1);
            if (cosAH < 1 && cosAH >= 1e-3) { slerpRef('SlerpMidpoint', t, q0, q1, r); }
        }
    });

    family.case('Slerp.throwParity', (io) => {
        const which = io.integer();
        io.integer();  // 'valid', which only steers the C++ generator
        const q0 = io.reals(4);
        const q1 = io.reals(4);
        const t = io.real();
        const c = io.real();
        const qh = io.reals(4);
        let r: number[];
        switch (which) {
            case 0: r = slerp(t, q0, q1); break;
            case 1: r = slerpUsingCosAngle(t, q0, q1, c); break;
            default: r = slerpUsingMidpoint(t, q0, q1, qh, c); break;
        }
        io.outReals(r);
    });

    // Within maxError(D) of the true slerp, per component scaled by the two
    // combined inputs (the estimate replaces both Chebyshev ratios).
    const slerpEstRef = (key: string, d: number, t: number, a: readonly number[],
        b: readonly number[], r: readonly number[]): void => {
        const truth = trueSlerp(t, a, b);
        const bound = getChebyshevRatioEstimateMaxError(d);
        for (let i = 0; i < r.length; ++i) {
            refCheck(`${key} degree ${d}`, (r[i] - truth[i]) / (Math.abs(a[i]) + Math.abs(b[i]) || 1),
                bound, 1e-12, () => `${key}(t = ${t}, a = ${a}, b = ${b})`);
        }
    };

    family.case('SlerpEstimate.slerpEstimate', (io) => {
        for (const n of DIMS) {
            const q0 = io.reals(n);
            const q1 = io.reals(n);
            const t = io.real();
            const c = dot(q0, q1);
            count('SlerpEstimate.slerpEstimate', c >= 0 ? 'angle <= pi/2' : 'angle > pi/2 (outside the documented range)');
            for (const d of DEG_1_16) {
                const r = slerpEstimate(t, q0, q1, d);
                io.outReals(r);
                if (c >= 0 && t >= 0 && t <= 1) { slerpEstRef('SlerpEstimate', d, t, q0, q1, r); }
            }
        }
    }, { exact: true });

    family.case('SlerpEstimate.slerpEstimateCosAngle', (io) => {
        for (const n of DIMS) {
            const q0 = io.reals(n);
            const q1 = io.reals(n);
            const t = io.real();
            const c = io.real();
            for (const d of DEG_1_16) {
                const r = slerpEstimateUsingCosAngle(t, q0, q1, c, d);
                io.outReals(r);
                if (c >= 0 && t >= 0 && t <= 1) { slerpEstRef('SlerpEstimateCosAngle', d, t, q0, q1, r); }
            }
        }
    }, { exact: true });

    family.case('SlerpEstimate.slerpEstimateMidpoint', (io) => {
        for (const n of DIMS) {
            const q0 = io.reals(n);
            const q1 = io.reals(n);
            const t = io.real();
            const qh = io.reals(n);
            const cosAH = io.real();
            count('SlerpEstimate.slerpEstimateMidpoint', 2 * t <= 1 ? '2t <= 1' : '2t > 1');
            for (const d of DEG_1_16) {
                const r = slerpEstimateUsingMidpoint(t, q0, q1, qh, cosAH, d);
                io.outReals(r);
                if (cosAH >= 1e-3 && t >= 0 && t <= 1) {
                    const first = 2 * t <= 1;
                    slerpEstRef('SlerpEstimateMidpoint', d, first ? 2 * t : 2 * t - 1,
                        first ? q0 : qh, first ? qh : q1, r);
                }
            }
        }
    }, { exact: true });

    family.finish();

    // ---- independent reference checks and branch coverage -----------------

    it('every estimate is within its published bound of the true function', () => {
        const lines: string[] = [];
        const violations: string[] = [];
        for (const [key, s] of [...refStats.entries()].sort()) {
            lines.push(`${key}: ${s.checked} checked, worst |error|/bound ${s.worst.toPrecision(3)}`);
            violations.push(...s.violations.map((v) => `${key}: ${v}`));
        }
        console.log(`v23 reference checks\n  ${lines.join('\n  ')}`);
        expect(refStats.size).toBeGreaterThan(0);
        expect(violations).toEqual([]);
    });

    it('every range-reduction and ratio branch is reached', () => {
        const lines: string[] = [];
        for (const [key, h] of [...histograms.entries()].sort()) {
            const parts = [...h.entries()].sort().map(([b, n]) => `${b}: ${n}`);
            lines.push(`${key}: ${parts.join(', ')}`);
        }
        console.log(`v23 branch histograms\n  ${lines.join('\n  ')}`);
        const required: Record<string, string[]> = {
            'Exp2Estimate.estimateRR': ['normal', 'normal y=0', 'subnormal', 'underflow to 0',
                'overflow', 'non-finite'],
            'ExpEstimate.estimateRR': ['normal', 'subnormal', 'underflow to 0', 'overflow', 'non-finite'],
            'Log2Estimate.estimateRR': ['p even', 'p odd', 'subnormal p even', 'subnormal p odd',
                '+0', '-0', 'non-finite'],
            'SqrtEstimate.estimateRR': ['p even', 'p odd', '+0', '-0'],
            'InvSqrtEstimate.estimateRR': ['p even', 'p odd'],
            'SinEstimate.estimateRR': ['r>t', 'r<-t', '|r|<=t', 'r>t tie', 'r<-t tie', 'non-finite'],
            'CosEstimate.estimateRR': ['r>t', 'r<-t', '|r|<=t', 'r>t tie', 'r<-t tie'],
            'TanEstimate.estimateRR': ['|y|<=pi/4', 'y>pi/4', 'y<-pi/4', 'y>pi/4 tie', 'y<-pi/4 tie'],
            'ATanEstimate.estimateRR': ['|x|<=1', 'x>1', 'x<-1', 'NaN (else arm)'],
            'ChebyshevRatio.ratio': ['angle == 0 (returns t)', 'angle in (0, pi)'],
            'ChebyshevRatio.ratioUsingCosAngle': ['cosA in (-1, 1)', 'cosA >= 1 (returns t)',
                'NaN (returns t)'],
            'Slerp.slerp': ['cosA in (-1, 1)', 'cosA >= 1 (returns t)'],
            'Slerp.slerpMidpoint': ['cosA in (-1, 1), 2t <= 1', 'cosA in (-1, 1), 2t > 1',
                'dot(q0, q1) < -0.99'],
            'SlerpEstimate.slerpEstimate': ['angle <= pi/2', 'angle > pi/2 (outside the documented range)'],
            'SlerpEstimate.slerpEstimateMidpoint': ['2t <= 1', '2t > 1']
        };
        const missing: string[] = [];
        for (const [key, branches] of Object.entries(required)) {
            for (const b of branches) {
                if (!(histograms.get(key)?.has(b) ?? false)) { missing.push(`${key}: ${b}`); }
            }
        }
        // TanEstimateRR's r > pi/2 / r < -pi/2 arms are dead (#57).
        const dead = [...(histograms.get('TanEstimate.estimateRR')?.keys() ?? [])]
            .filter((b) => b.includes('UNREACHABLE'));
        expect(dead).toEqual([]);
        expect(missing).toEqual([]);
    });
});
