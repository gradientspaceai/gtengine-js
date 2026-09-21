// Regression test for the back-substitution accumulation order of
// NaturalQuinticSpline (found by the C++ oracle, group 17 / v17-curves).
//
// Upstream NaturalQuinticSpline.h writes each polynomial coefficient as
//
//   poly[k] = invR(r,0)*B[j0] + invR(r,1)*B[j1] + invR(r,2)*B[j2]
//           + invR(r,3)*B[j3];
//
// which C++ evaluates strictly left to right, ((a+b)+c)+d. The port grouped
// the four terms pairwise, (a+b)+(c+d). The two differ in the last bits
// whenever B[j2] and B[j3] are both nonzero, which is the case for every
// closed and every clamped spline (a free spline has B[j2] = 0, which is why
// the defect was invisible there). The error is then amplified by the
// backward recursion over the segments: the C++ oracle measured scaled
// differences up to 2.3e-13 on the polynomial coefficients and the jets.
//
// The expected values below were produced by the real MSVC build of upstream
// (oracle/golden/v17-curves.txt, cases NaturalQuinticSpline.closed and
// NaturalQuinticSpline.clamped); 'PAIRWISE' holds what the port produced
// before the fix and is asserted to be different, so the test fails if the
// grouping is ever reintroduced.
import { describe, expect, it } from 'vitest';
import { NaturalQuinticSpline } from '../src/NaturalQuinticSpline.js';
import { Vector } from '../src/Vector.js';

function fromHex(hex: string): number {
    const view = new DataView(new ArrayBuffer(8));
    view.setBigUint64(0, BigInt('0x' + hex));
    return view.getFloat64(0);
}

// Every coefficient of every segment, as "<x>,<y>" hex pairs in the order
// getPolynomials() returns them.
const UPSTREAM_CLOSED = [
    '0000000000000000,0000000000000000', '3fe0000000000000,bff0000000000000',
    '4000452e1018a880,401ca9de1850d2f4', '4005bc884e159370,401a94a25546b8a0',
    'c014244d663a902c,c03049b7bcdffa88', '3ffc8dc8dc8dc8dc,4019e85e85e85e88',
    '4000000000000000,4008000000000000', '400e000000000000,3ff4000000000000',
    'c02b45231a8ab3cc,c0424dc0a4635156', '401a5679941ed290,4058768ab3cca0f6',
    '3fce3beee0523400,c0561874710447eb', '3fb07b4a7d00d900,403ae3283da53e80',
    'bff0000000000000,4010000000000000', 'c004000000000000,4004000000000000',
    '40205bd37a6f4de7,c011e1642c8590ae', '4020589d89d89d80,402a64ec4ec4ec52',
    'c03b225ac17f9260,c0443df096b05fe6', '402c90447eb73958,403941c4111fadc8'
];

const PAIRWISE_CLOSED = [
    '0000000000000000,0000000000000000', '3fe0000000000000,bff0000000000000',
    '4000452e1018a880,401ca9de1850d2f8', '4005bc884e159370,401a94a25546b8a0',
    'c014244d663a902c,c03049b7bcdffa84', '3ffc8dc8dc8dc8dc,4019e85e85e85e88',
    '4000000000000000,4008000000000000', '400e000000000000,3ff4000000000000',
    'c02b45231a8ab3cc,c0424dc0a4635155', '401a5679941ed290,4058768ab3cca0f5',
    '3fce3beee0523400,c0561874710447ea', '3fb07b4a7d00d900,403ae3283da53e7f',
    'bff0000000000000,4010000000000000', 'c004000000000000,4004000000000000',
    '40205bd37a6f4de7,c011e1642c8590ae', '4020589d89d89d80,402a64ec4ec4ec50',
    'c03b225ac17f925f,c0443df096b05fe6', '402c90447eb73958,403941c4111fadc9'
];

const UPSTREAM_CLAMPED = [
    '0000000000000000,0000000000000000', '3fe0000000000000,bff0000000000000',
    '3fc0000000000180,3fcffffffffffd80', '401fc55e2455e228,40399ae3ec2e3ed0',
    'c023855e2455e236,c040bae3ec2e3ec8', '40098abc48abc47c,402735c7d85c7d8c',
    '4000000000000000,4008000000000000', '400e000000000000,3ff4000000000000',
    'c030c39ce739ce6f,c047e3f7bdef7be1', '4029e1ce739ce71c,405e871ef7bdef7f',
    'c008b7bdef7bde80,c059d84a5294a52c', '3fc4d6b5ad6b5880,403ccc9ce739ce77',
    'bff0000000000000,4010000000000000', 'c004000000000000,4004000000000000',
    '401498421084210b,c02b76a5294a5296', '4001ae739ce739b6,c00b5042108420f6',
    'c0187739ce739cde,4034860842108422', '40020f7bdef7bdec,c023c15ad6b5ad6c'
];

const PAIRWISE_CLAMPED = [
    '0000000000000000,0000000000000000', '3fe0000000000000,bff0000000000000',
    '3fc0000000000180,3fcffffffffffe00', '401fc55e2455e228,40399ae3ec2e3ecc',
    'c023855e2455e236,c040bae3ec2e3ec8', '40098abc48abc47c,402735c7d85c7d8c',
    '4000000000000000,4008000000000000', '400e000000000000,3ff4000000000000',
    'c030c39ce739ce6f,c047e3f7bdef7be0', '4029e1ce739ce71c,405e871ef7bdef7f',
    'c008b7bdef7bde80,c059d84a5294a52c', '3fc4d6b5ad6b5880,403ccc9ce739ce77',
    'bff0000000000000,4010000000000000', 'c004000000000000,4004000000000000',
    '401498421084210b,c02b76a5294a5296', '4001ae739ce739b6,c00b5042108420f8',
    'c0187739ce739cde,4034860842108422', '40020f7bdef7bdec,c023c15ad6b5ad6c'
];

const f0 = [
    Vector.fromArray([0, 0]), Vector.fromArray([2, 3]),
    Vector.fromArray([-1, 4]), Vector.fromArray([0, 0])
];
const f1 = [
    Vector.fromArray([1, -2]), Vector.fromArray([3, 1]),
    Vector.fromArray([-2, 2]), Vector.fromArray([1, -2])
];
const times = [0, 0.5, 1.75, 3];

function coefficients(spline: NaturalQuinticSpline): number[][] {
    const out: number[][] = [];
    for (const poly of spline.getPolynomials()) {
        for (const c of poly) { out.push([c.values[0], c.values[1]]); }
    }
    return out;
}

function expected(hexes: readonly string[]): number[][] {
    return hexes.map((h) => h.split(',').map(fromHex));
}

describe('NaturalQuinticSpline back-substitution accumulation order', () => {
    it('matches upstream bit for bit on a closed spline', () => {
        const actual = coefficients(
            NaturalQuinticSpline.createClosed(f0, f1, times));
        const want = expected(UPSTREAM_CLOSED);
        expect(actual.length).toBe(want.length);
        for (let i = 0; i < want.length; ++i) {
            for (let k = 0; k < 2; ++k) {
                expect(Object.is(actual[i][k], want[i][k]),
                    `coefficient ${i} component ${k}`).toBe(true);
            }
        }
    });

    it('matches upstream bit for bit on a clamped spline', () => {
        const actual = coefficients(
            NaturalQuinticSpline.createClamped(f0, f1, times,
                Vector.fromArray([1, 2]), Vector.fromArray([-3, 0.5])));
        const want = expected(UPSTREAM_CLAMPED);
        expect(actual.length).toBe(want.length);
        for (let i = 0; i < want.length; ++i) {
            for (let k = 0; k < 2; ++k) {
                expect(Object.is(actual[i][k], want[i][k]),
                    `coefficient ${i} component ${k}`).toBe(true);
            }
        }
    });

    it('the pairwise grouping is a different computation', () => {
        // The two candidate groupings of the same four-term sum are not the
        // same function on these inputs, so the bit-identity above is
        // evidence and not a coincidence.
        for (const [upstream, pairwise] of
            [[UPSTREAM_CLOSED, PAIRWISE_CLOSED],
                [UPSTREAM_CLAMPED, PAIRWISE_CLAMPED]]) {
            const a = expected(upstream);
            const b = expected(pairwise);
            let differing = 0;
            for (let i = 0; i < a.length; ++i) {
                for (let k = 0; k < 2; ++k) {
                    if (!Object.is(a[i][k], b[i][k])) { ++differing; }
                }
            }
            expect(differing).toBeGreaterThan(0);
        }
    });
});
