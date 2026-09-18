/// <reference types="node" />
// Replays golden records produced by the C++ oracle (oracle/cpp, built
// against the upstream GTE headers with MSVC) through the TypeScript port and
// compares outputs. See ORACLE.md.
//
// A case body mirrors its C++ ORACLE_CASE: it pulls inputs from 'io' in the
// same order the C++ side generated them, runs the port, and pushes outputs
// in the same order the C++ side recorded them.

import { existsSync, mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { expect, it } from 'vitest';
import { Matrix } from '../../src/Matrix.js';
import { Vector } from '../../src/Vector.js';

const repoRoot = join(dirname(fileURLToPath(import.meta.url)), '..', '..');

// ORACLE_GOLDEN_DIR points the tests at a deep local run (oracle/out/deep)
// instead of the committed goldens.
const goldenDir = process.env['ORACLE_GOLDEN_DIR']
    ?? join(repoRoot, 'oracle', 'golden');

// ORACLE_STATS=1 writes per-family agreement statistics to oracle/out/stats.
const statsDir = process.env['ORACLE_STATS']
    ? join(repoRoot, 'oracle', 'out', 'stats') : undefined;

interface GoldenRecord {
    inputs: Float64Array;
    // null when the C++ code threw.
    outputs: Float64Array | null;
}

export interface CaseOptions {
    // Default tolerance for real outputs: |a - e| <= tol * max(1, |a|, |e|).
    // Default 1e-12. Use 0 with 'exact'.
    tol?: number;
    // Require every real output to be bit-identical to the C++ value
    // (+0 and -0 are distinguished, any NaN matches any NaN).
    exact?: boolean;
    // Skip the case with a documented reason. The reason is reported.
    skip?: string;
    // The port is EXPECTED to disagree with the C++ build on at least one
    // record because it deliberately fixes an upstream defect. The string
    // cites the record of that decision (issue number, UPSTREAM-FINDINGS
    // entry). The test fails when every record agrees, which means the
    // deviation no longer exists and the case should become a normal one.
    deviation?: string;
}

const DEFAULT_TOL = 1e-12;

function parseHex(tokens: string[], start: number): Float64Array {
    const out = new Float64Array(tokens.length - start);
    const view = new DataView(out.buffer);
    for (let i = start; i < tokens.length; ++i) {
        view.setBigUint64((i - start) * 8, BigInt('0x' + tokens[i]), true);
    }
    return out;
}

function toHex(x: number): string {
    const view = new DataView(new ArrayBuffer(8));
    view.setFloat64(0, x);
    return view.getBigUint64(0).toString(16).padStart(16, '0');
}

function sameBits(a: number, b: number): boolean {
    return Object.is(a, b) || (Number.isNaN(a) && Number.isNaN(b));
}

class Mismatch extends Error {}

export class OracleIO {
    private inPos = 0;
    private outPos = 0;
    readonly actual: number[] = [];
    // Value disagreements; the record keeps going so that every differing
    // output is reported. Structural disagreements throw Mismatch instead.
    readonly valueFailures: string[] = [];

    // Statistics over real outputs of this record.
    numReal = 0;
    numExact = 0;
    maxRelErr = 0;

    constructor(readonly index: number,
        private readonly record: GoldenRecord,
        private readonly options: CaseOptions) {}

    // ---- inputs ----
    real(): number {
        if (this.inPos >= this.record.inputs.length) {
            throw new Mismatch('the port consumed more inputs than C++ generated');
        }
        return this.record.inputs[this.inPos++];
    }

    integer(): number {
        const x = this.real();
        if (!Number.isInteger(x)) {
            throw new Mismatch(`input ${this.inPos - 1} is not an integer: ${x}`);
        }
        return x;
    }

    boolean(): boolean { return this.integer() !== 0; }

    vec(n: number): Vector {
        const v = new Vector(n);
        for (let i = 0; i < n; ++i) { v.set(i, this.real()); }
        return v;
    }

    reals(n: number): number[] {
        const a: number[] = [];
        for (let i = 0; i < n; ++i) { a.push(this.real()); }
        return a;
    }

    // ---- outputs ----
    private expected(what: string): number {
        const outputs = this.record.outputs;
        if (outputs === null) {
            throw new Mismatch(`C++ threw, but the port produced ${what}`);
        }
        if (this.outPos >= outputs.length) {
            throw new Mismatch(`the port produced more outputs than C++ (${outputs.length})`);
        }
        return outputs[this.outPos++];
    }

    outReal(x: number, tol?: number): void {
        const e = this.expected('a real');
        this.actual.push(x);
        ++this.numReal;
        if (sameBits(x, e)) {
            ++this.numExact;
            return;
        }
        const scale = Math.max(1, Math.abs(x), Math.abs(e));
        const err = Math.abs(x - e) / scale;
        if (Number.isFinite(err)) { this.maxRelErr = Math.max(this.maxRelErr, err); }
        const t = this.options.exact ? 0 : (tol ?? this.options.tol ?? DEFAULT_TOL);
        const ok = !this.options.exact && (x === e || err <= t);
        if (!ok) {
            this.valueFailures.push(`output ${this.outPos - 1}: port ${x} (${toHex(x)}) vs `
                + `C++ ${e} (${toHex(e)}), scaled error ${err.toExponential(3)}, `
                + (this.options.exact ? 'exact match required' : `tolerance ${t}`));
        }
    }

    // Discrete outputs always compare exactly.
    outInt(i: number): void {
        const e = this.expected('an integer');
        this.actual.push(i);
        if (i !== e) {
            throw new Mismatch(`output ${this.outPos - 1}: port integer ${i} vs C++ ${e}`);
        }
    }

    outBool(b: boolean): void {
        const e = this.expected('a boolean');
        this.actual.push(b ? 1 : 0);
        if ((b ? 1 : 0) !== e) {
            throw new Mismatch(`output ${this.outPos - 1}: port boolean ${b} vs C++ ${e !== 0}`);
        }
    }

    outVec(v: Vector, tol?: number): void {
        for (let i = 0; i < v.size; ++i) { this.outReal(v.get(i), tol); }
    }

    outReals(a: readonly number[], tol?: number): void {
        for (const x of a) { this.outReal(x, tol); }
    }

    // Row-major, matching Ctx::outMat.
    outMat(m: Matrix, tol?: number): void {
        for (let r = 0; r < m.numRows; ++r) {
            for (let c = 0; c < m.numCols; ++c) { this.outReal(m.get(r, c), tol); }
        }
    }

    // Size first, matching Ctx::outGVec.
    outGVec(v: Vector, tol?: number): void {
        this.outInt(v.size);
        this.outVec(v, tol);
    }

    // Rows, columns, elements, matching Ctx::outGMat.
    outGMat(m: Matrix, tol?: number): void {
        this.outInt(m.numRows);
        this.outInt(m.numCols);
        this.outMat(m, tol);
    }

    finish(): void {
        if (this.inPos !== this.record.inputs.length) {
            throw new Mismatch(`the port consumed ${this.inPos} of `
                + `${this.record.inputs.length} inputs`);
        }
        const outputs = this.record.outputs;
        if (outputs !== null && this.outPos !== outputs.length) {
            throw new Mismatch(`the port produced ${this.outPos} of ${outputs.length} outputs`);
        }
    }
}

interface CaseStats {
    records: number;
    threw: number;
    realOutputs: number;
    exactOutputs: number;
    maxRelErr: number;
    exactRequired: boolean;
    skipped?: string;
    deviation?: string;
    deviatingRecords?: number;
}

export class OracleFamily {
    private readonly cases = new Map<string, GoldenRecord[]>();
    private readonly claimed = new Set<string>();
    private readonly stats: Record<string, CaseStats> = {};
    readonly header: string[] = [];

    constructor(readonly family: string) {
        const path = join(goldenDir, `${family}.txt`);
        if (!existsSync(path)) {
            throw new Error(`golden file not found: ${path} (run npm run oracle:gen)`);
        }
        let current: GoldenRecord[] | undefined;
        let pending: Float64Array | undefined;
        for (const line of readFileSync(path, 'utf8').split('\n')) {
            if (line.length === 0) { continue; }
            if (line.startsWith('#')) { this.header.push(line); continue; }
            const tokens = line.split(' ');
            if (tokens[0] === 'case') {
                current = [];
                this.cases.set(line.substring(5), current);
            } else if (tokens[0] === 'i') {
                pending = parseHex(tokens, 1);
            } else if (tokens[0] === 'o') {
                current!.push({
                    inputs: pending!,
                    outputs: tokens[1] === '!' ? null : parseHex(tokens, 1)
                });
            } else {
                throw new Error(`${path}: unrecognized line '${line.substring(0, 40)}'`);
            }
        }
    }

    // Register a vitest test that replays every golden record of 'name'.
    case(name: string, body: (io: OracleIO) => void, options: CaseOptions = {}): void {
        const records = this.cases.get(name);
        this.claimed.add(name);
        const run = options.skip !== undefined ? it.skip : it;
        if (options.skip !== undefined) {
            this.stats[name] = {
                records: records?.length ?? 0, threw: 0, realOutputs: 0, exactOutputs: 0,
                maxRelErr: 0, exactRequired: false, skipped: options.skip
            };
        }
        run(`${name} matches the C++ oracle`, () => {
            expect(records, `no golden records for case '${name}'`).toBeDefined();
            const stats: CaseStats = {
                records: records!.length, threw: 0, realOutputs: 0, exactOutputs: 0,
                maxRelErr: 0, exactRequired: options.exact === true
            };
            const failures: string[] = [];
            records!.forEach((record, index) => {
                const io = new OracleIO(index, record, options);
                try {
                    let threw = false;
                    try {
                        body(io);
                    } catch (e) {
                        if (e instanceof Mismatch) { throw e; }
                        threw = true;
                        if (record.outputs !== null) {
                            throw new Mismatch(`the port threw (${(e as Error).message}) `
                                + 'but C++ returned normally');
                        }
                    }
                    if (record.outputs === null) {
                        ++stats.threw;
                        if (!threw) { throw new Mismatch('C++ threw but the port returned normally'); }
                    } else {
                        io.finish();
                    }
                    if (io.valueFailures.length > 0) {
                        throw new Mismatch(io.valueFailures.join('; '));
                    }
                } catch (e) {
                    if (!(e instanceof Mismatch)) { throw e; }
                    failures.push(`record ${index}: ${e.message}`);
                }
                stats.realOutputs += io.numReal;
                stats.exactOutputs += io.numExact;
                stats.maxRelErr = Math.max(stats.maxRelErr, io.maxRelErr);
            });
            this.stats[name] = stats;
            if (options.deviation !== undefined) {
                stats.deviation = options.deviation;
                stats.deviatingRecords = failures.length;
                expect(failures.length, `'${name}' is declared as a deliberate deviation `
                    + `(${options.deviation}) but every record agrees with C++`)
                    .toBeGreaterThan(0);
                return;
            }
            if (failures.length > 0) {
                const shown = failures.slice(0, 8).join('\n  ');
                expect.fail(`${failures.length} of ${records!.length} records disagree with `
                    + `the C++ oracle:\n  ${shown}`);
            }
        });
    }

    // Call last: fails when the golden file has cases no test claimed, and
    // writes the statistics when ORACLE_STATS is set.
    finish(): void {
        it('every golden case has a TypeScript replay', () => {
            const unclaimed = [...this.cases.keys()].filter((k) => !this.claimed.has(k));
            if (statsDir !== undefined) {
                mkdirSync(statsDir, { recursive: true });
                writeFileSync(join(statsDir, `${this.family}.json`),
                    JSON.stringify({ family: this.family, header: this.header,
                        cases: this.stats }, null, 1));
            }
            expect(unclaimed, 'golden cases without a replay').toEqual([]);
        });
    }
}
