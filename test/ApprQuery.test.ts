import { describe, it, expect, afterEach } from 'vitest';
import { ApprQuery } from '../src/ApprQuery';

// A miniature fitting query in the style of the Appr* files: fit the mean of
// a set of 1D observations. The model parameter is the mean; the model error
// for an observation is its squared deviation from the mean.
class ApprMean1 extends ApprQuery<number> {
    mean: number = 0;

    fitIndexed(observations: readonly number[],
        indices: readonly number[]): boolean {
        if (indices.length >= this.getMinimumRequired() &&
            this.validIndices(observations, indices)) {
            let sum = 0;
            for (const index of indices) {
                sum += observations[index];
            }
            this.mean = sum / indices.length;
            return true;
        }
        this.mean = 0;
        return false;
    }

    getMinimumRequired(): number {
        return 1;
    }

    error(observation: number): number {
        const diff = observation - this.mean;
        return diff * diff;
    }

    copyParameters(input: ApprQuery<number>): void {
        this.mean = (input as ApprMean1).mean;
    }
}

describe('ApprQuery', () => {
    afterEach(() => {
        ApprQuery.validateIndices = false;
    });

    describe('fit', () => {
        it('fits all observations', () => {
            const query = new ApprMean1();
            expect(query.fit([1, 2, 3, 4])).toBe(true);
            expect(query.mean).toBe(2.5);
        });

        it('fails with no observations', () => {
            const query = new ApprMean1();
            expect(query.fit([])).toBe(false);
        });

        it('fits a contiguous subset [imin, imax] inclusive', () => {
            const query = new ApprMean1();
            expect(query.fit([100, 2, 4, 6, 100], 1, 3)).toBe(true);
            expect(query.mean).toBe(4);

            // A single-element range.
            expect(query.fit([100, 2, 4, 6, 100], 4, 4)).toBe(true);
            expect(query.mean).toBe(100);
        });

        it('fails for a reversed contiguous range (imin > imax)', () => {
            const query = new ApprMean1();
            expect(query.fit([1, 2, 3], 2, 1)).toBe(false);
        });

        it('fits an indexed subset', () => {
            const query = new ApprMean1();
            expect(query.fit([10, 999, 20, 999, 30], [0, 2, 4])).toBe(true);
            expect(query.mean).toBe(20);
        });

        it('allows duplicate indices (weighted fit)', () => {
            const query = new ApprMean1();
            expect(query.fit([0, 6], [1, 1, 0])).toBe(true);
            expect(query.mean).toBe(4);
        });

        it('fits a prefix of the indices when numIndices is passed', () => {
            const query = new ApprMean1();
            expect(query.fit([10, 20, 30, 999], [0, 1, 2, 3], 3)).toBe(true);
            expect(query.mean).toBe(20);
        });

        it('clamps numIndices to indices.length', () => {
            const query = new ApprMean1();
            expect(query.fit([10, 20], [0, 1], 100)).toBe(true);
            expect(query.mean).toBe(15);
        });
    });

    describe('error', () => {
        it('is the squared deviation from the fitted model', () => {
            const query = new ApprMean1();
            query.fit([1, 3]);
            expect(query.mean).toBe(2);
            expect(query.error(2)).toBe(0);
            expect(query.error(5)).toBe(9);
        });
    });

    describe('validIndices', () => {
        it('accepts anything when validation is disabled (default)', () => {
            const query = new ApprMean1();
            expect(ApprQuery.validateIndices).toBe(false);
            expect(query.validIndices([1, 2], [5, -1, 7])).toBe(true);
        });

        it('checks bounds and counts when validation is enabled', () => {
            ApprQuery.validateIndices = true;
            const query = new ApprMean1();

            // Valid subsets.
            expect(query.validIndices([1, 2, 3], [0, 2])).toBe(true);
            expect(query.validIndices([1, 2, 3], [0, 1, 2])).toBe(true);

            // Out-of-range indices.
            expect(query.validIndices([1, 2, 3], [0, 3])).toBe(false);
            expect(query.validIndices([1, 2, 3], [-1, 0])).toBe(false);

            // Fewer indices than the minimum required.
            expect(query.validIndices([1, 2, 3], [])).toBe(false);

            // More indices than observations.
            expect(query.validIndices([1, 2], [0, 1, 0])).toBe(false);
        });

        it('makes fit reject invalid indices when validation is enabled', () => {
            ApprQuery.validateIndices = true;
            const query = new ApprMean1();
            expect(query.fit([1, 2, 3], [0, 7])).toBe(false);
            expect(query.fit([1, 2, 3], [0, 2])).toBe(true);
            expect(query.mean).toBe(2);
        });
    });

    describe('ransac', () => {
        // Twenty inliers near 10 and five gross outliers near 100.
        const inliers = Array.from({ length: 20 }, (_, i) => 10 + 0.01 * (i - 10));
        const outliers = [100, 101, 99, 102, 98];
        const observations = [...inliers, ...outliers];

        it('rejects the outliers and fits the inlier consensus set', () => {
            const candidateModel = new ApprMean1();
            const bestModel = new ApprMean1();
            const { success, bestConsensus } = ApprQuery.ransac(
                candidateModel, observations, 15, 1.0, 32, bestModel);

            expect(success).toBe(true);

            // The consensus set is exactly the inliers (every inlier is
            // within maxErrorForGoodFit of the inlier mean, and every
            // outlier is far away).
            expect(bestConsensus.length).toBe(inliers.length);
            expect([...bestConsensus].sort((a, b) => a - b)).toEqual(
                Array.from({ length: inliers.length }, (_, i) => i));

            // The best model is the least-squares fit of the consensus set.
            const inlierMean = inliers.reduce((a, b) => a + b, 0) / inliers.length;
            expect(bestModel.mean).toBeCloseTo(inlierMean, 12);
            expect(Math.abs(bestModel.mean - 10)).toBeLessThan(0.2);
        });

        it('is deterministic for fixed inputs', () => {
            const run = (): { mean: number, consensus: number[] } => {
                const candidateModel = new ApprMean1();
                const bestModel = new ApprMean1();
                const result = ApprQuery.ransac(
                    candidateModel, observations, 15, 1.0, 8, bestModel);
                expect(result.success).toBe(true);
                return { mean: bestModel.mean, consensus: result.bestConsensus };
            };
            const first = run();
            const second = run();
            expect(second.mean).toBe(first.mean);
            expect(second.consensus).toEqual(first.consensus);
        });

        it('fails when there are too few observations', () => {
            const { success, bestConsensus } = ApprQuery.ransac(
                new ApprMean1(), [], 1, 1.0, 10, new ApprMean1());
            expect(success).toBe(false);
            expect(bestConsensus).toEqual([]);
        });

        it('fits the entire set when it has exactly the minimum size', () => {
            const bestModel = new ApprMean1();
            const { success, bestConsensus } = ApprQuery.ransac(
                new ApprMean1(), [7], 1, 1.0, 10, bestModel);
            expect(success).toBe(true);
            expect(bestConsensus).toEqual([0]);
            expect(bestModel.mean).toBe(7);
        });

        it('fails when no consensus set is large enough', () => {
            // Spread-out data with a tight error tolerance: no subset of 6
            // observations agrees to within the tolerance.
            const spread = [0, 1000, 2000, 3000, 4000, 5000, 6000];
            const { success } = ApprQuery.ransac(
                new ApprMean1(), spread, 6, 1.0, 20, new ApprMean1());
            expect(success).toBe(false);
        });

        it('copies parameters into bestModel rather than aliasing', () => {
            const candidateModel = new ApprMean1();
            const bestModel = new ApprMean1();
            const { success } = ApprQuery.ransac(
                candidateModel, observations, 15, 1.0, 32, bestModel);
            expect(success).toBe(true);
            expect(bestModel).not.toBe(candidateModel);

            // Mutating the candidate model afterwards must not affect the
            // best model.
            const bestMean = bestModel.mean;
            candidateModel.fit([12345]);
            expect(bestModel.mean).toBe(bestMean);
        });
    });
});
