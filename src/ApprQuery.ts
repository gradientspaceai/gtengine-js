// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprQuery.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Base class support for least-squares fitting algorithms and for RANSAC
// algorithms.
//
// Upstream declares 'template <typename Real, typename ObservationType>
// class ApprQuery'. The Real template parameter is dropped (always 'number'
// in the port); ObservationType remains a type parameter. Every Appr*
// fitting class derives from ApprQuery and must implement this contract:
//
//   fitIndexed(observations, indices): boolean
//     Estimate the model parameters from the subset of observations selected
//     by 'indices' (indices into 'observations'; duplicates permitted by the
//     base class). Store the parameters in the derived-class object and
//     return whether the fit succeeded. All base-class fit(...) overloads
//     funnel into this function. Implementations should first call
//     this.validIndices(observations, indices) and fail when it returns
//     false.
//
//   getMinimumRequired(): number
//     The minimum number of observations required to fit the model (e.g. 2
//     for a line in 2D, 3 for a circle).
//
//   error(observation): number
//     The model error for the specified observation under the current model
//     parameters (e.g. the squared distance from the observation to the
//     fitted model). Used by RANSAC to classify inliers.
//
//   copyParameters(input): void
//     Copy the model parameters from 'input' (an object of the same derived
//     class) into 'this'. Used by RANSAC to copy the candidate-model
//     parameters to the current best-fit model. C++ value assignment copies;
//     TypeScript implementations must deep-copy any parameters stored in
//     objects or arrays.
//
// Port notes: The C++ fit functions take raw pointers with separate element
// counts (numObservations, numIndices); TypeScript arrays carry their
// length, so the pointer+count pairs collapse to arrays and the
// Fit(size_t, ObservationType const*) overload merges with
// Fit(std::vector const&) into fit(observations). The remaining C++ Fit
// overloads become TypeScript overloads of fit(...); the indexed overload
// Fit(observations, indices) is virtual upstream but is never overridden by
// a derived class -- the customization point in the port is fitIndexed.
// Upstream's compile-time index validation toggle
// GTE_APPR_QUERY_VALIDATE_INDICES becomes the mutable static
// ApprQuery.validateIndices (default false, matching the upstream default of
// the macro being undefined). Upstream RANSAC shuffles with a
// freshly-default-constructed std::default_random_engine each iteration; the
// port mirrors this with a fresh minstd-style (Lehmer) generator per
// iteration, so RANSAC is deterministic for fixed inputs in both versions,
// though the exact permutations are implementation-defined in C++ and
// therefore differ. The RANSAC output reference parameter bestConsensus
// becomes a field of the returned object literal.

// The multiplicative congruential generator x -> 16807 * x mod (2^31 - 1)
// with seed 1 (std::minstd_rand0, the common choice for
// std::default_random_engine). The products stay below 2^45, so the
// arithmetic is exact in IEEE doubles.
class DefaultRandomEngine {
    private mState = 1;

    // Returns a pseudo-random number in [1, 2^31 - 2].
    next(): number {
        this.mState = (16807 * this.mState) % 2147483647;
        return this.mState;
    }
}

// Randomly permute the array in place (the port of std::shuffle over the
// whole array): Fisher-Yates using the engine to draw an index in [0, i].
function shuffle(values: number[], engine: DefaultRandomEngine): void {
    for (let i = values.length - 1; i > 0; --i) {
        const j = engine.next() % (i + 1);
        const temp = values[i];
        values[i] = values[j];
        values[j] = temp;
    }
}

export abstract class ApprQuery<ObservationType> {
    // The port of the compile-time macro GTE_APPR_QUERY_VALIDATE_INDICES.
    // Set to true to make validIndices verify that the incoming indices to
    // the fitting functions are valid; when false (the default),
    // validIndices returns true and the caller is responsible for passing
    // correctly formed data.
    static validateIndices: boolean = false;

    // The base-class fit overloads are generic but need to call the indexed
    // fitting function for the specific derived class.
    abstract fitIndexed(observations: readonly ObservationType[],
        indices: readonly number[]): boolean;

    validIndices(observations: readonly ObservationType[],
        indices: readonly number[]): boolean {
        if (!ApprQuery.validateIndices) {
            // The caller is responsible for passing correctly formed data.
            return true;
        }

        if (this.getMinimumRequired() <= indices.length &&
            indices.length <= observations.length) {
            for (const index of indices) {
                if (index < 0 || index >= observations.length) {
                    return false;
                }
            }
            return true;
        }
        return false;
    }

    // Estimate the model parameters for all observations.
    fit(observations: readonly ObservationType[]): boolean;

    // Estimate the model parameters for the contiguous subset of
    // observations with indices imin through imax (inclusive).
    fit(observations: readonly ObservationType[],
        imin: number, imax: number): boolean;

    // Estimate the model parameters for an indexed subset of observations.
    fit(observations: readonly ObservationType[],
        indices: readonly number[]): boolean;

    // Estimate the model parameters for the subset of observations
    // specified by the first numIndices elements of indices, where
    // numIndices is possibly smaller than indices.length.
    fit(observations: readonly ObservationType[],
        indices: readonly number[], numIndices: number): boolean;

    fit(observations: readonly ObservationType[],
        arg1?: number | readonly number[], arg2?: number): boolean {
        if (arg1 === undefined) {
            // All observations.
            const indices = new Array<number>(observations.length);
            for (let i = 0; i < indices.length; ++i) {
                indices[i] = i;
            }
            return this.fitIndexed(observations, indices);
        }

        if (typeof arg1 === 'number') {
            // Contiguous subset [imin, imax].
            const imin = arg1;
            const imax = arg2 as number;
            if (imin <= imax) {
                const numIndices = imax - imin + 1;
                const indices = new Array<number>(numIndices);
                for (let i = 0; i < numIndices; ++i) {
                    indices[i] = imin + i;
                }
                return this.fitIndexed(observations, indices);
            }
            return false;
        }

        if (arg2 === undefined) {
            // Indexed subset.
            return this.fitIndexed(observations, arg1);
        }

        // Prefix of the indexed subset.
        const imax = Math.min(arg2, arg1.length);
        return this.fitIndexed(observations, arg1.slice(0, imax));
    }

    // Apply the RANdom SAmple Consensus algorithm for fitting a model to
    // observations. The algorithm requires three functions to be
    // implemented by the derived classes (in addition to fitIndexed).

    // The minimum number of observations required to fit the model.
    abstract getMinimumRequired(): number;

    // Compute the model error for the specified observation for the current
    // model parameters.
    abstract error(observation: ObservationType): number;

    // Copy the parameters between two models. This is used to copy the
    // candidate-model parameters to the current best-fit model.
    abstract copyParameters(input: ApprQuery<ObservationType>): void;

    // On success (a consensus set of at least numRequiredForGoodFit
    // observations was found), bestModel holds the best-fit parameters and
    // bestConsensus holds the indices of the observations in its consensus
    // set. candidateModel and bestModel must be distinct objects of the same
    // derived class.
    static ransac<ObservationType>(
        candidateModel: ApprQuery<ObservationType>,
        observations: readonly ObservationType[],
        numRequiredForGoodFit: number, maxErrorForGoodFit: number,
        numIterations: number, bestModel: ApprQuery<ObservationType>):
        { success: boolean, bestConsensus: number[] } {
        const numObservations = observations.length;
        const minRequired = candidateModel.getMinimumRequired();
        if (numObservations < minRequired) {
            // Too few observations for model fitting.
            return { success: false, bestConsensus: [] };
        }

        // The first part of the array will store the consensus set,
        // initially filled with the minimum number of indices that
        // correspond to the candidate inliers. The last part will store the
        // remaining indices. These points are tested against the model and
        // are added to the consensus set when they fit. All the index
        // manipulation is done in place. Initially, the candidates are the
        // identity permutation.
        const candidates = new Array<number>(numObservations);
        for (let i = 0; i < numObservations; ++i) {
            candidates[i] = i;
        }

        if (numObservations === minRequired) {
            // We have the minimum number of observations to generate the
            // model, so RANSAC cannot be used. Compute the model with the
            // entire set of observations.
            return { success: bestModel.fit(observations), bestConsensus: candidates };
        }

        let bestNumFittedObservations = minRequired;
        let bestConsensus: number[] = [];

        for (let i = 0; i < numIterations; ++i) {
            // Randomly permute the previous candidates, partitioning the
            // array into getMinimumRequired() indices (the candidate
            // inliers) followed by the remaining indices (candidates for
            // testing against the model).
            shuffle(candidates, new DefaultRandomEngine());

            // Fit the model to the inliers.
            if (candidateModel.fit(observations, candidates, minRequired)) {
                // Test each remaining observation whether it fits the
                // model. If it does, include it in the consensus set.
                let numFittedObservations = minRequired;
                for (let j = minRequired; j < numObservations; ++j) {
                    const error = candidateModel.error(observations[candidates[j]]);
                    if (error <= maxErrorForGoodFit) {
                        const temp = candidates[j];
                        candidates[j] = candidates[numFittedObservations];
                        candidates[numFittedObservations] = temp;
                        ++numFittedObservations;
                    }
                }

                if (numFittedObservations >= numRequiredForGoodFit) {
                    // We have observations that fit the model. Update the
                    // best model using the consensus set.
                    candidateModel.fit(observations, candidates, numFittedObservations);
                    if (numFittedObservations > bestNumFittedObservations) {
                        // The consensus set is larger than the previous
                        // consensus set, so its model becomes the best one.
                        bestModel.copyParameters(candidateModel);
                        bestConsensus = candidates.slice(0, numFittedObservations);
                        bestNumFittedObservations = numFittedObservations;
                    }
                }
            }
        }

        return {
            success: bestNumFittedObservations >= numRequiredForGoodFit,
            bestConsensus: bestConsensus
        };
    }
}
