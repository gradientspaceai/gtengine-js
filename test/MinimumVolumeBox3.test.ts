import { describe, it, expect } from 'vitest';
import {
    MVB3FloatingPoint,
    MVB3Rational,
    MVB3GPU,
    type MVB3ComputeType
} from '../src/MinimumVolumeBox3.js';
import { MinimumVolumeBox3FloatingPoint } from '../src/MinimumVolumeBox3FloatingPoint.js';
import { MinimumVolumeBox3Rational } from '../src/MinimumVolumeBox3Rational.js';
import { check, fc } from './helpers/arbitraries.js';

// MinimumVolumeBox3.h is a facade header: it defines only the compute-type
// tags (upstream std::integral_constant<std::size_t, N>) used to select a
// partial template specialization, plus an empty primary template. The port
// exposes the tags as marker constants; the implementations live in
// MinimumVolumeBox3FloatingPoint/Rational/GPU, ported separately.

describe('MinimumVolumeBox3 (facade)', () => {
    it('tags carry the upstream integral_constant values', () => {
        expect(MVB3FloatingPoint).toBe(0);
        expect(MVB3Rational).toBe(1);
        expect(MVB3GPU).toBe(2);
    });

    it('tags are pairwise distinct', () => {
        const tags = new Set<MVB3ComputeType>([MVB3FloatingPoint, MVB3Rational, MVB3GPU]);
        expect(tags.size).toBe(3);
    });

    it('MVB3ComputeType admits exactly the three tags', () => {
        // Compile-time check: each tag is assignable to the union type, and
        // a switch over the union is exhaustive.
        const describeTag = (tag: MVB3ComputeType): string => {
            switch (tag) {
                case MVB3FloatingPoint: return 'floating-point';
                case MVB3Rational: return 'rational';
                case MVB3GPU: return 'gpu';
            }
        };
        expect(describeTag(MVB3FloatingPoint)).toBe('floating-point');
        expect(describeTag(MVB3Rational)).toBe('rational');
        expect(describeTag(MVB3GPU)).toBe('gpu');
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). MinimumVolumeBox3.h contains
// no executable code: an empty primary template plus three
// std::integral_constant tags used for partial specialization. There is no
// numerical behaviour to cross-check, so the properties pin the tag values and
// the facade's documented mapping onto the ported implementations
// (MinimumVolumeBox3FloatingPoint / MinimumVolumeBox3Rational; the GPU
// specialization has no port).
// ---------------------------------------------------------------------------

describe('MinimumVolumeBox3 verification', () => {
    it('the tags are the upstream integral_constant values 0, 1, 2', () => {
        const tags: MVB3ComputeType[] = [MVB3FloatingPoint, MVB3Rational, MVB3GPU];
        check(fc.constantFrom(...tags), (tag) => {
            expect(Number.isInteger(tag)).toBe(true);
            expect(tags.indexOf(tag)).toBe(tag);
        });
    });

    it('a tag-keyed dispatch table selects the ported implementations', () => {
        // The facade exists only to name the specializations; this pins the
        // mapping recorded in the port notes.
        const implementations: Record<MVB3ComputeType, unknown> = {
            [MVB3FloatingPoint]: MinimumVolumeBox3FloatingPoint,
            [MVB3Rational]: MinimumVolumeBox3Rational,
            [MVB3GPU]: undefined
        };
        expect(typeof implementations[MVB3FloatingPoint]).toBe('function');
        expect(typeof implementations[MVB3Rational]).toBe('function');
        expect(implementations[MVB3GPU]).toBeUndefined();
    });
});
