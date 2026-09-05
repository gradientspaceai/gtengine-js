import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    include: ['test/**/*.test.ts'],
    passWithNoTests: true,
    // Property tests are CPU-bound and the suite runs alongside other
    // vitest processes (verifier agents) and on slow CI runners; the 5 s
    // default produced phantom timeouts in files that take ~1 s alone.
    testTimeout: 20000,
  },
});
