// Regenerates src/index.ts from the files present in src/. Run: npm run gen:index
import { readdirSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const files = readdirSync(join(root, 'src'))
  .filter(f => f.endsWith('.ts') && f !== 'index.ts')
  .sort((a, b) => a.localeCompare(b, 'en'));

const header = `// gtengine-js: TypeScript port of the Geometric Tools Engine (GTE) Mathematics library.
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// GENERATED FILE - run 'npm run gen:index' to regenerate. Do not edit by hand.

`;
const body = files.length === 0
  ? 'export {};\n'
  : files.map(f => `export * from './${f.replace(/\.ts$/, '')}';`).join('\n') + '\n';
writeFileSync(join(root, 'src', 'index.ts'), header + body);
console.log(`src/index.ts: ${files.length} modules exported`);
