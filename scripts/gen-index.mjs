// Regenerates src/index.ts from the files present in src/. Run: npm run gen:index
import { readdirSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const files = readdirSync(join(root, 'src'))
  .filter(f => f.endsWith('.ts') && f !== 'index.ts')
  .sort((a, b) => a.localeCompare(b, 'en'));

// Guard: with star-exports, duplicate exported names are silently dropped by
// TS/ESM instead of erroring. Detect them here and fail loudly.
import { readFileSync } from 'node:fs';
const owners = new Map();
let dupes = 0;
for (const f of files) {
  const text = readFileSync(join(root, 'src', f), 'utf8');
  for (const m of text.matchAll(/^export\s+(?:declare\s+)?(?:abstract\s+)?(?:function|class|interface|type|const|let|enum|namespace)\s+([A-Za-z_$][\w$]*)/gm)) {
    const name = m[1];
    if (owners.get(name) === f) continue; // overload signatures in the same file are fine
    if (owners.has(name)) {
      console.error(`DUPLICATE EXPORT: '${name}' in ${f} already exported by ${owners.get(name)}`);
      dupes++;
    } else {
      owners.set(name, f);
    }
  }
}
if (dupes > 0) {
  console.error(`${dupes} duplicate export name(s) — rename per PORTING.md before regenerating index.`);
  process.exit(1);
}

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
