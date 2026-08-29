// Marks porting-status.json entries 'ported' for every header whose .ts file
// exists in src/. Run by the orchestrator after merging port PRs.
import { readFileSync, writeFileSync, existsSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const path = join(root, 'porting-status.json');
const status = JSON.parse(readFileSync(path, 'utf8'));
let changed = 0;
for (const [header, entry] of Object.entries(status)) {
  if (entry.status === 'pending' && existsSync(join(root, 'src', header.replace(/\.h$/, '.ts')))) {
    entry.status = 'ported';
    changed++;
  }
}
writeFileSync(path, JSON.stringify(status, null, 1) + '\n');
const counts = {};
for (const e of Object.values(status)) counts[e.status] = (counts[e.status] ?? 0) + 1;
console.log(`updated ${changed};`, JSON.stringify(counts));
