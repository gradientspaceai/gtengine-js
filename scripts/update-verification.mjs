// Marks every file of a verification group as verified in porting-status.json.
// Usage: node scripts/update-verification.mjs <groupId> [<groupId> ...]
import { readFileSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const statusPath = join(root, 'porting-status.json');
const status = JSON.parse(readFileSync(statusPath, 'utf8'));
const plan = JSON.parse(readFileSync(join(root, 'plan', 'verify-groups.json'), 'utf8'));
const ids = process.argv.slice(2).map(Number);
if (ids.length === 0 || ids.some(Number.isNaN)) {
  console.error('usage: update-verification.mjs <groupId> ...');
  process.exit(1);
}
let changed = 0;
for (const id of ids) {
  const g = plan.groups.find(g => g.id === id);
  if (!g) { console.error(`no group ${id}`); process.exit(1); }
  for (const f of g.files) {
    const e = status[f.file];
    if (!e) { console.error(`no manifest entry ${f.file}`); process.exit(1); }
    if (e.verified !== id) { e.verified = id; changed++; }
  }
}
writeFileSync(statusPath, JSON.stringify(status, null, 1) + '\n');
const ported = Object.values(status).filter(e => e.status === 'ported');
const done = ported.filter(e => e.verified).length;
console.log(`marked ${changed}; verified ${done}/${ported.length}`);
