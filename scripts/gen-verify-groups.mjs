// Generates plan/verify-groups.json from plan/batches.json: every ported file,
// grouped by category (<= 5000 upstream lines and <= 24 files per group).
import { readFileSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const batches = JSON.parse(readFileSync(join(root, 'plan', 'batches.json'), 'utf8'));
const status = JSON.parse(readFileSync(join(root, 'porting-status.json'), 'utf8'));
const list = batches.batches ?? batches;
const byCat = {};
for (const b of list) (byCat[b.category] ??= []).push(b);
const groups = [];
for (const cat of Object.keys(byCat).sort()) {
  const bs = byCat[cat].sort((p, q) => p.depth - q.depth || p.id - q.id);
  let cur = null;
  for (const b of bs) {
    const files = b.files
      .filter(f => status[f.file]?.status === 'ported')
      .map(f => ({ file: f.file, lines: f.lines, batch: b.id }));
    if (files.length === 0) continue;
    const lines = files.reduce((n, f) => n + f.lines, 0);
    if (cur && cur.lines + lines <= 5000 && cur.files.length + files.length <= 24) {
      cur.files.push(...files);
      cur.lines += lines;
    } else {
      cur = { id: 0, category: cat, lines, files };
      groups.push(cur);
    }
  }
}
groups.forEach((g, i) => { g.id = i + 1; });
const plan = {
  description: 'Verification wave work queue (see VERIFYING.md): every ported file grouped by category, <=5000 upstream lines and <=24 files per group. One verify-batch issue and one PR per group.',
  groups,
};
writeFileSync(join(root, 'plan', 'verify-groups.json'), JSON.stringify(plan, null, 1) + '\n');
console.log(`groups ${groups.length}; files ${groups.reduce((n, g) => n + g.files.length, 0)}`);
