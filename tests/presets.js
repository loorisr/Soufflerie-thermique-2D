// Fait tourner chaque cas exemple de index.src.html sur le solveur processeur et vérifie la stabilité.
const fs = require('fs'), path = require('path');
const root = path.join(__dirname, '..');
const { Solver } = require(path.join(root, 'solver.js'));
const src = fs.readFileSync(path.join(root, 'index.src.html'), 'utf8');
const a = src.indexOf('const PRESETS = '), b = src.indexOf('\n};\n', a);
const PRESETS = new Function('return ' + src.slice(a + 'const PRESETS = '.length, b + 2))();
const only = process.argv[2];
let fail = 0;
// eau froide : de l'eau à 2 °C ou à 6 °C dans de l'eau à 4 °C doit monter dans les deux cas
{
  const s = new Solver({ fluid: 'eau4', L: 0.1, H: 0.1, nx: 16, Tamb: 4, gravity: true, depth: 0.1,
    sides: { gauche: 'paroi', droite: 'paroi', bas: 'paroi', haut: 'paroi' } }, []);
  const ok = s.gb === 0 && s.gb2 > 0;
  if (!ok) fail++;
  console.log(`${ok ? 'ok   ' : 'ÉCHEC'} eau4      flottabilité symétrique autour de 4 °C`);
}
for (const [name, p] of Object.entries(PRESETS)) {
  if (only && name !== only) continue;
  let id = 1;
  const s = new Solver(p.cfg, p.comps.map(c => ({ id: id++, ...c })));
  const t0 = Date.now();
  let line = '', bad = 0;
  for (let k = 0; k < 4; k++) {
    s.run(1500);
    let um = 0, Tmin = Infinity, Tmax = -Infinity; bad = 0;
    for (let i = 0; i < s.N; i++) {
      const u = Math.hypot(s.ux[i], s.uy[i]);
      if (!isFinite(u) || !isFinite(s.T[i])) { bad++; continue; }
      um = Math.max(um, u);
      if (s.type[i] === 0) { Tmin = Math.min(Tmin, s.T[i]); Tmax = Math.max(Tmax, s.T[i]); }
    }
    line = `t=${s.simTime().toFixed(1)} s  umax=${(um * s.cu).toFixed(2)} m/s  T=[${Tmin.toFixed(1)}, ${Tmax.toFixed(1)}] °C`;
    if (s.fans.length) line += '  Q=' + s.fans.map(f => s.fanReport(f).Q.toFixed(0)).join('/') + ' m³/h';
    if (s.exch.length) line += '  P=' + s.exch.map(e => e.P.toFixed(0)).join('/') + ' W';
  }
  if (bad) fail++;
  console.log(`${bad ? 'ÉCHEC' : 'ok   '} ${name.padEnd(9)} ${s.nx}x${s.ny}  ${line}  (${((Date.now() - t0) / 1000).toFixed(1)} s)`);
}
process.exit(fail ? 1 : 0);
