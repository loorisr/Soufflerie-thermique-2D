const { Solver } = require('../solver.js');
const cfg = { fluid: 'air', L: 2.0, H: 1.0, nx: 300, Tamb: 7, gravity: true, depth: 0.5,
  sides: { gauche: 'ouvert', droite: 'ouvert', bas: 'paroi', haut: 'ouvert' } };
const comps = [
  { id: 1, type: 'solid', x: 0.6, y: 0.2, w: 0.8, h: 0.03, fixedT: false, T: 20 },
  { id: 2, type: 'solid', x: 0.6, y: 0.77, w: 0.8, h: 0.03, fixedT: false, T: 20 },
  { id: 3, type: 'exchanger', x: 0.62, y: 0.23, w: 0.05, h: 0.54, mode: 'T', Tc: -3, ha: 30000, P: 3000, dp: 40, vnom: 2 },
  { id: 4, type: 'fan', x: 1.3, y: 0.23, w: 0.06, h: 0.54, dir: 'droite', Q: 1800 },
  { id: 5, type: 'solid', x: 0.0, y: 0.0, w: 2.0, h: 0.2, fixedT: false },
];
const s = new Solver(cfg, comps);
console.log('grid', s.nx, s.ny, 'dt', s.dt, 'tau0', s.tau0, s.tauPhys, 'taug0', s.taug0, 'Uref', s.Uref);
const t0 = Date.now();
for (let b = 0; b < 6; b++) {
  s.run(1000);
  let umax = 0, bad = 0;
  for (let i = 0; i < s.N; i++) { const u = Math.hypot(s.ux[i], s.uy[i]); if (!isFinite(u)) bad++; umax = Math.max(umax, u); }
  const fr = s.fanReport(s.fans[0]);
  console.log(s.step, 't=', s.simTime().toFixed(2), 's umax', (umax * s.cu).toFixed(2), 'Q', fr.Q.toFixed(0), 'dp', fr.dp.toFixed(1), 'P', s.exch[0].P.toFixed(0), 'Tm', s.exch[0].Tm.toFixed(2), 'bad', bad);
}
console.log('steps/s', (6000 / ((Date.now() - t0) / 1000)).toFixed(0));
let m=0,mi=0;for(let i=0;i<s.N;i++){const u=Math.hypot(s.ux[i],s.uy[i]);if(u>m){m=u;mi=i}}
console.log('max at', mi%s.nx, Math.floor(mi/s.nx), 'type', s.type[mi], 'fan', s.fan[mi], 'zone', s.zone[mi], 'rho', s.rho[mi]);
let cnt=0;for(let i=0;i<s.N;i++){if(Math.hypot(s.ux[i],s.uy[i])>0.25)cnt++}console.log('cells>0.25',cnt);
