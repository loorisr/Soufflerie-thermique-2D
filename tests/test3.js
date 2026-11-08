const { Solver } = require('../solver.js');
const floor = parseFloat(process.argv[2]);
const cfg = { fluid: 'air', L: 3.0, H: 0.6, nx: 320, Tamb: 15, gravity: false, depth: 0.4,
      sides: { gauche: 'ouvert', droite: 'ouvert', bas: 'paroi', haut: 'paroi' } };
const comps= [
      { id:1,type: 'fan', name: 'Ventilateur', x: 0.35, y: 0.0, w: 0.08, h: 0.6, dir: 'droite', Q: 1200 },
      { id:2,type: 'exchanger', name: 'Batterie chaude', x: 1.0, y: 0.0, w: 0.08, h: 0.6, mode: 'P', Tc: 50, ha: 30000, P: 4000, dp: 30, vnom: 1.5 },
      { id:3,type: 'solid', name: 'Obstacle', x: 1.8, y: 0.22, w: 0.16, h: 0.16, fixedT: false, T: 20 } ];
const s = new Solver(cfg, comps);
if (floor) { s.taug0 = Math.max(s.taug0, floor); }
for (let b=0;b<10;b++){ s.run(500); let lo=1e9,hi=-1e9,li=0,hiI=0; for(let i=0;i<s.N;i++) if(s.type[i]==0){ if(s.T[i]<lo){lo=s.T[i];li=i} if(s.T[i]>hi){hi=s.T[i];hiI=i}}
 console.log(s.step, 'T', lo.toFixed(1), hi.toFixed(1), 'hiat', hiI%s.nx, Math.floor(hiI/s.nx), s.zone[hiI], s.fan[hiI], 'P', s.exch[0].P.toFixed(0), 'Q', s.fanReport(s.fans[0]).Q.toFixed(0)); }
