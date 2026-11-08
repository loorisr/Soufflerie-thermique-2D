const { Solver } = require('../solver.js');
for (const fluid of ['air','eau']) {
const cfg = { fluid, L: 3.0, H: 2.5, nx: 200, Tamb: 18, gravity: true, depth: 1,
  sides: { gauche: 'paroi', droite: 'paroi', bas: 'paroi', haut: 'paroi' } };
const comps = [
  { id: 1, type: 'solid', x: 0.1, y: 0.2, w: 0.06, h: 0.6, fixedT: true, T: 60 },
  { id: 2, type: 'solid', x: 2.9, y: 0.5, w: 0.1, h: 1.5, fixedT: true, T: 8 },
];
const s = new Solver(cfg, comps);
console.log(fluid,'grid', s.nx, s.ny, 'dt', s.dt.toExponential(2), 'Uref', s.Uref.toFixed(3));
for (let b = 0; b < 4; b++) { s.run(2000);
  let umax=0,bad=0,Tm=0;for(let i=0;i<s.N;i++){const u=Math.hypot(s.ux[i],s.uy[i]);if(!isFinite(u)||!isFinite(s.T[i]))bad++;umax=Math.max(umax,u);if(s.type[i]==0)Tm+=s.T[i]}
  console.log(s.step,'t',s.simTime().toFixed(1),'umax',(umax*s.cu).toFixed(3),'Tmean',(Tm/s.N).toFixed(2),'bad',bad);}
}
