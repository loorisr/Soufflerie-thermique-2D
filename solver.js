// ===== Solveur LBM thermique 2D — D2Q9 MRT (écoulement) + D2Q5 BGK (température, diffusivité adaptée au Péclet de maille) =====
// Boussinesq, LES Smagorinsky, zones poreuses avec source de chaleur (échangeurs),
// zones de force asservies au débit (ventilateurs).
// Stockage : f~ = f - w (écart au repos) et θ = T - Tamb, en simple précision (identique au GPU).

const FLUIDS = {
  air: { name: 'Air', rho: 1.204, nu: 1.50e-5, Pr: 0.71, beta: 3.43e-3, cp: 1006 },
  eau: { name: 'Eau', rho: 998.2, nu: 1.00e-6, Pr: 7.0, beta: 2.07e-4, cp: 4182 },
};

const CX = [0, 1, 0, -1, 0, 1, -1, -1, 1];
const CY = [0, 0, 1, 0, -1, 1, 1, -1, -1];
const W9 = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36];
const OPP = [0, 3, 4, 1, 2, 7, 8, 5, 6];
const W5 = [1 / 3, 1 / 6, 1 / 6, 1 / 6, 1 / 6];
// Matrice MRT de Lallemand & Luo : (ρ, e, ε, jx, qx, jy, qy, pxx, pxy)
const MM = [
  [1, 1, 1, 1, 1, 1, 1, 1, 1],
  [-4, -1, -1, -1, -1, 2, 2, 2, 2],
  [4, -2, -2, -2, -2, 1, 1, 1, 1],
  [0, 1, 0, -1, 0, 1, -1, -1, 1],
  [0, -2, 0, 2, 0, 1, -1, -1, 1],
  [0, 0, 1, 0, -1, 1, 1, -1, -1],
  [0, 0, -2, 0, 2, 1, 1, -1, -1],
  [0, 1, -1, 1, -1, 0, 0, 0, 0],
  [0, 0, 0, 0, 0, 1, -1, 1, -1],
];
const MNORM = [9, 36, 36, 6, 12, 6, 12, 4, 4];
const S_E = 1.19, S_EPS = 1.4;   // taux de relaxation des moments non hydrodynamiques
const CS = 0.16, PRT = 0.85;     // Smagorinsky, Prandtl turbulent
const TAU_MIN = 0.503, TAUG_MIN = 0.503;
const PE_CELL = 6;               // Péclet de maille max. : la diffusivité thermique n'est relevée que là où l'écoulement est rapide
const MAXF = 16, MAXE = 16;

const DIRS = { droite: [1, 0], gauche: [-1, 0], haut: [0, 1], bas: [0, -1] };

class Solver {
  constructor(cfg, comps) {
    this.cfg = cfg;
    this.fl = FLUIDS[cfg.fluid];
    const nx = cfg.nx, ny = Math.max(8, Math.round(cfg.nx * cfg.H / cfg.L));
    this.nx = nx; this.ny = ny;
    const N = nx * ny; this.N = N;
    this.dx = cfg.L / nx;
    this.f0 = new Float32Array(9 * N); this.f1 = new Float32Array(9 * N);
    this.g0 = new Float32Array(5 * N); this.g1 = new Float32Array(5 * N);
    this.rho = new Float32Array(N).fill(1);
    this.ux = new Float32Array(N); this.uy = new Float32Array(N);
    this.T = new Float32Array(N).fill(cfg.Tamb);    // température réelle (affichage)
    this.th = new Float32Array(N);                   // θ = T - Tamb (calcul)
    this.type = new Uint8Array(N);                   // 0 fluide, 1 solide adiabatique, 2 solide à T imposée
    this.Tw = new Float32Array(N);                   // θ de paroi
    this.nb = new Int32Array(9 * N);
    this.sponge = new Float32Array(N);
    this.drag = new Float32Array(N);
    this.hk = new Float32Array(N);
    this.hT = new Float32Array(N);                   // θ batterie
    this.hq = new Float32Array(N);
    this.zone = new Int32Array(N).fill(-1);
    this.fan = new Int32Array(N).fill(-1);
    this.reset = new Uint32Array(N);                 // cellules à réinitialiser (GPU)
    this.step = 0;
    this.fans = []; this.exch = [];
    this.setScales(comps);
    this.setComponents(comps);
  }

  setScales(comps) {
    const cfg = this.cfg, fl = this.fl, dx = this.dx;
    let U = 0.05, dTmax = 2;
    for (const c of comps) {
      if (c.type === 'fan') {
        const d = DIRS[c.dir];
        const wperp = d[0] !== 0 ? c.h : c.w;
        U = Math.max(U, (c.Q / 3600) / cfg.depth / Math.max(wperp, dx));
      }
      if (c.type === 'exchanger') dTmax = Math.max(dTmax, c.mode === 'T' ? Math.abs(c.Tc - cfg.Tamb) : 15);
      if (c.type === 'solid' && c.fixedT) dTmax = Math.max(dTmax, Math.abs(c.T - cfg.Tamb));
    }
    if (cfg.gravity) U = Math.max(U, 0.5 * Math.sqrt(9.81 * fl.beta * dTmax * cfg.H));
    this.Uref = U;
    this.ulat = 0.07;
    this.dt = this.ulat * dx / U;
    this.cu = dx / this.dt;
    this.ps = fl.rho * this.cu * this.cu;
    const nul = fl.nu * this.dt / (dx * dx);
    const all = fl.nu / fl.Pr * this.dt / (dx * dx);
    this.tauPhys = 0.5 + 3 * nul;
    this.tau0 = Math.max(this.tauPhys, TAU_MIN);
    this.taug0 = Math.max(0.5 + 3 * all, TAUG_MIN);
    this.nuEff = (this.tau0 - 0.5) / 3 * dx * dx / this.dt;
    this.alphaEff = (this.taug0 - 0.5) / 3 * dx * dx / this.dt;
    this.gb = cfg.gravity ? 9.81 * fl.beta * this.dt * this.dt / dx : 0;
  }

  setComponents(comps) {
    const { nx, ny, N, dx, cfg, fl } = this;
    const prev = this.type.slice();
    const prevFans = this.fans;
    this.type.fill(0); this.drag.fill(0); this.hk.fill(0); this.hq.fill(0); this.hT.fill(0);
    this.zone.fill(-1); this.fan.fill(-1); this.reset.fill(0);
    this.exch = []; this.fans = [];
    const cells = (c) => {
      const out = [];
      const i0 = Math.max(0, Math.floor(c.x / dx)), i1 = Math.min(nx - 1, Math.ceil((c.x + c.w) / dx));
      const j0 = Math.max(0, Math.floor(c.y / dx)), j1 = Math.min(ny - 1, Math.ceil((c.y + c.h) / dx));
      for (let j = j0; j <= j1; j++) for (let i = i0; i <= i1; i++) {
        const x = (i + 0.5) * dx, y = (j + 0.5) * dx;
        if (x >= c.x && x < c.x + c.w && y >= c.y && y < c.y + c.h) out.push(j * nx + i);
      }
      return out;
    };
    for (const c of comps) {
      if (c.type !== 'fan' || this.fans.length >= MAXF) continue;
      const id = this.fans.length;
      const d = DIRS[c.dir];
      const cl = cells(c);
      const thick = Math.max(1, Math.round((d[0] !== 0 ? c.w : c.h) / dx));
      const wperp = d[0] !== 0 ? c.h : c.w;
      const old = prevFans.find(f => f.comp.id === c.id);
      this.fans.push({ comp: c, dx: d[0], dy: d[1], cells: cl, n: cl.length, thick, wperp,
        F: old ? old.F : 0, meas: old ? old.meas : 0 });
      for (const i of cl) this.fan[i] = id;
    }
    for (const c of comps) {
      if (c.type !== 'exchanger' || this.exch.length >= MAXE) continue;
      const id = this.exch.length;
      const cl = cells(c);
      const e = Math.max(dx, Math.min(c.w, c.h));
      const C2 = c.dp / (e * fl.rho * Math.max(0.01, c.vnom) ** 2);
      const area = cl.length * dx * dx;
      let k = 0, q = 0;
      if (c.mode === 'T') k = 1 - Math.exp(-c.ha / (fl.rho * fl.cp) * this.dt);
      else q = (c.P / cfg.depth) / Math.max(area, 1e-9) / (fl.rho * fl.cp) * this.dt;
      this.exch.push({ comp: c, cells: cl, n: cl.length, P: 0, Tm: cfg.Tamb });
      for (const i of cl) {
        this.zone[i] = id; this.drag[i] = C2 * dx;
        this.hk[i] = k; this.hT[i] = c.Tc - cfg.Tamb; this.hq[i] = q;
        this.fan[i] = -1;
      }
    }
    for (const c of comps) {
      if (c.type !== 'solid') continue;
      for (const i of cells(c)) {
        this.type[i] = c.fixedT ? 2 : 1; this.Tw[i] = c.fixedT ? c.T - cfg.Tamb : 0;
        this.zone[i] = -1; this.fan[i] = -1; this.drag[i] = 0; this.hk[i] = 0; this.hq[i] = 0;
      }
    }
    for (let i = 0; i < N; i++) {
      if (this.type[i] === 0 && prev[i] !== 0) { this.initCell(i); this.reset[i] = 1; }
      if (this.type[i] !== 0) {
        this.ux[i] = 0; this.uy[i] = 0; this.rho[i] = 1;
        this.T[i] = this.type[i] === 2 ? this.Tw[i] + cfg.Tamb : cfg.Tamb;
      }
    }
    this.buildNeighbors();
  }

  initCell(i) {
    const N = this.N;
    for (let k = 0; k < 9; k++) { this.f0[k * N + i] = 0; this.f1[k * N + i] = 0; }
    for (let k = 0; k < 5; k++) { this.g0[k * N + i] = 0; this.g1[k * N + i] = 0; }
    this.rho[i] = 1; this.ux[i] = 0; this.uy[i] = 0; this.th[i] = 0; this.T[i] = this.cfg.Tamb;
  }

  // voisins : >=0 source fluide ; -1 rebond adiabatique ; -2 bord ouvert ; <=-3 paroi à T imposée (-3-idx)
  buildNeighbors() {
    const { nx, ny, N, type, nb } = this;
    const s = this.cfg.sides;
    const Ls = Math.max(6, Math.round(Math.min(nx, ny) * 0.08));
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const idx = j * nx + i;
      for (let k = 0; k < 9; k++) {
        const si = i - CX[k], sj = j - CY[k];
        let v;
        if (si < 0 || si >= nx || sj < 0 || sj >= ny) {
          let open = true;
          if (si < 0) open = open && s.gauche === 'ouvert';
          if (si >= nx) open = open && s.droite === 'ouvert';
          if (sj < 0) open = open && s.bas === 'ouvert';
          if (sj >= ny) open = open && s.haut === 'ouvert';
          v = open ? -2 : -1;
        } else {
          const sidx = sj * nx + si;
          v = type[sidx] === 0 ? sidx : (type[sidx] === 2 ? -3 - sidx : -1);
        }
        nb[k * N + idx] = v;
      }
      let d = 1e9;
      if (s.gauche === 'ouvert') d = Math.min(d, i);
      if (s.droite === 'ouvert') d = Math.min(d, nx - 1 - i);
      if (s.bas === 'ouvert') d = Math.min(d, j);
      if (s.haut === 'ouvert') d = Math.min(d, ny - 1 - j);
      const r = d < Ls ? 1 - d / Ls : 0;
      this.sponge[idx] = 1.2 * r * r;
    }
  }

  ramp() { return Math.min(1, this.step / 3000); }
  fanTargetLat(fn) {
    const c = fn.comp;
    return (c.Q / 3600) / this.cfg.depth / Math.max(fn.wperp, this.dx) / this.cu;
  }
  fanKi(fn) { return 2e-3 / fn.thick; }

  run(n) { for (let s = 0; s < n; s++) this.stepOnce(); this.sync(); }

  stepOnce() {
    const { N, nb, type, drag, hk, hT, hq, zone, fan, Tw, sponge } = this;
    const fA = this.f0, fB = this.f1, gA = this.g0, gB = this.g1;
    const rhoA = this.rho, uxA = this.ux, uyA = this.uy, thA = this.th;
    const tau0 = this.tau0, taug0 = this.taug0, gb = this.gb;
    const Cs2x = 18 * Math.SQRT2 * CS * CS;
    const fans = this.fans, exch = this.exch;
    const ramp = this.ramp();
    const nf = fans.length, ne = exch.length;
    const fanF = new Float64Array(nf), fdx = new Float64Array(nf), fdy = new Float64Array(nf);
    for (let q = 0; q < nf; q++) { fanF[q] = fans[q].F; fdx[q] = fans[q].dx; fdy[q] = fans[q].dy; }
    const fanAcc = new Float64Array(nf), exAcc = new Float64Array(ne), exT = new Float64Array(ne);
    const f = new Float64Array(9), g = new Float64Array(5), m = new Float64Array(9);

    for (let idx = 0; idx < N; idx++) {
      if (type[idx] !== 0) continue;
      const u0x = uxA[idx], u0y = uyA[idx], u0sq = u0x * u0x + u0y * u0y;
      for (let k = 0; k < 9; k++) {
        const v = nb[k * N + idx];
        if (v >= 0) f[k] = fA[k * N + v];
        else if (v === -2) {
          const cu = CX[k] * u0x + CY[k] * u0y;
          f[k] = -fA[OPP[k] * N + idx] + 2 * W9[k] * (4.5 * cu * cu - 1.5 * u0sq);
        } else f[k] = fA[OPP[k] * N + idx];
      }
      const thc = thA[idx];
      for (let k = 0; k < 5; k++) {
        const v = nb[k * N + idx];
        if (v >= 0) g[k] = gA[k * N + v];
        else if (v === -2) {
          const cu = CX[k] * u0x + CY[k] * u0y;     // >0 : l'écoulement entre par ce bord
          const tin = cu > 0 ? 0 : thc;
          g[k] = W5[k] * tin * (1 + 3 * cu);
        } else if (v === -1) g[k] = gA[OPP[k] * N + idx];
        else g[k] = -gA[OPP[k] * N + idx] + 2 * W5[k] * Tw[-3 - v];
      }
      // moments
      {
        const f0 = f[0], f1 = f[1], f2 = f[2], f3 = f[3], f4 = f[4], f5 = f[5], f6 = f[6], f7 = f[7], f8 = f[8];
        const sa = f1 + f2 + f3 + f4, sd = f5 + f6 + f7 + f8;
        m[0] = f0 + sa + sd; m[1] = -4 * f0 - sa + 2 * sd; m[2] = 4 * f0 - 2 * sa + sd;
        m[3] = f1 - f3 + f5 - f6 - f7 + f8; m[4] = -2 * f1 + 2 * f3 + f5 - f6 - f7 + f8;
        m[5] = f2 - f4 + f5 + f6 - f7 - f8; m[6] = -2 * f2 + 2 * f4 + f5 + f6 - f7 - f8;
        m[7] = f1 - f2 + f3 - f4; m[8] = f5 - f6 + f7 - f8;
      }
      let drho = m[0];
      let rho = 1 + drho;
      if (!(rho > 0.2 && rho < 5)) { rho = 1; drho = 0; }
      const mx = m[3], my = m[5];
      const th = g[0] + g[1] + g[2] + g[3] + g[4];
      // forces
      let Fx = 0, Fy = rho * gb * th;
      const fq = fan[idx];
      if (fq >= 0) { Fx += rho * fanF[fq] * fdx[fq]; Fy += rho * fanF[fq] * fdy[fq]; }
      let ux = (mx + 0.5 * Fx) / rho, uy = (my + 0.5 * Fy) / rho;
      const C = drag[idx];
      if (C > 0) {
        const s = 0.5 * C * Math.sqrt(ux * ux + uy * uy) + 0.002;
        ux /= (1 + s); uy /= (1 + s);
        Fx -= 2 * rho * s * ux; Fy -= 2 * rho * s * uy;
      }
      const um = Math.sqrt(ux * ux + uy * uy);
      if (um > 0.3) { ux *= 0.3 / um; uy *= 0.3 / um; }
      if (fq >= 0) fanAcc[fq] += ux * fdx[fq] + uy * fdy[fq];
      const usq = ux * ux + uy * uy;
      // équilibres (écarts au repos)
      const eEq = -2 * drho + 3 * rho * usq;
      const epsEq = drho - 3 * rho * usq;
      const pxxEq = rho * (ux * ux - uy * uy), pxyEq = rho * ux * uy;
      const qxEq = -rho * ux, qyEq = -rho * uy;
      // LES : norme du tenseur hors équilibre (moments du second ordre)
      const nE = m[1] - eEq, nXX = m[7] - pxxEq, nXY = m[8] - pxyEq;
      const Pxx = (nE / 3 + nXX) * 0.5, Pyy = (nE / 3 - nXX) * 0.5, Pxy = nXY;
      const Q = Math.sqrt(Pxx * Pxx + Pyy * Pyy + 2 * Pxy * Pxy);
      const t0 = tau0 + sponge[idx];
      const tau = 0.5 * (t0 + Math.sqrt(t0 * t0 + Cs2x * Q / rho));
      const sv = 1 / tau;
      const sq = 8 * (2 - sv) / (8 - sv);
      // termes de force en moments (Guo)
      const uF = ux * Fx + uy * Fy;
      // collision
      m[1] = m[1] - S_E * (m[1] - eEq) + (1 - 0.5 * S_E) * 6 * uF;
      m[2] = m[2] - S_EPS * (m[2] - epsEq) + (1 - 0.5 * S_EPS) * (-6 * uF);
      m[3] = mx + Fx;
      m[4] = m[4] - sq * (m[4] - qxEq) + (1 - 0.5 * sq) * (-Fx);
      m[5] = my + Fy;
      m[6] = m[6] - sq * (m[6] - qyEq) + (1 - 0.5 * sq) * (-Fy);
      m[7] = m[7] - sv * (m[7] - pxxEq) + (1 - 0.5 * sv) * 2 * (ux * Fx - uy * Fy);
      m[8] = m[8] - sv * (m[8] - pxyEq) + (1 - 0.5 * sv) * (uy * Fx + ux * Fy);
      m[0] = drho;
      {
        const a0 = m[0] / 9, a1 = m[1] / 36, a2 = m[2] / 36, a3 = m[3] / 6, a4 = m[4] / 12, a5 = m[5] / 6, a6 = m[6] / 12, a7 = m[7] / 4, a8 = m[8] / 4;
        const b = a0 - a1 - 2 * a2, d = a0 + 2 * a1 + a2;
        fB[idx] = a0 - 4 * a1 + 4 * a2;
        fB[N + idx] = b + a3 - 2 * a4 + a7;
        fB[2 * N + idx] = b + a5 - 2 * a6 - a7;
        fB[3 * N + idx] = b - a3 + 2 * a4 + a7;
        fB[4 * N + idx] = b - a5 + 2 * a6 - a7;
        fB[5 * N + idx] = d + a3 + a4 + a5 + a6 + a8;
        fB[6 * N + idx] = d - a3 - a4 + a5 + a6 - a8;
        fB[7 * N + idx] = d - a3 - a4 - a5 - a6 + a8;
        fB[8 * N + idx] = d + a3 + a4 - a5 - a6 - a8;
      }
      // thermique TRT
      let src = 0;
      const z = zone[idx];
      if (z >= 0) {
        src = hk[idx] > 0 ? hk[idx] * (hT[idx] - th) : hq[idx];
        exAcc[z] += src; exT[z] += th;
      }
      const tg = Math.max(taug0 + (tau - tau0) / PRT, 0.5 + 3 * um / PE_CELL);
      const og = 1 / tg;
      for (let k = 0; k < 5; k++) {
        const ge = W5[k] * th * (1 + 3 * (CX[k] * ux + CY[k] * uy));
        gB[k * N + idx] = g[k] - og * (g[k] - ge) + W5[k] * src;
      }
      rhoA[idx] = rho; uxA[idx] = ux; uyA[idx] = uy; thA[idx] = th + src;
    }
    this.f0 = fB; this.f1 = fA; this.g0 = gB; this.g1 = gA;
    this.step++;
    for (let q = 0; q < nf; q++) {
      const fn = fans[q];
      fn.meas = fanAcc[q] / Math.max(1, fn.n);
      const target = ramp * this.fanTargetLat(fn);
      fn.F = Math.max(-0.02, Math.min(0.02, fn.F + this.fanKi(fn) * (target - fn.meas)));
    }
    const kP = this.fl.rho * this.fl.cp * this.dx * this.dx * this.cfg.depth / this.dt;
    for (let z = 0; z < ne; z++) {
      const e = exch[z];
      e.P = e.P * 0.98 + exAcc[z] * kP * 0.02;
      e.Tm = exT[z] / Math.max(1, e.n) + this.cfg.Tamb;
    }
  }

  // met à jour T réel pour l'affichage
  sync() { const Ta = this.cfg.Tamb; for (let i = 0; i < this.N; i++) if (this.type[i] === 0) this.T[i] = this.th[i] + Ta; }

  fanReport(fn) {
    return { Q: fn.meas * this.cu * fn.wperp * this.cfg.depth * 3600, dp: fn.F * fn.thick * this.ps };
  }
  simTime() { return this.step * this.dt; }
}

if (typeof module !== 'undefined') module.exports = { Solver, FLUIDS, DIRS };
