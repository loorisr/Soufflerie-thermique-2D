// ===== Portage WebGPU du solveur (mêmes équations que solver.js) =====
const WGSL = /* wgsl */`
struct Params {
  N: u32, nx: u32, nf: u32, ne: u32,
  tau0: f32, taug0: f32, gb: f32, cs2x: f32,
  ramp: f32, peCell: f32, sE: f32, sEps: f32,
  gb2: f32, pad0: f32, pad1: f32, pad2: f32,
};
@group(0) @binding(0) var<uniform> P: Params;
// fg : populations f (0..8) puis g (9..13), en double tampon
@group(0) @binding(1) var<storage, read> A: array<f32>;
@group(0) @binding(2) var<storage, read_write> B: array<f32>;
@group(0) @binding(3) var<storage, read_write> mac: array<vec4<f32>>;
@group(0) @binding(4) var<storage, read> nb: array<i32>;
// cell : [drag, hk, hT, hq], [sponge, Tw, -, -], [type, zone, fan, reset]
@group(0) @binding(5) var<storage, read> cell: array<vec4<f32>>;
// st : fanCfg 0..31 (2 par ventilateur), fanState 32..47, exState 48..63
@group(0) @binding(6) var<storage, read_write> st: array<vec4<f32>>;
@group(0) @binding(7) var<storage, read_write> acc: array<atomic<i32>, 48>;

const FAN_SCALE: f32 = 100000.0;
const HEAT_SCALE: f32 = 20000.0;
const T_SCALE: f32 = 100.0;

var<private> CX: array<f32, 9> = array<f32, 9>(0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0);
var<private> CY: array<f32, 9> = array<f32, 9>(0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0);
var<private> W9: array<f32, 9> = array<f32, 9>(0.444444444, 0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.0277777778, 0.0277777778, 0.0277777778, 0.0277777778);
var<private> OPP: array<u32, 9> = array<u32, 9>(0u, 3u, 4u, 1u, 2u, 7u, 8u, 5u, 6u);
var<private> W5: array<f32, 5> = array<f32, 5>(0.333333333, 0.166666667, 0.166666667, 0.166666667, 0.166666667);

@compute @workgroup_size(64)
fn step(@builtin(global_invocation_id) gid: vec3<u32>) {
  let idx = gid.x;
  let N = P.N;
  if (idx >= N) { return; }
  let ci = vec4<i32>(cell[3u * idx + 2u]);
  if (ci.x != 0) { return; }
  let c0 = cell[3u * idx];
  let c1 = cell[3u * idx + 1u];
  let mc = mac[idx];
  let u0x = mc.y; let u0y = mc.z; let thc = mc.w;
  let u0sq = u0x * u0x + u0y * u0y;
  var f: array<f32, 9>;
  var g: array<f32, 5>;
  for (var k = 0u; k < 9u; k++) {
    let v = nb[k * N + idx];
    if (v >= 0) { f[k] = A[k * N + u32(v)]; }
    else if (v == -2) {
      let cu = CX[k] * u0x + CY[k] * u0y;
      f[k] = -A[OPP[k] * N + idx] + 2.0 * W9[k] * (4.5 * cu * cu - 1.5 * u0sq);
    } else { f[k] = A[OPP[k] * N + idx]; }
  }
  for (var k = 0u; k < 5u; k++) {
    let v = nb[k * N + idx];
    if (v >= 0) { g[k] = A[(9u + k) * N + u32(v)]; }
    else if (v == -2) {
      let cu = CX[k] * u0x + CY[k] * u0y;
      var tin = thc;
      if (cu > 0.0) { tin = 0.0; }
      g[k] = W5[k] * tin * (1.0 + 3.0 * cu);
    } else if (v == -1) { g[k] = A[(9u + OPP[k]) * N + idx]; }
    else {
      let s = u32(-3 - v);
      g[k] = -A[(9u + OPP[k]) * N + idx] + 2.0 * W5[k] * cell[3u * s + 1u].y;
    }
  }
  // moments MRT
  let sa = f[1] + f[2] + f[3] + f[4];
  let sd = f[5] + f[6] + f[7] + f[8];
  var m0 = f[0] + sa + sd;
  var m1 = -4.0 * f[0] - sa + 2.0 * sd;
  var m2 = 4.0 * f[0] - 2.0 * sa + sd;
  let mx = f[1] - f[3] + f[5] - f[6] - f[7] + f[8];
  var m4 = -2.0 * f[1] + 2.0 * f[3] + f[5] - f[6] - f[7] + f[8];
  let my = f[2] - f[4] + f[5] + f[6] - f[7] - f[8];
  var m6 = -2.0 * f[2] + 2.0 * f[4] + f[5] + f[6] - f[7] - f[8];
  var m7 = f[1] - f[2] + f[3] - f[4];
  var m8 = f[5] - f[6] + f[7] - f[8];
  var drho = m0;
  var rho = 1.0 + drho;
  if (!(rho > 0.2 && rho < 5.0)) { rho = 1.0; drho = 0.0; }
  let th = g[0] + g[1] + g[2] + g[3] + g[4];
  // forces
  var Fx = 0.0;
  var Fy = rho * th * (P.gb + P.gb2 * th);
  let fq = ci.z;
  var fd = vec2<f32>(0.0, 0.0);
  if (fq >= 0) {
    let fc = st[2u * u32(fq)];
    fd = fc.xy;
    let F = st[32u + u32(fq)].x;
    Fx += rho * F * fd.x; Fy += rho * F * fd.y;
  }
  var ux = (mx + 0.5 * Fx) / rho;
  var uy = (my + 0.5 * Fy) / rho;
  let C = c0.x;
  if (C > 0.0) {
    let s = 0.5 * C * sqrt(ux * ux + uy * uy) + 0.002;
    ux = ux / (1.0 + s); uy = uy / (1.0 + s);
    Fx -= 2.0 * rho * s * ux; Fy -= 2.0 * rho * s * uy;
  }
  let um0 = sqrt(ux * ux + uy * uy);
  if (um0 > 0.3) { ux *= 0.3 / um0; uy *= 0.3 / um0; }
  let um = min(um0, 0.3);
  if (fq >= 0) {
    atomicAdd(&acc[u32(fq)], i32(round((ux * fd.x + uy * fd.y) * FAN_SCALE)));
  }
  let usq = ux * ux + uy * uy;
  let eEq = -2.0 * drho + 3.0 * rho * usq;
  let epsEq = drho - 3.0 * rho * usq;
  let pxxEq = rho * (ux * ux - uy * uy);
  let pxyEq = rho * ux * uy;
  let nE = m1 - eEq; let nXX = m7 - pxxEq; let nXY = m8 - pxyEq;
  let Pxx = (nE / 3.0 + nXX) * 0.5; let Pyy = (nE / 3.0 - nXX) * 0.5;
  let Q = sqrt(Pxx * Pxx + Pyy * Pyy + 2.0 * nXY * nXY);
  let t0 = P.tau0 + c1.x;
  let tau = 0.5 * (t0 + sqrt(t0 * t0 + P.cs2x * Q / rho));
  let sv = 1.0 / tau;
  let sq = 8.0 * (2.0 - sv) / (8.0 - sv);
  let uF = ux * Fx + uy * Fy;
  m1 = m1 - P.sE * (m1 - eEq) + (1.0 - 0.5 * P.sE) * 6.0 * uF;
  m2 = m2 - P.sEps * (m2 - epsEq) + (1.0 - 0.5 * P.sEps) * (-6.0 * uF);
  let m3 = mx + Fx;
  m4 = m4 - sq * (m4 + rho * ux) + (1.0 - 0.5 * sq) * (-Fx);
  let m5 = my + Fy;
  m6 = m6 - sq * (m6 + rho * uy) + (1.0 - 0.5 * sq) * (-Fy);
  m7 = m7 - sv * (m7 - pxxEq) + (1.0 - 0.5 * sv) * 2.0 * (ux * Fx - uy * Fy);
  m8 = m8 - sv * (m8 - pxyEq) + (1.0 - 0.5 * sv) * (uy * Fx + ux * Fy);
  let a0 = drho / 9.0; let a1 = m1 / 36.0; let a2 = m2 / 36.0; let a3 = m3 / 6.0; let a4 = m4 / 12.0;
  let a5 = m5 / 6.0; let a6 = m6 / 12.0; let a7 = m7 / 4.0; let a8 = m8 / 4.0;
  let b = a0 - a1 - 2.0 * a2; let d = a0 + 2.0 * a1 + a2;
  B[idx] = a0 - 4.0 * a1 + 4.0 * a2;
  B[N + idx] = b + a3 - 2.0 * a4 + a7;
  B[2u * N + idx] = b + a5 - 2.0 * a6 - a7;
  B[3u * N + idx] = b - a3 + 2.0 * a4 + a7;
  B[4u * N + idx] = b - a5 + 2.0 * a6 - a7;
  B[5u * N + idx] = d + a3 + a4 + a5 + a6 + a8;
  B[6u * N + idx] = d - a3 - a4 + a5 + a6 - a8;
  B[7u * N + idx] = d - a3 - a4 - a5 - a6 + a8;
  B[8u * N + idx] = d + a3 + a4 - a5 - a6 - a8;
  // thermique
  var src = 0.0;
  let z = ci.y;
  if (z >= 0) {
    if (c0.y > 0.0) { src = c0.y * (c0.z - th); } else { src = c0.w; }
    atomicAdd(&acc[16u + u32(z)], i32(round(src * HEAT_SCALE)));
    atomicAdd(&acc[32u + u32(z)], i32(round(th * T_SCALE)));
  }
  let tg = max(P.taug0 + (tau - P.tau0) / 0.85, 0.5 + 3.0 * um / P.peCell);
  let og = 1.0 / tg;
  for (var k = 0u; k < 5u; k++) {
    let ge = W5[k] * th * (1.0 + 3.0 * (CX[k] * ux + CY[k] * uy));
    B[(9u + k) * N + idx] = g[k] - og * (g[k] - ge) + W5[k] * src;
  }
  mac[idx] = vec4<f32>(rho, ux, uy, th + src);
}

@compute @workgroup_size(32)
fn control(@builtin(local_invocation_id) lid: vec3<u32>) {
  let t = lid.x;
  let stepc = st[63].x;
  let ramp = min(1.0, stepc / 3000.0);
  if (t < 16u) {
    if (t < P.nf) {
      let c0 = st[2u * t]; let c1 = st[2u * t + 1u];
      let meas = f32(atomicExchange(&acc[t], 0)) / FAN_SCALE * c0.z;
      var s = st[32u + t];
      s.x = clamp(s.x + c0.w * (ramp * c1.x - meas), -0.02, 0.02);
      s.y = meas;
      st[32u + t] = s;
    }
  } else {
    let e = t - 16u;
    if (e < P.ne) {
      let v = f32(atomicExchange(&acc[16u + e], 0)) / HEAT_SCALE;
      let tt = f32(atomicExchange(&acc[32u + e], 0)) / T_SCALE;
      var s = st[48u + e];
      s.x = s.x * 0.98 + v * 0.02;
      s.y = tt * s.z;
      st[48u + e] = s;
    }
  }
  workgroupBarrier();
  if (t == 31u) { st[63].x = stepc + 1.0; }
}

@compute @workgroup_size(64)
fn resetCells(@builtin(global_invocation_id) gid: vec3<u32>) {
  let idx = gid.x;
  let N = P.N;
  if (idx >= N) { return; }
  if (cell[3u * idx + 2u].w == 0.0) { return; }
  for (var k = 0u; k < 14u; k++) { B[k * N + idx] = 0.0; }
  mac[idx] = vec4<f32>(1.0, 0.0, 0.0, 0.0);
}
`;

class GpuSolver extends Solver {
  static async create(cfg, comps) {
    if (!navigator.gpu) throw new Error('WebGPU indisponible');
    const adapter = await navigator.gpu.requestAdapter({ powerPreference: 'high-performance' });
    if (!adapter) throw new Error('Aucun adaptateur WebGPU');
    const need = 14 * 4 * cfg.nx * Math.max(8, Math.round(cfg.nx * cfg.H / cfg.L));
    const lim = adapter.limits;
    if (need > lim.maxStorageBufferBindingSize) throw new Error('Grille trop grande pour ce GPU');
    const dev = await adapter.requestDevice({ requiredLimits: {
      maxStorageBufferBindingSize: Math.max(need, Math.min(lim.maxStorageBufferBindingSize, 134217728)),
      maxBufferSize: Math.max(need, Math.min(lim.maxBufferSize, 268435456)),
    } });
    const s = new GpuSolver(cfg, comps);
    s.initGpu(dev);
    return s;
  }

  initGpu(dev) {
    this.dev = dev;
    const N = this.N, B = GPUBufferUsage;
    const mk = (size, usage) => dev.createBuffer({ size: Math.max(16, Math.ceil(size / 16) * 16), usage });
    const S = B.STORAGE | B.COPY_DST | B.COPY_SRC;
    this.bufs = {
      params: mk(64, B.UNIFORM | B.COPY_DST),
      fg: [mk(14 * N * 4, S), mk(14 * N * 4, S)],
      mac: mk(4 * N * 4, S),
      nb: mk(9 * N * 4, S),
      cell: mk(12 * N * 4, S),
      st: mk(64 * 16, S),
      acc: mk(48 * 4, S),
      readMac: mk(4 * N * 4, B.MAP_READ | B.COPY_DST),
      readStat: mk(64 * 16, B.MAP_READ | B.COPY_DST),
    };
    const mod = dev.createShaderModule({ code: WGSL });
    this.shaderInfo = mod.getCompilationInfo ? mod.getCompilationInfo() : null;
    const layout = dev.createBindGroupLayout({ entries: [
      { binding: 0, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'uniform' } },
      { binding: 1, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
      { binding: 2, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },
      { binding: 3, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },
      { binding: 4, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
      { binding: 5, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'read-only-storage' } },
      { binding: 6, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },
      { binding: 7, visibility: GPUShaderStage.COMPUTE, buffer: { type: 'storage' } },
    ] });
    const pl = dev.createPipelineLayout({ bindGroupLayouts: [layout] });
    this.pipeStep = dev.createComputePipeline({ layout: pl, compute: { module: mod, entryPoint: 'step' } });
    this.pipeCtrl = dev.createComputePipeline({ layout: pl, compute: { module: mod, entryPoint: 'control' } });
    this.pipeReset = dev.createComputePipeline({ layout: pl, compute: { module: mod, entryPoint: 'resetCells' } });
    const b = this.bufs;
    this.bg = [0, 1].map(p => dev.createBindGroup({ layout, entries: [
      { binding: 0, resource: { buffer: b.params } },
      { binding: 1, resource: { buffer: b.fg[p] } }, { binding: 2, resource: { buffer: b.fg[1 - p] } },
      { binding: 3, resource: { buffer: b.mac } }, { binding: 4, resource: { buffer: b.nb } },
      { binding: 5, resource: { buffer: b.cell } }, { binding: 6, resource: { buffer: b.st } },
      { binding: 7, resource: { buffer: b.acc } },
    ] }));
    this.parity = 0;
    this.busy = false;
    this.lost = false;
    dev.lost.then(() => { this.lost = true; });
    this.uploadAll(true);
  }

  // super.setComponents est appelé par le constructeur avant l'init GPU
  setComponents(comps) {
    super.setComponents(comps);
    if (this.dev) this.uploadAll(false);
  }

  uploadAll(first) {
    const q = this.dev.queue, b = this.bufs, N = this.N;
    q.writeBuffer(b.nb, 0, this.nb);
    const cf = new Float32Array(12 * N);
    for (let i = 0; i < N; i++) {
      const o = 12 * i;
      cf[o] = this.drag[i]; cf[o + 1] = this.hk[i]; cf[o + 2] = this.hT[i]; cf[o + 3] = this.hq[i];
      cf[o + 4] = this.sponge[i]; cf[o + 5] = this.Tw[i];
      cf[o + 8] = this.type[i]; cf[o + 9] = this.zone[i]; cf[o + 10] = this.fan[i]; cf[o + 11] = first ? 1 : this.reset[i];
    }
    q.writeBuffer(b.cell, 0, cf);
    const ss = new Float32Array(256);
    this.fans.forEach((fn, i) => { ss[4 * (32 + i)] = fn.F; ss[4 * (32 + i) + 1] = fn.meas; });
    this.exch.forEach((e, i) => { ss[4 * (48 + i)] = e.Plat || 0; ss[4 * (48 + i) + 2] = 1 / Math.max(1, e.n); });
    ss[4 * 63] = this.step;
    q.writeBuffer(b.st, 0, ss);
    q.writeBuffer(b.acc, 0, new Int32Array(48));
    this.writeFanCfg();
    this.writeParams();
    // remise à zéro des cellules redevenues fluides (dans les deux tampons)
    const enc = this.dev.createCommandEncoder();
    const pass = enc.beginComputePass();
    pass.setPipeline(this.pipeReset);
    const wg = Math.ceil(N / 64);
    pass.setBindGroup(0, this.bg[0]); pass.dispatchWorkgroups(wg);
    pass.setBindGroup(0, this.bg[1]); pass.dispatchWorkgroups(wg);
    pass.end();
    q.submit([enc.finish()]);
  }

  writeFanCfg() {
    const fc = new Float32Array(128);
    this.fans.forEach((fn, i) => {
      fc[8 * i] = fn.dx; fc[8 * i + 1] = fn.dy; fc[8 * i + 2] = 1 / Math.max(1, fn.n); fc[8 * i + 3] = this.fanKi(fn);
      fc[8 * i + 4] = this.fanTargetLat(fn); fc[8 * i + 5] = fn.thick;
    });
    this.dev.queue.writeBuffer(this.bufs.st, 0, fc);
  }

  writeParams() {
    const buf = new ArrayBuffer(64), u = new Uint32Array(buf), f = new Float32Array(buf);
    u[0] = this.N; u[1] = this.nx; u[2] = this.fans.length; u[3] = this.exch.length;
    f[4] = this.tau0; f[5] = this.taug0; f[6] = this.gb; f[7] = 18 * Math.SQRT2 * CS * CS;
    f[8] = this.ramp(); f[9] = PE_CELL; f[10] = S_E; f[11] = S_EPS; f[12] = this.gb2;
    this.dev.queue.writeBuffer(this.bufs.params, 0, buf);
  }

  // n pas sur le GPU puis relecture des champs ; résout quand les données sont prêtes
  async runAsync(n) {
    if (this.busy || this.lost) return false;
    this.busy = true;
    const dev = this.dev, b = this.bufs, N = this.N;
    this.writeParams();
    this.writeFanCfg();
    const enc = dev.createCommandEncoder();
    const pass = enc.beginComputePass();
    const wg = Math.ceil(N / 64);
    for (let s = 0; s < n; s++) {
      pass.setPipeline(this.pipeStep);
      pass.setBindGroup(0, this.bg[this.parity]);
      pass.dispatchWorkgroups(wg);
      pass.setPipeline(this.pipeCtrl);
      pass.dispatchWorkgroups(1);
      this.parity ^= 1;
    }
    pass.end();
    enc.copyBufferToBuffer(b.mac, 0, b.readMac, 0, 16 * N);
    enc.copyBufferToBuffer(b.st, 0, b.readStat, 0, 1024);
    dev.queue.submit([enc.finish()]);
    this.step += n;
    try {
      await Promise.all([b.readMac.mapAsync(GPUMapMode.READ), b.readStat.mapAsync(GPUMapMode.READ)]);
    } catch (e) { this.busy = false; this.lost = true; return false; }
    const m = new Float32Array(b.readMac.getMappedRange());
    const Ta = this.cfg.Tamb;
    for (let i = 0; i < N; i++) {
      if (this.type[i] !== 0) continue;
      this.rho[i] = m[4 * i]; this.ux[i] = m[4 * i + 1]; this.uy[i] = m[4 * i + 2];
      this.th[i] = m[4 * i + 3]; this.T[i] = m[4 * i + 3] + Ta;
    }
    b.readMac.unmap();
    const st = new Float32Array(b.readStat.getMappedRange());
    const kP = this.fl.rho * this.fl.cp * this.dx * this.dx * this.cfg.depth / this.dt;
    this.fans.forEach((fn, i) => { fn.F = st[4 * (32 + i)]; fn.meas = st[4 * (32 + i) + 1]; });
    this.exch.forEach((e, i) => { e.Plat = st[4 * (48 + i)]; e.P = e.Plat * kP; e.Tm = st[4 * (48 + i) + 1] + Ta; });
    b.readStat.unmap();
    this.busy = false;
    return true;
  }

  run(n) { return this.runAsync(n); }
  destroy() { try { Object.values(this.bufs).flat().forEach(x => x.destroy()); this.dev.destroy(); } catch (e) {} }
}
