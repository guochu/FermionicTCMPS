# FermionicTCMPS 复活瘦身计划

**日期**: 2026-11-19 · **状态**: rev 4 待实施（D1''/D5 已由用户裁决关闭）
**目标**: 复活 FermionicTCMPS 为 GTEMPO 的独立对照实现——星型 thermofield AIM 的可分离初态直接实时演化，TensorKit U₁×U₁ 对称张量栈。

---

## 0. 需求裁决史

**rev 1 五条原需求**: ① 零温/有限温（thermofield）实时演化、初态 = 杂质⊗热库热平衡态、H 不含时 ② U₁×U₁ 对称、全 JW 玻色表象 ③ 只 TDVP2 ④ 单杂质(up/down)+1/2 bath（1 bath 杂质 site 1 / 2 bath 居中）⑤ TensorKit + ImpurityModelBase。

**rev 2 四点裁决**: ① chainmapping 全删（有限温默认 thermofield 星型）② 不求整体基态（`ground_state`/`gs_*` 删）③ quench 全删 ④ `thermal_state` 删——只做可分离初态直接实时演化 ⑤ MPO 构造弃 GeneralHamiltonians 字符串求和，改 **FSM 通道式直构**。

**rev 3 两点裁决**: ⑥ 引擎包只取极小子集，"多余的都去掉" ⑦ 离散化完全复用 ImpurityModelBase `DiscreteSpectrum`。

**rev 4 终局裁决（2026-11-19）**: ⑧ **"只参考 FiniteMPSAlgorithms 的思想，本包用 U₁×U₁ 对称性张量"**——FMA 稠密栈整体出局（不 vendor、不依赖），只继承其 FSM/通道布局思想与 TDVP 驱动模式；对称版 `SchurMPOTensor` 的实现参考 = `backup/GTEMPO/src/mpo/mpohamiltonian/schurmpotensor.jl`。

---

## 1. 架构（终态）

### 1.1 依赖栈

```
FermionicTCMPS
  ├── TensorKit          （U₁×U₁ 对称张量）
  ├── DMRG               （TensorKit MPS 引擎——极小使用面，见 §1.2）
  ├── ImpurityModelBase  （离散浴 + ED 参照；dev ../../ImpurityModelBase）
  └── LinearAlgebra
```

退役：GeneralHamiltonians、TEBD、InfiniteDMRG、FiniteMPSAlgorithms（连同 rev 3 设想的 vendor）、QuenchScheme/chainmapping/thermal_state/ground_state 全套功能。

### 1.2 DMRG 包的极小使用面（不裁剪仓库，只约束调用面）

DMRG 包**已导出全部所需容器与算法**（tmp/DMRG 逐条核实）：

| 符号 | 用途 |
|------|------|
| `SchurMPOTensor` / `MPOHamiltonian` / `SchurMPO` / `fromABCD` | **FSM MPO 的目标容器**：上三角算符矩阵 + 逐行/列 `ElementarySpace` 键空间（通道各带 U₁×U₁ 扇区）；`fromABCD(;B, C, D, A)` 即角块式 FSM 构造入口（B=字串起点向量、C=闭合向量、D=on-site、A=传播矩阵） |
| `environments` / `ExpectationCache` / `TDVP2` / `sweep!` | TDVP2 演化（`TDVPCache` 模式：`sweep!(environments(h, ψ), TDVP2(stepsize=-im*δt, trunc=...))`；rev-1 代码已验证 `MPOHamiltonian` 直入此链） |
| `prodmps` / `prodmpo` / `MPO` | 积态构造、可观测量小 MPO |
| `truncdimcutoff` / `expectation` / `expectationvalue` | 截断与观测 |

**明确不用**：`DMRG2`/`ground_state`（无基态需求）、`schurmpo/exponentialexpansion` 族（TEMPO 遗产）、`SparseGMPS` 等。

### 1.3 浴配置（rev 3 裁决不变）

```julia
bath = discretefermionicbath(ws, fs; β, μ)      # β=∞: discretefermionicvacuum(ws, fs; μ)
model = SISBD(bath; U, εd, J)                   # 有限温内部默认 thermofield 倍增
```

自家 discretization.jl 删；`StarConfiguration`/`StarThermofieldConfiguration` 类塌缩——模型直持 `AbstractDiscreteNormalBath`，`thermofield` 只剩热占据劈裂纯函数（`V → V√(1−n)` 上 a₁ 支、`V√n` 上 a₂ 支，a₁/a2 位点交错排布）。

### 1.4 FSM MPO 直构（核心新增：`src/siteops.jl` + `src/mpoconstruction.jl`）

**FMA 思想 + DMRG 容器 + GTEMPO 参考**的落点：

```
H = (εd·n + U·n↑n↓)|杂质 + Σ_j ω_j n_j + Σ_j V_j' (c_imp† c_j + h.c.)    V_j' = J·V_j（有限温经 thermofield 劈裂）
```

1. **局域算符**（`siteops.jl`，自写 ~80 行）：物理空间 `Rep[U₁×U₁]((0,0),(1,0),(0,1),(1,1))`；`d↑`/`d↓`（intra-site JW 相位折入，约定沿 GeneralHamiltonians `sp_u1u1`——rev-1 全套测试曾机器精度验证）、`P`（宇称内蕴算符，charge (0,0)）、`n↑`/`n↓`/`n↑n↓`；
2. **通道布局**（FMA `PairOpTerm` 思想 + 用户裁决的系数放置规则）：**系数恒放非共享端**（杂质是所有跳变项的共享端）——
   - 右浴 `c_imp†c_j`：杂质开通道，浴位 j 闭合、系数 `C[1,j]·d_j` 在闭合张量（FMA 原生约定）；
   - 右浴 h.c. `c_j†c_imp` 与左浴两方向：闭合端都是杂质（共享）→ 系数移到开启端 `V_j·d_j^{(†)}`（每浴位自己的开启块），闭合统一为杂质的 `d^{(†)}`；
   - **每方向一个通道**：星型耦合系数只依赖非杂质端（可分离），通道内线性累积 `Σ_j V_j·(字串)·d_j`，闭合算符全同——单通道合并成立（FMA `PairOpTerm` 之所以每起点一通道，是其任意耦合矩阵 C[i,j] 依赖两端、系数被迫放闭合端，本包不用该约定）；
   - 整链 MPO 键维 ~6 恒定（真空 + 4 方向通道 + 闭合），与浴大小无关；
3. **落点类型**：每站点一个 `DMRG.SchurMPOTensor`（通道键空间：真空/闭合/通道 = 平凡扇区——宇称串 charge-0 不改键荷）→ `SchurMPO(链)` = `MPOHamiltonian`；
4. **双浴**（杂质居中）：左右浴同一构造器，仅系数放置端不同；**无镜像布局**（方案 b = 系数放置规则，不是布局翻转）；
5. **实现参考**：`backup/GTEMPO` 同款 `SchurMPOTensor`（`MPSBondTensor` 恒等/宇称通道的 `isomorphism` 提升技巧、`isid` 标量化、`leftspaces/rightspaces` 记账）；
6. 构造 O(N)：对比旧 GeneralHamiltonians 路线（N 项 × O(N) 字符串 QTerm + simplify = O(N²)）。

### 1.5 演化与观测

- **初态**：`prodmps`（杂质位指定占据 + a₁ 位全 0 + a₂ 位全 1——热积态；a₂ 支的 PHT 约定在 G2 用稠密锚钉死）；
- **演化**：`evolve!(ψ, h, tspan; δt, trunc)` 薄驱动 ~30 行（TDVPCache + sweep! 循环，替代旧 TEBD 胶水）；
- **观测**：n_d/J/单点算符 → 单点张量嵌入恒等 MPO 走 `expectation(ψ, h)`；实时 GF（greater/lesser）= 演化到 t → 局域算符单点插入（键维不变）→ 续演化 → 单点 overlap；
- **U₁×U₁**：张量级扇区分块全程生效（需求 2 完整恢复）；测试另断言 ⟨N↑⟩/⟨N↓⟩ 守恒。

### 1.6 目标公开 API 面（export 草案，~8 符号）

```julia
export SISBD, SIDBD, mpo, separable_state, evolve!
export occupation_op, current_op, gf_greater_lesser_t
```

### 1.7 目标规模

src 2741 → **~600–900 行**（无 vendor）；导出 33 → 8。

---

## 2. 文件级判决表（rev 4 终版）

| 文件 | 判决 |
|------|------|
| `src/auxiliary/discretization.jl` | **删**（IMB `DiscreteSpectrum` 复用） |
| `src/auxiliary/freegf.jl`、`rungekutta.jl` | **删** |
| `src/operators/*` | **删**（`FermionicHamiltonian`/Term 体系被 FSM 直构替代） |
| `src/tcmps/star_to_chain.jl`、`quenchscheme.jl` | **删** |
| `src/tcmps/transformations.jl` | **重写**：只留 thermofield 劈裂纯函数 |
| `src/tcmps/abstractmodels.jl`、`current.jl`、`separablestate.jl` | **重写**（骨架/小 MPO 观测/prodmps 积态） |
| `src/tcmps/{one,two}bath/{oneimpurity,twoimpurity}.jl` | **重写**：SISBD/SIDBD |
| `{one,two}bath/{thermalstate,freefermions}.jl`、`gf/*` | **删** |
| `src/siteops.jl` | **新增**：U₁×U₁ 局域算符表 |
| `src/mpoconstruction.jl` | **新增**：FSM 直构 → `DMRG.SchurMPO` |
| `src/timeevo.jl` | **新增**：TDVP2 薄驱动 |
| 测试 570 行 | **重写**：G2/G3 三重锚 |

---

## 3. 实施步骤（五 Phase，~5–7 天）

### Phase 0 — 环境底座（~0.5 天）
1. `git checkout -b revival-slim`；`Project.toml` 重写：deps = `TensorKit` + `DMRG`(dev `tmp/DMRG`) + `ImpurityModelBase`(dev `../../ImpurityModelBase`) + `LinearAlgebra`；`test/Project.toml` 同源（`--project=test` 跑法）。
2. `tmp/GeneralHamiltonians`、`tmp/TEBD` 退役（G3 后删）；`tmp/FiniteMPSAlgorithms` 保留仅作 FMA 思想 diff 参照（无代码拷贝）。

**闸门 G0**: `--project=test` 全量加载成功。

### Phase 1 — 瘦身手术（~0.5 天）
按 §2 执行；`FermionicTCMPS.jl` 重写（export 8 符号）。

**闸门 G1**: 加载成功、旧符号不可见。

### Phase 2 — siteops + FSM MPO（~1.5 天）
1. `siteops.jl`：U₁×U₁ 局域算符（约定沿 `sp_u1u1`）；`mpoconstruction.jl`：统一 FSM 构造器（系数放非共享端），SISBD/SIDBD → `SchurMPOTensor` 链 → `SchurMPO`。
2. **稠密锚**（红线，防自证）：FSM MPO 的 Schur 链手工收缩展开成 4^N 稠密矩阵 vs 测试内逐项 Kron 装配的哈密顿量——**逐位一致**（JW 相位锚）；四配置覆盖（单/双浴 × 真空/thermofield）；双浴算例须含左浴两方向的系数-开启端路径；**thermofield a₂ 支的 PHT 约定在此钉死**。

**闸门 G2**: 稠密锚全绿。

### Phase 3 — 初态、演化、观测（~1.5 天）
`separable_state` + `evolve!` + n_d/J/GF；**自由模型锚**（U=0 vs IMB `freefermions_greater_lesser`/`freefermions_Gt` 机器精度）+ δt² 收敛阶 + ⟨N↑⟩/⟨N↓⟩ 守恒 + β→∞ 退化（a₂ 支退耦不 NaN）。

**闸门 G3**: 全套绿（验收闸门，不过全停）。

### Phase 4 — GTEMPO 对拍（~2–3 天，继承 rev 1–3 设计）
参数桥双侧同源字面量 → n_d(t)/J(t)/G^>_d/G^<_d；零温+有限温（β∈{0.5,1}）× 单浴+双浴；容差 ≤1e-3（积分量）/≤1e-2（逐点 GF）。`test/gtempo_crosscheck.jl`（秒级）+ `scripts/gtempo_benchmark.jl`（CSV 扫描）。

**闸门 G4**: 四组合各一算例过容差。

### Phase 5 — 文档收尾（~0.5 天）
README（极小范围/依赖/跑法/引用 2012.01424）；docstring 统一；推送 `revival-slim` 开 PR。

---

## 4. 已关闭裁决与遗留裁决

| # | 裁决 | 状态 |
|---|------|------|
| D1'' | 后端 = TensorKit U₁×U₁ 对称张量；FMA 只作思想参考 | **✅ 用户已裁（rev 4）** |
| D5 | vendor 机制 | **✅ 作废**（不 vendor；DMRG 容器直接可用） |
| D2' | 左 bath 通道方案 | **✅ 用户已裁：方案 b**——系数恒放非共享端，每方向单通道，整链 MPO 键维 ~6 恒定（无镜像布局） |
| rev2 ①–⑤、rev3 ⑥–⑦ | chainmapping/基态/quench/热态全删；极小使用面；IMB DiscreteSpectrum | **✅ 已裁** |
| D4' | GF 面 = 实时 greater/lesser + n_d + J | 待确认（默认采纳） |

---

## 5. 风险与开放问题

1. **PHT 约定**：thermofield a₂ 支的粒子- hole 变换（rev-1 `couplings` docstring 提及"particle-hole transformation on top of thermofield"）的符号/产生湮灭互换细节，必须在 G2 稠密锚逐位钉死，不得臆测。
2. **DMRG `MPOHamiltonian`×TDVP2 链**：rev-1 代码已跑通（`environments(h, ψ)` + `sweep!`），但当时 h 来自 `consolidate`；FSM 直构的 `SchurMPOTensor` 链在 Phase 2 冒烟中先行验证通道键空间记账（`leftspaces`/`rightspaces` 逐链一致）。
3. **系数-开启端路径的自证风险**：方案 b（系数放开启端、单通道线性累积）是 FMA `PairOpTerm` 约定（系数放闭合端）之外的路线，无现成实现对照——正确性完全依赖 G2 稠密锚；双浴四方向（左右浴 × 跳变/h.c.）逐一覆盖，不许只测右浴。
4. **U₁×U₁ 局域算符**：`d↑`/`d↓` 的 intra-site JW 相位与 `P` 的扇区对角值，沿 `sp_u1u1` 定义抄写并在 G2 锚定；`P` 为 charge-(0,0) 内蕴算符（宇称串不改变通道键荷）。
5. **β→∞ 退化**：`√(1−n_k)`→1、`√n_k`→0，a₂ 支退耦——覆盖测试。
6. **DMRG 版本**：dev 固定 `tmp/DMRG` 快照（f141e43 时代的 HEAD `a0f239c`）；上游变更手动同步。
