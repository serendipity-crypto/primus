# 算法选型与未来改进方向

本文档说明 `primus_gcd` 当前使用的算法、它们的取舍，以及如果将来需要优化时可以参考的更先进算法。不公开在 rustdoc 中，是给维护者的内部备忘。

---

## 当前实现

| 方法 | 算法 | 复杂度 | 来源 |
|---|---|---|---|
| `gcd` | Stein's binary GCD（shift-and-subtract，无除法）| O(BITS²) | 标准实现 |
| `xgcd` | FLINT 风格 Euclidean xgcd（quot=1/2/3 快路径 + 通用 `x / y`）| O(BITS²) | [FLINT `n_xgcd`](https://flintlib.org/doc/ulong_extras.html#c.n_xgcd) |
| `gcdinv` | 同上，只追踪一个 cofactor | O(BITS²) | [FLINT `n_gcdinv`](https://flintlib.org/doc/ulong_extras.html#c.n_gcdinv) |

**这个组合的取舍**：
- 实现简单、移植自经过工程验证的 FLINT 代码、对 u64 在现代 CPU 上性能可接受。
- 但 `xgcd`/`gcdinv` 仍然包含 `x / y` 这种除法：u64 没问题（硬件 div），u128 上会落到 `__udivti3` 软件除法。
- 控制流依赖输入值（quot=1/2/3 vs 通用分支），**不是常数时间**。如果将来涉及 secret-dependent inversion，需要替换。

---

## 更优算法（按场景分）

### 1. Bernstein–Yang "divsteps2"（2019）

固定迭代次数的 shift-and-subtract 风格 xgcd / modular inversion，**完全 branch-free、constant-time、无除法**。

- **何时考虑**：`gcdinv` 进入 FHE 热路径（key switching / NTT primitive root setup / Barrett `mu` 预计算之外的运行时求逆），或将来需要 side-channel 抗性。
- **现代 CPU 上对 u64 比 Euclidean 风格快约 1.5–3×**；对 u128 优势更大（消掉 `__udivti3`）。
- **参考**：Bernstein, Yang. *Fast constant-time gcd computation and modular inversion*. IACR TCHES 2019. <https://gcd.cr.yp.to/papers.html#safegcd>
- **Rust 参考实现**：[`crypto-bigint::modular::safegcd`](https://github.com/RustCrypto/crypto-bigint/blob/master/src/modular/safegcd.rs) —— 可以直接对照成熟的 Rust 实现。

### 2. Pornin "Optimized Binary GCD for Modular Inversion"（2020）

在 divsteps2 上做了 signed-digit + 多步合并优化。`libsecp256k1` / `blst` 等密码库的 modular inverse 走的就是这个变体。

- **何时考虑**：和 #1 类似，但实现更紧凑；对单字 u64 仍然比 Euclidean 风格快。
- **参考**：Pornin. *Optimized Binary GCD for Modular Inversion*. 2020. <https://eprint.iacr.org/2020/972>
- **参考实现**：[`libsecp256k1/src/modinv64_impl.h`](https://github.com/bitcoin-core/secp256k1/blob/master/src/modinv64_impl.h) —— C 代码，注释非常详细。

### 3. Lehmer / Half-GCD（HGCD）

**只对多字大整数（≥ 4 limbs）有意义**：取被除数 / 除数的高位 limb 当 u64 跑一段 xgcd 预测多步商，再一次性回代到完整大整数，把 O(n²) 降到 O(M(n) log n)。

- **何时考虑**：将来若 `BigUint<S>` 上要实现 xgcd / 模逆。单字 `Xgcd` trait 用不上。
- **参考**：GMP `mpn/generic/hgcd*.c`、FLINT `fmpz_gcd` / `fmpz_xgcd`。
- 注意：HGCD 的实现复杂度远高于单字算法，迁移成本不低。

### 4. 小幅优化（不换骨架）

如果不想换算法但想压一压性能：

- **扩展 quot 快路径**：当前 `quot=1/2/3` 短路；可以加到 `quot ∈ [1, 16]` 用查表 / 分支预测来命中更高商分布（注意 u128 下高商 case 比 u64 更多）。
- **u128 单独走 Lehmer kernel**：在 u128 xgcd 内部用 u64 预测两步商，再 fall back，避免 `__udivti3`。
- **`coeff_sub_mul` 内部用 widening mul**：当前 `factor.wrapping_mul(rhs)` 可能在 i64/i128 上没充分利用 `mulx` —— 不一定有提升，需要 benchmark。

---

## 决策路标

| 触发条件 | 推荐方向 |
|---|---|
| FHE 热路径出现频繁的运行时 `gcdinv`（非编译期 / 非启动预计算） | 切到 **divsteps2 / Pornin** |
| 安全审计要求 modular inverse 常数时间（secret 模数或 secret 输入） | 切到 **divsteps2 / Pornin** |
| u128 `gcdinv` 成为 profile 热点 | 先试 **扩展 quot 快路径** 或 **Lehmer kernel**；不行再上 divsteps2 |
| `BigUint<S>` 需要 xgcd / 模逆 | **HGCD**（参考 GMP `mpn_hgcd`）|
| 只想多挤 10–20% 性能、不改架构 | 扩展 quot 快路径 + bench 验证 |

---

## 不必改的情况

- 当前 `xgcd`/`gcdinv` 主要用于 BarrettModulus / RNS 建造期 `gcdinv` 计算 modular inverse，这些是**一次性 / 编译期 / 构造期**调用，性能并非瓶颈。
- 没有 secret-dependent inversion 的需求时，constant-time 不是必需。
- 在确认热点 profile 之前，**不要预先换算法**——FLINT 风格已经够好。
