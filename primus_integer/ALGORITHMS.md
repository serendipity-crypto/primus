# 算法选型与未来改进方向

本文档说明 `primus_integer` 当前使用的算法、它们的取舍，以及如果将来需要优化时可以参考的更先进算法。不公开在 rustdoc 中，是给维护者的内部备忘。

---

## 当前实现

| 模块 / 方法 | 算法 | 关键路径下沉 | 来源 |
|---|---|---|---|
| `CarryingAdd` / `BorrowingSub` (u8–u64, usize) | 委托 `core::*::carrying_add`/`borrowing_sub` | x86_64 `adc`/`sbb`、aarch64 `adds`/`subs` | stdlib |
| `CarryingAdd` / `BorrowingSub` (u128) | 委托 stdlib（编译器拆成两条 `adc`） | 同上 | stdlib |
| `WideningMul` / `CarryingMul` (u8–u64, usize) | 升宽到 `2*T` 后单次 `mul` | x86_64 `mul`/`mulx` (BMI2)、aarch64 `umulh`+`mul` | stdlib + cast |
| `WideningMul` / `CarryingMul` (u128) | 手写 schoolbook：4× `wrapping_mul(u128, u128)` + carry chain | 编译器把每个 `wrapping_mul` 视为 64×64→128，下沉到 `mulx`/`umulh`，但中间 `add` 不一定形成 ADX 链 | 自写 |
| `WideningMul` / `CarryingMul`（SIMD u8/u16/u32） | `lhs.cast::<W>() * rhs.cast::<W>()` + 拆位 | nightly `core::simd`：AVX2/AVX-512 / NEON | stdlib portable_simd |
| `WideningMul` / `CarryingMul`（SIMD u64） | 手写 schoolbook：4× 32×32→64 `Simd` mul + carry chain | nightly `core::simd`；**已知 LLVM 对此 Rust 翻译 codegen 次优**（[mod.rs 注释](src/integer_traits/widening/widening_mul/simd.rs)）| LLVM `__mulddi3` 的 Rust 移植 |
| `DivWide` (u8–u64) | 升宽到 `2*T` 后单次 `div` | x86_64 `div` (RDX:RAX)、aarch64 `udiv` | stdlib |
| `DivWide` (u128) | 委托 `div_rem_scalar(&[lo, hi])` → Knuth D | 见下 | 自写 |
| `DivRemScalar` (u8–u64) | 每 limb `div_half` 或 `div_wide` 单步硬件除法 | x86_64 `div`、aarch64 `udiv` | 自写 |
| `DivRemScalar` (u128) | Knuth Algorithm D（normalize + 2 个 q-hat 估算 + 修正循环），高位 64 走 native u128/u64 div | x86_64 `div`、`mulx`；correction loop 是普通 sub | TAOCP Vol. 2 §4.3.1，[LLVM compiler-rt `udivti3`](https://github.com/llvm/llvm-project/blob/main/compiler-rt/lib/builtins/udivmodti4.c) |
| `BigUint::add_assign` / `sub_assign` | 串行 limb 循环，单条进/借位标志（`bool`） | LLVM 通常能合并 `adc`，但**不保证 ADCX/ADOX 并行链** | 自写 |
| `BigUint::mul_value_assign` | 串行 `carrying_mul` 链 | 同上 | 自写 |
| `BigUint::add_modulo_assign` / `sub_modulo_assign` | add + 条件 sub modulus（常规 modular add） | 常规 limb 循环 | 自写 |
| `multiply_many_values[_except[_inplace]]` | 串行 `mul_value_assign` 累乘 | 同 `mul_value_assign` | 自写 |

**这个组合的取舍**：

- 标量 `u8`–`u64` 已经下沉到硬件 `mul` / `div` / `adc` / `sbb`，几乎没有进一步空间。
- `u128` 系列**没有走硬件 128-bit 路径**——Rust stdlib 没有提供 `u128::widening_mul` 这种直达 `mulx`+`umulh` 链的形式，自写的 schoolbook 会被 LLVM 拆成 4 条 64×64 + 一堆 64-bit add，**没有形成 ADCX/ADOX 并行链**。是 u128 路径上最大的可量化损失。
- SIMD u64 widening mul 走"四块 32-bit 乘 + 拆位"，在 x86_64 上对应到 `vpmuludq` 但不利用 AVX-512 `vpmullq`/`vpmadd52luq`（HEXL 路径）；在 AArch64 上也没有走 NEON 的 `umullh`。
- 多 limb `BigUint` 操作走纯 Rust for 循环，**没有 SIMD、没有 ADCX/ADOX 双链**——属于"功能完整，性能尚未压榨"。
- Knuth D 是经典实现，但每位商有一次 q-hat 估算 + 最多两次修正循环；Möller–Granlund (2011) 的改良能把修正循环平均次数降到 < 0.5。

---

## 更优算法（按场景分）

### 1. u128 widening / carrying mul 走 `mulx` + `umulh` intrinsic 链

**适用**：`primus_modulus` 的 Barrett / Montgomery / Shoup 路径里出现 u128 mul-hi（高 128 位）；`primus_rns` / `primus_decompose` 的大模数运算；`u128` 作为底层 word 的 FHE 参数集。

**做法**：x86_64 用 [`core::arch::x86_64::_mulx_u64`](https://doc.rust-lang.org/core/arch/x86_64/fn._mulx_u64.html)（BMI2，需要 `target_feature = "bmi2"`）做 64×64→128 乘法、用 [`_addcarry_u64`](https://doc.rust-lang.org/core/arch/x86_64/fn._addcarry_u64.html) + [`_addcarryx_u64`](https://doc.rust-lang.org/core/arch/x86_64/fn._addcarryx_u64.html)（ADX）形成两条并行 carry chain；aarch64 直接用 `umulh` 内嵌汇编或 `core::arch::aarch64::__umulh`（如有）。

- **预期收益**：u128 `widening_mul` 从 ~14 µops（4× mulx + 5× add + masks）压到 ~6 µops（2× mulx + 1× mul + ADX 链），大约 2–2.5×。
- **参考实现**：
  - [`crypto-bigint::Limb::widening_mul`](https://github.com/RustCrypto/crypto-bigint/blob/master/src/uint/mul.rs)（多 limb 乘法 inner loop）。
  - GMP `mpn/x86_64/coreihwl/mul_1.asm`（ADX + MULX 经典编排）。
  - BoringSSL `crypto/bn/asm/x86_64-mont.pl`（Montgomery 内层）。
- **风险**：`target_feature = "bmi2"` 不开启时回落到 schoolbook；得给 `cfg(target_feature = "bmi2")` + runtime fallback 都写一份；不开启 `target-cpu=native` 时无效。

### 2. SIMD u64 widening mul 走 AVX-512 IFMA / HEXL 路径

**适用**：`primus_ntt` / `primus_poly` 的 element-wise modmul（关键 hot path）；批量同模数乘 reduce。

**做法**：
- **AVX-512 IFMA（`vpmadd52luq` / `vpmadd52huq`）**：把 52-bit ×52-bit→104-bit 乘加压成单条指令。需要把模数限制在 52-bit 内（很多 FHE 参数集已经这样）。
- **AVX-512 普通路径（`vpmullq` + `vpmullq` 高半合成）**：直接 64×64→64 低，高半用 schoolbook 拆位。
- **AArch64 SVE / NEON**：NEON 没有原生 64×64→128，但 SVE2 有 `mul[ah]h`；当前最佳折中是 vpmuludq 风格的 32×32→64 拆位 + carry chain。

- **预期收益**：u64×u64 模乘在 AVX-512 IFMA 上比当前 schoolbook 快 3–4×（每 lane 1 条指令 vs 当前 ~10 条）。
- **参考实现**：
  - [Intel HEXL `EltwiseFMAModAVX512`](https://github.com/intel/hexl/blob/main/hexl/eltwise/eltwise-fma-mod-avx512.cpp)。
  - HEXL 论文：[Intel HEXL: Accelerating Homomorphic Encryption with Intel AVX512-IFMA52](https://arxiv.org/abs/2103.16400)。
  - `concrete-ntt` 已经在 `primus_ntt` 中用作 backend，可作对照。
- **风险**：IFMA 不是所有 x86_64 CPU 都有（Ice Lake-X / Tiger Lake 及更新）；CPUID 检测 + 多版本派发是标配；52-bit 模数限制需要上层 API 配合。

### 3. Knuth D 的 q-hat 改良（Möller–Granlund 2011）

**适用**：`udiv256_by_128_to_128` 已经是 hot path（u128 `DivWide` / `DivRemScalar`）；`primus_rns` 的多模数构造与 base conversion 走 u128 division 链。

**做法**：用预计算的归一化倒数 `v_recip = ⌊(B²−1)/v⌋ − B`（B = 2^64）替代 Knuth 的"高位试除 + 修正循环"。每位商只需 1 次 mul + 1 次 sub + ≤ 1 次修正，平均修正次数 < 0.5（Knuth D 是 ≤ 2）。

- **预期收益**：u128 单 limb 除法 1.3–1.7×（取决于 divisor 分布；高商分布下更明显）。
- **参考实现**：
  - 论文：Möller & Granlund. *Improved division by invariant integers*. IEEE TC 60(2), 2011.
  - [GMP `mpn/generic/sbpi1_div_q.c`](https://gmplib.org/repo/gmp/file/tip/mpn/generic/sbpi1_div_q.c)（`sbpi1_div_qr` 是这个算法的成熟实现）。
  - 也可参考 [`libdivide`](https://github.com/ridiculousfish/libdivide) 的 branchful 变种。
- **风险**：实现复杂度高于 Knuth D；correctness proof 不平凡；需要全位宽（u8/u16/u32/u64/u128）的独立倒数预计算。可优先只替换 u128 路径。

### 4. `BigUint::add_assign` / `sub_assign` / `mul_value_assign` 用 ADCX/ADOX 双链

**适用**：`primus_rns::compose` / `primus_lattice` 的多 limb 累加；`multiply_many_values` 在 `BigUint` 上的累乘内层循环。

**做法**：x86_64 上 Intel ADX 扩展（Broadwell+）允许两条并行的 carry chain：[`_addcarry_u64`](https://doc.rust-lang.org/core/arch/x86_64/fn._addcarry_u64.html)（用 CF）与 [`_addcarryx_u64`](https://doc.rust-lang.org/core/arch/x86_64/fn._addcarryx_u64.html)（用 OF）。两条链可并行执行，**multi-limb add/mul 内部 throughput 直接翻倍**。GMP、OpenSSL、BoringSSL、`crypto-bigint` 全在用。

- **预期收益**：BigUint add/sub 4-limb 以上 ~1.8×；`mul_value_assign` 内层 carry chain ~1.5×（因为 mul 本身比 add 慢，串行依赖较短）。
- **参考实现**：
  - GMP `mpn/x86_64/coreihwl/{add_n,sub_n,mul_1}.asm`。
  - BoringSSL `crypto/fipsmodule/bn/asm/x86_64-mont*.pl`。
  - [`crypto-bigint` `mul_inner`](https://github.com/RustCrypto/crypto-bigint/blob/master/src/uint/mul.rs)（用 `mulx_u64` + `_addcarryx_u64`）。
- **风险**：`adcx`/`adox` 在编译器内部表达需要 nightly 的 `core::arch::x86_64` intrinsic；ABI 限制使得手写 ASM 还更靠谱；不开启 `target_feature = "adx,bmi2"` 时无效。

### 5. `multiply_many_values` 走 Karatsuba（小幅）

**适用**：`primus_rns` 构造期算多模数乘积（一次性，但当模数 ≥ 16 时仍有效）；不在运行时 hot path。

**当前**：schoolbook 累乘，O(n²) limb 乘法（n = limb 数）。
**Karatsuba**：阈值 n ≥ 16–32 limb 时 O(n^1.585)。

- **预期收益**：n = 32 时约 1.5×，n = 128 时约 3×。但 `multiply_many_values` 输入通常 ≤ 16 limb，**收益微小，不必现在做**。
- **参考实现**：GMP `mpn/generic/{mul_basecase,toom22_mul,toom33_mul,mul_fft}.c`。
- 远期若 `BigUint × BigUint` 进入 API（不只是 mul-by-scalar），再考虑。

### 6. 小幅优化（不换骨架）

- **`mul_value_assign(1)` 早返回**：当前只对 `value == 0` 短路；`value == 1` 走全 loop（每 limb 一次 `carrying_mul(0, ele, 0) = (ele, 0)`，浪费）。加 1 行 `if value.is_one() { return T::ZERO; }`。
- **`BigUint::bits_count` 批量 OR**：当前 `iter().rev().find()` 顺序找最高非零 limb；超大 limb 数（≥ 64）下可一次 OR 4 个 u64 chunk + 末尾再细化。基本不在 hot path。
- **`BigUint::PartialEq::eq` 用 `slice::eq`**：当前 `iter().zip().all`；改用 `self.0.as_slice() == other.0.as_slice()` 让 stdlib SIMD memcmp 接管。
- **`udiv256_by_128_to_128` 把 `un128 / vn1_u64` 的 `vn1_u64 as u128` 改成保留 u64**：当前 `un128 / (vn1_u64 as u128)` 让 LLVM 偶尔回落 u128 路径；用 `(un128 as u64).into() / vn1_u64` 手动 hint hardware `div` 可能稳一点。Needs codegen check.

---

## 决策路标

| 触发条件 | 推荐方向 |
|---|---|
| u128 mul 进入 profile 热点（典型：FHE 大模数 Barrett / Montgomery 路径） | **#1 mulx + ADX 链** |
| NTT / element-wise modmul 是瓶颈 | **#2 AVX-512 IFMA / HEXL** |
| u128 division 在 RNS / Barrett 预计算外的运行时也频繁 | **#3 Möller–Granlund 改良** |
| `BigUint` 多 limb add/sub/mul 在 RNS / lattice 上 profile 热 | **#4 ADCX/ADOX 双链** |
| `multiply_many_values` 输入 ≥ 16 limb 且构造期成瓶颈 | **#5 Karatsuba**（很少触发）|
| 只想多挤 5–10% 性能、不改架构 | **#6 mul_value(1) 早返回 + slice eq + codegen 微调** |
| stable 用户也想要 SIMD | 写 stable `core::arch::x86_64` + `core::arch::aarch64` 路径（不依赖 `portable_simd`），通过 `cfg(target_feature)` 选择 |

---

## 不必改的情况

- 标量 `u8`–`u64` 的 `widening_mul` / `carrying_mul` / `div_wide` / `carrying_add` 已经全部走硬件 `mul`/`mulx`/`div`/`adc`，**不要再"优化"**——LLVM 已经把它们 lower 到最优指令。手写 intrinsic 反而会限制 inline 与寄存器分配。
- `Bits` / `ByteCount` / `ConstBounded` / `Checked*` / `Overflowing*` / `Wrapping*` 全是 stdlib delegate，对应到 `popcnt`/`lzcnt`/`tzcnt`/`bswap` 等硬件指令，**没有进一步压榨空间**。
- u128 schoolbook 写法虽然不最优，但在**没有 BMI2/ADX 的目标**（如 SSE-only / 老 Atom / 部分 aarch64）上 LLVM 仍会下沉到合理指令，删掉它前要做 fallback。
- SIMD 路径目前 nightly-only：在 stable 化 `portable_simd` 之前，**用户量极小**，深度优化的 ROI 低。先等 stable 落地或评估手写 `core::arch` 路径的工程成本。
- 在确认 profile 热点之前，**不要预先换算法**——`primus_integer` 目前作为底层 trait crate，主要是给上层 crate 提供干净的抽象，性能瓶颈多半在更上层（NTT/RNS/FHE），先 profile 再优化。

---

## 文献

- Knuth. *The Art of Computer Programming, Vol. 2: Seminumerical Algorithms*. §4.3.1 Algorithm D（多精度除法）.
- Möller, Granlund. *Improved division by invariant integers*. IEEE Transactions on Computers 60(2), 2011. <https://gmplib.org/~tege/division-paper.pdf>
- Boemer et al. *Intel HEXL: Accelerating Homomorphic Encryption with Intel AVX512-IFMA52*. WAHC 2021. <https://arxiv.org/abs/2103.16400>
- Granlund, Möller. *Division by Invariant Integers using Multiplication*. PLDI 1994. <https://gmplib.org/~tege/divcnst-pldi94.pdf>
- GMP 源码：<https://gmplib.org/repo/gmp/> （`mpn/x86_64/coreihwl/*.asm` 与 `mpn/generic/*.c`）
- `crypto-bigint` 源码：<https://github.com/RustCrypto/crypto-bigint>
- Intel SDM Vol. 2A & 2B：`mulx`、`adcx`、`adox`、`vpmadd52luq`、`vpmullq` 指令说明。
