/// Zip an arbitrary number of iterators into one yielding flat tuples.
///
/// `izip!(a, b, c)` is equivalent to `a.zip(b).zip(c).map(|((x, y), z)| (x, y, z))`
/// but produces a flat `(x, y, z)` directly. Supports 1–16 input iterators.
/// Useful when iterating several big-integer limb buffers in lock-step.
///
/// # Source
///
/// Copied verbatim from [tfhe-fft](https://github.com/zama-ai/tfhe-rs/blob/4a73b7bb4b2a3e4209f8210c64521a96f0f0b0c1/tfhe-fft/src/lib.rs#L89).
/// The `@ __closure @` internal arm is structured this way to work around
/// rust-analyzer issue [#13526](https://github.com/rust-lang/rust-analyzer/issues/13526).
#[macro_export]
macro_rules! izip {
    // implemented this way to avoid a bug with type hints in rust-analyzer
    // https://github.com/rust-lang/rust-analyzer/issues/13526
    (@ __closure @ ($a:expr)) => { |a| (a,) };
    (@ __closure @ ($a:expr, $b:expr)) => { |(a, b)| (a, b) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr)) => { |((a, b), c)| (a, b, c) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr)) => { |(((a, b), c), d)| (a, b, c, d) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr)) => { |((((a, b), c), d), e)| (a, b, c, d, e) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr)) => { |(((((a, b), c), d), e), f)| (a, b, c, d, e, f) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr)) => { |((((((a, b), c), d), e), f), g)| (a, b, c, d, e, f, g) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr)) => { |(((((((a, b), c), d), e), f), g), h)| (a, b, c, d, e, f, g, h) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr)) => { |((((((((a, b), c), d), e), f), g), h), i)| (a, b, c, d, e, f, g, h, i) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr, $j: expr)) => { |(((((((((a, b), c), d), e), f), g), h), i), j)| (a, b, c, d, e, f, g, h, i, j) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr, $j: expr, $k: expr)) => { |((((((((((a, b), c), d), e), f), g), h), i), j), k)| (a, b, c, d, e, f, g, h, i, j, k) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr, $j: expr, $k: expr, $l: expr)) => { |(((((((((((a, b), c), d), e), f), g), h), i), j), k), l)| (a, b, c, d, e, f, g, h, i, j, k, l) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr, $j: expr, $k: expr, $l: expr, $m:expr)) => { |((((((((((((a, b), c), d), e), f), g), h), i), j), k), l), m)| (a, b, c, d, e, f, g, h, i, j, k, l, m) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr, $j: expr, $k: expr, $l: expr, $m:expr, $n:expr)) => { |(((((((((((((a, b), c), d), e), f), g), h), i), j), k), l), m), n)| (a, b, c, d, e, f, g, h, i, j, k, l, m, n) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr, $j: expr, $k: expr, $l: expr, $m:expr, $n:expr, $o:expr)) => { |((((((((((((((a, b), c), d), e), f), g), h), i), j), k), l), m), n), o)| (a, b, c, d, e, f, g, h, i, j, k, l, m, n, o) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr, $j: expr, $k: expr, $l: expr, $m:expr, $n:expr, $o:expr, $p: expr)) => { |(((((((((((((((a, b), c), d), e), f), g), h), i), j), k), l), m), n), o), p)| (a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) };

    ( $first:expr $(,)?) => {
        {
            ::core::iter::IntoIterator::into_iter($first)
        }
    };
    ( $first:expr, $($rest:expr),+ $(,)?) => {
        {
            ::core::iter::IntoIterator::into_iter($first)
                $(.zip($rest))*
                .map(izip!(@ __closure @ ($first, $($rest),*)))
        }
    };
}

/// Generate `<Poly>Iter` / `<Poly>IterMut` chunk-iterators for a wrapper type.
///
/// Given a wrapper `$poly` of the form `pub struct $poly<S>(S)` over slice
/// storage, expands to two iterator types — `<$poly>Iter` and
/// `<$poly>IterMut` — that walk the underlying buffer in fixed-size chunks of
/// `[<$short_name _len>]` elements and yield `$poly<&[T]>` / `$poly<&mut [T]>`
/// items. Used to give big-integer/polynomial wrappers row-iteration without
/// duplicating boilerplate per wrapper type.
///
/// # Arguments
///
/// * `$poly` — the wrapper type name (e.g. `BigUint`).
/// * `$short_name` — short label used to derive the chunk-length parameter
///   name (e.g. `bit_uint` → `bit_uint_len`).
#[macro_export]
macro_rules! impl_iters {
    ($poly:ident, $short_name:ident) => {
        paste::paste! {
            pub struct [<$poly Iter>]<'a, T>
            where
                T: UnsignedInteger,
            {
                pub iter: core::slice::ChunksExact<'a, T>
            }

            impl<'a, T: UnsignedInteger> [<$poly Iter>]<'a, T> {
                #[inline]
                pub fn new(data:&'a [T], [<$short_name _len>]:usize) -> Self{
                    Self {
                        iter: data.chunks_exact([<$short_name _len>])
                    }
                }
            }

            impl<'a, T: UnsignedInteger> Iterator for [<$poly Iter>]<'a, T> {
                type Item = $poly<&'a [T]>;

                #[inline]
                fn next(&mut self) -> Option<Self::Item> {
                    self.iter.next().map(|slice| $poly(slice))
                }
            }

            impl<'a, T: UnsignedInteger> core::iter::FusedIterator for [<$poly Iter>]<'a, T> {}
        }

        paste::paste! {
            pub struct [<$poly IterMut>]<'a, T>
            where
                T: UnsignedInteger,
            {
                pub iter: core::slice::ChunksExactMut<'a, T>
            }

            impl<'a, T: UnsignedInteger> [<$poly IterMut>]<'a, T> {
                #[inline]
                pub fn new(data:&'a mut [T], [<$short_name _len>]:usize) -> Self{
                    Self {
                        iter: data.chunks_exact_mut([<$short_name _len>])
                    }
                }
            }

            impl<'a, T: UnsignedInteger> Iterator for [<$poly IterMut>]<'a, T> {
                type Item = $poly<&'a mut [T]>;

                #[inline]
                fn next(&mut self) -> Option<Self::Item> {
                    self.iter.next().map(|slice| $poly(slice))
                }
            }

            impl<'a, T: UnsignedInteger> core::iter::FusedIterator for [<$poly IterMut>]<'a, T> {}
        }
    };
}
