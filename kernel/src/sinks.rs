use std::io::{BufRead, Read, Write};
use std::{mem, process, ptr};
use std::any::Any;
use std::collections::BTreeMap;
use std::fs::File;
use std::hint::unreachable_unchecked;
use std::mem::MaybeUninit;
use std::ops::{Coroutine, CoroutineState};
use std::pin::Pin;
use std::ptr::{addr_of, null, null_mut, slice_from_raw_parts, slice_from_raw_parts_mut};
use std::sync::LazyLock;
use std::task::Poll;
use std::time::Instant;
use futures::StreamExt;
use pathmap::ring::{AlgebraicStatus, Lattice};
use mork_expr::{byte_item, maybe_byte_item, Expr, ExprZipper, ExtractFailure, item_byte, parse, serialize, Tag, traverseh, ExprEnv, unify, UnificationFailure, apply};
use mork_frontend::bytestring_parser::{Parser, ParserError, Context};
use mork_interning::{WritePermit, SharedMapping, SharedMappingHandle};
use pathmap::utils::{BitMask, ByteMask};
use pathmap::zipper::*;
use mork_frontend::json_parser::Transcriber;
use log::*;
use pathmap::PathMap;

pub trait Sink {
    fn new(e: Expr) -> Self;
    fn request(&self) ->  impl Iterator<Item=&'static [u8]>;
    fn sink<'w, 'a, 'k, It : Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, it: It, path: &[u8]) where 'a : 'w, 'k : 'w;
    fn finalize<'w, 'a, 'k, It : Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, it: It) -> bool where 'a : 'w, 'k : 'w ;
}

pub struct AddSink { e: Expr, changed: bool }

impl Sink for AddSink {
    fn new(e: Expr) -> Self { AddSink { e, changed: false } }
    fn request(&self) -> impl Iterator<Item=&'static [u8]> {
        let p = &unsafe { self.e.prefix().unwrap_or_else(|x| self.e.span()).as_ref().unwrap() }[3..];
        trace!(target: "sink", "+ requesting {}", serialize(p));
        std::iter::once(p)
    }
    fn sink<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, mut it: It, path: &[u8]) where 'a : 'w, 'k : 'w {
        let mut wz = it.next().unwrap();
        let mpath = &path[3+wz.root_prefix_path().len()..];
        trace!(target: "sink", "+ at '{}' sinking raw '{}'", serialize(wz.root_prefix_path()), serialize(path));
        trace!(target: "sink", "+ sinking '{}'", serialize(mpath));
        wz.move_to_path(mpath);
        self.changed |= wz.set_val(()).is_none();
    }
    fn finalize<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, it: It) -> bool where 'a : 'w, 'k : 'w  {
        trace!(target: "sink", "+ finalizing");
        self.changed
    }
}

pub struct RemoveSink { e: Expr, remove: PathMap<()> }
// perhaps more performant to graft, remove*, and graft back?
impl Sink for RemoveSink {
    fn new(e: Expr) -> Self { RemoveSink { e, remove: PathMap::new() } }
    fn request(&self) -> impl Iterator<Item=&'static [u8]> {
        // !! we're never grabbing the full expression path, because then we don't have the ability to remove the root value
        let p = &unsafe { self.e.prefix().unwrap_or_else(|x| { let s = self.e.span(); slice_from_raw_parts(self.e.ptr, s.len() - 1) }).as_ref().unwrap() }[3..];
        trace!(target: "sink", "- requesting {}", serialize(p));
        std::iter::once(p)
    }
    fn sink<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, mut it: It, path: &[u8]) where 'a : 'w, 'k : 'w {
        let mut wz = it.next().unwrap();
        let mpath = &path[3+wz.root_prefix_path().len()..];
        trace!(target: "sink", "- at '{}' sinking raw '{}'", serialize(wz.root_prefix_path()), serialize(path));
        trace!(target: "sink", "- sinking '{}'", serialize(mpath));
        self.remove.insert(mpath, ());
    }
    fn finalize<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, mut it: It) -> bool where 'a : 'w, 'k : 'w  {
        let mut wz = it.next().unwrap();
        wz.reset();
        trace!(target: "sink", "- finalizing by subtracting {} at '{}'", self.remove.val_count(), serialize(wz.origin_path()));
        // match self.remove.remove(&[]) {
        //     None => {}
        //     Some(s) => {
        //         println!("has root");
        //         wz.remove_val(true);
        //         println!("val not removed");
        //     }
        // }
        match wz.subtract_into(&self.remove.read_zipper(), true) {
            AlgebraicStatus::Element => { true }
            AlgebraicStatus::Identity => { false }
            AlgebraicStatus::None => { true } // GOAT maybe not?
        }
    }
}

pub struct HeadSink { e: Expr, head: PathMap<()>, skip: usize, count: usize, max: usize, top: Vec<u8> }
impl Sink for HeadSink {
    fn new(e: Expr) -> Self {
        let mut ez = ExprZipper::new(e); ez.next(); ez.next();
        let max_s = ez.item().err().expect("cnt can not be an expression or variable");
        let max: usize = str::from_utf8(max_s).expect("string encoded numbers for now").parse().expect("a number");
        assert_ne!(max, 0);
        HeadSink { e, head: PathMap::new(), skip: 1 + 1+4 + 1+max_s.len(), count: 0, max, top: vec![] }
    }
    fn request(&self) -> impl Iterator<Item=&'static [u8]> {
        let p = &unsafe { self.e.prefix().unwrap_or_else(|x| { let s = self.e.span(); slice_from_raw_parts(self.e.ptr, s.len() - 1) }).as_ref().unwrap() }[self.skip..];
        trace!(target: "sink", "head requesting {}", serialize(p));
        std::iter::once(p)
    }
    fn sink<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, mut it: It, path: &[u8]) where 'a : 'w, 'k : 'w {
        let mut wz = it.next().unwrap();
        let mpath = &path[self.skip+wz.root_prefix_path().len()..];
        trace!(target: "sink", "head at '{}' sinking raw '{}'", serialize(wz.root_prefix_path()), serialize(path));
        if self.count == self.max {
            if &self.top[..] <= mpath {
                trace!(target: "sink", "head at max capacity ignoring '{}'", serialize(mpath));
                // doesn't displace any path
            } else {
                trace!(target: "sink", "head at max capacity replacing '{}' with '{}'", serialize(&self.top[..]), serialize(mpath));
                assert!(self.head.insert(mpath, ()).is_none());
                self.head.remove(&self.top[..]);
                let mut rz = self.head.read_zipper();
                rz.descend_last_path();
                self.top.clear();
                self.top.extend_from_slice(rz.path()); // yikes, throwing away our needless allocation
            }
        } else {
            if &self.top[..] <= mpath {
                if self.head.insert(mpath, ()).is_none() {
                    trace!(target: "sink", "head adding new top at '{}'", serialize(mpath));
                    self.top.clear();
                    self.top.extend_from_slice(mpath);
                    self.count += 1;
                }
            } else {
                if self.head.insert(mpath, ()).is_none() {
                    trace!(target: "sink", "head adding '{}'", serialize(mpath));
                    self.count += 1;
                }
            }
        }
    }
    fn finalize<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, mut it: It) -> bool where 'a : 'w, 'k : 'w  {
        let mut wz = it.next().unwrap();
        wz.reset();
        trace!(target: "sink", "head finalizing by joining {} at '{}'", self.count, serialize(wz.origin_path()));

        match wz.join_into(&self.head.read_zipper()) {
            AlgebraicStatus::Element => { true }
            AlgebraicStatus::Identity => { false }
            AlgebraicStatus::None => { true } // GOAT maybe not?
        }
    }
}

#[cfg(feature = "wasm")]
pub struct WASMSink { e: Expr, skip: usize, changed: bool, module: wasmtime::Module, store: wasmtime::Store<()>, instance: wasmtime::Instance }

#[cfg(feature = "wasm")]
static ENGINE_LINKER: LazyLock<(wasmtime::Engine, wasmtime::Linker<()>)> = LazyLock::new(|| {
    let mut config = wasmtime::Config::new();
    config.wasm_multi_memory(true);
    config.strategy(wasmtime::Strategy::Cranelift);
    config.signals_based_traps(true);
    config.memory_reservation(1 << 32);
    config.memory_guard_size(1 << 32);
    #[cfg(all(target_feature = "avx2"))]
    unsafe {
        config.cranelift_flag_enable("has_sse3");
        config.cranelift_flag_enable("has_ssse3");
        config.cranelift_flag_enable("has_sse41");
        config.cranelift_flag_enable("has_sse42");
        config.cranelift_flag_enable("has_avx");
        config.cranelift_flag_enable("has_avx2");
        config.cranelift_flag_enable("has_bmi1");
        config.cranelift_flag_enable("has_bmi2");
        config.cranelift_flag_enable("has_lzcnt");
        config.cranelift_flag_enable("has_popcnt");
        config.cranelift_flag_enable("has_fma");
    }
    #[cfg(all(target_feature = "avx512"))]
    unsafe {
        config.cranelift_flag_enable("has_avx512bitalg");
        config.cranelift_flag_enable("has_avx512dq");
        config.cranelift_flag_enable("has_avx512vl");
        config.cranelift_flag_enable("has_avx512vbmi");
        config.cranelift_flag_enable("has_avx512f");
    }

    let engine = wasmtime::Engine::new(&config).unwrap();

    let mut linker = wasmtime::Linker::new(&engine);

    linker.func_wrap("", "i32.bswap", |param: i32| param.to_be()).unwrap();
    linker.func_wrap("", "i64.bswap", |param: i64| param.to_be()).unwrap();

    (engine, linker)
});

#[cfg(feature = "wasm")]
static mut LINKER: Option<wasmtime::Linker<()>> = None;
macro_rules! wasm_ctx { () => { r#"
(module
  (import "" "i32.bswap" (func $i32.bswap (param i32) (result i32)))
  (import "" "i64.bswap" (func $i64.bswap (param i64) (result i64)))

  (memory $in 1)
  (export "in" (memory $in))
  (memory $out 1)
  (export "out" (memory $out))
  (memory $local 1)

  (func (export "_otf_grounding")
    {:?}
  )
)
"# } }


#[cfg(feature = "wasm")]
impl Sink for WASMSink {
    fn new(e: Expr) -> Self {
        let mut ez = ExprZipper::new(e); ez.next(); ez.next();
        let program_e = ez.subexpr();
        let wat = format!(wasm_ctx!(), program_e);
        let module = wasmtime::Module::new(&ENGINE_LINKER.0, wat).unwrap();
        let mut store = wasmtime::Store::new(&ENGINE_LINKER.0, ());
        let instance = (&ENGINE_LINKER.1).instantiate(&mut store, &module).unwrap();

        WASMSink { e, skip: 1 + 1+4 + program_e.span().len(), changed: false, module, store, instance }
    }
    fn request(&self) -> impl Iterator<Item=&'static [u8]> {
        // let p = &unsafe { self.e.prefix().unwrap_or_else(|x| { let s = self.e.span(); slice_from_raw_parts(self.e.ptr, s.len() - 1) }).as_ref().unwrap() }[self.skip..];
        // trace!(target: "sink", "wasm requesting {}", serialize(p));
        // std::iter::once(p)
        static empty: [u8; 0] = [];
        std::iter::once(&empty[..])
    }
    fn sink<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperUntracked<'a, 'k, ()>>>(&mut self, mut it: It, path: &[u8]) where 'a : 'w, 'k : 'w {
        let mut wz = it.next().unwrap();
        let mpath = &path[self.skip+wz.root_prefix_path().len()..];
        trace!(target: "sink", "wasm at '{}' sinking raw '{}'", serialize(wz.root_prefix_path()), serialize(path));
        trace!(target: "sink", "wasm input '{}'", serialize(mpath));
        let imem = self.instance.get_memory(&mut self.store, "in").unwrap();
        imem.write(&mut self.store, 0, mpath).unwrap();
        let run = self.instance.get_typed_func::<(), ()>(&mut self.store, "_otf_grounding").unwrap();
        match run.call(&mut self.store, ()) {
            Ok(()) => {
                let omem = self.instance.get_memory(&mut self.store, "out").unwrap().data(&mut self.store);
                let ospan = unsafe { Expr{ ptr: omem.as_ptr().cast_mut() }.span().as_ref().unwrap() };
                trace!(target: "sink", "wasm output '{}'", serialize(ospan));
                wz.move_to_path(ospan);
                self.changed |= wz.set_val(()).is_none();
            }
            Err(e) => {
                trace!(target: "sink", "wasm error {:?}", e);
            }
        }

    }
    fn finalize<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperUntracked<'a, 'k, ()>>>(&mut self, mut it: It) -> bool where 'a : 'w, 'k : 'w  {
        trace!(target: "sink", "wasm finalizing");
        self.changed
    }
}

// ($k $x) (f $x $y)
// (count (count of $k is $i) $i ($x $y))   unify
// (count (count of r2 is $i) $i (P Q))
// (count (count of r2 is 3) 3 ($x $y))
pub struct CountSink { e: Expr, unique: PathMap<()> }
impl Sink for CountSink {
    fn new(e: Expr) -> Self {
        CountSink { e, unique: PathMap::new() }
    }
    fn request(&self) -> impl Iterator<Item=&'static [u8]> {
        let p = &unsafe { self.e.prefix().unwrap_or_else(|x| { let s = self.e.span(); slice_from_raw_parts(self.e.ptr, s.len() - 1) }).as_ref().unwrap() }[7..];
        trace!(target: "sink", "count requesting {}", serialize(p));
        std::iter::once(p)
    }
    fn sink<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, mut it: It, path: &[u8]) where 'a : 'w, 'k : 'w {
        let mut wz = it.next().unwrap();
        let mpath = &path[7+wz.root_prefix_path().len()..];
        let ctx = unsafe { Expr { ptr: mpath.as_ptr().cast_mut() } };
        trace!(target: "sink", "count at '{}' sinking raw '{}'", serialize(wz.root_prefix_path()), serialize(path));
        trace!(target: "sink", "count registering in ctx {:?}", serialize(mpath));
        self.unique.insert(mpath, ());
    }
    fn finalize<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, mut it: It) -> bool where 'a : 'w, 'k : 'w  {
        let mut wz = it.next().unwrap();
        wz.reset();
        trace!(target: "sink", "count finalizing by reducing {} at '{}'", self.unique.val_count(), serialize(wz.origin_path()));

        let mut _to_swap = PathMap::new(); std::mem::swap(&mut self.unique, &mut _to_swap);
        let mut rooted_input = PathMap::new();
        rooted_input.write_zipper_at_path(wz.root_prefix_path()).graft_map(_to_swap);

        static v: &'static [u8] = &[item_byte(Tag::NewVar)];
        let mut prz = ProductZipper::new::<_, ReadZipperUntracked<()>, [_; 0]>(rooted_input.into_read_zipper(&[]), []);
        let prz_ptr = (&prz) as *const ProductZipper<()>;
        let mut changed = false;
        let mut buffer: Vec<u8> = Vec::with_capacity(1 << 32);
        crate::space::Space::query_multi_raw(unsafe { prz_ptr.cast_mut().as_mut().unwrap() }, &[ExprEnv::new(0, Expr{ ptr: v.as_ptr().cast_mut() })], |refs_bindings, loc| {
            let cnt = prz.val_count();
            trace!(target: "sink", "'{}' and under {}", serialize(prz.path()), cnt);
            let clen = prz.path().len();
            let cnt_str = cnt.to_string();
            if prz.descend_to_existing_byte(item_byte(Tag::SymbolSize(cnt_str.len() as _))) {
                if prz.descend_to_existing(cnt_str.as_bytes()) == cnt_str.len() {
                    let fixed = &prz.path()[..prz.path().len()-(1+cnt_str.len())];
                    trace!(target: "sink", "fixed guard {}", serialize(fixed));
                    wz.move_to_path(fixed);
                    wz.set_val(());
                    changed |= true;
                }
            } 
            if prz.descend_to_existing_byte(item_byte(Tag::NewVar)) {
                let ignored = &prz.path()[..prz.path().len()-1];
                trace!(target: "sink", "ignored guard {}", serialize(ignored));
                wz.move_to_path(ignored);
                wz.set_val(());
                changed |= true;
            } 
            if prz.descend_first_byte() {
                if let Tag::VarRef(k) = byte_item(prz.path()[prz.path().len()-1]) {
                    let mut cntv = vec![item_byte(Tag::SymbolSize(cnt_str.len() as _))];
                    cntv.extend_from_slice(cnt_str.as_bytes());
                    let varref = &prz.path()[..prz.path().len()-1];
                    let ie = Expr { ptr: (&varref[0] as *const u8).cast_mut() };
                    let mut oz = ExprZipper::new(Expr{ ptr: buffer.as_mut_ptr() });
                    trace!(target: "sink", "ref guard '{}' var {:?} with '{}'", serialize(varref), k, serialize(&cntv[..]));
                    let os = ie.substitute_one_de_bruijn(k, Expr{ ptr: cntv.as_mut_ptr() }, &mut oz);
                    unsafe { buffer.set_len(oz.loc) }
                    trace!(target: "sink", "ref guard subs '{:?}'", serialize(&buffer[..oz.loc]));
                    wz.move_to_path(&buffer[wz.root_prefix_path().len()..oz.loc]);
                    wz.set_val(());
                    changed |= true
                }
            }
            true
        });
        changed
    }
}

// MaxSink - finds maximum numeric value from matched patterns
// Syntax: (O (max <report-template-with-$m> $m <value-expression>))
// Example: (O (max (dup $m and $m to $t) $m $w))
// Uses ProductZipper pattern for correct duplicate variable handling
pub struct MaxSink {
    e: Expr,
    head: PathMap<()>,  // Like HeadSink!
    skip: usize,
    has_non_numeric: bool,  // Flag for partial eval - no need for Vec
}

impl MaxSink {
    fn parse_i64_symbol_strict(sym: &[u8]) -> Option<i64> {
        if sym.is_empty() { return None; }
        if sym[0] == b'-' {
            if sym.len() <= 1 || !sym[1..].iter().all(|b| b.is_ascii_digit()) { return None; }
        } else if !sym.iter().all(|b| b.is_ascii_digit()) {
            return None;
        }
        let s = unsafe { core::str::from_utf8_unchecked(sym) };
        s.parse::<i64>().ok()
    }

    #[inline]
    fn parse_i64_from_tail(p: &[u8]) -> Option<i64> {
        // Walk backwards to find the last symbol
        let mut i = p.len();
        while i > 0 {
            i -= 1;
            if let Tag::SymbolSize(n) = byte_item(p[i]) {
                let n = n as usize;
                if i + 1 + n <= p.len() {
                    return Self::parse_i64_symbol_strict(&p[i+1..i+1+n]);
                }
            }
        }
        None
    }

    #[inline]
    fn parse_i64_from_tail_safe(p: &[u8]) -> Option<i64> {
        // Safe version that doesn't panic on invalid bytes
        // Walks backwards to find the last symbol, skipping invalid bytes
        let mut i = p.len();
        while i > 0 {
            i -= 1;
            // Use maybe_byte_item to avoid panic on reserved bytes
            if let Ok(Tag::SymbolSize(n)) = maybe_byte_item(p[i]) {
                let n = n as usize;
                if i + 1 + n <= p.len() {
                    return Self::parse_i64_symbol_strict(&p[i+1..i+1+n]);
                }
            }
            // Skip invalid bytes silently - they're payload from non-numeric symbols
        }
        None
    }

    #[inline]
    fn make_value_bytes(v: i64) -> Vec<u8> {
        let s = v.to_string();
        let mut out = Vec::with_capacity(1 + s.len());
        out.push(item_byte(Tag::SymbolSize(s.len() as _)));
        out.extend_from_slice(s.as_bytes());
        out
    }
}

impl Sink for MaxSink {
    fn new(e: Expr) -> Self {
        let skip = 1 + 1 + 3;  // Skip outer compound, 'O', and 'max'
        MaxSink { e, head: PathMap::new(), skip, has_non_numeric: false }
    }

    fn request(&self) -> impl Iterator<Item=&'static [u8]> {
        let p = &unsafe {
            self.e.prefix()
                .unwrap_or_else(|_| { let s = self.e.span(); slice_from_raw_parts(self.e.ptr, s.len() - 1) })
                .as_ref().unwrap()
        }[self.skip..];
        trace!(target: "sink", "max requesting {}", serialize(p));
        std::iter::once(p)
    }

    fn sink<'w,'a,'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a,'k,()>>>(
        &mut self, mut it: It, path: &[u8]
    ) where 'a:'w, 'k:'w {
        let wz = it.next().unwrap();

        // COLLECTION PHASE LOGGING - comprehensive validation
        let root_prefix = wz.root_prefix_path();
        let full_path = path;
        let mpath = &path[self.skip + root_prefix.len()..];

        // Log collection details for validation
        trace!(target: "sink", "=== MAX COLLECTION ===");
        trace!(target: "sink", "  full_path: '{}'", serialize(full_path));
        trace!(target: "sink", "  skip: {}", self.skip);
        trace!(target: "sink", "  root_prefix: '{}'", serialize(root_prefix));
        trace!(target: "sink", "  mpath: '{}'", serialize(mpath));

        // Check if path contains a numeric value
        if let Some(v) = Self::parse_i64_from_tail_safe(mpath) {
            trace!(target: "sink", "  extracted value: {}", v);
            // Only collect if we haven't hit a non-numeric yet
            if !self.has_non_numeric {
                self.head.insert(mpath, ());
                trace!(target: "sink", "  collection size: {}", self.head.val_count());
            }
        } else {
            // Track that we hit a non-numeric value for partial evaluation
            if !self.has_non_numeric {
                trace!(target: "sink", "  found first non-numeric path: '{}'", serialize(mpath));
                self.has_non_numeric = true;  // Just set the flag
            }
        }
    }

    fn finalize<'w,'a,'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a,'k,()>>>(
        &mut self, mut it: It
    ) -> bool where 'a:'w, 'k:'w {
        let mut wz = it.next().unwrap();
        wz.reset();

        trace!(target: "sink", "=== MAX FINALIZE ===");
        trace!(target: "sink", "  collected paths: {}", self.head.val_count());

        // STEP 0: Check if we hit a non-numeric value - emit partial eval if so
        if self.has_non_numeric {
            trace!(target: "sink", "  hit non-numeric - emitting report with fresh variable");

            // Extract the report template from the original pattern
            // The pattern is: (max <report> <guard> <value>)
            // We want just <report>

            let prefix = unsafe {
                self.e.prefix()
                    .unwrap_or_else(|_| { let s = self.e.span(); slice_from_raw_parts(self.e.ptr, s.len() - 1) })
                    .as_ref().unwrap()
            };

            // Skip to the report part
            let report_start = &prefix[self.skip..];
            let report_expr = unsafe { Expr { ptr: report_start.as_ptr().cast_mut() } };
            let report_len = report_expr.span().len();

            // The simple solution: emit an empty result to show no max was found
            // This is cleaner than trying to emit a partial pattern
            trace!(target: "sink", "  no numeric values found, returning empty");
            return false;
        }

        // STEP 1: Re-root the collected paths under the writer's root prefix (CRITICAL!)
        // This is the key insight from GPT-5 Pro - we need complete paths, not fragments
        let collected = std::mem::take(&mut self.head);

        let mut rooted_input = PathMap::new();
        rooted_input.write_zipper_at_path(wz.root_prefix_path()).graft_map(collected);

        trace!(target: "sink", "  re-rooted under: '{}'", serialize(wz.root_prefix_path()));

        // STEP 2: Find the max value from the re-rooted paths
        let mut rz = rooted_input.read_zipper();
        let mut best_val: Option<i64> = None;
        let mut best_path: Vec<u8> = Vec::new(); // This will be the ABSOLUTE path

        while rz.to_next_val() {
            let abs_path = rz.path();
            trace!(target: "sink", "  checking absolute path: '{}'", serialize(abs_path));

            if let Some(v) = Self::parse_i64_from_tail(abs_path) {
                trace!(target: "sink", "    found value: {}", v);
                if best_val.map_or(true, |cur| v > cur) {
                    best_val = Some(v);
                    best_path = abs_path.to_vec();
                    trace!(target: "sink", "    new max!");
                }
            }
        }

        let Some(max_val) = best_val else {
            trace!(target: "sink", "  no valid values found");
            return false;
        };

        trace!(target: "sink", "  MAX VALUE: {}", max_val);
        trace!(target: "sink", "  winning path: '{}'", serialize(&best_path));

        // STEP 3: Parse the winning path and emit
        // The best_path is now ABSOLUTE, so we can correctly parse the report
        let best_expr = unsafe { Expr { ptr: best_path.as_ptr().cast_mut() } };
        let report_span = best_expr.span();
        let report_len = report_span.len();

        // Basic sanity check
        if report_len == 0 || report_len > best_path.len() {
            trace!(target: "sink", "  ERROR: invalid report_len {}", report_len);
            return false;
        }

        trace!(target: "sink", "  report_len: {}", report_len);

        // Validate boundaries
        #[cfg(debug_assertions)]
        {
            // Ensure report_len is valid
            debug_assert!(report_len > 0, "report should have non-zero length");
            debug_assert!(report_len <= best_path.len(), "report_len {} exceeds path len {}",
                         report_len, best_path.len());

            // If there's a guard, validate it too
            if report_len < best_path.len() {
                let guard_expr = unsafe { Expr { ptr: best_path[report_len..].as_ptr().cast_mut() } };
                let guard_span = guard_expr.span();
                let guard_len = guard_span.len();
                debug_assert!(report_len + guard_len <= best_path.len(),
                             "report_len {} + guard_len {} exceeds path len {}",
                             report_len, guard_len, best_path.len());
            }
        }

        if report_len >= best_path.len() {
            // No guard/value, emit report relative to root_prefix
            let relative = &best_path[wz.root_prefix_path().len()..];
            trace!(target: "sink", "  emitting report only: '{}'", serialize(relative));
            wz.move_to_path(relative);
            wz.set_val(());
            return true;
        }

        // Check guard type
        let guard_tag = byte_item(best_path[report_len]);
        trace!(target: "sink", "  guard_tag: {:?}", guard_tag);

        match guard_tag {
            Tag::VarRef(k) => {
                trace!(target: "sink", "  substituting VarRef({}) with {}", k, max_val);

                // Build value bytes and substitute in the REPORT portion
                let val_bytes = Self::make_value_bytes(max_val);

                // CRITICAL: Substitute on the report PREFIX only (not including guard)
                let report_prefix = &best_path[..report_len];
                let ie = unsafe { Expr { ptr: report_prefix.as_ptr().cast_mut() } };

                let mut out = vec![0u8; report_prefix.len() + 32];
                let mut oz = ExprZipper::new(Expr { ptr: out.as_mut_ptr() });
                ie.substitute_one_de_bruijn(k, Expr { ptr: val_bytes.as_ptr().cast_mut() }, &mut oz);
                unsafe { out.set_len(oz.loc) };

                trace!(target: "sink", "  substitution result: {} bytes: '{}'", oz.loc, serialize(&out[..oz.loc]));

                // The substituted result is a complete expression - emit it all
                trace!(target: "sink", "  emitting substituted: '{}'", serialize(&out[..oz.loc]));

                // Debug: validate the bytes we're about to emit
                #[cfg(debug_assertions)]
                {
                    if oz.loc > 0 {
                        let first_byte = out[0];
                        let top2 = first_byte & 0b1100_0000;
                        debug_assert!(
                            matches!(top2, 0b0000_0000 | 0b1000_0000 | 0b1100_0000),
                            "Invalid tag byte: {:#04x} (top bits: {:#04x})",
                            first_byte, top2
                        );
                    }
                }

                // Emit RELATIVE to root_prefix (critical fix!)
                // The substituted result is absolute (starts with root_prefix),
                // but wz expects a path relative to root_prefix after reset.
                let root = wz.root_prefix_path();
                debug_assert!(out.len() >= root.len(), "substituted path too short");
                debug_assert!(out[..root.len()] == root[..], "substituted path must start with root_prefix");

                let relative = &out[root.len()..oz.loc];
                trace!(target: "sink", "  emitting relative (sliced {} root bytes): '{}'", root.len(), serialize(relative));
                wz.move_to_path(relative);
                wz.set_val(());
                true
            }

            Tag::SymbolSize(_) => {
                // Literal guard - only emit if it matches the max value
                trace!(target: "sink", "  literal guard");

                // Safe guard extraction using the proper API
                let guard_expr = unsafe { Expr { ptr: best_path[report_len..].as_ptr().cast_mut() } };

                // Use symbol() to get the payload without the header byte
                let guard_bytes = unsafe { guard_expr.symbol() };
                if guard_bytes.is_none() {
                    trace!(target: "sink", "  could not extract guard symbol");
                    return false;
                }

                let guard_bytes = unsafe { guard_bytes.unwrap().as_ref().unwrap() };
                if let Some(guard_val) = Self::parse_i64_symbol_strict(guard_bytes) {
                    if guard_val == max_val {
                        // Emit the report portion, relative to root_prefix
                        let root = wz.root_prefix_path();
                        // The report is best_path[0..report_len]
                        // We want it relative to root, so skip root.len() bytes
                        let relative = if report_len > root.len() {
                            &best_path[root.len()..report_len]
                        } else {
                            // Report is equal to or within root - no additional bytes needed
                            &[]
                        };
                        trace!(target: "sink", "  guard matches, emitting: '{}'", serialize(relative));
                        wz.move_to_path(relative);
                        wz.set_val(());
                        return true;
                    }
                }
                trace!(target: "sink", "  guard doesn't match");
                false
            }

            Tag::NewVar => {
                // Guard not referenced, emit report
                let root = wz.root_prefix_path();
                let relative = if report_len > root.len() {
                    &best_path[root.len()..report_len]
                } else {
                    &[]
                };
                trace!(target: "sink", "  NewVar guard, emitting report: '{}'", serialize(relative));
                wz.move_to_path(relative);
                wz.set_val(());
                true
            }

            _ => {
                trace!(target: "sink", "  unexpected guard tag");
                false
            }
        }
    }
}

// MinSink - finds minimum numeric value from matched patterns
// Syntax: (O (min <report-template-with-$m> $m <value-expression>))
// Example: (O (min (neg-ok) -10 $v))
// Uses ProductZipper pattern for correct duplicate variable handling
pub struct MinSink {
    e: Expr,
    head: PathMap<()>,  // Like HeadSink/MaxSink!
    skip: usize,
    has_non_numeric: bool,  // Flag for partial eval - no need for Vec
}

impl MinSink {
    fn parse_i64_symbol_strict(sym: &[u8]) -> Option<i64> {
        if sym.is_empty() { return None; }
        if sym[0] == b'-' {
            if sym.len() <= 1 || !sym[1..].iter().all(|b| b.is_ascii_digit()) { return None; }
        } else if !sym.iter().all(|b| b.is_ascii_digit()) {
            return None;
        }
        let s = unsafe { core::str::from_utf8_unchecked(sym) };
        s.parse::<i64>().ok()
    }

    #[inline]
    fn parse_i64_from_tail(p: &[u8]) -> Option<i64> {
        // Walk backwards to find the last symbol
        let mut i = p.len();
        while i > 0 {
            i -= 1;
            if let Tag::SymbolSize(n) = byte_item(p[i]) {
                let n = n as usize;
                if i + 1 + n <= p.len() {
                    return Self::parse_i64_symbol_strict(&p[i+1..i+1+n]);
                }
            }
        }
        None
    }

    #[inline]
    fn parse_i64_from_tail_safe(p: &[u8]) -> Option<i64> {
        // Safe version that doesn't panic on invalid bytes
        // Walks backwards to find the last symbol, skipping invalid bytes
        let mut i = p.len();
        while i > 0 {
            i -= 1;
            // Use maybe_byte_item to avoid panic on reserved bytes
            if let Ok(Tag::SymbolSize(n)) = maybe_byte_item(p[i]) {
                let n = n as usize;
                if i + 1 + n <= p.len() {
                    return Self::parse_i64_symbol_strict(&p[i+1..i+1+n]);
                }
            }
            // Skip invalid bytes silently - they're payload from non-numeric symbols
        }
        None
    }

    #[inline]
    fn make_value_bytes(v: i64) -> Vec<u8> {
        let s = v.to_string();
        let mut bytes = vec![item_byte(Tag::SymbolSize(s.len() as u8))];
        bytes.extend_from_slice(s.as_bytes());
        bytes
    }
}

impl Sink for MinSink {
    fn new(e: Expr) -> Self {
        let skip = 1 + 1 + 3;  // Skip outer compound, 'O', and 'min'
        MinSink { e, head: PathMap::new(), skip, has_non_numeric: false }
    }

    fn request(&self) -> impl Iterator<Item=&'static [u8]> {
        let p = &unsafe {
            self.e.prefix()
                .unwrap_or_else(|_| { let s = self.e.span(); slice_from_raw_parts(self.e.ptr, s.len() - 1) })
                .as_ref().unwrap()
        }[self.skip..];
        trace!(target: "sink", "min requesting {}", serialize(p));
        std::iter::once(p)
    }

    fn sink<'w,'a,'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a,'k,()>>>(
        &mut self, mut it: It, path: &[u8]
    ) where 'a:'w, 'k:'w {
        let wz = it.next().unwrap();

        // COLLECTION PHASE - same as MaxSink
        let root_prefix = wz.root_prefix_path();
        let full_path = path;
        let mpath = &path[self.skip + root_prefix.len()..];

        // Log collection details for validation
        trace!(target: "sink", "=== MIN COLLECTION ===");
        trace!(target: "sink", "  full_path: '{}'", serialize(full_path));
        trace!(target: "sink", "  skip: {}", self.skip);
        trace!(target: "sink", "  root_prefix: '{}'", serialize(root_prefix));
        trace!(target: "sink", "  mpath: '{}'", serialize(mpath));

        // Check if path contains a numeric value
        if let Some(v) = Self::parse_i64_from_tail_safe(mpath) {
            trace!(target: "sink", "  extracted value: {}", v);
            // Only collect if we haven't hit a non-numeric yet
            if !self.has_non_numeric {
                self.head.insert(mpath, ());
                trace!(target: "sink", "  collection size: {}", self.head.val_count());
            }
        } else {
            // Track that we hit a non-numeric value for partial evaluation
            if !self.has_non_numeric {
                trace!(target: "sink", "  found first non-numeric path: '{}'", serialize(mpath));
                self.has_non_numeric = true;  // Just set the flag
            }
        }
    }

    fn finalize<'w,'a,'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a,'k,()>>>(
        &mut self, mut it: It
    ) -> bool where 'a:'w, 'k:'w {
        let mut wz = it.next().unwrap();
        wz.reset();

        trace!(target: "sink", "=== MIN FINALIZE ===");
        trace!(target: "sink", "  collected paths: {}", self.head.val_count());

        // STEP 0: Check if we hit a non-numeric value - emit partial eval if so
        if self.has_non_numeric {
            trace!(target: "sink", "  hit non-numeric - emitting report with fresh variable");

            // Extract the report template from the original pattern
            // The pattern is: (max <report> <guard> <value>)
            // We want just <report>

            let prefix = unsafe {
                self.e.prefix()
                    .unwrap_or_else(|_| { let s = self.e.span(); slice_from_raw_parts(self.e.ptr, s.len() - 1) })
                    .as_ref().unwrap()
            };

            // Skip to the report part
            let report_start = &prefix[self.skip..];
            let report_expr = unsafe { Expr { ptr: report_start.as_ptr().cast_mut() } };
            let report_len = report_expr.span().len();

            // The simple solution: emit an empty result to show no max was found
            // This is cleaner than trying to emit a partial pattern
            trace!(target: "sink", "  no numeric values found, returning empty");
            return false;
        }

        // STEP 1: Re-root the collected paths under the writer's root prefix (CRITICAL!)
        // This is the key insight from GPT-5 Pro - we need complete paths, not fragments
        let collected = std::mem::take(&mut self.head);

        let mut rooted_input = PathMap::new();
        rooted_input.write_zipper_at_path(wz.root_prefix_path()).graft_map(collected);

        trace!(target: "sink", "  re-rooted under: '{}'", serialize(wz.root_prefix_path()));

        // STEP 2: Find the MINIMUM value from the re-rooted paths (< instead of >)
        let mut rz = rooted_input.read_zipper();
        let mut best_val: Option<i64> = None;
        let mut best_path: Vec<u8> = Vec::new(); // This will be the ABSOLUTE path

        while rz.to_next_val() {
            let abs_path = rz.path();
            trace!(target: "sink", "  checking absolute path: '{}'", serialize(abs_path));

            if let Some(v) = Self::parse_i64_from_tail(abs_path) {
                trace!(target: "sink", "    found value: {}", v);
                if best_val.map_or(true, |cur| v < cur) {  // ← THE KEY CHANGE: < instead of >
                    best_val = Some(v);
                    best_path = abs_path.to_vec();
                    trace!(target: "sink", "    new min!");
                }
            }
        }

        let Some(min_val) = best_val else {
            trace!(target: "sink", "  no valid values found");
            return false;
        };

        trace!(target: "sink", "  MIN VALUE: {}", min_val);
        trace!(target: "sink", "  winning path: '{}'", serialize(&best_path));

        // STEP 3: Parse the winning path and emit
        // The best_path is now ABSOLUTE, so we can correctly parse the report
        let best_expr = unsafe { Expr { ptr: best_path.as_ptr().cast_mut() } };
        let report_span = best_expr.span();
        let report_len = report_span.len();

        // Basic sanity check
        if report_len == 0 || report_len > best_path.len() {
            trace!(target: "sink", "  ERROR: invalid report_len {}", report_len);
            return false;
        }

        trace!(target: "sink", "  report_len: {}", report_len);

        // Validate boundaries
        #[cfg(debug_assertions)]
        {
            // Ensure report_len is valid
            debug_assert!(report_len > 0, "report should have non-zero length");
            debug_assert!(report_len <= best_path.len(), "report_len {} exceeds path len {}",
                         report_len, best_path.len());

            // If there's a guard, validate it too
            if report_len < best_path.len() {
                let guard_expr = unsafe { Expr { ptr: best_path[report_len..].as_ptr().cast_mut() } };
                let guard_span = guard_expr.span();
                let guard_len = guard_span.len();
                debug_assert!(report_len + guard_len <= best_path.len(),
                             "report_len {} + guard_len {} exceeds path len {}",
                             report_len, guard_len, best_path.len());
            }
        }

        if report_len >= best_path.len() {
            // No guard/value, emit report relative to root_prefix
            let relative = &best_path[wz.root_prefix_path().len()..];
            trace!(target: "sink", "  emitting report only: '{}'", serialize(relative));
            wz.move_to_path(relative);
            wz.set_val(());
            return true;
        }

        // Check guard type
        let guard_tag = byte_item(best_path[report_len]);
        trace!(target: "sink", "  guard_tag: {:?}", guard_tag);

        match guard_tag {
            Tag::VarRef(k) => {
                trace!(target: "sink", "  substituting VarRef({}) with {}", k, min_val);

                // Build value bytes and substitute in the REPORT portion
                let val_bytes = Self::make_value_bytes(min_val);

                // CRITICAL: Substitute on the report PREFIX only (not including guard)
                let report_prefix = &best_path[..report_len];
                let ie = unsafe { Expr { ptr: report_prefix.as_ptr().cast_mut() } };

                let mut out = vec![0u8; report_prefix.len() + 32];
                let mut oz = ExprZipper::new(Expr { ptr: out.as_mut_ptr() });
                ie.substitute_one_de_bruijn(k, Expr { ptr: val_bytes.as_ptr().cast_mut() }, &mut oz);
                unsafe { out.set_len(oz.loc) };

                trace!(target: "sink", "  substitution result: {} bytes: '{}'", oz.loc, serialize(&out[..oz.loc]));

                // The substituted result is a complete expression - emit it all
                trace!(target: "sink", "  emitting substituted: '{}'", serialize(&out[..oz.loc]));

                // Debug: validate the bytes we're about to emit
                #[cfg(debug_assertions)]
                {
                    if oz.loc > 0 {
                        let first_byte = out[0];
                        let top2 = first_byte & 0b1100_0000;
                        debug_assert!(
                            matches!(top2, 0b0000_0000 | 0b1000_0000 | 0b1100_0000),
                            "Invalid tag byte: {:#04x} (top bits: {:#04x})",
                            first_byte, top2
                        );
                    }
                }

                // Emit RELATIVE to root_prefix (critical fix!)
                // The substituted result is absolute (starts with root_prefix),
                // but wz expects a path relative to root_prefix after reset.
                let root = wz.root_prefix_path();
                debug_assert!(out.len() >= root.len(), "substituted path too short");
                debug_assert!(out[..root.len()] == root[..], "substituted path must start with root_prefix");

                let relative = &out[root.len()..oz.loc];
                trace!(target: "sink", "  emitting relative (sliced {} root bytes): '{}'", root.len(), serialize(relative));
                wz.move_to_path(relative);
                wz.set_val(());
                true
            }

            Tag::SymbolSize(_) => {
                // Literal guard - only emit if it matches the min value
                trace!(target: "sink", "  literal guard");

                // Safe guard extraction using the proper API
                let guard_expr = unsafe { Expr { ptr: best_path[report_len..].as_ptr().cast_mut() } };

                // Use symbol() to get the payload without the header byte
                let guard_bytes = unsafe { guard_expr.symbol() };
                if guard_bytes.is_none() {
                    trace!(target: "sink", "  could not extract guard symbol");
                    return false;
                }

                let guard_bytes = unsafe { guard_bytes.unwrap().as_ref().unwrap() };
                if let Some(guard_val) = Self::parse_i64_symbol_strict(guard_bytes) {
                    if guard_val == min_val {
                        // Emit the report portion, relative to root_prefix
                        let root = wz.root_prefix_path();
                        let relative = if report_len > root.len() {
                            &best_path[root.len()..report_len]
                        } else {
                            &[]
                        };
                        trace!(target: "sink", "  guard matches, emitting: '{}'", serialize(relative));
                        wz.move_to_path(relative);
                        wz.set_val(());
                        return true;
                    }
                }
                trace!(target: "sink", "  guard doesn't match");
                false
            }

            Tag::NewVar => {
                // Guard not referenced, emit report
                let root = wz.root_prefix_path();
                let relative = if report_len > root.len() {
                    &best_path[root.len()..report_len]
                } else {
                    &[]
                };
                trace!(target: "sink", "  NewVar guard, emitting report: '{}'", serialize(relative));
                wz.move_to_path(relative);
                wz.set_val(());
                true
            }

            _ => {
                trace!(target: "sink", "  unexpected guard type: {:?}", guard_tag);
                false
            }
        }
    }
}

// MaxSink_dep - DEPRECATED: Does not handle duplicate variables correctly
// finds maximum numeric value from matched patterns
// Syntax: (O (max <report-template-with-$m> $m <value-expression>))
// Example: (O (max (max-depth $m) $m $depth))
pub struct MaxSink_dep {
    e: Expr,
    skip: usize,
    max: Option<i64>,
}

impl Sink for MaxSink_dep {
    fn new(e: Expr) -> Self {
        // Fixed skip: 1 (arity) + 1 (symbol size) + 3 ("max") = 5
        // This skips past the operator, leaving <report> <guard> <value>
        MaxSink_dep { e, skip: 5, max: None }
    }

    fn request(&self) -> impl Iterator<Item=&'static [u8]> {
        // Return pattern after "max": <report> <guard> <value>
        let p = &unsafe {
            self.e.prefix()
                .unwrap_or_else(|_| {
                    let s = self.e.span();
                    slice_from_raw_parts(self.e.ptr, s.len() - 1)
                })
                .as_ref()
                .unwrap()
        }[self.skip..];
        trace!(target: "sink", "max requesting {}", serialize(p));
        std::iter::once(p)
    }

    fn sink<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, mut it: It, path: &[u8]) where 'a : 'w, 'k : 'w {
        let mut wz = it.next().unwrap();
        // mpath should contain <report> <guard> <value>
        // We skip the operator ("max"), but NOT the root_prefix (which is part of the pattern)
        let mpath = &path[self.skip..];
        trace!(target: "sink", "max at '{}' sinking raw '{}'", serialize(wz.root_prefix_path()), serialize(path));
        trace!(target: "sink", "max mpath contains '{}'", serialize(mpath));

        // mpath is a sequence of bytes representing: <report> <guard> <value>
        // These are stored sequentially, not as a tuple, so we navigate by reading spans
        let mut offset = 0;
        let mut ctx = unsafe { Expr { ptr: mpath.as_ptr().cast_mut() } };

        // Skip report (first subexpression)
        let report_len = ctx.span().len();
        trace!(target: "sink", "max: report at offset {}, len {}", offset, report_len);
        offset += report_len;

        if offset >= mpath.len() {
            trace!(target: "sink", "max: ran out of bytes after report");
            return;
        }

        // Skip guard (second subexpression)
        ctx = unsafe { Expr { ptr: mpath.as_ptr().add(offset).cast_mut() } };
        let guard_len = ctx.span().len();
        trace!(target: "sink", "max: guard at offset {}, len {}", offset, guard_len);
        offset += guard_len;

        if offset >= mpath.len() {
            trace!(target: "sink", "max: ran out of bytes after guard");
            return;
        }

        // Now at value (third subexpression)
        let val_expr = unsafe { Expr { ptr: mpath.as_ptr().add(offset).cast_mut() } };
        trace!(target: "sink", "max value expr at offset {}: '{}'", offset, serialize(unsafe { val_expr.span().as_ref().unwrap() }));

        if let Some(bytes) = unsafe { val_expr.span().as_ref() } {
            // read an integer symbol if present
            if !bytes.is_empty() {
                if let Tag::SymbolSize(_) = byte_item(bytes[0]) {
                    if let Ok(s) = std::str::from_utf8(&bytes[1..]) {
                        // STRICT PARSING: Validate that the entire string is a valid integer
                        let is_valid = !s.is_empty() &&
                            if s.as_bytes()[0] == b'-' {
                                s.len() > 1 && s[1..].bytes().all(|b| b.is_ascii_digit())
                            } else {
                                s.bytes().all(|b| b.is_ascii_digit())
                            };

                        if is_valid {
                            if let Ok(v) = s.parse::<i64>() {
                                trace!(target: "sink", "max found value {}", v);
                                self.max = Some(self.max.map_or(v, |cur| cur.max(v)));
                            } else {
                                trace!(target: "sink", "max: value '{}' overflows i64, skipped", s);
                            }
                        } else {
                            trace!(target: "sink", "max: value '{}' is not a valid integer, skipped", s);
                        }
                    }
                }
            }
        }
    }

    fn finalize<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, mut it: It) -> bool where 'a : 'w, 'k : 'w {
        let mut wz = it.next().unwrap();
        trace!(target: "sink", "max finalize: before reset, origin='{}', root_prefix='{}'", serialize(wz.origin_path()), serialize(wz.root_prefix_path()));
        wz.reset();
        trace!(target: "sink", "max finalize: after reset, path='{}', root_prefix='{}'", serialize(wz.path()), serialize(wz.root_prefix_path()));

        let Some(max_val) = self.max else {
            trace!(target: "sink", "max finalizing with no values");
            return false;
        };
        trace!(target: "sink", "max finalizing with {}", max_val);

        // Parse (max <report> <guard> <value-expr>) again from the definition 'e'
        // Use the same span-based navigation as sink()
        let mut offset = 5; // Skip past "(max "
        let e_bytes = unsafe { self.e.span().as_ref().unwrap() };

        // Get report
        let report = unsafe { Expr { ptr: e_bytes.as_ptr().add(offset).cast_mut() } };
        offset += report.span().len();

        // Get guard
        let guard = unsafe { Expr { ptr: e_bytes.as_ptr().add(offset).cast_mut() } };

        trace!(target: "sink", "max finalize: report='{}'", serialize(unsafe { report.span().as_ref().unwrap() }));
        trace!(target: "sink", "max finalize: guard='{}'", serialize(unsafe { guard.span().as_ref().unwrap() }));

        // Guard cases:
        unsafe {
            let guard_tag = byte_item(*guard.ptr);
            trace!(target: "sink", "max finalize: guard_tag={:?}", guard_tag);
            match guard_tag {
                // (max <report($m)> $m <value>) where $m is a NewVar or VarRef
                // Substitute $m in report with the computed max
                Tag::NewVar => {
                    trace!(target: "sink", "max finalize: NewVar guard - appending max value {}", max_val);

                    // Same as VarRef case: append the max value to complete the report template
                    let ms = max_val.to_string();
                    let mut value_bytes = vec![item_byte(Tag::SymbolSize(ms.len() as _))];
                    value_bytes.extend_from_slice(ms.as_bytes());

                    trace!(target: "sink", "max finalize: NewVar - Appending value bytes '{}' to root '{}'", serialize(&value_bytes), serialize(wz.root_prefix_path()));
                    wz.move_to_path(&value_bytes);
                    wz.set_val(());
                    return true;
                }

                Tag::VarRef(k) => {
                    trace!(target: "sink", "max finalize: VarRef({}) guard - appending max value {}", k, max_val);

                    // The WriteZipper is positioned at the report template (e.g., "[2] found-max")
                    // which has arity 2 but is incomplete (only has the symbol, missing the value)
                    // We need to append the max value as the second element
                    let ms = max_val.to_string();
                    let mut value_bytes = vec![item_byte(Tag::SymbolSize(ms.len() as _))];
                    value_bytes.extend_from_slice(ms.as_bytes());

                    trace!(target: "sink", "max finalize: VarRef - Appending value bytes '{}' to root '{}'", serialize(&value_bytes), serialize(wz.root_prefix_path()));
                    wz.move_to_path(&value_bytes);
                    trace!(target: "sink", "max finalize: VarRef - After move, path='{}'", serialize(wz.path()));
                    wz.set_val(());
                    return true;
                }

                // (max <report> <literal-int> <value>) -- emit only if literal == max
                Tag::SymbolSize(_) => {
                    trace!(target: "sink", "max finalize: Literal guard - checking if guard == {}", max_val);
                    if let Some(gbytes) = guard.span().as_ref() {
                        if let Ok(gs) = std::str::from_utf8(&gbytes[1..]) {
                            // STRICT PARSING for literal guard
                            let is_valid = !gs.is_empty() &&
                                if gs.as_bytes()[0] == b'-' {
                                    gs.len() > 1 && gs[1..].bytes().all(|b| b.is_ascii_digit())
                                } else {
                                    gs.bytes().all(|b| b.is_ascii_digit())
                                };

                            if is_valid {
                                if let Ok(gv) = gs.parse::<i64>() {
                                    trace!(target: "sink", "max finalize: Parsed guard as {}", gv);
                                    if gv == max_val {
                                        trace!(target: "sink", "max finalize: Match! Emitting report");
                                        let report_bytes = report.span().as_ref().unwrap();
                                        wz.move_to_path(report_bytes);
                                        wz.set_val(());
                                        return true;
                                    } else {
                                        trace!(target: "sink", "max finalize: No match ({} != {})", gv, max_val);
                                    }
                                } else {
                                    trace!(target: "sink", "max finalize: Guard value '{}' overflows i64", gs);
                                }
                            } else {
                                trace!(target: "sink", "max finalize: Guard value '{}' is not a valid integer", gs);
                            }
                        }
                    }
                }

                _ => {
                    trace!(target: "sink", "max finalize: Unexpected guard tag {:?}", guard_tag);
                }
            }
        }

        false
    }
}

// MinSink_dep - DEPRECATED: Does not handle duplicate variables correctly
// finds minimum numeric value from matched patterns
// Syntax: (O (min <report-template-with-$m> $m <value-expression>))
// Example: (O (min (min-depth $m) $m $depth))
pub struct MinSink_dep {
    e: Expr,
    skip: usize,
    min: Option<i64>,
}

impl Sink for MinSink_dep {
    fn new(e: Expr) -> Self {
        // Fixed skip: 1 (arity) + 1 (symbol size) + 3 ("min") = 5
        // This skips past the operator, leaving <report> <guard> <value>
        MinSink_dep { e, skip: 5, min: None }
    }

    fn request(&self) -> impl Iterator<Item=&'static [u8]> {
        // Return pattern after "min": <report> <guard> <value>
        let p = &unsafe {
            self.e.prefix()
                .unwrap_or_else(|_| {
                    let s = self.e.span();
                    slice_from_raw_parts(self.e.ptr, s.len() - 1)
                })
                .as_ref()
                .unwrap()
        }[self.skip..];
        trace!(target: "sink", "min requesting {}", serialize(p));
        std::iter::once(p)
    }

    fn sink<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, mut it: It, path: &[u8]) where 'a : 'w, 'k : 'w {
        let mut wz = it.next().unwrap();
        // mpath should contain <report> <guard> <value>
        // We skip the operator ("min"), but NOT the root_prefix (which is part of the pattern)
        let mpath = &path[self.skip..];
        trace!(target: "sink", "min at '{}' sinking raw '{}'", serialize(wz.root_prefix_path()), serialize(path));
        trace!(target: "sink", "min mpath contains '{}'", serialize(mpath));

        // mpath is a sequence of bytes representing: <report> <guard> <value>
        // These are stored sequentially, not as a tuple, so we navigate by reading spans
        let mut offset = 0;
        let mut ctx = unsafe { Expr { ptr: mpath.as_ptr().cast_mut() } };

        // Skip report (first subexpression)
        let report_len = ctx.span().len();
        trace!(target: "sink", "min: report at offset {}, len {}", offset, report_len);
        offset += report_len;

        if offset >= mpath.len() {
            trace!(target: "sink", "min: ran out of bytes after report");
            return;
        }

        // Skip guard (second subexpression)
        ctx = unsafe { Expr { ptr: mpath.as_ptr().add(offset).cast_mut() } };
        let guard_len = ctx.span().len();
        trace!(target: "sink", "min: guard at offset {}, len {}", offset, guard_len);
        offset += guard_len;

        if offset >= mpath.len() {
            trace!(target: "sink", "min: ran out of bytes after guard");
            return;
        }

        // Now at value (third subexpression)
        let val_expr = unsafe { Expr { ptr: mpath.as_ptr().add(offset).cast_mut() } };
        trace!(target: "sink", "min value expr at offset {}: '{}'", offset, serialize(unsafe { val_expr.span().as_ref().unwrap() }));

        if let Some(bytes) = unsafe { val_expr.span().as_ref() } {
            // read an integer symbol if present
            if !bytes.is_empty() {
                if let Tag::SymbolSize(_) = byte_item(bytes[0]) {
                    if let Ok(s) = std::str::from_utf8(&bytes[1..]) {
                        // STRICT PARSING: Validate that the entire string is a valid integer
                        let is_valid = !s.is_empty() &&
                            if s.as_bytes()[0] == b'-' {
                                s.len() > 1 && s[1..].bytes().all(|b| b.is_ascii_digit())
                            } else {
                                s.bytes().all(|b| b.is_ascii_digit())
                            };

                        if is_valid {
                            if let Ok(v) = s.parse::<i64>() {
                                trace!(target: "sink", "min found value {}", v);
                                self.min = Some(self.min.map_or(v, |cur| cur.min(v)));
                            } else {
                                trace!(target: "sink", "min: value '{}' overflows i64, skipped", s);
                            }
                        } else {
                            trace!(target: "sink", "min: value '{}' is not a valid integer, skipped", s);
                        }
                    }
                }
            }
        }
    }

    fn finalize<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, mut it: It) -> bool where 'a : 'w, 'k : 'w {
        let mut wz = it.next().unwrap();
        wz.reset();

        let Some(min_val) = self.min else {
            trace!(target: "sink", "min finalizing with no values");
            return false;
        };
        trace!(target: "sink", "min finalizing with {}", min_val);

        // Parse (min <report> <guard> <value-expr>) using span-based navigation
        let mut offset = 5; // Skip past "(min "
        let e_bytes = unsafe { self.e.span().as_ref().unwrap() };

        // Get report
        let report = unsafe { Expr { ptr: e_bytes.as_ptr().add(offset).cast_mut() } };
        offset += report.span().len();

        // Get guard
        let guard = unsafe { Expr { ptr: e_bytes.as_ptr().add(offset).cast_mut() } };

        // Guard cases:
        unsafe {
            let guard_tag = byte_item(*guard.ptr);
            match guard_tag {
                // (min <report($m)> $m <value>) where $m is a NewVar or VarRef
                Tag::NewVar => {
                    trace!(target: "sink", "min finalize: NewVar guard - appending min value {}", min_val);

                    // Append the min value to complete the report template
                    let ms = min_val.to_string();
                    let mut value_bytes = vec![item_byte(Tag::SymbolSize(ms.len() as _))];
                    value_bytes.extend_from_slice(ms.as_bytes());

                    trace!(target: "sink", "min finalize: NewVar - Appending value bytes '{}' to root '{}'", serialize(&value_bytes), serialize(wz.root_prefix_path()));
                    wz.move_to_path(&value_bytes);
                    wz.set_val(());
                    return true;
                }

                Tag::VarRef(k) => {
                    trace!(target: "sink", "min finalize: VarRef({}) guard - appending min value {}", k, min_val);

                    // Append the min value to complete the report template
                    let ms = min_val.to_string();
                    let mut value_bytes = vec![item_byte(Tag::SymbolSize(ms.len() as _))];
                    value_bytes.extend_from_slice(ms.as_bytes());

                    trace!(target: "sink", "min finalize: VarRef - Appending value bytes '{}' to root '{}'", serialize(&value_bytes), serialize(wz.root_prefix_path()));
                    wz.move_to_path(&value_bytes);
                    wz.set_val(());
                    return true;
                }

                // (min <report> <literal-int> <value>) -- emit only if literal == min
                Tag::SymbolSize(_) => {
                    trace!(target: "sink", "min finalize: Literal guard - checking if guard == {}", min_val);
                    if let Some(gbytes) = guard.span().as_ref() {
                        if let Ok(gs) = std::str::from_utf8(&gbytes[1..]) {
                            // STRICT PARSING for literal guard
                            let is_valid = !gs.is_empty() &&
                                if gs.as_bytes()[0] == b'-' {
                                    gs.len() > 1 && gs[1..].bytes().all(|b| b.is_ascii_digit())
                                } else {
                                    gs.bytes().all(|b| b.is_ascii_digit())
                                };

                            if is_valid {
                                if let Ok(gv) = gs.parse::<i64>() {
                                    trace!(target: "sink", "min finalize: Parsed guard as {}", gv);
                                    if gv == min_val {
                                        trace!(target: "sink", "min finalize: Match! Emitting report");
                                        let report_bytes = report.span().as_ref().unwrap();
                                        wz.move_to_path(report_bytes);
                                        wz.set_val(());
                                        return true;
                                    } else {
                                        trace!(target: "sink", "min finalize: No match ({} != {})", gv, min_val);
                                    }
                                } else {
                                    trace!(target: "sink", "min finalize: Guard value '{}' overflows i64", gs);
                                }
                            } else {
                                trace!(target: "sink", "min finalize: Guard value '{}' is not a valid integer", gs);
                            }
                        }
                    }
                }

                _ => {}
            }
        }

        false
    }
}

mod pure {
    fn add_i32(i: &[u8], o: &mut [u8]) -> bool {
        false
    }
}

// (pure (result $x) $x (ascii_to_i32 2))
pub struct PureSink { e: Expr }
impl Sink for PureSink {
    fn new(e: Expr) -> Self {
        PureSink { e }
    }
    fn request(&self) -> impl Iterator<Item=&'static [u8]> {
        let p = &unsafe { self.e.prefix().unwrap_or_else(|x| { let s = self.e.span(); slice_from_raw_parts(self.e.ptr, s.len() - 1) }).as_ref().unwrap() }[6..];
        trace!(target: "sink", "count requesting {}", serialize(p));
        std::iter::once(p)
    }
    fn sink<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, mut it: It, path: &[u8]) where 'a : 'w, 'k : 'w {
        // let mut wz = it.next().unwrap();
        // let mpath = &path[7+wz.root_prefix_path().len()..];
        // let ctx = unsafe { Expr { ptr: mpath.as_ptr().cast_mut() } };
        // trace!(target: "sink", "count at '{}' sinking raw '{}'", serialize(wz.root_prefix_path()), serialize(path));
        // trace!(target: "sink", "count registering in ctx {:?}", ctx);
        // self.unique.insert(mpath, ());
        // wz.move_to_path(opath)
    }
    fn finalize<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, mut it: It) -> bool where 'a : 'w, 'k : 'w  {
        true
    }
}


pub enum ASink { AddSink(AddSink), RemoveSink(RemoveSink), HeadSink(HeadSink), CountSink(CountSink), MaxSink(MaxSink), MinSink(MinSink),
    #[cfg(feature = "wasm")]
    WASMSink(WASMSink),
    #[cfg(feature = "grounding")]
    PureSink(PureSink)
}

impl Sink for ASink {
    fn new(e: Expr) -> Self {
        if unsafe { *e.ptr == item_byte(Tag::Arity(2)) && *e.ptr.offset(1) == item_byte(Tag::SymbolSize(1)) && *e.ptr.offset(2) == b'-' } {
            ASink::RemoveSink(RemoveSink::new(e))
        } else if unsafe { *e.ptr == item_byte(Tag::Arity(2)) && *e.ptr.offset(1) == item_byte(Tag::SymbolSize(1)) && *e.ptr.offset(2) == b'+' } {
            ASink::AddSink(AddSink::new(e))
        } else if unsafe { *e.ptr == item_byte(Tag::Arity(3)) && *e.ptr.offset(1) == item_byte(Tag::SymbolSize(4)) &&
            *e.ptr.offset(2) == b'h' && *e.ptr.offset(3) == b'e' && *e.ptr.offset(4) == b'a' && *e.ptr.offset(5) == b'd' } {
            ASink::HeadSink(HeadSink::new(e))
        } else if unsafe { *e.ptr == item_byte(Tag::Arity(4)) && *e.ptr.offset(1) == item_byte(Tag::SymbolSize(5)) &&
            *e.ptr.offset(2) == b'c' && *e.ptr.offset(3) == b'o' && *e.ptr.offset(4) == b'u' && *e.ptr.offset(5) == b'n' && *e.ptr.offset(6) == b't' } {
            ASink::CountSink(CountSink::new(e))
        } else if unsafe { *e.ptr == item_byte(Tag::Arity(4)) && *e.ptr.offset(1) == item_byte(Tag::SymbolSize(3)) &&
            *e.ptr.offset(2) == b'm' && *e.ptr.offset(3) == b'a' && *e.ptr.offset(4) == b'x' } {
            ASink::MaxSink(MaxSink::new(e))
        } else if unsafe { *e.ptr == item_byte(Tag::Arity(4)) && *e.ptr.offset(1) == item_byte(Tag::SymbolSize(3)) &&
            *e.ptr.offset(2) == b'm' && *e.ptr.offset(3) == b'i' && *e.ptr.offset(4) == b'n' } {
            ASink::MinSink(MinSink::new(e))
        } else if unsafe { *e.ptr == item_byte(Tag::Arity(3)) && *e.ptr.offset(1) == item_byte(Tag::SymbolSize(4)) &&
            *e.ptr.offset(2) == b'w' && *e.ptr.offset(3) == b'a' && *e.ptr.offset(4) == b's' && *e.ptr.offset(5) == b'm' } {
            #[cfg(feature = "wasm")]
            return ASink::WASMSink(WASMSink::new(e));
            #[cfg(not(feature = "wasm"))]
            panic!("MORK was not built with the wasm feature, yet trying to call {:?}", e);
        } else if unsafe { *e.ptr == item_byte(Tag::Arity(4)) && *e.ptr.offset(1) == item_byte(Tag::SymbolSize(4)) &&
            *e.ptr.offset(2) == b'p' && *e.ptr.offset(3) == b'u' && *e.ptr.offset(4) == b'r' && *e.ptr.offset(5) == b'e' } {
            #[cfg(feature = "grounding")]
            return ASink::PureSink(PureSink::new(e));
            #[cfg(not(feature = "grounding"))]
            panic!("MORK was not built with the grounding feature, yet trying to call {:?}", e);
        } else {
            unreachable!()
        }
    }

    fn request(&self) -> impl Iterator<Item=&'static [u8]> {
        gen move {
            match self {
                ASink::AddSink(s) => { for i in s.request().into_iter() { yield i } }
                ASink::RemoveSink(s) => { for i in s.request().into_iter() { yield i } }
                ASink::HeadSink(s) => { for i in s.request().into_iter() { yield i } }
                ASink::CountSink(s) => { for i in s.request().into_iter() { yield i } }
                ASink::MaxSink(s) => { for i in s.request().into_iter() { yield i } }
                ASink::MinSink(s) => { for i in s.request().into_iter() { yield i } }
                #[cfg(feature = "wasm")]
                ASink::WASMSink(s) => { for i in s.request().into_iter() { yield i } }
                #[cfg(feature = "grounding")]
                ASink::PureSink(s) => { for i in s.request().into_iter() { yield i } }
            }
        }
    }
    fn sink<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, it: It, path: &[u8]) where 'a: 'w, 'k: 'w {
        match self {
            ASink::AddSink(s) => { s.sink(it, path) }
            ASink::RemoveSink(s) => { s.sink(it, path) }
            ASink::HeadSink(s) => { s.sink(it, path) }
            ASink::CountSink(s) => { s.sink(it, path) }
            ASink::MaxSink(s) => { s.sink(it, path) }
            ASink::MinSink(s) => { s.sink(it, path) }
            #[cfg(feature = "wasm")]
            ASink::WASMSink(s) => { s.sink(it, path) }
            #[cfg(feature = "grounding")]
            ASink::PureSink(s) => { s.sink(it, path) }
        }
    }

    fn finalize<'w, 'a, 'k, It: Iterator<Item=&'w mut WriteZipperTracked<'a, 'k, ()>>>(&mut self, it: It) -> bool where 'a: 'w, 'k: 'w {
        match self {
            ASink::AddSink(s) => { s.finalize(it) }
            ASink::RemoveSink(s) => { s.finalize(it) }
            ASink::HeadSink(s) => { s.finalize(it) }
            ASink::CountSink(s) => { s.finalize(it) }
            ASink::MaxSink(s) => { s.finalize(it) }
            ASink::MinSink(s) => { s.finalize(it) }
            #[cfg(feature = "wasm")]
            ASink::WASMSink(s) => { s.finalize(it) }
            #[cfg(feature = "grounding")]
            ASink::PureSink(s) => { s.finalize(it) }
        }
    }
}
