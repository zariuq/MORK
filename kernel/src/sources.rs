use log::trace;
use pathmap::arena_compact::{ACTMmapZipper};
use pathmap::PathMap;
use pathmap::zipper::{PolyZipper, PrefixZipper, ReadZipperUntracked, Zipper, ZipperAbsolutePath, ZipperIteration, ZipperMoving, DependentProductZipperG, ReadZipperOwned, ZipperSubtries};
use mork_expr::{byte_item, destruct, item_byte, serialize, Expr, Tag};
use mork_expr::macros::SerializableExpr;

// Note: FilterZipper implementation is more complex than expected
// For now, we'll keep the simpler approach where predicates don't work as intended
// The max/min non-numeric handling is the main improvement achieved

pub(crate) enum ResourceRequest {
    BTM(&'static [u8]),
    ACT(&'static str)
}

pub(crate) enum Resource<'trie, 'path> {
    BTM(ReadZipperUntracked<'trie, 'path, ()>),
    ACT(ACTMmapZipper<'trie, ()>)
}

pub trait Source {
    // step 1: parsing the source
    fn new(e: Expr) -> Self;
    // step 2: request access to resources before running
    fn request(&self) -> impl Iterator<Item=ResourceRequest>;
    // step 3: create the factor in the product/the (virtual) zipper for the source
    fn source<'trie, 'path, It : Iterator<Item=Resource<'trie, 'path>>>(&self, it: It) -> AFactor<'trie, ()> where 'path : 'trie;
}

struct BTMSource {
    e: Expr
}
impl Source for BTMSource {
    fn new(e: Expr) -> Self {
        BTMSource { e }
    }

    fn request(&self) -> impl Iterator<Item=ResourceRequest> {
        std::iter::once(ResourceRequest::BTM([].as_slice()))
    }

    fn source<'trie, 'path, It: Iterator<Item=Resource<'trie, 'path>>>(&self, mut it: It) -> AFactor<'trie, ()> where 'path : 'trie {
        // (I (BTM <pat1>) (ACT <filename> <pat2>)
        //    --factor1--  -----factor2---------
        // prefix: '[2] BTM'
        static PREFIX: [u8; 5] = [item_byte(Tag::Arity(2)), item_byte(Tag::SymbolSize(3)), b'B', b'T', b'M'];
        let Resource::BTM(rz) = it.next().unwrap() else { unreachable!() };
        let rz = PrefixZipper::new(&PREFIX[..], rz);
        AFactor::PosSource(rz)
    }
}

struct ACTSource {
    e: Expr,
    act: &'static str
}
impl Source for ACTSource {
    fn new(e: Expr) -> Self {
        destruct!(e, ("ACT" {act: &str} se), {
            return ACTSource{ e, act }
        }, _err => { panic!("act not the right shape") });
    }

    fn request(&self) -> impl Iterator<Item=ResourceRequest> {
        std::iter::once(ResourceRequest::ACT(self.act))
    }

    fn source<'trie, 'path, It: Iterator<Item=Resource<'trie, 'path>>>(&self, mut it: It) -> AFactor<'trie, ()> where 'path : 'trie {
        // prefix: '[3] ACT <filename>'
        static CONSTANT_PREFIX: [u8; 5] = [item_byte(Tag::Arity(3)), item_byte(Tag::SymbolSize(3)), b'A', b'C', b'T'];
        let Resource::ACT(rz) = it.next().unwrap() else { unreachable!() };
        let mut prefix = vec![];
        prefix.extend_from_slice(&CONSTANT_PREFIX[..]);
        prefix.push(item_byte(Tag::SymbolSize( (self.act.size() as u8) - 1)));
        prefix.extend_from_slice(self.act.as_bytes());
        trace!(target: "source", "prefix {}", serialize(&prefix[..]));
        let rz = PrefixZipper::new(prefix, rz);
        AFactor::ACTSource(rz)
    }
}


struct CmpSource {
    e: Expr,
    cmp: usize
}

impl CmpSource {
    fn policy(ctx: (usize, PathMap<()>), p: &[u8], c: usize) -> ((usize, PathMap<()>), Option<ReadZipperOwned<()>>) {
        let (cmp, map) = ctx;
        if c == 0 {
            if cmp == 0 {
                trace!(target: "source", "== enrolling at {}", serialize(p));
                // bug: de bruijn levels broken, easy fix: shift the copy of p by introductions(p)
                let e = Expr{ ptr: p.as_ptr().cast_mut() };
                let mut qv = p.to_vec();
                e.shift(e.newvars() as _, &mut mork_expr::ExprZipper::new(Expr{ ptr: qv.as_mut_ptr() }));
                ((cmp, map), Some(PathMap::single(&qv[..], ()).into_read_zipper(&[])))
            } else if cmp == 1 {
                let mut cloned = map.clone();
                let present = cloned.remove(p).is_some();
                trace!(target: "source", "!= enrolling (present {:?}) at {}", present, serialize(p));
                ((cmp, map), Some(cloned.into_read_zipper(&[])))
            } else {
                unreachable!()
            }
        } else {
            ((cmp, map), None)
        }
    }
}

impl Source for CmpSource {
    fn new(e: Expr) -> Self {
        let cmp = if unsafe { *e.ptr.offset(2) == b'=' } {
            assert!(unsafe { *e.ptr.offset(3) == b'=' });
            0
        } else if unsafe { *e.ptr.offset(2) == b'!' } {
            assert!(unsafe { *e.ptr.offset(3) == b'=' });
            1
        } else {
            // todo < <= #=
            panic!("comparator not implemented")
        };
        // trace!(target: "source", "cmp {cmp} source");
        CmpSource { e, cmp }
    }

    fn request(&self) -> impl Iterator<Item=ResourceRequest> {
        std::iter::once(ResourceRequest::BTM([].as_slice()))
    }

    fn source<'trie, 'path, It: Iterator<Item=Resource<'trie, 'path>>>(&self, mut it: It) -> AFactor<'trie, ()> where 'path : 'trie {
        static EQ_PREFIX: [u8; 4] = [item_byte(Tag::Arity(3)), item_byte(Tag::SymbolSize(2)), b'=', b'='];
        static NE_PREFIX: [u8; 4] = [item_byte(Tag::Arity(3)), item_byte(Tag::SymbolSize(2)), b'!', b'='];
        let Resource::BTM(rz) = it.next().unwrap() else { unreachable!() };
        let map = rz.make_map().unwrap();
        let rz = DependentProductZipperG::new_enroll(rz, (self.cmp, map),
            CmpSource::policy as for<'a> fn((usize, PathMap<()>), &'a [u8], usize) -> ((usize, PathMap<()>), Option<ReadZipperOwned<()>>));
        let rz = PrefixZipper::new(
            if self.cmp == 0 { &EQ_PREFIX[..] }
            else if self.cmp == 1 { &NE_PREFIX[..] }
            else { unreachable!() }, rz);
        AFactor::CmpSource(rz)
    }
}

/* // IsNumberSource - checks if an expression is a numeric value
// Syntax: (is-number <expr>)
// NOTE: Implementation deferred - needs complex FilterZipper work
struct IsNumberSource {
    e: Expr
}

impl IsNumberSource {
    // Prefix for (is-number $x) pattern
    const PREFIX: [u8; 11] = [
        item_byte(Tag::Arity(2)),
        item_byte(Tag::SymbolSize(9)),
        b'i', b's', b'-', b'n', b'u', b'm', b'b', b'e', b'r'
    ];

    fn arg_is_number(full_path: &[u8]) -> bool {
        // Skip the prefix to get to the argument
        if full_path.len() <= Self::PREFIX.len() {
            return false;
        }

        let arg = &full_path[Self::PREFIX.len()..];

        // Check if the argument is a numeric symbol
        if arg.is_empty() {
            return false;
        }

        // The argument should be a single symbol
        if let Ok(Tag::SymbolSize(n)) = mork_expr::maybe_byte_item(arg[0]) {
            let n = n as usize;
            if arg.len() >= n + 1 {
                let sym = &arg[1..n + 1];
                // Check if it's a valid number
                if sym.is_empty() { return false; }
                if sym[0] == b'-' {
                    if sym.len() <= 1 || !sym[1..].iter().all(|b| b.is_ascii_digit()) {
                        return false;
                    }
                } else if !sym.iter().all(|b| b.is_ascii_digit()) {
                    return false;
                }
                // Try to parse as i64
                let s = unsafe { core::str::from_utf8_unchecked(sym) };
                return s.parse::<i64>().is_ok();
            }
        }
        false
    }
}

impl Source for IsNumberSource {
    fn new(e: Expr) -> Self {
        IsNumberSource { e }
    }

    fn request(&self) -> impl Iterator<Item=ResourceRequest> {
        // Request the entire BTM space to filter
        std::iter::once(ResourceRequest::BTM([].as_slice()))
    }

    fn source<'trie, 'path, It: Iterator<Item=Resource<'trie, 'path>>>(&self, mut it: It) -> AFactor<'trie, ()> where 'path : 'trie {
        let Resource::BTM(rz) = it.next().unwrap() else { unreachable!() };

        // Simple prefix approach for now - filtering would need more complex zipper work
        let prefixed = PrefixZipper::new(&Self::PREFIX[..], rz);

        AFactor::IsNumberSource(prefixed)
    }
}

// IsVarSource - checks if an expression is a variable
// Syntax: (is-var <expr>)
struct IsVarSource {
    e: Expr
}

impl IsVarSource {
    // Prefix for (is-var $x) pattern
    const PREFIX: [u8; 8] = [
        item_byte(Tag::Arity(2)),
        item_byte(Tag::SymbolSize(6)),
        b'i', b's', b'-', b'v', b'a', b'r'
    ];

    fn arg_is_var(full_path: &[u8]) -> bool {
        // Skip the prefix to get to the argument
        if full_path.len() <= Self::PREFIX.len() {
            return false;
        }

        let arg = &full_path[Self::PREFIX.len()..];

        // Check if the argument is a variable
        if arg.is_empty() {
            return false;
        }

        // Check the first byte of the argument
        match mork_expr::maybe_byte_item(arg[0]) {
            Ok(Tag::NewVar) => true,      // $
            Ok(Tag::VarRef(_)) => true,    // _1, _2, etc
            _ => false
        }
    }
}

impl Source for IsVarSource {
    fn new(e: Expr) -> Self {
        IsVarSource { e }
    }

    fn request(&self) -> impl Iterator<Item=ResourceRequest> {
        std::iter::once(ResourceRequest::BTM([].as_slice()))
    }

    fn source<'trie, 'path, It: Iterator<Item=Resource<'trie, 'path>>>(&self, mut it: It) -> AFactor<'trie, ()> where 'path : 'trie {
        let Resource::BTM(rz) = it.next().unwrap() else { unreachable!() };

        // Simple prefix approach for now - filtering would need more complex zipper work
        let prefixed = PrefixZipper::new(&Self::PREFIX[..], rz);

        AFactor::IsVarSource(prefixed)
    }
}
*/

pub enum ASource { PosSource(BTMSource), ACTSource(ACTSource), CmpSource(CmpSource) }

#[derive(PolyZipper)]
pub enum AFactor<'trie, V: Clone + Send + Sync + Unpin + 'static = ()> {
    PosSource(PrefixZipper<'trie, ReadZipperUntracked<'trie, 'trie, V>>),
    ACTSource(PrefixZipper<'trie, ACTMmapZipper<'trie, V>>),
    CmpSource(PrefixZipper<'trie, DependentProductZipperG<'trie, ReadZipperUntracked<'trie, 'trie, V>,
        ReadZipperOwned<V>, V, (usize, PathMap<()>), for<'a> fn((usize, PathMap<()>), &'a [u8], usize) -> ((usize, PathMap<()>), Option<ReadZipperOwned<V>>)>>),
}

impl Source for ASource {
    fn new(e: Expr) -> Self {
        if unsafe { *e.ptr == item_byte(Tag::Arity(2)) && *e.ptr.offset(1) == item_byte(Tag::SymbolSize(3)) && *e.ptr.offset(2) == b'B' && *e.ptr.offset(3) == b'T' && *e.ptr.offset(4) == b'M' } {
            ASource::PosSource(BTMSource::new(e))
        } else if unsafe { *e.ptr == item_byte(Tag::Arity(3)) && *e.ptr.offset(1) == item_byte(Tag::SymbolSize(3)) && *e.ptr.offset(2) == b'A' && *e.ptr.offset(3) == b'C' && *e.ptr.offset(4) == b'T' } {
            ASource::ACTSource(ACTSource::new(e))
        } else if unsafe { *e.ptr == item_byte(Tag::Arity(3)) && *e.ptr.offset(1) == item_byte(Tag::SymbolSize(2)) && (*e.ptr.offset(2) == b'=' || *e.ptr.offset(2) == b'!') && *e.ptr.offset(3) == b'=' } {
            ASource::CmpSource(CmpSource::new(e))
        } else {
            unreachable!()
        }
    }

    fn request(&self) -> impl Iterator<Item=ResourceRequest> {
        gen move {
            match self {
                ASource::PosSource(s) => { for i in s.request().into_iter() { yield i } }
                ASource::ACTSource(s) => { for i in s.request().into_iter() { yield i } }
                ASource::CmpSource(s) => { for i in s.request().into_iter() { yield i } }
            }
        }
    }

    fn source<'trie, 'path, It: Iterator<Item=Resource<'trie, 'path>>>(&self, mut it: It) -> AFactor<'trie, ()> where 'path : 'trie {
        match self {
            ASource::PosSource(s) => { s.source(it) }
            ASource::ACTSource(s) => { s.source(it) }
            ASource::CmpSource(s) => { s.source(it) }
        }
    }
}
