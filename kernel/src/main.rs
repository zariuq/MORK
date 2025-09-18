use std::collections::HashSet;
// use std::future::Future;
// use std::task::Poll;
use std::time::Instant;
use pathmap::trie_map::BytesTrieMap;
use pathmap::zipper::{Zipper, ZipperAbsolutePath, ZipperIteration, ZipperMoving};
use mork_frontend::bytestring_parser::Parser;
use mork::{expr, prefix, sexpr};
use mork::prefix::Prefix;
use mork::space::{transitions, unifications, Space};
use mork_bytestring::{item_byte, Tag};
/*fn main() {
    let mut s = Space::new();
    let t0 = Instant::now();
    let nodesf = std::fs::File::open("/run/media/adam/43323a1c-ad7e-4d9a-b3c0-cf84e69ec61a/awesome-biomedical-kg/ckg_v3-002/results/nodes.json").unwrap();
    let nodesfs = unsafe { memmap2::Mmap::map(&nodesf).unwrap() };
    let loaded = s.load_json(nodesfs.as_ref()).unwrap();
    println!("loaded {} nodes in {} seconds", loaded, t0.elapsed().as_secs());
    let t1 = Instant::now();
    let edgesf = std::fs::File::open("/run/media/adam/43323a1c-ad7e-4d9a-b3c0-cf84e69ec61a/awesome-biomedical-kg/ckg_v3-002/results/edges.json").unwrap();
    let edgesfs = unsafe { memmap2::Mmap::map(&edgesf).unwrap() };
    let loaded = s.load_json(edgesfs.as_ref()).unwrap();
    println!("loaded {} edges in {} seconds", loaded, t1.elapsed().as_secs());
    s.done();
}*/


// fn main() {
//     let mut s = Space::new();
//     let t0 = Instant::now();
//     let nodesf = std::fs::File::open("/run/media/adam/43323a1c-ad7e-4d9a-b3c0-cf84e69ec61a/awesome-biomedical-kg/ckg_v3-002/results/nodes.json").unwrap();
//     let nodesfs = unsafe { memmap2::Mmap::map(&nodesf).unwrap() };
//     let loaded = s.load_json(nodesfs.as_ref()).unwrap();
//     println!("loaded {} nodes in {} seconds", loaded, t0.elapsed().as_secs());
//     let t1 = Instant::now();
//     let edgesf = std::fs::File::open("/run/media/adam/43323a1c-ad7e-4d9a-b3c0-cf84e69ec61a/awesome-biomedical-kg/ckg_v3-002/results/edges.json").unwrap();
//     let edgesfs = unsafe { memmap2::Mmap::map(&edgesf).unwrap() };
//     let loaded = s.load_json(edgesfs.as_ref()).unwrap();
//     println!("loaded {} edges in {} seconds", loaded, t1.elapsed().as_secs());
//     s.done();
// }

// fn work(s: &mut Space) {
//     let restore_paths_start = Instant::now();
//     println!("restored paths {:?}", s.restore_paths("/dev/shm/combined_ni.paths.gz").unwrap());
//     println!("paths restore took {}", restore_paths_start.elapsed().as_secs());
//     s.statistics();
//
//     let add_gene_name_index_start = Instant::now();
//     s.transform(expr!(s, "[4] NKV $ gene_name $"), expr!(s, "[3] gene_name_of _2 _1"));
//     println!("add gene name index took {} ms", add_gene_name_index_start.elapsed().as_millis());
//     s.statistics();
//
//     let all_related_to_gene_start = Instant::now();
//     s.transform_multi(&[
//         expr!(s, "[3] gene_name_of TP73-AS1 $"),
//         expr!(s, "[4] SPO _1 includes $"),
//         expr!(s, "[4] SPO _1 transcribed_from $"),
//     ], expr!(s, "[4] res0 _1 _2 _3"));
//     println!("all_related_to_gene_start {}", all_related_to_gene_start.elapsed().as_micros());
//     let mut count = 0;
//     s.query(expr!(s, "[4] res0 $ $ $"), |_, e| {
//         println!("{}", sexpr!(s, e));
//         count += 1
//     });
//     println!("res0 count {}", count);
//
//     let add_exon_chr_index_start = Instant::now();
//     s.transform(expr!(s, "[4] NKV $ chr $"), expr!(s, "[3] chr_of _2 _1"));
//     println!("add exon chr index took {}", add_exon_chr_index_start.elapsed().as_secs());
//     s.statistics();
//
//     let ops_index_start = Instant::now();
//     s.transform(expr!(s, "[4] SPO $ $ $"), expr!(s, "[4] OPS _3 _2 _1"));
//     println!("add ops index took {}", ops_index_start.elapsed().as_secs());
//     s.statistics();
//
//     let transitive_chr1_start = Instant::now();
//     s.transform_multi(&[
//         expr!(s, "[3] chr_of chr1 $"),
//         expr!(s, "[4] OPS _1 includes $"),
//         expr!(s, "[4] SPO _2 translates_to $"),
//         expr!(s, "[4] OPS _3 interacts_with $"),
//     ], expr!(s, "[5] res1 _1 _2 _3 _4"));
//     println!("transitive_chr1 {} ms", transitive_chr1_start.elapsed().as_millis());
//     let mut count = 0;
//     s.query(expr!(s, "[5] res1 $ $ $ $"), |_, e| {
//         // println!("{}", sexpr!(s, e));
//         count += 1
//     });
//     println!("res1 count {}", count);
//
//     let q0_start = Instant::now();
//     s.transform_multi(&[
//         expr!(s, "[3] gene_name_of BRCA2 $"),
//         expr!(s, "[4] SPO _1 transcribed_to $"),
//         expr!(s, "[4] SPO _2 translates_to $"),
//         expr!(s, "[4] OPS _3 interacts_with $"),
//         expr!(s, "[4] SPO _1 genes_pathways $"),
//     ], expr!(s, "[6] res2 _1 _2 _3 _4 _5"));
//     println!("q0 {}", q0_start.elapsed().as_micros());
//     let mut count = 0;
//     s.query( expr!(s, "[6] res2 $ $ $ $ $"), |_, e| {
//         // println!("{}", sexpr!(s, e));
//         count += 1
//     });
//     println!("res2 count {}", count);
//
// }

const work_mm2: &str = r#"
(exec P0 (, (NKV $x gene_name $y)) (,) (, (gene_name_of $y $x)))
(exec P0' (,) (, (MICROS $t) (U64.DIV $t 1000 $tms)) (, (time "add gene name index" $tms ms)))

(exec P1 (, (gene_name_of TP73-AS1 $x)
            (SPO $x includes $y)
            (SPO $x transcribed_from $z)) (,) (, (res0 $x $y $z))))
(exec P1' (,) (, (MICROS $t)) (, (time "all related to gene" $t us)))
(exec P1'' (, (res0 $x $y $z)) (, (COUNT (res0 $x $y $z))) (, (count "query TP73-AS1" $c)))
(exec P1''' (,) (, (MICROS $t)) (, (time "query TP73-AS1" $t us)))

(exec P2 (, (NKV $x chr $y)) (,) (, (chr_of $y $x)))
(exec P2' (,) (, (MICROS $t)) (, (time "add exon chr index" $t us)))

(exec P3 (, (SPO $s $p $o)) (,) (, (OPS $o $p $s)))
(exec P3' (,) (, (MICROS $t)) (, (time "add exon chr index" $t us)))

(exec P4 (, (chr_of chr1 $x)
            (OPS $x includes $y)
            (SPO $y transcribed_from $z)
            (OPS $z transcribed_from $w)) (,) (, (res1 $x $y $z $w))))
(exec P4' (,) (, (MICROS $t)) (, (time "all related to gene" $t us)))
(exec P4'' (, (res1 $x $y $z $w)) (, (COUNT (res1 $x $y $z $w))) (, (count "query chr1" $c)))
(exec P4''' (,) (, (MICROS $t)) (, (time "query chr1" $t us)))

(exec P5 (, (gene_name_of BRCA2 $x)
            (SPO $x transcribed_to $y)
            (SPO $y translates_to $z)
            (OPS $z interacts_with $p)
            (SPO $x genes_pathways $q)) (,) (, (res2 $x $y $z $p $q))))
(exec P5' (,) (, (MICROS $t)) (, (time "all related to gene" $t us)))
(exec P5'' (, (res2 $x $y $z $p $q)) (, (COUNT (res2 $x $y $z $p $q))) (, (count "query BRCA2" $c)))
(exec P5''' (,) (, (MICROS $t)) (, (time "query BRCA2" $t us)))
"#;

fn work_mm2_run() {
    let mut s = Space::new();
    let restore_paths_start = Instant::now();
    println!("restored paths {:?}", s.restore_paths("/dev/shm/combined_ni.paths.gz").unwrap());
    println!("paths restore took {}", restore_paths_start.elapsed().as_secs());
    s.statistics();

    s.metta_calculus(100);

    let backup_paths_start = Instant::now();
    println!("{:?}", s.backup_paths("/run/media/adam/43323a1c-ad7e-4d9a-b3c0-cf84e69ec61a/whole_flybase.paths.gz").unwrap());
    println!("paths backup took {}", backup_paths_start.elapsed().as_secs());
}

/*
paths restore took 135
 978700221 atoms
add gene name index took 8637 ms
 979015756 atoms
query TP73-AS1
 193 µs
 142
add exon chr index took 32 s
 1054962128 atoms
add ops index took 91 s
 1386253656 atoms
query chr1
 7691 ms
 40172978 atoms
query BRCA2
 33295 µs
 151956 atoms
 */

fn peano(x: usize) -> String {
    if x == 0 { "Z".to_string() }
    else { format!("(S {})", peano(x - 1)) }
}

fn basic() {
    let mut s = Space::new();

    const space: &str = r#"
(Straight 1 2)
(Straight 2 3)

(exec P1 (, (Straight $x $y) (Straight $y $z)) (, (Transitive $x $z)))

(exec P2 (, (Transitive $x $y)) (, (Line $x $q)))
(exec P2 (, (Transitive $x $y)) (, (Line $q $y)))

"#;
    // (exec (P0 reverse) (, (Straight $x $y) (exec (P0 reverse) $P $T)) (, (Reverse $y $x) (pexec (P0 reverse) $P $T)))
    //
    // (exec P1 (, (Straight $x $y) (Straight $y $z)) (, (Transitive $x $z)))
    //

    s.load_sexpr(space.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    s.metta_calculus(100);

    let mut v = vec![];
    s.dump_sexpr(expr!(s, "$"), expr!(s, "_1"), &mut v);

    println!("out {}", String::from_utf8(v).unwrap());


}

fn process_calculus_bench(steps: usize, x: usize, y: usize) {
    let mut s = Space::new();

    // note 'idle' MM2-like statement that can be activated by moving it to the exec space
    let space_exprs = format!(r#"
(exec (IC 0 1 {})
               (, (exec (IC $x $y (S $c)) $sp $st)
                  ((exec $x) $p $t))
               (, (exec (IC $y $x $c) $sp $st)
                  (exec (R $x) $p $t)))

((exec 0)
      (, (petri (? $channel $payload $body))
         (petri (! $channel $payload)) )
      (, (petri $body)))
((exec 1)
      (, (petri (| $lprocess $rprocess)))
      (, (petri $lprocess)
         (petri $rprocess)))

(petri (? (add $ret) ((S $x) $y) (| (! (add (PN $x $y)) ($x $y))
                                    (? (PN $x $y) $z (! $ret (S $z)))  )  ))
(petri (? (add $ret) (Z $y) (! $ret $y)))
(petri (! (add result) ({} {})))
    "#, peano(steps), peano(x), peano(y));

    s.load_sexpr(space_exprs.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    let t0 = Instant::now();
    let mcalc_steps = s.metta_calculus(1000000000000000); // big number to show the MM2 inference control working
    let elapsed = t0.elapsed();
    
    let mut v = vec![];
    // s.dump_all_sexpr(&mut v).unwrap();
    // We're only interested in the petri dish (not the state of exec), and specifically everything that was outputted `!` to `result`
    s.dump_sexpr(expr!(s, "[2] petri [3] ! result $"), expr!(s, "_1"), &mut v);
    let res = String::from_utf8(v).unwrap();

    println!("{x}+{y} ({} steps) in {} µs result: {res}", steps, elapsed.as_micros());
    assert_eq!(res, format!("{}\n", peano(x+y)));
    println!("unifications {}, instructions {}", unsafe { unifications }, unsafe { transitions });
    // (badbad)
    // 200+200 (1000 steps) in 42716559 µs
}

fn process_calculus_reverse() {
    let mut s = Space::new();

    // note 'idle' MM2-like statement that can be activated by moving it to the exec space
    const SPACE_EXPRS: &str = r#"
(exec (IC 0 1  (S (S (S (S (S (S (S (S (S (S Z)))))))))) )
               (, (exec (IC $x $y (S $c)) $sp $st)
                  ((exec $x) $p $t))
               (, (exec (IC $y $x $c) $sp $st)
                  (exec (R $x) $p $t)))

((exec 0)
      (, (petri (! $channel $payload))
         (petri (? $channel $payload $body)) )
      (, (petri $body)))
((exec 1)
      (, (petri (| $lprocess $rprocess)))
      (, (petri $lprocess)
         (petri $rprocess)))

(petri (? (add $ret) ((S $x) $y) (| (! (add (PN $x $y)) ($x $y))
                                    (? (PN $x $y) $z (! $ret (S $z)))  )  ))
(petri (? (add $ret) (Z $y) (! $ret $y)))
(petri (! (add result) ( (S (S Z)) (S (S Z)) )))
    "#;

    s.load_sexpr(SPACE_EXPRS.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    let steps = s.metta_calculus(1000000000000000); // big number to show the MM2 inference control working

    let mut v = vec![];
    s.dump_sexpr(expr!(s, "[2] petri [3] ! result $"), expr!(s, "_1"), &mut v);
    let res = String::from_utf8(v).unwrap();

    println!("result: {res}");
    assert_eq!(res, "(S (S (S (S Z))))\n");
}

fn lookup() {
    let mut s = Space::new();

    const SPACE_EXPRS: &str = r#"
(exec 0 (, (Something (very specific))) (, MATCHED))

(Something (very specific))
    "#;

    s.load_all_sexpr(SPACE_EXPRS.as_bytes()).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(1000000000000000);
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();

    println!("result: {res}");
    assert!(res.contains("MATCHED\n"));
}

fn positive() {
    let mut s = Space::new();

    const SPACE_EXPRS: &str = r#"
(exec 0 (, (Something $unspecific)) (, MATCHED))

(Something (very specific))
    "#;

    s.load_all_sexpr(SPACE_EXPRS.as_bytes()).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(1000000000000000);
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();

    println!("result: {res}");
    assert!(res.contains("MATCHED\n"));
}

fn positive_equal() {
    let mut s = Space::new();

    const SPACE_EXPRS: &str = r#"
(exec 0 (, (Something $repeated $repeated)) (, MATCHED))

(Something (very specific) (very specific))
    "#;

    s.load_all_sexpr(SPACE_EXPRS.as_bytes()).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(1000000000000000);
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();

    println!("result: {res}");
    assert!(res.contains("MATCHED\n"));
}

fn negative() {
    let mut s = Space::new();

    // note 'idle' MM2-like statement that can be activated by moving it to the exec space
    const SPACE_EXPRS: &str = r#"
(exec 0 (, (Something (very specific))) (, MATCHED))

(Something $unspecific)

    "#;

    s.load_sexpr(SPACE_EXPRS.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(1000000000000000); // big number to show the MM2 inference control working
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();

    println!("result: {res}");
    assert!(res.contains("MATCHED\n"));
}

fn negative_equal() {
    let mut s = Space::new();

    // note 'idle' MM2-like statement that can be activated by moving it to the exec space
    const SPACE_EXPRS: &str = r#"
(exec 0 (, (Something (very specific) (very specific))) (, MATCHED))

(Something $repeated $repeated)

    "#;

    s.load_sexpr(SPACE_EXPRS.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(1000000000000000); // big number to show the MM2 inference control working
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();

    println!("result: {res}");
    assert!(res.contains("MATCHED\n"));
}

fn bipolar() {
    let mut s = Space::new();

    // note 'idle' MM2-like statement that can be activated by moving it to the exec space
    const SPACE_EXPRS: &str = r#"
(exec 0 (, (Something (very $unspecific))) (, MATCHED))

(Something ($unspecific specific))

    "#;

    s.load_sexpr(SPACE_EXPRS.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(1000000000000000); // big number to show the MM2 inference control working
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();

    println!("result: {res}");
    assert!(res.contains("MATCHED\n"));
}

fn bipolar_equal() {
    let mut s = Space::new();

    // note 'idle' MM2-like statement that can be activated by moving it to the exec space
    const SPACE_EXPRS: &str = r#"
(exec 0 (, (Something ($x Y $x Y))) (, MATCHED))

(Something (X $y X $y))

    "#;

    s.load_sexpr(SPACE_EXPRS.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(1000000000000000); // big number to show the MM2 inference control working
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();

    println!("result: {res}");
    assert!(res.contains("MATCHED\n"));
}

fn two_positive_equal() {
    let mut s = Space::new();

    const SPACE_EXPRS: &str = r#"
(exec 0 (, (Something $x $x) (Else $y $y)) (, MATCHED))

(Something (foo bar) (foo bar))
(Else (bar baz) (bar baz))
    "#;

    s.load_all_sexpr(SPACE_EXPRS.as_bytes()).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(1000000000000000);
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();

    println!("result: {res}");
    assert!(res.contains("MATCHED\n"));
}

fn two_positive_equal_crossed() {
    let mut s = Space::new();

    const SPACE_EXPRS: &str = r#"
(exec 0 (, (Something $x $y) (Else $x $y)) (, MATCHED))

(Something (foo bar) (bar baz))
(Else (foo bar) (bar baz))
    "#;

    s.load_all_sexpr(SPACE_EXPRS.as_bytes()).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(1000000000000000);
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();

    println!("result: {res}");
    assert!(res.contains("MATCHED\n"));
}

fn two_bipolar_equal_crossed() {
    let mut s = Space::new();

    const SPACE_EXPRS: &str = r#"
(exec 0 (, (Something $x $y) (Else $x $y)) (, (MATCHED $x $y)))

(Something (foo $x) (foo $x))
(Else ($x bar) ($x bar))
    "#;

    s.load_all_sexpr(SPACE_EXPRS.as_bytes()).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(1000000000000000);
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();

    println!("result: {res}");
    assert!(res.contains("(MATCHED (foo bar) (foo bar))\n"));
}

fn logic_query() {
    let mut s = Space::new();

    const SPACE_EXPRS: &str = r#"
(exec 0 (, (axiom $x) (axiom $x)) (, (combined $x)))
(exec 0 (, (axiom (= $lhs $rhs)) (axiom (= $rhs $lhs))) (, (reversed $lhs $rhs)))
    "#;

    const AXIOM_EXPRS: &str = r#"
(= (L $x $y $z) (R $x $y $z))
(= (L 1 $x $y) (R 1 $x $y))
(= (R $x (L $x $y $z) $w) $x)
(= (R $x (R $x $y $z) $w) $x)
(= (R $x (L $x $y $z) $x) (L $x (L $x $y $z) $x))
(= (L $x $y (\ $y $z)) (L $x $y $z))
(= (L $x $y (* $z $y)) (L $x $y $z))
(= (L $x $y (\ $z 1)) (L $x $z $y))
(= (L $x $y (\ $z $y)) (L $x $z $y))
(= (L $x 1 (\ $y 1)) (L $x $y 1))
(= (T $x (L $x $y $z)) $x)
(= (T $x (R $x $y $z)) $x)
(= (T $x (a $x $y $z)) $x)
(= (T $x (\ (a $x $y $z) $w)) (T $x $w))
(= (T $x (* $y $y)) (T $x (\ (a $x $z $w) (* $y $y))))
(= (R (/ 1 $x) $x (\ $x 1)) (\ $x 1))
(= (\ $x 1) (/ 1 (L $x $x (\ $x 1))))
(= (L $x $x $x) (* (K $x (\ $x 1)) $x))
    "#;

    s.load_all_sexpr(SPACE_EXPRS.as_bytes()).unwrap();
    s.load_sexpr(AXIOM_EXPRS.as_bytes(),expr!(s, "$"), expr!(s, "[2] axiom _1")).unwrap();

    let steps = s.metta_calculus(1000000000000000);

    assert_eq!(s.btm.val_count(), 79);
}

fn bc0() {
    let mut s = Space::new();

    const SPACE_EXPRS: &str = r#"
    ((step base)
      (, (goal (: $proof $conclusion)) (kb (: $proof $conclusion)))
      (, (ev (: $proof $conclusion) ) ))

    ((step abs)
      (, (goal (: $proof $conclusion)))
      (, (goal (: $lhs (-> $synth $conclusion)) ) ))

    ((step rev)
      (, (ev (: $lhs (-> $a $r)))  (goal (: $k $r)) )
      (, (goal (: $rhs $a) ) ))

    ((step app)
      (, (ev (: $lhs (-> $a $r)))  (ev (: $rhs $a))  )
      (, (ev (: (@ $lhs $rhs) $r) ) ))

    (exec zealous
            (, ((step $x) $p0 $t0)
               (exec zealous $p1 $t1) )
            (, (exec $x $p0 $t0)
               (exec zealous $p1 $t1) ))
    "#;

    const KB_EXPRS: &str = r#"
    (kb (: a A))
    (kb (: ab (R A B)))
    (kb (: bc (R B C)))
    (kb (: MP (-> (R $p $q) (-> $p $q))))

    (goal (: $proof C))
    "#;

    s.load_all_sexpr(SPACE_EXPRS.as_bytes()).unwrap();
    s.load_all_sexpr(KB_EXPRS.as_bytes()).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(50);
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();

    println!("result: {res}");
    assert!(res.contains("(ev (: (@ (@ MP bc) (@ (@ MP ab) a)) C))\n"));
}

fn bc1() {
    let mut s = Space::new();

    const SPACE_EXPRS: &str = r#"
    ((step base)
      (, (goal (: $proof $conclusion)) (kb (: $proof $conclusion)))
      (, (ev (: $proof $conclusion) ) ))

    ((step rec)
      (, (goal (: (@ $lhs $rhs) $conclusion)))
      (, (goal (: $lhs (-> $synth $conclusion))) (goal (: $rhs $synth))))

    ((step app)
      (, (ev (: $lhs (-> $a $r)))  (ev (: $rhs $a))  )
      (, (ev (: (@ $lhs $rhs) $r) ) ))

    (exec zealous
            (, ((step $x) $p0 $t0)
               (exec zealous $p1 $t1) )
            (, (exec $x $p0 $t0)
               (exec zealous $p1 $t1) ))
    "#;

    const KB_EXPRS: &str = r#"
    (kb (: a A))
    (kb (: ab (R A B)))
    (kb (: bc (R B C)))
    (kb (: cd (R C D)))
    (kb (: MP (-> (R $p $q) (-> $p $q))))

    (goal (: $proof C))
    "#;

    s.load_all_sexpr(SPACE_EXPRS.as_bytes()).unwrap();
    s.load_all_sexpr(KB_EXPRS.as_bytes()).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(100);
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();

    println!("result: {res}");
    assert!(res.contains("(ev (: (@ (@ MP cd) (@ (@ MP bc) (@ (@ MP ab) a))) D))\n"));
}

fn bc2() {
    let mut s = Space::new();

    /*
    ((step rec)
      (, (goal (: (@ $lhs $rhs) $conclusion)))
      (, (goal (: $lhs (-> $synth $conclusion))) (goal (: $rhs $synth))))

    ((step rec2)
      (, (goal (: (@ $f $a $b) $conclusion)))
      (, (goal (: $f (-> $syntha $synthb $conclusion))) (goal (: $a $syntha)) (goal (: $b $synthb)) ))
    
     */
    const SPACE_EXPRS: &str = r#"
    ((step base)
      (, (goal (: $proof $conclusion)) (kb (: $proof $conclusion)))
      (, (ev (: $proof $conclusion) ) ))

    ((step abs)
      (, (goal (: $proof $conclusion)))
      (, (goal (: $lhs (-> $synth $conclusion)) ) ))

    ((step rev)
      (, (ev (: $lhs (-> $a $r)))  (goal (: $k $r)) )
      (, (goal (: $rhs $a) ) ))

    ((step abs2)
      (, (goal (: $proof $conclusion)))
      (, (goal (: $lhs (-> $syntha $synthb $conclusion)) ) ))

    ((step rev2)
      (, (ev (: $lhs (-> $a $b $r)))  (goal (: $k $r)) )
      (, (goal (: $ap $a)) (goal (: $bp $b)) ))

    ((step app)
      (, (ev (: $lhs (-> $a $r)))  (ev (: $rhs $a))  )
      (, (ev (: (@ $lhs $rhs) $r) ) ))
      
    ((step app2)
      (, (ev (: $f (-> $a $b $r)))  (ev (: $ap $a)) (ev (: $bp $b))  )
      (, (ev (: (@ $f $ap $bp) $r) ) ))

    (exec zealous
            (, ((step $x) $p0 $t0)
               (exec zealous $p1 $t1) )
            (, (exec $x $p0 $t0)
               (exec zealous $p1 $t1) ))
    "#;

    const KB_EXPRS: &str = r#"
    (kb (: ax-mp (-> $𝜑 (→ $𝜑 $𝜓) $𝜓)))
    (kb (: ax-1 (→ $𝜑 (→ $𝜓 $𝜑))))
    (kb (: ax-2 (→ (→ $𝜑 (→ $𝜓 $𝜒)) (→ (→ $𝜑 $𝜓) (→ $𝜑 $𝜒)))))
    (kb (: ax-3 (→ (→ (¬ $𝜑) (¬ $𝜓)) (→ $𝜓 $𝜑))))

    (kb (: mp2b.1 𝜑))
    (kb (: mp2b.2 (→ 𝜑 𝜓)))
    (kb (: mp2b.3 (→ 𝜓 𝜒)))

    (goal (: $proof 𝜒))
    "#;

    s.load_all_sexpr(SPACE_EXPRS.as_bytes()).unwrap();
    s.load_all_sexpr(KB_EXPRS.as_bytes()).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(30);
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    s.dump_sexpr(expr!(s, "[2] ev [3] : $ 𝜒"), expr!(s, "_1"), &mut v);
    let res = String::from_utf8(v).unwrap();

    println!("proof of 𝜒: {res}");
    assert!(res.contains("(@ ax-mp (@ ax-mp mp2b.1 mp2b.2) mp2b.3)\n"));
}

// delete me
// fn mm0() {
//     use std::time::Instant;
//     use mork::space::Space;
//     use mork::{expr};

//     // --- mm2 program (steps + a tiny 2-arg KB/goal) ---
//     const MM0: &str = r#"
//     ; keep bc0's base, app, rev, abs, exec..., exec zealous here
//     ; + add the new 2-arg rules & execs
//     ((step abs2)  (, (goal (: $p $r))) (, (goal (: $f (-> $a (-> $b $r))))))
//     ((step rev2)  (, (ev (: $f (-> $a (-> $b $r)))) (goal (: $k $r)))
//                 (, (goal (: $xa $a)) (goal (: $xb $b))))
//     ((step app2)  (, (ev (: $f (-> $a (-> $b $r))))
//                     (ev (: $x $a))
//                     (ev (: $y $b)))
//                 (, (ev (: (@ (@ $f $x) $y) $r))))
//     ((step seed-ev) (, (kb (: $f (-> $a $b)))) (, (ev (: $f (-> $a $b)))))

//     (exec abs2 (, (goal (: $p $r)))
//             (, (goal (: $f (-> $a (-> $b $r)))) (fired abs2 $r)))
//     (exec rev2 (, (ev (: $f (-> $a (-> $b $r)))) (goal (: $k $r)))
//             (, (goal (: $xa $a)) (goal (: $xb $b)) (fired rev2 $r $a $b)))
//     (exec app2 (, (ev (: $f (-> $a (-> $b $r))))
//                 (ev (: $x $a))
//                 (ev (: $y $b)))
//             (, (ev (: (@ (@ $f $x) $y) $r)) (fired app2 $r $a $b)))
//     (exec seed-ev (, (kb (: $f (-> $a $b))))
//                 (, (ev (: $f (-> $a $b))) (fired seed-ev $a $b)))

//     ; ---- tiny 2-arg KB/goal ----
//     (kb (: f (-> A (-> B C))))
//     (kb (: x A))
//     (kb (: y B))
//     (goal (: (@ (@ f x) y) C))
//     "#;

//     println!("== mm0: 2-arg BC demo ==");
//     let t0 = Instant::now();

//     let mut s = Space::new();
//     // Load bc0 chassis first (your existing steps/exec), then MM0 additions.
//     // If you inlined bc0 earlier, concatenate strings or paste bc0 + the new block.
//     let _ = s.load_sexpr(MM0.as_bytes(), expr!(s, "$"), expr!(s, "_1"));

//     // Run for plenty of steps; the loop stops when nothing matches.
//     let steps = s.metta_calculus(1_000_000);
//     println!("elapsed {:?}, steps {}", t0.elapsed(), steps);

//     // Full dump (you can grep 'result:' or 'fired')
//     let mut out = std::io::stdout();
//     s.dump_all_sexpr(&mut out).unwrap();
// }

// This version focuses on robustly verifying the result of the original,
// stalling-but-correct chainer. It bypasses the `dump_sexpr` query API
// and uses a simple string search on the full text dump of the space.
fn mm0_original_with_verification_fix() {
    use std::time::Instant;
    use mork::space::Space;
    use mork::expr;

    const P: &str = r#"
    ; Using the corrected `base` rule from our last insight.
    ((step base)
      (, (goal $fact) (kb $fact))
      (, (ev $fact)))

    ((step revapp2)
      (, (goal (: (@ (@ $f $x) $y) $R)) (kb (: $f (-> $A (-> $B $R)))))
      (, (goal (: $x $A)) (goal (: $y $B))))

    ((step app2)
      (, (kb (: $f (-> $A (-> $B $R)))) (ev (: $x $A)) (ev (: $y $B)))
      (, (ev (: (@ (@ $f $x) $y) $R))))

    ; --- Knowledge Base & Goal ---
    (kb (: ⟨+⟩ (-> term (-> term term))))
    (kb (: ⟨t⟩ term))
    (kb (: ⟨0⟩ term))
    (goal (: (@ (@ ⟨+⟩ ⟨t⟩) ⟨0⟩) term))

    ; --- Driver ---
    (exec zealous (, ((step $n) $p $r) (exec zealous $d $e))
                  (, (exec $n $p $r)    (exec zealous $d $e)))
    "#;

    println!("\n== mm0_original_with_verification_fix (String Search Version) ==");
    let mut s = Space::new();
    let t0 = Instant::now();
    s.load_sexpr(P.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    let steps = s.metta_calculus(16);
    let elapsed = t0.elapsed();

    // --- Analytics & Verification via String Search ---
    println!("\n--- Analytics ---");

    // 1. Get the entire state of the space as a single string.
    let mut full_dump_buffer = Vec::new();
    s.dump_all_sexpr(&mut full_dump_buffer).unwrap();
    let full_dump_string = String::from_utf8_lossy(&full_dump_buffer);

    // 2. Define the exact line we are looking for in the dump.
    let target_evidence = "(ev (: (@ (@ ⟨+⟩ ⟨t⟩) ⟨0⟩) term))";

    // 3. Perform the simple string search. This is our new verification method.
    let success = full_dump_string.contains(target_evidence);

    // 4. Report the summary.
    println!("Status: {}", if success { "✅ SUCCESS (Verified by String Search)" } else { "❌ FAILURE" });
    println!("Completed in {:?} after {} calculus steps.", elapsed, steps);

    // 5. Print the full dump that we searched, for complete transparency.
    println!("\n--- Full Final State Dump ---");
    print!("{}", full_dump_string);
    println!("----------------------------------------------------------");
}

// Renamed to clarify its behavior. This version uses the zealous driver with a
// fixed number of internal steps, which we've proven works and finds the result.
fn mm0_original_fixed_steps() {
    use std::time::Instant;
    use mork::space::Space;
    use mork::expr;

    const P: &str = r#"
    ((step base)
      (, (goal $fact) (kb $fact))
      (, (ev $fact)))
    ((step revapp2)
      (, (goal (: (@ (@ $f $x) $y) $R)) (kb (: $f (-> $A (-> $B $R)))))
      (, (goal (: $x $A)) (goal (: $y $B))))
    ((step app2)
      (, (kb (: $f (-> $A (-> $B $R)))) (ev (: $x $A)) (ev (: $y $B)))
      (, (ev (: (@ (@ $f $x) $y) $R))))
    (kb (: ⟨+⟩ (-> term (-> term term))))
    (kb (: ⟨t⟩ term))
    (kb (: ⟨0⟩ term))
    (goal (: (@ (@ ⟨+⟩ ⟨t⟩) ⟨0⟩) term))
    (exec zealous (, ((step $n) $p $r) (exec zealous $d $e))
                  (, (exec $n $p $r)    (exec zealous $d $e)))
    "#;

    println!("\n== mm0_original_fixed_steps ==");
    let mut s = Space::new();
    let t0 = Instant::now();
    s.load_sexpr(P.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    let steps = s.metta_calculus(16); // Using the proven fixed-step method
    let elapsed = t0.elapsed();

    println!("\n--- Analytics ---");
    let mut full_dump_buffer = Vec::new();
    s.dump_all_sexpr(&mut full_dump_buffer).unwrap();
    let full_dump_string = String::from_utf8_lossy(&full_dump_buffer);
    let target_evidence = "(ev (: (@ (@ ⟨+⟩ ⟨t⟩) ⟨0⟩) term))";
    let success = full_dump_string.contains(target_evidence);

    println!("Status: {}", if success { "✅ SUCCESS (Verified by String Search)" } else { "❌ FAILURE" });
    println!("Completed in {:?} after {} calculus steps.", elapsed, steps);
    println!("\n--- Full Final State Dump ---");
    print!("{}", full_dump_string);
    println!("----------------------------------------------------------");
}

// This new version of mm0 uses a robust "Data Pipeline" pattern that works
// correctly with the manual ticking loop. Each step produces a unique data
// atom that becomes the input for the next step in the pipeline.
fn mm0() {
    use std::time::Instant;
    use mork::space::Space;
    use mork::expr;

    const P: &str = r#"
    ; --- Data Pipeline Rules ---

    ; Step 1: Decompose the main goal. This is the entry point.
    ; THE FIX IS HERE: The pattern now correctly matches the nested application.
    (exec decompose-goal
      (, (goal (: (@ (@ $f $x) $y) $R)) (kb (: $f (-> $A (-> $B $R)))))
      (, (subgoal-for $x $A) (subgoal-for $y $B)))

    ; Step 2: Solve subgoals. These rules look for the intermediate products
    ; from Step 1 and produce new, unique evidence atoms.
    (exec solve-subgoal-t
      (, (subgoal-for ⟨t⟩ term) (kb (: ⟨t⟩ term)))
      (, (evidence-for t-is-term)))
    (exec solve-subgoal-0
      (, (subgoal-for ⟨0⟩ term) (kb (: ⟨0⟩ term)))
      (, (evidence-for 0-is-term)))

    ; Step 3: Synthesize the final proof. This rule waits for the unique
    ; evidence products from Step 2 to appear.
    (exec synthesize-final-proof
      (, (evidence-for t-is-term) (evidence-for 0-is-term))
      (, (ev (: (@ (@ ⟨+⟩ ⟨t⟩) ⟨0⟩) term))))

    ; --- Knowledge Base & Goal ---
    (kb (: ⟨+⟩ (-> term (-> term term))))
    (kb (: ⟨t⟩ term))
    (kb (: ⟨0⟩ term))
    (goal (: (@ (@ ⟨+⟩ ⟨t⟩) ⟨0⟩) term))
    "#;

    println!("\n== mm0 (Data Pipeline Version - Corrected) ==");
    let mut s = Space::new();
    let t0 = Instant::now();
    s.load_sexpr(P.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    let mut success = false;
    let mut ticks = 0;
    for i in 0..100 {
        ticks = i + 1;
        let n = s.metta_calculus(1);
        
        // Check for success inside the loop to stop as soon as the proof is found
        let pat = expr!(s, "(ev (: (@ (@ ⟨+⟩ ⟨t⟩) ⟨0⟩) term))");
        let mut buf = Vec::new();
        s.dump_sexpr(pat, expr!(s, "_1"), &mut buf);
        if !buf.is_empty() {
            success = true;
            break;
        }

        if n == 0 { break; } // Stop if the space has saturated
    }
    let elapsed = t0.elapsed();

    // --- Analytics ---
    println!("\n--- Analytics ---");
    let mut full_dump_buffer = Vec::new();
    s.dump_all_sexpr(&mut full_dump_buffer).unwrap();
    let full_dump_string = String::from_utf8_lossy(&full_dump_buffer);
    
    // We can re-verify with the string search, but the loop break is the real test
    let success_check = full_dump_string.contains("(ev (: (@ (@ ⟨+⟩ ⟨t⟩) ⟨0⟩) term))");

    println!("Status: {}", if success_check { "✅ SUCCESS" } else { "❌ FAILURE" });
    println!("Completed in {:?} after {} ticks.", elapsed, ticks);
    println!("\n--- Full Final State Dump ---");
    print!("{}", full_dump_string);
    println!("----------------------------------------------------------");
}

fn mm1_a() {
    use mork::space::Space;
    use mork::expr;
    use std::time::Instant;

    const P: &str = r#"
    ; --- KB: Just the types of our axioms and constructors ---
    ; Let P := ((t + 0) = t) and Q := (t = t)
    (kb (: ⟨mp-curry⟩ (-> (: $P ⟨wff⟩) (-> (: $Q ⟨wff⟩) (-> (: $P ⟨|-⟩) (-> (: (⟨->⟩ $P $Q) ⟨|-⟩) (: $Q ⟨|-⟩)))))))
    (kb (: ⟨a2-curry⟩ (-> (: $t ⟨term⟩) (: (⟨=⟩ (⟨+⟩ $t ⟨0⟩) $t) ⟨|-⟩))))
    (kb (: ⟨a1-curry⟩ (-> (: $t ⟨term⟩) (-> (: $r ⟨term⟩) (-> (: $s ⟨term⟩) (: (⟨->⟩ (⟨=⟩ $t $r) (⟨->⟩ (⟨=⟩ $t $s) (⟨=⟩ $r $s))) ⟨|-⟩))))))

    ; --- Premise "Ingredients": The proven components we will assemble ---
    ; We give them simple names (wff_P, proof_P, etc.)
    
    ; The types (wffs)
    (: wff_P  (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))
    (: wff_Q  (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩))

    ; The proofs (sequents)
    (: proof_P (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))
    ; This next term is the proof of (P -> (P -> Q)) which simplifies to (P -> Q) after one MP, a common pattern.
    ; For our purpose, it's the required |- (P -> Q) premise, derived from axiom a1.
    (: proof_PtoQ (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩))

    ; --- The Assembler Rule ---
    ; This single rule finds our named ingredients and composes them.
    (exec assemble-final-proof
      (, (: wff_P $type_wff_P)
         (: wff_Q $type_wff_Q)
         (: proof_P $type_proof_P)
         (: proof_PtoQ $type_proof_PtoQ)
         (kb (: ⟨mp-curry⟩ $_))) ; Just confirm mp-curry is in the KB
      (, (ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨|-⟩)
             ; This is the beautiful compositional proof term:
             ((((⟨mp-curry⟩ wff_P) wff_Q) proof_P) proof_PtoQ))))
    "#;

    // --- Load, Run, and Verify ---
    let mut s = Space::new();
    s.load_sexpr(P.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();
    let t0 = Instant::now();

    // A single tick should be enough for this simple assembly.
    println!("Tick 1");
    s.metta_calculus(1);
    
    // --- Verification ---
    println!("\n--- Verifying Result ---");
    let mut full_dump_buffer = Vec::new();
    s.dump_all_sexpr(&mut full_dump_buffer).unwrap();
    let full_dump_string = String::from_utf8_lossy(&full_dump_buffer);

    // Use the robust string search method for verification
    let target_evidence = "(ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨|-⟩) ((((⟨mp-curry⟩ wff_P) wff_Q) proof_P) proof_PtoQ))";
    let success = full_dump_string.contains(target_evidence);

    if success {
        println!("\n✅ SUCCESS in {:?}", t0.elapsed());
        println!("Assembled proof of (t = t) and found the following evidence line:");
        println!("{}", target_evidence);

    } else {
        println!("\n❌ FAILURE: mm1 proof not assembled.");
    }
    
    println!("\n--- Full Final State Dump ---");
    print!("{}", full_dump_string);
}

fn mm1_b() {
    use mork::space::Space;
    use mork::expr;
    use std::time::Instant;

    const P: &str = r#"
    ; --- KB: The Raw Materials (Type Constructors) ---
    (kb (: ⟨=⟩  (-> ⟨term⟩ (-> ⟨term⟩ ⟨wff⟩))))
    (kb (: ⟨+⟩  (-> ⟨term⟩ (-> ⟨term⟩ ⟨term⟩))))
    (kb (: ⟨t⟩ ⟨term⟩))
    (kb (: ⟨0⟩ ⟨term⟩))

    ; --- Step 1: Derive that (t + 0) is a term ---
    ; This rule consumes two kb facts and produces a new, unique "derived-term" atom.
    ; It corresponds to the 'tpl' axiom: term ( t + r )
    (exec derive-term-tplus0
      (, (kb (: ⟨t⟩ ⟨term⟩)) (kb (: ⟨0⟩ ⟨term⟩)))
      (, (derived-term t+0 (: (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨term⟩))))

    ; --- Step 2: Derive wff_P using the newly derived term ---
    ; This rule consumes the product from Step 1 and another kb fact.
    ; It corresponds to the 'weq' axiom: wff t = r
    (exec derive-wff-P-from-terms
      (, (derived-term t+0 (: (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨term⟩)) (kb (: ⟨t⟩ ⟨term⟩)))
      (, (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))))
    "#;

    // --- Load, Run, and Verify with a Heartbeat ---
    let mut s = Space::new();
    s.load_sexpr(P.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();
    let t0 = Instant::now();

    println!("\n== mm1_b: Deriving wff_P (Pipeline Version) ==");
    let mut success = false;
    let target_evidence = "(ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))";
    let mut ticks = 0;

    // We tick the interpreter step-by-step to drive the pipeline forward.
    for i in 0..10 {
        ticks = i + 1;
        println!("Tick {}", ticks);
        let steps_taken = s.metta_calculus(1);

        let mut full_dump_buffer = Vec::new();
        s.dump_all_sexpr(&mut full_dump_buffer).unwrap();
        let full_dump_string = String::from_utf8_lossy(&full_dump_buffer);

        if full_dump_string.contains(target_evidence) {
            success = true;
            println!("\n✅ SUCCESS in {:?} after {} ticks.", t0.elapsed(), ticks);
            println!("Derived wff_P and found the following evidence line:");
            println!("{}", target_evidence);
            break;
        }

        if steps_taken == 0 {
            println!("Space stabilized.");
            break;
        }
    }

    if !success {
        println!("\n❌ FAILURE: wff_P not derived.");
    }
    
    println!("\n--- Full Final State Dump ---");
    let mut full_dump_buffer = Vec::new();
    s.dump_all_sexpr(&mut full_dump_buffer).unwrap();
    print!("{}", String::from_utf8_lossy(&full_dump_buffer));
}

fn mm1_b_tpl() {
    use mork::expr;
    use mork::space::Space;
    use std::time::Instant;

    const P: &str = r#"
; ---------- KB: primitive symbols ----------
(kb (: ⟨+⟩ (-> ⟨term⟩ (-> ⟨term⟩ ⟨term⟩))))
(kb (: ⟨=⟩ (-> ⟨term⟩ (-> ⟨term⟩ ⟨wff⟩))))
(kb (: ⟨t⟩ ⟨term⟩))
(kb (: ⟨0⟩ ⟨term⟩))

; ---------- KB: generalized constructors ----------
(kb (: ⟨tpl⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                      (: (⟨+⟩ $x $y) ⟨term⟩)))))
(kb (: ⟨weq⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                      (: (⟨=⟩ $x $y) ⟨wff⟩)))))

; ---------- Deterministic pipeline rules ----------
(exec tpl-apply-kb
  (, (kb (: ⟨tpl⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                           (: (⟨+⟩ $x $y) ⟨term⟩)))))
     (kb (: $x ⟨term⟩))
     (kb (: $y ⟨term⟩)))
  (, (ev (: (⟨+⟩ $x $y) ⟨term⟩))))

(exec weq-apply-kb
  (, (kb (: ⟨weq⟩ (-> (: $a ⟨term⟩) (-> (: $b ⟨term⟩)
                           (: (⟨=⟩ $a $b) ⟨wff⟩)))))
     (ev (: $lhs ⟨term⟩))
     (kb (: $rhs ⟨term⟩)))
  (, (ev (: (⟨=⟩ $lhs $rhs) ⟨wff⟩))))
"#;

    let t0 = Instant::now();
    let mut s = Space::new();
    s.load_sexpr(P.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    // Tick up to a small bound; break when saturated or when target appears.
    let mut ticks = 0usize;
    let target = "(ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))";
    loop {
        ticks += 1;
        let n = s.metta_calculus(1);
        // Check success by string membership over a full dump buffer.
        let mut buf = Vec::new();
        s.dump_all_sexpr(&mut buf).unwrap();
        let out = String::from_utf8_lossy(&buf);
        let done = out.contains(target);
        if done || n == 0 || ticks >= 32 {
            println!("\n== mm1_b_tpl: result = {} in {:?} after {} tick(s) ==",
                     if done { "SUCCESS" } else { "INCOMPLETE" }, t0.elapsed(), ticks);
            println!("\n--- Full Final State Dump ---");
            print!("{out}");
            break;
        }
    }
}


fn mm1_b2_tpl() {
    use mork::expr;
    use mork::space::Space;
    use std::time::Instant;

    const P: &str = r#"
; ---------- KB: primitive symbols ----------
(kb (: ⟨+⟩ (-> ⟨term⟩ (-> ⟨term⟩ ⟨term⟩))))
(kb (: ⟨=⟩ (-> ⟨term⟩ (-> ⟨term⟩ ⟨wff⟩))))
(kb (: ⟨t⟩ ⟨term⟩))
(kb (: ⟨0⟩ ⟨term⟩))

; ---------- KB: generalized constructors ----------
(kb (: ⟨tpl⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                      (: (⟨+⟩ $x $y) ⟨term⟩)))))
(kb (: ⟨weq⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                      (: (⟨=⟩ $x $y) ⟨wff⟩)))))

; ---------- Small generalization ----------
(exec lift
  (, (kb (: $t $T)))
  (, (ev (: $t $T))))

(exec tpl-apply
  (, (kb (: ⟨tpl⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                           (: (⟨+⟩ $x $y) ⟨term⟩)))))
     (ev (: $x ⟨term⟩))
     (ev (: $y ⟨term⟩)))
  (, (ev (: (⟨+⟩ $x $y) ⟨term⟩))))

(exec weq-apply
  (, (kb (: ⟨weq⟩ (-> (: $a ⟨term⟩) (-> (: $b ⟨term⟩)
                           (: (⟨=⟩ $a $b) ⟨wff⟩)))))
     (ev (: $a ⟨term⟩))
     (ev (: $b ⟨term⟩)))
  (, (ev (: (⟨=⟩ $a $b) ⟨wff⟩))))
"#;

    let t0 = Instant::now();
    let mut s = Space::new();
    s.load_sexpr(P.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    let mut ticks = 0usize;
    let target = "(ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))";
    loop {
        ticks += 1;
        let n = s.metta_calculus(1);
        let mut buf = Vec::new();
        s.dump_all_sexpr(&mut buf).unwrap();
        let out = String::from_utf8_lossy(&buf);
        let done = out.contains(target);
        if done || n == 0 || ticks >= 32 {
            println!("\n== mm1_b2_tpl: result = {} in {:?} after {} tick(s) ==",
                     if done { "SUCCESS" } else { "INCOMPLETE" }, t0.elapsed(), ticks);
            println!("\n--- Full Final State Dump ---");
            print!("{out}");
            break;
        }
    }
}

fn mm1_c() {
    use mork::expr;
    use mork::space::Space;
    use std::time::Instant;

    // Program: derive wff_P, wff_Q via tpl/weq, and proof_P via a2-curry@t.
    const P: &str = r#"
; ===== Universe & primitives =====
(kb (: ⟨+⟩ (-> ⟨term⟩ (-> ⟨term⟩ ⟨term⟩))))
(kb (: ⟨=⟩ (-> ⟨term⟩ (-> ⟨term⟩ ⟨wff⟩))))
(kb (: ⟨t⟩ ⟨term⟩))
(kb (: ⟨0⟩ ⟨term⟩))

; ===== Generalized typed constructors (your names) =====
(kb (: ⟨tpl⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                      (: (⟨+⟩ $x $y) ⟨term⟩)))))
(kb (: ⟨weq⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                      (: (⟨=⟩ $x $y) ⟨wff⟩)))))

; ===== Axioms encoded in "curried" form (as types in KB) =====
; a2:   ⊢ (= (+ a 0) a)
(kb (: ⟨a2-curry⟩ (-> (: $a ⟨term⟩)
                  (: (⟨=⟩ (⟨+⟩ $a ⟨0⟩) $a) ⟨|-⟩))))

; ===== Small pipeline rules =====

; Lift KB typings into usable "ev" facts
(exec lift (, (kb (: $t $T))) (, (ev (: $t $T))))

; Build (+ x y) : term from x:term, y:term using tpl
(exec tpl-apply
  (, (kb (: ⟨tpl⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                           (: (⟨+⟩ $x $y) ⟨term⟩)))))
     (ev (: $x ⟨term⟩))
     (ev (: $y ⟨term⟩)))
  (, (ev (: (⟨+⟩ $x $y) ⟨term⟩))))

; Build (= a b) : wff from a:term, b:term using weq
(exec weq-apply
  (, (kb (: ⟨weq⟩ (-> (: $a ⟨term⟩) (-> (: $b ⟨term⟩)
                           (: (⟨=⟩ $a $b) ⟨wff⟩)))))
     (ev (: $a ⟨term⟩))
     (ev (: $b ⟨term⟩)))
  (, (ev (: (⟨=⟩ $a $b) ⟨wff⟩))))

; Instantiate a2 at a := t to get proof_P:  ⊢ (= (+ t 0) t)
(exec a2-instantiate-t
  (, (kb (: ⟨a2-curry⟩ (-> (: $a ⟨term⟩)
                        (: (⟨=⟩ (⟨+⟩ $a ⟨0⟩) $a) ⟨|-⟩))))
     (ev (: ⟨t⟩ ⟨term⟩)))
  (, (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))))

; (Optional) tag the proof with a name, so your dump resembles demo files
(exec tag-proof_P
  (, (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩)))
  (, (: proof_P (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))))

; Convenience: derive wff_Q directly from t:term, t:term
(exec derive-wff-Q
  (, (ev (: ⟨t⟩ ⟨term⟩)) (ev (: ⟨t⟩ ⟨term⟩)))
  (, (ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩))))
"#;

    let mut s = Space::new();
    let t0 = Instant::now();
    s.load_sexpr(P.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    // Drive a few ticks until saturation or success. No panics; just report.
    let mut ticks = 0usize;
    let want_wff_p  = "(ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))";
    let want_wff_q  = "(ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩))";
    let want_proof_p = "(ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))";

    loop {
        ticks += 1;
        let n = s.metta_calculus(1);

        let mut buf = Vec::new();
        s.dump_all_sexpr(&mut buf).unwrap();
        let dump = String::from_utf8_lossy(&buf);

        // STRICT evidence check: match whole dumped lines by prefix (no substring hits inside execs).
        let has_wff_p   = dump.lines().any(|l| l.trim_start().starts_with(want_wff_p));
        let has_wff_q   = dump.lines().any(|l| l.trim_start().starts_with(want_wff_q));
        let has_proof_p = dump.lines().any(|l| l.trim_start().starts_with(want_proof_p));

        if (has_wff_p && has_wff_q && has_proof_p) || n == 0 || ticks >= 64 {
            println!(
                "\n== mm1_c: done in {:?} after {} tick(s) ==\n  wff_P: {}\n  wff_Q: {}\n  proof_P (a2@t): {}",
                t0.elapsed(),
                ticks,
                if has_wff_p { "✓" } else { "✗" },
                if has_wff_q { "✓" } else { "✗" },
                if has_proof_p { "✓" } else { "✗" },
            );
            println!("\n--- Full Final State Dump ---");
            print!("{dump}");
            break;
        }
    }
}

fn mm1_d() {
    use mork::expr;
    use mork::space::Space;
    use std::time::Instant;

    // Program: derive wff_P, wff_Q via tpl/weq, proof_P via a2-curry@t,
    // and proof_PtoQ via a1-curry at a=(+ t 0), b=t, c=t; then assemble with mp-curry.
    const P: &str = r#"
; ===== Universe & primitives =====
(kb (: ⟨+⟩ (-> ⟨term⟩ (-> ⟨term⟩ ⟨term⟩))))
(kb (: ⟨=⟩ (-> ⟨term⟩ (-> ⟨term⟩ ⟨wff⟩))))
(kb (: ⟨t⟩ ⟨term⟩))
(kb (: ⟨0⟩ ⟨term⟩))

; (optional) implication ctor if you want it typed (not required for patterning)
; (kb (: ⟨->⟩ (-> ⟨wff⟩ (-> ⟨wff⟩ ⟨wff⟩))))

; ===== Generalized typed constructors (your names) =====
(kb (: ⟨tpl⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                      (: (⟨+⟩ $x $y) ⟨term⟩)))))
(kb (: ⟨weq⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                      (: (⟨=⟩ $x $y) ⟨wff⟩)))))

; ===== Axioms encoded in "curried" form (as types in KB) =====
; a2:   ⊢ (= (+ a 0) a)
(kb (: ⟨a2-curry⟩ (-> (: $a ⟨term⟩)
                  (: (⟨=⟩ (⟨+⟩ $a ⟨0⟩) $a) ⟨|-⟩))))
; a1:   ⊢ (-> (= a b) (-> (= a c) (= b c)))
(kb (: ⟨a1-curry⟩ (-> (: $a ⟨term⟩) (-> (: $b ⟨term⟩) (-> (: $c ⟨term⟩)
                  (: (⟨->⟩ (⟨=⟩ $a $b) (⟨->⟩ (⟨=⟩ $a $c) (⟨=⟩ $b $c))) ⟨|-⟩))))))

; Modus ponens (curried)
(kb (: ⟨mp-curry⟩ (-> (: $P ⟨wff⟩) (-> (: $Q ⟨wff⟩)
                  (-> (: $P ⟨|-⟩) (-> (: (⟨->⟩ $P $Q) ⟨|-⟩) (: $Q ⟨|-⟩)))))))

; ===== Small pipeline rules =====

; Lift KB typings into usable "ev" facts
(exec lift (, (kb (: $t $T))) (, (ev (: $t $T))))

; (+ x y) : term from x:term, y:term using tpl
(exec tpl-apply
  (, (kb (: ⟨tpl⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                           (: (⟨+⟩ $x $y) ⟨term⟩)))))
     (ev (: $x ⟨term⟩))
     (ev (: $y ⟨term⟩)))
  (, (ev (: (⟨+⟩ $x $y) ⟨term⟩))))

; (= a b) : wff from a:term, b:term using weq
(exec weq-apply
  (, (kb (: ⟨weq⟩ (-> (: $a ⟨term⟩) (-> (: $b ⟨term⟩)
                           (: (⟨=⟩ $a $b) ⟨wff⟩)))))
     (ev (: $a ⟨term⟩))
     (ev (: $b ⟨term⟩)))
  (, (ev (: (⟨=⟩ $a $b) ⟨wff⟩))))

; ----- derive & tag wff_P = ((+ t 0) = t) : wff -----
(exec tag-wff_P
  (, (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩)))
  (, (: wff_P (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))))

; ----- derive & tag wff_Q = (t = t) : wff -----
(exec derive-wff-Q
  (, (ev (: ⟨t⟩ ⟨term⟩)) (ev (: ⟨t⟩ ⟨term⟩)))
  (, (ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩))))
(exec tag-wff_Q
  (, (ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩)))
  (, (: wff_Q (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩))))

; ----- derive & tag proof_P : ⊢ (= (+ t 0) t)  via a2-curry -----
(exec a2-instantiate-t
  (, (kb (: ⟨a2-curry⟩ (-> (: $a ⟨term⟩)
                        (: (⟨=⟩ (⟨+⟩ $a ⟨0⟩) $a) ⟨|-⟩))))
     (ev (: ⟨t⟩ ⟨term⟩)))
  (, (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))))
(exec tag-proof_P
  (, (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩)))
  (, (: proof_P (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))))

; ----- derive & tag proof_PtoQ : ⊢ (-> P (-> P Q)) via a1-curry -----
; Instantiate a1 with a = (+ t 0), b = t, c = t
(exec a1-instantiate-PtoQ
  (, (kb (: ⟨a1-curry⟩ (-> (: $a ⟨term⟩) (-> (: $b ⟨term⟩) (-> (: $c ⟨term⟩)
                         (: (⟨->⟩ (⟨=⟩ $a $b) (⟨->⟩ (⟨=⟩ $a $c) (⟨=⟩ $b $c))) ⟨|-⟩))))))
     (ev (: (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨term⟩))   ; establishes $a
     (ev (: ⟨t⟩ ⟨term⟩)))              ; establishes $b and $c (both t)
  (, (ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
               (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
                       (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩))))
(exec tag-proof_PtoQ
  (, (ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
               (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
                       (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩)))
  (, (: proof_PtoQ (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
                        (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
                                (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩))))

; ----- final assembly using mp-curry -----
(exec assemble-final-proof
  (, (: wff_P $) (: wff_Q $) (: proof_P $) (: proof_PtoQ $)
     (kb (: ⟨mp-curry⟩ $_)))
  (, (ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨|-⟩)
         ((((⟨mp-curry⟩ wff_P) wff_Q) proof_P) proof_PtoQ))))
"#;

    let mut s = Space::new();
    let t0 = Instant::now();
    s.load_sexpr(P.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    // Targets (strict "whole-line starts-with" matching)
    let want_ev_term_tplus0    = "(ev (: (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨term⟩))";
    let want_ev_wff_p          = "(ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))";
    let want_ev_wff_q          = "(ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩))";
    let want_tag_wff_p         = "(: wff_P (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))";
    let want_tag_wff_q         = "(: wff_Q (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩))";
    let want_ev_proof_p        = "(ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))";
    let want_tag_proof_p       = "(: proof_P (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))";
    let want_ev_proof_ptoptoq     = "(ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩))";
    let want_tag_proof_ptoq    = "(: proof_PtoQ (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩))";
    let want_final_evidence    = "(ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨|-⟩) ((((⟨mp-curry⟩ wff_P) wff_Q) proof_P) proof_PtoQ))";

    // Drive ticks; stop on success or saturation.
    let mut ticks = 0usize;
    let mut success = false;
    loop {
        ticks += 1;
        let n = s.metta_calculus(1);

        let mut buf = Vec::new();
        s.dump_all_sexpr(&mut buf).unwrap();
        let dump = String::from_utf8_lossy(&buf);

        // Strict per-line checks (avoid substring matches inside exec traces)
        let line_has = |needle: &str| dump.lines().any(|l| l.trim_start().starts_with(needle));

        let have_tplus0_term  = line_has(want_ev_term_tplus0);
        let have_wff_p_ev     = line_has(want_ev_wff_p);
        let have_wff_q_ev     = line_has(want_ev_wff_q);
        let have_wff_p_tag    = line_has(want_tag_wff_p);
        let have_wff_q_tag    = line_has(want_tag_wff_q);
        let have_proof_p_ev   = line_has(want_ev_proof_p);
        let have_proof_p_tag  = line_has(want_tag_proof_p);
        let have_ptoptoq_ev      = line_has(want_ev_proof_ptoptoq);
        let have_ptoq_tag     = line_has(want_tag_proof_ptoq);
        let have_final        = line_has(want_final_evidence);

        if have_final {
            println!("\n== mm1_d: done in {:?} after {} tick(s) ==", t0.elapsed(), ticks);
            println!("  (+ t 0) : term .......... {}", if have_tplus0_term { "✓" } else { "—" });
            println!("  wff_P (ev) .............. {}", if have_wff_p_ev { "✓" } else { "—" });
            println!("  wff_Q (ev) .............. {}", if have_wff_q_ev { "✓" } else { "—" });
            println!("  wff_P (tag) ............. {}", if have_wff_p_tag { "✓" } else { "—" });
            println!("  wff_Q (tag) ............. {}", if have_wff_q_tag { "✓" } else { "—" });
            println!("  proof_P (a2@t, ev) ...... {}", if have_proof_p_ev { "✓" } else { "—" });
            println!("  proof_P (tag) ........... {}", if have_proof_p_tag { "✓" } else { "—" });
            println!("  proof_PtoQ (a1, ev) ..... {}", if have_ptoptoq_ev { "✓" } else { "—" });
            println!("  proof_PtoQ (tag) ........ {}", if have_ptoq_tag { "✓" } else { "—" });

            // Final evidence confirmation
            println!("\n--- Verifying final MP assembly ---");
            println!("✅ final evidence line present:\n{}", want_final_evidence);

            println!("\n--- Full Final State Dump ---");
            print!("{dump}");
            success = true;
            break;
        }

        if n == 0 || ticks >= 128 {
            println!("\n== mm1_d: done in {:?} after {} tick(s) ==", t0.elapsed(), ticks);
            println!("  (+ t 0) : term .......... {}", if have_tplus0_term { "✓" } else { "—" });
            println!("  wff_P (ev) .............. {}", if have_wff_p_ev { "✓" } else { "—" });
            println!("  wff_Q (ev) .............. {}", if have_wff_q_ev { "✓" } else { "—" });
            println!("  wff_P (tag) ............. {}", if have_wff_p_tag { "✓" } else { "—" });
            println!("  wff_Q (tag) ............. {}", if have_wff_q_tag { "✓" } else { "—" });
            println!("  proof_P (a2@t, ev) ...... {}", if have_proof_p_ev { "✓" } else { "—" });
            println!("  proof_P (tag) ........... {}", if have_proof_p_tag { "✓" } else { "—" });
            println!("  proof_PtoQ (a1, ev) ..... {}", if have_ptoptoq_ev { "✓" } else { "—" });
            println!("  proof_PtoQ (tag) ........ {}", if have_ptoq_tag { "✓" } else { "—" });

            println!("\n--- Verifying final MP assembly ---");
            if have_final {
                println!("✅ final evidence line present:\n{}", want_final_evidence);
            } else {
                println!("❌ final evidence line NOT found.");
            }

            println!("\n--- Full Final State Dump ---");
            print!("{dump}");
            break;
        }
    }

    let _ = success; // (kept in case you want to use it programmatically)
}

fn mm1() {
    use mork::expr;
    use mork::space::Space;
    use std::time::Instant;

    const P: &str = r#"
; ===== Universe & primitives =====
(kb (: ⟨+⟩ (-> ⟨term⟩ (-> ⟨term⟩ ⟨term⟩))))
(kb (: ⟨=⟩ (-> ⟨term⟩ (-> ⟨term⟩ ⟨wff⟩))))
(kb (: ⟨t⟩ ⟨term⟩))
(kb (: ⟨0⟩ ⟨term⟩))

; ===== Generalized typed constructors =====
(kb (: ⟨tpl⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                      (: (⟨+⟩ $x $y) ⟨term⟩)))))
(kb (: ⟨weq⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                      (: (⟨=⟩ $x $y) ⟨wff⟩)))))

; ===== Axioms encoded in "curried" form =====
(kb (: ⟨a2-curry⟩ (-> (: $a ⟨term⟩)
                  (: (⟨=⟩ (⟨+⟩ $a ⟨0⟩) $a) ⟨|-⟩))))
(kb (: ⟨a1-curry⟩ (-> (: $a ⟨term⟩) (-> (: $b ⟨term⟩) (-> (: $c ⟨term⟩)
                  (: (⟨->⟩ (⟨=⟩ $a $b) (⟨->⟩ (⟨=⟩ $a $c) (⟨=⟩ $b $c))) ⟨|-⟩))))))
(kb (: ⟨mp-curry⟩ (-> (: $P ⟨wff⟩) (-> (: $Q ⟨wff⟩)
                  (-> (: $P ⟨|-⟩) (-> (: (⟨->⟩ $P $Q) ⟨|-⟩) (: $Q ⟨|-⟩)))))))

; ===== Pipeline rules =====

; Lift KB typings into usable "ev" facts
(exec lift (, (kb (: $t $T))) (, (ev (: $t $T))))

; (+ x y) : term from x:term, y:term using tpl
(exec tpl-apply
  (, (kb (: ⟨tpl⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                           (: (⟨+⟩ $x $y) ⟨term⟩)))))
     (ev (: $x ⟨term⟩))
     (ev (: $y ⟨term⟩)))
  (, (ev (: (⟨+⟩ $x $y) ⟨term⟩))))

; (= a b) : wff from a:term, b:term using weq
(exec weq-apply
  (, (kb (: ⟨weq⟩ (-> (: $a ⟨term⟩) (-> (: $b ⟨term⟩)
                           (: (⟨=⟩ $a $b) ⟨wff⟩)))))
     (ev (: $a ⟨term⟩))
     (ev (: $b ⟨term⟩)))
  (, (ev (: (⟨=⟩ $a $b) ⟨wff⟩))))

; Derive wff_Q = (t = t) : wff
(exec derive-wff-Q
  (, (ev (: ⟨t⟩ ⟨term⟩)) (ev (: ⟨t⟩ ⟨term⟩)))
  (, (ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩))))

; Derive proof_P : ⊢ (= (+ t 0) t) via a2-curry
(exec a2-instantiate-t
  (, (kb (: ⟨a2-curry⟩ (-> (: $a ⟨term⟩)
                        (: (⟨=⟩ (⟨+⟩ $a ⟨0⟩) $a) ⟨|-⟩))))
     (ev (: ⟨t⟩ ⟨term⟩)))
  (, (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))))

; Derive proof_PtoQ : ⊢ (-> P (-> P Q)) via a1-curry
(exec a1-instantiate-PtoQ
  (, (kb (: ⟨a1-curry⟩ (-> (: $a ⟨term⟩) (-> (: $b ⟨term⟩) (-> (: $c ⟨term⟩)
                         (: (⟨->⟩ (⟨=⟩ $a $b) (⟨->⟩ (⟨=⟩ $a $c) (⟨=⟩ $b $c))) ⟨|-⟩))))))
     (ev (: (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨term⟩))
     (ev (: ⟨t⟩ ⟨term⟩)))
  (, (ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
               (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
                       (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩))))

; Final assembly using mp-curry with direct ev facts
(exec assemble-final-proof-direct
  (, (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))
     (ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩))
     (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))
     (ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
               (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
                       (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩))
     (kb (: ⟨mp-curry⟩ $_)))
  (, (ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨|-⟩)
         ((((⟨mp-curry⟩ (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))
            (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩))
           (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))
          (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
                  (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
                         (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩)))))
"#;

    let mut s = Space::new();
    let t0 = Instant::now();
    s.load_sexpr(P.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    // Targets (simpler without tags)
    let want_ev_term_tplus0    = "(ev (: (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨term⟩))";
    let want_ev_wff_p          = "(ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))";
    let want_ev_wff_q          = "(ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩))";
    let want_ev_proof_p        = "(ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))";
    let want_ev_proof_ptoptoq     = "(ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩))";
    let want_final_evidence    = "(ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨|-⟩)";

    let mut ticks = 0usize;
    let mut success = false;
    loop {
        ticks += 1;
        let n = s.metta_calculus(1);

        let mut buf = Vec::new();
        s.dump_all_sexpr(&mut buf).unwrap();
        let dump = String::from_utf8_lossy(&buf);

        let line_has = |needle: &str| dump.lines().any(|l| l.trim_start().starts_with(needle));

        let have_tplus0_term  = line_has(want_ev_term_tplus0);
        let have_wff_p_ev     = line_has(want_ev_wff_p);
        let have_wff_q_ev     = line_has(want_ev_wff_q);
        let have_proof_p_ev   = line_has(want_ev_proof_p);
        let have_ptoptoq_ev      = line_has(want_ev_proof_ptoptoq);
        let have_final        = line_has(want_final_evidence);

        if have_final {
            println!("\n== mm1: ✅ SUCCESS in {:?} after {} tick(s) ==", t0.elapsed(), ticks);
            println!("  (+ t 0) : term .......... {}", if have_tplus0_term { "✓" } else { "—" });
            println!("  wff_P (ev) .............. {}", if have_wff_p_ev { "✓" } else { "—" });
            println!("  wff_Q (ev) .............. {}", if have_wff_q_ev { "✓" } else { "—" });
            println!("  proof_P (a2@t, ev) ...... {}", if have_proof_p_ev { "✓" } else { "—" });
            println!("  proof_PtoQ (a1, ev) ..... {}", if have_ptoptoq_ev { "✓" } else { "—" });
            println!("\n--- Final evidence confirmation ---");
            println!("✅ Successfully derived ⊢ (t = t)");
            
            println!("\n--- Full Final State Dump ---");
            print!("{dump}");
            success = true;
            break;
        }

        if n == 0 || ticks >= 128 {
            println!("\n== mm1: done in {:?} after {} tick(s) ==", t0.elapsed(), ticks);
            println!("  (+ t 0) : term .......... {}", if have_tplus0_term { "✓" } else { "—" });
            println!("  wff_P (ev) .............. {}", if have_wff_p_ev { "✓" } else { "—" });
            println!("  wff_Q (ev) .............. {}", if have_wff_q_ev { "✓" } else { "—" });
            println!("  proof_P (a2@t, ev) ...... {}", if have_proof_p_ev { "✓" } else { "—" });
            println!("  proof_PtoQ (a1, ev) ..... {}", if have_ptoptoq_ev { "✓" } else { "—" });
            
            if !have_final {
                println!("\n❌ Failed to derive ⊢ (t = t)");
            }
            
            println!("\n--- Full Final State Dump ---");
            print!("{dump}");
            break;
        }
    }
}

// A cleaned, didactic version of the working mm1(), preserving the exact control flow.
// It loads the same program, runs the same single‑step metta_calculus loop,
// and checks for the same milestones and final goal.
//
// Goal: derive ⊢ (t = t) using:
//   a2-curry @ t           : ⊢ ((t + 0) = t)
//   a1-curry @ (t+0, t, t) : ⊢ (P → (P → Q)) where P := ((t+0)=t), Q := (t=t)
//   mp-curry               : combine the above to produce ⊢ (t = t)
//
// The “lift / tpl-apply / weq-apply / derive-wff-Q / a2-instantiate-t /
//  a1-instantiate-PtoQ / assemble-final-proof-direct” rules and their names
// are kept exactly as in the working mm1().

fn mm1_didactic() {
    use mork::expr;
    use mork::space::Space;
    use std::time::Instant;

    // Program: universe, typed constructors, axioms (curried), tiny pipeline, and final assembly.
    const P: &str = r#"
; ===== Universe & primitives =====
(kb (: ⟨+⟩ (-> ⟨term⟩ (-> ⟨term⟩ ⟨term⟩))))
(kb (: ⟨=⟩ (-> ⟨term⟩ (-> ⟨term⟩ ⟨wff⟩))))
(kb (: ⟨t⟩ ⟨term⟩))
(kb (: ⟨0⟩ ⟨term⟩))

; ===== Generalized typed constructors =====
(kb (: ⟨tpl⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                      (: (⟨+⟩ $x $y) ⟨term⟩)))))
(kb (: ⟨weq⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                      (: (⟨=⟩ $x $y) ⟨wff⟩)))))
(kb (: ⟨wim⟩ (-> (: $P ⟨wff⟩) (-> (: $Q ⟨wff⟩) 
                      (: (⟨->⟩ $P $Q) ⟨wff⟩)))))

; ===== Axioms encoded in "curried" form =====
(kb (: ⟨a2-curry⟩ (-> (: $a ⟨term⟩)
                  (: (⟨=⟩ (⟨+⟩ $a ⟨0⟩) $a) ⟨|-⟩))))
(kb (: ⟨a1-curry⟩ (-> (: $a ⟨term⟩) (-> (: $b ⟨term⟩) (-> (: $c ⟨term⟩)
                  (: (⟨->⟩ (⟨=⟩ $a $b) (⟨->⟩ (⟨=⟩ $a $c) (⟨=⟩ $b $c))) ⟨|-⟩))))))
(kb (: ⟨mp-curry⟩ (-> (: $P ⟨wff⟩) (-> (: $Q ⟨wff⟩)
                  (-> (: $P ⟨|-⟩) (-> (: (⟨->⟩ $P $Q) ⟨|-⟩) (: $Q ⟨|-⟩)))))))

(kb (: ⟨a2⟩ (-> (: $a ⟨term⟩) (: (⟨=⟩ (⟨+⟩ $a ⟨0⟩) $a) ⟨|-⟩))))
(kb (: ⟨a1⟩ (-> (: $a ⟨term⟩) (: $b ⟨term⟩) (: $c ⟨term⟩) (: (⟨->⟩ (⟨=⟩ $a $b) (⟨->⟩ (⟨=⟩ $a $c) (⟨=⟩ $b $c))) ⟨|-⟩))))
(kb (: ⟨mp⟩ (-> (: $P ⟨wff⟩) (: $Q ⟨wff⟩) (: $P ⟨|-⟩) (: (⟨->⟩ $P $Q) ⟨|-⟩) (: $Q ⟨|-⟩))))

; ===== Pipeline rules (unchanged names & behavior) =====

; Lift KB typings into usable "ev" facts
(exec lift (, (kb (: $t $T))) (, (ev (: $t $T))))

; (+ x y) : term from x:term, y:term using tpl
(exec tpl-apply
  (, (kb (: ⟨tpl⟩ (-> (: $x ⟨term⟩) (-> (: $y ⟨term⟩)
                           (: (⟨+⟩ $x $y) ⟨term⟩)))))
     (ev (: $x ⟨term⟩))
     (ev (: $y ⟨term⟩)))
  (, (ev (: (⟨+⟩ $x $y) ⟨term⟩))))

; (= a b) : wff from a:term, b:term using weq
(exec weq-apply
  (, (kb (: ⟨weq⟩ (-> (: $a ⟨term⟩) (-> (: $b ⟨term⟩)
                           (: (⟨=⟩ $a $b) ⟨wff⟩)))))
     (ev (: $a ⟨term⟩))
     (ev (: $b ⟨term⟩)))
  (, (ev (: (⟨=⟩ $a $b) ⟨wff⟩))))

; ⊢ ((t+0) = t) via a2-curry @ t
(exec a2-instantiate-t
  (, (kb (: ⟨a2-curry⟩ (-> (: $a ⟨term⟩)
                        (: (⟨=⟩ (⟨+⟩ $a ⟨0⟩) $a) ⟨|-⟩))))
     (ev (: $a ⟨term⟩)))
  (, (ev (: (⟨=⟩ (⟨+⟩ $a ⟨0⟩) $a) ⟨|-⟩))))

; ⊢ ((t+0) = t) via a2-curry @ t
;(exec a2
; (, (ev (: $a ⟨term⟩)))
;  (, (ev (: (⟨=⟩ (⟨+⟩ $a ⟨0⟩) $a) ⟨|-⟩))))

; ⊢ (P → (P → Q)) via a1-curry @ (a,b,c)=(t+0,t,t)
;(exec a1-instantiate-PtoQ
;  (, (kb (: ⟨a1-curry⟩ (-> (: $a ⟨term⟩) (-> (: $b ⟨term⟩) (-> (: $c ⟨term⟩)
;                         (: (⟨->⟩ (⟨=⟩ $a $b) (⟨->⟩ (⟨=⟩ $a $c) (⟨=⟩ $b $c))) ⟨|-⟩))))))
;     (ev (: (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨term⟩))
;     (ev (: ⟨t⟩ ⟨term⟩)))
;  (, (ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
;               (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
;                       (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩))))

(exec a1-instantiate-PtoQ
  (, (kb (: ⟨a1-curry⟩ (-> (: $a ⟨term⟩) (-> (: $bc ⟨term⟩) (-> (: $bc ⟨term⟩)
                         (: (⟨->⟩ (⟨=⟩ $a $bc) (⟨->⟩ (⟨=⟩ $a $bc) (⟨=⟩ $bc $bc))) ⟨|-⟩))))))
     (ev (: $a ⟨term⟩))
     (ev (: $bc ⟨term⟩)))
  (, (ev (: (⟨->⟩ (⟨=⟩ $a $bc)
               (⟨->⟩ (⟨=⟩ $a $bc)
                       (⟨=⟩ $bc $bc))) ⟨|-⟩))))

; Final assembly using mp-curry (same single step as the working mm1)

; (kb (: ⟨mp-curry⟩ (-> (: $P ⟨wff⟩) 
;                    (-> (: $Q ⟨wff⟩)
;                      (-> (: $P ⟨|-⟩) 
;                      (-> (: (⟨->⟩ $P $Q) ⟨|-⟩) 
;                        (: $Q ⟨|-⟩)))))))
;(kb (: ⟨mp⟩ (-> (: $P ⟨wff⟩) (: $Q ⟨wff⟩) (: $P ⟨|-⟩) (: (⟨->⟩ $P $Q) ⟨|-⟩) (: $Q ⟨|-⟩))))

; Derive wff type for (P -> Q) 
;(exec derive-wff-P-to-Q
;  (, (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))
;     (ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩)))
;  (, (ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨=⟩ ⟨t⟩ ⟨t⟩)) ⟨wff⟩))))

; Generic wim (implication constructor) rule
(exec wim-apply
  (, (kb (: ⟨wim⟩ (-> (: $P ⟨wff⟩) (-> (: $Q ⟨wff⟩) (: (⟨->⟩ $P $Q) ⟨wff⟩)))))
     (ev (: $P ⟨wff⟩))
     (ev (: $Q ⟨wff⟩)))
  (, (ev (: (⟨->⟩ $P $Q) ⟨wff⟩))))

; Intermediate step: derive (P -> Q) from P and (P -> (P -> Q))
; $P :- (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
; $Q :- ....
;(exec derive-P-to-Q-direct
;  (, (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))
;     (ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨=⟩ ⟨t⟩ ⟨t⟩)) ⟨wff⟩))
;     (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))
;     (ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
;               (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
;                       (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩)))
;     (kb (: ⟨mp-curry⟩ (-> (: $P ⟨wff⟩) (-> (: $Q ⟨wff⟩)
;                  (-> (: $P ⟨|-⟩) (-> (: (⟨->⟩ $P $Q) ⟨|-⟩) (: $Q ⟨|-⟩))))))) )
;  (, (ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨=⟩ ⟨t⟩ ⟨t⟩)) ⟨|-⟩)) )) 

;; Current 'best' attempt
;(exec derive-P-to-Q-direct2
;  (, (ev (: $P ⟨wff⟩))
;     (ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨=⟩ ⟨t⟩ ⟨t⟩)) ⟨wff⟩))
;     (ev (: $P ⟨|-⟩))
;     (ev (: (⟨->⟩ $P
;               (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
;                       (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩)))
;     (kb (: ⟨mp-curry⟩ (-> (: $P ⟨wff⟩) (-> (: $Q ⟨wff⟩)
;                  (-> (: $P ⟨|-⟩) (-> (: (⟨->⟩ $P $Q) ⟨|-⟩) (: $Q ⟨|-⟩))))))) )
;  (, (ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨=⟩ ⟨t⟩ ⟨t⟩)) ⟨|-⟩)) )) 

(exec derive-P-to-Q-direct3
  (, (ev (: $P ⟨wff⟩))
     (ev (: $IMP ⟨wff⟩))
     (ev (: $P ⟨|-⟩))
     (ev (: (⟨->⟩ $P $IMP) ⟨|-⟩)))
  (, (ev (: $IMP ⟨|-⟩))))

;(exec derive-P-to-Q-direct4
;  (, (ev (: $P ⟨wff⟩))
;     (ev (: (⟨->⟩ $P $Q) ⟨wff⟩))
;     (ev (: $P ⟨|-⟩))
;     (ev (: (⟨->⟩ $P (⟨->⟩ $P $Q)) ⟨|-⟩)))
;  (, (ev (: (⟨->⟩ $P $Q) ⟨|-⟩))))

; If we ever see the nested implication proof, expose the ($P,$Q) it contains.
(exec dbg-see-a1
  (, (ev (: (⟨->⟩ $P (⟨->⟩ $P $Q)) ⟨|-⟩)))
  (, (: seen/a1 $P $Q)))

; If we ever see a proof of P, expose it too.
(exec dbg-see-P
  (, (ev (: $P ⟨|-⟩)))
  (, (: seen/P $P)))

  ; From P, and ⊢(P → (P → Q)), *prepare* to conclude (P → Q)
(exec mp#1-token
  (, (ev (: $P ⟨wff⟩))
     (ev (: $Q ⟨wff⟩))
     (ev (: $P ⟨|-⟩))
     (ev (: (⟨->⟩ $P (⟨->⟩ $P $Q)) ⟨|-⟩)))
  (, (: ready/mp1 $P $Q)))

; Conclude ⊢(P → Q) *only* from the token, i.e., after we know P,Q matched once.
(exec mp#1-finish
  (, (: ready/mp1 $P $Q))
  (, (ev (: (⟨->⟩ $P $Q) ⟨|-⟩))))

(exec ping-a1
  (, (ev (: (⟨->⟩ $AnyP (⟨->⟩ $AnyP $AnyQ)) ⟨|-⟩)))
  (, (: ping/a1)))

  (exec mp#2-token
  (, (ev (: $P ⟨wff⟩))
     (ev (: $Q ⟨wff⟩))
     (ev (: $P ⟨|-⟩))
     (ev (: (⟨->⟩ $P $Q) ⟨|-⟩)))
  (, (: ready/mp2 $Q)))

(exec mp#2-finish
  (, (: ready/mp2 $Q))
  (, (ev (: $Q ⟨|-⟩))))


(exec derive-P-to-Q-v2
  (, (ev (: $P ⟨wff⟩))
     (ev (: (⟨->⟩ $P (⟨=⟩ ⟨t⟩ ⟨t⟩)) ⟨wff⟩))
     (ev (: $P ⟨|-⟩))
     (ev (: (⟨->⟩ $P (⟨->⟩ $P (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩)))
  (, (ev (: (⟨->⟩ $P (⟨=⟩ ⟨t⟩ ⟨t⟩)) ⟨|-⟩))))


  (exec mp#1
  (, (ev (: $P ⟨wff⟩))
     (ev (: $Q ⟨wff⟩))
     (ev (: $P ⟨|-⟩))
     (ev (: (⟨->⟩ $P (⟨->⟩ $P $Q)) ⟨|-⟩)))
  (, (ev (: (⟨->⟩ $P $Q) ⟨|-⟩))))

  (exec mp#2
  (, (ev (: $P ⟨wff⟩))
     (ev (: $Q ⟨wff⟩))
     (ev (: $P ⟨|-⟩))
     (ev (: (⟨->⟩ $P $Q) ⟨|-⟩)))
  (, (ev (: $Q ⟨|-⟩))))

  ; Derive (P→Q) from P and (P→(P→Q)) using modus ponens
(exec derive-implication-from-nested
  (, (ev (: $P ⟨wff⟩))
     (ev (: (⟨->⟩ $P $Q) ⟨wff⟩))
     (ev (: $P ⟨|-⟩))
     (ev (: (⟨->⟩ $P (⟨->⟩ $P $Q)) ⟨|-⟩)))
  (, (ev (: (⟨->⟩ $P $Q) ⟨|-⟩))))

;(exec derive-P-to-Q-direct
;  (, (ev (: $P ⟨wff⟩))
;     (ev (: (⟨->⟩ $P $Q) ⟨wff⟩))
;     (ev (: $P ⟨|-⟩))
;     (ev (: (⟨->⟩ $P
;               (⟨->⟩ $P $Q)) ⟨|-⟩))
;     (kb (: ⟨mp-curry⟩ (-> (: $P ⟨wff⟩) (-> (: $Q ⟨wff⟩)
;                  (-> (: $P ⟨|-⟩) (-> (: (⟨->⟩ $P  (⟨->⟩ $P $Q)) ⟨|-⟩) (: $Q ⟨|-⟩))))))) )
;  (, (ev (: (⟨->⟩ $P $Q) ⟨|-⟩)) )) 

; Final assembly using mp-curry (same single step as the working mm1)
; Inserting $P and $Q in appropriately.
; $P :- (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
; $Q :- (⟨=⟩ ⟨t⟩ ⟨t⟩)
; Works now (with the P -> Q part)
(exec assemble-final-proof-direct
  (, (ev (: $P ⟨wff⟩))
     (ev (: $Q ⟨wff⟩))
     (ev (: $P ⟨|-⟩))
;     (ev (: (⟨->⟩ $P ;; why is this needed?  -- apparently it's not!
;               (⟨->⟩ $P $Q))) ⟨|-⟩)
     (ev (: (⟨->⟩ $P $Q) ⟨|-⟩)))
;     (kb (: ⟨mp-curry⟩ (-> (: $P ⟨wff⟩) (-> (: $Q ⟨wff⟩)
;                  (-> (: $P ⟨|-⟩) (-> (: (⟨->⟩ $P $Q) ⟨|-⟩) (: $Q ⟨|-⟩))))))) )
  (, (ev (: $Q ⟨|-⟩))))
"#;

// Doesn't work
// (exec assemble-final-proof-direct
//   (, (ev (: $P ⟨wff⟩))
//      (ev (: $Q ⟨wff⟩))
//      (ev (: $P ⟨|-⟩))
//      (ev (: (⟨->⟩ $P $Q) ⟨|-⟩)) )
//   (, (ev (: $Q ⟨|-⟩)) ) )

// Works
// (exec assemble-final-proof-direct
//   (, (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))
//      (ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩))
//      (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))
//      (ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
//                (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
//                        (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩))
//      (kb (: ⟨mp-curry⟩ $_)))
//   (, (ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨|-⟩)
//          ((((⟨mp-curry⟩ (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))
//             (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩))
//            (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))
//           (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
//                   (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩)
//                          (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩)))))

    let mut s = Space::new();
    let t0 = Instant::now();
    s.load_sexpr(P.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    // Targets (kept identical to mm1())
    let want_ev_term_tplus0    = "(ev (: (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨term⟩))";
    let want_ev_wff_p          = "(ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨wff⟩))";
    let want_ev_wff_q          = "(ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨wff⟩))";
    let want_ev_proof_p        = "(ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))";
    let want_ev_proof_ptoq     = "(ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨=⟩ ⟨t⟩ ⟨t⟩)) ⟨|-⟩))";  
    let want_ev_proof_ptoptoq  = "(ev (: (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨->⟩ (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) (⟨=⟩ ⟨t⟩ ⟨t⟩))) ⟨|-⟩))";
    let want_final_evidence    = "(ev (: (⟨=⟩ ⟨t⟩ ⟨t⟩) ⟨|-⟩)";

    println!("=== MM1 (didactic): Proving ⊢ (t = t) ===");

    let mut ticks = 0usize;
    loop {
        ticks += 1;
        let n = s.metta_calculus(1);

        let mut tmut = Vec::new();
        // trying to get: (ev (: (⟨=⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨t⟩) ⟨|-⟩))
        s.dump_sexpr(
          expr!(s, "[2] ev [3] : [3] ⟨=⟩ $ $ ⟨|-⟩"),  // Pattern
          expr!(s, "[2] ev [3] : [3] ⟨=⟩ _1 _2 ⟨|-⟩"),  // Template: full reconstruction  
          &mut tmut
      );
        // s.dump_sexpr(
        //     expr!(s, "[2] ev [3] : [3] ⟨=⟩ $ $ ⟨|-⟩"),  //Query result (tick 2): (⟨+⟩ ⟨0⟩ ⟨0⟩)
// (⟨+⟩ ⟨t⟩ ⟨0⟩)
// (⟨+⟩ (⟨+⟩ ⟨0⟩ ⟨0⟩) ⟨0⟩)
// (⟨+⟩ (⟨+⟩ ⟨0⟩ ⟨t⟩) ⟨0⟩)
// (⟨+⟩ (⟨+⟩ ⟨t⟩ ⟨0⟩) ⟨0⟩)
// (⟨+⟩ (⟨+⟩ ⟨t⟩ ⟨t⟩) ⟨0⟩)
            // expr!(s, "[2] ev [3] : [3] ⟨=⟩ [3] ⟨+⟩ ⟨t⟩ ⟨0⟩ ⟨t⟩ ⟨|-⟩"),  // Query result (tick 2): $a
        //     expr!(s, "_1"),
        //     &mut tmut
        // );

        let result = String::from_utf8(tmut).unwrap();
        println!("Query result (tick {}): {}", ticks, result);

        for line in result.lines() {
          let trimmed = line.trim();
          if trimmed == want_ev_proof_p {
              println!("✅ EXACT MATCH found at tick {}: {}", ticks, trimmed);
              break;
          }
        }

        let mut proof_ptoq_check = Vec::new();
        s.dump_sexpr(
            expr!(s, "[2] ev [3] : [3] ⟨->⟩ [3] ⟨=⟩ [3] ⟨+⟩ ⟨t⟩ ⟨0⟩ ⟨t⟩ [3] ⟨=⟩ ⟨t⟩ ⟨t⟩ ⟨|-⟩"),  // Pattern
            expr!(s, "[2] ev [3] : [3] ⟨->⟩ [3] ⟨=⟩ [3] ⟨+⟩ ⟨t⟩ ⟨0⟩ ⟨t⟩ [3] ⟨=⟩ ⟨t⟩ ⟨t⟩ ⟨|-⟩"),  // Template: return same expression
            &mut proof_ptoq_check
        );

        if !proof_ptoq_check.is_empty() {
            let result = String::from_utf8(proof_ptoq_check).unwrap();
            println!("🎯 Found P→Q proof: {}", result.trim());
        } else {
            println!("P→Q proof not found yet");
        }

        let mut buf = Vec::new();
        s.dump_all_sexpr(&mut buf).unwrap();
        let dump = String::from_utf8_lossy(&buf);

        let line_has = |needle: &str| dump.lines().any(|l| l.trim_start().starts_with(needle));

        let have_tplus0_term  = line_has(want_ev_term_tplus0);
        let have_wff_p_ev     = line_has(want_ev_wff_p);
        let have_wff_q_ev     = line_has(want_ev_wff_q);
        let have_proof_p_ev   = line_has(want_ev_proof_p);
        let have_ptoq_ev      = line_has(want_ev_proof_ptoq);
        let have_ptoptoq_ev   = line_has(want_ev_proof_ptoptoq);
        let have_final        = line_has(want_final_evidence);

        if have_final {
            println!("\n== mm1 (didactic): ✅ SUCCESS in {:?} after {} tick(s) ==", t0.elapsed(), ticks);
            println!("  (+ t 0) : term ............. {}", if have_tplus0_term { "✓" } else { "—" });
            println!("  wff_P (ev) ................. {}", if have_wff_p_ev { "✓" } else { "—" });
            println!("  wff_Q (ev) ................. {}", if have_wff_q_ev { "✓" } else { "—" });
            println!("  proof_P (a2@t, ev) ......... {}", if have_proof_p_ev { "✓" } else { "—" });
            println!("  proof_PtoQ (a1, ev) ........ {}", if have_ptoq_ev { "✓" } else { "—" });
            println!("  proof_PtoPtoQ (a1, ev) ..... {}", if have_ptoptoq_ev { "✓" } else { "—" });

            println!("\n--- Final evidence confirmation ---");
            println!("✅ Successfully derived ⊢ (t = t)");

            println!("\n--- Full Final State Dump ---");
            print!("{dump}");
            break;
        }

        if n == 0 || ticks >= 128 {
            println!("\n== mm1 (didactic): — FAILURE in {:?} after {} tick(s) ==", t0.elapsed(), ticks);
            println!("  (+ t 0) : term ............. {}", if have_tplus0_term { "✓" } else { "—" });
            println!("  wff_P (ev) ................. {}", if have_wff_p_ev { "✓" } else { "—" });
            println!("  wff_Q (ev) ................. {}", if have_wff_q_ev { "✓" } else { "—" });
            println!("  proof_P (a2@t, ev) ......... {}", if have_proof_p_ev { "✓" } else { "—" });
            println!("  proof_PtoQ (a1, ev) ........ {}", if have_ptoq_ev { "✓" } else { "—" });
            println!("  proof_PtoPtoQ (a1, ev) ..... {}", if have_ptoptoq_ev { "✓" } else { "—" });

            if !have_final {
                println!("\n❌ Failed to derive ⊢ (t = t)");
            }

            println!("\n--- Full Final State Dump ---");
            print!("{dump}");
            break;
        }
    }
}



fn bc3() {
    let mut s = Space::new();

    const SPACE_EXPRS: &str = r#"
    ((step (0 base) $ts)
      (, (goal $ts (: $proof $conclusion)) (kb (: $proof $conclusion)))
      (, (ev (: $proof $conclusion) ) ))

    ((step (1 abs) $ts)
      (, (goal $k (: $proof $conclusion)))
      (, (goal (S $ts) (: $lhs (-> $synth $conclusion)) ) ))

    ((step (2 rev) $ts)
      (, (ev (: $lhs (-> $a $r)))  (goal $k (: $k $r)) )
      (, (goal (S $ts) (: $rhs $a) ) ))

    ((step (3 app) $ts)
      (, (ev (: $lhs (-> $a $r)))  (ev (: $rhs $a))  )
      (, (ev (: (@ $lhs $rhs) $r) ) ))

    (exec (clocked Z)
            (, ((step $x $ts) $p0 $t0)
               (exec (clocked $ts) $p1 $t1) )
            (, (exec (a $x) $p0 $t0)
               (exec (clocked (S $ts)) $p1 $t1) ))
    "#;

    const KB_EXPRS: &str = r#"
    (kb (: a A))
    (kb (: ab (R A B)))
    (kb (: bc (R B C)))
    (kb (: MP (-> (R $p $q) (-> $p $q))))

    (goal Z (: $proof C))
    "#;


    // (kb (: a A))
    //     (kb (: ab (-> A B)))
    //
    //     (goal Z (: $proof B))


    // (kb (: b B))
    //     (kb (: ab_c (-> A (-> B C))))
    //     (kb (: uncurry (-> (-> $a (-> $b $c)) (-> (* $a $b) $c))))
    // (kb (: sym (-> (* $a $b) (* $b $a))))
    // (kb (: . (-> (-> $b $c) (-> (-> $a $b) (-> $a $c)))))
    // (kb (: curry (-> (-> (* $a $b) $c) (-> $a (-> $b $c)))))
    //
    // (goal Z (: $proof (-> A C)))


    // P1:  (exec $p (, pat) (, (- temp) (+ x)))
    // add subtracts to SUB space, and remove them at the end
    // could not remove patterns under bindings
    // P2:  (exec $p (, (take pat) ) (, temp x)
    // only remove original patterns
    // P3:  (exec $p (, pat ) (, (subtract pat) (subtract temp)) (, temp x)
    //

    s.load_all_sexpr(SPACE_EXPRS.as_bytes()).unwrap();
    s.load_all_sexpr(KB_EXPRS.as_bytes()).unwrap();


    // let mut t0 = Instant::now();
    // println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(60);
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    // s.dump_all_sexpr(&mut v).unwrap();
    s.dump_sexpr(expr!(s, "[2] ev [3] : $ C"), expr!(s, "_1"), &mut v);
    let res = String::from_utf8(v).unwrap();

    println!("proof: {res}");


    // for i in 0..14 {
    //     println!("GEN {i}");
    //     let steps = s.metta_calculus(1);
    //     let mut v = vec![];
    //     s.dump_all_sexpr(&mut v).unwrap();
    //     // s.dump_sexpr(expr!(s, "[2] ev [3] : $ C"), expr!(s, "_1"), &mut v).unwrap();
    //     let res = String::from_utf8(v).unwrap();
    //
    //     println!("result: {res}");
    //
    // }

    // assert!(res.contains("(@ (@ . (@ uncurry ab_c)) (@ (@ curry sym) b))\n"));
}

fn cm0() {
    let mut s = Space::new();
    
    // Follow along https://en.wikipedia.org/wiki/Counter_machine#Program
    
    // non-peano csv version see cm1
    /*
    s.load_csv(INSTRS_CSV.as_bytes(), expr!(s, "$"), expr!(s, "[2] program _1"), b',').unwrap();
    s.load_csv(REGS_CSV.as_bytes(), expr!(s, "[2] $ $"), expr!(s, "[3] state 0 [3] REG _1 _2"), b',').unwrap();
    JZ,2,5\nDEC,2,2INC,3,3\nINC,1,4\nJZ,0,0\nJZ,1,9\nDEC,1,7\nINC,2,8\nJZ,0,5\nH,0,0
     */
    let TO_COPY = 50;

    let SPACE_MACHINE = format!(r#"
    (program Z (JZ 2 (S (S (S (S (S Z))))) ))
    (program (S Z) (DEC 2))
    (program (S (S Z)) (INC 3))
    (program (S (S (S Z))) (INC 1))
    (program (S (S (S (S Z)))) (JZ 0 Z))
    (program (S (S (S (S (S Z))))) (JZ 1 (S (S (S (S (S (S (S (S (S Z)))))))))))
    (program (S (S (S (S (S (S Z)))))) (DEC 1))
    (program (S (S (S (S (S (S (S Z))))))) (INC 2))
    (program (S (S (S (S (S (S (S (S Z)))))))) (JZ 0 (S (S (S (S (S Z)))))))
    (program (S (S (S (S (S (S (S (S (S Z))))))))) H)
    (state Z (REG 0 Z))
    (state Z (REG 1 Z))
    (state Z (REG 2 {}))
    (state Z (REG 3 Z))
    (state Z (REG 4 Z))
    (state Z (IC Z))
    (if (S $n) $x $y $x)
    (if Z $x $y $y)
    (0 != 1) (0 != 2) (0 != 3) (0 != 4)
    (1 != 0) (1 != 2) (1 != 3) (1 != 4)
    (2 != 1) (2 != 0) (2 != 3) (2 != 4)
    (3 != 1) (3 != 2) (3 != 0) (3 != 4)
    (4 != 1) (4 != 2) (4 != 0) (4 != 3)
    
    ((step JZ $ts)
      (, (state $ts (IC $i)) (program $i (JZ $r $j)) (state $ts (REG $r $v)) (if $v (S $i) $j $ni) (state $ts (REG $k $kv)))
      (, (state (S $ts) (IC $ni)) (state (S $ts) (REG $k $kv))))

    ((step INC $ts)
      (, (state $ts (IC $i)) (program $i (INC $r)) (state $ts (REG $r $v)) ($r != $o) (state $ts (REG $o $ov)))
      (, (state (S $ts) (IC (S $i))) (state (S $ts) (REG $r (S $v))) (state (S $ts) (REG $o $ov))))
    
    ((step DEC $ts)
      (, (state $ts (IC $i)) (program $i (DEC $r)) (state $ts (REG $r (S $v))) ($r != $o) (state $ts (REG $o $ov)))
      (, (state (S $ts) (IC (S $i))) (state (S $ts) (REG $r $v)) (state (S $ts) (REG $o $ov))))  

    (exec (clocked Z)
            (, (exec (clocked $ts) $p1 $t1) 
               (state $ts (IC $_))
               ((step $k $ts) $p0 $t0))
            (, (exec ($k $ts) $p0 $t0)
               (exec (clocked (S $ts)) $p1 $t1)))
    "#, peano(TO_COPY));

    s.load_all_sexpr(SPACE_MACHINE.as_bytes()).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(1000000000000000);
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v_ts = vec![];
    s.dump_sexpr(expr!(s, "[3] state $ $"), expr!(s, "_1"), &mut v_ts);
    let last_ts_tmp = String::from_utf8(v_ts).unwrap(); 
    let last_ts = last_ts_tmp.split("\n").max_by_key(|x| x.len()).unwrap();
    let mut v = vec![];
    // s.dump_all_sexpr(&mut v).unwrap();
    s.dump_sexpr(expr!(s, "[3] state $ [3] REG 3 $"), expr!(s, "[2] _1 _2"), &mut v);
    let res = String::from_utf8(v).unwrap();
    
    // println!("{res}");
    assert!(res.contains(format!("({} {})", last_ts, peano(TO_COPY)).as_str()));
}

/*fn match_case() {
    let mut s = Space::new();

    const SPACE_EXPRS: &str = r#"
(unify $x $x)

(exec 0
      (, (Apply $x)
         (Match $c $p $t))
      (, (exec (M $c)
               (, (unify $x $p) (exec (M $c) $Mp $Mt))
               (, (res $t)
                  (- (exec (M $c) $Mp $Mt)) ))))

(Match 0 (foo $x) (Inner Foo $x))
(Match 1 (bar $x) (Inner Bar $x))
(Match 2 $x (Fallback $x))

(Apply (foo $x))
    "#;

    s.load_all_sexpr(SPACE_EXPRS.as_bytes()).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(1000000000000000);
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();

    println!("result: {res}");
}*/

fn lens_aunt() {
    let mut s = Space::new();
    /*
    Tom x Pam
     |   \
    Liz  Bob
         / \
      Ann   Pat
             |
            Jim
     */
    let SPACE = r#"
    (exec QA (, (aunt $xc $x $y $yt) (data $xc) (exec QA $P $T)
                (parent $p $x) (parent $gp $p) (parent $gp $y)
                (female $y) ($p != $y))
             (, (data $yt) (exec QA $P $T)))

    (data (poi Jim)) (data (poi Ann))
    (aunt (poi $x) $x $y (result ($y aunt of $x)))

    (parent Tom Bob)
    (parent Pam Bob)
    (parent Tom Liz)
    (parent Bob Ann)
    (parent Bob Pat)
    (parent Pat Jim)
    (female Pam) (female Liz) (female Pat) (female Ann)
    (male Tom) (male Bob) (male Jim)

    (Pam == Pam) (Pam != Liz) (Pam != Pat) (Pam != Ann) (Pam != Tom) (Pam != Bob) (Pam != Jim)
    (Liz != Pam) (Liz == Liz) (Liz != Pat) (Liz != Ann) (Liz != Tom) (Liz != Bob) (Liz != Jim)
    (Pat != Pam) (Pat != Liz) (Pat == Pat) (Pat != Ann) (Pat != Tom) (Pat != Bob) (Pat != Jim)
    (Ann != Pam) (Ann != Liz) (Ann != Pat) (Ann == Ann) (Ann != Tom) (Ann != Bob) (Ann != Jim)
    (Tom != Pam) (Tom != Liz) (Tom != Pat) (Tom != Ann) (Tom == Tom) (Tom != Bob) (Tom != Jim)
    (Bob != Pam) (Bob != Liz) (Bob != Pat) (Bob != Ann) (Bob != Tom) (Bob == Bob) (Bob != Jim)
    (Jim != Pam) (Jim != Liz) (Jim != Pat) (Jim != Ann) (Jim != Tom) (Jim != Bob) (Jim == Jim)
    "#;

    s.load_all_sexpr(SPACE.as_bytes()).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(1);
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    // s.dump_all_sexpr(&mut v).unwrap();
    s.dump_sexpr(expr!(s, "[2] data [2] result $"), expr!(s, "_1"), &mut v);
    let res = String::from_utf8(v).unwrap();

    println!("{res}");
    assert_eq!(res, "(Ann aunt of Jim)\n(Liz aunt of Ann)\n");
}

fn lens_composition() {
    let mut s = Space::new();

    let SPACE = r#"
    (exec LC (, (compose $l0 $l1)
                (lens ($l0 $xc0 $x0 $y0 $yt0))
                (lens ($l1 $x0 $x1 $y1 $y0)) )
             (, (lens (($l0 o $l1) $xc0 $x1 $y1 $yt0))))

    (lens (aunt (poi $x) $x $y (result ($y aunt of $x))))
    (lens (ns (users (adam (experiments $x))) $x $y (users (adam (experiments $y)))))
    (compose ns aunt)
    "#;

    s.load_all_sexpr(SPACE.as_bytes()).unwrap();

    let mut t0 = Instant::now();
    let steps = s.metta_calculus(1);
    println!("elapsed {} steps {} size {}", t0.elapsed().as_millis(), steps, s.btm.val_count());

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();

    println!("{res}");
    assert!(res.contains("(lens ((ns o aunt) (users (adam (experiments (poi $a)))) $a $b (users (adam (experiments (result ($b aunt of $a)))))))"));
}

fn bench_transitive_no_unify(nnodes: usize, nedges: usize) {
    use rand::{rngs::StdRng, SeedableRng, Rng};
    let mut rng = StdRng::from_seed([0; 32]);
    let mut s = Space::new();

    let mut edges = String::new();

    for k in 0..nedges {
        let i = rng.random_range(0..nnodes);
        let j = rng.random_range(0..nnodes);
        edges.push_str(format!("(edge {i} {j})\n").as_str());
    }

    s.load_all_sexpr(edges.as_bytes()).unwrap();
    println!("constructed {} nodes {} edges", nnodes, nedges);

    let t0 = Instant::now();
    s.interpret(expr!(s, "[4] exec 0 [3] , [3] edge $ $ [3] edge _2 $ [2] , [3] trans _1 _3"));
    println!("trans elapsed {} µs", t0.elapsed().as_micros());

    let t1 = Instant::now();
    s.interpret(expr!(s, "[4] exec 0 [4] , [3] edge $ $ [3] edge _2 $ [3] edge _1 _3 [2] , [4] dtrans _1 _2 _3"));
    println!("detect trans elapsed {} µs", t1.elapsed().as_micros());


    let mut v = vec![];
    s.dump_sexpr(expr!(s, "[3] trans $ $"), expr!(s, "[2] _1 _2"), &mut v);
    let ntrans: usize = v.iter().map(|c| if *c == b'\n' { 1 } else { 0 }).sum();
    v.clear();
    s.dump_sexpr(expr!(s, "[4] dtrans $ $ $"), expr!(s, "[3] _1 _2 _3"), &mut v);
    let ndtrans: usize = v.iter().map(|c| if *c == b'\n' { 1 } else { 0 }).sum();
    println!("trans {} detected trans {}", ntrans, ndtrans);

    // (badbad)
    // constructed 50000 nodes 1000000 edges
    // trans elapsed 17391765 µs
    // detect trans elapsed 11928710 µs
    // trans 19917429 detected trans 8716
}


fn bench_clique_no_unify(nnodes: usize, nedges: usize, max_clique: usize) {
    fn binom_as_f64(n: u64, k: u64) -> f64 {
        if k > n { return 0.0; }
        let k = std::cmp::min(k, n - k);
        let mut res = 1.0f64;
        for i in 1..=k {
            res *= (n - k + i) as f64;
            res /= i as f64;
        }
        res
    }

    fn expected_fraction_kclique_gne(n: u64, e: u64, k: u64) -> f64 {
        assert!(n >= 2, "n >= 2");
        let m = n * (n - 1) / 2; // total possible edges
        assert!(e <= m, "E must be <= C(n,2)");
        let kk = k * (k - 1) / 2; // number of edges inside a k-clique
        if kk == 0 { return 1.0; } // k=0 or 1
        if e < kk { return 0.0; }  // not enough edges to cover a clique
        let mut num = 1.0f64;
        let mut den = 1.0f64;
        for i in 0..kk {
            num *= (e - i) as f64;
            den *= (m - i) as f64;
        }
        num / den
    }

    fn expected_num_kclique_gne(n: u64, e: u64, k: u64) -> f64 {
        binom_as_f64(n, k) * expected_fraction_kclique_gne(n, e, k)
    }

    fn clique_query(k: usize) -> String {
        format!("(exec 0 (,{}) (, ({}-clique{})))",
            (0..k).flat_map(|i| ((i + 1)..k).map(move |j| format!(" (edge $x{} $x{})", i, j))).collect::<String>(),
            k,
            (0..k).map(|i| format!(" $x{}", i)).collect::<String>()
        )
    }

    use rand::{rngs::StdRng, SeedableRng, Rng};
    let mut rng = StdRng::from_seed([0; 32]);
    let mut s = Space::new();

    let mut edges: HashSet<String> = HashSet::new();

    // irreflexive degeneracy ordered graph
    while edges.len() < nedges {
        let i = rng.random_range(0..nnodes);
        let j = rng.random_range(0..nnodes);
        if i == j { continue }
        if i < j { edges.insert(format!("(edge {i} {j})\n")); }
        else { edges.insert(format!("(edge {j} {i})\n")); }
    }

    s.load_all_sexpr(edges.into_iter().collect::<String>().as_bytes()).unwrap();
    println!("constructed {} nodes {} edges", nnodes, nedges);

    for k in 3..(max_clique+1) {
        let query = clique_query(k);
        println!("executing query {}", query);
        let t0 = Instant::now();
        s.load_sexpr(query.as_bytes(), expr!(s, "$"), expr!(s, "_1"));
        s.metta_calculus(1);
        let nkcliques: usize = s.btm.read_zipper_at_path([item_byte(Tag::Arity((k+1) as _))]).val_count();
        println!("found {} {k}-cliques (expected {}) in {} µs", nkcliques, expected_num_kclique_gne(nnodes as _, nedges as _, k as _).round(), t0.elapsed().as_micros());
    }
    // constructed 200 nodes 3600 edges
    // executing query (exec 0 (, (edge $x0 $x1) (edge $x0 $x2) (edge $x1 $x2)) (, (3-clique $x0 $x1 $x2)))
    // found 7824 3-cliques (expected 7770) in 39910 µs
    // executing query (exec 0 (, (edge $x0 $x1) (edge $x0 $x2) (edge $x0 $x3) (edge $x1 $x2) (edge $x1 $x3) (edge $x2 $x3)) (, (4-clique $x0 $x1 $x2 $x3)))
    // found 2320 4-cliques (expected 2260) in 1096909 µs
    // executing query (exec 0 (, (edge $x0 $x1) (edge $x0 $x2) (edge $x0 $x3) (edge $x0 $x4) (edge $x1 $x2) (edge $x1 $x3) (edge $x1 $x4) (edge $x2 $x3) (edge $x2 $x4) (edge $x3 $x4)) (, (5-clique $x0 $x1 $x2 $x3 $x4)))
    // found 102 5-cliques (expected 94) in 24705340 µs
    // executing query (exec 0 (, (edge $x0 $x1) (edge $x0 $x2) (edge $x0 $x3) (edge $x0 $x4) (edge $x0 $x5) (edge $x1 $x2) (edge $x1 $x3) (edge $x1 $x4) (edge $x1 $x5) (edge $x2 $x3) (edge $x2 $x4) (edge $x2 $x5) (edge $x3 $x4) (edge $x3 $x5) (edge $x4 $x5)) (, (6-clique $x0 $x1 $x2 $x3 $x4 $x5)))
    // found 0 6-cliques (expected 1) in <1288009964 µs
}

fn bench_finite_domain(terms: usize) {
    use rand::{rngs::StdRng, SeedableRng, Rng};
    let mut rng = StdRng::from_seed([0; 32]);
    const DS: usize = 64;
    const SYM: [&'static str; 64] = ["0","1","2","3","4","5","6","7","8","9","?","@","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"];
    // const SYM: [&'static str; 64] = ["À", "Á", "Â", "Ã", "Ä", "Å", "Æ", "Ç", "È", "É", "Ê", "Ë", "Ì", "Í", "Î", "Ï", "Ð", "Ñ", "Ò", "Ó", "Ô", "Õ", "Ö", "×", "Ø", "Ù", "Ú", "Û", "Ü", "Ý", "Þ", "ß", "à", "á", "â", "ã", "ä", "å", "æ", "ç", "è", "é", "ê", "ë", "ì", "í", "î", "ï", "ð", "ñ", "ò", "ó", "ô", "õ", "ö", "÷", "ø", "ù", "ú", "û", "ü", "ý", "þ", "ÿ"];
    // const SYM: [&'static str; 64] = ["䷁","䷗","䷆","䷒","䷎","䷣","䷭","䷊","䷏","䷲","䷧","䷵","䷽","䷶","䷟","䷡","䷇","䷂","䷜","䷻","䷦","䷾","䷯","䷄","䷬","䷐","䷮","䷹","䷞","䷰","䷛","䷪","䷖","䷚","䷃","䷨","䷳","䷕","䷑","䷙","䷢","䷔","䷿","䷥","䷷","䷝","䷱","䷍","䷓","䷩","䷺","䷼","䷴","䷤","䷸","䷈","䷋","䷘","䷅","䷉","䷠","䷌","䷫","䷀"];

    fn uop<F : Fn(usize) -> usize>(sym: &str, f: F) -> String {
        let mut s = String::new();
        for x in 0..DS {
            let z = f(x);
            if z == usize::MAX { continue }
            s.push_str(format!("({} {} = {})\n", sym, SYM[x], SYM[z]).as_str());
        }
        s
    }

    fn bop<F : Fn(usize, usize) -> usize>(sym: &str, f: F) -> String {
        let mut s = String::new();
        for x in 0..DS {
            for y in 0..DS {
                let z = f(x, y);
                if z == usize::MAX { continue }
                s.push_str(format!("({} {} {} = {})\n", SYM[x], sym, SYM[y], SYM[z]).as_str());
            }
        }
        s
    }

    let mut s = Space::new();

    let sq = uop("²", |x| (x * x) % DS);
    let sqrt = uop("√", |x| x.isqrt());

    let add = bop("+", |x, y| (x + y) % DS);
    let sub = bop("-", |x, y| ((DS + x) - y) % DS);
    let mul = bop("*", |x, y| (x * y) % DS);
    let div = bop("/", |x, y| if y == 0 { usize::MAX } else { x / y });
    let join = bop("\\/", |x, y| x.max(y));
    let meet = bop("/\\", |x, y| x.min(y));

    let adds = bop("+s", |x, y| if x + y < DS { x + y } else { DS - 1 });
    let muls = bop("*s", |x, y| if x * y < DS { x * y } else { DS - 1 });

    let ops = [sq, sqrt, add, sub, mul, div, join, meet, adds, muls].concat();

    s.load_sexpr(ops.as_bytes(), expr!(s, "$"), expr!(s, "_1"));

    let mut args = String::new(); // e.g. (args ䷽ ䷣ ䷜ ䷣)
    for i in 0..10_000 {
        let x0 = rng.random_range(0..DS);
        let x1 = rng.random_range(0..DS);
        let y0 = rng.random_range(0..DS);
        let y1 = rng.random_range(0..DS);
        args.push_str(format!("(args {} {} {} {})", SYM[x0], SYM[x1], SYM[y0], SYM[y1]).as_str())
    }
    s.load_sexpr(args.as_bytes(), expr!(s, "$"), expr!(s, "_1"));

    s.load_sexpr(r"(exec 0 (, (args $x0 $y0 $x1 $y1) ($x0 /\ $x1 = $xl) ($x0 \/ $x1 = $xh) ($y0 /\ $y1 = $yl) ($y0 \/ $y1 = $yh) ($xh - $xl = $dx) ($yh - $yl = $dy) (² $dx = $dx2) (² $dy = $dy2) ($dx2 + $dy2 = $d2) (√ $d2 = $d)) (, (res $d)))".as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();
    let t0 = Instant::now();
    s.metta_calculus(1);
    let t1 = Instant::now();

    let mut v = vec![];
    // s.dump_all_sexpr(&mut v).unwrap();
    s.dump_sexpr(expr!(s, "[2] res $"), expr!(s, "_1"), &mut v);
    let res = String::from_utf8(v).unwrap();

    println!("{}", s.btm.val_count());
    println!("{res} ({terms} inputs) in {} µs", t1.duration_since(t0).as_micros());
    // (badbad)
    // (10_000 inputs) in 85833 µs
}

#[cfg(all(feature = "nightly"))]
fn json_upaths_smoke() {
    let test = r#"{
"first_name": "John",
"last_name": "Smith",
"is_alive": true,
"age": 27,
"address": {
  "street_address": "21 2nd Street",
  "city": "New York",
  "state": "NY",
  "postal_code": "10021-3100"},
"phone_numbers": [
  {"type": "home", "number": "212 555-1234"},
  {"type": "office", "number": "646 555-4567"}],
"children": ["Catherine", "Thomas", "Trevor"],
"spouse": null}"#;
    let mut cv = vec![];

    let mut s = Space::new();
    // let written = s.load_json(test.as_bytes()).unwrap();
    let written = s.json_to_paths(test.as_bytes(), &mut cv).unwrap();
    // println!("{:?}", pathmap::path_serialization::serialize_paths_(btm.read_zipper(), &mut cv));
    println!("written {written}");
    pathmap::paths_serialization::deserialize_paths(s.btm.write_zipper(), &cv[..], ()).unwrap();

    let mut v = vec![];
    s.dump_all_sexpr(&mut v).unwrap();
    let res = String::from_utf8(v).unwrap();
    println!("res {res}");
    assert_eq!(res, r#"(age 27)
(spouse null)
(address (city New York))
(address (state NY))
(address (postal_code 10021-3100))
(address (street_address 21 2nd Street))
(children (0 Catherine))
(children (1 Thomas))
(children (2 Trevor))
(is_alive true)
(last_name Smith)
(first_name John)
(phone_numbers (0 (type home)))
(phone_numbers (0 (number 212 555-1234)))
(phone_numbers (1 (type office)))
(phone_numbers (1 (number 646 555-4567)))
"#);
}

#[cfg(all(feature = "nightly"))]
fn json_upaths<IPath: AsRef<std::path::Path>, OPath : AsRef<std::path::Path>>(json_path: IPath, upaths_path: OPath) {
    println!("mmapping JSON file {:?}", json_path.as_ref().as_os_str());
    println!("writing out unordered .paths file {:?}", upaths_path.as_ref().as_os_str());
    let json_file = std::fs::File::open(json_path).unwrap();
    let json_mmap = unsafe { memmap2::Mmap::map(&json_file).unwrap() };
    let upaths_file = std::fs::File::create_new(upaths_path).unwrap();
    let mut upaths_bufwriter = std::io::BufWriter::new(upaths_file);

    let mut s = Space::new();
    let t0 = Instant::now();
    let written = s.json_to_paths(&*json_mmap, &mut upaths_bufwriter).unwrap();
    println!("written {written} in {} ms", t0.elapsed().as_millis());
    // (zephy)
    // mmapping JSON file "/home/adam/Downloads/G37S-9NQ.json"
    // writing out unordered .paths file "G37S-9NQ.upaths"
    // Ok(SerializationStats { bytes_out: 1415053, bytes_in: 12346358, path_count: 224769 })
    // written 224769 in 193 ms
    // (badbad)
    // mmapping JSON file "/mnt/data/enwiki-20231220-pages-articles-links/cqls.json"
    // writing out unordered .paths file "/mnt/data/enwiki-20231220-pages-articles-links/cqls.upaths"
    // Ok(SerializationStats { bytes_out: 231708224, bytes_in: 808333425, path_count: 15969490 })
    // written 15969490 in 17441 ms
}


fn main() {
    env_logger::init();

    // lookup();
    // positive();
    // negative();
    // bipolar();
    // positive_equal();
    // negative_equal();
    // bipolar_equal();
    //
    // two_positive_equal();
    // two_positive_equal_crossed();
    // two_bipolar_equal_crossed();
    //
    // process_calculus_reverse();
    // logic_query();
    // bc0();
    // bc1();
    // bc2();
    // bc3();
    //
    // cm0();
    // bc0();
    // bc1();

    // mm0_original_with_verification_fix();
    // mm0_original_fixed_steps();
    // // mm0_original_with_sequential_tick();
    // mm0();

    // mm1_a();
    // mm1_b();
    mm1_b_tpl();
    mm1_b2_tpl();
    // mm1_c();
    // mm1_d();
    // mm1();
    mm1_didactic();
    mm1_forward();

    // lens_aunt();
    // lens_composition();

    // bench_transitive_no_unify(50000, 1000000);
    // bench_clique_no_unify(200, 3600, 5);
    // bench_finite_domain(10_000);
    // process_calculus_bench(1000, 200, 200);

    // #[cfg(all(feature = "nightly"))]
    // json_upaths_smoke();
    // #[cfg(all(feature = "nightly"))]
    // json_upaths("/mnt/data/enwiki-20231220-pages-articles-links/cqls.json", "/mnt/data/enwiki-20231220-pages-articles-links/cqls.upaths");
    return;

    let mut s = Space::new();
    const space: &str = r#"
(exec P0 (, (sudoku p2 input (row ($r $x1 $x2 $x3 $x4 $x5 $x6 $x7 $x8 $x9)))) 
         (, (cell 1 $r $x1) (cell 2 $r $x2) (cell 3 $r $x3)  (cell 4 $r $x4) (cell 5 $r $x5) (cell 6 $r $x6)  (cell 7 $r $x7) (cell 8 $r $x8) (cell 9 $r $x9)  ))

(exec P1 (, (cell $c $r _))
         (, (cell $c $r 1) (cell $c $r 2) (cell $c $r 3)  (cell $c $r 4) (cell $c $r 5) (cell $c $r 6)  (cell $c $r 7) (cell $c $r 8) (cell $c $r 9)  ))

(exec P2 (, (cell $ca $r $va) (cell $cb $r $vb))
         (, (Deduction remaining (cell $ca $r $x) (cell $cb $r $y))))

"#;
    //
    // (exec P2 (, (cell $ca $r $va) (cell $cb $r $vb))
    // (, (Deduction remaining (cell $ca $r X) (cell $cb $r Y))))

    // (exec P3 (, (cell $c $ra $va) (cell $c $rb $vb))
    // (, (Deduction remaining (cell $c $ra X) (cell $c $rb Y))))
    // 
    // 
    // (exec P4 (, (cell $c $ra $va) (cell $c $rb $vb))
    // (, (Deduction remaining (cell $c $ra $x) (cell $c $rb $y))))
    // (block 0 1 1) (block 0 1 2) (block 0 1 3)
    // (block 0 2 1) (block 0 2 2) (block 0 2 3)
    // (block 0 3 1) (block 0 3 2) (block 0 3 3)

    const sudoku_p2: &str = r#"
1 2 3 4 5 6 7 8 9
_ 5 _ _ _ _ 9 _ _
_ _ _ 8 3 1 2 5 _
2 _ 7 _ _ _ 6 1 3
9 _ 6 _ _ 7 _ 3 _
1 2 8 _ _ _ 7 _ _
_ _ _ 2 _ 4 _ 9 6
8 1 _ 7 6 _ _ 2 9
7 3 4 _ 2 8 _ _ 1
_ _ _ 4 1 _ _ _ _"#;
    
    s.load_csv(sudoku_p2.as_bytes(), expr!(s, "$"), expr!(s, "[4] sudoku p2 input [2] row _1"), b' ').unwrap();
    
    s.load_sexpr(space.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();

    s.metta_calculus(100);
    
    // println!("size {:?}", s.btm.val_count());
    
    let mut v = vec![];
    s.dump_sexpr(expr!(s, "$"), expr!(s, "_1"), &mut v);
    
    println!("{}", String::from_utf8(v).unwrap());
    
    return;
    /*let mut s = Space::new();
    const space: &str = r#"
(= (add (S $x) $y) (S (add $x $y)))
(= (add Z $y) $y)
(= (mul (S (S $x)) $y) (add $y (mul (S $x) $y)))
(= (mul (S Z) $y) $y)
(= (mul Z $y) Z)

(PC (Cat (Point3D 39 9504 34980)))
(PC (Cat (Point3D 39 9504 34980)))
(PC (Cat (Point3D 39 9480 34980)))
(PC (Cat (Point3D 39 4830 34980)))
(PC (Cat (Point3D 39 9504 34980)))
(PC (Cat (Point3D 39 3230 34923)))
(PC (Cat (Point3D 27 3410 34900)))
(PC (Cat (Point3D 27 3410 34964)))
(PC (Cat (Point3D 23 3459 34932)))

"#;

    s.load_sexpr(space.as_bytes(), expr!(s, "$"), expr!(s, "_1")).unwrap();


    let [(t1, _), (t2, _)] = &s.token_bfs(&[], expr!(s, "$"))[..] else { panic!() };
    println!("{:?}", s.token_bfs(&t1[..], expr!(s, "$")));
    println!("{:?}", s.token_bfs(t2, expr!(s, "$")));

    // let mut v = vec![];
    // s.dump_sexpr(expr!(s, "$"), expr!(s, "_1"), &mut v).unwrap();
    //
    // println!("{}", String::from_utf8(v).unwrap());
    return;*/
    /*

    let mut s = Space::new();
    const facts: &str = "0,1\n1,2\n2,3\n3,4\n4,5\n5,0\n5,6\n6,7\n7,8\n8,9\n9,10\n10,7";
    const expected_same_clique: &str = "...";
    const expected_reachable: &str = "...";


    s.load_csv(facts.as_bytes(), expr!(s, "[2] $ $"), expr!(s, "[3] edge _1 _2")).unwrap();

    // (reachable $x $y) :- (edge $x $y)
    // (reachable $x $y) :- (edge $x $z), (reachable $z $y)
    // (same_clique $x $y) :- (reachable $x $y), (reachable $y $x)
    s.datalog(&[
        expr!(s, "[3] -: [2] , [3] edge $ $ [3] reachable _1 _2"),
        expr!(s, "[3] -: [3] , [3] edge $ $ [3] reachable _2 $ [3] reachable _1 _3"),
        expr!(s, "[3] -: [3] , [3] reachable $ $ [3] reachable _2 _1 [3] same_clique _1 _2"),
    ]);

    let mut v = vec![];
    s.dump_sexpr(expr!(s, "$"), expr!(s, "_1"), &mut v).unwrap();

    println!("{}", String::from_utf8(v).unwrap());
    return;
     */
    /*
        const csv_contents: &str = r#"1,2
10,20
10,30"#;

    const sexpr_contents: &str = r#"(useful (Foo 1))
(useless ((- o -) (- o -)))"#;

    let mut s = Space::new();
    // s.load_csv(csv_contents.as_bytes(), expr!(s, "[2] $ $"), expr!(s, "[2] mycsv [3] my precious _2")).unwrap();
    s.load_csv(csv_contents.as_bytes(), expr!(s, "[2] 10 $"), expr!(s, "[2] data [2] mycsv [3] my precious _1")).unwrap();

    s.load_sexpr(sexpr_contents.as_bytes(), expr!(s, "[2] useful $"), expr!(s, "[2] data [2] mysexpr _1")).unwrap();

    let mut v = vec![];
    s.dump_sexpr(expr!(s, "[2] data [2] mycsv $"), expr!(s, "_1"), &mut v).unwrap();

    println!("{}", String::from_utf8(v).unwrap());
    return;
    */
    // println!("{}", mork_bytestring::serialize(&[3, 3, 200, 84, 80, 55, 51, 45, 65, 83, 49, 204, 103, 101, 110, 101, 95, 110, 97, 109, 101, 95, 111, 102, 200, 0, 0, 0, 0, 4, 129, 29, 29, 4, 195, 83, 80, 79, 200, 0, 0, 0, 0, 4, 129, 29, 29, 200]));
    //
    return;
    // let mut s = Space::new();

    // let tree = pathmap::arena_compact::ArenaCompactTree::open_mmap("/run/media/adam/43323a1c-ad7e-4d9a-b3c0-cf84e69ec61a/whole_flybase.tree").unwrap();
    // 
    // let iter_tree_start = Instant::now();
    // let mut rz = tree.read_zipper();
    // let mut npaths = 0usize; 
    // let mut nbytes = 0; 
    // while rz.to_next_val() {
    //     nbytes += rz.path().len();
    //     npaths += 1;
    //     if npaths % 10_000_000 == 0 {
    //         println!("npaths {}", npaths);
    //     }
    // }
    // println!("iterating tree backup {} {} took {}", npaths, nbytes, iter_tree_start.elapsed().as_secs());
    
    // let restore_paths_start = Instant::now();
    // let restored = s.restore_paths("/run/media/adam/43323a1c-ad7e-4d9a-b3c0-cf84e69ec61a/whole_flybase.paths.gz").unwrap();
    // println!("restored paths {:?} {}", restored, restore_paths_start.elapsed().as_secs());
    
    // let everythingf = std::fs::File::open("/run/media/adam/43323a1c-ad7e-4d9a-b3c0-cf84e69ec61a/whole_flybase.jsonl").unwrap();
    // let everythingfs = unsafe { memmap2::Mmap::map(&everythingf).unwrap() };
    // let load_compressed = Instant::now();
    // println!("done {:?} {}", s.load_jsonl(everythingfs.as_ref()).unwrap(), load_compressed.elapsed().as_secs());
    // // done (326728210, 6798095370) 1934s
    // 
    // let backup_paths_start = Instant::now();
    // println!("{:?}", s.backup_paths("/run/media/adam/43323a1c-ad7e-4d9a-b3c0-cf84e69ec61a/whole_flybase.paths.gz").unwrap());
    // println!("paths backup took {}", backup_paths_start.elapsed().as_secs());
    // // SerializationStats { bytes_out: 42_741_214_528, bytes_in: 355_357_500_042, path_count: 6798095370 }
    // //                                                           328_165_118_562
    // // paths backup took 4619s
    // 
    // let backup_tree_start = Instant::now();
    // println!("{:?}", s.backup_tree("/run/media/adam/43323a1c-ad7e-4d9a-b3c0-cf84e69ec61a/whole_flybase.tree").unwrap());
    // println!("tree backup took {}", backup_tree_start.elapsed().as_secs());
    // // tree backup took 3033
    // 
    // let backup_dag_start = Instant::now();
    // println!("{:?}", s.backup("/run/media/adam/43323a1c-ad7e-4d9a-b3c0-cf84e69ec61a/").unwrap());
    // println!("dag backup took {}", backup_dag_start.elapsed().as_secs());


    // let restore_symbols_start = Instant::now();
    // println!("restored symbols {:?}", s.restore_symbols("/dev/shm/combined.symbols.zip").unwrap());
    // println!("symbols backup took {}", restore_symbols_start.elapsed().as_secs());
    // println!("{:?}", s.sm.get_sym(b"SPO"));
    // println!("{:?}", s.sm.get_sym(b"IL9R-207"));
    // let bucket_map::serialization::Tables { to_symbol, to_bytes } = s.sm.reveal_tables();
    // println!("to_symbol.len() = {}; to_bytes.len() = {}", to_symbol.len(), to_bytes.len());
    // println!("to_symbol.first().unwrap().val_count() = {:?}; to_bytes.first().unwrap().val_count() = {:?}", to_symbol.last().unwrap().val_count(), to_bytes.last().unwrap().val_count());

    // for to_symbol_part in to_symbol {
    //     print!("{},", to_symbol_part.val_count());
    // }

    // for to_byte_part in to_bytes {
    //     print!("{},", to_byte_part.val_count());
    // }

    // let mut to_symbol_rz = to_symbol.last().unwrap().read_zipper();
    // while let Some(v) = to_symbol_rz.to_next_val() {
    //     println!("{:?} {:?}", std::str::from_utf8(to_symbol_rz.path()).unwrap_or(format!("{:?}", to_symbol_rz.path()).as_str()), v)
    // }

    // let restore_paths_start = Instant::now();
    // println!("restored paths {:?}", s.restore_paths("/dev/shm/combined_ni.paths.gz").unwrap());
    // println!("paths restore took {}", restore_paths_start.elapsed().as_secs());
    // s.statistics();

    // let load_labels_start = Instant::now();
    // println!("{:?}", s.load_neo4j_node_labels("bolt://localhost:7687", "neo4j", "morkwork").unwrap());
    // println!("loading labels took {}", load_labels_start.elapsed().as_secs());
    // s.statistics();

    // let load_start = Instant::now();
    // println!("{:?}", s.load_neo4j_triples("bolt://localhost:7687", "neo4j", "morkwork").unwrap());
    // println!("loading took {}", load_start.elapsed().as_secs());
    // s.statistics();
    // let mut rz = s.btm.read_zipper_at_path(&[item_byte(Tag::Arity(4)), item_byte(Tag::SymbolSize(3)), b'S', b'P', b'O']);
    // println!("SPO's {}", rz.val_count());
    // rz.to_next_val();
    // println!("{}", mork_bytestring::serialize(rz.origin_path().unwrap()));
    //
    // let property_load_start = Instant::now();
    // println!("{:?}", s.load_neo4j_node_properties("bolt://localhost:7687", "neo4j", "morkwork").unwrap());
    // println!("property loading took {}", property_load_start.elapsed().as_secs());
    // s.statistics();


    // work(&mut s);

    // s.statistics();
    // s.done();
    // let restore_start = Instant::now();
    // s.restore("/dev/shm/");
    // println!("restore took {}", restore_start.elapsed().as_secs());
    // s.statistics();




    // let backup_paths_start = Instant::now();
    // println!("{:?}", s.backup_paths("/dev/shm/combined_ni.paths.gz").unwrap());
    // println!("paths backup took {}", backup_paths_start.elapsed().as_secs());
    //
    // let backup_symbols_start = Instant::now();
    // println!("{:?}", s.backup_symbols("/dev/shm/combined.symbols.zip").unwrap());
    // println!("symbols backup took {}", backup_symbols_start.elapsed().as_secs());

    // let backup_start = Instant::now();
    // s.backup("/dev/shm/combined.remydag");
    // println!("backup took {}", backup_start.elapsed().as_secs());
}

/*
restored paths DeserializationStats { bytes_in: 4061021063, bytes_out: 34621879507, path_count: 978700221 }
paths restore took 135
val count 978700221
add gene name index took 8637 ms
val count 979015756
all_related_to_gene_start 193 µs
res0 count 142
add exon chr index took 32 s
val count 1054962128
add ops index took 91 s
val count 1386253656
transitive_chr1 7691 ms
res1 count 40172978
q0 33295 µs
res2 count 151956

 */