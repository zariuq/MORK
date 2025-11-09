# MM2 Examples

Collection of MM2 (Minimal MeTTa 2) example programs demonstrating various patterns and techniques.

## Documentation

- **`../HOWTO_MM2.md`** - Complete MM2 tutorial and reference guide
- **`docs/LIST_REPRESENTATIONS.md`** - Indexed lists vs Cons/Nil comparison
- **`docs/MM2_LIST_PATTERNS.md`** - List processing patterns

## Directory Structure

- **`i_wrapper_tests/`** - Complete test suite for all I-wrapper (input) patterns (BTM, ACT, ==, !=)
- **`sinks/`** - Complete test suite for O-wrapper (output/sink) patterns (+, -)
- **`fibonacci_unfold_fold/`** - Recursion, unfold/fold, and memoization examples
- **`map_fold_swap/`** - List processing, loops, and SIMD patterns
- **`if_cond/`** - Conditional branching examples
- **`metamath/`** - Metamath formal proof databases in MM2 format
- **`misc/`** - Miscellaneous patterns and experiments

## Unfold/Fold Recursion Pattern

Examples demonstrating how recursive functions work in MM2 through unfold and fold phases:

### `unfold_fold_demo.mm2`
**Simplest example** - Demonstrates the basic unfold/fold pattern without deep recursion.
```bash
../target/release/mork run --steps 50 unfold_fold_demo.mm2
# Result: (final task 10)
```
Computes `double(5) = 10` by unfolding to `add(5,5)` then folding the result.

### `unfold_fold_simple.mm2`
Non-recursive version with multiple tasks showing parallel computation.
```bash
../target/release/mork run --steps 100 unfold_fold_simple.mm2
# Results: (final-result task1 6), (final-result task2 4), (final-result task3 2)
```

### `fibonacci_clean.mm2`
**Fibonacci with recursion** - Computes fib(2) = 1
```bash
../target/release/mork run --steps 100 fibonacci_clean.mm2
# Result: (fib-result main 1)
```
Shows unfold: `fib(2) → fib(1) + fib(0)` then fold: `1 + 0 = 1`

### `fibonacci_complete.mm2`
**Full Fibonacci example** - Computes fib(3) = 2 with complete documentation
```bash
../target/release/mork run --steps 200 fibonacci_complete.mm2
# Result: (fib-result main 2)
```
Demonstrates nested recursion:
```
fib(3) → fib(2) + fib(1)
  fib(2) → fib(1) + fib(0) → 1
  fib(1) → 1
Result: 2
```

### `fibonacci_table.mm2`
**Generic table-based Fibonacci** - Uses preloaded arithmetic tables
Works for any N within table bounds (up to fib(12) = 144).
```bash
../target/release/mork run --steps 2000 fibonacci_table.mm2
# Results: fib(8)=21, fib(10)=55, fib(12)=144
```
**Key optimization**: Uses only S(n, n+1) bidirectionally via pattern matching:
- `S(4, $s)` matches to get successor: $s=5
- `S($p, 5)` matches to get predecessor: $p=4

No separate predecessor table needed! Addition table A(a, b, sum) handles the fold phase.

### `UNFOLD_FOLD_TUTORIAL.md`
Complete tutorial explaining the unfold/fold pattern with implementation details.

## Memoization Pattern

Demonstrates caching computed values to avoid redundant computation in recursive functions.

### `fibonacci_memoized.mm2`
**Basic memoization** - Small example computing fib(5)
```bash
../target/release/mork run --steps 250 fibonacci_memoized.mm2
# Shows cache hits and cached values
```
Uses priority-0 rule to check cache before computing.

### `fibonacci_table_memoized.mm2`
**Full memoized Fibonacci** - With large tables up to fib(12)
```bash
../target/release/mork run --steps 1000 fibonacci_table_memoized.mm2
# Performance: 30% faster than non-memoized version
```

**Key technique**: Priority-based cache checking
```metta
; Priority 0 (HIGHEST): Check cache first
((step (0 fib-from-cache))
  (, (fib $n $id)
     (fib-value $n $val))
  (O (+ (fib-result $id $val))
     (+ (cache-hit $n))
     (- (fib $n $id))))

; Priority 4: Fold and ADD TO CACHE
((step (4 fib-fold))
  (, (wait-add $id)
     (fib-result (left $id) $r1)
     (fib-result (right $id) $r2)
     (A $r1 $r2 $sum)
     (computing $n $id))      ; Links ID to original n
  (O (+ (fib-result $id $sum))
     (+ (fib-value $n $sum))  ; Cache the result
     ...))
```

**Performance for fib(10):**
- Non-memoized: 1632 unifications, 86ms
- Memoized: 903 unifications (-45%), 60ms (-30%)

### `fibonacci_table_19.mm2` & `fibonacci_table_memoized_19.mm2`
**Large-scale Fibonacci** - Up to fib(19) = 4181
```bash
# Non-memoized (slower)
../target/release/mork run --steps 5000 fibonacci_table_19.mm2

# Memoized (7× faster!)
../target/release/mork run --steps 5000 fibonacci_table_memoized_19.mm2
```

**Files:** 108MB each, 6.2M+ expressions loaded

**Test cases (same in both):**
- fib(8) = 21
- fib(10) = 55
- fib(12) = 144
- fib(19) = 4181

**Performance for fib(19):**
- Non-memoized: 25,359 unifications, 8295ms
- Memoized: 6,499 unifications (-74%), 1185ms (-86% - **7× faster!**)

The performance gap widens dramatically as N increases! See `FIBONACCI_19_RESULTS.md`.

### Documentation

- `MEMOIZATION_SUMMARY.md` - Technical deep-dive on memoization pattern
- `MEMOIZATION_DEMO.md` - User-friendly demonstration
- `FIBONACCI_19_RESULTS.md` - Performance analysis for fib(19)

## I-Wrapper Tests

Complete test suite covering ALL I-wrapper patterns found in MORK!

See `i_wrapper_tests/README.md` for detailed documentation.

**Test coverage:**
- Basic patterns: `(, (pattern))` comma-wrapper queries
- BTM patterns: `(I (BTM ...) (BTM ...))` cross-matching main space
- ACT patterns: `(I (ACT ...) (ACT ...))` cross-matching external files
- Mixed patterns: `(I (BTM ...) (ACT ...))` joining internal + external data
- Comparisons: `(I (== ...))` expression construction, `(I (!= ...))` all-different
- Constraints: Multiple `!=` for constraint satisfaction

**Quick test:**
```bash
cd i_wrapper_tests
mork run --steps 10 test_btm_cross_match.mm2
mork run --steps 20 test_not_equal.mm2
```

## O-Wrapper (Sink) Tests

Complete test suite covering output operations (sinks) in MM2!

See `sinks/README.md` for detailed documentation.

**Test coverage:**
- ADD: `(O (+ <template>))` add facts to the space
- REMOVE: `(O (- <template>))` remove facts with exact match semantics
- COUNT: `(O (count <report> $k <template>))` count aggregator
- Pattern-based removal: Using bound variables from input patterns
- Combined add/remove operations

**Quick test:**
```bash
cd sinks
mork run --steps 10 test_add_simple.mm2
mork run --steps 10 test_count_var.mm2
```

**NOT implemented:**
- `(O (ACT <file> <template>))` - ACT writing (use shell scripts or Rust API)
- `(O (PATHS ...))`, `(O (SUM ...))` - Not yet available

See `sinks/NOT_IMPLEMENTED.md` for workarounds.

## Control Flow Examples

### `if_cond/if_simple_true.mm2`
Basic conditional branching using exec priority ordering.

### `if_cond/if_macro_codegen.mm2`
Macro-based conditional code generation.

### `map_fold_swap/loop_dynamic.mm2`
Dynamic loop construction.

### `misc/loop_batch_codegen.mm2`
Batch code generation for loops.

## List Processing Examples

See `map_fold_swap/` directory for comprehensive list processing examples:
- Sequential processing (map, fold, swap)
- SIMD-style parallel processing
- Loop templates and patterns

### `misc/list_walkers.mm2`
Pattern for walking and processing list structures.

## ACT Files - External Data Sources

ACT files are serialized space fragments that MM2 programs can query without loading into memory. Perfect for large datasets and shared knowledge bases!

**Quick start:**
```bash
# Create ACT file and query it
./examples/misc/make_fruit_colors_act.sh
../target/release/mork run --steps 20 misc/act_demo_query.mm2
```

See `misc/ACT_QUICKSTART.md` for a 30-second demo!

### `misc/act_demo_query.mm2`
Query fruit color facts from an external ACT file.
```bash
../target/release/mork run --steps 50 misc/act_demo_query.mm2
# Result: (result apple-is red), (yellow-fruit banana), (yellow-fruit lemon)
```

### `misc/act_math_demo.mm2`
Use an external addition table as a lookup database.
```bash
../target/release/mork run --steps 20 misc/act_math_demo.mm2
# Result: (result two-plus-two-equals 4), (ways-to-make-four 2 and 2), ...
```

### `misc/act_demo_same_color.mm2`
Two-stage approach for cross-matching ACT data (finds fruits with same color).
```bash
../target/release/mork run --steps 50 misc/act_demo_same_color.mm2
# Result: (same-color apple strawberry red), (same-color banana lemon yellow), ...
```

**Key benefits:**
- Data stays on disk (memory-mapped)
- Fast prefix indexing for efficient queries
- Persistent across program runs
- Can be shared between multiple programs

**Documentation:**
- `misc/ACT_QUICKSTART.md` - **Start here!** 30-second demo
- `misc/ACT_TUTORIAL.md` - Learning guide with hands-on examples
- `misc/ACT_PRACTICE_GUIDE.md` - Exercises and practice patterns
- `misc/ACT_EXAMPLES_SUMMARY.md` - Summary of working patterns
- `misc/ACT_README.md` - Technical reference and implementation details
- `misc/CREATE_ACT_FILES.md` - How to create ACT files

**Scripts:**
- `misc/make_fruit_colors_act.sh` - Create fruit_colors.act
- `misc/make_math_facts_act.sh` - Create math_facts.act

## Metamath Examples

The `metamath/` directory contains Metamath formal proof databases converted to MM2 format:

### `metamath/demo0.mm2` (100 expressions)
Simple demonstration database with basic arithmetic and logic. Proves `|- t = t` (reflexivity).

### `metamath/disjoint2.mm2` (93 expressions)
Demonstrates disjoint variable constraints.

### `metamath/miu.mm2` (95 expressions)
The MIU formal system from Hofstadter's "Gödel, Escher, Bach".

**Loading examples:**
```bash
../target/release/mork run --steps 0 metamath/demo0.mm2
# Result: loaded 100 expressions, dumped successfully
```

See `metamath/README.md` for details on contents and how to work with these files.

## Key Techniques

1. **Broadcast Chain Pattern** - Keep exec rules persistent
   ```metta
   (exec bc
     (, ((step ($priority $name)) $premises0 $conclusions0)
        (exec bc $premises1 $conclusions1))
     (O (+ (exec ($priority $name) $premises0 $conclusions0))
        (+ (exec bc $premises1 $conclusions1))
        (- ((step ($priority $name)) $premises0 $conclusions0))
        (+ ((step ((S ($priority)) $name)) $premises0 $conclusions0))))
   ```

2. **Wait States** - Coordinate multi-step computations
3. **Phase Gates** - Control execution flow through phases

## Resources

- [MM2 Wiki](https://github.com/trueagi-io/MORK/wiki/Minimal-MeTTa-2-(MM2))
- MORK Repository: https://github.com/trueagi-io/MORK
