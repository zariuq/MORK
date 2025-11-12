# MORK Aggregator Tests

This directory contains consolidated tests for MORK's aggregator sinks (max, min, count).

## Core Files

### Working Examples
- **aggregators_core.mm2** - Basic MAX, MIN, COUNT usage
- **aggregators_practical.mm2** - Real-world use cases
- **aggregators_advanced.mm2** - Edge cases (ties, negatives, single values)

### Bug Documentation
- **negative_test_debruijn_bug.mm2** - Demonstrates known de Bruijn index bug
- **BUG_MULTIPLE_AGGREGATORS.mm2** - Comprehensive bug test cases
- **BUG_ANALYSIS_DEBRUIJN.md** - Complete technical analysis

### Documentation
- **README_CONSOLIDATED.md** - Migration guide from old tests
- **CONSOLIDATION_SUMMARY.md** - Organization notes
- **VARIABLE_BINDING_GUIDE.md** - How variables work in aggregators

## Known Bug: Multiple Aggregators with Different Variables

When using multiple aggregators in one exec with DIFFERENT variables, only the first works:

```metta
; BROKEN
(O (max (r $m) $m X)
   (min (s $n) $n Y))
; Output: (r 42), (s $a)  ← Second fails!

; WORKS
(O (max (r $m) $m X)
   (min (s $m) $m Y))
; Output: (r 42), (s 42)  ← Both work!
```

**Workaround**: Use same variable for all aggregators, or split into separate execs.

See `negative_test_debruijn_bug.mm2` and `BUG_ANALYSIS_DEBRUIJN.md` for details.

## Archived Tests

50 original test files have been moved to `archive/` for reference.

## Running Tests

```bash
# Run core examples (all should work)
./target/release/mork run --steps 10 examples/sinks/aggregators_core.mm2

# Run negative test (demonstrates bugs)
./target/release/mork run --steps 10 examples/sinks/negative_test_debruijn_bug.mm2
```
