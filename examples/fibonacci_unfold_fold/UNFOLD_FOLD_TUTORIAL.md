# Unfold/Fold Recursion Pattern in MM2

## Overview

Recursive functions in functional programming follow an **unfold/fold** pattern:
1. **UNFOLD**: Break the problem into smaller subproblems
2. **BASE CASES**: Stop recursion when reaching the simplest cases
3. **FOLD**: Combine results from subproblems back together

## Example: Fibonacci in MM2

### Conceptual Flow for fib(3)
```
fib(3)
  ↓ UNFOLD
fib(2) + fib(1)
  ↓ UNFOLD        ↓ BASE
fib(1) + fib(0)   1
  ↓ BASE  ↓ BASE
  1       0
  ↓ FOLD
  1
  ↓ FOLD (final)
  2 ✓
```

### MM2 Implementation

**Key Pattern**: Use broadcast chain to keep exec rules persistent

```metta
; Broadcast chain - extracts execs from step declarations
(exec bc
  (, ((step ($priority $name)) $premises0 $conclusions0)
     (exec bc $premises1 $conclusions1))
  (O (+ (exec ($priority $name) $premises0 $conclusions0))
     (+ (exec bc $premises1 $conclusions1))
     (- ((step ($priority $name)) $premises0 $conclusions0))
     (+ ((step ((S ($priority)) $name)) $premises0 $conclusions0))))
```

**Why needed?** MM2 exec rules are consumed after firing, even if they don't match anything. The broadcast chain pattern keeps regenerating them.

### Rule Structure

```metta
; BASE CASE
((step (0 fib-base-0))
  (, (fib 0 $id))
  (O (+ (fib-result $id 0))
     (- (fib 0 $id))))

; UNFOLD: Break into subproblems
((step (2 fib-2-unfold))
  (, (fib 2 $id))
  (O (+ (fib 1 (left-2 $id)))      ; Create subproblem 1
     (+ (fib 0 (right-2 $id)))      ; Create subproblem 2
     (+ (wait-2 $id))                ; Create wait state
     (- (fib 2 $id))))               ; Remove original

; FOLD: Combine results
((step (3 fib-2-fold))
  (, (wait-2 $id)
     (fib-result (left-2 $id) 1)    ; Got result from left
     (fib-result (right-2 $id) 0))  ; Got result from right
  (O (+ (fib-result $id 1))         ; Produce combined result
     (- (wait-2 $id))                ; Clean up
     (- (fib-result (left-2 $id) 1))
     (- (fib-result (right-2 $id) 0))))
```

## Running the Examples

```bash
# Simple unfold/fold (no recursion)
../target/release/mork run --steps 100 ../../metamath/mmverify/mm2/unfold_fold_simple.mm2

# Fibonacci fib(2) = 1
../target/release/mork run --steps 100 ../../metamath/mmverify/mm2/fibonacci_clean.mm2

# Fibonacci fib(3) = 2 (full recursion)
../target/release/mork run --steps 200 ../../metamath/mmverify/mm2/fibonacci_complete.mm2
```

## Results

```
# fib(2)
executing 100 steps took 10 ms (unifications 88, writes 348, transitions 8576)
(fib-result main 1) ✓

# fib(3)
executing 200 steps took 20 ms (unifications 181, writes 718, transitions 18212)
(fib-result main 2) ✓
```

## Key Lessons

1. **Exec rules are consumed** even when they don't match
2. **Broadcast chain pattern** keeps rules persistent
3. **Wait states** (like `wait-2`) coordinate fold phase
4. **Priority ordering** matters but is complex (see MM2 wiki)
5. Unfold/fold maps naturally to MM2's rewrite semantics
6. **Bidirectional pattern matching** - tables work both ways:
   - S(4, $s) → $s=5 (successor)
   - S($p, 5) → $p=4 (predecessor)
   - No need for separate predecessor table!

## Files

- `unfold_fold_simple.mm2` - Basic pattern without recursion
- `fibonacci_clean.mm2` - fib(2) with full recursion
- `fibonacci_complete.mm2` - fib(3) with documentation
- `fibonacci_table.mm2` - Generic version using arithmetic tables (up to fib(12)=144)
