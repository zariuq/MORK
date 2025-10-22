# Memoization in MM2 Fibonacci

## Problem: Redundant Computation

Without memoization, recursive Fibonacci has exponential time complexity because it recomputes the same values many times:

```
fib(5)
  calls fib(4) and fib(3)
    fib(4) calls fib(3) and fib(2)  ← fib(3) computed AGAIN
    fib(3) calls fib(2) and fib(1)  ← fib(2) computed AGAIN
```

For fib(10), this results in **hundreds of redundant computations**.

## Solution: Global Cache

The memoized version uses a global cache to store computed values:

```metta
; Cache structure: (fib-value $n $result)
; Examples:
(fib-value 0 0)
(fib-value 1 1)
(fib-value 2 1)
(fib-value 5 5)
```

### Priority-Based Cache Checking

**Priority 0 (HIGHEST): Check cache first**
```metta
((step (0 fib-from-cache))
  (, (fib $n $id)
     (fib-value $n $val))    ; ← Look for cached value
  (O (+ (fib-result $id $val))
     (+ (cache-hit $n))       ; ← Record cache hit
     (- (fib $n $id))))
```

**Priority 1-2: Base cases populate cache**
```metta
((step (1 fib-base-0))
  (, (fib 0 $id))
  (O (+ (fib-result $id 0))
     (+ (fib-value 0 0))      ; ← Add to cache
     (- (fib 0 $id))))
```

**Priority 3: Unfold (only if not cached)**
```metta
((step (3 fib-unfold))
  (, (fib $n $id)
     (S $n1 $n)
     (S $n2 $n1))
  (O (+ (fib $n1 (left $id)))
     (+ (fib $n2 (right $id)))
     (+ (wait-add $id))
     (+ (computing $n $id))   ; ← Track what we're computing
     (- (fib $n $id))))
```

**Priority 4: Fold and ADD TO CACHE**
```metta
((step (4 fib-fold))
  (, (wait-add $id)
     (fib-result (left $id) $r1)
     (fib-result (right $id) $r2)
     (A $r1 $r2 $sum)
     (computing $n $id))      ; ← Know which n this is for
  (O (+ (fib-result $id $sum))
     (+ (fib-value $n $sum))  ; ← Add computed value to cache
     (- (wait-add $id))
     (- (computing $n $id))
     (- (fib-result (left $id) $r1))
     (- (fib-result (right $id) $r2))))
```

## Key Insights

1. **Priority ordering is crucial**
   - Cache check must be priority 0 (highest) to fire before unfold
   - Without this, unfold would happen first and we'd lose memoization benefit

2. **Track computation context**
   - `(computing $n $id)` links the computation ID to the original n value
   - This allows fold rule to know which cache entry to create

3. **Cache hits prevent unfold**
   - When cache hit fires, it removes `(fib $n $id)` fact
   - This prevents unfold rule from firing for that value

## Performance Results: fib(10)

### Without Memoization
```
executing 715 steps took 86 ms
  unifications: 1632
  writes: 5816
  transitions: 307720
```

### With Memoization
```
executing 855 steps took 60 ms
  unifications: 903    (↓ 45%)
  writes: 3739         (↓ 36%)
  transitions: 126421  (↓ 59%)

Cache hits: 5
Cached values: fib(0) through fib(10)
```

### Speedup: 30% faster! ⚡

## Files

- **fibonacci_memoized.mm2** - Small example with fib(5)
- **fibonacci_table_memoized.mm2** - Full version with tables up to 144
- **fibonacci_table.mm2** - Non-memoized version for comparison

## Usage

```bash
# Run memoized version
cd /home/zar/claude/hyperon/MORK/kernel
../target/release/mork run --steps 1000 ../examples/fibonacci_table_memoized.mm2

# Compare performance
bash /tmp/compare_fib_10.sh
```

## Pattern: Memoization in MM2

This pattern can be applied to any recursive computation:

1. Define cache structure: `(cache-name $input $output)`
2. Add priority-0 rule to check cache
3. Add `(computing $input $id)` fact during unfold
4. Populate cache in fold rule using `(computing $input $id)` context
5. Record `(cache-hit $input)` for diagnostics

The key innovation is using `(computing $n $id)` to link the abstract computation ID to the concrete input value, allowing the fold rule to create the right cache entry.
