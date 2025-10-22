# Fibonacci up to F₁₉ = 4181

## Files Created

- **fibonacci_table_19.mm2** - Non-memoized version
- **fibonacci_table_memoized_19.mm2** - Memoized version
- Size: 108MB each (6.2M+ expressions loaded)

## Fibonacci Sequence
```
F0=0, F1=1, F2=1, F3=2, F4=3, F5=5, F6=8, F7=13, F8=21, F9=34, F10=55
F11=89, F12=144, F13=233, F14=377, F15=610, F16=987, F17=1597, F18=2584, F19=4181
```

## Test Cases
Both files test the same values:
```metta
(fib 8 test-8)    ; Expected: 21
(fib 10 test-10)  ; Expected: 55
(fib 12 test-12)  ; Expected: 144
(fib 19 test-19)  ; Expected: 4181
```

## Performance Results

### Non-Memoized Version
```
Loaded: 6,265,690 expressions
Time: 8295 ms
Unifications: 25,359
Writes: 87,194
Transitions: 120,843,026
```

### Memoized Version
```
Loaded: 6,265,691 expressions
Time: 1185 ms (↓ 86% - 7× FASTER!)
Unifications: 6,499 (↓ 74%)
Writes: 27,542 (↓ 68%)
Transitions: 13,628,132 (↓ 89%)

Cache hits: 9 (fib values 0-8 reused)
```

## Key Insights

**Massive speedup for larger Fibonacci numbers!**
- fib(10): 30% faster with memoization
- fib(19): **86% faster with memoization (7× speedup!)**

The performance gap widens dramatically as N increases, because:
1. More intermediate values get reused
2. Exponential redundancy without caching
3. fib(19) requires computing fib(18) + fib(17), which both need many overlapping subproblems

## Table Sizes

**Successor table:** 4,182 entries (S 0 1) through (S 4181 4182)

**Addition table:** 6,261,499 entries
- Covers A(a, b, sum) for a, b ∈ [0, 2599]
- Maximum sum: 4200

Yes, the pre-loading is space-intensive, but the pattern matching on tables makes the actual computation very clean!

## Usage

```bash
cd /home/zar/claude/hyperon/MORK/kernel

# Non-memoized
../target/release/mork run --steps 5000 ../examples/fibonacci_table_19.mm2

# Memoized (much faster!)
../target/release/mork run --steps 5000 ../examples/fibonacci_table_memoized_19.mm2

# Compare performance
bash /tmp/compare_fib_19.sh
```

## Results
```
✓ (fib-result test-8 21)
✓ (fib-result test-10 55)
✓ (fib-result test-12 144)
✓ (fib-result test-19 4181)
```

All tests pass in both versions!
