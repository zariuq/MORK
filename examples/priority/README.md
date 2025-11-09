# Priority Ordering in MM2/MORK

This directory contains tests exploring how MM2's priority system works and how it differs from standard lexicographic ordering.

## Key Findings

### 1. Basic Principle: SHORTLEX Ordering

MM2 uses **shortlex ordering** (short lexicographic), NOT standard lexicographic ordering.

**Shortlex**: Compare by length first, then by value
**Lexicographic**: Compare character-by-character (dictionary style)

Example difference:
```
Shortlex:      0, 1, 2, 3, 10, 20, 100   (length first)
Lexicographic: 0, 1, 10, 100, 2, 20, 3   (digit by digit)
```

### 2. Numeric Priorities

Lower numbers = higher priority (execute first)

```metta
(exec 0 ...) ; Executes first
(exec 1 ...) ; Executes second
(exec 2 ...) ; Executes third
```

Shortlex applies to numbers:
- Single digits (0-9) execute before double digits (10-99)
- Double digits execute before triple digits (100-999)
- Within same length, numeric value determines order

See: `test_priority_basic.mm2`, `test_priority_shortlex.mm2`

### 3. Symbol Priorities

Symbols also use shortlex ordering:
- Shorter symbols execute first
- Within same length, alphabetical ordering applies

```metta
; Execution order:
(exec a ...)   ; 1 char
(exec b ...)   ; 1 char (b > a)
(exec aa ...)  ; 2 chars
(exec bb ...)  ; 2 chars (bb > aa)
(exec aaa ...) ; 3 chars
```

See: `test_priority_symbols.mm2`

### 4. Numbers vs Symbols

When priorities have the same byte length, **numbers come before symbols**:

```metta
; Execution order with 2-byte priorities:
(exec 42 ...) ; Number, 2 digits
(exec 99 ...) ; Number, 2 digits
(exec ab ...) ; Symbol, 2 chars - comes AFTER numbers
```

See: `test_priority_boundaries.mm2`

### 5. Structured Priorities

Complex expressions can be used as priorities. They're compared by their byte representation:

```metta
; Execution order:
(exec (0 0) ...)   ; Tuple comparison
(exec (0 1) ...)   ; Element-wise: (0 0) < (0 1)
(exec (1 0) ...)   ; Element-wise: (0 1) < (1 0)

; Also works:
(exec (phase 0) ...)
(exec (phase 1) ...)
(exec (phase 2) ...)
```

The tuples are compared element by element, applying shortlex recursively.

See: `test_priority_structured.mm2`

### 6. Nested Structures

Nesting depth affects the byte representation and thus ordering:

```metta
(exec ((0))   ...) ; 2 levels of nesting
(exec (0)     ...) ; 1 level of nesting
(exec (((0))) ...) ; 3 levels of nesting
```

The execution order depends on the internal byte encoding of these structures.

See: `test_priority_boundaries.mm2`

### 7. Unsupported Priority Values

**Negative numbers**: Cause panic (unreachable code)
**Empty tuples ()**: Cause panic (unreachable code)

These are not valid priority values in MM2.

## Test Files

| File | Purpose |
|------|---------|
| `test_priority_basic.mm2` | Simple numeric priorities (0, 1, 2, 3) |
| `test_priority_shortlex.mm2` | Demonstrates shortlex: 0, 1, 2, 3, 10, 20, 100 |
| `test_priority_symbols.mm2` | Symbol ordering: a, b, aa, bb, aaa |
| `test_priority_structured.mm2` | Tuples and complex priorities |
| `test_priority_boundaries.mm2` | Edge cases: large numbers, nesting, numbers vs symbols |

## Running Tests

From the MORK root directory:

```bash
# Basic numeric test
./target/release/mork run --steps 10 examples/priority/test_priority_basic.mm2

# Shortlex demonstration
./target/release/mork run --steps 20 examples/priority/test_priority_shortlex.mm2

# Symbol ordering
./target/release/mork run --steps 20 examples/priority/test_priority_symbols.mm2

# Structured priorities
./target/release/mork run --steps 20 examples/priority/test_priority_structured.mm2

# Boundary cases
./target/release/mork run --steps 20 examples/priority/test_priority_boundaries.mm2
```

## Implementation Details

Priority ordering is implemented in the PathMap data structure, which stores atoms in sorted order. When exec rules are selected for execution, they're retrieved via `to_next_val()` which walks the PathMap in its natural byte-order, giving shortlex behavior.

See: `kernel/src/space.rs:1651-1680` (metta_calculus function)

## Practical Usage

### Phased Execution

Use numeric priorities for clear phase ordering:

```metta
; Phase 1: Initialization
(exec 0 ...)

; Phase 2: Processing
(exec 1 ...)

; Phase 3: Cleanup
(exec 2 ...)
```

### Named Phases

Use structured priorities for self-documenting code:

```metta
(exec (init 0) ...)
(exec (process 1) ...)
(exec (cleanup 2) ...)
```

### Fine-Grained Control

Use shortlex to control execution order within phases:

```metta
; All execute in phase 1, but in specific order:
(exec 1 ...)   ; First
(exec 2 ...)   ; Second
(exec 10 ...)  ; Third (after all single digits)
(exec 11 ...)  ; Fourth
```

## Related Documentation

- `HOWTO_MM2.md` - Basic MM2 syntax and exec rules
- `examples/sources/` - I-wrapper (input) patterns
- `examples/sinks/` - O-wrapper (output) operations
