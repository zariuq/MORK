# I-Wrapper Tests

Complete test suite demonstrating all I-wrapper patterns found in MORK.

## Basic Patterns

### `test_basic_query.mm2`
Basic pattern matching using comma-wrapper.
```bash
mork run --steps 5 test_basic_query.mm2
```
**Demonstrates:** `(, (pattern))` for querying main space

### `test_add_remove.mm2`
Add and remove operations in O-wrapper.
```bash
mork run --steps 5 test_add_remove.mm2
```
**Demonstrates:** `(O (+ ...) (- ...))` for adding/removing facts

### `test_multi_pattern.mm2`
Conjunction of multiple patterns.
```bash
mork run --steps 5 test_multi_pattern.mm2
```
**Demonstrates:** `(, (pattern1) (pattern2) (pattern3))` for conjunctions

## BTM (Main Space) Patterns

### `test_btm_cross_match.mm2`
Cross-matching two patterns from main space.
```bash
mork run --steps 10 test_btm_cross_match.mm2
```
**Demonstrates:** `(I (BTM <pat1>) (BTM <pat2>))` for finding consistent bindings

**Example:**
```metta
(Something (foo bar) (foo bar))
(Else (foo bar) (foo bar))

(I (BTM (Something $x $y))
   (BTM (Else $x $y)))

; Output: (MATCHED (foo bar) (foo bar))
```

## Comparison Operators

### `test_bind_var.mm2`
How `==` constructs and binds expressions.
```bash
mork run --steps 10 test_bind_var.mm2
```
**Demonstrates:** `(I (BTM ...) (== <expr> $var))` for expression construction

**Key insight:** `==` doesn't "check equality" - it **builds** expressions!

### `test_bind_only.mm2`
Using `==` alone for pure construction.
```bash
mork run --steps 5 test_bind_only.mm2
```
**Demonstrates:** `(I (== <expr> $var))` without other sources

### `test_not_equal.mm2`
How `!=` generates all-different pairs.
```bash
mork run --steps 20 test_not_equal.mm2
```
**Demonstrates:** `(I (!= <expr1> <expr2>))` for all-different enumeration

**Example:**
```metta
(VAL X) (VAL Y) (VAL Z)

(I (!= (VAL $x) (VAL $y)))

; Output: All 6 different pairs
; (OUT (X != Y)), (OUT (X != Z)), (OUT (Y != X))
; (OUT (Y != Z)), (OUT (Z != X)), (OUT (Z != Y))
```

### `test_not_equal_constraint.mm2`
Using multiple `!=` for constraint satisfaction.
```bash
mork run --steps 20 test_not_equal_constraint.mm2
```
**Demonstrates:** Multiple `!=` constraints for all-different requirements

**Example:**
```metta
(I (!= (color $c1) (color $c2))
   (!= (color $c2) (color $c3))
   (!= (color $c3) (color $c1)))

; Output: All 6 permutations of 3 colors
```

## ACT (External File) Patterns

**Prerequisites:** Create ACT files first:
```bash
../misc/make_fruit_colors_act.sh
```

### `test_act_cross_match.mm2`
Cross-matching two ACT queries.
```bash
mork run --steps 10 test_act_cross_match.mm2
```
**Demonstrates:** `(I (ACT file <pat1>) (ACT file <pat2>))` for cross-matching external data

### `test_btm_act_mix.mm2`
Mixing BTM (main space) and ACT (external file).
```bash
mork run --steps 10 test_btm_act_mix.mm2
```
**Demonstrates:** `(I (BTM <pat>) (ACT file <pat>))` for joining internal and external data

**Example:**
```metta
; Main space
(likes Alice apple)

; ACT file
(color apple red)

(I (BTM (likes $person $fruit))
   (ACT fruit_colors (color $fruit $c)))

; Output: (preference Alice likes apple which-is red)
```

**Key insight:** Like a SQL JOIN between in-memory and external data!

## Pattern Summary

All I-wrapper patterns found in MORK codebase:

| Pattern | File | Purpose |
|---------|------|---------|
| `(, (pattern))` | test_basic_query.mm2 | Basic query |
| `(I (BTM ...) (BTM ...))` | test_btm_cross_match.mm2 | BTM cross-match |
| `(I (BTM ...) (ACT ...))` | test_btm_act_mix.mm2 | BTM + ACT join |
| `(I (BTM ...) (== ...))` | test_bind_var.mm2 | Query + construct |
| `(I (== ...))` | test_bind_only.mm2 | Pure construction |
| `(I (!= ...))` | test_not_equal.mm2 | All-different |
| `(I (ACT ...) (ACT ...))` | test_act_cross_match.mm2 | ACT cross-match |

## Key Concepts

**BTM** = "BoTtoM" = Main space (the PathMap in MORK)

**ACT** = External serialized space fragment (file-based)

**`==`** = Expression constructor/unifier (not equality check!)

**`!=`** = All-different generator (not inequality check!)

## What's NOT Implemented

- **Negative patterns:** `(I (- <pattern>))` - Not supported
- **ACT writing from MM2:** Can only read ACT files, not write
- **PATHS/COUNT/SUM in I-wrapper:** These are O-wrapper operations

## See Also

- `../misc/ACT_*.md` - ACT file documentation
- `../README.md` - Main examples documentation
