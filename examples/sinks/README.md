# O-Wrapper (Sink) Tests

Complete test suite for **output operations** (sinks) in MM2.

## What are Sinks?

Sinks are **output operations** in the `(O ...)` wrapper that modify the space or external storage:

- `(+ <template>)` - **ADD** to the space
- `(- <template>)` - **REMOVE** from the space (exact match)
- `(count <report> $k <template>)` - **COUNT** aggregator (count unique instances)
- `(ACT <file> <template>)` - **SERIALIZE** to ACT file (NOT IMPLEMENTED)

## Test Files

### Basic ADD/REMOVE Operations

**test_add_simple.mm2**
- Demonstrates `(O (+ <fact>))` to add facts
- Adds three color facts to the space
```bash
mork run --steps 10 test_add_simple.mm2
# Result: (color apple red), (color banana yellow), (color grape purple)
```

**test_remove_simple.mm2**
- Demonstrates `(O (- <fact>))` to remove facts
- Shows exact match semantics
```bash
mork run --steps 10 test_remove_simple.mm2
# Result: (color grape purple) remains, others removed
```

**test_remove_pattern.mm2**
- Shows how `(- <pattern>)` uses bound variable values
- Pattern matching in input, exact removal in output
```bash
mork run --steps 10 test_remove_pattern.mm2
# Result: (removed buy-milk), other tasks remain
```

**test_add_remove.mm2**
- Combined example showing both operations
- Marks tasks as completed and removes them
```bash
mork run --steps 10 test_add_remove.mm2
# Result: One task completed and removed
```

### COUNT Aggregator

**test_count_simple.mm2**
- Basic COUNT usage to count matching patterns
- Counts 5 color facts
```bash
mork run --steps 10 test_count_simple.mm2
# Result: (total-colors 5)
```

**test_count_var.mm2**
- Using $k variable in the report template
- Counts 18 combinations (3×2×3)
```bash
mork run --steps 10 test_count_var.mm2
# Result: (all 18)
```

**test_count_guard.mm2**
- Using a constant as guard to filter by expected count
- Only outputs if count matches the guard value
```bash
mork run --steps 10 test_count_guard.mm2
# Result: (all-eighteen) but NOT (all-sixteen)
```

**test_count_explained.mm2**
- Detailed walkthrough of how COUNT works internally
- Shows template instantiation step-by-step
```bash
mork run --steps 10 test_count_explained.mm2
# Result: (total 3) plus the matched facts
```

**test_count_template_meaning.mm2**
- Demonstrates that template symbols (cux, result, etc.) are arbitrary
- All three tests with different symbols produce the same count
```bash
mork run --steps 10 test_count_template_meaning.mm2
# Result: (with-cux 4), (with-result 4), (with-arbitrary 4)
```

## Key Concepts

### ADD: `(+ <template>)`
Adds a new fact to the space. Variables in the template use their bound values.

```metta
(exec 0
  (, (trigger))
  (O (+ (new-fact 42))
     (+ (another-fact xyz))))
```

### REMOVE: `(- <template>)`
Removes facts with **exact match semantics**. If you use variables, they must be bound from the input pattern.

```metta
(exec 0
  (, (task $name))
  (O (- (task $name))))  ; Removes the matched task
```

**Important:** `(- <pattern>)` is NOT pattern matching! It only removes exact matches.

### COUNT: `(count <report-template-with-$k> $k <template-to-count>)`

Counts **unique instantiations** of a template and binds the count to variable `$k`.

**How it works:**
```metta
(exec 0
  (, (color $fruit $c))
  (O (count (total $k) $k (result $fruit $c))))
```

**Step by step:**
1. Input pattern `(color $fruit $c)` matches multiple times with different bindings
2. For each match, COUNT instantiates the template `(result $fruit $c)`:
   - `(result apple red)`
   - `(result banana yellow)`
   - `(result grape purple)`
3. COUNT collects all **unique** instantiations in an internal PathMap
4. COUNT returns the number of unique instantiations (e.g., 3)
5. Report template `(total $k)` becomes `(total 3)`

**Important:**
- The template symbol (`result`, `cux`, `counted`, etc.) is **arbitrary** - use any symbol you want!
- Template instantiations are **never added to the space** - they only exist inside COUNT
- COUNT counts **unique combinations** of the bound variables

**Using $k in the report:**
```metta
(O (count (all $k) $k (cux $z $y $x)))
```
The variable `$k` is bound to the count and can appear in the report template: `(all 18)`.

**Guard mode (constant instead of $k):**
```metta
(O (count (all-eighteen) 18 (cux $z $y $x)))
```
Only produces output if the count equals 18. Acts as a filter/assertion.

**Three modes:**
1. **Variable in report**: `(count (total $k) $k ...)` → `(total 18)`
2. **Guard/filter**: `(count (only-if-18) 18 ...)` → `(only-if-18)` or nothing
3. **Constant in report**: `(count (found) $k ...)` → `(found)` with count checked

### ACT Writing: NOT IMPLEMENTED

The syntax `(O (ACT <file> <template>))` is documented but **not yet implemented** in MORK.

**Workaround:** Create ACT files using:
1. Shell scripts (see `examples/misc/make_fruit_colors_act.sh`)
2. Rust API: `space.backup_tree("filename.act")`

See `NOT_IMPLEMENTED.md` for details.

## NOT_IMPLEMENTED.md

Documents sink operations that are mentioned in Otter's docs but not yet implemented:
- `(O (ACT <file> <template>))` - Write to ACT file
- `(O (PATHS <file> <template>))` - Dump paths
- `(O (COUNT ...))` - Count operations
- `(O (SUM ...))` - Sum operations

## Running Tests

```bash
cd examples/sinks

# Test ADD
mork run --steps 10 test_add_simple.mm2

# Test REMOVE
mork run --steps 10 test_remove_simple.mm2

# Test pattern-based removal
mork run --steps 10 test_remove_pattern.mm2

# Test combined ADD/REMOVE
mork run --steps 10 test_add_remove.mm2
```

## See Also

- **`../i_source_tests/`** - I-wrapper (input) patterns
- **`../misc/ACT_README.md`** - How to create ACT files
- **`kernel/src/sinks.rs`** - Sink implementation (if exists)
- **`kernel/src/space.rs`** - Space modification operations
