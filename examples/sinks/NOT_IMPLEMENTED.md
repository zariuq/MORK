# O-Wrapper (Sink) Features NOT Implemented

Based on analysis of MORK source code and documentation, here are sink operations mentioned in Otter's documentation that are **NOT yet implemented**:

## ACT Writing - `(O (ACT <file> <template>))` (NOT IMPLEMENTED)

**From Otter's docs:**
> `(O (ACT <filename> <template>))` should serialize expressions to an ACT file

**Status:** ❌ **Not implemented in MM2**

**Evidence:** See `examples/misc/ACT_README.md`:
> "The MM2 syntax `(O (ACT "filename" <template>))` is **not yet implemented** in MORK."

**Current workaround:**
ACT files must be created via:
1. **Shell scripts:** See `examples/misc/make_fruit_colors_act.sh`
2. **Rust API:** `space.backup_tree(format!("{ACT_PATH}filename.act"))`

**Example (theoretical, will not work):**
```metta
; This DOESN'T WORK - ACT writing not implemented!
(color apple red)
(color banana yellow)

(exec 0
  (, (color $fruit $c))
  (O (ACT my_colors (color $fruit $c))))  ; Would write to my_colors.act
```

**Working alternative (shell script):**
```bash
#!/bin/bash
# make_colors_act.sh
cat > /tmp/colors_space.mm2 << 'EOF'
(color apple red)
(color banana yellow)
EOF

mork run --steps 0 /tmp/colors_space.mm2
# Then use Rust API to call backup_tree()
```

## PATHS - `(O (PATHS <file> <template>))` (NOT IMPLEMENTED)

**From Otter's docs:**
> `(O (PATHS <filename> <template>))` should dump path information to a file

**Status:** ❌ **Not implemented**

**What it might do (if implemented):**
- Dump internal path representations to file
- Useful for debugging PathMap structure
- Export space topology information

**Example (theoretical):**
```metta
; Doesn't work
(exec 0
  (, (color $fruit $c))
  (O (PATHS debug.txt (color $fruit $c))))
```

## COUNT - ✅ **IMPLEMENTED!**

**Status:** ✅ **Fully implemented and working**

**Syntax:** `(count <report-template-with-$k> $k <template-to-count>)`

**Working example:**
```metta
; This WORKS!
(exec 0
  (, (color $fruit $c))
  (O (count (total-colors $k) $k (result $fruit $c))))
; Produces: (total-colors 5)
```

**See test files:**
- `test_count_simple.mm2` - Basic counting
- `test_count_var.mm2` - Using $k in report template
- `test_count_guard.mm2` - Guard mode (filter by expected count)

**Implementation:** `kernel/src/sinks.rs:267` - `CountSink` struct

## SUM - `(O (SUM ...))` (NOT IMPLEMENTED)

**Status:** ❌ **Not checked, likely not implemented**

**What it might do:**
- Sum numeric values
- Aggregate computations

**Example (theoretical):**
```metta
; Theoretical usage
(exec 0
  (, (price $item $cost))
  (O (SUM total-cost $cost)))
```

## What IS Implemented

For reference, these O-wrapper operations **do work**:

✅ **`(+ <template>)`** - Add facts to the space
✅ **`(- <template>)`** - Remove facts from the space (exact match)

See the test files in this directory for working examples!

## Future Work

If you want to implement ACT writing or other sink features, you would need to:

1. **Examine `kernel/src/sinks.rs`** (or equivalent) for sink implementations
2. **Add new sink types** to the sink enum
3. **Add parsing** for the new sink syntax
4. **Implement the sink operation** (file I/O, aggregation, etc.)
5. **Add tests** in `kernel/src/main.rs`

The pattern to follow is in the existing `AddSink` and `RemoveSink` implementations.

## See Also

- **`examples/i_source_tests/NOT_IMPLEMENTED.md`** - I-wrapper features not implemented
- **`examples/misc/ACT_README.md`** - Current ACT file limitations
- **`examples/misc/CREATE_ACT_FILES.md`** - How to create ACT files today
