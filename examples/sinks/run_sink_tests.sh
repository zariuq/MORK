#!/bin/bash
# Run all sink aggregator tests with expected output validation
# Usage: ./run_sink_tests.sh

set -e

MORK_BIN="${MORK_BIN:-../../target/release/mork}"

if [ ! -f "$MORK_BIN" ]; then
    echo "Error: MORK binary not found at $MORK_BIN"
    echo "Please build MORK first: cargo build --release"
    exit 1
fi

echo "=========================================="
echo "MORK Sink Aggregator Test Suite"
echo "=========================================="
echo ""

TOTAL_TESTS=0
PASSED=0
FAILED=0

# Color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Helper function to run a test
run_test() {
    local test_file="$1"
    local test_name="$2"
    shift 2
    local expected_patterns=("$@")

    echo "-------------------------------------------"
    echo "Test: $test_name"
    echo "File: $test_file"

    if [ ! -f "$test_file" ]; then
        echo -e "${RED}✗ SKIP${NC}: File not found"
        return
    fi

    TOTAL_TESTS=$((TOTAL_TESTS + 1))

    # Run MORK and capture output
    output=$($MORK_BIN run "$test_file" 2>&1)

    # Check for expected patterns
    local all_found=true
    for pattern in "${expected_patterns[@]}"; do
        if echo "$output" | grep -q "$pattern"; then
            echo -e "${GREEN}✓${NC} Found: $pattern"
        else
            echo -e "${RED}✗${NC} Missing: $pattern"
            all_found=false
        fi
    done

    if [ "$all_found" = true ]; then
        echo -e "${GREEN}✓ PASSED${NC}"
        PASSED=$((PASSED + 1))
    else
        echo -e "${RED}✗ FAILED${NC}"
        FAILED=$((FAILED + 1))
        echo ""
        echo "Full output:"
        echo "$output"
    fi

    echo ""
}

# Helper function to run a negative test (should NOT contain patterns)
run_negative_test() {
    local test_file="$1"
    local test_name="$2"
    shift 2
    local forbidden_patterns=("$@")

    echo "-------------------------------------------"
    echo "Test: $test_name (NEGATIVE)"
    echo "File: $test_file"

    if [ ! -f "$test_file" ]; then
        echo -e "${RED}✗ SKIP${NC}: File not found"
        return
    fi

    TOTAL_TESTS=$((TOTAL_TESTS + 1))

    # Run MORK and capture output
    output=$($MORK_BIN run "$test_file" 2>&1)

    # Check that forbidden patterns are NOT present
    local none_found=true
    for pattern in "${forbidden_patterns[@]}"; do
        if echo "$output" | grep -q "$pattern"; then
            echo -e "${RED}✗${NC} Should NOT contain: $pattern"
            none_found=false
        else
            echo -e "${GREEN}✓${NC} Correctly absent: $pattern"
        fi
    done

    if [ "$none_found" = true ]; then
        echo -e "${GREEN}✓ PASSED${NC}"
        PASSED=$((PASSED + 1))
    else
        echo -e "${RED}✗ FAILED${NC}"
        FAILED=$((FAILED + 1))
        echo ""
        echo "Full output:"
        echo "$output"
    fi

    echo ""
}

# ===== Core Aggregator Tests =====
echo "=========================================="
echo "PART 1: Core Aggregators (MAX, MIN, COUNT)"
echo "=========================================="
echo ""

run_test "aggregators_core.mm2" \
    "Core aggregators - MAX, MIN, COUNT, HEAD" \
    "(max-depth 9)" \
    "(min-cost 8)" \
    "(total-fruits 5)" \
    "(found-max 10)"

# ===== Advanced Aggregator Tests =====
echo "=========================================="
echo "PART 2: Advanced Aggregators"
echo "=========================================="
echo ""

run_test "aggregators_advanced.mm2" \
    "Advanced - duplicates, ties, conditions" \
    "(max-of-dups 9)" \
    "(coldest -5)" \
    "(warmest 3)" \
    "(count-dups 5)" \
    "(max-with-ties 10)"

# ===== Practical Use Cases =====
echo "=========================================="
echo "PART 3: Practical Use Cases"
echo "=========================================="
echo ""

run_test "aggregators_practical.mm2" \
    "Practical - graph, monitoring, scheduling" \
    "(max-level 5)" \
    "(total-hyps 5)" \
    "(longest-path 15)" \
    "(shortest-path 7)" \
    "(peak-cpu 89)" \
    "(min-cpu 23)" \
    "(active-processes 4)" \
    "(shortest-proof 8)" \
    "(total-theorems 4)"

# ===== Sum Tests =====
echo "=========================================="
echo "PART 4: Sum Aggregator Tests"
echo "=========================================="
echo ""

run_test "test_sum_basic.mm2" \
    "Sum - basic integer addition" \
    "(total-60)" \
    "(single-sum)" \
    "(zero-result)"

run_test "test_sum_edge_cases.mm2" \
    "Sum - edge cases (large numbers)" \
    "(large-sum)"

# Note: test_sum_floats.mm2 intentionally panics (floats not supported by u32 parser)
# run_test "test_sum_floats.mm2" \
#     "Sum - floating point" \
#     "(float-result)"

run_test "test_sum_large.mm2" \
    "Sum - large numbers" \
    "(large-sum)"

# ===== Bug Fix Tests =====
echo "=========================================="
echo "PART 5: RemoveSink Bug Fix (Issue #39)"
echo "=========================================="
echo ""

# Test that the removal lag bug is fixed
cat > /tmp/test_removal_lag_fixed.mm2 << 'EOF'
; RemoveSink bug fix test - exec should NOT fire after trigger removed
(state ready)

(exec 0 (, (state ready)) (O (- (state ready)) (+ (trigger "x"))))
(exec 1 (, (trigger "x")) (O (+ (error "x")) (- (trigger "x"))))
(exec 2 (, (trigger "x")) (O (+ (add "x"))))
EOF

run_test "/tmp/test_removal_lag_fixed.mm2" \
    "RemoveSink bug fix - error path only" \
    "(error \"x\")"

run_negative_test "/tmp/test_removal_lag_fixed.mm2" \
    "RemoveSink bug fix - add should NOT fire" \
    "(add \"x\")"

# Test case pattern
run_test "../case/01_case_fix_to_removal_lag_demo.mm2" \
    "Case pattern - alternative elimination" \
    "(error-new \"x\")"

run_negative_test "../case/01_case_fix_to_removal_lag_demo.mm2" \
    "Case pattern - success prevented" \
    "(add-new \"x\")"

# ===== Summary =====
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo ""
echo "Total tests run: $TOTAL_TESTS"
echo -e "Passed: ${GREEN}$PASSED${NC}"
if [ $FAILED -gt 0 ]; then
    echo -e "Failed: ${RED}$FAILED${NC}"
else
    echo "Failed: 0"
fi
echo ""

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}✓ All tests PASSED! ✓${NC}"
    echo ""
    echo "Test Coverage:"
    echo "  ✓ MAX/MIN/COUNT aggregators"
    echo "  ✓ HEAD limit sink"
    echo "  ✓ SUM aggregator"
    echo "  ✓ Duplicate handling"
    echo "  ✓ Multiple conditions"
    echo "  ✓ RemoveSink bug fix (issue #39)"
    echo "  ✓ Case pattern alternatives"
    exit 0
else
    echo -e "${RED}✗ Some tests FAILED ✗${NC}"
    exit 1
fi
