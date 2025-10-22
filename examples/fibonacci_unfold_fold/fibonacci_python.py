#!/usr/bin/env python3

def fib_fast_doubling(n: int) -> int:
    """Compute F(n) using fast doubling in O(log n)."""
    def _pair(k: int):
        if k == 0:
            return (0, 1)
        a, b = _pair(k >> 1)
        c = a * (2 * b - a)
        d = a * a + b * b
        if k & 1:
            return (d, c + d)
        else:
            return (c, d)

    return _pair(n)[0]


def main():
    # Match the format seen in fibonacci_table.mm2 output
    cases = [("test-8", 8), ("test-10", 10), ("test-12", 12), ("test-19", 19), ("test-1000", 1000)]
    for label, n in cases:
        value = fib_fast_doubling(n)
        print(f"(fib-result {label} {value})")


if __name__ == "__main__":
    main()

