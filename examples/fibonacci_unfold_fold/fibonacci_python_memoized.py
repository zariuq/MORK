#!/usr/bin/env python3

from functools import lru_cache


@lru_cache(maxsize=None)
def fib(n: int) -> int:
    if n < 0:
        raise ValueError("n must be non-negative")
    if n < 2:
        return n
    return fib(n - 1) + fib(n - 2)


def main():
    # Compute in ascending order so later calls hit cache.
    cases = [("test-8", 8), ("test-10", 10), ("test-12", 12), ("test-19", 19)]
    for label, n in cases:
        print(f"(fib-result {label} {fib(n)})")


if __name__ == "__main__":
    main()

