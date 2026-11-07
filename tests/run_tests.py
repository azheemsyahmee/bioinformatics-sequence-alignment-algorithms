"""
Comprehensive Test Runner for Sequence Alignment Algorithms.
"""

import unittest
import sys
import os

# Ensure package root is in path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from tests.test_needleman_wunsch import TestNeedlemanWunsch
from tests.test_smith_waterman import TestSmithWaterman


def run_all_tests():
    """Discover and execute all test cases in the test suite."""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    suite.addTests(loader.loadTestsFromTestCase(TestNeedlemanWunsch))
    suite.addTests(loader.loadTestsFromTestCase(TestSmithWaterman))

    runner = unittest.TextTestRunner(verbosity=2)
    print("=" * 70)
    print("RUNNING BIOINFORMATICS SEQUENCE ALIGNMENT TEST SUITE")
    print("=" * 70)
    result = runner.run(suite)
    print("=" * 70)

    if result.wasSuccessful():
        print(f"SUCCESS: All {result.testsRun} tests passed successfully!")
        return 0
    else:
        print(f"FAILED: {len(result.failures)} failures, {len(result.errors)} errors.")
        return 1


if __name__ == "__main__":
    sys.exit(run_all_tests())
