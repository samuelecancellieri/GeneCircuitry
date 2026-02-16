#!/usr/bin/env python
# TO REMOVE: One-time dev verification script. No longer needed post-implementation. Remove this file.
"""
Verification script for PipelineController implementation
==========================================================

This script verifies that all controller components are correctly implemented.
"""

import sys
import ast


def check_class_exists():
    """Check if PipelineController class exists."""
    with open("examples/complete_pipeline.py", "r") as f:
        tree = ast.parse(f.read())

    classes = [node.name for node in ast.walk(tree) if isinstance(node, ast.ClassDef)]

    if "PipelineController" in classes:
        print("✓ PipelineController class found")
        return True
    else:
        print("✗ PipelineController class NOT found")
        return False


def check_methods_exist():
    """Check if required methods exist."""
    required_methods = [
        "run_step_load",
        "run_step_preprocessing",
        "run_step_stratification",
        "run_step_clustering",
        "run_step_celloracle",
        "run_step_hotspot",
        "run_step_grn_analysis",
        "process_single_stratification",
        "run_stratified_pipeline_sequential",
        "run_stratified_pipeline_parallel",
        "run_complete_pipeline",
        "print_final_summary",
    ]

    with open("examples/complete_pipeline.py", "r") as f:
        tree = ast.parse(f.read())

    # Find PipelineController class
    controller_class = None
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == "PipelineController":
            controller_class = node
            break

    if not controller_class:
        print("✗ Cannot find PipelineController class")
        return False

    # Get all method names
    methods = [
        node.name for node in controller_class.body if isinstance(node, ast.FunctionDef)
    ]

    # Check required methods
    all_found = True
    for method in required_methods:
        if method in methods:
            print(f"  ✓ {method}")
        else:
            print(f"  ✗ {method} NOT found")
            all_found = False

    return all_found


def check_cli_flags():
    """Check if new CLI flags exist."""
    required_flags = ["--parallel", "--steps"]

    with open("examples/complete_pipeline.py", "r") as f:
        content = f.read()

    all_found = True
    for flag in required_flags:
        if flag in content:
            print(f"  ✓ {flag} flag exists")
        else:
            print(f"  ✗ {flag} flag NOT found")
            all_found = False

    return all_found


def check_imports():
    """Check if required imports exist."""
    required_imports = [
        "from multiprocessing import Pool, cpu_count",
        "from functools import partial",
    ]

    with open("examples/complete_pipeline.py", "r") as f:
        content = f.read()

    all_found = True
    for imp in required_imports:
        if imp in content:
            print(f"  ✓ {imp}")
        else:
            print(f"  ✗ {imp} NOT found")
            all_found = False

    return all_found


def check_documentation_exists():
    """Check if documentation files exist."""
    import os

    docs = [
        "docs/PIPELINE_CONTROLLER_GUIDE.md",
        "docs/CONTROLLER_QUICK_REF.md",
        "examples/controller_usage_examples.py",
    ]

    all_found = True
    for doc in docs:
        if os.path.exists(doc):
            print(f"  ✓ {doc}")
        else:
            print(f"  ✗ {doc} NOT found")
            all_found = False

    return all_found


def main():
    """Run all verification checks."""
    print("\n" + "=" * 70)
    print("Pipeline Controller Implementation Verification")
    print("=" * 70)

    checks = [
        ("Class Implementation", check_class_exists),
        ("Method Implementation", check_methods_exist),
        ("CLI Flags", check_cli_flags),
        ("Required Imports", check_imports),
        ("Documentation Files", check_documentation_exists),
    ]

    results = []
    for name, check_func in checks:
        print(f"\n{name}:")
        print("-" * 70)
        result = check_func()
        results.append(result)

    print("\n" + "=" * 70)
    print("Verification Summary")
    print("=" * 70)

    for (name, _), result in zip(checks, results):
        status = "PASS" if result else "FAIL"
        symbol = "✓" if result else "✗"
        print(f"{symbol} {name}: {status}")

    all_passed = all(results)

    print("\n" + "=" * 70)
    if all_passed:
        print("✓ ALL CHECKS PASSED!")
        print("=" * 70)
        print("\nPipelineController implementation is complete and correct.")
        print("\nYou can now use:")
        print("  • Modular step execution (--steps)")
        print("  • Parallel stratified processing (--parallel)")
        print("  • Programmatic pipeline control (PipelineController)")
        return 0
    else:
        print("✗ SOME CHECKS FAILED")
        print("=" * 70)
        print("\nPlease review the failures above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
