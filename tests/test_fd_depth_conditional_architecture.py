"""Synthetic-only tests of archived FD algebra; no data or protected execution."""
from pathlib import Path
import runpy

def test_depth_conditional_exact_identities():
    p=Path(__file__).resolve().parents[1]/"shared_db_research/scripts/fd_depth_conditional_architecture.py"
    runpy.run_path(str(p))["self_test"]()
