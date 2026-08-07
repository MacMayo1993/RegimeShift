"""Provenance capture (``regimeshift.manifest``).

A manifest recording only ``commit`` plus ``dirty: true`` cannot identify the
source state a run came from -- any account of *how* the tree differed lives in
prose beside the results, not in the record. These tests pin the two mechanisms
that close that gap: a dirty tree is described rather than merely flagged, and a
release run can refuse to be certified from one at all.
"""

from __future__ import annotations

import json

import pytest

from regimeshift import manifest as manifest_mod
from regimeshift.manifest import write_manifest


SPEC = {"groups": (3,), "scenarios": ("shared_orbit",), "effects": (0.5,),
        "segment_lengths": (100,), "n_alt": 2, "n_null": 2}


def _write(tmp_path, **kwargs):
    return write_manifest(
        tmp_path, grid="quick", spec=SPEC, base_seed=1,
        n_configs=1, n_datasets=4, workers=1, **kwargs
    )


def test_manifest_records_what_made_the_tree_dirty(tmp_path, monkeypatch):
    monkeypatch.setattr(manifest_mod, "git_commit", lambda: {
        "commit": "abc123",
        "dirty": True,
        "tracked_diff_sha256": "deadbeef",
        "status_lines": [" M regimeshift/detectors.py", "?? notes.txt"],
        "untracked": ["notes.txt"],
    })
    path = _write(tmp_path)
    git = json.loads(path.read_text())["git"]
    assert git["dirty"] is True
    assert git["tracked_diff_sha256"] == "deadbeef"
    assert git["untracked"] == ["notes.txt"]


def test_require_clean_refuses_a_dirty_tree(tmp_path, monkeypatch):
    monkeypatch.setattr(manifest_mod, "git_commit", lambda: {
        "commit": "abc123", "dirty": True,
        "status_lines": [" M regimeshift/detectors.py"], "untracked": [],
    })
    with pytest.raises(RuntimeError, match="dirty working tree"):
        _write(tmp_path, require_clean=True)
    assert not (tmp_path / "run_manifest.json").exists(), (
        "a refused release run must leave no manifest behind"
    )


def test_require_clean_accepts_a_clean_tree(tmp_path, monkeypatch):
    monkeypatch.setattr(manifest_mod, "git_commit",
                        lambda: {"commit": "abc123", "dirty": False})
    path = _write(tmp_path, require_clean=True)
    assert json.loads(path.read_text())["git"]["commit"] == "abc123"


def test_git_commit_reports_this_checkout():
    """The real helper runs against the repository the tests live in."""
    record = manifest_mod.git_commit()
    assert set(record) >= {"commit"}
    if record.get("dirty"):
        assert "tracked_diff_sha256" in record
        assert record["status_lines"]


def test_write_manifest_accepts_a_captured_git_state(tmp_path, monkeypatch):
    """The CLI captures provenance before the run; the manifest must use it.

    Recording the state at the end would describe a tree the run did not start
    from -- and refusing a release run only after it finished would be useless.
    """
    monkeypatch.setattr(manifest_mod, "git_commit", lambda: {
        "commit": "end-state", "dirty": True, "status_lines": ["?? results/x.csv"],
    })
    captured = {"commit": "start-state", "dirty": False}
    path = _write(tmp_path, require_clean=True, git=captured)
    assert json.loads(path.read_text())["git"]["commit"] == "start-state"
