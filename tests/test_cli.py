"""The CLI surface. These exist because every defect here shipped: implemented
subcommands were refused with a GPU message, and a missing catalogue surfaced as a raw
sqlite3 traceback."""

from __future__ import annotations

import pytest

from dhfrcamp.cli import GPU_GATED, build_parser, main


def test_every_subcommand_is_routed():
    """A subcommand that parses but falls through to 'unknown command' is a dead entry."""
    parser = build_parser()
    actions = [a for a in parser._actions if hasattr(a, "choices") and a.choices]
    names = set(actions[0].choices)
    assert names == {"prepare", "decoys", "screen", "generate", "campaign", "evaluate"}


@pytest.mark.parametrize("command", sorted(GPU_GATED))
def test_gpu_gated_commands_say_why_and_point_somewhere(command):
    with pytest.raises(SystemExit) as excinfo:
        main([command])
    message = str(excinfo.value)
    assert "GPU" in message
    assert "campaign" in message


@pytest.mark.parametrize("command", ["decoys", "evaluate"])
def test_cpu_commands_are_not_refused_as_gpu_work(command):
    """decoys and evaluate are implemented and CPU-only; they used to be refused."""
    assert command not in GPU_GATED


def test_missing_catalogue_names_the_command_that_builds_it(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        main(
            [
                "--results-dir",
                str(tmp_path),
                "campaign",
                "--catalog",
                str(tmp_path / "nope.sqlite"),
            ]
        )
    message = str(excinfo.value)
    assert "catalogue not found" in message
    assert "build_catalog.py" in message


def test_evaluate_without_findings_says_what_to_run(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        main(["--results-dir", str(tmp_path), "evaluate"])
    assert "campaign" in str(excinfo.value)
