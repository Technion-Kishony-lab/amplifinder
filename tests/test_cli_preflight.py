import subprocess
from pathlib import Path

import pytest

from amplifinder.config import Config, ISDetectionMethod
from amplifinder.exceptions import ToolNotFoundError


def _minimal_config(tmp_path: Path, *, is_method: ISDetectionMethod = ISDetectionMethod.GENBANK) -> Config:
    iso_dir = tmp_path / "iso"
    iso_dir.mkdir()
    return Config(
        iso_fastq_path=iso_dir,
        ref_name="U00096",
        is_detection_method=is_method,
    )


def test_preflight_missing_bowtie2(monkeypatch, tmp_path):
    import amplifinder.cli as cli

    cfg = _minimal_config(tmp_path)

    def fake_get_tool_path(name, config_path=None, required=True):
        if name == "bowtie2":
            raise ToolNotFoundError("bowtie2 not found", include_help=False)
        return Path(f"/fake/{name}")

    monkeypatch.setattr(cli, "BRESEQ_DOCKER", False)
    monkeypatch.setattr(cli, "get_tool_path", fake_get_tool_path)

    with pytest.raises(ToolNotFoundError) as exc:
        cli._preflight_dependencies(cfg)

    assert "bowtie2 not found" in str(exc.value)


def test_preflight_docker_missing(monkeypatch, tmp_path):
    import amplifinder.cli as cli

    cfg = _minimal_config(tmp_path)

    monkeypatch.setattr(cli, "BRESEQ_DOCKER", True)
    monkeypatch.setattr(cli.shutil, "which", lambda _: None)

    # Tools other than docker are present
    monkeypatch.setattr(cli, "get_tool_path", lambda name, config_path=None, required=True: Path(f"/fake/{name}"))

    with pytest.raises(ToolNotFoundError) as exc:
        cli._preflight_dependencies(cfg)

    assert "docker not found" in str(exc.value).lower()


def test_preflight_docker_daemon_not_running(monkeypatch, tmp_path):
    import amplifinder.cli as cli

    cfg = _minimal_config(tmp_path)

    monkeypatch.setattr(cli, "BRESEQ_DOCKER", True)
    monkeypatch.setattr(cli.shutil, "which", lambda _: "/usr/bin/docker")

    def fake_run(*args, **kwargs):
        raise subprocess.CalledProcessError(returncode=1, cmd=["docker", "info"])

    monkeypatch.setattr(cli.subprocess, "run", fake_run)
    monkeypatch.setattr(cli, "get_tool_path", lambda name, config_path=None, required=True: Path(f"/fake/{name}"))

    with pytest.raises(ToolNotFoundError) as exc:
        cli._preflight_dependencies(cfg)

    msg = str(exc.value).lower()
    assert "docker is installed but not running" in msg
    assert "docker info" in msg


def test_preflight_local_breseq_missing(monkeypatch, tmp_path):
    import amplifinder.cli as cli

    cfg = _minimal_config(tmp_path)

    def fake_get_tool_path(name, config_path=None, required=True):
        if name == "breseq":
            raise ToolNotFoundError("breseq not found", include_help=False)
        return Path(f"/fake/{name}")

    monkeypatch.setattr(cli, "BRESEQ_DOCKER", False)
    monkeypatch.setattr(cli, "get_tool_path", fake_get_tool_path)

    with pytest.raises(ToolNotFoundError) as exc:
        cli._preflight_dependencies(cfg)

    assert "breseq not found" in str(exc.value)


def test_preflight_isfinder_requires_blast(monkeypatch, tmp_path):
    import amplifinder.cli as cli

    cfg = _minimal_config(tmp_path, is_method=ISDetectionMethod.ISFINDER)

    def fake_get_tool_path(name, config_path=None, required=True):
        if name == "blastn":
            raise ToolNotFoundError("blastn not found", include_help=False)
        return Path(f"/fake/{name}")

    monkeypatch.setattr(cli, "BRESEQ_DOCKER", False)
    monkeypatch.setattr(cli, "get_tool_path", fake_get_tool_path)

    with pytest.raises(ToolNotFoundError) as exc:
        cli._preflight_dependencies(cfg)

    assert "blastn not found" in str(exc.value)


def test_create_config_bypasses_preflight(monkeypatch, tmp_path):
    import amplifinder.cli as cli

    iso_dir = tmp_path / "iso"
    iso_dir.mkdir()
    out_dir = tmp_path / "out"

    # If called, this test should fail
    monkeypatch.setattr(cli, "_preflight_dependencies", lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError()))

    cfg_kwargs = {
        "iso_fastq_path": iso_dir,
        "ref_name": "U00096",
        "output_dir": out_dir,
    }
    create_path = tmp_path / "template.yaml"

    cli._run_single(
        config_kwargs=cfg_kwargs,
        create_config_path=create_path,
        breseq_only=False,
        verbose=False,
    )

    assert create_path.exists()

