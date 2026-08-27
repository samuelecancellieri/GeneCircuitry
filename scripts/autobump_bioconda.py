#!/usr/bin/env python3
"""
Autobump Bioconda Recipe Script for GeneCircuitry

Updates meta.yaml for genecircuitry with the new version and sha256 checksum from PyPI,
and optionally forks bioconda/bioconda-recipes, commits the change, and opens a Pull Request.

Usage:
    python scripts/autobump_bioconda.py --version 0.2.2 --dry-run
    python scripts/autobump_bioconda.py --version 0.2.2 --token "$BIOCONDA_TOKEN"
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional, Tuple


def get_pypi_release_info(
    package_name: str, version: str, max_retries: int = 12, retry_delay: int = 30
) -> Tuple[str, str]:
    """
    Poll PyPI for release metadata and compute the SHA256 checksum of the sdist.

    Returns:
        Tuple of (tarball_url, sha256_checksum)
    """
    url = f"https://pypi.org/pypi/{package_name}/{version}/json"
    print(f"Polling PyPI metadata at {url} ...")

    data = None
    for attempt in range(1, max_retries + 1):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "GeneCircuitry-Bioconda-Autobump/1.0"}
            )
            with urllib.request.urlopen(req, timeout=15) as resp:
                if resp.status == 200:
                    data = json.loads(resp.read().decode("utf-8"))
                    print(f"PyPI package {package_name} v{version} is available.")
                    break
        except urllib.error.HTTPError as e:
            if e.code == 404:
                print(
                    f"Attempt {attempt}/{max_retries}: Version {version} not yet on PyPI (HTTP 404). Retrying in {retry_delay}s..."
                )
            else:
                print(
                    f"Attempt {attempt}/{max_retries}: HTTP error {e.code}. Retrying in {retry_delay}s..."
                )
        except Exception as e:
            print(
                f"Attempt {attempt}/{max_retries}: Network error ({e}). Retrying in {retry_delay}s..."
            )

        if attempt < max_retries:
            time.sleep(retry_delay)

    if not data:
        raise RuntimeError(
            f"Package {package_name} v{version} did not appear on PyPI within {max_retries * retry_delay}s."
        )

    # Find sdist (.tar.gz)
    sdist_file = None
    for item in data.get("urls", []):
        if item.get("packagetype") == "sdist" or item.get("filename", "").endswith(
            ".tar.gz"
        ):
            sdist_file = item
            break

    if sdist_file and "digests" in sdist_file and "sha256" in sdist_file["digests"]:
        sdist_url = sdist_file["url"]
        sha256_hash = sdist_file["digests"]["sha256"]
        print(f"Found sdist from PyPI API: {sdist_url}")
        print(f"SHA256 from PyPI API: {sha256_hash}")
        return sdist_url, sha256_hash

    # Fallback to direct download & hash computation
    sdist_url = f"https://pypi.io/packages/source/{package_name[0]}/{package_name}/{package_name}-{version}.tar.gz"
    print(f"Downloading sdist to compute SHA256: {sdist_url}")
    req = urllib.request.Request(
        sdist_url, headers={"User-Agent": "GeneCircuitry-Bioconda-Autobump/1.0"}
    )
    with urllib.request.urlopen(req, timeout=30) as resp:
        content = resp.read()
        sha256_hash = hashlib.sha256(content).hexdigest()
    print(f"Computed SHA256: {sha256_hash}")
    return sdist_url, sha256_hash


def update_meta_yaml_content(content: str, version: str, sha256: str) -> str:
    """
    Update version, sha256, and reset build number in meta.yaml content string.
    """
    # 1. Update version: {% set version = "..." %}
    content = re.sub(
        r'^(\s*\{%\s*set\s+version\s*=\s*")[^"]*("\s*%\})',
        r"\g<1>" + version + r"\g<2>",
        content,
        flags=re.MULTILINE,
    )

    # 2. Update sha256: handles both commented and uncommented forms
    content = re.sub(
        r"^(\s*)(?:#\s*)?sha256:\s*.*$",
        r"\g<1>sha256: " + sha256,
        content,
        flags=re.MULTILINE,
    )

    # 3. Reset build number: number: 0
    content = re.sub(
        r"^(\s*number:\s*)\d+",
        r"\g<1>0",
        content,
        flags=re.MULTILINE,
    )

    return content


def update_local_recipe(
    recipe_path: Path, version: str, sha256: str, dry_run: bool = False
) -> bool:
    """
    Updates the local recipe file. Returns True if file was changed.
    """
    if not recipe_path.exists():
        raise FileNotFoundError(f"Recipe file not found: {recipe_path}")

    original = recipe_path.read_text(encoding="utf-8")
    updated = update_meta_yaml_content(original, version, sha256)

    if original == updated:
        print(f"Local recipe at {recipe_path} is already up to date.")
        return False

    if dry_run:
        print(f"[DRY-RUN] Would update {recipe_path} to version {version}, sha256 {sha256}")
    else:
        recipe_path.write_text(updated, encoding="utf-8")
        print(f"Updated {recipe_path} to version {version}, sha256 {sha256}")

    return True


def run_cmd(
    cmd: list[str],
    cwd: Optional[Path] = None,
    env: Optional[dict] = None,
    check: bool = True,
    capture_output: bool = True,
) -> subprocess.CompletedProcess:
    """Helper to run a shell command and print outputs."""
    print(f"+ {' '.join(cmd)}" + (f" (in {cwd})" if cwd else ""))
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            env=env,
            check=check,
            capture_output=capture_output,
            text=True,
        )
        return result
    except subprocess.CalledProcessError as exc:
        if exc.stdout:
            print(f"[Command stdout]\n{exc.stdout.strip()}")
        if exc.stderr:
            print(f"[Command stderr]\n{exc.stderr.strip()}", file=sys.stderr)
        raise


def submit_bioconda_pr(
    package_name: str,
    version: str,
    sha256: str,
    meta_yaml_source_path: Path,
    bioconda_repo: str = "bioconda/bioconda-recipes",
    token: Optional[str] = None,
    dry_run: bool = False,
) -> Optional[str]:
    """
    Forks and clones bioconda/bioconda-recipes, updates recipes/<package>/meta.yaml,
    commits, pushes to branch, and creates a Pull Request.

    Returns:
        The PR URL if created, or None.
    """
    # Only use explicit BIOCONDA_TOKEN or GH_TOKEN PAT. Do not fall back to GITHUB_TOKEN,
    # as the repository-scoped GITHUB_TOKEN cannot fork or push to user forks.
    token = token or os.environ.get("BIOCONDA_TOKEN") or os.environ.get("GH_TOKEN")
    if not token and not dry_run:
        print(
            "::warning::No BIOCONDA_TOKEN or GH_TOKEN secret found.\n"
            "::warning::Cross-repository Bioconda PR submission requires a Personal Access Token (PAT) "
            "with 'public_repo' (or 'repo') scope.\n"
            "::warning::The default GITHUB_TOKEN is scoped to this repository and cannot push to user forks.\n"
            "::warning::To enable automated Bioconda PRs:\n"
            "::warning::  1. Create a GitHub Personal Access Token (classic) with 'public_repo' scope.\n"
            "::warning::  2. Add it as repository secret 'BIOCONDA_TOKEN' under Settings > Secrets and variables > Actions.\n"
            "::warning::Skipping Bioconda PR submission."
        )
        return None

    branch_name = f"bump-{package_name}-v{version}"
    pr_title = f"{package_name} v{version}"
    pr_body = (
        f"Automated recipe update for `{package_name}` version **{version}** from PyPI.\n\n"
        f"- **Package**: `{package_name}`\n"
        f"- **Version**: `{version}`\n"
        f"- **SHA256**: `{sha256}`\n"
        f"- **Source**: https://pypi.org/project/{package_name}/{version}/\n\n"
        f"---\n"
        f"*Triggered by release of {package_name} v{version} via GitHub Actions.*"
    )

    if dry_run:
        print(f"[DRY-RUN] Target upstream repo: {bioconda_repo}")
        print(f"[DRY-RUN] Branch: {branch_name}")
        print(f"[DRY-RUN] PR Title: {pr_title}")
        print(f"[DRY-RUN] PR Body:\n{pr_body}")
        return None

    # Check for GitHub CLI
    if not shutil.which("gh"):
        raise RuntimeError("GitHub CLI (`gh`) is required to fork, push, and open PRs.")

    env = os.environ.copy()
    if token:
        env["GH_TOKEN"] = token
        env["GITHUB_TOKEN"] = token

    # Check authentication
    auth_check = run_cmd(["gh", "auth", "status"], env=env, check=False)
    if auth_check.returncode != 0 and not token:
        print("::warning::GitHub CLI is not authenticated. Please set BIOCONDA_TOKEN.")
        return None

    # Get authenticated username
    user_res = run_cmd(["gh", "api", "user", "-q", ".login"], env=env, check=False)
    gh_user = user_res.stdout.strip() if user_res.returncode == 0 else ""
    if not gh_user:
        gh_user = os.environ.get("GITHUB_REPOSITORY_OWNER") or "samuelecancellieri"
    print(f"Authenticated as GitHub user: {gh_user}")

    temp_dir = Path(tempfile.mkdtemp(prefix="bioconda_bump_"))
    try:
        print(f"Working in temporary directory: {temp_dir}")
        print(f"Cloning {bioconda_repo} (master branch)...")
        run_cmd(
            ["git", "clone", "--depth", "10", "--branch", "master", f"https://github.com/{bioconda_repo}.git", "bioconda-recipes"],
            cwd=temp_dir,
            env=env,
        )
        work_tree = temp_dir / "bioconda-recipes"

        # Configure git user
        run_cmd(["git", "config", "user.name", "github-actions[bot]"], cwd=work_tree, env=env)
        run_cmd(["git", "config", "user.email", "github-actions[bot]@users.noreply.github.com"], cwd=work_tree, env=env)

        # Create new branch
        run_cmd(["git", "checkout", "-b", branch_name], cwd=work_tree, env=env)

        # Locate recipe
        recipe_dir = work_tree / "recipes" / package_name
        recipe_file = recipe_dir / "meta.yaml"
        recipe_dir.mkdir(parents=True, exist_ok=True)

        if recipe_file.exists():
            original_content = recipe_file.read_text(encoding="utf-8")
            updated_content = update_meta_yaml_content(original_content, version, sha256)
        else:
            # If package recipe does not exist in bioconda yet, copy our local recipe
            print(f"Recipe not found in bioconda-recipes for {package_name}. Copying from source...")
            original_content = meta_yaml_source_path.read_text(encoding="utf-8")
            updated_content = update_meta_yaml_content(original_content, version, sha256)

        recipe_file.write_text(updated_content, encoding="utf-8")

        # Check diff
        diff_res = run_cmd(["git", "diff", "--stat"], cwd=work_tree, env=env)
        if not diff_res.stdout.strip():
            print(f"No changes detected in recipes/{package_name}/meta.yaml.")
            return None

        # Commit
        commit_msg = f"{package_name} v{version}"
        run_cmd(["git", "add", f"recipes/{package_name}/meta.yaml"], cwd=work_tree, env=env)
        run_cmd(["git", "commit", "-m", commit_msg], cwd=work_tree, env=env)

        # Set up push remote: push to user fork
        print(f"Ensuring fork exists for user {gh_user}...")
        fork_res = run_cmd(["gh", "repo", "fork", bioconda_repo, "--clone=false"], env=env, check=False)
        if fork_res.returncode != 0:
            err_msg = (fork_res.stderr or "").strip()
            if err_msg and "already exists" not in err_msg.lower():
                print(f"gh repo fork note: {err_msg}")

        fork_remote_url = f"https://x-access-token:{token}@github.com/{gh_user}/bioconda-recipes.git"
        run_cmd(["git", "remote", "add", "user_fork", fork_remote_url], cwd=work_tree, env=env)

        print(f"Pushing branch {branch_name} to fork {gh_user}/bioconda-recipes...")
        try:
            run_cmd(["git", "push", "-u", "user_fork", branch_name, "--force"], cwd=work_tree, env=env)
        except subprocess.CalledProcessError:
            print(
                f"::error::Failed to push branch '{branch_name}' to fork '{gh_user}/bioconda-recipes'.\n"
                f"::error::Please check that:\n"
                f"::error::  1. The secret BIOCONDA_TOKEN is a valid GitHub Personal Access Token (classic) with 'public_repo' (or 'repo') scope.\n"
                f"::error::  2. The fork repository https://github.com/{gh_user}/bioconda-recipes exists and the token has write access to it.\n"
                f"::error::  3. The token has not expired.",
                file=sys.stderr,
            )
            raise

        # Check if PR already exists
        print("Checking if PR already exists...")
        pr_list = run_cmd(
            ["gh", "pr", "list", "--repo", bioconda_repo, "--head", f"{gh_user}:{branch_name}", "--json", "url", "-q", ".[0].url"],
            env=env,
            check=False,
        )
        existing_pr = pr_list.stdout.strip()
        if existing_pr:
            print(f"PR already exists: {existing_pr}")
            return existing_pr

        # Create PR
        print(f"Opening Pull Request to {bioconda_repo}...")
        pr_create = run_cmd(
            [
                "gh", "pr", "create",
                "--repo", bioconda_repo,
                "--base", "master",
                "--head", f"{gh_user}:{branch_name}",
                "--title", pr_title,
                "--body", pr_body,
            ],
            cwd=work_tree,
            env=env,
        )
        pr_url = pr_create.stdout.strip()
        print(f"Successfully created Bioconda PR: {pr_url}")
        return pr_url

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def main() -> int:
    parser = argparse.ArgumentParser(description="Autobump Bioconda recipe for GeneCircuitry")
    parser.add_argument("--version", type=str, help="Target release version (e.g. 0.2.3)")
    parser.add_argument("--sha256", type=str, default=None, help="Pre-computed SHA256 hash (optional)")
    parser.add_argument("--package", type=str, default="genecircuitry", help="Package name (default: genecircuitry)")
    parser.add_argument("--recipe-path", type=str, default="conda-recipe/meta.yaml", help="Path to local meta.yaml")
    parser.add_argument("--bioconda-repo", type=str, default="bioconda/bioconda-recipes", help="Upstream Bioconda repo")
    parser.add_argument("--token", type=str, default=None, help="GitHub token with repo permissions")
    parser.add_argument("--local-only", action="store_true", help="Only update local recipe, do not open Bioconda PR")
    parser.add_argument("--dry-run", action="store_true", help="Simulate run without writing files or making network mutations")

    args = parser.parse_args()

    # Determine version if not given
    version = args.version
    if not version:
        tag = os.environ.get("GITHUB_REF_NAME", "")
        if tag.startswith("v"):
            version = tag[1:]
        elif tag:
            version = tag

    if not version:
        # Try finding version from local package
        try:
            import genecircuitry
            version = getattr(genecircuitry, "__version__", None)
        except Exception:
            pass

    if not version:
        print("Error: Version could not be determined. Please specify --version.")
        return 1

    # Clean version string (remove any leading 'v')
    version = version.lstrip("v")
    print(f"Starting autobump for {args.package} v{version}...")

    # Compute or fetch SHA256 from PyPI
    sha256 = args.sha256
    if not sha256:
        try:
            _, sha256 = get_pypi_release_info(args.package, version)
        except Exception as e:
            print(f"Error obtaining PyPI release info: {e}")
            return 1

    # 1. Update local conda-recipe/meta.yaml
    local_recipe = Path(args.recipe_path)
    if local_recipe.exists():
        update_local_recipe(local_recipe, version, sha256, dry_run=args.dry_run)
    else:
        print(f"Warning: Local recipe {local_recipe} not found.")

    # 2. Submit PR to bioconda/bioconda-recipes
    if not args.local_only:
        try:
            pr_url = submit_bioconda_pr(
                package_name=args.package,
                version=version,
                sha256=sha256,
                meta_yaml_source_path=local_recipe,
                bioconda_repo=args.bioconda_repo,
                token=args.token,
                dry_run=args.dry_run,
            )
            if pr_url:
                print(f"Bioconda Pull Request: {pr_url}")
        except Exception as e:
            print(f"Error submitting Bioconda PR: {e}")
            return 1

    print(f"Autobump completed successfully for {args.package} v{version}!")
    return 0


if __name__ == "__main__":
    sys.exit(main())

