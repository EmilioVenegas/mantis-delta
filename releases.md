# Releasing mantis-delta

How a new version of `mantis-delta` gets to PyPI.

## How it works (read once)

- **Versioning is automatic** — [setuptools-scm](https://setuptools-scm.readthedocs.io)
  derives the version from the latest **git tag**. There is no `version =` field in
  `pyproject.toml` and no `__version__` string to bump. The tag *is* the version.
- **Publishing is automatic** — pushing a **GitHub Release** triggers
  `.github/workflows/publish.yml`, which builds the wheel + sdist and uploads them to
  PyPI via **trusted publishing (OIDC)**. No API tokens, no `twine upload` by hand.
- **CI builds are always clean** — the workflow checks out the exact tagged commit in a
  fresh clone, so released versions never carry a `dev`/`dirty` suffix (see below).
- **The changelog is generated** — [git-cliff](https://git-cliff.org) builds `CHANGELOG.md`
  from the Conventional-Commits history (config in `cliff.toml`). The same tool produces
  the GitHub Release notes, so `CHANGELOG.md` is the single source of truth.
- **The version number is derived, not typed** — `scripts/release.sh` asks git-cliff to
  compute the next version from your commit types since the last tag (`feat:` → MINOR,
  `fix:` → PATCH, per the `[bump]` policy in `cliff.toml`). This is the guardrail against
  picking a number by hand and accidentally shipping features as a patch release.
  **Never pass `git-cliff --tag vX.Y.Z` yourself.**

So a release is just: **run `scripts/release.sh` → push → cut a GitHub Release.**

Build tooling below uses [`uv`](https://docs.astral.sh/uv/) (`uv build`, `uvx`), so there
is nothing to `pip install` — `uvx` runs `twine` and `git-cliff` in throwaway environments.

## Release checklist (in order)

```fish
cd ~/Documents/hairpin/mantis

# 1. Start from an up-to-date main
git checkout main
git pull

# 2. Make sure the working tree is clean and the tests pass
git status                      # must say "nothing to commit, working tree clean"
pytest -q                       # run the suite in your dev env

# 3. Cut the release. ONE command does it all: it computes the next version from
#    your commits, shows you what will be released, then (after you confirm)
#    regenerates CHANGELOG.md, commits "chore(release): vX.Y.Z", and tags it.
#    It does NOT push or publish — those stay manual (steps 5–6).
scripts/release.sh

# 4. (optional) Validate the build locally. The tag now exists, so setuptools-scm
#    stamps the clean version. PyPI is IMMUTABLE — catch errors here, not after upload.
uv build                        # writes wheel + sdist to dist/
uvx twine check dist/*          # metadata / README render check
rm -rf dist                     # clean up the throwaway artifacts

# 5. Push the release commit AND its tag together.
git push origin main --follow-tags
set TAG (git describe --tags)   # the version the script just created, e.g. v0.3.0

# 6. Cut the GitHub Release using the changelog section as the body.
#    THIS is the event that publishes to PyPI.
uvx git-cliff --latest --strip header -o /tmp/notes.md
gh release create $TAG --title "$TAG" --notes-file /tmp/notes.md
#    ...or in the web UI: Releases -> Draft a new release -> pick the tag ->
#    paste that version's section from CHANGELOG.md -> Publish.

# 7. Watch the publish job
gh run watch                    # or the Actions tab on GitHub
```

After the workflow goes green, `pip install mantis-delta==<that version>` works.

> **Doing it by hand instead?** If you ever bypass `release.sh`, the version must still
> come from git-cliff — `set TAG (uvx git-cliff --bumped-version)` — never a number you
> guessed. Then `uvx git-cliff --bump -o CHANGELOG.md`, commit, and `git tag $TAG`.

## How the version number is chosen

You don't pick it — `scripts/release.sh` derives it from your commits via the `[bump]`
policy in `cliff.toml`, following [semantic versioning](https://semver.org)
`vMAJOR.MINOR.PATCH`:

- **PATCH** (`v0.2.0 -> v0.2.1`) — only `fix:`/`perf:` etc. since the last tag.
- **MINOR** (`v0.2.0 -> v0.3.0`) — any `feat:` present (`features_always_bump_minor`).
- **MAJOR** (`-> v1.0.0`) — a breaking change. **Held to MINOR while pre-1.0**
  (`breaking_always_bump_major = false`); flip that in `cliff.toml` when you're ready for 1.0.

## Commit messages drive the changelog

git-cliff groups entries by the [Conventional Commits](https://www.conventionalcommits.org)
prefix on each commit subject. Keep using them:

| Prefix      | Shows up under | Example                                  |
| ----------- | -------------- | ---------------------------------------- |
| `feat:`     | Features       | `feat: add τ-leap adaptive step control` |
| `fix:`      | Bug Fixes      | `fix: correct propensity recalculation`  |
| `perf:`     | Performance    | `perf: vectorize the Gillespie inner loop` |
| `refactor:` | Refactor       | `refactor: split the ODE/SSA backends`   |
| `docs:`     | Documentation  | `docs: refresh README examples`          |
| `test:`     | Testing        | `test: cover the τ-leap sampler`         |

`chore:`, `build:`, `ci:`, `style:` are intentionally **excluded** from the changelog
(housekeeping noise). Commits with **no recognized prefix are dropped entirely** — that is
the "N commits skipped due to parse errors" warning git-cliff prints. Usually harmless, but
mind it: a real `feat`/`fix` with a forgotten prefix vanishes silently from the changelog.

**Write a commit body, not just a subject.** `cliff.toml` renders the commit *body*
(indented, under the bullet) and *scope* (bold prefix) in the changelog and release notes.
A subject-only commit produces a one-line, undetailed entry. To make this easy, the
`prepare-commit-msg` git hook drafts a Conventional-Commits message (subject + body) from
your staged diff and pre-fills the editor when you run `git commit` — review and edit it,
then save. (It's in `.git/hooks/`, so it's local to your clone and not version-controlled.)

Regenerate the full file from history at any time with `uvx git-cliff -o CHANGELOG.md`.
Adjust grouping/exclusions in `cliff.toml`.

## Understanding the version string

setuptools-scm produces a clean `X.Y.Z` **only** when HEAD is exactly on a tag **and**
the working tree is clean. Otherwise it adds a suffix that tells you why this build is not
a real release:

| Version you see                 | Meaning                                              |
| ------------------------------- | ---------------------------------------------------- |
| `0.2.0`                         | clean release — HEAD is on tag `v0.2.0`              |
| `0.2.1.dev4+g211454a`           | 4 commits **past** tag `v0.2.0` (untagged work)      |
| `0.2.1.dev4+g211454a.d20260613` | ...and the working tree is **dirty** (uncommitted)   |

These suffixes only appear in local builds — that is exactly their purpose: marking an
artifact as "not an official release." CI builds from a clean tagged checkout, so it never
sees them.

## Gotchas

- **PyPI is immutable.** Once `0.2.0` is uploaded it can never be overwritten or reused.
  A mistake means burning the number and releasing the next one. Hence the local `uv build`
  + `uvx twine check` in step 4.
- **Only a published Release triggers the upload.** Pushing a tag alone does *not* publish
  — the `release: published` event does.
- **Never delete + recreate a tag to "redo" a release** that already hit PyPI; the version
  is gone. Bump to the next number instead.
