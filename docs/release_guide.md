# Release Guide (GitFlow + GitHub Releases → (Test)PyPI)

This repository uses **GitFlow** and **GitHub Releases** as the single source of truth for publishing:
- **Pre-releases** publish automatically to **TestPyPI**
- **Releases** publish automatically to **PyPI**

Publishing is triggered by creating a **GitHub Release** (UI) that points to the correct tag/commit.

---

## Branching standard

### Long-lived branches
- `master`: production-ready history, tagged releases
- `develop`: integration branch for upcoming changes

### Supporting branches
- `feature/<name>`: new features (branch from `develop`, merge back into `develop`)
- `release/<version>`: stabilization & release preparation (branch from `develop`, merge into `master` and `develop`)
- `hotfix/<version>`: urgent fixes (branch from `master`, merge into `master` and `develop`)

---

## Versioning standard

Use **SemVer**: `MAJOR.MINOR.PATCH`
- `PATCH`: bugfix
- `MINOR`: new backwards-compatible features
- `MAJOR`: breaking changes

### Tag naming
- Final release tag: `vX.Y.Z` (example: `v0.3.2`)
- Pre-release tag: `vX.Y.Zrc1` (rc = release candidate)

---

## Required repo state before releasing

1. **All desired changes are on `develop`** (via PRs).
2. CI on `develop` is green. (OPTIONAL: WHEN WE HAVE TESTS WE CAN AUTOMATE)
3. Your release workflows exist on the branch/commit you will release from:
   - If you create a tag pointing to a commit **that predates** the workflow files, the workflow **will not run**.
   - Therefore: ensure the workflows are already merged into `develop`/`master` before tagging and releasing.

---

## Normal release flow 

### 1) Create a release branch
From `develop`:
- Create: `release/X.Y.Z`

On `release/X.Y.Z`, do only release prep:
- Bump version in `pyproject.toml` (and anywhere else the project defines its version)
- Update `CHANGELOG.md` (OPTIONAL: NEED TO THINK WHETHER WE NEED THIS)
- Ensure packaging metadata is correct (README, license, classifiers, etc.)
- Run tests locally (OPTIONAL: WHEN WE HAVE PROPER TESTS)

Push the branch and open a PR:
- `release/X.Y.Z` → `master`

Use **Squash & merge** or **Merge commit** consistently (team choice).
- If you use **Squash & merge**, GitHub will create a new commit on `master`. Tags/releases should point to that commit.

### 2) Merge release branch into `master`
After merging, ensure `master` is green.

### 3) Create the GitHub Release (UI) → triggers PyPI publish
In GitHub:
- Go to **Releases** → **Draft a new release**
- **Tag**: `vX.Y.Z` (create new tag)
- **Target**: `master` (important)
- Title: `vX.Y.Z`
- Description: paste changelog notes
- Ensure **“Set as a pre-release”** is **unchecked**
- Publish release

Result: Your workflow should publish to **PyPI** automatically.

### 4) Merge back into `develop`
Open a PR:
- `master` → `develop`
This keeps `develop` aligned with the released state (including version bumps and changelog).

---

## Pre-release flow (TestPyPI)

Use this when you want to test packaging and installation before a final release.

### Option A (common): pre-release from `develop`
1. Ensure `develop` is green and contains the publishing workflows.
2. Create GitHub Release:
   - Tag: `vX.Y.Z-rc.1` (or `-beta.1`)
   - Target: `develop`
   - Check **“Set as a pre-release”**
   - Publish release

Result: workflow should publish to **TestPyPI**.

### Option B: pre-release from a `release/X.Y.Z` branch
Same as above, but target the `release/X.Y.Z` branch/commit.

**Important:** The workflow must exist on that commit.

---

## Hotfix flow (urgent fix)

1. Branch from `master`:
   - `hotfix/X.Y.(Z+1)`
2. Apply fix + bump version + update changelog
3. PR `hotfix/...` → `master` and merge
4. Create GitHub Release on `master` with tag `vX.Y.(Z+1)` (not pre-release)
   - triggers PyPI publish
5. PR `master` → `develop` to sync

---

## GitHub Release checklist

Before clicking **Publish release**, verify:
- [ ] Tag name is correct (`vX.Y.Z` or `vX.Y.Z-rc.N`)
- [ ] Target branch/commit is correct (`master` for final releases)
- [ ] Workflows exist in the target commit (otherwise no publish)
- [ ] Pre-release checkbox matches intent
- [ ] Release notes are included
- [ ] CI is green

---

## Troubleshooting

### “The publish workflow didn’t run”
Common causes:
- The tag points to a commit where the workflow file didn’t exist yet.
- Release was created against an unexpected target branch.
- The workflow trigger is configured for a different release event/type.

Fix pattern:
1. Confirm the workflow YAML is present in the commit referenced by the tag.
2. If not, create a new tag on a newer commit that includes the workflow, and create a new GitHub Release pointing to it.

### “GitKraken doesn’t show a merge commit / PR relationship”
If you used **Squash & merge**, GitHub creates a single new commit on the base branch.
- That’s expected: no merge commit exists.
- The PR is still visible in GitHub history; GitKraken may show it as a normal commit.

---

## Team conventions (recommended defaults)

- PR merges:
  - Prefer **Squash & merge** for feature branches to keep history clean
  - Prefer **Merge commit** for `release/` and `hotfix/` if you want explicit branch structure
  - Pick one approach and document it (consistency matters more than choice)

- Always create releases via **GitHub UI** (not `git tag` alone), because publishing is coupled to GitHub Releases.

---

## Quick reference

### Final release to PyPI
1. `develop` → `release/X.Y.Z`
2. PR `release/X.Y.Z` → `master`
3. GitHub Release:
   - Tag `vX.Y.Z`, target `master`, not pre-release
4. PR `master` → `develop`

### Pre-release to TestPyPI
- GitHub Release:
  - Tag `vX.Y.Z-rc.1` (or `beta.1`)
  - Target `develop` (or `release/X.Y.Z`)
  - pre-release
