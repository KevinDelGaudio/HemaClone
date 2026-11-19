# Contributing to HemaClone

Rules for changing the code, configuration and documentation. Read
[`docs/CONTEXT.md`](docs/CONTEXT.md) before proposing a change: it describes
the design, the steps, the criteria, the outputs and the open items.

## Proposing a change

- Open an issue or describe the change before implementing it, with the reason and the expected effect on calls or outputs.
- One branch and one pull request per change; nothing is committed directly to `main`.
- Pull requests contain complete, working files, not partial edits.

## Code style

- Keep code simple and readable: straightforward constructs over clever one-liners or deep nesting.
- Preserve existing code unless a rewrite is the purpose of the change.
- Code, comments, documentation and commit messages in English.

## Workflow rules

- Every setting is configurable from the launch command; defaults live only in `config/config.yaml`.
- Launch examples are always complete, with every setting written out.
- Fail loudly: preflight checks on paths and resources, explicit errors on ID or join mismatches, row counts printed.
- Strict control of samples per batch and hard CPU and memory limits on the container.
- Steps are idempotent; a failed sample never stops the others and is listed.
- Never use a panel of normals built from the cohort being analysed: recurrent CH hotspots would be removed as artefacts.

## Versioning

- Every change to code or defaults bumps `VERSION` and adds a `CHANGELOG.md` entry.
  - Calls or classification change: at least a MINOR bump.
  - Every output identical: PATCH.
- `VERSION`, the git tag `vX.Y.Z` and the Docker image `hemaclone:X.Y.Z` always match. Release tags are created by the maintainer as GitHub releases after merge.
- Update `docs/CONTEXT.md` in the same pull request whenever steps, criteria, defaults or outputs change, and `README.md` when usage changes.
- The Docker tools layer stays stable: `workflow/docker/environment.yaml` changes only to change a tool version, and every package is pinned.

## Publishability

- No lab, cohort, patient or sample identifiers, collaborator names or site-specific paths in code, configuration, documentation, commit messages or pull requests. Sample lists and local paths stay on the analysis machine and are passed at run time.

## Performance claims

- No claim of better performance without a benchmark on shared data against the tools in [`docs/COMPARATORS.md`](docs/COMPARATORS.md).
- Current gaps are stated until they are closed: single caller; no strand-bias, repeat-region or within-cohort recurrence filters; sensitivity not yet measured.
- Sensitivity and precision are reported by VAF bin, with the truth set, data source and versions stated.

## Interpretation

- Calls are reported as somatic blood variants consistent with CH, not as CH outright.
- Technical causes are checked first (contamination, recurrent positions, indels in low-complexity regions, depth, library quality), then biological ones (tumour-derived DNA, therapy-related clones, mosaicism, rare germline variants).
