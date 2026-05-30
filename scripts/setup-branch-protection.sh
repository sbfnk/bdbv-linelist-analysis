#!/usr/bin/env bash
#
# Enable branch protection on `main` for epiforecasts/bdbv-linelist-analysis.
#
# Requires the `gh` CLI authenticated with a token that has `repo` scope and
# admin rights on the repository. Run once (or after changing required status
# checks). Re-running is idempotent — the PUT call replaces the existing
# protection configuration.
#
# Usage:
#   ./scripts/setup-branch-protection.sh
#
# Protections enabled:
#   * Require pull requests before merging (1 approving review, dismiss stale).
#   * Require status checks to pass and branches to be up to date.
#     - "Documentation" (from .github/workflows/docs.yml).
#   * Require conversation resolution before merging.
#   * Block force-pushes and branch deletion.
#   * Enforce all rules on administrators.

set -euo pipefail

REPO="epiforecasts/bdbv-linelist-analysis"
BRANCH="main"

if ! command -v gh >/dev/null 2>&1; then
    echo "Error: gh CLI not found. Install from https://cli.github.com/." >&2
    exit 1
fi

if ! gh auth status >/dev/null 2>&1; then
    echo "Error: gh CLI not authenticated. Run 'gh auth login' first." >&2
    exit 1
fi

echo "Applying branch protection to ${REPO}@${BRANCH}..."

gh api \
    --method PUT \
    -H "Accept: application/vnd.github+json" \
    -H "X-GitHub-Api-Version: 2022-11-28" \
    "repos/${REPO}/branches/${BRANCH}/protection" \
    --input - <<'JSON'
{
  "required_status_checks": {
    "strict": true,
    "contexts": ["Documentation"]
  },
  "enforce_admins": true,
  "required_pull_request_reviews": {
    "dismiss_stale_reviews": true,
    "require_code_owner_reviews": false,
    "required_approving_review_count": 1,
    "require_last_push_approval": false
  },
  "restrictions": null,
  "required_linear_history": false,
  "allow_force_pushes": false,
  "allow_deletions": false,
  "block_creations": false,
  "required_conversation_resolution": true,
  "lock_branch": false,
  "allow_fork_syncing": false
}
JSON

echo "Branch protection applied successfully."
