#!/bin/bash
# Module 5: put the project under version control so it can go on GitHub (for the CV)
# Only code, metadata, software versions and small results are tracked.
# Raw data, the reference and Salmon outputs are far too large for GitHub
# and can be regenerated from the scripts, so they are excluded.
# Run once:  bash ~/rnaseq_qNSC/scripts/11_github.sh
# Time: seconds (pushing to GitHub is a separate step, see the end)
set -euo pipefail
cd ~/rnaseq_qNSC

cat > .gitignore << 'EOF'
data/
ref/
results/salmon/*/
results/deseq2/*.rds
*.fastq.gz
*.sra
EOF

if [ ! -d .git ]; then git init; fi
git add .gitignore scripts metadata environment*.yml results README.md 2>/dev/null || \
git add .gitignore scripts metadata environment*.yml results
git status --short
echo
echo "Next, once: set your name and email for git, then commit:"
echo '  git config --global user.name  "Romi Marina George"'
echo '  git config --global user.email "your-email@example.com"'
echo '  git commit -m "Bulk RNA-seq analysis of GSE116997 Experiment B"'
echo "Then create an empty repository on github.com and follow its 'push an existing repository' instructions."
