# Setup checklist

Everything that still needs a human. Delete this file once the repository is live.

## 1. Fill in the placeholders

| File | Placeholder | What to put there |
|---|---|---|
| `CITATION.cff` | `date-released` | the date you cut the GitHub release |
| `CITATION.cff` | `doi` | the Zenodo DOI (step 4) |
| `README.md` | "How to cite" block | the same DOI, and the manuscript citation once it has one |

Owner (`KarlisPleiko`), repository name (`differential-homing`), author (Karlis Pleiko, ORCID
0000-0002-1073-197X) and the MIT licence are already filled in. `CITATION.cff` also
cross-references phader as a related work — add its DOI there once both are minted.

## 2. Check before the first push

- [ ] Confirm `app.R` is the exact file used for the manuscript analysis. Its sha256 as
      published here is
      `f87d96a595c25aac377c9582f4601d311000f551f8e0174b0cf241a51a77cbb5`
      (`shasum -a 256 app.R`).
- [ ] Read the two launcher changes in `run.R` and confirm you want them. They are the
      `127.0.0.1` bind and the non-upgrading Bioconductor install; both are documented in
      KNOWN_ISSUES H12 and H13. If you would rather publish `run.R` byte-identical too, revert
      them and move both entries back into the issue list.
- [ ] **Decide how to handle H1 for the manuscript.** Depleted peptides are drawn as "Not
      significant" in both volcanoes. If a figure from this version is in the paper, either
      state in the caption that only enrichment is highlighted, or annotate the depleted
      peptides from the results CSV. This is the one item that affects what a reader concludes.
- [ ] Run `Rscript scripts/capture_renv.R` on the machine that produced the published figures
      and commit `renv.manuscript.lock`. The shipped `renv.lock` describes the validation
      environment, not yours.
- [ ] Confirm no data is staged: `git status --short` should list only text files, the example
      count tables, the reference sheet and the four example PNGs. `.gitignore` excludes
      `*.fastq`, `*_results.txt`, `*.xlsx`, `data/` and `results/`, with explicit exceptions
      for the files under `examples/`.

## 3. Create the repository

On github.com: **New repository**, owner `KarlisPleiko`, name `differential-homing`, Public, and **do
not** tick "Add a README", ".gitignore" or "license" — this bundle already contains all three.

```bash
cd differential-homing
git init
git add .
git commit -m "differential-homing v1.0.0: differential homing analysis for phage display peptide counts"
git branch -M main
git remote add origin git@github.com:KarlisPleiko/differential-homing.git
git push -u origin main
git tag -a v1.0.0 -m "Version used for the manuscript analysis"
git push origin v1.0.0
```

Set the repository description and topics (`phage-display`, `biopanning`, `edger`,
`differential-abundance`, `shiny`, `peptide-library`) so the tool is findable.

If you upload through the browser instead: Finder hides `.gitignore`, so add it afterwards with
**Add file → Create new file**; and drag the `docs`, `examples`, `scripts` folders themselves
onto the drop zone so their paths are preserved.

## 4. Mint a DOI

1. In your Zenodo GitHub settings, switch `differential-homing` **On**. Do this *before* publishing the
   release — a release published while the toggle is off is never archived.
2. On GitHub: **Releases → Draft a new release**, tag `v1.0.0`, title `differential-homing v1.0.0`,
   release label **None** (not pre-release), publish.
3. Wait a minute or two, reload your Zenodo GitHub settings, and copy the DOI.
4. Paste it into `CITATION.cff` and `README.md`, commit, push.

Cite the **version** DOI in the manuscript, not the concept DOI, so readers land on exactly
this code.

## 5. Cross-link the two tools

Once both repositories are live and both DOIs exist:

- [ ] Add the differential-homing link and DOI to phader's README, and the phader link and DOI to this
      README, so a reader arriving at either finds the other.
- [ ] Add phader's DOI to the `references` block in `CITATION.cff`.
- [ ] In the manuscript's code-availability statement, cite both version DOIs and say which
      tool did which part: phader counts peptides from the FASTQ files, differential-homing performs
      the differential abundance analysis on those counts.

## 6. Next release

The two tools share a file format, so their v2 work is coupled — see the Planned section of
[CHANGELOG.md](CHANGELOG.md). Do not change phader's output format without changing this app's
reader in the same release; the five-line header skip is a hard dependency (H8).
