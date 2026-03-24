# nf-core/crossmap portability report (for running on another machine)

## Executive summary

This pipeline is **not yet production-runnable end-to-end** on a new machine in normal mode because a required module (`BAM_TO_GFF`) is still a stub that exits with error.

You can still run a **smoke test** on another machine using `-stub-run` (dry run with fake outputs), and you can run many upstream real steps.  

---

## Current readiness snapshot

### Working components

- Main DSL2 orchestration is in place (`main.nf`, `workflows/crossmap.nf`).
- Most local preprocessing modules are implemented:
  - `FILTER_FEATURES`
  - `AGAT_LONGEST_ISOFORM`
  - `GET_NON_OVERLAPPING`
  - `GET_INTERGENIC_INTERVALS`
  - `AGAT_EXTRACT_SEQUENCES`
- nf-core modules for `minimap2`, `gffcompare`, `gffread`, `busco`, `multiqc` are wired.

### Blocking / incomplete components

1. **Hard blocker: `modules/local/bam_to_gff.nf`**
   - Contains placeholder logic and `exit 1`.
   - This step is mandatory in workflow, so full execution fails.

2. **Optional blocker: `modules/local/relocate_loci.nf`**
   - Placeholder logic and `exit 1`.
   - Only used when `--skip_decoy false` (default is false).  
   - Workaround: run with `--skip_decoy`.

3. **Optional future components: annocli**
   - `annocli_download.nf` and `annocli_alias.nf` are placeholders.
   - Currently not active in main path because `PREPARE_GENOME` is commented out in `workflows/crossmap.nf`.

4. **Test profile data mismatch**
   - `conf/test.config` points to `assets/samplesheet.csv`.
   - `assets/samplesheet.csv` contains placeholder paths (`/path/to/...`) that will not exist on another machine.

5. **Docs still template-like**
   - `README.md` and `docs/usage.md` still contain nf-core template TODO sections and non-crossmap examples.

---

## What to do on another machine right now

## 1) Environment prerequisites

- Java 17+ (required by Nextflow)
- Nextflow `>= 24.04.2`
- One runtime profile:
  - recommended: `docker`
  - alternatives: `singularity` / `apptainer` / `podman` / `conda`
- Internet access for container/data pulls (unless pre-cached)

## 2) Clone and pin code

```bash
git clone https://github.com/nf-core/crossmap.git
cd crossmap
git checkout <commit-or-tag-you-trust>
```

## 3) Create a real samplesheet (required)

Use real, existing files (absolute paths recommended):

```csv
sample,taxid,role,assembly,annotation
Species_A,12345,both,/abs/path/species_a.fna,/abs/path/species_a.gff3
Species_B,67890,target,/abs/path/species_b.fna,
```

Rules:
- `role` must be `source`, `target`, or `both`
- `assembly` required for all rows
- `annotation` required for `source`/`both`

## 4) Run a smoke test (portable sanity check)

Because of current stubs, use:

```bash
nextflow run main.nf \
  -profile docker \
  --input /abs/path/samplesheet.csv \
  --outdir /abs/path/results_stub \
  --skip_decoy \
  --skip_busco \
  --skip_gffcompare \
  -stub-run
```

This validates orchestration portability (channels, module wiring, publish paths) but not biological correctness.

---

## What must be implemented before real end-to-end runs

1. **Implement `BAM_TO_GFF` (highest priority)**
   - Convert spliced BAM alignments to valid GFF3 transcripts/features.
   - This is the main gate to operational pipeline usage.

2. **Either implement `RELOCATE_LOCI` or keep `--skip_decoy` default true**
   - If decoy controls are required scientifically, implement relocation logic.

3. **Decide annocli strategy**
   - If download/alias normalization is needed, implement annocli modules and re-enable `PREPARE_GENOME`.
   - If not needed, keep download disabled and rely on user-provided files.

4. **Fix `conf/test.config` and assets test data**
   - Replace placeholder samplesheet paths with usable tiny test dataset.

5. **Clean docs for transferability**
   - Update `README.md` and `docs/usage.md` with real crossmap input/output instructions.

---

## Recommended short-term run policy for portability

Until blockers are fixed, use this policy on any machine:

- Always provide your own samplesheet with real paths.
- Always run with:
  - `--skip_decoy`
  - `--skip_busco` (or set real `--busco_lineage`)
  - `--skip_gffcompare` for simpler first portability checks
- Use `-stub-run` for CI/smoke portability checks.

---

## Suggested params file for another machine

Create `params_portable.yaml`:

```yaml
input: "/abs/path/samplesheet.csv"
outdir: "/abs/path/results"
feature_types: "lnc_RNA,mRNA"
minimap2_preset: "splice"
skip_decoy: true
skip_busco: true
skip_gffcompare: true
skip_download: true
```

Run:

```bash
nextflow run main.nf -profile docker -params-file params_portable.yaml -stub-run
```

Remove `-stub-run` only after implementing `BAM_TO_GFF` (and any other required stubs).
