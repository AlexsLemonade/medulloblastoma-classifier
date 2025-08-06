# medulloblastoma-classifier

## Input data and models

The following table summarizes how different data types are handled by default.
If a method applies additional transformations, that is captured in the next table.

| Dataset type | Baseline normalization/processing |
|--------------|-----------------------------------|
| Microarray | [refine.bio processed with Single Channel Array Normalization](https://docs.refine.bio/en/latest/main_text.html#microarray-pipelines), quantile normalization is skipped |
| Bulk RNA-seq | TPM | 
| Smart-seq2 scRNA-seq | TPM |
| 10X scRNA-seq | Counts | 

This table summarizes the models used in this work, the packages from which they originate, and any transformations applied to the input gene expression measures.

| Model | Package | Additional transformations (if applicable) |
|-------|---------|-------|
| k-Top Scoring Pairs (kTSP) | [`multiclassPairs`](https://cran.r-project.org/web/packages/multiclassPairs/index.html) | N/A |  
| Random Forest (RF) | [`multiclassPairs`](https://cran.r-project.org/web/packages/multiclassPairs/index.html) | N/A | 
| MM2S ([Gendoo and Haibe-Kains. 2016.](https://doi.org/10.1186/s13029-016-0053-y)) | [`MM2S`](https://cran.r-project.org/src/contrib/Archive/MM2S/) | N/A |
| medulloPackage ([Rathi et al. 2020.](https://doi.org/10.1371/journal.pcbi.1008263)) | [`medulloPackage`](https://github.com/d3b-center/medullo-classifier-package/) | All RNA-seq data is log2-transformed |
| LASSO Logistic Regression | [`glmnet`](https://cran.r-project.org/web/packages/glmnet/index.html) | Each sample is scaled to sum to 1 | 


## Internal development guidelines

These guidelines are intended to be used by Data Lab members and collaborators.

### Dependency management

#### Docker

We expect development to primarily occur within the project Docker container.
We use renv and conda as part of the build process, so please make use of those approaches when updating the `Dockerfile` (see sections below).

A GitHub Actions workflow builds and pushes the Docker image to the GitHub Container Registry any time the relevant environment files or `Dockerfile` are updated.
It also checks on pull requests that alter relevant files that the image can be built.

To pull the most recent copy of the Docker image, use the following command:

```sh
docker pull ghcr.io/alexslemonade/medulloblastoma-classifier:latest
```

To run the container, use the following command from the root of this repository:

```sh
docker run \
  --mount type=bind,target=/home/rstudio/medulloblastoma-classifier,source=$PWD \
  -e PASSWORD={PASSWORD} \
  -p 8787:8787 \
  ghcr.io/alexslemonade/medulloblastoma-classifier:latest
```

Be sure to replace `{PASSWORD}`, including the curly braces, with a password of your choice.

You can then access the RStudio at <http://localhost:8787> using the username `rstudio` and the password you just set.

For Apple Silicon users, include the `--platform linux/amd64` tag in the `docker pull` and `docker run` commands.

#### Managing R packages with renv

We manage R package dependencies using [renv](https://rstudio.github.io/renv/articles/renv.html).

When you install additional packages, please update the lockfile with the following command:

```r
renv::snapshot()
```

When prompted, respond `y` to save the new packages in your `renv.lock` file.
Commit the changes to the `renv.lock` file.

To pin any packages that are not automatically captured in the lockfile, you can add loading them to the `dependencies.R` file in the root of the repository.

#### Managing command-line tools and Python packages with Conda

We use Conda to manage command-line tools and Python packages.

To create and activate the environment, run the following from the root of the repository (requires conda-lock to be installed):

```sh
conda-lock install --name medulloblastoma-classifier conda-lock.yml
conda activate medulloblastoma-classifier
```

To add new packages to the Conda environment, add them to `environment.yml`, and then update the `conda-lock.yml` file:

```sh
conda-lock --file environment.yml
```

#### Pre-commit

We use [pre-commit](https://pre-commit.com/) to make sure large files or secrets are not committed to the repository.
The Conda environment contains pre-commit.

To setup the pre-commit hooks for this project, run the following from the root of the repository:

```sh
pre-commit install
```

##### Additional hooks for local development

If you would like to add additional hooks to use locally (e.g., to style and lint R files), you can by creating and using a `.pre-commit-local.yaml` file like so:

```sh
# make and activate a local pre-commit configuration
cp .pre-commit-config.yaml .pre-commit-local.yaml
pre-commit install --config .pre-commit-local.yaml
```

`.pre-commit-local.yaml` is ignored by Git, so you can modify that file without affecting other contributors.

### Data and model management

We use an S3 bucket (`s3://data-lab-mb-ssp`) with versioning enabled to manage the files in the following directories:

- `data`
- `models`
- `processed_data`
- `plots/data`
- `results`

Which are all present in the `.gitignore` file.

To push files to S3, use the following command from the root of the repository:

```sh
aws s3 sync {directory} s3://data-lab-mb-ssp/{directory}
```

Where `{directory}` should be one of: `data`, `models`, `processed_data`, `plots/data`, or `results`.

To pull files locally, use the following command from the root of the repository:

```sh
aws s3 sync s3://data-lab-mb-ssp/{directory} {directory}
```

A non-exhaustive list of [`aws s3 sync` flags](https://awscli.amazonaws.com/v2/documentation/api/latest/reference/s3/sync.html) that may be useful:

- `--delete`: Delete files that exist in the destination that are not in the source.
- `--dryrun`: Performs a dry run without running the command.
- `--profile`: A profile from your credential file.
- `--exclude`: Exclude objects or files that match this pattern.
- `--include`: Don't exclude objects or files that match this pattern.
