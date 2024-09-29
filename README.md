# medulloblastoma-classifier

## Internal Development Guidelines

These guidelines are intended to be used by Data Lab members and collaborators.

### Docker

We expect development to primarily occur within the project Docker container.
We use renv and conda as part of the build process, so please make use of those approaches when updating the `Dockerfile` (see sections below).

A GitHub Actions workflow builds and pushes the Docker image to the GitHub Container Registry any time the relevant environment files or `Dockerfile` are updated.
It also checks on pull requests that alter relevant files that the image can be built.

To pull the most recent copy of the Docker image, use the following command:

```sh
docker pull ghcr.io/AlexsLemonade/medulloblastoma-classifier:latest
```

To run the container, use the following command from the root of this repository:

```sh
docker run \
  --mount type=bind,target=/home/rstudio/medulloblastoma-classifier,source=$PWD \
  -e PASSWORD={PASSWORD} \
  -p 8787:8787 \
  ghcr.io/AlexsLemonade/medulloblastoma-classifier:latest
```

Be sure to replace `{PASSWORD}`, including the curly braces, with a password of your choice.

You can then access the RStudio at <http://localhost:8787> using the username `rstudio` and the password you just set.

### Managing R packages with renv

We manage R package dependencies using [renv](https://rstudio.github.io/renv/articles/renv.html).

When you install additional packages, please update the lockfile with the following command:

```r
renv::snapshot()
```

When prompted, respond `y` to save the new packages in your `renv.lock` file.
Commit the changes to the `renv.lock` file.

To pin any packages that are not automatically captured in the lockfile, you can add loading them to the `dependencies.R` file in the root of the repository.

### Managing command-line tools and Python packages with Conda

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

### Pre-commit

We use [pre-commit](https://pre-commit.com/) to make sure large files or secrets are not committed to the repository.
The Conda environment contains pre-commit.

To setup the pre-commit hooks for this project, run the following from the root of the repository:

```sh
pre-commit install
```

#### Additional hooks for local development

If you would like to add additional hooks to use locally (e.g., to style and lint R files), you can by creating and using a `.pre-commit-local.yaml` file like so:

```sh
# make and activate a local pre-commit configuration
cp .pre-commit-config.yaml .pre-commit-local.yaml
pre-commit install --config .pre-commit-local.yaml
```

`.pre-commit-local.yaml` is ignored by Git, so you can modify that file without affecting other contributors.
