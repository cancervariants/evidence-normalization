# Evidence Normalization

Service for normalizing evidence

## Developer instructions

The following sections include instructions specifically for developers.

### Installation

#### Pipenv
For a development install, we recommend using Pipenv. See the 
[pipenv docs](https://pipenv-fork.readthedocs.io/en/latest/#install-pipenv-today) 
for direction on installing pipenv in your compute environment.
 
Once installed, from the project root dir, just run:

```commandline
pipenv shell
pipenv lock && pipenv sync
```

#### Pip

If you wish to install developer dependencies for `evidence.dev`:
```commandline
pip install evidence-normalizer[dev]
```

If you do not need the extra dependencies:
```commandline
pip install evidence-normalizer
```

### Backend Services

Evidence Normalization relies on [Variation Normalization](https://github.com/cancervariants/variation-normalization) for normalizing Cancer Hotspots data. You will need to setup backend services and set the appropriate environment variables. See the [README](https://github.com/cancervariants/variation-normalization#variation-normalization) for more information.


### Starting the Evidence Normalization Service Locally

To start the service, run the following:

```commandline
uvicorn evidence.main:app --reload
```

Next, view the OpenAPI docs on your local machine:
http://127.0.0.1:8000/evidence

### Init coding style tests

Code style is managed by [flake8](https://github.com/PyCQA/flake8) and checked prior to commit.

We use [pre-commit](https://pre-commit.com/#usage) to run conformance tests.

This ensures:

* Check code style
* Check for added large files
* Detect AWS Credentials
* Detect Private Key

Before first commit run:

```commandline
pre-commit install
```


### Running unit tests

Running unit tests is as easy as pytest.

```commandline
pipenv run pytest
```
