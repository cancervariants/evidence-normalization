# Evidence Normalization

Service for normalizing evidence

## Installation

Install from [PyPI](https://pypi.org/project/evidence-normalizer):

```shell
python3 -m pip install evidence-normalizer
```

To install ETL dependencies:

```shell
python3 -m pip install evidence-normalizer[etl]
```

## Development

Clone the repo and create a virtual environment:

```shell
git clone https://github.com/cancervariants/evidence-normalization
cd evidence-normalization
python3 -m virtualenv venv
source venv/bin/activate
```

Install development dependencies and `pre-commit`:

```shell
python3 -m pip install -e '.[dev,etl,tests]'
pre-commit install
```

### Backend Services

Evidence Normalization relies on [Variation Normalization](https://github.com/cancervariants/variation-normalization) for normalizing Cancer Hotspots data. You will need to setup backend services and set the appropriate environment variables. See the [README](https://github.com/cancervariants/variation-normalization#variation-normalization) for more information.
