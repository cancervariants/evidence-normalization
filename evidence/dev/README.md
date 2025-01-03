# Dev Module

Intended for developer purposes

## CLI

### Transforming Cancer Hotspots data

Requires Variation Normalization to be configured. See the [README](https://github.com/cancervariants/variation-normalization#variation-normalization) for more information.
```commandline
python3 -m evidence.dev.cli --transform_cancer_hotspots
```

### Transforming cBioPortal data

```commandline
python3 -m evidence.dev.cli --transform_cbioportal
```

### Transform all source data

```commandline
python3 -m evidence.dev.cli --transform_all
```
