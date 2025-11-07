# SeqStoreEstimator
## run app

```bash
poetry run python -m shiny run --reload
```

## Dev

Update requirements.txt if packages are updated/changed:

```bash
poetry export -f requirements.txt --output requirements.txt
```

With each commit update the manifest

```bash
poetry run rsconnect write-manifest shiny .
```