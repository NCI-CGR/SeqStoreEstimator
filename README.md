# SeqStoreEstimator

## Live

See the live version here:

https://appshare-dev.cancer.gov/content/2e91bacb-7334-416f-8159-03f8d366ea3a


## run app

```bash
poetry run python -m shiny run --reload
```

## Dev

If you haven't already, install the export plugin:

```bash
poetry self add poetry-plugin-export
```

Update requirements.txt if packages are updated/changed:

```bash
poetry export -f requirements.txt --output requirements.txt --without-hashes
```

With each commit update the manifest

```bash
poetry run rsconnect write-manifest shiny --overwrite .
```