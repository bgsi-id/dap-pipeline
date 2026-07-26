# DAP Pipelines

Approved workflows executed by DAP through WES. A repository may contain many
pipelines; each pipeline is isolated in its own directory and registered in
`catalog.yaml`.

```text
pipelines/
  <pipeline-id>/
    pipeline.yaml
    main.nf
    nextflow.config
```

Runtime submissions always pin this repository to a full Git commit SHA and a
safe pipeline directory.
