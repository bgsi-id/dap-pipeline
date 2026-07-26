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

Auth synchronises `catalog.yaml` from the configured branch into its own
persistent registry. Repository `active` controls whether an entry is eligible
for import; a separate platform enable/disable flag controls whether users can
launch it. Synchronisation updates metadata without overriding a platform
disable.
