# Start here

This folder is the first working project scaffold, not merely an audit.

## What you need to do now

1. Open `viralquantsim.Rproj` in RStudio.
2. Replace the placeholder email in `DESCRIPTION`.
3. Install the development dependencies:

```r
install.packages(c(
  "devtools", "testthat", "shiny", "bslib",
  "ggplot2", "plotly", "DT"
))
```

4. Run:

```r
devtools::load_all()
devtools::test()
```

5. Launch the interface:

```r
run_app()
```

The original ZIP remains reference material only. Do not copy its hard-coded vectors into this public repository.
