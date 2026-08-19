# AGENTS.md — Configuration for AI Coding Agents

This file is the **canonical source of truth** for all AI coding agents
(including Google Jules and GitHub Copilot) working on the `simcyto` repository.
Copilot-specific files must defer to this file rather than becoming a separate
source of repository guidance.

> **Instruction update rule for AI agents:**
> When a shared coding pattern, guideline or best practice needs to change,
> update `AGENTS.md`. Do not duplicate or maintain competing shared instructions
> in Copilot-specific configuration.

## 1. Project context

`simcyto` is an R package for simulating cytometry data for testing and
benchmarking. Keep the package focused on reusable simulation machinery rather
than downstream statistical evaluation or benchmarking logic.

The main simulation hierarchy is:

`simCytExperiment()` -> `simCytSample()` -> `simCytCondition()` -> `simCytCluster()`.

Scenario builders provide the marker-expression, cluster and response settings
consumed by the simulation engine. The package also supports transformations,
perturbations and ratio correction used by the simulator.

## 2. Repository structure

- `R/`: package implementation.
- `tests/testthat/`: package tests.
- `DESCRIPTION`: package metadata and dependencies.
- `NAMESPACE`: generated/exported namespace information.
- `.github/workflows/copilot-setup-steps.yml`: GitHub Copilot cloud-agent setup.
- `.github/copilot-instructions.md`: thin Copilot pointer back to this file.

## 3. Environment setup and dependency management

Run package commands from the repository root, the directory containing
`DESCRIPTION`.

For ordinary development, install the package dependencies and development
tools, then load the package with:

```r
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::local_install_deps(dependencies = TRUE)

devtools::load_all()
```

The package depends on Bioconductor cytometry packages, including `flowCore`
and `Biobase`. If dependency installation fails, fix the dependency or system
setup rather than working around the failure in package code.

### GitHub Copilot cloud agent

The Copilot cloud agent is prepared by
`.github/workflows/copilot-setup-steps.yml`. The workflow installs R, system and
package dependencies, and the development tools needed for package work. It
uses Posit Public Package Manager where possible and caches successfully built
package libraries so dependencies do not need to be rebuilt from source on
every run.

In the Copilot agent, use the prepared development library directly:

```r
devtools::load_all()
devtools::test()
```

If expected packages are unavailable there, treat that as an infrastructure
setup failure. Fix `.github/workflows/copilot-setup-steps.yml` rather than
adding package-code workarounds.

## 4. Build, test and quality checks

Useful commands from the repository root are:

```r
devtools::load_all()
devtools::test()
devtools::document()
styler::style_pkg()
lintr::lint_package()
devtools::check()
```

Before opening or updating a pull request, run the checks relevant to the files
you changed. For package-code changes, this should normally include at least
`devtools::test()` and `devtools::check()`. Regenerate roxygen documentation
with `devtools::document()` when documentation or exports change.

Do not edit generated `.Rd` files manually.

## 5. Coding conventions

- Preserve the existing `simCyt*` naming convention for exported simulation
  functions.
- Keep internal helpers internal unless there is a clear package API reason to
  export them.
- Keep simulation machinery in `simcyto`; do not move downstream method
  evaluation or benchmarking analysis into the package merely for convenience.
- Add or update tests for behavioural changes.
- Prefer explicit namespace qualification for external package calls unless the
  package intentionally imports a symbol.
- Keep files free of trailing whitespace and end text files with a final newline.

## 6. Instruction authority

`AGENTS.md` is authoritative for shared AI-agent instructions. The
`.github/copilot-instructions.md` file exists only so GitHub Copilot is directed
to this file. If the two ever appear to disagree, follow `AGENTS.md` and update
the Copilot file so that it simply points here.
