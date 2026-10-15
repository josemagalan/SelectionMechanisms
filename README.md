![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
[![Shiny App](https://img.shields.io/badge/Shiny-Live_App-brightgreen)](https://josemagalan.shinyapps.io/SelectionMechanisms/)
![R](https://img.shields.io/badge/R-4.0%2B-blue)
![Purpose](https://img.shields.io/badge/Purpose-Teaching_Tool-informational)

# SelectionMechanisms

Interactive Shiny application to explore classical selection mechanisms
in Genetic Algorithms and analyse their impact on convergence dynamics
and genetic diversity.

🌐 **Live application:**\
https://josemagalan.shinyapps.io/SelectionMechanisms/


## Authors

-   José Manuel Galán
-   Silvia Díaz-de La Fuente
-   Virginia Ahedo
-   María Pereda
-   José Ignacio Santos


## Overview

Selection pressure is one of the key drivers of evolutionary dynamics in
Genetic Algorithms.
Different selection operators induce different convergence speeds and
diversity loss patterns.

This interactive Shiny application allows users to:

-   Experiment with multiple classical selection mechanisms
-   Adjust selection parameters
-   Run multiple stochastic simulations
-   Visualise convergence and diversity trajectories
-   Compare models side by side

The tool is designed for teaching, experimentation and conceptual
understanding of evolutionary search processes.


## Implemented Selection Methods

-   **Tournament selection**
    -   Adjustable tournament size (k)
    -   Optional sampling without replacement
-   **Roulette Wheel selection**
    -   Optional linear scaling (f' = a·f + b)
-   **Stochastic Universal Sampling (SUS)**
    -   Optional linear scaling
-   **Rank selection (Baker linear ranking)**
    -   Adjustable selection pressure parameter
-   **Truncation selection**
    -   Adjustable top proportion
-   **Boltzmann selection**
    -   Adjustable temperature parameter


## Metrics and Visualisations

For each simulation and generation, the app computes:

-   Number of unique individuals (richness)
-   Shannon diversity index
-   Mean fitness
-   Maximum fitness

The interface includes:

-   Diversity plots (mean with percentile bands)
-   Fitness dynamics plots
-   Summary statistics
-   Multi-model comparison panel


## Pedagogical Purpose

The application is designed as a classroom laboratory to:

-   Illustrate the concept of selection pressure
-   Compare stochastic sampling variance (Roulette vs SUS)
-   Analyse trade-offs between exploitation and diversity preservation
-   Understand how parameterisation affects convergence dynamics
-   Connect evolutionary computation concepts with probabilistic
    modelling

It is suitable for undergraduate and graduate courses in:

-   Evolutionary Computation
-   Metaheuristics
-   Artificial Intelligence
-   Operations Research
-   Industrial Engineering
-   Stochastic Optimization


## Running Locally

Clone the repository and run:

``` r
shiny::runApp()
```

Required R packages:

-   shiny
-   dplyr
-   purrr
-   tibble
-   ggplot2
-   rlang


## Deployment

The application is deployed using **shinyapps.io** and is publicly
accessible at:

https://josemagalan.shinyapps.io/SelectionMechanisms/


<!-- ga-series-links:start -->
## Teaching resources / Materiales docentes

[Portal en español](https://josemagalan.github.io/genetic-algorithms-lab/?lang=es) · [Portal in English](https://josemagalan.github.io/genetic-algorithms-lab/?lang=en)

| Resource / Recurso | Use or read / Usar o leer | Repository / Repositorio |
| --- | --- | --- |
| Portal: Genetic Algorithms Lab / Laboratorio de algoritmos genéticos | [Open / Abrir](https://josemagalan.github.io/genetic-algorithms-lab/) | [GitHub](https://github.com/josemagalan/genetic-algorithms-lab) |
| Representation / Representación | [Open / Abrir](https://josemagalan.github.io/representation-ga/) | [GitHub](https://github.com/josemagalan/representation-ga) |
| Initialisation and evaluation / Inicialización y evaluación | [Open / Abrir](https://josemagalan.github.io/initialization-and-evaluation-ga/) | [GitHub](https://github.com/josemagalan/initialization-and-evaluation-ga) |
| Selection / Selección | [Open / Abrir](https://josemagalan.github.io/selection-ga/) | [GitHub](https://github.com/josemagalan/selection-ga) |
| Crossover / Cruce | [Open / Abrir](https://josemagalan.github.io/crossover-ga/) | [GitHub](https://github.com/josemagalan/crossover-ga) |
| Mutation / Mutación | [Open / Abrir](https://josemagalan.github.io/mutation-ga/) | [GitHub](https://github.com/josemagalan/mutation-ga) |
| EvoTraveller | [Open / Abrir](https://evotraveller.streamlit.app/) | [GitHub](https://github.com/jismartin/evotraveller) |
| SelectionMechanisms | [Open / Abrir](https://josemagalan.shinyapps.io/SelectionMechanisms/) | [GitHub](https://github.com/josemagalan/SelectionMechanisms) |
| Book / Libro: *Evolución a la vista* | [Read / Leer](https://josemagalan.github.io/genetic-algorithms-lab-book/) | [GitHub](https://github.com/josemagalan/genetic-algorithms-lab-book) |

El libro está en construcción: su contenido y estructura pueden cambiar. / The book, in Spanish, is a work in progress: its content and structure may change.
<!-- ga-series-links:end -->

## License

This project is released under the MIT License.



