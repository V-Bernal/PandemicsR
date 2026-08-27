# PandemicsR - R package for Epidemic Network Simulation
An R package for simulating epidemic dynamics on networks.
The simulation includes 3 layers: an opinion model, a Schelling model, and an Epidemic model.

## Documentation and Vignettes

The package includes detailed documentation in the form of R Markdown vignettes:

- **Introduction** – overview of the model and its components.
- **Model description** – mathematical and conceptual description of the simulation.
- **Parameters** – description of all user-defined and internal parameters.
   
After installing the package, the vignettes can be accessed in R with:

```r
browseVignettes("PandemicsR")
```
---
## Install R package from GitHub 
```{r}
install.packages("devtools")
library(devtools)
install_github("V-Bernal/PandemicsR", subdir="PandemicsR-main")
# Run the script or Shiny for testing. 
library(PandemicsR)
shiny::runApp()
```
## Architecture
Pipeline:
1. Generate a random network (structure, opinion, and initial epidemic states)
2. Generate network
3. Assign individual attributes/opinions
4. Run epidemic simulation
5. Analyze outcomes

Output:
Graphs: PDF (`.pdf`, `.jpg`, `.jpeg`)

## Project Structure
```text
PandemicsR-main/
├── R/
│   ├── network/
│   ├── opinions/
│   ├── groups/
│   ├── epidemic/
│   ├── stability.R
│   ├── gillespie.R
│   ├── run_simulation.R
│   └── ...
│
├── vignettes/
│   ├── introduction.Rmd
│   ├── model_description.Rmd
│   └── Parameters.Rmd
│
├── app.R
├── DESCRIPTION
├── NAMESPACE
└── README.md
```
## Model interaction
```text
                 ┌─────────────────┐
                 │     Network     │
                 └────────┬────────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │     Opinions    │
                 │   Voter model   │
                 └────────┬────────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │      Groups     │
                 │ Schelling model │
                 └────────┬────────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │     Epidemic    │
                 │    SIR model    │
                 └─────────────────┘

       Social state ───────────────► Epidemic transmission
       Epidemic state ─────────────► Social dynamics

```
## Requirements
- `R >= 4.3.0`
- `igraph`
- `Matrix`
- `shiny`
- `zip`
---
