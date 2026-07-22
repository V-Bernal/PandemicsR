# PandemicsR - R package for Epidemic Network Simulation
An R package for simulating epidemic dynamics on networks.
The simulation includes 3 layers: an opinions model, a Schelling model, and an Epidemic model.

---
## Architecture
Pipeline:
1. Generate a random network (structure, opinion, and initial epidemic states)
2. Generate network
3. Assign individual attributes/opinions
4. Run epidemic simulation
5. Analyze outcomes

Output:
1.  Graphs: PDF (`.pdf`, `.jpg`, `.jpeg`)

---
## Install R package from GitHub 
```{r}
install.packages("devtools")
library(devtools)
install_github("V-Bernal/PandemicsR", subdir="PandemicsR-main")

```
Run the script or shinny for testing. 

---
## Requirements
- R version
- Shinny
- Matrix
- 

