This is a simple R function to do general data psychometric analysis based on the SEM framework.

Upload the file to your R session with
``` 
source("Your Drive\\psicomAgh.R")
``` 
Simple usage instructions:
``` 
 psicomA <- function(model, bd, ordered = FALSE, stdlv = TRUE, vd = TRUE, sw = NULL, ...) {
``` 
   
   ### Function arguments
   
   model   : lavaan model
   
   bd      : data.base (data.frame)
   
   ordered : logical — Ordinal variables?
   
   stdlv   : logical — standardize lanten variables?
   
   vd      : locical  - Calculated discriminant/convergent validity?
   
   sw      : Sampling weights (opcional)
 
  
### Example:

### Lavaan Model
``` 
mod <- "
         F1 =~ INT1 + INT2 + INT3 + PER1 + PER2
         F2 =~ ECO2 + ECO3 + PROF1 + PROF2 + PER3
          "
``` 
### Select items from mod object
``` 
items <- c("INT1", "INT2", "INT3", "PER1", "PER2", "ECO2", "ECO3", "PROF1", "PROF2", "PER3")
``` 
### Upload function
``` 
source("D:\\OneDrive - Ensino Lusófona\\psicomA.R")          # Replace the "D:\\OneDrive - Ensino Lusófona" by your own
``` 
### Fit the model
``` 
fitM <- psicomA(model=mod,bd=bd[,items], ordered=T,vd=T, sw=NULL, stdlv=T)
``` 
### See results
``` 
fitM       ### Works best in R Studio > Knitr > to HTML
``` 
### To retrieve the lavaan object (for other analysis
``` 
fitM[[1]]
``` 
### To retrieve individual results objects
 ``` 
fit$Sensibilidade             ##### Descriptive statistics for items
fit$PesosFatoriais            ##### Std factor loadings table
fit$IndicesGOF                #### GOF tabel
fit$Fiabilidade               #### Reliability table 
fit$ValidadeConvergente       #### Convergent Validity table
fit$ValidadeDiscriminante     ### Discriminant validity table
fit$PathPlot
```

Questions, bugs, improvements: use my email jpmaroco at gmail dot com

To cite (APA format):

Marôco, J. (2024). psicomA: Psychometric analysis tools in R (Version 1.0) [Software]. GitHub. https://github.com/joaomaroco/psicomA

