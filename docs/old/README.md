# Create PDF from PIPELINE.rmd

PIPELINE.rmd is R notebook file. To create PDF from it you will need to install **R** and **rmarkdown** library for it.

After installation enter R console and folder containing PIPELINE.rmd and related files.
Run following code in R console.

```
library("rmarkdown")
render("PIPELINE.rmd", "all")
```
 
 This code will convert PIPELINE.rmd into all formats (pdf only at the moment) listed in it.
