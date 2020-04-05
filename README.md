A list of r package r requiered for run this app in a desk:

| Package         | Version |
| --------------- | :-----: |
| data.table      | 1.12.8  |
| DESeq2          | 1.22.2  |
| DT              |  0.11   |
| EnhancedVolcano |  1.5.4  |
| plotly          |  4.9.1  |
| reshape2        |  1.4.3  |
| shiny           |  1.4.0  |
| shinydashboard  |  0.7.1  |
| tidyverse       |  1.3.0  |

1. A easy way to check/install this package in your R session, run:

```r
.cran_packages <- c('shiny', 'plotly', 'data.table', 'reshape2', 'DT', 'tidyverse', 'shinydashboard')

.bioc_packages <- c('DESeq2', 'EnhancedVolcano')

.inst <- .cran_packages %in% installed.packages()

if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}

.inst <- .bioc_packages %in% installed.packages()

if(any(!.inst)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(.bioc_packages[!.inst], ask = F)
}

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)


# Check version
ip <- installed.packages()[,c(1,3)]
rownames(ip) <- NULL
ip <- ip[ip[,1] %in% c(.cran_packages,.bioc_packages),]
print(ip, row.names=FALSE)



```

2. After check the package installation downolad the `app.R` file and run from an r session
3. A web-app version can be used from [here](https://rjhgore-dc.shinyapps.io/Exploratory_analisys/). You can make a _fake_ dataset and run the app on the web:

```r
wd <- '~/'

dds2 <- makeExampleDESeqDataSet(n=10000, m=12)

# Create fake experimental design

m <- data.frame(colData(dds2))
m[,2] <- rep(paste0('rep',1:3), each = 1, times = 4)

design <- rep(0:3, each = 3, times = 1)
design <- paste("G", design, sep="") 

m[,3] <- design
m[,4] <- colnames(dds2)

colnames(m) <- c('condition', 'Rep', 'group', 'names')


saveRDS(dds2, file = paste0(wd, "dds_example.rds"))
write.table(m, file = paste0(wd, 'colData_example.csv'))

#file.size(paste0(wd, "dds_example.rds")) / 10^6
#as.numeric(object.size(dds2) / 10^6)

```



