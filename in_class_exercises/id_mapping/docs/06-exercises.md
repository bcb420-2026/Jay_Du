# Exercises

These are intended to be done **after** completing the worked examples.

## Exercise 1 — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119732

Using **GSE119732**, confirm whether the ID column contains Ensembl IDs with version suffixes.

1. Extract the first 20 IDs.
2. Count how many contain a `.`.
3. Create a new column with versions stripped.
4. Map the identifiers to HGNC symbols.

### Helper functions


``` r
# Load Library
library(tibble)
library(readr)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

``` r
# Print first n from head
kable_head <- function(x, n = 5, caption = NULL) {
  knitr::kable(utils::head(x, n), caption = caption)
}

# Safe read files
safe_read <- function(file) {
  # First attempt: read as TSV
  df <- tryCatch(
    readr::read_tsv(file, show_col_types = FALSE),
    error = function(e) NULL   # catch fatal errors
  )
  
  # If read_tsv failed entirely:
  if (is.null(df)) {
    message("TSV read failed — reading as space-delimited file instead.")
    return(readr::read_table(file, show_col_types = FALSE))
  }
  
  # If read_tsv returned but with parsing issues:
  probs <- problems(df)
  if (nrow(probs) > 0) {
    message("Parsing issues detected in TSV — reading as space-delimited file instead.")
    return(readr::read_table(file, show_col_types = FALSE))
  }
  
  # If everything was fine:
  return(df)
}

# Strip version number
strip_ensembl_version <- function(ids) sub("\\..*$", "", ids)


# Load fetching function
source("./fetch_geo_supp.R")
```


### Question 1 

``` r
gse <- "GSE119732"
pattern <- "GSE119732.*\\.gz"

path <- file.path("data", gse)
files <- list.files(path, pattern = pattern, 
                    full.names = TRUE, recursive = TRUE)

kable_head(tibble(file = basename(files)), min(10, length(files)), paste(gse,": extracted file list (first 10)"))
```



Table: (\#tab:unnamed-chunk-2)GSE119732 : extracted file list (first 10)

|file                                 |
|:------------------------------------|
|GSE119732_count_table_RNA_seq.txt.gz |

``` r
x <- safe_read(files[1])

kable_head(x[, 1:min(6, ncol(x))], 20, paste(gse,": raw table preview"))
```



Table: (\#tab:unnamed-chunk-2)GSE119732 : raw table preview

|gene_id           |  A1|  A2|  A3|  A4|  B1|
|:-----------------|---:|---:|---:|---:|---:|
|ENSG00000223972.5 |   0|   0|   0|   0|   0|
|ENSG00000227232.5 |  79| 119|  84|  50|  80|
|ENSG00000278267.1 |  17|  10|  22|  19|  19|
|ENSG00000243485.4 |   0|   0|   0|   0|   0|
|ENSG00000237613.2 |   0|   0|   0|   0|   0|
|ENSG00000268020.3 |   0|   0|   0|   0|   0|
|ENSG00000240361.1 |   0|   0|   0|   0|   0|
|ENSG00000186092.4 |   0|   0|   0|   0|   0|
|ENSG00000238009.6 |   0|   0|   0|   1|   1|
|ENSG00000239945.1 |   0|   0|   0|   0|   0|
|ENSG00000233750.3 |   0|   0|   0|   0|   0|
|ENSG00000268903.1 |   7|  10|   9|   6|  29|
|ENSG00000269981.1 |  17|  11|  14|   4|  27|
|ENSG00000239906.1 |   0|   0|   1|   0|   0|
|ENSG00000241860.6 |   6|  19|  12|   8|  60|
|ENSG00000222623.1 |   0|   0|   0|   0|   0|
|ENSG00000241599.1 |   0|   0|   0|   0|   0|
|ENSG00000279928.1 |   8|   5|   9|   6|   0|
|ENSG00000279457.3 | 282| 298| 334| 205| 156|
|ENSG00000273874.1 |   0|   0|   0|   0|   0|


### Question 2

``` r
# Extract the first 20 IDs
first_20_ids <- head(x$gene_id, 20)

# Count how many of the first 20 contain a "."
dot_count <- sum(grepl("\\.", first_20_ids))
print(paste("Number of IDs with version suffixes in first 20:", dot_count))
```

```
## [1] "Number of IDs with version suffixes in first 20: 20"
```

``` r
# Count how many total contain a "."
dot_count <- sum(grepl("\\.", x$gene_id))
print(paste("Total number of IDs with version suffixes:", dot_count))
```

```
## [1] "Total number of IDs with version suffixes: 58037"
```

All 20 from the 1st 20 have a ".", all 58037 ids of this whole set
have a ".".


### Question 3 

``` r
# Apply to your tibble 'x'
id_col <- names(x)[1]
ids <- x[[1]] |> as.character()
kable_head(tibble(raw_id = head(ids, 10), 
                  stripped = strip_ensembl_version(head(ids, 10))), 
           10, paste(gse,": ID preview"))
```



Table: (\#tab:unnamed-chunk-4)GSE119732 : ID preview

|raw_id            |stripped        |
|:-----------------|:---------------|
|ENSG00000223972.5 |ENSG00000223972 |
|ENSG00000227232.5 |ENSG00000227232 |
|ENSG00000278267.1 |ENSG00000278267 |
|ENSG00000243485.4 |ENSG00000243485 |
|ENSG00000237613.2 |ENSG00000237613 |
|ENSG00000268020.3 |ENSG00000268020 |
|ENSG00000240361.1 |ENSG00000240361 |
|ENSG00000186092.4 |ENSG00000186092 |
|ENSG00000238009.6 |ENSG00000238009 |
|ENSG00000239945.1 |ENSG00000239945 |


### Question 4

``` r
if (any(grepl("^ENSG", strip_ensembl_version(ids)))) {
  library(biomaRt)
  ensembl_ids <- unique(strip_ensembl_version(ids))
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  map <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
  
  kable_head(map, 10, paste(gse,": mapping preview"))
  
  expr_mapped <- x %>%
    mutate(ensembl_gene_id = strip_ensembl_version(.data[[id_col]])) %>%
    left_join(map, by = "ensembl_gene_id") %>%
    dplyr::select( ensembl_gene_id,hgnc_symbol, everything())
  
  
  kable_head(expr_mapped[, 1:min(8, ncol(expr_mapped))], 5, paste(gse,": mapped table preview"))
}
```

```
## Warning in left_join(., map, by = "ensembl_gene_id"): Detected an unexpected many-to-many relationship between `x` and `y`.
## ℹ Row 29117 of `x` matches multiple rows in `y`.
## ℹ Row 53969 of `y` matches multiple rows in `x`.
## ℹ If a many-to-many relationship is expected, set `relationship =
##   "many-to-many"` to silence this warning.
```



Table: (\#tab:unnamed-chunk-5)GSE119732 : mapped table preview

|ensembl_gene_id |hgnc_symbol |gene_id           | A1|  A2| A3| A4| B1|
|:---------------|:-----------|:-----------------|--:|---:|--:|--:|--:|
|ENSG00000223972 |DDX11L1     |ENSG00000223972.5 |  0|   0|  0|  0|  0|
|ENSG00000227232 |WASH7P      |ENSG00000227232.5 | 79| 119| 84| 50| 80|
|ENSG00000278267 |MIR6859-1   |ENSG00000278267.1 | 17|  10| 22| 19| 19|
|ENSG00000243485 |MIR1302-2HG |ENSG00000243485.4 |  0|   0|  0|  0|  0|
|ENSG00000237613 |FAM138A     |ENSG00000237613.2 |  0|   0|  0|  0|  0|

## Exercise 2 — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122380

Using **GSE122380**, confirm whether the ID column contains Ensembl IDs with version suffixes.

1. Extract the first 20 IDs.
2. Create a new column with versions stripped.
3. Map the identifiers to HGNC symbols.
4. What is different about this file? 


### Question 1 

``` r
gse <- "GSE122380"
pattern <- "GSE122380.*\\.gz"

path <- file.path("data", gse)
files <- list.files(path, pattern = pattern, 
                    full.names = TRUE, recursive = TRUE)

kable_head(tibble(file = basename(files)), min(10, length(files)), paste(gse,": extracted file list (first 10)"))
```



Table: (\#tab:unnamed-chunk-6)GSE122380 : extracted file list (first 10)

|file                        |
|:---------------------------|
|GSE122380_raw_counts.txt.gz |

``` r
x <- safe_read(files[1])
```

```
## Warning: One or more parsing issues, call `problems()` on your data frame for details,
## e.g.:
##   dat <- vroom(...)
##   problems(dat)
```

```
## Parsing issues detected in TSV — reading as space-delimited file instead.
```

``` r
kable_head(x[, 1:min(6, ncol(x))], 20, paste(gse,": raw table preview"))
```



Table: (\#tab:unnamed-chunk-6)GSE122380 : raw table preview

|Gene_id         | 18489_0| 18489_10| 18489_11| 18489_12| 18489_13|
|:---------------|-------:|--------:|--------:|--------:|--------:|
|ENSG00000000419 |     825|      549|      576|      599|      607|
|ENSG00000000457 |     175|      238|      279|      254|      327|
|ENSG00000000460 |     337|      103|      120|      124|      154|
|ENSG00000000938 |       7|        1|        0|        0|        0|
|ENSG00000000971 |       0|        3|        5|       10|        4|
|ENSG00000001036 |    2251|      860|      819|      732|     1071|
|ENSG00000001084 |     980|      833|      694|      615|      959|
|ENSG00000001167 |     886|      991|      889|     1164|     1000|
|ENSG00000001460 |     239|      215|      142|      171|      225|
|ENSG00000001461 |     402|      798|      664|      777|      907|
|ENSG00000001561 |     336|      293|      207|      195|      309|
|ENSG00000001617 |    1850|      955|      472|     1174|      737|
|ENSG00000001626 |     100|       55|       46|      101|       52|
|ENSG00000001629 |    1266|     1055|     1037|     1255|     1248|
|ENSG00000001630 |     157|      176|      150|      130|      179|
|ENSG00000001631 |     641|      843|      812|      861|      984|
|ENSG00000002016 |     403|      543|      354|      428|      464|
|ENSG00000002330 |     131|      299|      301|      298|      374|
|ENSG00000002549 |     742|      869|      902|      713|      973|
|ENSG00000002587 |      11|       50|       51|       54|       30|


### Question 2

``` r
# Extract the first 20 IDs
first_20_ids <- head(x$gene_id, 20)
```

```
## Warning: Unknown or uninitialised column: `gene_id`.
```

``` r
# Count how many of the first 20 contain a "."
dot_count <- sum(grepl("\\.", first_20_ids))
print(paste("Number of IDs with version suffixes in first 20:", dot_count))
```

```
## [1] "Number of IDs with version suffixes in first 20: 0"
```

``` r
# Count how many total contain a "."
dot_count <- sum(grepl("\\.", x$gene_id))
```

```
## Warning: Unknown or uninitialised column: `gene_id`.
```

``` r
print(paste("Total number of IDs with version suffixes:", dot_count))
```

```
## [1] "Total number of IDs with version suffixes: 0"
```

0 from the 1st 20 have a ".", all ids of this whole set
dpn't have a ".".


### Question 3 

``` r
# Define the ID column name
id_col <- names(x)[1]
ids <- x[[1]] |> as.character()

# Since IDs don't have dots, we use them directly
if (any(grepl("^ENSG", ids))) {
  library(biomaRt)
  library(dplyr)
  
  ensembl_ids <- unique(ids)
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Fetch the mapping
  map <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
  
  kable_head(map, 10, paste(gse, ": mapping preview (cleaned)"))
  
  # Join using the original IDs directly
  expr_mapped <- x %>%
    mutate(ensembl_gene_id = .data[[id_col]]) %>%
    left_join(map, by = "ensembl_gene_id") %>%
    dplyr::select(ensembl_gene_id,hgnc_symbol, everything())
  
  # Preview the final table
  kable_head(expr_mapped[, 1:min(8, ncol(expr_mapped))], 10, paste(gse, ": mapped table preview"))
}
```



Table: (\#tab:unnamed-chunk-8)GSE122380 : mapped table preview

|ensembl_gene_id |hgnc_symbol |Gene_id         | 18489_0| 18489_10| 18489_11| 18489_12| 18489_13|
|:---------------|:-----------|:---------------|-------:|--------:|--------:|--------:|--------:|
|ENSG00000000419 |DPM1        |ENSG00000000419 |     825|      549|      576|      599|      607|
|ENSG00000000457 |SCYL3       |ENSG00000000457 |     175|      238|      279|      254|      327|
|ENSG00000000460 |FIRRM       |ENSG00000000460 |     337|      103|      120|      124|      154|
|ENSG00000000938 |FGR         |ENSG00000000938 |       7|        1|        0|        0|        0|
|ENSG00000000971 |CFH         |ENSG00000000971 |       0|        3|        5|       10|        4|
|ENSG00000001036 |FUCA2       |ENSG00000001036 |    2251|      860|      819|      732|     1071|
|ENSG00000001084 |GCLC        |ENSG00000001084 |     980|      833|      694|      615|      959|
|ENSG00000001167 |NFYA        |ENSG00000001167 |     886|      991|      889|     1164|     1000|
|ENSG00000001460 |STPG1       |ENSG00000001460 |     239|      215|      142|      171|      225|
|ENSG00000001461 |NIPAL3      |ENSG00000001461 |     402|      798|      664|      777|      907|


### Question 4
GSE122380 uses stable Ensembl IDs without version suffixes, whereas your 
previous file used versioned IDs (e.g., .5) that required stripping, 
which made is simpler to map to HGNC.


## Exercise 3 - 
 
Can you use the worked example to process the above two GEO records?  How?

Yes — you can process both GEO records using the exact same workflow but with
minor change in variables and remove the stripping version number step if not
needed.
