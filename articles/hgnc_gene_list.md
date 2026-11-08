# Clean a gene list with HGNC

## The problem

Gene lists arrive from papers, spreadsheets and colleagues with whatever
symbols their authors used: aliases (`HER2`), old symbols (`SEPT1`),
symbols HGNC has withdrawn, and names that are simply not symbols at
all. Before such a list can be joined to DrugBank targets, HPO
phenotypes or anything else, each entry has to be mapped to an approved
HGNC symbol, and the ambiguous ones have to be seen, not guessed.

This article does that with
[`parseHGNC()`](https://docs.ropensci.org/dbparser/reference/parseHGNC.md)
and `dplyr` on the HGNC sample that ships with `dbparser`. It runs
offline in a few seconds.

``` r

library(dbparser)
library(dplyr)
```

## The list

``` r

genes <- c("PD1", "PDL1", "HER2", "SEPT1", "MARCH1", "FOX1",
           "CBBM", "A12M1", "TP53")
```

The bundled sample holds only 24 genes and a few withdrawn entries, so
every input above was chosen to be in it (or, for `FOX1`, to show what a
miss looks like). On HGNC’s full files the same code covers all approved
genes; see [Full data](#full-data) below.

## Parse the HGNC sample

``` r

hgnc_dir <- system.file("extdata", "hgnc", package = "dbparser")
hgnc <- parseHGNC(file.path(hgnc_dir, "hgnc_complete_set.txt"),
                  withdrawn_path = file.path(hgnc_dir, "withdrawn.txt"))

names(hgnc)
#> [1] "gene"      "gene_name" "withdrawn"
nrow(hgnc$gene)
#> [1] 24
```

Every `dvobject` records where its data came from, including the licence
and the citation HGNC asks for:

``` r

glimpse(get_source_metadata(hgnc))
#> Rows: 1
#> Columns: 12
#> $ metadata_version <chr> "1"
#> $ source           <chr> "hgnc"
#> $ source_name      <chr> "HUGO Gene Nomenclature Committee (HGNC)"
#> $ release          <chr> NA
#> $ release_date     <date> NA
#> $ retrieved_at     <date> NA
#> $ source_url       <chr> "https://www.genenames.org/download/statistics-and-fi…
#> $ licence          <chr> "CC0-1.0"
#> $ licence_url      <chr> "https://www.genenames.org/about/license/"
#> $ attribution      <chr> "Seal RL, Braschi B, Gray K, McClay J, Tweedie S, Bru…
#> $ files            <list> [<data.frame[2 x 2]>]
#> $ parser           <chr> "dbparser 2.3.0"
```

HGNC publishes no release number for its current files, so `release`,
`release_date` and `retrieved_at` are `NA` here; pass them to
[`parseHGNC()`](https://docs.ropensci.org/dbparser/reference/parseHGNC.md)
when you parse your own download.

## Match each input

`gene_name` has one row per approved symbol, alias and previous symbol
(and the names, which we do not need here). Each row carries its gene’s
`hgnc_id`, so a match gives the approved symbol directly.

``` r

approved <- hgnc$gene |>
  select(hgnc_id, matched = symbol)

name_hits <- hgnc$gene_name |>
  filter(kind %in% c("symbol", "alias_symbol", "prev_symbol")) |>
  inner_join(approved, by = "hgnc_id") |>
  transmute(input = name,
            matched,
            how = recode(kind,
                         symbol       = "approved",
                         alias_symbol = "alias",
                         prev_symbol  = "previous"))
```

Withdrawn symbols are in their own table, with the gene or genes HGNC
merged or split them into. An entry withdrawn outright points nowhere.

``` r

withdrawn_hits <- hgnc$withdrawn |>
  transmute(input   = withdrawn_symbol,
            matched = merged_into_symbol,
            how     = "withdrawn")
```

Join the list to both, keeping every match: an input that fits several
genes gets several rows here.

``` r

hits <- tibble(input = genes) |>
  left_join(bind_rows(name_hits, withdrawn_hits),
            by = "input", relationship = "one-to-many") |>
  distinct()

hits
#> # A tibble: 14 × 3
#>    input  matched how      
#>    <chr>  <chr>   <chr>    
#>  1 PD1    PDCD1   alias    
#>  2 PD1    SNCA    alias    
#>  3 PD1    SPATA2  alias    
#>  4 PDL1   CD274   alias    
#>  5 HER2   ERBB2   alias    
#>  6 SEPT1  SEPTIN1 previous 
#>  7 MARCH1 MARCHF1 previous 
#>  8 FOX1   NA      NA       
#>  9 CBBM   OPN1LW  previous 
#> 10 CBBM   OPN1MW  previous 
#> 11 CBBM   OPN1MW  withdrawn
#> 12 CBBM   OPN1LW  withdrawn
#> 13 A12M1  NA      withdrawn
#> 14 TP53   TP53    approved
```

## One row per input

Collapse to one row per input. `needs_choice` is `TRUE` when an input
fits more than one approved gene; these are listed, not resolved.

``` r

clean <- hits |>
  group_by(input) |>
  summarise(
    needs_choice = n_distinct(na.omit(matched)) > 1,
    how          = if (all(is.na(how))) "no match"
                   else paste(unique(na.omit(how)), collapse = ", "),
    matched      = paste(sort(unique(na.omit(matched))), collapse = ", "),
    .groups      = "drop") |>
  select(input, matched, how, needs_choice) |>
  arrange(match(input, genes))

knitr::kable(clean)
```

| input  | matched             | how                 | needs_choice |
|:-------|:--------------------|:--------------------|:-------------|
| PD1    | PDCD1, SNCA, SPATA2 | alias               | TRUE         |
| PDL1   | CD274               | alias               | FALSE        |
| HER2   | ERBB2               | alias               | FALSE        |
| SEPT1  | SEPTIN1             | previous            | FALSE        |
| MARCH1 | MARCHF1             | previous            | FALSE        |
| FOX1   |                     | no match            | FALSE        |
| CBBM   | OPN1LW, OPN1MW      | previous, withdrawn | TRUE         |
| A12M1  |                     | withdrawn           | FALSE        |
| TP53   | TP53                | approved            | FALSE        |

What the table shows:

- **`TP53`** is an approved symbol; nothing to do.
- **`HER2`** and **`PDL1`** are aliases of `ERBB2` and `CD274`.
- **`SEPT1`** and **`MARCH1`** are previous symbols: HGNC renamed them
  `SEPTIN1` and `MARCHF1`. The `gene` table has the date of each change.
- **`PD1`** is an alias of three genes (`PDCD1`, `SNCA` and `SPATA2`).
  Most readers mean `PDCD1`, but the list does not say so, and only its
  author can.
- **`CBBM`** is a previous symbol of both `OPN1LW` and `OPN1MW`, and
  HGNC’s withdrawn file records it as split into the same two genes, so
  it needs a choice too.
- **`A12M1`** was withdrawn and points to no gene.
- **`FOX1`** finds nothing. In this sample `RBFOX1` has the alias
  `FOX-1`, with a hyphen; exact matching does not guess that the two are
  the same.

``` r

hgnc$gene |>
  filter(symbol %in% c("SEPTIN1", "MARCHF1")) |>
  select(symbol, prev_symbol, date_symbol_changed)
#> # A tibble: 2 × 3
#>   symbol  prev_symbol date_symbol_changed
#>   <chr>   <chr>       <date>             
#> 1 MARCHF1 MARCH1      2019-07-09         
#> 2 SEPTIN1 DIFF6|SEPT1 2019-04-04
```

## Withdrawn symbols

The withdrawn table keeps HGNC’s own status for the withdrawn ID and for
each gene it points to:

``` r

hgnc$withdrawn |>
  filter(withdrawn_symbol %in% genes) |>
  select(withdrawn_symbol, status, merged_into_symbol, merged_into_status)
#> # A tibble: 3 × 4
#>   withdrawn_symbol status          merged_into_symbol merged_into_status
#>   <chr>            <chr>           <chr>              <chr>             
#> 1 A12M1            Entry Withdrawn NA                 NA                
#> 2 CBBM             Merged/Split    OPN1MW             Approved          
#> 3 CBBM             Merged/Split    OPN1LW             Approved
```

## Next steps

Inputs with `needs_choice = TRUE`, a `withdrawn` match with no target,
or no match at all go back to whoever wrote the list. The rest map
one-to-one to an approved symbol and can be joined to other sources, for
example DrugBank targets or HPO gene-phenotype associations.

## Full data

The sample is enough to show the method, not to clean a real list. The
script
[`announce/hgnc_full.R`](https://github.com/ropensci/dbparser/blob/cran/announce/hgnc_full.R)
in the package repository runs the same steps on your own download of
`hgnc_complete_set.txt` and `withdrawn.txt` from
[genenames.org](https://www.genenames.org/download/statistics-and-files/).
