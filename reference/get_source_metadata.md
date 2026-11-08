# Get the source metadata of a dvobject

Returns where each source in a dvobject came from: release, dates, URLs,
licence, attribution, the parsed files with their SHA-256 hashes, and
the dbparser version that parsed them. A merged dvobject has one row per
source. Unknown values are `NA`.

## Usage

``` r
get_source_metadata(dvobject)
```

## Arguments

- dvobject:

  a dvobject returned by a parser or a merge function

## Value

a tibble with one row per source and columns `metadata_version`,
`source`, `source_name`, `release`, `release_date`, `retrieved_at`,
`source_url`, `licence`, `licence_url`, `attribution`, `files` (a list
column of data.frames with `name` and `sha256`) and `parser`. It has no
rows when the dvobject carries no source metadata (e.g. one saved by an
older dbparser version).

## See also

Other utility:
[`add_database_info()`](https://docs.ropensci.org/dbparser/reference/add_database_info.md),
[`get_genes_for_phenotype()`](https://docs.ropensci.org/dbparser/reference/get_genes_for_phenotype.md),
[`search_phenotypes()`](https://docs.ropensci.org/dbparser/reference/search_phenotypes.md),
[`show_dvobject_metadata()`](https://docs.ropensci.org/dbparser/reference/show_dvobject_metadata.md),
[`subset_drugbank_dvobject()`](https://docs.ropensci.org/dbparser/reference/subset_drugbank_dvobject.md),
[`subset_onsides_dvobject()`](https://docs.ropensci.org/dbparser/reference/subset_onsides_dvobject.md)

## Examples

``` r
drugbank <- parseDrugBank(
  system.file("extdata", "drugbank_record_biotech.xml",
              package = "dbparser"))
#> Loading DrugBank DB from path: /github/home/R/x86_64-pc-linux-gnu-library/4.6/dbparser/extdata/drugbank_record_biotech.xml
#> Completed loading DrugBank DB into memory
#> ...........................................
#> parsing drugs elements
#> drugs information parsing has started
#> parsing drugs general information
get_source_metadata(drugbank)
#> # A tibble: 1 × 12
#>   metadata_version source   source_name release release_date retrieved_at
#>   <chr>            <chr>    <chr>       <chr>   <date>       <date>      
#> 1 1                drugbank DrugBank    5.1     2018-07-03   NA          
#> # ℹ 6 more variables: source_url <chr>, licence <chr>, licence_url <chr>,
#> #   attribution <chr>, files <list>, parser <chr>
```
