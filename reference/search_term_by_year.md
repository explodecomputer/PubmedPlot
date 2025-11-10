# Search for a term in the database, split up by year.

This function searches for a given term in the database and returns the
results. It splits the search into years to avoid hitting the NCBI
server limits.

## Usage

``` r
search_term_by_year(term, years)
```

## Arguments

- term:

  The term to search for.

- years:

  The years to search for.

## Value

The search results.
