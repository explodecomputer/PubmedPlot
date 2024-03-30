#' Search for a term in the database.
#'
#' This function searches for a given term in the database and returns the results.
#'
#' @param term The term to search for.
#' @return The search results.
#' @export
search <- function(term) {
    search_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_params <- list(
        db = "pubmed",
        term = term,
        retmode = "json",
        usehistory = "y",
        retmax = 20000
    )

    search_response <- httr::GET(url = search_url, query = search_params)
    search_content <- httr::content(search_response, "text")
    search_result <- jsonlite::fromJSON(search_content)

    pmids <- search_result$esearchresult$idlist
    # length(pmids)

    summary_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

    # Prepare the body of the POST request
    # This includes specifying db, id, and retmode as before, but formatted for a POST request
    body_list <- list(
        db = "pubmed",
        id = paste(pmids, collapse = ","),
        retmode = "xml"
    )
    body <- paste(names(body_list), body_list, sep = "=", collapse = "&")

    # Make the POST request
    summary_response <- httr::POST(url = summary_url, body = body, encode = "form")
    summary_content <- httr::content(summary_response, "text", encoding = "UTF-8")

    doc <- XML::xmlParse(summary_content)

    pub_dates <- XML::xpathApply(doc, "//DocSum", function(node) {
        pmid <- XML::xmlValue(node[["Id"]])
        pub_date <- XML::xmlValue(node[[2]])
        dplyr::tibble(pmid = pmid, pub_date = pub_date)
    }) %>% dplyr::bind_rows() %>% dplyr::mutate(pub_date = lubridate::ymd(pub_date))

    max_pubmed_date <- max(pub_dates$pub_date, na.rm=TRUE)

    time_levels <- c("weeks","months","years")
    cap_verb="Data from PubMed search for Mendelian randomi[s/z]ation, title only. "
    xlab="Date"
    ylab_root="New PubMed Entries"
    mycol="dodgerblue3"
    lwidth=.4

    res <- lapply(c(time_levels), function(z){  
        ## assign time groups then remove record duplications from e.g. pub comments
        dat <- pub_dates %>% 
            dplyr::mutate(tend = lubridate::ceiling_date(pub_date, unit = z) - 1) %>%
            dplyr::group_by(pmid) %>%
            dplyr::filter(pub_date==min(pub_date)) %>%
            dplyr::ungroup()
        
        ## ensure all date rows exist even as zero.
        agg <- dplyr::left_join(
            tidyr::crossing(
                tend=seq.Date(from=min(dat$tend),to=max(dat$tend),by=paste(1,gsub("s$","",z)))) %>%
                dplyr::mutate(tend=lubridate::ceiling_date(tend,unit=z)-1),
            dat %>%
                dplyr::group_by(tend) %>%
                dplyr::summarise(n = as.numeric(dplyr::n()), .groups='drop'),by = dplyr::join_by(tend)
        ) %>%
            dplyr::mutate(n=dplyr::if_else(!is.finite(n),as.numeric(0),n)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(time_level=z) 
    }) %>%
        dplyr::bind_rows() %>%
        dplyr::select(pubmed_date=tend,n_publications=n,time_level) %>%
        dplyr::filter(pubmed_date < lubridate::today()) %>%
        dplyr::mutate(time_level=stringr::str_to_title(time_level)) %>%
        dplyr::mutate(last_pubmed_date=max_pubmed_date)
    return(list(pub_dates=pub_dates, res=res))
}

