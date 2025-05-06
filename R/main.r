#' Search for a term in the database.
#'
#' This function searches for a given term in the database and returns the results.
#'
#' @param term The term to search for.
#' @return The search results.
#' @export
search_term <- function(term) {
  Sys.sleep(1)
  search_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
  search_params <- list(
    db = "pubmed",
    term = term,
    retmode = "json",
    usehistory = "y",
    retmax = 0  # We only want the count and search history, not the IDs
  )

  search_response <- httr::GET(url = search_url, query = search_params)
  search_content <- httr::content(search_response, "text")
  search_result <- jsonlite::fromJSON(search_content)

  # Get the WebEnv, query_key, and count
  web_env <- search_result$esearchresult$webenv
  query_key <- search_result$esearchresult$querykey
  count <- as.numeric(search_result$esearchresult$count)

  cat("Found", count, "articles matching the search term.\n")

  if(count == 0) {
    return(NULL)
  }

  # Step 2: Use EFetch to download records in batches
  fetch_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  batch_size <- 500  # NCBI recommends no more than 500 records per request
  all_pub_dates <- list()

  for (start in seq(0, count - 1, by = batch_size)) {
    cat("Fetching records", start + 1, "to", min(start + batch_size, count), "\n")
  
    fetch_params <- list(
      db = "pubmed",
      query_key = query_key,
      WebEnv = web_env,
      retstart = start,
      retmax = batch_size,
      retmode = "xml",
      rettype = "abstract"
    )
    
    # Add a delay to avoid overloading NCBI servers (3-10 requests per second limit)
    if (start > 0) {
      Sys.sleep(0.5)  # 500ms delay between requests
    }

    fetch_response <- httr::GET(url = fetch_url, query = fetch_params)
    fetch_content <- httr::content(fetch_response, "text", encoding = "UTF-8")

    if (httr::status_code(fetch_response) == 200) {
      doc <- XML::xmlParse(fetch_content)
      xmltop <- XML::xmlRoot(doc)
      # xmlSize(xmltop)
      # xmlName(xmltop[[1]][[1]][[1]])
      # xmlValue(xmltop[[1]][[]][["PMID"]])

      pub_dates <- XML::xpathApply(doc, '//PubmedArticle', \(x) {
      dplyr::tibble(
        pmid = XML::xmlValue(x[[1]][["PMID"]]),
        ab = XML::xmlValue(x[[1]][["Article"]][["Abstract"]]),
        pub_date = lubridate::ymd(
        paste(
          XML::xmlValue(x[["PubmedData"]][["History"]][["PubMedPubDate"]][["Year"]]),
          XML::xmlValue(x[["PubmedData"]][["History"]][["PubMedPubDate"]][["Month"]]),
          XML::xmlValue(x[["PubmedData"]][["History"]][["PubMedPubDate"]][["Day"]])
        )
        ),
        title = XML::xmlValue(x[[1]][["Article"]][["ArticleTitle"]]),
        journal_issn = XML::xmlValue(x[[1]][["Article"]][["Journal"]][["ISSN"]]),
        journal = XML::xmlValue(x[[1]][["Article"]][["Journal"]][["Title"]]),
        author_affil = XML::xmlValue(x[[1]][["Article"]][["AuthorList"]][[1]][["AffiliationInfo"]])
      )
      }) %>% dplyr::bind_rows()
      all_pub_dates[[start+1]] <- pub_dates
    } else {
      cat("Error fetching records:", httr::status_code(fetch_response), "\n")
      cat("Response content:", fetch_content, "\n")
    }
  }
  return(all_pub_dates)
}


#' Search for a term in the database, split up by year.
#'
#' This function searches for a given term in the database and returns the results.
#' It splits the search into years to avoid hitting the NCBI server limits.
#'
#' @param term The term to search for.
#' @param years The years to search for.
#' @return The search results.
#' @export
search_term_by_year <- function(term, years) {
  pub_dates <- list()
  for (year in years) {
    cat("Searching for year:", year, "\n")
    sterm <- paste0(term, ' AND ("', year, '/01/01"[dp] : "', year, '/12/31"[dp])')
    cat("Search term:", sterm, "\n")
    pub_dates[[as.character(year)]] <- search_term(sterm)
  }
  pub_dates <- dplyr::bind_rows(pub_dates) %>% 
    dplyr::filter(!duplicated(pmid))
  return(pub_dates)
}




#' Organise results into counts by time interval
#'
#' Provide a weekly, monthly or yearly count of the number of publications
#'
#' @param pub_dates Output from search_term or search_term_by_year
#' @return The search results.
#' @export
group_by_time_interval <- function(pub_dates) {
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


# Document this function
#' Plot the results of the search
#'
#' This function plots the results of the search.
#' It uses ggplot2 to create a line plot of the number of publications over time.
#' 
#' @param x The results group_by_time_interval
#' @param interval The time interval to plot. Default is "Weeks". Choose "Weeks", "Months", "Years".
#' @return A ggplot object.
#' @importFrom ggplot2 '%+replace%'
#' @export
plot_time_interval <- function(x, interval = "Weeks") {

  lab <- tolower(interval) %>% gsub("s$", "", .) %>% paste0("by ", .)
  cap_verb="Data from PubMed search for Mendelian randomi[s/z]ation, title only. "
  xlab="Date"
  ylab_root="New PubMed Entries"
  mycol="dodgerblue3"
  lwidth=.4

  '%+replace%' <- ggplot2::'%+replace%' # nolint
  theme_mrlit <- function(...){
    ggplot2::theme_classic(base_size=14) %+replace%
      ggplot2::theme(legend.position="bottom",
            strip.background=ggplot2::element_rect(colour='grey99',fill='grey98'),
            panel.grid.major=ggplot2::element_line(linewidth=lwidth/3,colour='grey'),
            panel.grid.minor=ggplot2::element_line(linewidth=lwidth/5,colour='grey'),
            plot.background=ggplot2::element_rect(fill=NA,linewidth=lwidth/2,colour='grey'),
            panel.border=ggplot2::element_rect(fill=NA,linewidth=lwidth/2,colour='grey'),
            axis.line=ggplot2::element_line(linewidth=.2),
            axis.ticks=ggplot2::element_line(linewidth=.2))
  }

  ## plot ####
  p <- ggplot2::ggplot(x$res |> 
            dplyr::filter(time_level == interval & pubmed_date != max(pubmed_date,na.rm=F)), ggplot2::aes(pubmed_date,n_publications)) +
    ggplot2::geom_point(size = .25, colour = mycol, alpha = 1) +
    ggplot2::geom_line(alpha = .75, colour = mycol) +
    ggplot2::geom_line(data = x$res %>% dplyr::filter(time_level == interval), linetype=2, alpha=.75, colour=mycol) +
    ggplot2::scale_x_date(date_minor_breaks = "1 year") +
    
    ggplot2::labs(x="Date",
        y=paste(ylab_root, lab),
        caption=paste("Most recent publication date:", max(subset(x$res, time_level==interval)$pubmed_date)),
        x="Date") +
    theme_mrlit()

  p

}
