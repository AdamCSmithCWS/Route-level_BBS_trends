cv_folds <- function(raw_data,
                     k = 10,  # number of folds
                     fold_groups = "route",
                     omit_singles = TRUE,
                     quiet = FALSE){
  
  
  # fold_groups is the critical grouping factor for the leave future out
  # cross-validation. Each fold-k identifies a test-set of observations that
  # include all future observations from a randomly selected set of groups in
  # fold_group (e.g., all future observations for a set of observers)
  # balanced across all routes and strata
  
  # ID training and testing
  # Each training-fold (inverse of those in fold-k) must include all values of
  # fold_groups. This is challenging for fold_groups with only 1-year of data -
  # they can't be used in any of the testing sets
  if(is.null(raw_data)){stop("Data file is missing")
    return(NULL)}
  # Split observers in to K groups, based on fold_groups
  cv <- dplyr::mutate(raw_data,
                      grouping = .data[[fold_groups]])
  
  n_groups <- length(unique(cv$grouping))
  
  # Make a dataframe that assigns each group to one of the K folds
  grps <- cv %>%
    dplyr::distinct(.data$grouping) %>%
    dplyr::mutate(fold = rep(1:k, length.out = n_groups))
  
  cv <- cv %>%
    dplyr::left_join(grps, by = "grouping") %>%
    dplyr::group_by(grouping) %>%
    # removes first-year observations so they are always in training set
    dplyr::mutate(fold = dplyr::if_else(.data$firstyr ==1, NA, .data$fold))
  
  # Optional: Exclude all groups with no replication from testing sets
  # if(omit_singles){
  #   n_events_by_group <- cv %>%
  #     dplyr::group_by(.data$grouping, .data$first_year) %>%
  #     dplyr::summarise(n_events = dplyr::n(), .groups = "drop") %>%
  #     dplyr::group_by(.data$grouping) %>%
  #     dplyr::summarise(n_cats = dplyr::n(),
  #                      only_first_year = dplyr::if_else(.data$n_cats == 1,
  #                                                       TRUE,
  #                                                       FALSE),
  #                      n_events_cv = sum(.data$n_events),
  #                      .groups = "drop")
  #   
  #   cv <- cv %>%
  #     dplyr::left_join(n_events_by_group, by = "grouping") %>%
  #     dplyr::left_join(grps, by = "grouping") %>%
  #     dplyr::group_by(grouping) %>%
  #     dplyr::mutate(
  #       # bumps first-year observations to the previous fold
  #       fold2 = .data$fold - .data$first_year,
  #       # wrap around the list of K values for fold == 1 and first_year == 1
  #       fold = dplyr::if_else(.data$fold2 < 1,
  #                             .env$k,
  #                             .data$fold2),
  #       # remove groups with only first-year observer route combinations
  #       fold = dplyr::if_else(.data$only_first_year, NA_real_, .data$fold))
  #   
  # } else {
    # cv <- cv %>%
    #   dplyr::left_join(grps, by = "grouping") %>%
    #   dplyr::group_by(grouping) %>%
    #   # bumps first-year observations to the next fold
    #   dplyr::mutate(fold2 = .data$fold + .data$first_year,
    #                 fold = dplyr::if_else(.data$fold2 > .env$k, 1, .data$fold2))
 # }
  
  return(dplyr::pull(cv, .data$fold))
}
