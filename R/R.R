# sensiverse: all core code in one script
# =============================================================================
# Handbook order:
# 1) Model Space Estimation (OLS universe)
# 2) Filtering
# 3) Sign shares + plotting
# 4) Uncertainty sources: Multinomial Logit (MLOGIT)
# 5) Uncertainty sources: Nearest Neighbour (1NN)
# 6) Uncertainty sources: Neural Network (Keras + SHAP)
# =============================================================================

# --- Package-level imports ----------------------------------------------------
#' @title sensiverse
#' @description
#' Tools to estimate, filter, and interpret large model spaces for specification
#' sensitivity. The package provides (i) fast OLS grid estimation, (ii) filters
#' for curated universes, (iii) sign-share diagnostics with plots, and (iv) three
#' complementary methods to pinpoint uncertainty sources (MLOGIT, 1NN,
#' and NN+SHAP).
#'
#' @keywords internal
#' @importFrom stats AIC logLik predict reorder lm
#' @importFrom utils globalVariables
#' @importFrom lmtest coeftest
#' @importFrom sandwich vcovHC
#' @importFrom multiwayvcov cluster.vcov
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr rename mutate filter select case_when group_by summarise
#' @importFrom dplyr bind_rows left_join arrange
#' @importFrom purrr list_rbind map
#' @importFrom ggplot2 ggplot aes geom_col coord_flip theme_minimal labs
#' @importFrom ggplot2 scale_fill_manual element_blank element_text
"_PACKAGE"

# Silence CMD check for NSE
utils::globalVariables(c(
  ".", ".data", "IV", "MODELID", "coef", "se", "pval", "set", "sig_sign",
  "contset", "fes", "dep", "sample", "r2", "aic", "loglh", "nobs",
  "tval", "dimension"
))


# =============================================================================
# Helpers (internal)
# =============================================================================

#' Map fixed-effects label to a formula fragment
#' @param fes_label Character scalar. One of "none","time","group","time+group".
#' @return Character string (RHS fragment) or error on unknown label.
#' @keywords internal
.fes_to_formula <- function(fes_label) {
  switch(
    fes_label,
    "none"       = "",
    "time"       = "factor(FEtime)",
    "group"      = "factor(FEgroup)",
    "time+group" = "factor(FEgroup)+factor(FEtime)",
    stop("Unknown fixed-effect label: ", fes_label)
  )
}

#' Build an lm() formula from components
#' @param dep Character scalar (dependent variable).
#' @param contset Character scalar like "CONTa+CONTb".
#' @param fes_label Character scalar FE label (see \code{.fes_to_formula}).
#' @return A formula object \code{dep ~ ...}.
#' @keywords internal
.build_formula <- function(dep, contset, fes_label) {
  rhs <- c(contset, .fes_to_formula(fes_label))
  rhs <- rhs[rhs != ""]
  as.formula(paste(dep, "~", paste(rhs, collapse = " + ")))
}


# =============================================================================
# 1) Model Space Estimation (OLS universe)
# =============================================================================

#' Estimate a specification universe (OLS) for a focus variable
#'
#' @description
#' Runs OLS regressions across a user-defined grid (dependent variables,
#' control sets, fixed effects, SE type, and sample split) and collects results
#' only for \code{focus_var}.
#'
#' @param df Data frame containing \code{DEP*}, \code{CONT*}, \code{FEgroup},
#'   \code{FEtime}, and \code{SAMPLE*} columns.
#' @param focus_var Character. IV of interest (e.g., \code{"CONTa"}).
#' @param specs Named list with components:
#' \describe{
#'   \item{dep}{Character vector of dependent variables present in \code{df}.}
#'   \item{cont}{Character vector of candidate controls (e.g., \code{CONT*}).}
#'   \item{fes}{Subset of \code{c("none","time","group","time+group")}.}
#'   \item{set}{Subset of \code{c("simple","robust","clustered")}.}
#'   \item{sample}{Character vector of sample-split dummies (cols in \code{df}).}
#' }
#' @param space_n Optional integer. If the full grid is larger, a random sample
#'   of at most \code{space_n} specifications is estimated (for speed).
#'
#' @return A tibble with coefficient, SE, p-value, significance class, and
#'   full specification metadata for \code{focus_var}.
#' @export
#' @family model-space
#' @examples
#' out <- estimate_model_space(df, "CONTa", specs)
#' }
estimate_model_space <- function(df, focus_var, specs, space_n = 150000) {
  # (code unchanged)

  # Expand grid manually from specs
  contsets <- unlist(lapply(1:length(specs$cont), function(k) {
    cmb <- utils::combn(specs$cont, k, simplify = FALSE)
    cmb <- Filter(function(v) focus_var %in% v, cmb)
    vapply(cmb, function(v) paste(v, collapse = "+"), character(1))
  }), use.names = FALSE)

  grid <- expand.grid(
    dep     = specs$dep,
    contset = contsets,
    fes     = specs$fes,
    set     = specs$set,
    sample  = specs$sample,
    stringsAsFactors = FALSE
  )

  if (!is.null(space_n) && nrow(grid) > space_n) {
    grid <- dplyr::slice_sample(grid, n = space_n)
  }
  grid$MODELID <- sprintf("%015d", seq_len(nrow(grid)))

  # Internal helper
  .est_one <- function(spec_row) {
    data_temp <- subset(df, base::get(spec_row$sample) == 1)
    form_temp <- .build_formula(spec_row$dep, spec_row$contset, spec_row$fes)
    reg <- stats::lm(formula = form_temp, data = data_temp)

    coef_tbl <- function(vcov_mat, set_label) {
      out <- lmtest::coeftest(reg, vcov = vcov_mat)
      class(out) <- "matrix"
      as.data.frame(out) |>
        tibble::rownames_to_column("IV") |>
        dplyr::rename(coef = 2, se = 3, tval = 4, pval = 5) |>
        dplyr::filter(.data$IV == focus_var) |>
        dplyr::mutate(set = set_label)
    }

    est <- list()
    if ("simple" %in% spec_row$set)
      est[["simple"]] <- coef_tbl(sandwich::vcovHC(reg, type = "const"), "simple")
    if ("robust" %in% spec_row$set)
      est[["robust"]] <- coef_tbl(sandwich::vcovHC(reg, type = "HC1"), "robust")
    if ("clustered" %in% spec_row$set && "FEgroup" %in% names(df))
      est[["clustered"]] <- coef_tbl(multiwayvcov::cluster.vcov(reg, data_temp$FEgroup), "clustered")

    est <- purrr::list_rbind(est)
    if (nrow(est)) {
      est <- est |>
        dplyr::mutate(
          sig_sign = dplyr::case_when(
            pval < 0.05 & coef > 0 ~ "significant positive",
            pval < 0.05 & coef < 0 ~ "significant negative",
            TRUE                   ~ "not significant"
          ),
          MODELID = spec_row$MODELID,
          contset = spec_row$contset,
          fes     = spec_row$fes,
          dep     = spec_row$dep,
          sample  = spec_row$sample
        )
    }

    meta <- tibble::tibble(
      MODELID = spec_row$MODELID,
      dep     = spec_row$dep,
      contset = spec_row$contset,
      fes     = spec_row$fes,
      sample  = spec_row$sample,
      r2      = summary(reg)[["r.squared"]],
      aic     = AIC(reg),
      loglh   = as.numeric(logLik(reg)),
      nobs    = length(reg$residuals)
    )

    list(est = est, meta = meta)
  }

  # Run grid
  res_list <- lapply(seq_len(nrow(grid)), function(i) .est_one(grid[i, ]))
  est_all  <- purrr::map(res_list, "est")  |> purrr::list_rbind()
  meta_all <- purrr::map(res_list, "meta") |> purrr::list_rbind()

  dplyr::left_join(est_all, meta_all,
                   by = c("MODELID","dep","contset","fes","sample"))
}


# =============================================================================
# 2) Filtering
# =============================================================================

#' Filter a model-space result
#'
#' @description
#' Post-process the output of \code{estimate_model_space()} to retain
#' a subset based on AIC quantiles, forbidden pairs, or an arbitrary logical
#' expression evaluated on the result.
#'
#' @param res Tibble from \code{estimate_model_space()}.
#' @param top_aic_pct Numeric in [0,1]. Keep the best (lowest) \code{top_aic_pct}
#'   share by AIC. If \code{NULL}, do not filter by AIC.
#' @param drop_expr Character string, an expression evaluated in the context of
#'   \code{res} (e.g., \code{"dep == 'DEP1' & grepl('CONTb', contset)"}).
#'
#' @return Filtered tibble (same schema as \code{res}).
#' @export
#' @family model-space
#' @examples
filter_model_space <- function(res, top_aic_pct = NULL, drop_expr = NULL) {
  # (code unchanged)
  keep_id <- rep(TRUE, nrow(res))
  if (!is.null(top_aic_pct)) {
    thr <- stats::quantile(res$aic, probs = top_aic_pct, na.rm = TRUE)
    keep_id <- keep_id & (res$aic <= thr)
  }
  if (!is.null(drop_expr)) {
    to_drop <- with(res, eval(parse(text = drop_expr)))
    keep_id <- keep_id & !to_drop
  }
  res[keep_id, , drop = FALSE]
}


# =============================================================================
# 3) Sign shares + plotting
# =============================================================================

#' Compute sign shares for the focus variable
#'
#' @description
#' Aggregates the distribution of \code{sig_sign} across all specs or by a
#' given dimension (e.g., \code{"fes"}, \code{"set"}, \code{"sample"}, \code{"dep"}).
#'
#' @param est_tbl Tibble from \code{estimate_model_space()} (or a filtered subset).
#' @param dimension Optional character. If provided, compute shares by this
#'   column; if \code{NULL}, returns overall shares.
#'
#' @return Tibble with counts, total, and signed percentages (\code{perc}).
#'   Negative values reflect negative significant shares.
#' @export
#' @family model-space
#' @examples
#' @export
calculate_sign_shares <- function(est_tbl, dimension = NULL) {
  # (code unchanged)
  if (is.null(est_tbl) || !nrow(est_tbl)) return(tibble::tibble())
  est_tbl <- est_tbl |> dplyr::filter(!is.na(IV) & nzchar(IV))
  if (is.null(dimension)) {
    tmp <- est_tbl |>
      dplyr::group_by(sig_sign) |>
      dplyr::summarise(freq = dplyr::n(), .groups = "drop")
    tmp$freq_total <- sum(tmp$freq)
  } else {
    tmp <- est_tbl |>
      dplyr::group_by(!!rlang::sym(dimension), sig_sign) |>
      dplyr::summarise(freq = dplyr::n(), .groups = "drop") |>
      dplyr::group_by(!!rlang::sym(dimension)) |>
      dplyr::mutate(freq_total = sum(freq)) |>
      dplyr::ungroup()
    names(tmp)[1] <- "dimension"
  }
  tmp |>
    dplyr::mutate(share = freq / freq_total,
                  perc  = dplyr::case_when(
                    sig_sign == "significant negative" ~ -share,
                    sig_sign == "significant positive" ~  share,
                    TRUE ~ 0))
}

#' Plot sign shares of the model space
#'
#' @description
#' Visualizes signed shares from \code{calculate_sign_shares()} either overall
#' or by dimension
#'
#' @param shares Output of \code{calculate_sign_shares()}.
#' @param dimension Character or \code{NULL}. If grouped by dimension, pass the
#'   same value used in \code{calculate_sign_shares()} for consistent axis labels.
#'
#' @return A \pkg{ggplot2} object.
#' @export
#' @family model-space
#' @examples
plot_sign_share <- function(shares, dimension = NULL) {
  shares <- as.data.frame(shares) %>% subset(., sig_sign != 'not significant')

  # Prettier axis labels
  pretty_labels <- c(
    "set"     = "Standard Error Type",
    "fes"     = "Fixed Effects",
    "dep"     = "Dependent Variable",
    "sample"  = "Sample Split",
    "contset" = "Control Set"
  )

  # Ensure there's ALWAYS a column named `dimension` for the x-axis
  if (is.null(dimension)) {
    shares$dimension <- factor("All")
    xlab <- ""
  } else {
    xlab <- ifelse(dimension %in% names(pretty_labels),
                   pretty_labels[[dimension]],
                   tools::toTitleCase(dimension))
  }

  # Order legend nicely
  if ("sig_sign" %in% names(shares)) {
    shares$sig_sign <- factor(
      shares$sig_sign,
      levels = c("significant negative", "significant positive")
    )
  }

  # ---- order by positive share only ----
  if ("perc" %in% names(shares) && "sig_sign" %in% names(shares)) {
    pos_order <- shares |>
      dplyr::filter(sig_sign == "significant positive") |>
      dplyr::group_by(dimension) |>
      dplyr::summarise(pos_value = sum(perc, na.rm = TRUE), .groups = "drop")

    shares <- dplyr::left_join(shares, pos_order, by = "dimension") |>
      dplyr::mutate(
        pos_value = dplyr::if_else(is.na(pos_value), 0, pos_value),
        dimension = reorder(dimension, pos_value, decreasing = TRUE)
      )
  }

  ggplot2::ggplot(shares, ggplot2::aes(x = dimension, y = perc, fill = sig_sign)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = xlab, y = "Share of significant estimates") +
    ggplot2::scale_fill_manual(
      values = c(
        "significant positive" = "darkgreen",
        "significant negative" = "firebrick"
      ),
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(breaks = seq(-1, 1, 0.25), limits = c(-1, 1)) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
    )
}


# =============================================================================
# 4) Uncertainty sources — MLOGIT source method
# =============================================================================

#' Build the MLOGIT source dataset from model-space results
#'
#' @keywords internal
#' @description
#' Converts \code{res} into a design matrix with factors for \code{dep}, \code{fes},
#' \code{set}, \code{sample} and dummies for every \code{CONT*} ever used
#' in \code{contset}. The outcome is the three-class factor \code{sig_sign}.
#'
#' @param res Tibble from \code{estimate_model_space()} (or filtered).
#' @param focus Character. The focus variable (matches \code{IV} column).
#' @return Data frame with predictors and \code{outcome} factor.
#' @export
#' @family sources
#' @examples
make_source_dataset <- function(res, focus) {
  # (code unchanged)
  df <- res[res$IV == focus, , drop = FALSE]
  if (!nrow(df)) stop("No rows for focus variable in `res`.")

  all_controls <- unique(unlist(strsplit(df$contset, "\\+")))
  all_controls <- sort(all_controls[!is.na(all_controls) & nzchar(all_controls) & grepl("^CONT", all_controls)])

  if (length(all_controls)) {
    ctrl_mat <- sapply(
      all_controls,
      function(ct) as.integer(grepl(paste0("(^|\\+)", ct, "($|\\+)"), df$contset)),
      simplify = TRUE
    )
    ctrl_df <- as.data.frame(ctrl_mat, check.names = FALSE)
  } else ctrl_df <- data.frame()

  spec_df <- data.frame(
    dep    = factor(df$dep),
    fes    = factor(df$fes, levels = c("none","group","time","time+group")),
    set    = factor(df$set, levels = c("simple","robust","clustered")),
    sample = factor(df$sample),
    stringsAsFactors = FALSE
  )

  out <- factor(
    df$sig_sign,
    levels = c("not significant","significant negative","significant positive")
  )

  out_df <- cbind(spec_df, ctrl_df)
  out_df$outcome <- droplevels(out)

  keep <- vapply(out_df, function(x) length(unique(x)) > 1, logical(1))
  keep["outcome"] <- TRUE
  out_df <- out_df[, keep, drop = FALSE]
  rownames(out_df) <- NULL
  out_df
}

#' Safe probability predictions for multinom objects
#'
#' @description
#' A defensive wrapper around \code{predict(type = "probs")} that guarantees a
#' named vector across possible shapes returned by \pkg{nnet}.
#'
#' @param model A fitted \code{nnet::multinom} object.
#' @param newdata Single-row data frame of predictors.
#' @param outcome_levels Character vector of all outcome level names.
#' @return Named numeric vector of probabilities for \code{outcome_levels}.
#' @keywords internal
safe_predict_probs <- function(model, newdata, outcome_levels) {
  # (code unchanged)
  raw <- predict(model, newdata = newdata, type = "probs")
  if (is.vector(raw) && !is.matrix(raw)) {
    vec <- as.numeric(raw)
    names(vec) <- outcome_levels[seq_along(vec)]
    out <- setNames(rep(0, length(outcome_levels)), outcome_levels)
    out[names(vec)] <- vec
    return(out)
  }
  if (is.matrix(raw)) {
    vec <- as.numeric(raw[1, ])
    nms <- colnames(raw)
    if (is.null(nms)) nms <- outcome_levels[seq_along(vec)]
    out <- setNames(rep(0, length(outcome_levels)), outcome_levels)
    out[nms] <- vec
    return(out)
  }
  stop("Unexpected predict() output type.")
}

#' Partial-probability changes by specification feature
#'
#' @description
#' Computes the change in predicted probabilities for each outcome class when
#' toggling one feature at a time away from a baseline spec.
#'
#' @param model A fitted \code{nnet::multinom} object.
#' @param df A dataset produced by \code{make_source_dataset()}.
#' @param ref_outcome Character. Outcome level to reference (must exist in \code{df$outcome}).
#' @return Data frame of probability levels and differences by feature/value.
#' @keywords internal
pp_change_by_feature <- function(model, df, ref_outcome = "not significant") {
  # (code unchanged)
  outcome_levels <- levels(df$outcome)
  if (!(ref_outcome %in% outcome_levels)) stop("Reference outcome not found.")
  base <- df[1, , drop = FALSE]
  for (f in c("dep","fes","set","sample")) {
    if (f %in% names(df)) base[[f]] <- factor(levels(df[[f]])[1], levels = levels(df[[f]]))
  }
  cont_cols <- grep("^CONT", names(df), value = TRUE)
  for (cc in cont_cols) base[[cc]] <- 0L
  pred_cols <- setdiff(names(df), "outcome")
  p0 <- safe_predict_probs(model, base[, pred_cols, drop = FALSE], outcome_levels)
  res_list <- list()
  # vary factors
  for (f in c("dep","fes","set","sample")) {
    if (f %in% names(df)) {
      for (lv in levels(df[[f]])) {
        tmp <- base; tmp[[f]] <- factor(lv, levels = levels(df[[f]]))
        p1 <- safe_predict_probs(model, tmp[, pred_cols, drop = FALSE], outcome_levels)
        for (outc in outcome_levels) {
          res_list[[length(res_list)+1]] <- data.frame(
            SPEC=f, val=lv, outcome_level=outc,
            prop=p1[outc], prop_baseline=p0[outc],
            prop_diff=p1[outc]-p0[outc]
          )
        }
      }
    }
  }
  # vary CONT
  for (cc in cont_cols) {
    tmp <- base; tmp[[cc]] <- 1L
    p1 <- safe_predict_probs(model, tmp[, pred_cols, drop = FALSE], outcome_levels)
    for (outc in outcome_levels) {
      res_list[[length(res_list)+1]] <- data.frame(
        SPEC=cc, val="1", outcome_level=outc,
        prop=p1[outc], prop_baseline=p0[outc],
        prop_diff=p1[outc]-p0[outc]
      )
    }
  }
  do.call(rbind, res_list)
}

#' Identify uncertainty sources via mlogit
#'
#' @description
#' Trains \code{nnet::multinom} on the MLOGIT source dataset and evaluates
#' how changing each spec feature shifts predicted probabilities across the
#' three outcome classes.
#'
#' @param res Tibble from \code{estimate_model_space()}.
#' @param focus Character. Focus variable (matches \code{IV}).
#'
#' @return A list with elements:
#' \describe{
#'   \item{method}{\code{"MLOGIT"}}\item{focus}{Focus variable}
#'   \item{model}{Fitted \code{nnet::multinom}} \item{pp_change}{Data frame of probability changes}
#' }
#' @export
#' @family estimate-sources
find_uncertainty_source_mlogit <- function(res, focus) {
  # (code unchanged)
  dat <- make_source_dataset(res, focus) |> na.omit()
  if (!nrow(dat)) stop("No usable rows.")
  model <- nnet::multinom(outcome ~ ., data = dat, trace = FALSE)
  pp_out <- pp_change_by_feature(model, dat)
  list(method="MLOGIT", focus=focus, model=model, pp_change=pp_out)
}

#' Plot feature importance for mlogit approach
#'
#' @description
#' Visual summary of average absolute probability changes by dimension or by
#' specific feature/value when using the MLOGIT method.
#'
#' @param out Output list from \code{find_uncertainty_source_mlogit()}.
#' @param aggregate Logical. If \code{TRUE}, aggregates to dimensions.
#'
#' @return A \pkg{ggplot2} object.
#' @export
#' @family plot-sources
plot_importance_mlogit <- function(out, aggregate=TRUE) {
  # (code unchanged)
  df <- out$pp_change
  df$importance <- abs(df$prop_diff)

  # map dimensions
  df$dimension <- dplyr::case_when(
    df$SPEC == "dep"    ~ "Dependent variable",
    df$SPEC == "fes"    ~ "Fixed effects",
    df$SPEC == "set"    ~ "SE type",
    df$SPEC == "sample" ~ "Sample split",
    grepl("^CONT", df$SPEC) ~ "Control set",
    TRUE ~ df$SPEC
  )

  if (aggregate) {
    agg <- df |>
      dplyr::group_by(dimension, outcome_level) |>
      dplyr::summarise(importance = mean(importance, na.rm = TRUE), .groups = "drop")

    ggplot2::ggplot(agg, ggplot2::aes(x = reorder(dimension, importance),
                                      y = importance,
                                      color = outcome_level,
                                      shape = outcome_level)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Specification dimension", y = "Avg. change in predicted prob.") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom",
                     legend.title = ggplot2::element_blank())

  } else {
    # Custom labels for CONT variables
    df$spec_label <- ifelse(grepl("^CONT", df$SPEC),
                            paste0("cont: ", df$SPEC),
                            paste(df$SPEC, df$val, sep = ": "))

    ggplot2::ggplot(df, ggplot2::aes(x = reorder(spec_label, importance),
                                     y = importance,
                                     color = outcome_level)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Spec decision", y = "change in predicted prob.") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom",
                     legend.title = ggplot2::element_blank())
  }
}


# =============================================================================
# 5) Uncertainty sources — Nearest Neighbour (1NN)
# =============================================================================

#' Identify uncertainty sources via nearest neighbour approach
#'
#' @description
#' Monte-Carlo heuristic: sample a baseline specification, flip exactly one
#' dimension (dep/fes/set/sample/cont) to its nearest neighbour, re-estimate,
#' and record whether the significance class changes for \code{focus}.
#'
#' @param specs The same \code{specs} list used in \code{estimate_model_space()}.
#' @param df Original data used for estimation.
#' @param focus Character. Focus variable.
#' @param n_draws Integer. Number of perturbation draws.
#' @param seed Random seed.
#' @param max_iter Safety cap on attempts to complete \code{n_draws}.
#'
#' @return Data frame with one row per draw and \code{flip} indicator.
#' @export
#' @family estimate-sources
find_uncertainty_source_neigh <- function(specs, df, focus, n_draws = 100, seed = 123, max_iter = 10 * n_draws) {
  set.seed(seed)
  res_list <- vector("list", n_draws)

  # ----- helper: does a sample have any rows? -----
  has_rows <- function(sample_var) {
    sample_var %in% names(df) && any(df[[sample_var]] == 1, na.rm = TRUE)
  }
  valid_samples <- specs$sample[vapply(specs$sample, has_rows, logical(1))]
  if (length(valid_samples) == 0) stop("All sample dummies are empty (no rows == 1).")

  # ----- estimation helper (same as yours, shortened here) -----
  run_single_estimation <- function(dep, contset, fes, set, sample_var, changeSpec = "baseline") {
    data_temp <- df[df[[sample_var]] == 1, , drop = FALSE]
    if (!nrow(data_temp)) return(NULL)
    fes_term <- switch(fes, "none"=NULL, "time"="FEtime", "group"="FEgroup", "time+group"="FEtime + FEgroup")
    form <- paste(dep, "~", contset, if (!is.null(fes_term)) paste("+", fes_term) else "")
    reg <- tryCatch(lm(as.formula(form), data = data_temp), error = function(e) NULL)
    if (is.null(reg)) return(NULL)
    temp <- tryCatch({
      if (set == "robust" || grepl("robust", changeSpec)) {
        lmtest::coeftest(reg, vcov = sandwich::vcovHC(reg, "HC1"))
      } else if (set == "clustered" || grepl("clustered", changeSpec)) {
        lmtest::coeftest(reg, vcov = multiwayvcov::cluster.vcov(reg, data_temp$FEgroup))
      } else {
        lmtest::coeftest(reg, vcov = sandwich::vcovHC(reg, "const"))
      }
    }, error = function(e) NULL)
    if (is.null(temp)) return(NULL)

    ct <- as.data.frame(unclass(temp))
    ct <- tibble::rownames_to_column(ct, "IV")
    nm <- names(ct)
    # normalize column names
    names(ct)[match(c("Estimate","Std. Error","t value","Pr(>|t|)"), nm, nomatch = 0)] <- c("Estimate","Std.Error","t.value","p.value")
    if (!all(c("Estimate","Std.Error","t.value","p.value") %in% names(ct))) return(NULL)

    est <- ct |>
      dplyr::filter(IV == focus) |>
      dplyr::mutate(
        outcome = dplyr::case_when(
          Estimate >= 0 & p.value <= 0.1 ~ "significant positive",
          Estimate <  0 & p.value <= 0.1 ~ "significant negative",
          TRUE ~ "not significant"
        )
      )
    if (!nrow(est)) return(NULL)
    est
  }

  draws_done <- 0
  attempts <- 0
  while (draws_done < n_draws && attempts < max_iter) {
    attempts <- attempts + 1

    # Baseline random spec
    dep        <- sample(specs$dep, 1)
    other_cont <- setdiff(specs$cont, focus)
    if (length(other_cont) == 0) next
    extra_cont <- sample(other_cont, 1)
    contset    <- paste(c(focus, extra_cont), collapse = "+")
    fes        <- sample(specs$fes, 1)
    set        <- sample(specs$set, 1)
    sample_var <- sample(valid_samples, 1)

    base <- run_single_estimation(dep, contset, fes, set, sample_var, "baseline")
    if (is.null(base)) next
    out_base <- base$outcome[1]

    # Pick a dimension with at least one alternative
    dims <- c("dep","fes","set","sample","CONT")
    dims <- dims[
      c(length(setdiff(specs$dep, dep)) > 0,
        length(setdiff(specs$fes, fes)) > 0,
        length(setdiff(specs$set, set)) > 0,
        length(setdiff(valid_samples, sample_var)) > 0,
        length(specs$cont) > 0)
    ]
    if (length(dims) == 0) next
    dim_type <- sample(dims, 1)

    change_spec <- NULL
    neigh <- NULL
    if (dim_type == "dep") {
      change_spec <- sample(setdiff(specs$dep, dep), 1)
      neigh <- run_single_estimation(change_spec, contset, fes, set, sample_var, change_spec)
    } else if (dim_type == "fes") {
      change_spec <- sample(setdiff(specs$fes, fes), 1)
      neigh <- run_single_estimation(dep, contset, change_spec, set, sample_var, change_spec)
    } else if (dim_type == "set") {
      change_spec <- sample(setdiff(specs$set, set), 1)
      neigh <- run_single_estimation(dep, contset, fes, change_spec, sample_var, change_spec)
    } else if (dim_type == "sample") {
      change_spec <- sample(setdiff(valid_samples, sample_var), 1)
      neigh <- run_single_estimation(dep, contset, fes, set, change_spec, change_spec)
    } else if (dim_type == "CONT") {
      change_spec <- sample(specs$cont, 1)
      contset_alt <- if (grepl(change_spec, contset)) {
        gsub(paste0("(^|\\+)", change_spec, "($|\\+)"), "", contset)
      } else {
        paste0(contset, "+", change_spec)
      }
      contset_alt <- gsub("\\+\\+", "+", contset_alt)
      contset_alt <- gsub("^\\+|\\+$", "", contset_alt)
      if (contset_alt == "") contset_alt <- focus
      neigh <- run_single_estimation(dep, contset_alt, fes, set, sample_var, change_spec)
    }

    if (is.null(neigh)) next
    out_alt <- neigh$outcome[1]
    flip <- as.integer(out_base != out_alt)

    res_list[[draws_done + 1]] <- data.frame(
      dimension      = dim_type,
      val            = change_spec,
      baseline_class = out_base,
      altered_class  = out_alt,
      flip           = flip
    )
    draws_done <- draws_done + 1
  }

  if (draws_done < n_draws) {
    warning(sprintf("Only %d/%d draws completed (insufficient valid alternatives).", draws_done, n_draws))
  }
  dplyr::bind_rows(res_list)
}

#' Plot flip rates by dimension using nearest-neighbour approach
#'
#' @description
#' Shows how often a one-step spec perturbation flips the significance class of
#' the focus variable—summarized by dimension.
#'
#' @param out Data frame from \code{find_uncertainty_source_neigh()}.
#'
#' @return A \pkg{ggplot2} object.
#' @export
#' @family plot-sources
plot_importance_neigh <- function(out) {
  # (code unchanged)
  agg <- out %>%
    dplyr::group_by(dimension) %>%
    dplyr::summarise(flip_rate = mean(flip), .groups = "drop") %>%
    dplyr::ungroup()

  # Pretty labels
  agg$dimension <- dplyr::recode(agg$dimension,
                                 "dep"    = "Dependent variable",
                                 "fes"    = "Fixed effects",
                                 "set"    = "SE type",
                                 "sample" = "Sample split",
                                 "CONT"   = "Control set"
  )

  ggplot2::ggplot(agg, ggplot2::aes(x = reorder(dimension, flip_rate), y = flip_rate)) +
    ggplot2::geom_point(size = 3, color = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::labs(x = "Specification Dimension", y = "Flip Rate (%)") +
    ggplot2::theme_minimal()
}


# =============================================================================
# 6) Uncertainty sources — Neural Network (Keras + SHAP)
# =============================================================================

#' Identify uncertainty sources via neural network approach
#'
#' @description
#' Trains separate one-vs-rest neural networks for each outcome class
#' (\code{"not significant"}, \code{"significant negative"}, \code{"significant positive"})
#' and computes SHAP importances (via \pkg{vip}), and
#' aggregates to both features and dimensions.
#'
#' @param res Tibble from \code{estimate_model_space()}.
#' @param focus Character. Focus variable.
#' @param grid List of hyperparameters with elements \code{dropout}, \code{layers},
#'   \code{units}, \code{learning_rate}. A small grid is recommended.
#' @param split Named numeric vector with \code{train}, \code{val}, \code{test} shares.
#' @param epochs Integer training epochs.
#' @param batch_size Mini-batch size.
#' @param seed Random seed.
#' @param verbose Keras verbosity (0/1/2).
#'
#' @return List with method label, focus, and per-class SHAP summaries:
#' \describe{
#'   \item{all_signs}{Named list (per outcome) with \code{shap_feature} and \code{shap_dimension}.}
#' }
#' @family estimate-sources
find_uncertainty_source_neuronet <- function(
    res, focus,
    grid = list(dropout = c(0, 0.3),
                layers = c(3, 5),
                units  = c(32, 64),
                learning_rate = c(1e-3)),
    split = c(train = 0.6, val = 0.2, test = 0.2),
    epochs = 15, batch_size = 256, seed = 123, verbose = 0
) {
  # (code unchanged)
  set.seed(seed)
  signs <- c("not significant","significant negative","significant positive")
  results <- list()

  # --- helper: build dataset
  make_source_dataset <- function(res, focus) {
    df <- res[res$IV == focus, , drop = FALSE]
    if (!nrow(df)) stop("No rows for focus variable in `res`.")
    all_controls <- unique(unlist(strsplit(df$contset, "\\+")))
    all_controls <- sort(all_controls[!is.na(all_controls) & nzchar(all_controls) & grepl("^CONT", all_controls)])
    if (length(all_controls)) {
      ctrl_mat <- sapply(
        all_controls,
        function(ct) as.integer(grepl(paste0("(^|\\+)", ct, "($|\\+)"), df$contset)),
        simplify = TRUE
      )
      ctrl_df <- as.data.frame(ctrl_mat, check.names = FALSE)
    } else ctrl_df <- data.frame()
    spec_df <- data.frame(
      dep    = factor(df$dep),
      fes    = factor(df$fes, levels = c("none","time","group","time+group")),
      set    = factor(df$set, levels = c("simple","robust","clustered")),
      sample = factor(df$sample),
      stringsAsFactors = FALSE
    )
    out <- factor(
      df$sig_sign,
      levels = c("not significant","significant negative","significant positive")
    )
    out_df <- cbind(spec_df, ctrl_df)
    out_df$outcome <- droplevels(out)
    keep <- vapply(out_df, function(x) length(unique(x)) > 1, logical(1))
    keep["outcome"] <- TRUE
    out_df <- out_df[, keep, drop = FALSE]
    rownames(out_df) <- NULL
    out_df
  }

  # --- one-vs-rest training per sign
  for (sign_label in signs) {
    message("Training NN for outcome: ", sign_label)

    dat <- make_source_dataset(res, focus) %>% stats::na.omit()
    if (!nrow(dat)) next
    dat[] <- lapply(dat, function(x) if (is.factor(x)) droplevels(x) else x)

    y_bin <- as.integer(dat$outcome == sign_label)
    if (length(unique(y_bin)) < 2) next

    X_mm <- model.matrix(~ . - outcome, data = dat)
    colnames(X_mm) <- make.names(colnames(X_mm), unique = TRUE)

    # split train/val/test
    idx <- sample(seq_len(nrow(X_mm)))
    n <- nrow(X_mm)
    n_tr <- floor(split["train"] * n); n_val <- floor(split["val"] * n)
    tr_id <- idx[1:n_tr]; val_id <- idx[(n_tr+1):(n_tr+n_val)]; te_id <- idx[(n_tr+n_val+1):n]

    x_tr <- X_mm[tr_id,,drop=FALSE]; y_tr <- y_bin[tr_id]
    x_val <- X_mm[val_id,,drop=FALSE]; y_val <- y_bin[val_id]
    x_te  <- X_mm[te_id,,drop=FALSE];  y_te  <- y_bin[te_id]

    to_categorical <- function(v) keras::to_categorical(as.integer(v), num_classes = 2)
    y_tr_m <- to_categorical(y_tr); y_val_m <- to_categorical(y_val); y_te_m <- to_categorical(y_te)

    # grid search
    grid_df <- expand.grid(grid, stringsAsFactors = FALSE)
    best <- list(score=-Inf, params=NULL, model=NULL)

    for (i in seq_len(nrow(grid_df))) {
      p <- as.list(grid_df[i,])
      model <- keras::keras_model_sequential() %>%
        keras::layer_dense(units=p$units, activation="relu", input_shape=ncol(x_tr)) %>%
        keras::layer_dropout(rate=p$dropout)
      if (p$layers > 1) {
        for (u in seq_len(p$layers-1)) {
          model %>% keras::layer_dense(units=p$units, activation="relu") %>%
            keras::layer_dropout(rate=p$dropout)
        }
      }
      model %>% keras::layer_dense(units=2, activation="softmax") %>%
        keras::compile(
          optimizer=keras::optimizer_rmsprop(learning_rate=p$learning_rate),
          loss="binary_crossentropy",
          metrics=c("accuracy", keras::metric_auc(name="AUC"))
        )
      hist <- keras::fit(model, x=x_tr, y=y_tr_m,
                         epochs=epochs, batch_size=batch_size, verbose=verbose,
                         validation_data=list(x_val,y_val_m))
      last <- tail(as.data.frame(hist$metrics),1)
      score <- if (!is.null(last$val_AUC)) last$val_AUC else last$val_accuracy
      if (!is.na(score) && score > best$score) {
        best <- list(score=score, params=p, model=model)
      }
    }
    network <- best$model

    pred_wrapper <- function(object, newdata) as.numeric(predict(object, x=as.matrix(newdata))[,2])
    train_df <- as.data.frame(x_tr)
    fip <- vip::vi(object=network, method="shap", train=train_df,
                   feature_names=colnames(train_df), pred_wrapper=pred_wrapper, scale=TRUE)
    fip$Importance <- abs(fip$Importance)

    # map to dimensions (drop "Other")
    fip$dimension <- dplyr::case_when(
      grepl("^dep", fip$Variable)    ~ "Dependent Variable",
      grepl("^fes", fip$Variable)    ~ "Fixed Effects",
      grepl("^set", fip$Variable)    ~ "Standard Error Type",
      grepl("^sample", fip$Variable) ~ "Sample Split",
      grepl("^CONT", fip$Variable)   ~ "Control Set",
      TRUE ~ NA_character_
    )
    fip <- dplyr::filter(fip, !is.na(dimension))

    fip$Variable_clean <- dplyr::case_when(
      grepl("^dep", fip$Variable)    ~ sub("^dep","dep: ", fip$Variable),
      grepl("^fes", fip$Variable)    ~ sub("^fes","fes: ", fip$Variable),
      grepl("^set", fip$Variable)    ~ sub("^set","set: ", fip$Variable),
      grepl("^sample", fip$Variable) ~ sub("^sample","sample: ", fip$Variable),
      grepl("^CONT", fip$Variable)   ~ sub("^CONT","cont: ", fip$Variable),
      TRUE ~ fip$Variable
    )

    fip$share <- 100 * fip$Importance / sum(fip$Importance, na.rm=TRUE)
    agg <- fip %>% dplyr::group_by(dimension) %>%
      dplyr::summarise(share=mean(share, na.rm=TRUE), .groups="drop")

    results[[sign_label]] <- list(sign_label=sign_label, shap_feature=fip, shap_dimension=agg)
  }

  list(method="NN_SHAP", focus=focus, all_signs=results)
}

#' Plot importances (SHAP) of neural network approach
#'
#' @description
#' Displays either dimension-level average SHAP shares (\code{aggregate = TRUE})
#' or feature-level SHAP shares (\code{aggregate = FALSE}), split by outcome class.
#'
#' @param out Output of \code{find_uncertainty_source_neuronet()}.
#' @param aggregate Logical. If \code{TRUE}, plot dimension aggregates.
#'
#' @return A \pkg{ggplot2} object.
#' @export
#' @family plot-sources
plot_importance_neuronet <- function(out, aggregate=TRUE) {
  # (code unchanged)
  if (is.null(out$all_signs) || length(out$all_signs)==0) stop("No NN results found.")
  df_list <- lapply(out$all_signs, function(x) {
    if (aggregate) {
      dplyr::mutate(x$shap_dimension, sign_label=x$sign_label)
    } else {
      dplyr::mutate(x$shap_feature, sign_label=x$sign_label)
    }
  })
  df <- dplyr::bind_rows(df_list)

  if (aggregate) {
    ggplot2::ggplot(df, ggplot2::aes(x=reorder(dimension, share), y=share,
                                     color=sign_label, shape=sign_label)) +
      ggplot2::geom_point(size=3) + ggplot2::coord_flip() +
      ggplot2::labs(x="Specification Dimension", y="Avg. SHAP Importance") +
      ggplot2::theme_minimal() + ggplot2::theme(legend.position="bottom",
                                                legend.title=ggplot2::element_blank())
  } else {
    ggplot2::ggplot(df, ggplot2::aes(x=reorder(Variable_clean, share), y=share,
                                     color=sign_label, shape=sign_label)) +
      ggplot2::geom_point(size=2) + ggplot2::coord_flip() +
      ggplot2::labs(x="Specification Feature", y="SHAP Importance") +
      ggplot2::theme_minimal() + ggplot2::theme(legend.position="bottom",
                                                legend.title=ggplot2::element_blank())
  }
}
