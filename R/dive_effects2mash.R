#' @title Wrapper to run mash given a phenotype data frame
#'
#' @description Though step-by-step GWAS, preparation of mash inputs, and mash
#'     allows you the most flexibility and opportunities to check your results
#'     for errors, once those sanity checks are complete, this function allows
#'     you to go from a phenotype data.frame of a few phenotypes you want to
#'     compare to a mash result. Some exception handling has been built into
#'     this function, but the user should stay cautious and skeptical of any
#'     results that seem 'too good to be true'.
#'
#' @param effects fbm created using 'dive_phe2effects' or 'dive_phe2mash'.
#'     Saved under the name "gwas_effects_{suffix}.rds" and can be loaded into
#'     R using the bigstatsr function "big_attach".
#' @param snp A "bigSNP" object; load with \code{snp_attach()}.
#' @param metadata Metadata created using 'dive_phe2effects' or 'dive_phe2mash'.
#'     Saved under the name "gwas_effects_{suffix}_associated_metadata.csv".
#' @param suffix Optional character vector to give saved files a unique search string/name.
#' @param outputdir Optional file path to save output files.
#' @param thr.r2 Value between 0 and 1. Threshold of r2 measure of linkage
#'     disequilibrium. Markers in higher LD than this will be subset using clumping.
#' @param thr.m "sum" or "max". Type of threshold to use to clump values for
#'     mash inputs. "sum" sums the -log10pvalues for each phenotype and uses
#'     the maximum of this value as the threshold. "max" uses the maximum
#'     -log10pvalue for each SNP across all of the univariate GWAS.
#' @param num.strong Integer. Number of SNPs used to derive data-driven covariance
#'     matrix patterns, using markers with strong effects on phenotypes.
#' @param num.random Integer. Number of SNPs used to derive the correlation structure
#'     of the null tests, and the mash fit on the null tests.
#' @param scale.phe Logical. Should effects for each phenotype be scaled to fall
#'     between -1 and 1? Default is TRUE.
#' @param U.ed Mash data-driven covariance matrices. Specify these as a list or a path
#'     to a file saved as an .rds. Creating these can be time-consuming, and
#'     generating these once and reusing them for multiple mash runs can save time.
#' @param U.hyp Other covariance matrices for mash. Specify these as a list. These
#'     matrices must have dimensions that match the number of phenotypes where
#'     univariate GWAS ran successfully.
#' @param verbose Output some information on the iterations? Default is `TRUE`.
#'
#' @return A mash object made up of all phenotypes where univariate GWAS ran
#'     successfully.
#'
#' @importFrom ashr get_fitted_g
#' @importFrom tibble tibble enframe add_row add_column rownames_to_column
#' @importFrom bigsnpr snp_autoSVD
#' @importFrom dplyr group_by summarise left_join select slice slice_max slice_sample mutate filter
#' @import bigstatsr
#' @import mashr
#' @importFrom cowplot save_plot
#' @importFrom tidyr replace_na
#' @importFrom matrixStats colMaxs rowMaxs
#' @importFrom stats predict
#' @importFrom bigassertr printf
#'
#' @export
dive_effects2mash <- function(effects, snp, metadata, suffix = "", outputdir = ".",
                          thr.r2 = 0.2,
                          thr.m = c("max", "sum"), num.strong = 1000,
                          num.random = NA,
                          scale.phe = TRUE, U.ed = NA,
                          U.hyp = NA, verbose = TRUE){
  # 1. Stop if not functions. ----
  if (attr(snp, "class") != "bigSNP") {
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  if (!dir.exists(outputdir)) {
    dir.create(outputdir)
  }
  if (!grepl("_$", suffix) & suffix != ""){
    suffix <- paste0("_", suffix)
  }

  ## 1a. Generate useful values ----
  G <- snp$genotypes
  #nSNP_M <- round(snp$genotypes$ncol/1000000, digits = 1)
  #nSNP <- paste0(nSNP_M, "_M")
  #if (nSNP_M < 1) {
  #  nSNP_K <- round(snp$genotypes$ncol/1000, digits = 1)
  #  nSNP <- paste0(nSNP_K, "_K")
  #}
  #   nInd <- snp$genotypes$nrow
  #plants <- snp$fam$sample.ID
  #bonferroni <- -log10(0.05/length(snp$map$physical.pos))
  markers <- tibble(CHR = snp$map$chromosome, POS = snp$map$physical.pos,
                    marker.ID = snp$map$marker.ID) %>%
    mutate(CHRN = as.numeric(as.factor(.data$CHR)),
           CHR = as.factor(.data$CHR))
  gwas_ok <- floor(effects$ncol / 3)

  printf2(verbose = verbose, "\nNow preparing gwas effects for use in mash.\n")
  # 4. mash input ----
  ## prioritize effects with max(log10p) or max(sum(log10p))
  ## make a random set of relatively unlinked SNPs
  ind_estim <- (1:(gwas_ok))*3 - 2
  ind_se <- (1:(gwas_ok))*3 - 1
  ind_p <- (1:(gwas_ok))*3
  colnames_fbm <- metadata$phe

  if(effects$ncol %% 3 == 0){
    if (thr.m[1] == "sum") {
      thr_log10p <- big_apply(effects,
                              a.FUN = function(X, ind) rowSums(X[, ind]),
                              ind = ind_p,
                              a.combine = 'plus')
    } else if(thr.m[1] == "max"){
      log10pmax_f <- function(X, ind) rowMaxs(as.matrix(X[, ind]))
      thr_log10p <- big_apply(effects,
                              a.FUN = log10pmax_f,
                              ind = ind_p, a.combine = 'c')
    }
    effects$add_columns(ncol_add = 1)
    colnames_fbm <- c(colnames_fbm, paste0(thr.m[1], "_thr_log10p"))
    effects[,(sum(gwas_ok)*3 + 1)] <- thr_log10p
    effects$save()
  } else if (effects$ncol %% 3 == 1){
    if (thr.m[1] == "sum") {
      thr_log10p <- big_apply(effects,
                              a.FUN = function(X, ind) rowSums(X[, ind]),
                              ind = ind_p,
                              a.combine = 'plus')
    } else if(thr.m[1] == "max"){
      log10pmax_f <- function(X, ind) rowMaxs(as.matrix(X[, ind]))
      thr_log10p <- big_apply(effects,
                              a.FUN = log10pmax_f,
                              ind = ind_p, a.combine = 'c')
    }
    colnames_fbm <- c(colnames_fbm, paste0(thr.m[1], "_thr_log10p"))
    effects[,(sum(gwas_ok)*3 + 1)] <- thr_log10p
    effects$save()
  } else {
    stop("Effect df should have a multiple of three columns, or can have one additional
         column for the p value threshold.")
  }

  ## replace NA or Nan values
  # Replace SE with 1's, estimates and p values with 0's.
  replace_na_1 <- function(X, ind) replace_na(X[, ind], 1)
  replace_na_0 <- function(X, ind) replace_na(X[, ind], 0)
  effects[, ind_se] <- big_apply(effects, a.FUN = replace_na_1, ind = ind_se,
                               a.combine = 'plus')
  effects[, ind_estim] <- big_apply(effects, a.FUN = replace_na_0, ind = ind_estim,
                                  a.combine = 'plus')
  effects[, ind_p] <- big_apply(effects, a.FUN = replace_na_0, ind = ind_p,
                              a.combine = 'plus')
  effects[, (sum(gwas_ok)*3+1)] <- big_apply(effects, a.FUN = replace_na_0,
                                           ind = (sum(gwas_ok)*3 + 1),
                                           a.combine = 'plus')
  effects$save()

  strong_clumps <- snp_clumping(G, infos.chr = markers$CHRN, thr.r2 = thr.r2,
                                infos.pos = markers$POS, S = thr_log10p)
  random_clumps <- snp_clumping(G, infos.chr = markers$CHRN, thr.r2 = thr.r2,
                                infos.pos = markers$POS)
  # this should be a top_n (slice_min/slice_max/slice_sample) with numSNPs, not a quantile
  strong_sample <- add_column(markers, thr_log10p) %>%
    rownames_to_column(var = "value") %>%
    mutate(value = as.numeric(.data$value)) %>%
    filter(.data$value %in% strong_clumps) %>%
    slice_max(order_by = .data$thr_log10p, n = num.strong) %>%
    arrange(.data$value)

  if (is.na(num.random)[1]) {
    num.random <- num.strong*2
  }
  random_sample <- add_column(markers, thr_log10p) %>%
    rownames_to_column(var = "value") %>%
    mutate(value = as.numeric(.data$value)) %>%
    filter(!is.na(.data$thr_log10p)) %>%
    filter(.data$value %in% random_clumps) %>%
    slice_sample(n = num.random) %>%
    arrange(.data$value)

  ## scaling
  if (scale.phe == TRUE) {
    colmaxes <- function(X, ind) colMaxs(abs(as.matrix(X[, ind])))
    scale.effects <- big_apply(effects, a.FUN = colmaxes,
                               ind = ind_estim, a.combine = 'c')
    colstand <- function(X, ind, v) X[,ind] / v
    for (j in seq_along(scale.effects)) {  # standardize one gwas at a time.
      effects[,c(ind_estim[j], ind_se[j])] <-
        big_apply(effects, a.FUN = colstand, ind = c(ind_estim[j], ind_se[j]),
                  v = scale.effects[j], a.combine = 'plus')
    }
    effects$save()
    gwas_metadata <- metadata %>% mutate(scaled = TRUE)
  } else {
    gwas_metadata <- metadata %>% mutate(scaled = FALSE)
  }

  write_csv(tibble(colnames_fbm), file.path(outputdir,
                                            paste0("gwas_effects_", suffix,
                                                   "_column_names.csv")))
  write_csv(gwas_metadata, file.path(outputdir,
                                     paste0("gwas_effects_", suffix,
                                            "_associated_metadata.csv")))
  ## make mash input data.frames (6x or more)

  Bhat_strong <- as.matrix(effects[strong_sample$value, ind_estim], )
  Shat_strong <- as.matrix(effects[strong_sample$value, ind_se])

  Bhat_random <- as.matrix(effects[random_sample$value, ind_estim])
  Shat_random <- as.matrix(effects[random_sample$value, ind_se])

  ## name the columns for these conditions (usually the phenotype)
  colnames(Bhat_strong) <- gwas_metadata$phe
  colnames(Shat_strong) <- gwas_metadata$phe
  colnames(Bhat_random) <- gwas_metadata$phe
  colnames(Shat_random) <- gwas_metadata$phe

  # 5. mash ----

  data_r <- mashr::mash_set_data(Bhat_random, Shat_random)
  printf2(verbose = verbose, "\nEstimating correlation structure in the null tests from a random sample of clumped data.\n")
  Vhat <- mashr::estimate_null_correlation_simple(data = data_r)

  data_strong <- mashr::mash_set_data(Bhat_strong, Shat_strong, V = Vhat)
  data_random <- mashr::mash_set_data(Bhat_random, Shat_random, V = Vhat)
  U_c <- mashr::cov_canonical(data_random)

  if (is.na(U.ed[1])) {
    printf2(verbose = verbose, "\nNow estimating data-driven covariances using
    the strong tests. NB: This step may take some time to complete.\n")
    if (length(ind_p) < 6) {
      cov_npc <- length(ind_p) - 1
    } else {
      cov_npc <- 5
    }
    U_pca = mashr::cov_pca(data_strong, npc = cov_npc)
    U_ed = mashr::cov_ed(data_strong, U_pca)
    saveRDS(U_ed, file = file.path(outputdir, paste0("Mash_U_ed", suffix,
                                                     ".rds")))
  } else if (typeof(U.ed) == "list") {
    U_ed <- U.ed
  } else if (typeof(U.ed) == "character") {
    U_ed <- readRDS(file = U.ed)
  } else {
    stop("U.ed should be NA, a list created using 'mashr::cov_ed', ",
         "or a file path of a U_ed saved as an .rds")
  }

  if (typeof(U.hyp) == "list") {
    m = mashr::mash(data_random, Ulist = c(U_ed, U_c, U.hyp), outputlevel = 1)
  } else if (typeof(U.hyp) == "character") {
    U_hyp <- readRDS(file = U.hyp)
    m = mashr::mash(data_random, Ulist = c(U_ed, U_c, U_hyp), outputlevel = 1)
  } else {
    m = mashr::mash(data_random, Ulist = c(U_ed, U_c), outputlevel = 1)
    printf2(verbose = verbose, "\nNo user-specified covariance matrices were included in the mash fit.")
  }

  printf2(verbose = verbose, "\nComputing posterior weights for all effects
  using the mash fit from the random tests.")
  ## Batch process SNPs through this, don't run on full set if > 20000 rows.
  ## Even for 1M SNPs, because computing posterior weights scales quadratically
  ## with the number of rows in Bhat and Shat. 10K  = 13s, 20K = 55s; 40K = 218s
  ## By my calc, this starts getting unwieldy between 4000 and 8000 rows.
  ## See mash issue: https://github.com/stephenslab/mashr/issues/87
  if(effects$nrow > 20000){
    subset_size <- 4000
    n_subsets <- ceiling(effects$nrow / subset_size)
    printf2(verbose = verbose, "\nSplitting data into %s sets of 4K markers to speed computation.\n",
            n_subsets)
    for (i in 1:n_subsets) {
      if(i < n_subsets){
        from <- (i*subset_size - (subset_size - 1))
        to <- i*subset_size
        row_subset <- from:to
      } else {
        from <- n_subsets*subset_size - (subset_size - 1)
        to <- effects$nrow
        row_subset <- from:to
      }
      Bhat_subset <- as.matrix(effects[row_subset, ind_estim])
      Shat_subset <- as.matrix(effects[row_subset, ind_se])
      colnames(Bhat_subset) <- gwas_metadata$phe
      colnames(Shat_subset) <- gwas_metadata$phe
      data_subset <- mashr::mash_set_data(Bhat_subset, Shat_subset, V = Vhat)
      m_subset = mashr::mash(data_subset, g = ashr::get_fitted_g(m), fixg = TRUE)

      if (i == 1){
        m2 <- m_subset
      } else {     # make a new mash object with the combined data.
        PosteriorMean = rbind(m2$result$PosteriorMean, m_subset$result$PosteriorMean)
        PosteriorSD = rbind(m2$result$PosteriorSD, m_subset$result$PosteriorSD)
        lfdr = rbind(m2$result$lfdr, m_subset$result$lfdr)
        NegativeProb = rbind(m2$result$NegativeProb, m_subset$result$NegativeProb)
        lfsr = rbind(m2$result$lfsr, m_subset$result$lfsr)
        posterior_matrices = list(PosteriorMean = PosteriorMean,
                                  PosteriorSD = PosteriorSD,
                                  lfdr = lfdr,
                                  NegativeProb = NegativeProb,
                                  lfsr = lfsr)
        loglik = m2$loglik # NB must recalculate from sum(vloglik) at end
        vloglik = rbind(m2$vloglik, m_subset$vloglik)
        null_loglik = c(m2$null_loglik, m_subset$null_loglik)
        alt_loglik = rbind(m2$alt_loglik, m_subset$alt_loglik)
        fitted_g = m2$fitted_g        # all four components are equal
        posterior_weights = rbind(m2$posterior_weights, m_subset$posterior_weights)
        alpha = m2$alpha  # equal
        m2 = list(result = posterior_matrices,
                  loglik = loglik,
                  vloglik = vloglik,
                  null_loglik = null_loglik,
                  alt_loglik = alt_loglik,
                  fitted_g = fitted_g,
                  posterior_weights = posterior_weights,
                  alpha = alpha)
        class(m2) = "mash"
      }
    }
    loglik = sum(m2$vloglik)
    m2$loglik <- loglik
    # total loglik in mash function is: loglik = sum(vloglik)
  } else {
    Bhat_full <- as.matrix(effects[, ind_estim])
    Shat_full <- as.matrix(effects[, ind_se])
    colnames(Bhat_full) <- gwas_metadata$phe
    colnames(Shat_full) <- gwas_metadata$phe
    data_full <- mashr::mash_set_data(Bhat_full, Shat_full, V = Vhat)
    m2 = mashr::mash(data_full, g = ashr::get_fitted_g(m), fixg = TRUE)
  }

  return(m2)
}


