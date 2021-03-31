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
#' @param df Dataframe containing phenotypes for mash where the first column is
#'     'sample.ID', which should match values in the snp$fam$sample.ID column.
#' @param snp A "bigSNP" object; load with \code{snp_attach()}.
#' @param type Character string, or a character vector the length of the number
#'     of phenotypes. Type of univarate regression to run for GWAS.
#'     Options are "linear" or "logistic".
#' @param svd A "big_SVD" object; Optional covariance matrix to use for
#'     population structure correction.
#' @param suffix Optional character vector to give saved files a unique search string/name.
#' @param outputdir Optional file path to save output files.
#' @param min.phe Integer. Minimum number of individuals phenotyped in order to
#'     include that phenotype in GWAS. Default is 200. Use lower values with
#'     caution.
#' @param save.plots Logical. Should Manhattan and QQ-plots be generated and
#'     saved to the working directory for univariate GWAS? Default is TRUE.
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
#' @param roll.size Integer. Used to create the svd for GWAS.
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
dive_phe2mash <- function(df, snp, type = "linear", svd = NULL, suffix = "",
                          outputdir = ".",
                          min.phe = 200, save.plots = TRUE, thr.r2 = 0.2,
                          thr.m = c("sum", "max"), num.strong = 1000,
                          num.random = NA,
                          scale.phe = TRUE, roll.size = 50, U.ed = NA,
                          U.hyp = NA, verbose = TRUE){
  # 1. Stop if not functions. ----
  if (attr(snp, "class") != "bigSNP") {
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  if (colnames(df)[1] != "sample.ID") {
    stop("First column of phenotype dataframe (df) must be 'sample.ID'.")
  }
  if (length(type) > 1) {
    if (length(type) != ncol(df) - 1) {
      stop(paste0("Specify either one GWAS type (type = 'linear' or type = ",
      "'logistic'), or one type for each phenotype in 'df'."))
    }
  } else {
    type <- rep(type, ncol(df) - 1)
  }
  if (!dir.exists(outputdir)) {
    dir.create(outputdir)
  }

  ## 1a. Generate useful values ----
  G <- snp$genotypes
  nSNP_M <- round(snp$genotypes$ncol/1000000, digits = 1)
  nSNP <- paste0(nSNP_M, "_M")
  if (nSNP_M < 1) {
    nSNP_K <- round(snp$genotypes$ncol/1000, digits = 1)
    nSNP <- paste0(nSNP_K, "_K")
  }
  #   nInd <- snp$genotypes$nrow
  plants <- snp$fam$sample.ID
  bonferroni <- -log10(0.05/length(snp$map$physical.pos))
  markers <- tibble(CHR = snp$map$chromosome, POS = snp$map$physical.pos,
                    marker.ID = snp$map$marker.ID) %>%
    mutate(CHRN = as.numeric(as.factor(.data$CHR)),
           CHR = as.factor(.data$CHR))

  # 2. Pop Structure Correction  ----
  if (is.null(svd)) {
    printf2(verbose = verbose, "\nCovariance matrix (svd) was not supplied - ")
    printf2(verbose = verbose, "\nthis will be generated using snp_autoSVD()")
    svd <- snp_autoSVD(G = G, infos.chr = markers$CHRN, infos.pos = markers$POS,
                       k = 10, thr.r2 = thr.r2, roll.size = roll.size)
    } else {
    stopifnot(attr(svd, "class") == "big_SVD")
      }
  pc_max <- ncol(svd$u)
  gwas_ok <- c()
  first_gwas_ok <- FALSE

  for (i in 2:ncol(df)) {
    df1 <- df %>%
      dplyr::select(.data$sample.ID, all_of(i))
    phename <- names(df1)[2]
    df1 <- df1 %>%
      group_by(.data$sample.ID) %>%
      filter(!is.na(.data[[phename]])) %>%
      summarise(phe = mean(.data[[phename]]), .groups = "drop_last")
    df1 <- plants %>%
      enframe(name = NULL, value = "sample.ID") %>%
      mutate(sample.ID = as.character(.data$sample.ID)) %>%
      left_join(df1, by = "sample.ID")
    nPhe <- length(which(!is.na(df1[,2])))
    nLev <- nrow(unique(df1[which(!is.na(df1[,2])),2]))

    # Checks for correct combinations of phenotypes and GWAS types.
    gwas_ok[i-1] <- check_gwas(df1 = df1, phename = phename, type = type[i-1],
                               nPhe = nPhe, minphe = min.phe, nLev = nLev)


    # Find best # PCs to correct for population structure for each phenotype.
    if(gwas_ok[i-1]){

      lambdagc_df <- div_lambda_GC(df = df1, type = type[i-1], snp = snp,
                                   svd = svd, npcs = c(0:pc_max))
      PC_df <- get_best_PC_df(lambdagc_df)
      PC_df <- PC_df[1,]

  # 3. GWAS  ----

    # run gwas using best npcs from step 2 (best pop structure correction)
    gwas <- div_gwas(df = df1, snp = snp, type = type[i - 1], svd = svd,
                     npcs = PC_df$NumPCs)
    gwas <- gwas %>% mutate(pvalue = predict(gwas, log10 = FALSE),
                            log10p = -log10(.data$pvalue))
    gwas_data <- tibble(phe = phename, type = type[i - 1],
                        nsnp = nSNP, npcs = PC_df$NumPCs, nphe = nPhe,
                        nlev = nLev, lambda_GC = PC_df$lambda_GC,
                        bonferroni = bonferroni)

    # plot  QQ if save.plots == TRUE
    if (save.plots == TRUE) {
      qqplot <- get_qqplot(ps = gwas$pvalue, lambdaGC = TRUE)
      }

    # save gwas outputs together in a fbm
    gwas <- gwas %>%
      select(.data[["estim"]], .data[["std.err"]], .data[["log10p"]])

    if(!first_gwas_ok){     # save .bk and .rds file the first time through the loop.
      if (!grepl("_$", suffix) & suffix != ""){
        suffix <- paste0("_", suffix)
      }
      first_gwas_ok <- TRUE
      fbm.name <- file.path(outputdir, paste0("gwas_effects", suffix))

      colnames_fbm <- c(paste0(phename, "_Effect"), paste0(phename, "_SE"),
                       paste0(phename, "_log10p"))
      as_FBM(gwas, backingfile = fbm.name)$save()
      gwas2 <- big_attach(paste0(fbm.name, ".rds"))
      gwas_metadata <- gwas_data

    } else {
      colnames_fbm <- c(colnames_fbm, paste0(phename, "_Effect"),
                        paste0(phename, "_SE"), paste0(phename, "_log10p"))
      gwas2$add_columns(ncol_add = 3)
      gwas2[, c(sum(gwas_ok)*3 - 2, sum(gwas_ok)*3 - 1,
                sum(gwas_ok)*3)] <- gwas
      gwas2$save()
      gwas_metadata <- add_row(gwas_metadata, phe = phename, type = type[i - 1],
                           nsnp = nSNP, npcs = PC_df$NumPCs, nphe = nPhe,
                           nlev = nLev, lambda_GC = PC_df$lambda_GC,
                           bonferroni = bonferroni)
    }
    # plot Manhattan and QQ if save.plots == TRUE
    if (save.plots == TRUE) {
      # set aspect ratio based on number of SNPs in snp file
      asp <- log10(snp$genotypes$ncol)/2
      if(asp < 1.1){
        asp <- 1.1
      }

      manhattan <- get_manhattan(X = gwas2, ind = sum(gwas_ok)*3, snp = snp,
                                 thresh = bonferroni)
      plotname <- paste0(gwas_data$phe, "_", gwas_data$type, "_model_",
                         gwas_data$nphe, "g_", gwas_data$nsnp, "_SNPs_",
                         gwas_data$npcs, "_PCs.png")
      save_plot(filename = file.path(outputdir, paste0("QQplot_", plotname)),
                plot = qqplot, base_asp = 1, base_height = 4)
      save_plot(filename = file.path(outputdir, paste0("Manhattan_", plotname)),
                plot = manhattan, base_asp = asp, base_height = 3.75)

    }
  rm(gwas)
  printf2(verbose = verbose, "\nFinished GWAS on phenotype %s. ",
          names(df)[i])
    } else {
      printf2(verbose = verbose, "\nSkipping GWAS on phenotype %s. ",
              names(df)[i])
    }
  }

  printf2(verbose = verbose, "\nNow preparing gwas effects for use in mash.\n")
  # 4. mash input ----
      ## prioritize effects with max(log10p) or max(sum(log10p))
      ## make a random set of relatively unlinked SNPs
  ind_estim <- (1:sum(gwas_ok))*3 - 2
  ind_se <- (1:sum(gwas_ok))*3 - 1
  ind_p <- (1:sum(gwas_ok))*3

  if (thr.m == "sum") {
  thr_log10p <- big_apply(gwas2,
                         a.FUN = function(X, ind) rowSums(X[, ind]),
                         ind = ind_p,
                         a.combine = 'plus')
  } else if(thr.m == "max"){
    log10pmax_f <- function(X, ind) rowMaxs(as.matrix(X[, ind]))
    thr_log10p <- big_apply(gwas2,
                           a.FUN = log10pmax_f,
                           ind = ind_p, a.combine = 'c')
  }
  gwas2$add_columns(ncol_add = 1)
  colnames_fbm <- c(colnames_fbm, paste0(thr.m, "_thr_log10p"))
  gwas2[,(sum(gwas_ok)*3 + 1)] <- thr_log10p
  gwas2$save()

    ## replace NA or Nan values
  # Replace SE with 1's, estimates and p values with 0's.
  replace_na_1 <- function(X, ind) replace_na(X[, ind], 1)
  replace_na_0 <- function(X, ind) replace_na(X[, ind], 0)
  gwas2[, ind_se] <- big_apply(gwas2, a.FUN = replace_na_1, ind = ind_se,
                              a.combine = 'plus')
  gwas2[, ind_estim] <- big_apply(gwas2, a.FUN = replace_na_0, ind = ind_estim,
                                 a.combine = 'plus')
  gwas2[, ind_p] <- big_apply(gwas2, a.FUN = replace_na_0, ind = ind_p,
                                 a.combine = 'plus')
  gwas2[, (sum(gwas_ok)*3+1)] <- big_apply(gwas2, a.FUN = replace_na_0,
                                          ind = (sum(gwas_ok)*3 + 1),
                                          a.combine = 'plus')
  gwas2$save()

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
    scale.effects <- big_apply(gwas2, a.FUN = colmaxes,
                               ind = ind_estim, a.combine = 'c')
    colstand <- function(X, ind, v) X[,ind] / v
    for (j in seq_along(scale.effects)) {  # standardize one gwas at a time.
      gwas2[,c(ind_estim[j], ind_se[j])] <-
        big_apply(gwas2, a.FUN = colstand, ind = c(ind_estim[j], ind_se[j]),
                  v = scale.effects[j], a.combine = 'plus')
      }
    gwas2$save()
    gwas_metadata <- gwas_metadata %>% mutate(scaled = TRUE)
  } else {
    gwas_metadata <- gwas_metadata %>% mutate(scaled = FALSE)
  }

  write_csv(tibble(colnames_fbm), file.path(outputdir,
                                            paste0("gwas_effects", suffix,
                                                   "_column_names.csv")))
  write_csv(gwas_metadata, file.path(outputdir,
                                     paste0("gwas_effects", suffix,
                                            "_associated_metadata.csv")))
    ## make mash input data.frames (6x or more)

    Bhat_strong <- as.matrix(gwas2[strong_sample$value, ind_estim], )
    Shat_strong <- as.matrix(gwas2[strong_sample$value, ind_se])

    Bhat_random <- as.matrix(gwas2[random_sample$value, ind_estim])
    Shat_random <- as.matrix(gwas2[random_sample$value, ind_se])

    ## Full data: Both Bhat and Shat are zero (or near zero) for some input data.
    ## Filter this data from the input, or set Shat to a positive number to
    ## avoid numerical issues. which rowSums are 0, filter these out or make +.
    ## Eventually want to batch process SNPs through this, not make a full set.
    Bhat_full <- as.matrix(gwas2[, ind_estim])
    Shat_full <- as.matrix(gwas2[, ind_se])

    ## name the columns for these conditions (usually the phenotype)
    colnames(Bhat_strong) <- gwas_metadata$phe
    colnames(Shat_strong) <- gwas_metadata$phe
    colnames(Bhat_random) <- gwas_metadata$phe
    colnames(Shat_random) <- gwas_metadata$phe
    colnames(Bhat_full) <- gwas_metadata$phe
    colnames(Shat_full) <- gwas_metadata$phe

  # 5. mash ----

  data_r <- mashr::mash_set_data(Bhat_random, Shat_random)
  printf2(verbose = verbose, "\nEstimating correlation structure in the null tests from a random sample of clumped data.")
  Vhat <- mashr::estimate_null_correlation_simple(data = data_r)

  data_strong <- mashr::mash_set_data(Bhat_strong, Shat_strong, V = Vhat)
  data_random <- mashr::mash_set_data(Bhat_random, Shat_random, V = Vhat)
  data_full <- mashr::mash_set_data(Bhat_full, Shat_full, V = Vhat)
  U_c <- mashr::cov_canonical(data_random)

  if (is.na(U.ed[1])) {
    printf2(verbose = verbose, "\nNow estimating data-driven covariances using
    the strong tests. NB: This step may take some time to complete.\n")
    if (length(ind_p) < 6) {
      cov_npc <- ind_p - 1
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
  m2 = mashr::mash(data_full, g = ashr::get_fitted_g(m), fixg = TRUE)

  return(m2)
}




#' Wrapper for bigsnpr for GWAS
#'
#' @description Given a dataframe of phenotypes associated with sample.IDs, this
#'     function is a wrapper around bigsnpr functions to conduct linear or
#'     logistic regression on wheat. The main advantages of this
#'     function over just using the bigsnpr functions is that it automatically
#'     removes individual genotypes with missing phenotypic data
#'     and that it can run GWAS on multiple phenotypes sequentially.
#'
#' @param df Dataframe of phenotypes where the first column is sample.ID
#' @param type Character string. Type of univarate regression to run for GWAS.
#'     Options are "linear" or "logistic".
#' @param snp Genomic information to include for wheat.
#' @param svd Optional covariance matrix to include in the regression. You
#'     can generate these using \code{bigsnpr::snp_autoSVD()}.
#' @param npcs Integer. Number of PCs to use for population structure correction.
#'
#' @import bigsnpr
#' @import bigstatsr
#' @importFrom dplyr mutate rename case_when
#' @importFrom purrr as_vector
#' @importFrom tibble as_tibble enframe
#' @importFrom rlang .data
#'
#' @return The gwas results for the last phenotype in the dataframe. That
#'     phenotype, as well as the remaining phenotypes, are saved as RDS objects
#'     in the working directory.
#'
#' @export
div_gwas <- function(df, snp, type, svd, npcs){
  stopifnot(type %in% c("linear", "logistic"))
  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  if(colnames(df)[1] != "sample.ID"){
    stop("First column of phenotype dataframe (df) must be 'sample.ID'.")
  }
  G <- snp$genotypes
  pc_max = ncol(svd$u)

  for(i in seq_along(names(df))[-1]){
    y1 <- as_vector(df[which(!is.na(df[,i])), i])
    ind_y <- which(!is.na(df[,i]))

    if(type == "linear"){
      if(npcs > 0){
        ind_u <- matrix(svd$u[which(!is.na(df[,i])),1:npcs], ncol = npcs)
        gwaspc <- big_univLinReg(G, y.train = y1, covar.train = ind_u,
                                 ind.train = ind_y, ncores = 1)
      } else {
        gwaspc <- big_univLinReg(G, y.train = y1, ind.train = ind_y,
                                 ncores = 1)
      }
    } else if(type == "logistic"){
      message(paste0("For logistic models, if convergence is not reached by ",
                     "the main algorithm for any SNP, the corresponding `niter` element ",
                     "is set to NA, and glm is used instead. If glm can't ",
                     "converge either, those SNP estimations are set to NA."))
      if(npcs > 0){
        ind_u <- matrix(svd$u[which(!is.na(df[,i])),1:npcs], ncol = npcs)
        gwaspc <- suppressMessages(big_univLogReg(G, y01.train = y1,
                                                  covar.train = ind_u,
                                                  ind.train = ind_y,
                                                  ncores = 1))
      } else {
        gwaspc <- suppressMessages(big_univLogReg(G, y01.train = y1,
                                                  ind.train = ind_y,
                                                  ncores = 1))
      }
    } else {
      stop(paste0("Type of GWAS not recognized: please choose one of 'linear'",
                  " or 'logistic'"))
    }

  }
  return(gwaspc)
}

#' Verbose?
#' @importFrom bigassertr printf
printf2 <- function(verbose, ...) if (verbose) { printf(...) }

#' Create a quantile-quantile plot with ggplot2.
#'
#' @description Assumptions for this quantile quantile plot:
#'     Expected P values are uniformly distributed.
#'     Confidence intervals assume independence between tests.
#'     We expect deviations past the confidence intervals if the tests are
#'     not independent.
#'     For example, in a genome-wide association study, the genotype at any
#'     position is correlated to nearby positions. Tests of nearby genotypes
#'     will result in similar test statistics.
#'
#' @param ps Numeric vector of p-values.
#' @param ci Numeric. Size of the confidence interval, 0.95 by default.
#' @param lambdaGC Logical. Add the Genomic Control coefficient as subtitle to
#'     the plot?
#'
#' @import ggplot2
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @importFrom stats qbeta ppoints
#' @param tol Numeric. Tolerance for optional Genomic Control coefficient.
#'
#' @return A ggplot2 plot.
#'
#' @export
get_qqplot <- function(ps, ci = 0.95, lambdaGC = FALSE, tol = 1e-8) {
  ps <- ps[which(!is.na(ps))]
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  df_round <- round_xy(df$expected, df$observed, cl = df$clower, cu = df$cupper)
  log10Pe <- expression(paste("Expected -log"[10], plain("("), italic(p-value),
                              plain(")")))
  log10Po <- expression(paste("Observed -log"[10], plain("("), italic(p-value),
                              plain(")")))
  p1 <- ggplot(as_tibble(df_round)) +
    geom_point(aes(.data$expected, .data$observed), shape = 1, size = 1) +
    geom_abline(intercept = 0, slope = 1, size = 1.5, color = "red") +
    geom_line(aes(.data$expected, .data$cupper), linetype = 2) +
    geom_line(aes(.data$expected, .data$clower), linetype = 2) +
    xlab(log10Pe) +
    ylab(log10Po) +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.line.x = element_line(size = 0.35, colour = 'grey50'),
          axis.line.y = element_line(size = 0.35, colour = 'grey50'),
          axis.ticks = element_line(size = 0.25, colour = 'grey50'),
          legend.justification = c(1, 0.75), legend.position = c(1, 0.9),
          legend.key.size = unit(0.35, 'cm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 9),
          legend.text.align = 0, legend.background = element_blank(),
          plot.subtitle = element_text(size = 10, vjust = 0),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0.5, size = 10 ,vjust = 0),
          strip.placement = 'outside', panel.spacing.x = unit(-0.4, 'cm'))

  if (lambdaGC) {
    lamGC <- get_lambdagc(ps = ps, tol = tol)
    expr <- substitute(expression(lambda[GC] == l), list(l = lamGC))
    p1 + labs(subtitle = eval(expr))
  } else {
    p1
  }
}


#' Return a number rounded to some number of digits
#'
#' @description Given some x, return the number rounded to some number of
#'     digits.
#'
#' @param x A number or vector of numbers
#' @param at Numeric. Rounding factor or size of the bin to round to.
#'
#' @return A number or vector of numbers
round2 <- function(x, at) ceiling(x / at) * at

#' Return a dataframe binned into 2-d bins by some x and y.
#'
#' @description Given a dataframe of x and y values (with some optional
#'     confidence intervals surrounding the y values), return only the unique
#'     values of x and y in some set of 2-d bins.
#'
#' @param x Numeric vector. The first vector for binning.
#' @param y Numeric vector. the second vector for binning
#' @param cl Numeric vector. Optional confidence interval for the y vector,
#'     lower bound.
#' @param cu Numeric vector. Optional confidence interval for the y vector,
#'     upper bound.
#' @param roundby Numeric. The amount to round the x and y vectors by for 2d
#'     binning.
#'
#' @return A dataframe containing the 2-d binned values for x and y, and their
#'     confidence intervals.
round_xy <- function(x, y, cl = NA, cu = NA, roundby = 0.001){
  expected <- round2(x, at = roundby)
  observed <- round2(y, at = roundby)
  if(!is.na(cl[1]) & !is.na(cu[1])){
    clower <- round2(cl, at = roundby)
    cupper <- round2(cu, at = roundby)
    tp <- cbind(expected, observed, clower, cupper)
    return(tp[!duplicated(tp),])
  } else {
    tp <- cbind(expected, observed)
    return(tp[!duplicated(tp),])
  }
}

get_manhattan <- function(X, ind, snp, thresh){
  roundFBM <- function(X, ind, at) ceiling(X[, ind] / at) * at
  observed <- big_apply(X, ind = ind, a.FUN = roundFBM, at = 0.01,
                        a.combine = 'plus')

  plot_data <- tibble(CHR = snp$map$chromosome, POS = snp$map$physical.pos,
                    marker.ID = snp$map$marker.ID, observed = observed)

  if (length(unique(snp$map$physical.pos)) >= 500000) {
    plot_data <- plot_data %>%
      mutate(POS = round2(.data$POS, at = 250000))
  }
  plot_data <- plot_data %>%
    group_by(.data$CHR, .data$POS, .data$observed) %>%
    slice(1) %>%
    mutate(CHR = as.factor(.data$CHR))

  nchr <- length(unique(plot_data$CHR))

  p1 <- plot_data %>%
    ggplot(aes(x = .data$POS, y = .data$observed)) +
    geom_point(aes(color = .data$CHR, fill = .data$CHR)) +
    geom_hline(yintercept = thresh, color = "black", linetype = 2,
               size = 1) +
    facet_wrap(~ .data$CHR, nrow = 1, scales = "free_x",
               strip.position = "bottom") +
    scale_color_manual(values = rep(c("#1B0C42FF", "#48347dFF",
                                      "#95919eFF"), ceiling(nchr/3)),
                       guide = FALSE) +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill=NA),
          legend.position = "none",
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.line.x = element_line(size = 0.35, colour = 'grey50'),
          axis.line.y = element_line(size = 0.35, colour = 'grey50'),
          axis.ticks = element_line(size = 0.25, colour = 'grey50'),
          legend.justification = c(1, 0.75),
          legend.key.size = unit(0.35, 'cm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 9),
          legend.text.align = 0, legend.background = element_blank(),
          plot.subtitle = element_text(size = 10, vjust = 0),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0.5, size = 10 ,vjust = 0),
          strip.placement = 'outside', panel.spacing.x = unit(-0.1, 'cm')) +
    labs(x = "Chromosome", y = "-log10(p value)") +
    scale_x_continuous(expand = c(0.15, 0.15))
  return(p1)
}

#' Return lambda_GC for different numbers of PCs for GWAS on Panicum virgatum.
#'
#' @description Given a dataframe of phenotypes associated with sample.IDs and
#'     output from a PCA to control for population structure, this function will
#'     return a .csv file of the lambda_GC values for the GWAS upon inclusion
#'     of different numbers of PCs. This allows the user to choose a number of
#'     PCs that returns a lambda_GC close to 1, and thus ensure that they have
#'     done adequate correction for population structure.
#'
#' @param df Dataframe of phenotypes where the first column is sample.ID and each
#'     sample.ID occurs only once in the dataframe.
#' @param type Character string. Type of univarate regression to run for GWAS.
#'     Options are "linear" or "logistic".
#' @param snp A bigSNP object with sample.IDs that match the df.
#' @param svd big_SVD object; Covariance matrix to include in the regression.
#'      Generate these using \code{bigsnpr::snp_autoSVD()}.
#' @param ncores Number of cores to use. Default is one.
#' @param npcs Integer vector of principle components to use.
#'     Defaults to c(0:10).
#' @param saveoutput Logical. Should output be saved as a csv to the
#'     working directory?
#'
#' @import bigsnpr
#' @import bigstatsr
#' @importFrom dplyr mutate rename case_when mutate_if
#' @importFrom purrr as_vector
#' @importFrom tibble as_tibble enframe
#' @importFrom rlang .data
#' @importFrom readr write_csv
#' @importFrom utils tail
#'
#' @return A dataframe containing the lambda_GC values for each number of PCs
#'     specified. This is also saved as a .csv file in the working directory.
#'
#' @export
div_lambda_GC <- function(df, type = c("linear", "logistic"), snp,
                            svd = NA, ncores = 1, npcs = c(0:10),
                            saveoutput = FALSE){
  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  if(colnames(df)[1] != "sample.ID"){
    stop("First column of phenotype dataframe (df) must be 'sample.ID'.")
  }
  if(length(svd) == 1){
    stop(paste0("Need to specify covariance matrix (svd) and a vector of",
                " PC #'s to test (npcs)."))
  }


  G <- snp$genotypes

  LambdaGC <- as_tibble(matrix(data =
                                 c(npcs, rep(NA, (ncol(df) - 1)*length(npcs))),
                               nrow = length(npcs), ncol = ncol(df),
                               dimnames = list(npcs, colnames(df))))
  LambdaGC <- LambdaGC %>%
    dplyr::rename("NumPCs" = .data$sample.ID) %>%
    mutate_if(is.integer, as.numeric)

  for (i in seq_along(names(df))[-1]) {

    for (k in c(1:length(npcs))) {

      if (type == "linear") {

        y1 <- as_vector(df[which(!is.na(df[,i])), i])
        ind_y <- which(!is.na(df[,i]))

        if (npcs[k] == 0) {

          gwaspc <- big_univLinReg(G, y.train = y1, ind.train = ind_y,
                                   ncores = ncores)
        } else {

          ind_u <- matrix(svd$u[which(!is.na(df[,i])),1:npcs[k]],
                          ncol = npcs[k])
          gwaspc <- big_univLinReg(G, y.train = y1, covar.train = ind_u,
                                   ind.train = ind_y, ncores = ncores)
        }
      } else if(type == "logistic"){
        message(paste0("For logistic models, if convergence is not reached by ",
                       "the main algorithm for some SNPs, the corresponding `niter` element ",
                       "is set to NA, and glm is used instead. If glm can't ",
                       "converge either, those SNP estimations are set to NA."))
        y1 <- as_vector(df[which(!is.na(df[,i])), i])
        ind_y <- which(!is.na(df[,i]))
        if(npcs[k] == 0){
          gwaspc <- suppressMessages(big_univLogReg(G, y01.train = y1,
                                                    ind.train = ind_y,
                                                    ncores = ncores))
        } else {
          ind_u <- matrix(svd$u[which(!is.na(df[,i])),1:npcs[k]],
                          ncol = npcs[k])
          gwaspc <- suppressMessages(big_univLogReg(G, y01.train = y1,
                                                    covar.train = ind_u,
                                                    ind.train = ind_y,
                                                    ncores = ncores))
        }
      }
      ps <- predict(gwaspc, log10 = FALSE)
      LambdaGC[k,i] <- get_lambdagc(ps = ps)
      #message(paste0("Finished Lambda_GC calculation for ", names(df)[i],
      #               " using ", npcs[k], " PCs."))
    }

    if(saveoutput == TRUE){
      write_csv(LambdaGC, path = paste0("Lambda_GC_", names(df)[i], ".csv"))
    }
    #message(paste0("Finished phenotype ", i-1, ": ", names(df)[i]))
  }
  if(saveoutput == TRUE){
    write_csv(LambdaGC, path = paste0("Lambda_GC_", names(df)[2], "_to_",
                                      tail(names(df), n = 1), "_Phenotypes_",
                                      npcs[1], "_to_", tail(npcs, n = 1),
                                      "_PCs.csv"))
    best_LambdaGC <- get_best_PC_df(df = LambdaGC)
    write_csv(best_LambdaGC, path = paste0("Best_Lambda_GC_", names(df)[2],
                                           "_to_", tail(names(df), n = 1),
                                           "_Phenotypes_", npcs[1], "_to_",
                                           tail(npcs, n = 1), "_PCs.csv"))
  }
  return(LambdaGC)
}

#' Find lambda_GC value for non-NA p-values
#'
#' @description Finds the lambda GC value for some vector of p-values.
#'
#' @param ps Numeric vector of p-values. Can have NA's.
#' @param tol Numeric. Tolerance for optional Genomic Control coefficient.
#'
#' @importFrom stats median uniroot
#'
#' @return A lambda GC value (some positive number, ideally ~1)
#'
#' @export
get_lambdagc <- function(ps, tol = 1e-8){
  ps <- ps[which(!is.na(ps))]
  xtr <- log10(ps)
  MEDIAN <- log10(0.5)
  f.opt <- function(x) (x - MEDIAN)
  xtr_p <- median(xtr) / uniroot(f.opt, interval = range(xtr),
                                 check.conv = TRUE,
                                 tol = tol)$root
  lamGC <- signif(xtr_p)
  return(lamGC)
}


#' Return best number of PCs in terms of lambda_GC
#'
#' @description Given a dataframe created using div_lambda_GC, this function
#'     returns the first lambda_GC less than 1.05, or the smallest lambda_GC,
#'     for each column in the dataframe.
#'
#' @param df Dataframe of phenotypes where the first column is NumPCs and
#'     subsequent column contains lambda_GC values for some phenotype.
#'
#' @importFrom dplyr filter top_n select full_join arrange
#' @importFrom tidyr gather
#' @importFrom rlang .data sym !!
#' @importFrom tidyselect all_of
#'
#' @return A dataframe containing the best lambda_GC value and number of PCs
#'     for each phenotype in the data frame.
get_best_PC_df <- function(df){
  column <- names(df)[ncol(df)]
  bestPCs <- df %>%
    filter(!! sym(column) < 1.05| !! sym(column) == min(!! sym(column))) %>%
    top_n(n = -1, wt = .data$NumPCs) %>%
    select(.data$NumPCs, all_of(column))

  if(ncol(df) > 2){
    for(i in c((ncol(df)-2):1)){
      column <- names(df)[i+1]

      bestPCs <- df %>%
        filter(!! sym(column) < 1.05 | !! sym(column) == min(!! sym(column))) %>%
        top_n(n = -1, wt = .data$NumPCs) %>%
        select(.data$NumPCs, all_of(column)) %>%
        full_join(bestPCs, by = c("NumPCs", (column)))
    }
  }

  bestPCdf <- bestPCs %>%
    arrange(.data$NumPCs) %>%
    gather(key = "trait", value = "lambda_GC", 2:ncol(bestPCs)) %>%
    filter(!is.na(.data$lambda_GC))

  return(bestPCdf)
}


div_mash <- function(){}


check_gwas <- function(df1, phename, type, nPhe, minphe, nLev){
  if(nPhe < minphe){
    message(paste0("The phenotype ", phename, " does not have the minimum ",
                   "number of phenotyped sample.ID's, (", minphe, ") and so ",
                   "will not be used for GWAS."))
    gwas_ok <- FALSE
  } else if(nLev < 2){
    message(paste0("The phenotype ", phename, " does not have two or more ",
                   "distinct non-NA values and will not be used for GWAS."))
    gwas_ok <- FALSE
  } else if(nLev > 2 & type == "logistic"){
    message(paste0("The phenotype ", phename, " has more than two distinct ",
                   "non-NA values and will not be used for GWAS with 'type=",
                   "logistic'."))
    gwas_ok <- FALSE
  } else if(!(unique(df1[which(!is.na(df1[,2])),2])[1,1] %in% c(0,1)) &
            !(unique(df1[which(!is.na(df1[,2])),2])[2,1] %in% c(0,1)) &
            type == "logistic"){
    message(paste0("The phenotype ", phename, " has non-NA values that are ",
                   "not 0 or 1 and will not be used for GWAS with 'type=",
                   "logistic'."))
    gwas_ok <- FALSE
  } else {
    gwas_ok <- TRUE
  }
  return(gwas_ok)
}


## @title Basic sanity check for covariance matrices
## @param X input matrix
check_covmat_basics = function(x) {
  label = substitute(x)
  if (!is.matrix(x))
    labelled_stop(label, "is not a matrix")
  if (!is.numeric(x))
    labelled_stop(label, "is not a numeric matrix")
  if (any(is.na(x)))
    labelled_stop(label, "cannot contain NA values")
  if (any(is.infinite(x)))
    labelled_stop(label, "cannot contain Inf values")
  if (any(is.nan(x)))
    labelled_stop(label, "cannot contain NaN values")
  if (nrow(x) != ncol(x))
    labelled_stop(label, "is not a square matrix")
  if (!isSymmetric(x, check.attributes = FALSE))
    labelled_stop(label, "is not a symmetric matrix")
  return(TRUE)
}

## @title check matrix for positive definitness
## @param X input matrix
check_positive_definite = function(x) {
  check_covmat_basics(x)
  tryCatch(chol(x),
           error = function(e) labelled_stop(substitute(x),
                                             "must be positive definite"))
  return(TRUE)
}
