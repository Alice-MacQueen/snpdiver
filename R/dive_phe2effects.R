#' @title Wrapper to run mash given a phenotype data frame
#'
#' @description This function allows
#'     you to go from a phenotype data.frame of a few phenotypes you want to
#'     compare to filebacked matrix of univariate GWAS effects, standard errors,
#'     and -log10pvalues. This output object can be used in "dive_effects2mash"
#'     function. Some exception handling has been built into
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
#' @param suffix Optional character vector to give saved files a unique search
#'     string/name.
#' @param outputdir Optional file path to save output files.
#' @param min.phe Integer. Minimum number of individuals phenotyped in order to
#'     include that phenotype in GWAS. Default is 200. Use lower values with
#'     caution.
#' @param ncores Optional integer to specify the number of cores to be used
#'     for parallelization. You can specify this with bigparallelr::nb_cores().
#' @param save.plots Logical. Should Manhattan and QQ-plots be generated and
#'     saved to the working directory for univariate GWAS? Default is TRUE.
#' @param thr.r2 Value between 0 and 1. Threshold of r2 measure of linkage
#'     disequilibrium. Markers in higher LD than this will be subset using clumping.
#' @param roll.size Integer. Used to create the svd for GWAS.
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
#' @importFrom bigparallelr nb_cores
#'
#' @export
dive_phe2effects <- function(df, snp, type = "linear", svd = NULL, suffix = "",
                          outputdir = ".",
                          min.phe = 200, ncores = NA,
                          save.plots = TRUE, thr.r2 = 0.2,
                          roll.size = 50, verbose = TRUE){
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
  if(is.na(ncores)){
    ncores <- nb_cores()
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
                       k = 10, thr.r2 = thr.r2, roll.size = roll.size,
                       ncores = ncores)
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
      summarise(phe = mean(.data[[phename]], na.rm = TRUE),
                .groups = "drop_last")
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
                                   svd = svd, npcs = c(0:pc_max),
                                   ncores = ncores)
      PC_df <- get_best_PC_df(lambdagc_df)
      PC_df <- PC_df[1,]

      # 3. GWAS  ----

      # run gwas using best npcs from step 2 (best pop structure correction)
      gwas <- div_gwas(df = df1, snp = snp, type = type[i - 1], svd = svd,
                       npcs = PC_df$NumPCs, ncores = ncores)
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
                                   thresh = bonferroni, ncores = ncores)
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

  printf2(verbose = verbose, "\nNow replacing NA gwas effects with 0's for use in mash.\n")
  # 4. mash input ----
  ## prioritize effects with max(log10p) or max(sum(log10p))
  ## make a random set of relatively unlinked SNPs
  ind_estim <- (1:sum(gwas_ok))*3 - 2
  ind_se <- (1:sum(gwas_ok))*3 - 1
  ind_p <- (1:sum(gwas_ok))*3

  ## replace NA or Nan values
  # Replace SE with 1's, estimates and p values with 0's.
  replace_na_1 <- function(X, ind) replace_na(X[, ind], 1)
  replace_na_0 <- function(X, ind) replace_na(X[, ind], 0)
  gwas2[, ind_se] <- big_apply(gwas2, a.FUN = replace_na_1, ind = ind_se,
                               a.combine = 'plus', ncores = ncores)
  gwas2[, ind_estim] <- big_apply(gwas2, a.FUN = replace_na_0, ind = ind_estim,
                                  a.combine = 'plus', ncores = ncores)
  gwas2[, ind_p] <- big_apply(gwas2, a.FUN = replace_na_0, ind = ind_p,
                              a.combine = 'plus', ncores = ncores)
  gwas2$save()

  ## No scaling in this function.
  gwas_metadata <- gwas_metadata %>% mutate(scaled = FALSE)

  ## Save column names and gwas metadata.
  write_csv(tibble(colnames_fbm), file.path(outputdir,
                                            paste0("gwas_effects", suffix,
                                                   "_column_names.csv")))
  write_csv(gwas_metadata, file.path(outputdir,
                                     paste0("gwas_effects", suffix,
                                            "_associated_metadata.csv")))
  return(list(effects = gwas2, effect_cols = colnames_fbm,
              metadata = gwas_metadata))
}
