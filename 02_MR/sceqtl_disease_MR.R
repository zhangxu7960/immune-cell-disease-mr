args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])
print('------------start--------------')

library(data.table)
library(readr)
library(utils)
library(TwoSampleMR)
library(foreach)
library(doParallel)
library(RadialMR)
library(ggplot2)
library(MRPRESSO)
library(dplyr)
library(mr.raps)
library(readxl)


namematch <- readxl::read_xlsx("/public/home/zhangxu02/immune_disease/6640-disease-information.xlsx")

immune_cells_basepath <- "/public/home/zhangxu02/immune_disease/data_fliter_1e5_0.001_10000kb"
immune_cells <- list.dirs(immune_cells_basepath,full.names = F, recursive = F)

exposure_files_list <- unlist(lapply(immune_cells, function(immune_cell) {
  exposure_adress <- file.path(immune_cells_basepath, immune_cell)
  exposure_files <- list.files(exposure_adress, full.names = TRUE)
  lapply(exposure_files, function(full_path) {
    read.table(full_path, sep = "\t",header = T,
               colClasses = c("A1" = "character","A2" = "character"))  # Keep allele codes such as T as character strings.

  })
}), recursive = FALSE)

blank_plot <- function(text) {

  plot <- ggplot() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    geom_text(aes(x = 0.5, y = 0.5, label = text), size = 6, hjust = 0.5, vjust = 0.5)+
    theme_bw()

  return(plot)
}

plotall <- function(file, mr_results, dat) {
  dat <- as.data.frame(dat)
  dat <- dat[dat$mr_keep %in% TRUE, , drop = FALSE]
  rownames(dat) <- NULL

  n <- nrow(dat)
  status <- list()

  # Record skipped or failed sensitivity steps without stopping the other steps.
  step <- function(name, min_n, expr) {
    state <- "ok"
    info <- ""
    value <- NULL

    if (n < min_n) {
      state <- "skipped"
      info <- paste0("Need >= ", min_n, " SNPs")
    } else {
      value <- tryCatch(force(expr), error = function(e) {
        state <<- "error"
        info <<- conditionMessage(e)
        NULL
      })
    }

    status[[length(status) + 1L]] <<- data.frame(
      step = name, status = state, nsnp = n, message = info
    )
    value
  }

  save_table <- function(x, name) {
    data.table::fwrite(x, file.path(file, name))
  }

  save_plot <- function(p, name, width = 8, height = 6) {
    grDevices::pdf(
      file.path(file, name), width = width, height = height
    )
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
  }

  step("pleiotropy", 3, {
    save_table(
      TwoSampleMR::mr_pleiotropy_test(dat), "ple.csv"
    )
  })

  step("heterogeneity", 2, {
    methods <- "mr_ivw"
    if (n >= 3) methods <- c(methods, "mr_egger_regression")

    save_table(
      TwoSampleMR::mr_heterogeneity(
        dat, method_list = methods
      ),
      "H.csv"
    )
  })

  step("mrpresso", 4, {
    nb <- max(10000L, 100L * n + 1L)

    ans <- MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      data = dat,
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      NbDistribution = nb,
      SignifThreshold = 0.05
    )

    saveRDS(ans, file.path(file, "mrpresso.rds"))
    capture.output(ans, file = file.path(file, "mrpresso.txt"))

    save_table(ans[["Main MR results"]], "mrpresso_main.csv")

    tests <- ans[["MR-PRESSO results"]]

    global <- as.data.frame(tests[["Global Test"]])
    global$nsnp <- n
    global$NbDistribution <- nb
    save_table(global, "mrpresso_global.csv")

    out <- tests[["Outlier Test"]]
    if (!is.null(out)) {
      out <- as.data.frame(out)
      idx <- match(rownames(out), rownames(dat))

      if (anyNA(idx)) {
        stop("Cannot match MR-PRESSO rows to input SNPs")
      }

      out <- data.frame(
        SNP = dat$SNP[idx], out, check.names = FALSE
      )
      save_table(out, "mrpresso_outliers.csv")
    }

    distortion <- tests[["Distortion Test"]]
    if (!is.null(distortion)) {
      save_table(
        data.frame(
          Outliers_Indices = paste(
            distortion[["Outliers Indices"]], collapse = ";"
          ),
          Distortion_Coefficient =
            distortion[["Distortion Coefficient"]],
          Pvalue = as.character(distortion[["Pvalue"]])
        ),
        "mrpresso_distortion.csv"
      )
    }
  })

  step("scatter", 2, {
    save_plot(
      TwoSampleMR::mr_scatter_plot(mr_results, dat)[[1]],
      "scatter.pdf"
    )
  })

  single <- step("single_snp", 1, {
    s <- TwoSampleMR::mr_singlesnp(dat)
    save_table(s, "single_snp.csv")
    s
  })

  if (!is.null(single)) {
    step("forest", 2, {
      save_plot(
        TwoSampleMR::mr_forest_plot(single)[[1]],
        "single_forest.pdf"
      )
    })

    step("funnel", 2, {
      save_plot(
        TwoSampleMR::mr_funnel_plot(single)[[1]],
        "funnel.pdf"
      )
    })
  }

  step("leave_one_out", 3, {
    loo <- TwoSampleMR::mr_leaveoneout(dat)
    save_table(loo, "leave_one_out.csv")

    save_plot(
      TwoSampleMR::mr_leaveoneout_plot(loo)[[1]],
      "leave-one-out.pdf"
    )
  })

  radial_ivw <- step("radial_ivw", 3, {
    fit <- RadialMR::ivw_radial(
      r_input = dat, alpha = 0.05, weights = 1,
      tol = 0.0001, summary = FALSE
    )
    saveRDS(fit, file.path(file, "radial_ivw.rds"))
    fit
  })

  radial_egger <- step("radial_egger", 3, {
    fit <- RadialMR::egger_radial(
      r_input = dat, alpha = 0.05,
      weights = 1, summary = FALSE
    )
    saveRDS(fit, file.path(file, "radial_egger.rds"))
    fit
  })

  if (!is.null(radial_ivw) || !is.null(radial_egger)) {
    step("radial_plot", 3, {
      obj <- if (!is.null(radial_ivw) &&
                 !is.null(radial_egger)) {
        c(radial_ivw, radial_egger)
      } else if (!is.null(radial_ivw)) {
        radial_ivw
      } else {
        radial_egger
      }

      save_plot(
        RadialMR::plot_radial(obj), "radial.pdf", 12, 8
      )
    })
  }

  save_table(
    data.table::rbindlist(status), "sensitivity_status.csv"
  )
  invisible(NULL)
}

MR_function = function(exposure_dat, outcome_id, df_batch, trait){

  exposure_name <- as.character(exposure_dat$GENE[1])
  celltype <- as.character(exposure_dat$CELL_ID[1])

  colnames(exposure_dat)[colnames(exposure_dat) == "P_VALUE"] <- "pval"
  colnames(exposure_dat)[colnames(exposure_dat) == "RSID"] <- "SNP"
  colnames(exposure_dat)[colnames(exposure_dat) == "A2_FREQ_ONEK1K"] <- "eaf"

  exposure_dat$F <- (exposure_dat$beta / exposure_dat$se)^2

  exposure_dat <- exposure_dat[
    is.finite(exposure_dat$beta) &
      is.finite(exposure_dat$se) &
      exposure_dat$se > 0 &
      is.finite(exposure_dat$F) &
      exposure_dat$F > 10,
    ,
    drop = FALSE
  ]

  if(nrow(exposure_dat) > 0){
    dat <- df_batch[df_batch$SNP %in% exposure_dat$SNP,]
    dat = as.data.frame(dat)

    if(nrow(dat)>0){
      dat$id = outcome_id
      dat$trait = trait

      outcome_dat <- TwoSampleMR::format_data(
        dat,
        snps = exposure_dat$SNP,
        type = "outcome",
        snp_col = "SNP",
        beta_col = "beta",
        se_col = "se",
        eaf_col = "eaf",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        pval_col = "pval"
      )

      outcome_dat$outcome = trait
      outcome_dat$id.outcome = outcome_id

      exposure_dat <- dplyr::rename(exposure_dat,
                                    effect_allele.exposure = A2,  # A2 is the effect allele.
                                    other_allele.exposure = A1,  # A1 is the other allele.
                                    beta.exposure = beta,
                                    eaf.exposure = eaf,
                                    se.exposure = se,
                                    pval.exposure = pval,
                                    exposure = GENE
      )
      exposure_dat$id.exposure <- exposure_dat$exposure[1]

      dat <- TwoSampleMR::harmonise_data(exposure_dat = exposure_dat,
                                         outcome_dat = outcome_dat
      )

      dat = dat[dat$pval.outcome >= 1e-02,]

      if(nrow(dat)>0){

        result = TwoSampleMR::mr(dat, method_list = c('mr_wald_ratio', 'mr_ivw','mr_ivw_mre','mr_ivw_fe','mr_simple_mode'
                                                      ,'mr_egger_regression','mr_simple_median','mr_raps','mr_weighted_median',
                                                      'mr_weighted_mode'))

        if(nrow(result) > 0){
          dat$outcome = result$outcome[1]
          dat$id.outcome = result$id.outcome[1]
          dat$exposure = result$exposure[1]
          dat$id.exposure = result$id.exposure[1]

          rslt <- file.path("/public/home/zhangxu02/immune_disease/1e5_0.001_10000kb",outcome_id,celltype,exposure_name)
          dir.create(rslt,recursive = TRUE,showWarnings = FALSE)

          result$b <- as.numeric(result$b)
          result$se <- as.numeric(result$se)

          result$lo_ci <- result$b - 1.96 * result$se
          result$up_ci <- result$b + 1.96 * result$se

          result$or <- exp(result$b)
          result$or_lci <- exp(result$lo_ci)
          result$or_uci <- exp(result$up_ci)

          result$cell_type <- celltype

          data.table::fwrite(result, file = paste0(rslt, '/mr_result.csv'),row.names = F)
          data.table::fwrite(dat, file = paste0(rslt, '/dat.csv'),row.names = F)

          tryCatch(
            plotall(rslt, result, dat),
            error = function(e) {
              message(
                "Sensitivity analysis failed: outcome=", outcome_id,
                "; cell_type=", celltype,
                "; gene=", exposure_name,
                "; error=", conditionMessage(e)
              )
            }
          )

          print('-------finish-------')
          return(result)
        }
      }
    }
  }
}

for (i in task_id) {
  print(paste0("-------------------",i,"-------------------"))

  outcome_file <- params[i,]$file
  outcome_id = params[i,]$id
  trait <- namematch$phenotype[which(namematch$GWASID == outcome_id)]
  df_batch <- data.table::fread(outcome_file, select =  c("SNP", "beta", "eaf","se","pval","effect_allele","other_allele"))

  results_list = data.frame()
  results_list <- parallel::mclapply(
    exposure_files_list,
    function(x) {
      tryCatch(
        MR_function(x, outcome_id, df_batch, trait),
        error = function(e) {
          message(
            "MR failed: outcome=", outcome_id,
            "; cell_type=", x$CELL_ID[1],
            "; gene=", x$GENE[1],
            "; call=", paste(deparse(conditionCall(e)), collapse = " "),
            "; error=", conditionMessage(e)
          )
          NULL
        }
      )
    },
    mc.cores = 2,
    mc.preschedule = TRUE
  )

  results_list <- Filter(
    function(x) is.data.frame(x) && nrow(x) > 0,
    results_list
  )

  if (length(results_list) > 0) {
    ivw <- data.table::rbindlist(
      results_list,
      use.names = TRUE,
      fill = TRUE
    )
  } else {

    ivw <- data.table::data.table(
      id.exposure = character(),
      id.outcome = character(),
      outcome = character(),
      exposure = character(),
      method = character(),
      nsnp = integer(),
      b = numeric(),
      se = numeric(),
      pval = numeric(),
      lo_ci = numeric(),
      up_ci = numeric(),
      or = numeric(),
      or_lci = numeric(),
      or_uci = numeric(),
      cell_type = character()
    )

    message("No valid MR results returned: ", outcome_id)
  }

  file_path <- paste0(
    "/public/home/zhangxu02/immune_disease/",
    "1e5ivw_0.001_10000kb/",
    outcome_id, ".csv"
  )

  dir.create(
    dirname(file_path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  data.table::fwrite(ivw, file = file_path, row.names = FALSE)
}
