
##### House Keeping Functions #####

#' @title Remove NAs
#'
#' @description Function to remove CpGs with NAs from test_result
#' @param test_result DF of summaryExtractTest() from DML summary
#' @return a modified data frame with CpGs containing NAs removed
#' @examples
#' betas <- remove_NAs(betas)
#' @export
remove_NAs <- function(test_result) {
  remove <- c()
  for (i in c(1:ncol(test_result))) {
    if (length(grep("FPval", colnames(test_result)[i])) == 0) {
      t <- which(is.na(test_result[, i]))
      remove <- c(remove, t)
    }
  }

  if (length(remove) > 0) {
    test_result <- test_result[-remove, ]
  }
  print(paste(length(remove), "NA CpGs removed from test result."))
  print(paste("Final test_result dimensions (", paste(dim(test_result), collapse = ","), ")", sep = ""))

  return(test_result)
}



#' @title Matrix Variance Calculation
#'
#' @description Function to calculate the z-scores for each row of a matrix
#' @param x A matrix where rows are measurements and columns are replicates
#' @return A modified matrix of same input dimensions but z-score values calculated for each row.
#' @examples
#' betas_zscore <- MatVar(betas)
#' @export
MatVar <- function (x) {
  a <- rowSds(x)
  b <- rowMeans(x)
  v <- (x - b)/a
  return(v)
}

###### Unit Tests #####

#' @title Test precense of a input condition
#'
#' @description For each of the conditions supplied by the formula input, this function will check if there is a colname present in the provided sample sheet that matches for further analysis.
#' @param condition A character string of one of the variables from formula input
#' @param sample_sheet_df Data frame of read-in sample sheet provided for pipeline
#' @return If condition column not found in sample sheet, will stop the pipeline and report error
#' @examples
#' test_condition_presence("Group", sample_sheet)
#' @export
test_condition_presence <- function(condition, sample_sheet_df) {
  l <- length(grep(condition, colnames(sample_sheet_df)))
  if (l == 0) {
    stop(paste("The following condition is not found in the sample sheet:", condition))
  }
}


#' @title Tests Output Architecture
#'
#' @description Tests for output directory architecture to make sure that pipeline can output the necessary documents/plots
#' @param out_dir Path to the output directory
#' @return If architecture incorrect, will stop the pipeline and report error
#' @examples
#' test_out_architecture("./output")
#' @export
test_output_architecture <- function(out_dir) {

  ## first level subdirs
  dir_list <- list.dirs(out_dir, full.names = F, recursive = F)
  fl_dirs <- c("QC", "PCA", "DML")
  for (i in fl_dirs) {
    if (!(i %in% dir_list)) {
      stop(paste("You are missing the following directory in the output directory:", i))
    }
  }

  ## second level subdirs
  DML_dirs <-  list.dirs(paste(out_dir, "/DML", sep = ""), full.names = F, recursive = F)
  sl_dirs <- c("DMR_Analysis", "GO_Enrichment", "Heatmaps", "testEnrichments", "Volcano_plots")
  for (j in sl_dirs) {
    if (!(j %in% DML_dirs)) {
      stop(paste("You are missing the following directory in the output directory: /DML/", j, sep = ""))
    }
  }
}


##### Read in data and preprocessing #####

#' @title Read in Sample Sheet
#'
#' @description Function to read in the provided SeSAMe STREET sample sheet in the correct formatting
#' @param sample_sheet Path to the sample sheet
#' @return Dataframe of the read in sample sheet data
#' @examples
#' sample_sheet_df <- SS_samplesheet("SeSAMe_STREET_Sample_Sheet.xlsx")
#' @export
SS_samplesheet <- function(sample_sheet) {
  sheet <- readxl::read_xlsx(sample_sheet, skip = 14)[, -c(1:5)]
  return(sheet)
}

#' @title Calculate Beta Values
#'
#' @description Function to read in Idat files and return the calculated beta values
#' @param Idat_dir Path to the methylation array data files
#' @param prep Character string of the preprocessing code
#' @return Data frame of beta values
#' @examples
#' betas <- get_betas("./Idat_dir", "TQCDPB")
#' @export
get_betas <- function(Idat_dir, prep) {
  betas <- openSesame(Idat_dir, func = getBetas, prep = prep)
  betas <- remove_NAs(betas)
  print(paste("Betas DF Dims:", dim(betas)))
  return(betas)
}

#' @title Retreive raw intensities
#'
#' @description Function to read in Idats and return intensity values
#' @param Idat_dir Path to the methylation array data files
#' @param prep Character string of the preprocessing code
#' @return list of data frames for each replicates reporting their intensity values
#' @examples
#' sdfs <- get_sdfs("./Idat_dir", "TQCDPB")
#' @export
get_sdfs <- function(Idat_dir, prep) {
  sdfs <- openSesame(Idat_dir, func = NULL, prep = prep)
  return(sdfs)
}



##### Quality control ####

#' @title Quality Metrics
#'
#' @description Function to calculate and write quality metics & ranked metrics
#' @param sdfs list of data frames of intensity values for each replicate; output of get_sdfs
#' @param prep Character string of the preprocessing code
#' @param out_dir Path to output directory
#' @return Outputs quality metrics to /QC subdir in out_dir
#' @examples
#' sdfs <- get_sdfs("./Idat_dir", "TQCDPB")
#' qual_metrics(sdfs, "TQCDPB", "output")
#' @export
qual_metrics <- function(sdfs, prep, out_dir) {

  message("Calculating quality metrics...")

  ## quality metrics
  qcs <- openSesame(sdfs, prep, func=sesameQC_calcStats, funs="detection")

  ## this crashes so I removed the output for now, will fix this later
  #qc_df <- head(do.call(rbind, lapply(qcs, as.data.frame)))
  out_file = paste(out_dir, "/QC/quality_metics.csv", sep = "")
  #write_csv(x = qc_df, file = out_file)


  ## ranked quality metrics
  ranked_qcs <- lapply(X = c(1:length(sdfs)), FUN = function(X){qc <- sesameQC_calcStats(sdfs[[X]], "intensity")
  sesameQC_rankStats(qc, platform="MM285")})

  out_file = paste(out_dir, "/QC/ranked_quality_metrics.csv", sep = "")
  #write_csv(x = ranked_qcs, file = out_file)
  #Problem with writing ranks: Error in as.data.frame.default(x[[i]], optional = TRUE) :
  #  cannot coerce class ‘structure("sesameQC", package = "sesame")’ to a data.frame
}


#' @title Plot Detection Rate
#'
#' @description Function to plot detection rate QC
#' @param sdfs list of data frames of intensity values for each replicate; output of get_sdfs
#' @param out_dir Path to output directory
#' @return Saves plot to /QC subdir
#' @examples
#' sdfs <- get_sdfs("./Idat_dir", "TQCDPB")
#' detection_rate(sdfs, "output)
#' @export
detection_rate <- function(sdfs, out_dir) {

  message("Calculating detection rates...")

  out_file <- paste(out_dir, "/QC/detection_rate.pdf", sep = "")

  pdf(out_file)
  plot <- sesameQC_plotBar(lapply(sdfs, sesameQC_calcStats, "detection"))
  print(plot)
  dev.off()

}


#' @title Green-Red QQ Plots
#'
#' @description Function to generate the Green-Red signal QQ plots for each replicate
#' @param sdfs list of data frames of intensity values for each replicate; output of get_sdfs
#' @param out_dir Path to output directory
#' @return Saves QQ plots to /QC subdir
#' @examples
#' sdfs <- get_sdfs("./Idat_dir", "TQCDPB")
#' GR_QQ_plots(sdfs, "output)
#' @export
GR_QQ_plots <- function(sdfs, out_dir) {

  message("Generating QQ Plots...")

  out_file <- paste(out_dir, "/QC/GR_QQ_plots.pdf", sep = "")
  pdf(out_file)
  plot_list <- lapply(X = c(1:length(sdfs)), FUN = function(X){sesameQC_plotRedGrnQQ(sdfs[[X]]) + title(names(sdfs)[X], adj = 0)})
  dev.off()

}

#' @title Beta Intensity Plots
#'
#' @description Function to plot the beta value intensity for each replicate
#' @param sdfs list of data frames of intensity values for each replicate; output of get_sdfs
#' @param out_dir Path to output directory
#' @return Saves intensity plots in /QC subdir
#' @examples
#' sdfs <- get_sdfs("./Idat_dir", "TQCDPB")
#' intensity_beta_plots(sdfs, "output)
#' @export
intensity_beta_plots <- function(sdfs, out_dir) {

  message("Generating intensity plots...")

  out_file <- paste(out_dir, "/QC/Intensity_plots.pdf", sep = "")
  pdf(out_file)
  plot_list <- lapply(X = c(1:length(sdfs)), FUN = function(X){sesameQC_plotIntensVsBetas(sdfs[[X]]) + title(names(sdfs)[X], adj = 0)})
  dev.off()

}



#' @title Extra SNP Freq Plot
#'
#' @description Function to plot the extra SNP allele frequencies
#' @param sdfs list of data frames of intensity values for each replicate; output of get_sdfs
#' @param out_dir Path to output directory
#' @return Saves extra SNP freq plots in /QC subdir
#' @examples
#' sdfs <- get_sdfs("./Idat_dir", "TQCDPB")
#' extra_SNP_freq(sdfs, "output)
#' @export
extra_SNP_freq <- function(sdfs, out_dir) {

  message("Calculating extra SNP freq...")

  out_file <- paste(out_dir, "/QC/extra_SNP_allele_freq.pdf", sep = "")
  pdf(file = out_file)


  afs <- openSesame(sdfs, func = getAFs, mask = FALSE)
  afs <- both.cluster(afs)$mat
  rg <- apply(afs, 1, function(x) {
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  })
  afs <- afs[rg > 0.3, ]

  if (nrow(afs) > 0) {
    plot <- sesameQC_plotHeatSNPs(sdfs)
    print(plot)
  }else{
    warning("nrow(afs) == 0; skipping extra SNP freq heatmap")
  }
  dev.off()
}

#' @title Model Infered Species
#'
#' @description Function to find the infered species for all replicates of a sesame run SDFobject
#' @param sdfs list of data frames of intensity values for each replicate; output of get_sdfs
#' @param out_dir Path to output directory
#' @return Writes csv of infered species to /QC subdir
#' @examples
#' sdfs <- get_sdfs("./Idat_dir", "TQCDPB")
#' all_inferedspecies(sdfs, "output)
#' @export
all_inferedspecies <- function(sdfs, out_dir) {
  message("Running infered species...")
  infspec <- unlist(lapply(c(1:length(sdfs)), function(X){inferSpecies(sdfs[[X]], return.species = TRUE)}))

  infspec_list <- list()
  for (i in unique(names(infspec))) {
    v <- unique(infspec[which(names(infspec) == i)])
    temp_list <- list(i = v)
    infspec_list <- c(infspec_list, temp_list)
  }
  names(infspec_list) <- unique(names(infspec))

  out_file <- paste(out_dir, "/QC/Infered_Species.csv", sep = "")
  arw <- paste(names(infspec_list), (infspec_list), sep = ": ")
  write_csv(as.data.frame(arw), file = out_file)
}

#' Model Infered Strains
#'
#' @description Plotting the infered mouse strain for each replicate
#' @param sdfs list of data frames of intensity values for each replicate; output of get_sdfs
#' @param out_dir Path to output directory
#' @return Plots probabilities of infered strains to /QC subdir
#' @examples
#' sdfs <- get_sdfs("./Idat_dir", "TQCDPB")
#' plot_inferedstrains(sdfs, "output)
#' @export
plot_inferedstrains <- function(sdf, out_dir) {
  message("Running infered strains...")
  infstr <- data.frame()
  for (i in c(1:length(sdf))) {
    p = inferStrain(sdf[[i]], return.probability = TRUE)
    df = data.frame(strain=names(p), probs=p, replicate = names(sdf)[i])
    infstr <- rbind(infstr, df)
  }

  plot <- ggplot(data = infstr,  aes(x = strain, y = probs, fill = replicate)) +
    geom_bar(stat = "identity", position = "dodge") +
    ggtitle("Strain Probabilities") +
    ylab("Probability") + xlab("") +
    scale_x_discrete(position = "top") +
    scale_fill_viridis_d() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0))

  out_file = paste(out_dir, "/QC/Infered_strain_prob.pdf", sep = "")

  ggsave(out_file, plot, device = "pdf")
}

#' @title QC Subassembly
#'
#' @description Function to run all QC functions
#' @param Idat_dir Path to the methylation array data files
#' @param prep Code for preprocessing steps from SeSAMe Street run input argument
#' @param out_dir Path to output direcrtory
#' @return Files and plots saved to /QC subdir
#' @examples
#' run_QC("path/Idat_dir", "TQCDPB", "path/out_dir")
#' @export
run_QC <- function(Idat_dir, prep, out_dir) {

  sdfs <- get_sdfs(Idat_dir, prep)
  out_file <- paste(out_dir, "/QC/sdfs.RData", sep = "")
  save(sdfs, file = out_file)

  qual_metrics(sdfs, prep, out_dir)
  detection_rate(sdfs, out_dir)
  GR_QQ_plots(sdfs, out_dir)
  intensity_beta_plots(sdfs, out_dir)
  extra_SNP_freq(sdfs, out_dir)

  ## add conditional stratement for if mouse species is selected to run or not. Version 1 of pipeline will only be compatable with MM285
  all_inferedspecies(sdfs, out_dir)
  plot_inferedstrains(sdfs, out_dir)
}



##### PCA #####

#' @title PC Variance Plot
#'
#' @description Function to plot the PC proportion of explained variance
#' @param pca PCA reults of prcomp from preprocessed beta values
#' @param out_dir Path to output directory
#' @return GGplit2 formatted Plot of the explained variance saved to /PCA subdir
#' @examples
#' PCA_variance_plot(pca, "path/out_dir")
#' @export
PCA_variance_plot <- function(pca, out_dir) {

  prop_var <- data.frame(t(summary(pca)$importance))
  names(prop_var) = c('sd', 'prop', 'cum')
  prop_var$num = 1:nrow(prop_var)

  out_file <- paste(out_dir, "/DimRed/Component_Variance.pdf", sep = "")

  pdf(file = out_file)

  var_plot <- ggplot(prop_var, aes(x=num, y=prop)) +
    geom_point(size=1.5) +
    geom_line() +
    scale_x_continuous(limits = c(1, 12), breaks = 1:12) +
    xlab("Principal Component") +
    ylab("Prop. of Variance") +
    ggtitle("") +
    theme_bw(10) +
    theme(
      axis.title.y = element_text(vjust=1),
      plot.margin = unit(c(0,0,0,6), "mm")
    )
  print(var_plot)

  dev.off()

  #ggsave(filename = out_file, var_plot, device = NULL)
}

#' @title Plot PC1 & PC2
#'
#' @description Function to plot the first two PCs and color by some condition value
#' @param pca PCA reults of prcomp from preprocessed beta values
#' @param condition Character value of a condition name which can be linked to a sample sheet column for coloring of the replicates by the levels of that condition
#' @param sample_sheet_df Dataframe of read in sample sheet
#' @return GGplot2 formatted dotplot of PC1 and 2 colored by condition variables
#' @examples
#' plot_PCA(pca, "Group", sample_sheet_df)
#' @export
plot_PCA <- function(pca, condition, sample_sheet_df) {

  pr_comps <- data.frame(pca$x)
  pr_comps$sample <- rownames(pr_comps)

  ss_n <- unlist(sample_sheet_df[,grep("Sample_Name", colnames(sample_sheet_df))[1]])
  cond <- sample_sheet_df[match(rownames(pr_comps), ss_n), which(colnames(sample_sheet_df) == condition)]
  pr_comps$Condition <- unlist(cond)
  #pr_comps$Condition <- sample_sheet_df[match(rownames(pr_comps), ss_n), which(colnames(sample_sheet_df) == condition)]
  pr_comps$Condition <- factor(pr_comps$Condition)

  plot <- ggplot(pr_comps, aes(x=PC1, y=PC2, color = Condition)) +
    geom_point(size=2)  +
    #theme_bw(12) +
    #scale_color_brewer(palette = "Set1") +
    #geom_text_repel(data = pca_for_plotting, aes(x = PC1, y = PC2, label = sample), max.overlaps = Inf) +
    theme_classic() +
    theme(legend.position = "right", aspect.ratio = 1, text = element_text(face = "bold", size = 15),
          axis.line = element_line(size = 1), plot.title = element_text(hjust = .5))
  return(plot)

}


#' @title Calulate tSNE data frame
#'
#' @description Function to calculate tSNE1 and 2 for the beta values by sample
#' @param betas Beta values matrix
#' @param perplexity Set to default of Rtsne; 30, but a test is set in place in case this is too high
#' @return data frame of tSNE1 and 2 with rownames set to sample names
#' @examples
#' calc_tSNE(betas)
#' @export
calc_tSNE <- function(betas, perplexity = 30) {

  ppl <- floor((nrow(t(betas)) - 1) / 3)
  if (perplexity > ppl) {
    message(paste("Perplexity is set too high, adjusting to highest allowed; ", ppl, sep = ""))
    perplexity = ppl
  }

  set.seed(1)
  tSNE_dims <- Rtsne(t(betas), perplexity = perplexity)
  tSNE_df <- as.data.frame(tSNE_dims$Y)
  rownames(tSNE_df) <- colnames(betas)

  return(tSNE_df)
}


#' @title Plot tSNE
#'
#' @description Function to plot tSNE1 and 2 from calcualted tSNE data frame
#' @param tSNE_df Data frame of tSNE1 and 2 with sample names as rownames (output of calc_tSNE)
#' @param condition Character value of a condition name which can be linked to a sample sheet column for coloring of the replicates by the levels of that condition
#' @param sample_sheet_df Dataframe of read in sample sheet
#' @return GGplot2 formatted dotplot of tSNE1 and 2 colored by condition variables
#' @examples
#' plot_tSNE(tSNE_df, "Group", sample_sheet_df)
#' @export
plot_tSNE <- function(tSNE_df, condition, sample_sheet_df) {

  #if (!(condition %in% colnames(sample_sheet_df))) {
  #  stop(paste("condition", condition, "not found in colnames(sample_sheet_df)"))
  #}
  tSNE_df$condition <- sample_sheet_df[[condition]][match(rownames(tSNE_df), sample_sheet_df$Sample_Name)]

  tSNE_plot <- ggplot(tSNE_df, aes(x=V1,y=V2,color=condition))+
    geom_point()+
    ylab("tSNE2") + xlab("tSNE1")+
    theme_classic() +
    theme(legend.position = "right", aspect.ratio = 1, text = element_text(face = "bold", size = 15),
          axis.line = element_line(size = 1), plot.title = element_text(hjust = .5))

  return(tSNE_plot)
}

#' @title Calculate UMAP
#'
#' @description Function to calculate UMASP embeddings of samples from betas matrix
#' @param betas matrix of beta values
#' @return Data frame of UMAP1 and 2
#' @examples
#' calc_UMAP(betas)
#' @export
calc_UMAP <- function(betas) {
  set.seed(1)
  nn <- ifelse(ncol(betas) < 15, ncol(betas)-1, 15)
  if (nn < 15) {
    message(paste("UMAP n_neighbors too high, adjusting to ", nn))
  }
  UMAP <- umap(t(betas), n_neighbors = nn)
  UMAP_df <- as.data.frame(UMAP)

  return(UMAP_df)
}

#' @title Plot UMAP
#'
#' @description Function to plot UMAP1 and 2 from calcualted tSNE data frame
#' @param UMAP_df Data frame of UMAP1 and 2 with sample names as rownames (output of calc_UMAP)
#' @param condition Character value of a condition name which can be linked to a sample sheet column for coloring of the replicates by the levels of that condition
#' @param sample_sheet_df Dataframe of read in sample sheet
#' @return GGplot2 formatted dotplot of UMAP1 and 2 colored by condition variables
#' @examples
#' plot_UMAP(UMAP_df, "Group", sample_sheet_df)
#' @export
plot_UMAP <- function(UMAP_df, condition, sample_sheet_df) {
  UMAP_df$condition <- sample_sheet_df[[condition]][match(rownames(UMAP_df), sample_sheet_df$Sample_Name)]

  UMAP_plot <- ggplot(UMAP_df, aes(x=V1,y=V2,color=condition))+
    geom_point()+
    ylab("UMAP2") + xlab("UMAP1")+
    theme_classic() +
    theme(legend.position = "right", aspect.ratio = 1, text = element_text(face = "bold", size = 15),
          axis.line = element_line(size = 1), plot.title = element_text(hjust = .5))

  return(UMAP_plot)
}







#' @title Dimensional Reduction Subassembly
#'
#' @description A function to run all PCA analysis
#' @param betas Preprocessed beta values
#' @param sample_sheet_df Data frame of read in sample sheet
#' @param formula Formula that will be used for DML, but is also used to select columns to color PCA for
#' @param out_dir Path to output directory
#' @param perplexity Value to use for tSNE perplexity
#' @return Files and plots saved in the /PCA subdir
#' @examples
#' run_PCA(betas, sample_sheet_df, ~ Group, "path/out_dir")
#' @export
run_DimRed <- function(betas, sample_sheet_df, formula, out_dir, perplexity) {

  pca <- prcomp(t(betas))

  PCA_variance_plot(pca, out_dir)

  out_file <- paste(out_dir, "/DimRed/tSNE_df.RData", sep = "")
  tSNE_df <- calc_tSNE(betas, perplexity = perplexity)
  save(tSNE_df, file = out_file)

  out_file <- paste(out_dir, "/DimRed/UMAP_df.RData", sep = "")
  UMAP_df <- calc_UMAP(betas)
  save(UMAP_df, file = out_file)

  message("Plotting DimReds by Sample Sheet Condition...")
  #conds <- strsplit(as.character(formula)[2], " + ")[[1]]
  conds <- all.vars(formula[[2]])
  for (i in conds) {

    ## unit test for condition presence
    test_condition_presence(i, sample_sheet_df)

    #pca plot
    out_file <- paste(out_dir, "/DimRed/PCA_", i, ".pdf", sep = "")
    pdf(out_file)
    plot <- plot_PCA(pca, condition = i, sample_sheet_df)
    print(plot)
    dev.off()

    #tSNE plot
    out_file <- paste(out_dir, "/DimRed/tSNE_", i, ".pdf", sep = "")
    pdf(out_file)
    plot <- plot_tSNE(tSNE_df, condition = i, sample_sheet_df)
    print(plot)
    dev.off()

    #UMAP PLOT
    out_file <- paste(out_dir, "/DimRed/UMAP_", i, ".pdf", sep = "")
    pdf(out_file)
    plot <- plot_UMAP(UMAP_df, condition = i, sample_sheet_df)
    print(plot)
    dev.off()
    #ggsave(filename = out_file, plot = plot, device = NULL)
  }

  ## saving PCA object
  message("Saving PCA object as RData...")
  out_file <- paste(out_dir, "/DimRed/PCA.RData", sep = "")
  save(pca, file = out_file)

}


##### DML #####

#' @title Create Summarized Experiment Object
#'
#' @description Function to properly create and format summarized experiement for the DML calculation
#' @param betas Preprocessed beta values
#' @param sample_sheet_df Data frame of read in sample sheet
#' @param formula Formula that will be used for DML, but is also used to select columns to color PCA for
#' @return Summarized Experiemnt object to be used for input of DML run
#' @examples
#' create_SE(betas, sample_sheet_df, ~ Group)
#' @export
create_SE <-function(betas, sample_sheet_df, formula) {

  ## setting formula columns as factors
  f <- all.vars(formula[[2]])
  for (i in f) {
    cn <- grep(i, colnames(sample_sheet_df))[1]
    sample_sheet_df[, cn] <- factor(unlist(sample_sheet_df[, cn]))
  }

  se <- SummarizedExperiment(assays = list(betas = betas),
                             colData = sample_sheet_df[match(colnames(betas), unlist(sample_sheet_df[, 1])), ])

  cd = as.data.frame(colData(se)); rownames(cd) = NULL
  #cd

  ## Removing NAs for differential methyltaion by levels of conditions
  f <- all.vars(formula[[2]])
  for (i in f) {
    ## unit test for condition presence
    test_condition_presence(i, sample_sheet_df)
    se_ok = (checkLevels(assay(se), colData(se)[[i]]))
    se = se[se_ok,]
  }

  #se_ok = (checkLevels(assay(se), colData(se)$Condition))
  #se = se[se_ok,]

  return(se)

}



#' @title Run DML
#'
#' @description Function to run the DML and summary stats calculation
#' @param betas Preprocessed beta values
#' @param sample_sheet_df Data frame of read in sample sheet
#' @param formula Formula that will be used for DML, but is also used to select columns to color PCA for
#' @param cores Number of cores to use for parallelization of DML run
#' @param out_dir Path to output directory
#' @return Summary statistics object to be used for downstream analysis
#' @examples
#' run_DML(betas, sample_sheet_df, ~ Group, 8, "path/out_dir")
#' @export
run_DML <- function(betas, sample_sheet_df, formula, cores, out_dir) {

  se <- create_SE(betas, sample_sheet_df, formula)
  message("Running DML...")
  smry <- DML(se, formula, BPPARAM = MulticoreParam(workers = cores))

  out_file <- paste(out_dir, "/DML/smry.RData", sep = "")
  save(smry, file = out_file)

  return(smry)

}



##### DML analysis/figures #####

#' @title Condition-Level DMLs
#'
#' @description Function to find the DMLs for a specific comparison of condition levels. To be used for plot generation for this specific level of a condition
#' @param test_result Test results of the DML summary statistics
#' @param CONDITION Comparisons to be me made from the colnames of the sample sheet
#' @param LEVEL Which level comparison of the condition DML are we looking at
#' @param pval Threshold for pvalue cutoff
#' @param effSize Threshold for effect size cutoff
#' @return Updated test_result with colum specific for CONDITION and LEVEL that shows Up, Down, and non-sig DMLs
#' @examples
#' condLEVEL_DMLs(test_result, "Group", "Treatment", pval = 0.05, effSize = 0.05)
#' @export
condLEVEL_DMLs <- function(test_result, CONDITION, LEVEL, pval = 0.05, effSize = 0.05) {

  pval_name <- paste("Pval_", CONDITION, LEVEL, sep = "")
  pval_col <- which(colnames(test_result) == pval_name)

  est_name <- paste("Est_", CONDITION, LEVEL, sep = "")
  est_col <- which(colnames(test_result) == est_name)

  eff_name <- paste("Eff_", CONDITION, sep = "")
  eff_col <- which(colnames(test_result) == eff_name)

  if (length(pval_col) == 0 | length(est_col) == 0) {
    test_result <- c(CONDITION, LEVEL)
  }else{
    DML_colname <- paste("DML", CONDITION, LEVEL, sep = "_")

    test_result[, DML_colname] <- ifelse(test_result[, pval_col] <= pval,
                                         ifelse(test_result[, eff_col] > effSize,
                                                ifelse(test_result[, est_col] > 0, "Up", "Down"),
                                                "Non Sig"),
                                         "Non Sig")
    #if (length(unique(test_result[[DML_colname]])) == 3) {
    #  test_result[[DML_colname]] <- factor(test_result[[DML_colname]], levels = c("Non Sig", "Up", "Down"))
    #}else{
    #  test_result[[DML_colname]] <- test_result[[DML_colname]]
    #}
    test_result[[DML_colname]] <- factor(test_result[[DML_colname]], levels = c("Non Sig", "Up", "Down"))
  }

  return(test_result)

}



#' @title Volcano Plot
#'
#' @description Function to plot CpG volcano plot
#' @param test_result Test results of the DML summary statistics
#' @param CONDITION Comparisons to be me made from the colnames of the sample sheet
#' @param LEVEL Which level comparison of the condition DML are we looking at
#' @param out_dir Path to output directory
#' @param base_condition_level Name of level used to compare DML for condition
#' @return Ggplot2 formatted volcano plot
#' @examples
#' volcano_plot(test_result, "Group", "Treatment", "path/out_dir")
#' @export
volcano_plot <- function(test_result, CONDITION, LEVEL, out_dir, base_condition_level) {

  x <- paste("Est_", CONDITION, LEVEL, sep = "")
  y <- paste("Pval_", CONDITION, LEVEL, sep = "")
  DML_colname <- paste("DML", CONDITION, LEVEL, sep = "_")

  plot <- ggplot(test_result %>% arrange(DML_colname), aes_string(x =x, y = sprintf("-log10(%s)", y))) + #y = -log10(y))) +
    geom_point(aes_string(color = DML_colname)) +
    theme_classic() +
    ggtitle(paste(CONDITION, "; ", LEVEL, " vs ", base_condition_level, sep = "")) +
    scale_color_manual(values = c("grey", "red", "blue")) +
    labs(color = "") +
    theme(legend.position = "right", aspect.ratio = 1, text = element_text(face = "bold", size = 15),
          axis.line = element_line(size = 1), plot.title = element_text(hjust = .5))


  out_file <- paste(out_dir, "/DML/Volcano_plots/", CONDITION, "_", LEVEL, "_vs_", base_condition_level,"_volcano.pdf", sep = "")
  pdf(file = out_file)
  print(plot)
  dev.off()

}

#' @title DML Heatmap
#'
#' @description Function to plot heatmap of DMLs frrom specific CONDITION-LEVEL comparison
#' @param betas Preprocessed beta values
#' @param sample_sheet_df Data frame of read in sample sheet
#' @param test_result Test results of DML summary stats
#' @param CONDITION Metavalue realted to colname of sample sheet
#' @param LEVEL Which comparison of condition levels to use DML results from
#' @param out_dir Path to output directory
#' @param base_condition_level Name of level used to compare DML for condition
#' @return Pheatmap formatted heatmap of DML beta values across replicates colored by selected conditions.
#' @examples
#' heatmap_DMLs(betas, sample_sheet_df, test_result, "Group", "Treatment", "path/out_dir")
#' @export
heatmap_DMLs <- function(betas, sample_sheet_df, test_result, CONDITION, LEVEL, out_dir, base_condition_level) {

  DML_colname <- paste("DML", CONDITION, LEVEL, sep = "_")
  DML_CpGs <- test_result$Probe_ID[which(test_result[[DML_colname]] %in% c("Up", "Down"))]

  if (length(DML_CpGs) == 0) {
    message(paste("The following comparison has 0 DMLs; ", DML_colname, sep = ""))
  }else{
    DML_betas <- betas[which(rownames(betas) %in% DML_CpGs), ]

    if (nrow(DML_betas) > 65500) {
      DML_betas <- DML_betas[sample(1:nrow(DML_betas), 65500),]
      message(">65500 DMLs, reducing number to 65500 DMLs for heatmap plotting")
    }

    ## custom heatmap
    hm_colors <- colorRampPalette(colors = c("dark blue","blue", "cyan", "green","yellow", "red", "dark red"))(1000)

    match_ss <- match(colnames(DML_betas), unlist(sample_sheet_df[, 2]))
    mat_row <- data.frame(Condition = factor(unlist(sample_sheet_df[, grep(CONDITION, colnames(sample_sheet_df))[1]])[match_ss]))
    rownames(mat_row) <- unlist(sample_sheet_df[,2])[match_ss]

    lc <- levels(mat_row$Condition)
    if (length(lc) <= 8) {
      row_colors <-  RColorBrewer::brewer.pal(length(lc), "Dark2")
    }else{
      row_colors <- rainbow(length(lc))
    }

    cond_array <- c()
    for (i in c(1:length(lc))) {
      ta <- row_colors[i]
      names(ta) <- lc[i]
      cond_array <- c(cond_array, ta)
    }
    ann_cols <- list(Condition = cond_array)



    ## raw betas
    out_file <- paste(out_dir, "/DML/Heatmaps/", CONDITION, "_", LEVEL, "_vs_", base_condition_level,"_heatmap.pdf", sep = "")
    pdf(out_file)
    plot <- pheatmap(t(DML_betas), color = hm_colors, annotation_row = mat_row, annotation_colors = ann_cols,
                     cluster_rows = T, cluster_cols = T, show_colnames = F, show_rownames = T, #cellwidth = .15, cellheight = 10,
                     annotation_names_row = F, fontsize = 5)
    print(plot)
    dev.off()


    ## z score betas
    out_file <- paste(out_dir, "/DML/Heatmaps/", CONDITION, "_", LEVEL, "_vs_", base_condition_level,"_zscore_heatmap.pdf", sep = "")
    pdf(out_file)
    is.na(DML_betas) %>% table()
    plot <- pheatmap(t(MatVar(DML_betas)), annotation_row = mat_row, annotation_colors = ann_cols,
                     cluster_rows = T, cluster_cols = T, show_colnames = F, show_rownames = T, #cellwidth = .15, cellheight = 10,
                     annotation_names_row = F, fontsize = 5)
    print(plot)
    dev.off()
  }
}


#' @title Plot Test Enrichments
#'
#' @description Function to calculate enrichment for different SeSAMe CpG databases from DML CpGs
#' @param test_result Test results of DML summary stats
#' @param CONDITION Metavalue realted to colname of sample sheet
#' @param LEVEL Which comparison of condition levels to use DML results from
#' @param out_dir Path to output directory
#' @param base_condition_level Name of level used to compare DML for condition
#' @return Pdf of test enrichment results for CONDITION-LEVEL DML CpGs
#' @examples
#' testEnrichment_DMLs(test_result, "Group", "Treatment", "path/out_dir")
#' @export
testEnrichment_DMLs <- function(test_result, CONDITION, LEVEL, out_dir, base_condition_level) {

  ## testing all CpGs together
  DML_colname <- paste("DML", CONDITION, LEVEL, sep = "_")
  DML_CpGs <- test_result$Probe_ID[which(test_result[[DML_colname]] %in% c("Up", "Down"))]

  tE_results <- testEnrichment(DML_CpGs)

  out_file <- paste(out_dir, "/DML/testEnrichments/", CONDITION, "_", LEVEL,"_vs_", base_condition_level, "_testEnrichment_All.pdf", sep = "")
  pdf(out_file)
  plot <- KYCG_plotEnrichAll(tE_results)
  print(plot)
  dev.off()

  ## testing only the hypermethylated CpGs
  DML_up <- test_result$Probe_ID[which(test_result[[DML_colname]] == "Up")]
  tE_up <- testEnrichment(DML_up)

  out_file <- paste(out_dir, "/DML/testEnrichments/", CONDITION, "_", LEVEL,"_vs_", base_condition_level, "_testEnrichment_Up.pdf", sep = "")
  pdf(out_file)
  plot <- KYCG_plotEnrichAll(tE_up)
  print(plot)
  dev.off()

  ## testing only the hypermethylated CpGs
  DML_down <- test_result$Probe_ID[which(test_result[[DML_colname]] == "Down")]
  tE_down <- testEnrichment(DML_down)

  out_file <- paste(out_dir, "/DML/testEnrichments/", CONDITION, "_", LEVEL,"_vs_", base_condition_level, "_testEnrichment_Down.pdf", sep = "")
  pdf(out_file)
  plot <- KYCG_plotEnrichAll(tE_down)
  print(plot)
  dev.off()


}




#' @title Plot Test Enrichments for a custom set of local CpG databases
#'
#' @description Function to calculate enrichment for different SeSAMe CpG .RData databases from DML CpGs
#' @param custom_set_paths Path to directory of .RData file(s) containing custom CpG sets
#' @param test_result Test results of DML summary stats
#' @param CONDITION Metavalue realted to colname of sample sheet
#' @param LEVEL Which comparison of condition levels to use DML results from
#' @param out_dir Path to output directory
#' @param base_condition_level Name of level used to compare DML for condition
#' @return Pdf of test enrichment results for CONDITION-LEVEL DML CpGs
#' @examples
#' testEnrichment_DMLs(test_result, "Group", "Treatment", "path/out_dir")
#' @export
test_enrichment_custom_sets <- function(custom_set_paths, test_result, CONDITION, LEVEL, out_dir, base_condition_level) {

  for (file in custom_set_paths) {
    load(file, verbose = T)
    db_name = load(file)

    for (db in db_name) {
      ## testing all CpGs together
      DML_colname <- paste("DML", CONDITION, LEVEL, sep = "_")
      DML_CpGs <- test_result$Probe_ID[which(test_result[[DML_colname]] %in% c("Up", "Down"))]

      tE_results <- testEnrichment(DML_CpGs, databases = get(db))

      out_file <- paste(out_dir, "/DML/testEnrichments/custom_sets/",db, "_", CONDITION, "_", LEVEL,"_vs_", base_condition_level, "_testEnrichment_All.pdf", sep = "")
      pdf(out_file)
      plot <- KYCG_plotEnrichAll(tE_results)
      print(plot)
      dev.off()

      ## testing only the hypermethylated CpGs
      DML_up <- test_result$Probe_ID[which(test_result[[DML_colname]] == "Up")]
      tE_up <- testEnrichment(DML_up, databases = get(db))

      out_file <- paste(out_dir, "/DML/testEnrichments/custom_sets/", db, "_", CONDITION, "_", LEVEL,"_vs_", base_condition_level, "_testEnrichment_Up.pdf", sep = "")
      pdf(out_file)
      plot <- KYCG_plotEnrichAll(tE_up)
      print(plot)
      dev.off()

      ## testing only the hypermethylated CpGs
      DML_down <- test_result$Probe_ID[which(test_result[[DML_colname]] == "Down")]
      tE_down <- testEnrichment(DML_down, databases = get(db))

      out_file <- paste(out_dir, "/DML/testEnrichments/custom_sets/",db, "_", CONDITION, "_", LEVEL,"_vs_", base_condition_level, "_testEnrichment_Down.pdf", sep = "")
      pdf(out_file)
      plot <- KYCG_plotEnrichAll(tE_down)
      print(plot)
      dev.off()
    }
  }

}



#' @title GO Analysis Plot
#'
#' @description Function to plot GO results
#' @param gostres GO results from SeSAMe package
#' @param num_to_plot Number of top terms per source to plot
#' @param plot_title Option to change the plot title
#' @return Ggplot2 formatted bar plot of the top terms enriched from each source
#' @examples
#' GO_plot(gostres, 5, "Top 5 GO Terms per Source")
#' @export
GO_plot <- function(gostres, num_to_plot = 5, plot_title = "Top 5 GO Terms per Source") {

  GO_res <- gostres$result
  GO_pal <- c(`GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618",
              KEGG = "#dd4477", REAC = "#3366cc", WP = "#0099c6",
              TF = "#5574a6", MIRNA = "#22aa99", HPA = "#6633cc",
              CORUM = "#66aa00", HP = "#990099")

  plot_df <- data.frame()
  for (i in unique(GO_res$source)) {
    s_df <- GO_res[which(GO_res$source == i), ]
    plot_df <- rbind(plot_df, s_df[c(1:num_to_plot), ])
  }

  plot_df$source <- factor(plot_df$source, levels = names(GO_pal))

  plot <- ggplot(plot_df, aes(y = reorder(term_name, -p_value), x = -log10(p_value), fill = source)) +
    geom_bar(stat = "identity") +
    facet_wrap(~source, switch = "y", drop = F, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = GO_pal) +
    ggtitle(plot_title) +
    ylab(label = NULL) +
    theme_classic() +
    theme(legend.position = "none", text = element_text(face = "bold"),
          axis.line = element_line(size = 1), plot.title = element_text(hjust = .5))

  return(plot)

}

#' @title GO Analysis of DMLs
#'
#' @description Function to calculate and plot the Gene Ontology results for CONDITION-LEVEL DML CpGs
#' @param test_result Test results of DML summary stats
#' @param CONDITION Metavalue realted to colname of sample sheet
#' @param LEVEL Which comparison of condition levels to use DML results from
#' @param out_dir Path to output directory
#' @param base_condition_level Name of level used to compare DML for condition
#' @return Saves a bar plot and Go results csv for up and down GO resepectively
#' @examples
#' GO_analysis_DMLs(test_result, "Group", "Treatment", "path/out_dir")
#' @export
GO_analysis_DMLs <- function(test_result, CONDITION, LEVEL, out_dir, base_condition_level) {

  DML_colname <- paste("DML", CONDITION, LEVEL, sep = "_")

  ## Running GO on Up and Down regulated CpGs independantly
  up_DML <- test_result$Probe_ID[which(test_result[[DML_colname]] == "Up")]
  down_DML <- test_result$Probe_ID[which(test_result[[DML_colname]] == "Down")]

  genes_up <- sesameData_getGenesByProbes(up_DML)
  genes_down <- sesameData_getGenesByProbes(down_DML)

  gostres_up <- gost(genes_up$gene_name, organism = "mmusculus")
  gostres_down <- gost(genes_down$gene_name, organism = "mmusculus")

  ## plotting top enrichments for each run

  if (is.null(gostres_up)) {
    print(paste("Hyper-Methylated CpG GO has no results for the following comparison:", DML_colname))
  }else{
    up_title <- paste("Top 5 GO Terms per Source \n", CONDITION, "; ", LEVEL," vs ", base_condition_level, " UP", sep = "")
    plot_up <- GO_plot(gostres_up, plot_title = up_title)
    up_plot_file <- paste(out_dir, "/DML/GO_Enrichment/", CONDITION,"_", LEVEL,"_vs_", base_condition_level, "_GO_UP.pdf", sep = "")

    pdf(file = up_plot_file)

    print(plot_up)

    dev.off()
    #ggsave(up_plot_file, plot_up)
    ## saving dataframes for each GO run
    up_outfile <- paste(out_dir, "/DML/GO_Enrichment/", CONDITION,"_", LEVEL, "_vs_", base_condition_level, "_GO_UP.csv", sep = "")
    write_csv(gostres_up$result, up_outfile)
  }


  if (is.null(gostres_down)) {
    print(paste("Hypo-Methylated CpG GO has no results for the following comparison:", DML_colname))
  }else{
    down_title <- paste("Top 5 GO Terms per Source \n", CONDITION, "; ", LEVEL," vs ", base_condition_level, " DOWN", sep = "")
    plot_down <- GO_plot(gostres_down, plot_title = down_title)
    down_plot_file <- paste(out_dir, "/DML/GO_Enrichment/", CONDITION,"_", LEVEL,"_vs_", base_condition_level, "_GO_DOWN.pdf", sep = "")

    pdf(file = down_plot_file)

    print(plot_down)

    dev.off()
    #ggsave(down_plot_file, plot_down)

    down_outfile <- paste(out_dir, "/DML/GO_Enrichment/", CONDITION, "_", LEVEL,"_vs_", base_condition_level, "_GO_DOWN.csv", sep = "")
    write_csv(gostres_down$result, down_outfile)
  }

}


#' @title Diff Methylated Region Analysis
#'
#' @description Function to run DMR analysis for CONDITION-LEVEL DML CpGs
#' @param se Summarized experiment object output from create_SE
#' @param smry Summary statistics from DML run
#' @param CONDITION Metavalue realted to colname of sample sheet
#' @param LEVEL Which comparison of condition levels to use DML results from
#' @param out_dir Path to output directory
#' @param base_condition_level Name of level used to compare DML for condition
#' @return Writes a csv of the DMR results for each CONDITION-LEVEL contrast
#' @examples
#' DMR_DMLs(se, smry, "Group", "Treatment", "path/out_dir")
#' @export
DMR_DMLs <- function(se, smry, CONDITION, LEVEL, out_dir, base_condition_level) {

  comp <- paste(CONDITION, LEVEL, sep = "")

  merged = DMR(se, smry, comp)
  sig_merged <- merged %>% dplyr::filter(Seg_Pval_adj < 0.05)

  out_file <- paste(out_dir, "/DML/DMR_Analysis/", CONDITION, "_", LEVEL, "_vs_", base_condition_level, "_DMR.csv", sep = "")
  write_csv(sig_merged, file = out_file)

}






#' @title DML Analysis Subassembly
#'
#' @description Sub-assembly of the SeSAMe street pipeline to analyse and plot basic DML results
#' @param betas Beta values of CpGs
#' @param sample_sheet_df Data frame of the sample sheet
#' @param smry Summary statistics of DML run
#' @param formula Formula that will be used for DML
#' @param out_dir Path to output directory
#' @param pval P value threshold for significant CpG DMLs
#' @param effSize Effect size threshold for signficant CpG DMLs
#' @param custom_set_paths Path to directory of .RData file(s) containing custom CpG sets
#' @return Saving the updated test_result object now containing Up, Down, not-sig DML results for each comparison
#' @examples
#' DML_analysis(betas, sample_sheet_df, smry, ~Group, "path/out_dir")
#' @export
DML_analysis <- function(betas, sample_sheet_df, smry, formula, out_dir, pval = pval, effSize = effSize, custom_set_paths=F) {

  test_result <- summaryExtractTest(smry)
  test_result <- remove_NAs(test_result)

  out_file <- paste(out_dir, "/DML/test_result.RData", sep = "")
  save(test_result, file = out_file)

  if (nrow(test_result) == 0) {
    stop("Problem in test results, NAs present in comparisons, please check the test_result.Rdata object to identify problem")
  }

  f <- all.vars(formula[[2]])

  ## loop to generate DML for conditionLEVEL plots and output
  for (i in f) {

    ## unit test for condition presence
    test_condition_presence(i, sample_sheet_df)

    lvs <- levels(factor(sample_sheet_df[[i]]))
    for (k in lvs[c(2:length(lvs))]) { ## removing first level since this is the comparison level

      CONDITION = i
      LEVEL = k

      ls <- strsplit(LEVEL, " |-")[[1]]
      if (length(ls) > 0) {LEVEL <- paste(ls, collapse = ".")}

      base_condition_level = lvs[1] ## used for naming titles and file scheme

      ## running analysis for each CONDITION-LEVEL Subsets
      cond_test_result <- condLEVEL_DMLs(test_result, CONDITION, LEVEL, pval = pval, effSize = effSize)

      if (length(cond_test_result) == 2 && class(cond_test_result[1]) == "character") {
        warning(paste("The following DML comparison not found in the test results; ", paste(cond_test_result, collapse = " "), sep = ""))
      }else{
        volcano_plot(cond_test_result, CONDITION, LEVEL, out_dir, base_condition_level)
        heatmap_DMLs(betas, sample_sheet_df, cond_test_result, CONDITION, LEVEL, out_dir, base_condition_level)
        testEnrichment_DMLs(cond_test_result, CONDITION, LEVEL, out_dir, base_condition_level)
        if (custom_sets_path != F){
          test_enrichment_custom_sets(custom_sets_path, cond_test_result, CONDITION, LEVEL, out_dir, base_condition_level)
        }
        GO_analysis_DMLs(cond_test_result, CONDITION, LEVEL, out_dir, base_condition_level)

        se <- create_SE(betas, sample_sheet_df, formula)
        DMR_DMLs(se, smry, CONDITION, LEVEL, out_dir, base_condition_level)
      }
    }
  }
}






##### Final Assembly #####


#' @title Final Assembly of SeSAMe STREET Subassemblies/Utilities
#'
#' @description Function to run the SeSAMe Street Pipeline
#' @param Idat_dir Path to directory containing raw Idat files
#' @param out_dir Path to output directory
#' @param sample_sheet Path to SeSAMe STREET sample sheet
#' @param prep Flag code for SeSAMe Preprocessing
#' @param formula Formula for DML calculation
#' @param cores Number of cores to run DML summary stats calculation
#' @param subsample Number of beta values to use for DML calc and analysis for test runs (leave blank for full scale run)
#' @param pval P value threshold for significant CpG DMLs
#' @param effSize Effect size threshold for signficant CpG DMLs
#' @param pipeline_alacarte Option to choose which steps of the pipeline to run. Must include a selection of only the following; c("QC", "PCA", "DML")
#' @param perplexity Value to be used for tSNE of beta values in run_PCA subassembly
#' @param custom_set_paths Path to directory of .RData file(s) containing custom CpG sets
#' @return All outputs in or a sub directory of the out_dir
#' @examples
#' SeSAMe_STREET("path/Idat_dir", "path/out_dir", "path/sample_sheet.xlsx", "TQCDPB", ~ Group, 16, NA)
#' @export
SeSAMeStr <- function(Idat_dir, out_dir, sample_sheet, prep, formula, cores, subsample = NA, pval = 0.05, effSize = 0.1, perplexity = 30, pipeline_alacarte = c("QC", "DimRed", "DML"), custom_set_paths = F) {

  for (order in pipeline_alacarte) {
    if (length(pipeline_alacarte) == 0) {
      stop("pipeline_alacarte must be atleast one step secected. Please select from the following menu; c('QC', 'DimRed', 'DML')")
    }
    if (!(order %in% c("QC", "DimRed", "DML"))) {
      stop(paste("The following value is not allowed in pipeline_alacarte; '", order,"'. Please select from the following menu; c('QC', 'DimRed', 'DML')", sep = ""))
    }
  }

  ## sinking output to a log file
  #log_file <- paste(out_dir, "/SeSAMe_STREET_log.txt", sep = "")
  #sink(log_file, append = F)

  ## unit test for output architecture
  test_output_architecture(out_dir)

  ## loading in sample sheet
  message("Loading in sample sheet")
  sample_sheet_df <- SS_samplesheet(sample_sheet)

  if ("QC" %in% pipeline_alacarte) {
    ## Preprocessing and QC analysis
    message("Running QC")
    run_QC(Idat_dir, prep, out_dir)
  }

  ## Calculating beta values
  if ("betas.RData" %in% list.files(paste(out_dir, "/DML", sep = ""))) {
    message("Loading in betas.RData")
    load(paste(out_dir, "/DML/betas.RData", sep = ""))
  }else{
    message("Calculating beta values")
    betas <- get_betas(Idat_dir, prep)
    out_file <- paste(out_dir, "/DML/betas.RData", sep = "")
    save(betas, file = out_file)
  }

  if ("DimRed" %in% pipeline_alacarte) {
    ## PCA analysis
    message("Running Dimensional Redcution")
    run_DimRed(betas, sample_sheet_df, formula, out_dir, perplexity)
  }

  ## Calculating summary stats
  message("Calculating summary statistics")

  if (!(is.na(subsample))) {
    betas <- betas[c(1:subsample), ]
  }

  if ("DML" %in% pipeline_alacarte) {
    if ("smry.RData" %in% list.files(paste(out_dir, "/DML", sep = ""))) {
      message("Loading in smry.RData")
      load(paste(out_dir, "/DML/smry.RData", sep = ""))
    }else{
      message("Calculating summary statistics")
      smry <- run_DML(betas, sample_sheet_df, formula, cores, out_dir)
    }

    ## First pass of basic analysis on DML summary statistics
    message("Running DML analysis")
    DML_analysis(betas, sample_sheet_df, smry, formula, out_dir, pval = pval, effSize = effSize, custom_set_paths = custom_set_paths)
  }

  ## closing sink log
  sink()
  ##closeAllConnections()
}



