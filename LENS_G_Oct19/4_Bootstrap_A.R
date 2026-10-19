# ============================================================================
# BOOTSTRAP: MODELS B, C, D
# Reconciled version: explicit inputs, full Bailey/LENS-G objective,
# stratified regional bootstrap for C/D, fixed result slots, QA/QC.
# ============================================================================

rm(list = ls())
#setwd("C:/Users/drewgron/Desktop/ilana_R/R_test_A")
#load("Oct11.RData")

load("Q:/Hydro/projects/indopacom_PRIVATE/DREW_ilana_R/R_test_A/Oct11.RData")
setwd("Q:/Hydro/projects/indopacom_PRIVATE/DREW_ilana_R/R_test_A")

library(DEoptim)


# ---- Required inputs -------------------------------------------------------
required_objects 	<- c("final_results", "regions_present", "island_width_k", "K_VALUES")
missing 			<- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]
if (length(missing) > 0) stop(paste("Missing required objects:", paste(missing, collapse=", ")))

required_cols <- c("Region", "Orientation", "Ensemble.Precip", "Avg.Width", "Z_TD", "Z_MAX_O")
missing_cols <- required_cols[!required_cols %in% names(final_results)]
if (length(missing_cols) > 0) stop(paste("final_results is missing:", paste(missing_cols, collapse=", ")))

n_regions <- length(regions_present)
if (n_regions < 1 || anyDuplicated(regions_present)) stop("Invalid regions_present.")
if (!all(regions_present %in% unique(final_results$Region))) stop("A region is missing from final_results.")

cat("\nINPUT QA/QC\n")
cat("Observations:", nrow(final_results), "\n")
cat("Regions:", n_regions, "\n")
print(table(factor(final_results$Region, levels=regions_present)))

# ---- Settings --------------------------------------------------------------
set.seed(2026)

n_boot 				<- 1000
bootstrap_itermax 	<- 1000		# Can even go up to 50,000 (for model B at least)
bootstrap_steptol 	<- 200		# Went up to 1,000 for model B
bootstrap_reltol 	<- 1e-5
bootstrap_CR 		<- 0.5		# Cross over constant; initially 0.9 (final thing to check; higher = faster)
bootstrap_F 		<- 0.8
batch_size 			<- 500

Fr_min				<- 0.001		# Model B is ok with 0, but Model C collapses onto 0
Fr_max				<- 0.999				

# NP is DEoptim population size, not parameter count.
# D=64 matches the previously successful Master workflow.
bootstrap_NP_B 	<- 20		#8
bootstrap_NP_C 	<- 100		#16
bootstrap_NP_D 	<- 200		#64


# ---- K-correction lookup ---------------------------------------------------
width_clamped 	<- pmin(pmax(final_results$Avg.Width, 150), 2500)
Kcorr_matrix 	<- sapply(K_VALUES, function(K_value) {
  approx(island_width_k$width,
         island_width_k[[paste0("k", K_value)]],
         xout=width_clamped, rule=2)$y
})
Kcorr_matrix <- as.matrix(Kcorr_matrix)
colnames(Kcorr_matrix) <- paste0(".Kcorr_", K_VALUES)
if (nrow(Kcorr_matrix) != nrow(final_results)) stop("Kcorr_matrix row mismatch.")

# ---- Full Bailey/LENS-G calculation ---------------------------------------
calculate_K_correction <- function(data, Kcorr_matrix_val, K_values_val, K_scale=NULL) {
  orientation 	<- as.character(data$Orientation)
  K_orientation <- rep(200, nrow(data))
  K_orientation[!is.na(orientation) & orientation == "windward"] 	<- 400
  K_orientation[!is.na(orientation) & orientation == "leeward"] 	<- 50
  K_orientation[!is.na(orientation) & orientation == "offset"] 		<- 200
  K_mday <- if (is.null(K_scale)) K_orientation else K_orientation * K_scale
  K_mday <- pmin(pmax(K_mday, 25), 400)
  K_mday <- pmin(pmax(K_mday, min(K_values_val)), max(K_values_val))
  lower_index <- findInterval(K_mday, K_values_val)
  lower_index <- pmin(pmax(lower_index, 1), length(K_values_val)-1)
  upper_index <- lower_index + 1
  upper_index[K_mday >= max(K_values_val)] <- length(K_values_val)
  lower_index[K_mday >= max(K_values_val)] <- length(K_values_val)
  row_index <- seq_len(nrow(data))
  lo <- Kcorr_matrix_val[cbind(row_index, lower_index)]
  hi <- Kcorr_matrix_val[cbind(row_index, upper_index)]
  Klo <- K_values_val[lower_index]
  Khi <- K_values_val[upper_index]
  ifelse(Khi == Klo, lo, lo + ((K_mday-Klo)/(Khi-Klo))*(hi-lo))
}

calculate_ZMAX <- function(params, data, regions_val, K_values_val, Kcorr_matrix_val, model=c("B","C","D")) {
  model <- match.arg(model)
  nr <- length(regions_val)
  if (model == "B") {
    if (length(params) != 1) stop("Model B requires 1 parameter.")
    rf <- setNames(rep(params[1], nr), regions_val)
    ks <- NULL
  } else if (model == "C") {
    if (length(params) != nr) stop("Model C parameter count mismatch.")
    rf <- setNames(params, regions_val)
    ks <- NULL
  } else {
    if (length(params) != 2*nr) stop("Model D parameter count mismatch.")
    rf <- setNames(params[seq_len(nr)], regions_val)
    ks <- setNames(params[(nr+1):(2*nr)], regions_val)
  }
  ri <- match(as.character(data$Region), regions_val)
  if (anyNA(ri)) return(rep(NA_real_, nrow(data)))
  recharge_frac <- rf[ri]
  if (is.null(ks)) {
    K_correction <- calculate_K_correction(data, Kcorr_matrix_val, K_values_val)
  } else {
    K_correction <- calculate_K_correction(data, Kcorr_matrix_val, K_values_val, ks[ri])
  }
  rainfall_m <- data$Ensemble.Precip
  width <- data$Avg.Width
  Z_TD <- data$Z_TD
  R <- rainfall_m * recharge_frac
  B <- 5 * exp(-0.00461 * width)
  C <- 1 - 0.15 * exp(-0.00288 * width)
  Y <- -20 * (1 - exp(-0.001386 * width))
  (((Z_TD-Y)*R)/(B+R) + Y) * K_correction * C
}

################# check on the sum(valid) < 10 criteria 

make_objective <- function(model_name) {
  force(model_name)
  function(params, data, regions_val, K_values_val, Kcorr_matrix_val) {
    z <- calculate_ZMAX(params, data, regions_val, K_values_val, Kcorr_matrix_val, model_name)
    valid <- is.finite(z) & is.finite(data$Z_MAX_O)
    if (sum(valid) < 10) return(1e10)
    sqrt(mean((z[valid]-data$Z_MAX_O[valid])^2))
  }
}
obj_fn_B <- make_objective("B")
obj_fn_C <- make_objective("C")
obj_fn_D <- make_objective("D")


###########################################################################
# ---- Regional data ---------------------------------------------------------
data_by_region 	<- split(	final_results, 
							factor(final_results$Region, 
							levels=regions_present))
regional_n 		<- vapply(data_by_region, nrow, integer(1))
if (any(regional_n < 1)) stop("A region has no observations.")
cat("\nRegional bootstrap sample sizes:\n")
print(regional_n)
###########################################################################


###########################################################################
# ---- Generic bootstrap helpers --------------------------------------------
run_DE <- function(fn, lower, upper, data, regions_val, Kcorr, np) {
  tryCatch(
    DEoptim(fn = fn, lower = lower, upper = upper, data = data,
            regions_val = regions_val, K_values_val = K_VALUES,
            Kcorr_matrix_val = Kcorr,
            control = 
				DEoptim.control(	itermax = bootstrap_itermax, 
									NP		= np,
                                    CR		= bootstrap_CR, 
									F		= bootstrap_F, 
									trace 	= FALSE,
                                    steptol 		= bootstrap_steptol, 
									reltol  		= bootstrap_reltol,
                                    storepopfrom	= bootstrap_itermax+1
								)
			),
    error=function(e) list(.error=conditionMessage(e))
  )
}

stratified_indices <- function(data, regions_val) 
	{
		pieces <- lapply(regions_val, function(r) {
			ii <- which(data$Region == r)
			ii[sample(seq_along(ii), length(ii), replace=TRUE)]
													}
						)
  unlist(pieces, use.names=FALSE)
}

check_ids <- function(results, n, model) {
  ids <- vapply(results, function(x) if (is.null(x)) NA_integer_ else x$replicate, integer(1))
  if (!identical(ids, seq_len(n))) stop(paste(model, "replicate IDs are missing, duplicated, or out of order."))
}



#######################################################################################
#######################################################################################
#######################################################################################
# ---- Model B ---------------------------------------------------------------
#######################################################################################
#######################################################################################
#######################################################################################


bootstrap_B <- function(i) {
  idx 	<- sample(	seq_len(nrow(final_results)), 
					nrow(final_results), 
					replace = TRUE)
  d 	<- final_results[idx,,drop=FALSE]; rownames(d) <- NULL
  k 	<- Kcorr_matrix[idx,,drop=FALSE]
  fit 	<- run_DE	(
			fn		= obj_fn_B, 
			lower 	= Fr_min, 
			upper	= Fr_max, 
			data	= d, 
			regions_val	= regions_present, 
			Kcorr	= k,						# I'm still not sure entirely what Kcorr does 
			np		= bootstrap_NP_B
					)
	if (".error" %in% names(fit)) 
		return(list(
				success		= FALSE, 
				replicate	= i, 
				params		= NULL, 
				objective	= NA_real_, 
				error		= fit$.error))
			
		list(	success		= TRUE, 
				replicate	= i, 
				# params		= c(rf_global = fit$optim$bestmem) #[1], 
				params		= fit$optim$bestmem,
				objective 	= fit$optim$bestval, 
				# n_generations = length(fit$member$bestvalit)
				error		= NA_character_
								
			)
}

#####################################################################
########### RUN MODEL B BOOTSRAP !!
#####################################################################

start_B <- Sys.time()
n_success <- 0L

boot_B_results <- lapply(seq_len(n_boot), function(i) {

  result <- bootstrap_B(i)

  if (isTRUE(result$success)) {
    n_success <<- n_success + 1L
  }

  elapsed_sec <- as.numeric(
    difftime(Sys.time(), start_B, units = "secs")
  )
  remaining_sec <- (elapsed_sec / i) * (n_boot - i)

  cat(sprintf(
    paste0(
      "\rCompleted: %d/%d (%3.0f%%) | ",
      "Successful: %d | Failed: %d | ",
      "Estimated remaining: %.1f min"
    ),
    i,
    n_boot,
    100 * i / n_boot,
    n_success,
    i - n_success,
    remaining_sec / 60
  ))
  flush.console()

  result
})

cat("\n")

elapsed_B <- difftime(Sys.time(), start_B, units = "mins")
print(elapsed_B)

#####################################################################
###############	END BOOTSTRAP RUN 
#####################################################################
#####################################################################

check_ids(	boot_B_results,
			n_boot,
			"Model B")
			
B_ok 	<- vapply(boot_B_results,function(x) isTRUE(x$success),logical(1)); 
B_ok

boot_B 	<- as.data.frame(do.call(rbind,lapply(boot_B_results[B_ok], function(x)x$params)))
	#if (any(boot_B$rf_global<.10 | boot_B$rf_global>.90)) stop("Model B bound violation.")

names(boot_B) <- "rf_global"

boot_B.obj 	<- as.data.frame(do.call(rbind,lapply(boot_B_results[B_ok],function(x)x$objective)))

write.csv(boot_B,"bootstrap_Model_B_joint.csv",row.names=FALSE)

hist(unlist(boot_B))
plot(unlist(boot_B))
plot(unlist(boot_B), unlist(boot_B.obj), type = "p", pch = 21)


## take best 100

b.step = cbind(boot_B, boot_B.obj)
b.best = b.step[order(b.step[,2])[1:100],]
b.best

hist(unlist(b.best[,1]), freq = FALSE)
plot(unlist(b.best[,1]))



#####################################################################
###############	END MODEL B
#####################################################################
#####################################################################



#####################################################################
#####################################################################
#####################################################################
# ---- Model C --------------------------------------------------
#####################################################################
#####################################################################
#####################################################################

# Check on lower and upper bounds on recharge fraction; set above

bootstrap_C <- function(i) {
  idx 	<- stratified_indices(final_results, regions_present)
  d 	<- final_results[idx,,drop = FALSE]; rownames(d) <- NULL
  k 	<- Kcorr_matrix[idx,,drop = FALSE]
  fit 	<- run_DE	(	obj_fn_C, 
						rep(Fr_min, n_regions), 
						rep(Fr_max, n_regions),    #changed from Fr_min to Fr_max
						d, regions_present, k, bootstrap_NP_C
					)
  if (".error" %in% names(fit)) 
	return(list(
			success		= FALSE, 
			replicate	= i, 
			params		= NULL, 
			objective	= NA_real_, 
			error		= fit$.error)
			)
  
  p <- fit$optim$bestmem; 
	names(p) <- paste0("rf_",regions_present)
	
	list(	success 	= TRUE, 
			replicate	= i, 
			params		= p, 
			objective	= fit$optim$bestval, 
			error		= NA_character_
		)
}

#####################################################################
########### BEGIN Model C bootstrap #################################
#####################################################################
# 
# start_C <- Sys.time()
# n_success <- 0L
# 
# boot_C_results <- lapply(seq_len(n_boot), function(i) {
# 
#   result <- bootstrap_C(i)
# 
#   if (isTRUE(result$success)) {
#     n_success <<- n_success + 1L
#   }
# 
#   elapsed_sec <- as.numeric(
#     difftime(Sys.time(), start_C, units = "secs")
#   )
#   remaining_sec <- (elapsed_sec / i) * (n_boot - i)
# 
#   cat(sprintf(
#     paste0(
#       "\rCompleted: %d/%d (%3.0f%%) | ",
#       "Successful: %d | Failed: %d | ",
#       "Estimated remaining: %.1f min"
#     ),
#     i,
#     n_boot,
#     100 * i / n_boot,
#     n_success,
#     i - n_success,
#     remaining_sec / 60
#   ))
#   flush.console()
# 
#   result
# })
# 
# cat("\n")
# 
# elapsed_C <- difftime(Sys.time(), start_C, units = "mins")
# print(elapsed_C)













start_C <- Sys.time()

boot_C_results <- vector("list", n_boot)
for (s in seq(1,n_boot,by=batch_size)) {
  e <- min(s+batch_size-1,n_boot); ii <- s:e
  boot_C_results[ii] <- lapply(ii, bootstrap_C)
  cat("Model C:", e, "/", n_boot, "\n")
}

#####################################################################
########### END Model C bootstrap #################################
#####################################################################






check_ids(boot_C_results,n_boot,"Model C")
C_ok 		<- vapply(boot_C_results,function(x)isTRUE(x$success),logical(1))
boot_C 		<- as.data.frame(do.call(rbind,lapply(boot_C_results[C_ok],function(x)x$params)))
expected_C 	<- paste0("rf_",regions_present)
if (!identical(colnames(boot_C),expected_C)) stop("Model C column mismatch.")
# if (any(as.matrix(boot_C)<.10 | as.matrix(boot_C)>.90)) stop("Model C bound violation.")
write.csv(boot_C,"bootstrap_Model_C_joint.csv", row.names = FALSE)

elapsed_C <- difftime(Sys.time(),start_C,units="hours")
cat("Model C complete:",sum(C_ok),"/",n_boot,"; elapsed",round(as.numeric(elapsed_C),2),"hours\n")


####################################################################################
####################################################################################
####################################################################################
# ---- Model D ---------------------------------------------------------------
####################################################################################
####################################################################################
####################################################################################


cat("\n################################################################################\nMODEL D: REGIONAL RF + K SCALAR\n################################################################################\n")
cat("NP =",bootstrap_NP_D,"(matches previous successful Master workflow)\n")
cat("Note: DEoptim may print its heuristic NP warning; NP=64 is retained deliberately.\n")
bootstrap_D <- function(i) {
  idx <- stratified_indices(final_results, regions_present)
  d <- final_results[idx,,drop=FALSE]; rownames(d) <- NULL
  k <- Kcorr_matrix[idx,,drop=FALSE]
  lower <- c(rep(Fr_min,n_regions),rep(.25,n_regions))
  upper <- c(rep(Fr_max,n_regions),rep(4.00,n_regions))
  fit <- run_DE(obj_fn_D, lower, upper, d, regions_present, k, bootstrap_NP_D)
  if (".error" %in% names(fit)) return(list(success=FALSE, replicate=i, params=NULL, objective=NA_real_, error=fit$.error))
  p <- fit$optim$bestmem
  names(p) <- c(paste0("rf_",regions_present),paste0("Ks_",regions_present))
  list(success=TRUE, replicate=i, params=p, objective=fit$optim$bestval, error=NA_character_)
}
start_D <- Sys.time(); boot_D_results <- vector("list", n_boot)
for (s in seq(1,n_boot,by=batch_size)) {
  e <- min(s+batch_size-1,n_boot); ii <- s:e
  # Fixed-slot assignment. No reset/reinitialization occurs here.
  batch_results <- lapply(ii, bootstrap_D)
  boot_D_results[ii] <- batch_results
  cat("Model D:", e, "/", n_boot, "\n")
}
check_ids(boot_D_results,n_boot,"Model D")
D_ok <- vapply(boot_D_results,function(x)isTRUE(x$success),logical(1))
if (!any(D_ok)) stop("Model D produced no successful replicates.")
boot_D <- as.data.frame(do.call(rbind,lapply(boot_D_results[D_ok],function(x)x$params)))
expected_D <- c(paste0("rf_",regions_present),paste0("Ks_",regions_present))
if (!identical(colnames(boot_D),expected_D)) stop("Model D column mismatch.")
rf_cols <- paste0("rf_",regions_present); ks_cols <- paste0("Ks_",regions_present)
if (any(as.matrix(boot_D[,rf_cols,drop=FALSE])<Fr_min | as.matrix(boot_D[,rf_cols,drop=FALSE])>Fr_max)) stop("Model D RF bound violation.")
if (any(as.matrix(boot_D[,ks_cols,drop=FALSE])<.25 | as.matrix(boot_D[,ks_cols,drop=FALSE])>4.00)) stop("Model D K bound violation.")
write.csv(boot_D,"bootstrap_Model_D_joint.csv",row.names=FALSE)
elapsed_D <- difftime(Sys.time(),start_D,units="hours")


n_plot <- min(450, nrow(boot_B), nrow(boot_C), nrow(boot_D))
idx_plot <- seq_len(n_plot)


# ==============================================================================
# ============================================================================
# PART 4: PLOTTING (All three models)
# ============================================================================
# ==============================================================================

# ==============================================================================
# MODEL B
# ==============================================================================


rf_B_values <- boot_B$rf_global

par(mfrow = c(1, 2))
par(mar = c(3, 3, 2, 1))
par(mgp = c(2, 0.7, 0))

plot(
  rf_B_values[500:1000],
  ylim = c(0, 1.1),
  type = "l",
  pch = 20,
  xlab = "Bootstrap replicate",
  ylab = "Global Recharge Fraction",
  main = "Model B: Trace Plot",
  las = 1,
  cex.main = 1.2
)

# Histogram
hist(
  rf_B_values,
  freq = FALSE,
  xlim = c(0, 1),
  xlab = "Global Recharge Fraction",
  ylab = "Probability Density",
  main = "Model B: Distribution",
  las = 1,
  cex.main = 1.2,
  breaks = 20
)



dev.print(
  png,
  filename = file.path(getwd(), "Model_B_trace_and_distribution.png"),
  width = 2400,
  height = 1200,
  res = 300
)


#CORRECTED PLOTTING CODE -------------------------------------------------------


# Correct the Model B parameter name
names(boot_B) <- "rf_global"
rf_B_values <- boot_B$rf_global

# Save the figure
png(
  file.path(getwd(), "Model_B_trace_and_distribution.png"),
  width = 2400, height = 1200, res = 300
)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Trace-style plot of all 100 bootstrap estimates
plot(
  seq_along(rf_B_values), rf_B_values,
  ylim = c(0, 1),
  type = "l",
  xlab = "Bootstrap replicate",
  ylab = "Global Recharge Fraction",
  main = "Model B: Bootstrap Estimates"
)

# Bootstrap distribution
hist(
  boot_B$rf_global,
  probability = TRUE,
  xlim = c(0, 1),
  xlab = "Global Recharge Fraction",
  ylab = "Probability Density",
  main = "Model B: Bootstrap Distribution",
  breaks = 15
)

dev.off()


# ==============================================================================
# MODEL C
# ==============================================================================


rf_cols_C <- names(boot_C)
region_names_C <- gsub("^rf_", "", rf_cols_C)


windows()

par(mfrow = c(4, 4))
par(mar = c(3.0, 3.0, 2.0, 1))
par(oma = c(3.5, 3.5, 3.5, 1))

# First 8 panels: trace plots
for (i in seq_along(rf_cols_C)) {
  
  rf_values <- boot_C[[rf_cols_C[i]]]
  
  plot(
    rf_values[idx_plot],
    ylim = c(0, 1.1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    type = "b",
    pch = 20,
    cex = 0.45,
    lwd = 0.5
  )
  
  box()
  
  title(
    main = region_names_C[i],
    cex.main = 0.9
  )
  
  axis(1, cex.axis = 0.8)
  axis(2, las = 1, cex.axis = 0.8)
}

# Second 8 panels: histograms
for (i in seq_along(rf_cols_C)) {
  
  rf_values <- boot_C[[rf_cols_C[i]]]
  
  hist(
    rf_values,
    main = region_names_C[i],
    freq = FALSE,
    xlim = c(0, 1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    breaks = 20,
    col = "lightgray",
    border = "gray30",
    cex.main = 0.9
  )
  
  box()
  axis(1, cex.axis = 0.8)
  axis(2, las = 1, cex.axis = 0.8)
}

mtext(
  "Bootstrap replicate / Recharge Fraction",
  side = 1,
  outer = TRUE,
  line = 2.0,
  cex = 1.1
)

mtext(
  "Model C: Regional Recharge Fraction",
  side = 3,
  outer = TRUE,
  line = 1.0,
  cex = 1.3,
  font = 2
)

dev.print(
  png,
  filename = file.path(getwd(), "Model_C_trace_and_histograms.png"),
  width = 3000,
  height = 3000,
  res = 300
)

cat("Model C plots created.\n\n")


# Summary statistics for Model C bootstrap estimates
for (col in names(boot_C)) {
  
  x <- boot_C[[col]]
  x <- x[is.finite(x)]
  
  cat("\n", col, "\n")
  cat("-----------------------------\n")
  print(summary(x))
  cat("Standard deviation:", sd(x), "\n")
  cat("2.5th percentile:", unname(quantile(x, 0.025)), "\n")
  cat("97.5th percentile:", unname(quantile(x, 0.975)), "\n")
  cat("Estimates >= 0.99:", sum(x >= 0.99), "of", length(x), "\n")
}

# ==============================================================================
# MODEL D
# ==============================================================================

cat("Creating Model D plots...\n")

# Model D has:
#   columns 1:8   = regional recharge fractions
#   columns 9:16  = regional K scalars

rf_cols_D <- paste0("rf_", regions_present)
K_cols_D  <- paste0("Ks_", regions_present)

# Confirm expected columns exist
missing_D_cols <- c(rf_cols_D, K_cols_D)[
  !c(rf_cols_D, K_cols_D) %in% names(boot_D)
]

if (length(missing_D_cols) > 0) {
  stop(
    paste(
      "Missing expected Model D bootstrap columns:",
      paste(missing_D_cols, collapse = ", ")
    )
  )
}

# Regional sample sizes used to support Model D parameters
region_count_table <- table(final_results$Region)

region_counts <- data.frame(
  Region = names(region_count_table),
  n = as.numeric(region_count_table),
  row.names = NULL
)

region_counts <- region_counts[
  match(regions_present, region_counts$Region),
  ,
  drop = FALSE
]


# ==============================================================================
# MODEL D: RF TRACES + HISTOGRAMS
# ==============================================================================

windows()

par(mfrow = c(4, 4))
par(mar = c(3.0, 3.0, 2.0, 1))
par(oma = c(3.5, 3.5, 3.5, 1))

# First 8 panels: trace plots
for (i in seq_along(rf_cols_D)) {
  
  rf_values <- boot_D[[rf_cols_D[i]]]
  
  plot(
    rf_values[idx_plot],
    ylim = c(0, 1.1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    type = "b",
    pch = 20,
    cex = 0.45,
    lwd = 0.5
  )
  
  box()
  
  title(
    main = paste0(
      region_counts$Region[i],
      " (n=",
      region_counts$n[i],
      ")"
    ),
    cex.main = 0.9
  )
  
  axis(1, cex.axis = 0.8)
  axis(2, las = 1, cex.axis = 0.8)
}


# Second 8 panels: histograms
for (i in seq_along(rf_cols_D)) {
  
  rf_values <- boot_D[[rf_cols_D[i]]]
  
  hist(
    rf_values,
    main = paste0(
      region_counts$Region[i],
      " (n=",
      region_counts$n[i],
      ")"
    ),
    freq = FALSE,
    xlim = c(0, 1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    breaks = 20,
    col = "lightgray",
    border = "gray30",
    cex.main = 0.9
  )
  
  box()
  axis(1, cex.axis = 0.8)
  axis(2, las = 1, cex.axis = 0.8)
}

mtext(
  "Bootstrap replicate / Recharge Fraction",
  side = 1,
  outer = TRUE,
  line = 2.0,
  cex = 1.1
)

mtext(
  "Model D: Regional Recharge Fraction",
  side = 3,
  outer = TRUE,
  line = 1.0,
  cex = 1.3,
  font = 2
)

dev.print(
  png,
  filename = file.path(getwd(), "Model_D_RF_traces_and_histograms.png"),
  width = 3000,
  height = 3000,
  res = 300
)


# Save PDF with maximum width of 15 cm
ggsave(
  filename = file.path(getwd(), "Model_D_RF_traces_and_histograms.pdf"),
  plot = Model_D_RF_traces_and_histograms,
  width = 15,
  height = 12,
  units = "cm",
  dpi = 600,
  bg = "white"
)


# ==============================================================================
# MODEL D: K-SCALAR TRACES + HISTOGRAMS
# ==============================================================================

windows()

par(mfrow = c(4, 4))
par(mar = c(3.0, 3.0, 2.0, 1))
par(oma = c(3.5, 3.5, 3.5, 1))

# First 8 panels: trace plots
for (i in seq_along(K_cols_D)) {
  
  K_values_region <- boot_D[[K_cols_D[i]]]
  
  plot(
    K_values_region[idx_plot],
    ylim = c(0, 4.2),
    axes = FALSE,
    xlab = "",
    ylab = "",
    type = "b",
    pch = 20,
    cex = 0.45,
    lwd = 0.5
  )
  
  box()
  
  title(
    main = paste0(
      region_counts$Region[i],
      " (n=",
      region_counts$n[i],
      ")"
    ),
    cex.main = 0.9
  )
  
  axis(1, cex.axis = 0.8)
  axis(2, las = 1, cex.axis = 0.8)
}


# Second 8 panels: histograms
for (i in seq_along(K_cols_D)) {
  
  K_values_region <- boot_D[[K_cols_D[i]]]
  
  hist(
    K_values_region,
    main = paste0(
      region_counts$Region[i],
      " (n=",
      region_counts$n[i],
      ")"
    ),
    freq = FALSE,
    xlim = c(0.25, 4.0),
    axes = FALSE,
    xlab = "",
    ylab = "",
    breaks = 20,
    col = "lightgray",
    border = "gray30",
    cex.main = 0.9
  )
  
  box()
  axis(1, cex.axis = 0.8)
  axis(2, las = 1, cex.axis = 0.8)
}

mtext(
  "Bootstrap replicate / K scalar",
  side = 1,
  outer = TRUE,
  line = 2.0,
  cex = 1.1
)

mtext(
  "Model D: Regional Hydraulic-Conductivity Scalar",
  side = 3,
  outer = TRUE,
  line = 1.0,
  cex = 1.3,
  font = 2
)

dev.print(
  png,
  filename = file.path(getwd(), "Model_D_K_traces_and_histograms.png"),
  width = 3000,
  height = 3000,
  res = 300
)

# Save PDF with maximum width of 15 cm
ggsave(
  filename = file.path(getwd(), "Model_D_K_traces_and_histograms.pdf"),
  plot = Model_D_K_traces_and_histograms,
  width = 15,
  height = 12,
  units = "cm",
  dpi = 600,
  bg = "white"
)

cat("Model D K-scalar plots created.\n\n")


# ============================================
# Model D: Bootstrap summary statistics
# Regional recharge fractions and K scalars
# ============================================

summarize_bootstrap <- function(boot_data, parameter_name,
                                lower_bound, upper_bound) {
  
  for (col in names(boot_data)) {
    
    x <- boot_data[[col]]
    x <- x[is.finite(x)]
    
    cat("\n", col, "\n")
    cat("-----------------------------------\n")
    
    print(summary(x))
    
    cat("Standard deviation:", sd(x), "\n")
    cat("2.5th percentile:",
        unname(quantile(x, 0.025)), "\n")
    cat("97.5th percentile:",
        unname(quantile(x, 0.975)), "\n")
    
    cat("Estimates at lower bound:",
        sum(x <= lower_bound + 1e-6), "of", length(x), "\n")
    
    cat("Estimates at upper bound:",
        sum(x >= upper_bound - 1e-6), "of", length(x), "\n")
    
    cat("Estimates >= 99% of upper bound:",
        sum(x >= upper_bound * 0.99), "of", length(x), "\n")
  }
}

# Regional recharge fractions
summarize_bootstrap(
  boot_D[, grep("^rf_", names(boot_D)), drop = FALSE],
  parameter_name = "Recharge fraction",
  lower_bound = Fr_min,
  upper_bound = Fr_max
)

# Regional hydraulic-conductivity scalars
ks_cols <- grep("^Ks_", names(boot_D), value = TRUE)

summarize_bootstrap(
  boot_D[, ks_cols, drop = FALSE],
  parameter_name = "K scalar",
  lower_bound = 0.25,
  upper_bound = 4.00
)

# ==============================================================================
# MODEL D: RF VS K SCATTERPLOTS
# ==============================================================================

cat("Creating Model D RF-vs-K scatterplots...\n")

windows()

par(mfrow = c(2, 4))
par(mar = c(0, 0, 0, 0))
par(oma = c(4, 4, 4, 2))

for (i in seq_along(regions_present)) {
  
  rf_values <- boot_D[[rf_cols_D[i]]]
  K_values_region <- boot_D[[K_cols_D[i]]]
  
  plot(
    x = K_values_region,
    y = rf_values,
    xlim = c(0, 4.5),
    ylim = c(0, 1.2),
    axes = FALSE,
    xlab = "",
    ylab = "",
    pch = 20,
    cex = 0.35
  )
  
  box()
  
  text(
    x = 0,
    y = 1.1,
    paste(
      region_counts$Region[i],
      paste0("(n=", region_counts$n[i], ")")
    ),
    pos = 4
  )
  
  if (i == 1) {
    axis(2)
    axis(3)
    mtext("K scalar", side = 3, line = 2.5)
    mtext("Rf", side = 2, line = 2.5)
  }
}

dev.print(
  png,
  filename = file.path(getwd(), "Model_D_RF_vs_K_scatterplots.png"),
  width = 3000,
  height = 1800,
  res = 300
)

cat("Model D RF-vs-K plots created.\n\n")


# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n============================================================\n")
cat("BOOTSTRAP PLOTTING COMPLETE!\n")
cat("============================================================\n\n")

cat("Successful bootstrap replicates:\n")
cat("  Model B: ", nrow(boot_B), "\n", sep = "")
cat("  Model C: ", nrow(boot_C), "\n", sep = "")
cat("  Model D: ", nrow(boot_D), "\n\n", sep = "")

cat("Output plots:\n")
cat("  - Model_B_trace_and_distribution.png\n")
cat("  - Model_C_trace_and_histograms.png\n")
cat("  - Model_D_RF_traces_and_histograms.png\n")
cat("  - Model_D_K_traces_and_histograms.png\n")
cat("  - Model_D_RF_vs_K_scatterplots.png\n\n")

cat("Trace plots use the final ", n_plot,
    " bootstrap replicates; histograms and RF-vs-K plots use all successful ",
    nrow(boot_D), " bootstrap replicates.\n", sep = "")

cat("============================================================\n\n")




# ==============================================================================
# PART 5: PARAMETER-UNCERTAINTY INTERVALS AND COVERAGE
# ==============================================================================

cat("\n\n")
cat("################################################################################\n")
cat("PART 5: PARAMETER-UNCERTAINTY INTERVALS AND COVERAGE\n")
cat("################################################################################\n\n")

# ------------------------------------------------------------------------------
# Prediction helper
# ------------------------------------------------------------------------------

predict_bootstrap_intervals <- function(boot_params,
                                        model_name,
                                        data,
                                        regions_val,
                                        K_values_val,
                                        Kcorr_matrix_val) {
  
  n_boot_local <- nrow(boot_params)
  n_obs <- nrow(data)
  
  prediction_matrix <- matrix(
    NA_real_,
    nrow = n_boot_local,
    ncol = n_obs
  )
  
  cat("Generating predictions for Model", model_name, "...\n")
  
  for (i in seq_len(n_boot_local)) {
    
    params_i <- as.numeric(boot_params[i, ])
    
    prediction_matrix[i, ] <- calculate_ZMAX(
      params = params_i,
      data = data,
      regions_val = regions_val,
      K_values_val = K_values_val,
      Kcorr_matrix_val = Kcorr_matrix_val,
      model = model_name
    )
  }
  
  lower <- apply(
    prediction_matrix,
    2,
    quantile,
    probs = 0.025,
    na.rm = TRUE
  )
  
  upper <- apply(
    prediction_matrix,
    2,
    quantile,
    probs = 0.975,
    na.rm = TRUE
  )
  
  median_pred <- apply(
    prediction_matrix,
    2,
    median,
    na.rm = TRUE
  )
  
  observed <- data$Z_MAX_O
  
  covered <- observed >= lower & observed <= upper
  
  interval_width <- upper - lower
  
  data.frame(
    Region = data$Region,
    Observed = observed,
    Lower_95 = lower,
    Median_Prediction = median_pred,
    Upper_95 = upper,
    Interval_Width = interval_width,
    Covered_95 = covered
  )
}


# ------------------------------------------------------------------------------
# Calculate intervals for Models B, C, and D
# ------------------------------------------------------------------------------

coverage_B <- predict_bootstrap_intervals(
  boot_params = boot_B,
  model_name = "B",
  data = final_results,
  regions_val = regions_present,
  K_values_val = K_VALUES,
  Kcorr_matrix_val = Kcorr_matrix
)

write.csv(coverage_B, "coverage_B.csv")

coverage_C <- predict_bootstrap_intervals(
  boot_params = boot_C,
  model_name = "C",
  data = final_results,
  regions_val = regions_present,
  K_values_val = K_VALUES,
  Kcorr_matrix_val = Kcorr_matrix
)

write.csv(coverage_C, "coverage_C.csv")

coverage_D <- predict_bootstrap_intervals(
  boot_params = boot_D,
  model_name = "D",
  data = final_results,
  regions_val = regions_present,
  K_values_val = K_VALUES,
  Kcorr_matrix_val = Kcorr_matrix
)

write.csv(coverage_D, "coverage_D.csv")


#-------------------------------------------------------------------------------

# ==============================================================================
# TABLE S8: FULL-DATA ESTIMATES + BOOTSTRAP 95% INTERVALS
# ==============================================================================

cat("\nCreating Table S8...\n")

# ---- 1. Regional sample sizes -----------------------------------------------

region_counts <- table(
  factor(final_results$Region, levels = regions_present)
)

# ---- 2. Fit Models C and D to the full calibration dataset -----------------
# These are full-data estimates, not bootstrap means.
# Uses the optimization settings defined in this script.

fit_C_full <- run_DE(
  fn = obj_fn_C,
  lower = rep(Fr_min, n_regions),
  upper = rep(Fr_max, n_regions),
  data = final_results,
  regions_val = regions_present,
  Kcorr = Kcorr_matrix,
  np = bootstrap_NP_C
)

if (".error" %in% names(fit_C_full)) {
  stop(paste("Full-data Model C failed:", fit_C_full$.error))
}

fit_D_full <- run_DE(
  fn = obj_fn_D,
  lower = c(rep(Fr_min, n_regions), rep(0.25, n_regions)),
  upper = c(rep(Fr_max, n_regions), rep(4.00, n_regions)),
  data = final_results,
  regions_val = regions_present,
  Kcorr = Kcorr_matrix,
  np = bootstrap_NP_D
)

if (".error" %in% names(fit_D_full)) {
  stop(paste("Full-data Model D failed:", fit_D_full$.error))
}

# Name the full-data parameter estimates
rf_C_full <- setNames(
  fit_C_full$optim$bestmem,
  paste0("rf_", regions_present)
)

rf_D_full <- setNames(
  fit_D_full$optim$bestmem[seq_len(n_regions)],
  paste0("rf_", regions_present)
)

ks_D_full <- setNames(
  fit_D_full$optim$bestmem[n_regions + seq_len(n_regions)],
  paste0("Ks_", regions_present)
)

# ---- 3. Helper function for bootstrap percentile intervals -----------------

get_boot_summary <- function(boot_data, parameter) {
  
  if (!parameter %in% names(boot_data)) {
    stop(paste("Parameter not found in bootstrap data:", parameter))
  }
  
  x <- boot_data[[parameter]]
  x <- x[is.finite(x)]
  
  if (length(x) < 2) {
    stop(paste("Too few valid bootstrap estimates for:", parameter))
  }
  
  c(
    lower_95 = unname(quantile(x, 0.025)),
    upper_95 = unname(quantile(x, 0.975)),
    n_successful = length(x)
  )
}

# ---- 4. Assemble Table S8 ----------------------------------------------------

table_S8 <- data.frame(
  Region = regions_present,
  n = as.integer(region_counts$n),
  Model_C_RF = as.numeric(rf_C_full),
  Model_D_RF = as.numeric(rf_D_full),
  Model_D_K_scalar = as.numeric(ks_D_full),
  Model_C_RF_lower_95 = NA_real_,
  Model_C_RF_upper_95 = NA_real_,
  Model_D_RF_lower_95 = NA_real_,
  Model_D_RF_upper_95 = NA_real_,
  Model_D_K_lower_95 = NA_real_,
  Model_D_K_upper_95 = NA_real_,
  row.names = NULL
)

for (i in seq_along(regions_present)) {
  
  region <- regions_present[i]
  
  c_summary <- get_boot_summary(
    boot_C,
    paste0("rf_", region)
  )
  
  d_rf_summary <- get_boot_summary(
    boot_D,
    paste0("rf_", region)
  )
  
  d_k_summary <- get_boot_summary(
    boot_D,
    paste0("Ks_", region)
  )
  
  table_S8$Model_C_RF_lower_95[i] <- c_summary["lower_95"]
  table_S8$Model_C_RF_upper_95[i] <- c_summary["upper_95"]
  
  table_S8$Model_D_RF_lower_95[i] <- d_rf_summary["lower_95"]
  table_S8$Model_D_RF_upper_95[i] <- d_rf_summary["upper_95"]
  
  table_S8$Model_D_K_lower_95[i] <- d_k_summary["lower_95"]
  table_S8$Model_D_K_upper_95[i] <- d_k_summary["upper_95"]
}

# ---- 5. Print results and save CSV -----------------------------------------

print(table_S8, digits = 4)

write.csv(
  table_S8,
  "Table_S8_regional_parameters_and_uncertainty.csv",
  row.names = FALSE
)

cat("\nSuccessful bootstrap replicates:\n")
cat("Model C:", nrow(boot_C), "of", n_boot, "\n")
cat("Model D:", nrow(boot_D), "of", n_boot, "\n")

cat("\nFull-data Model C RMSE:", fit_C_full$optim$bestval, "\n")
cat("Full-data Model D RMSE:", fit_D_full$optim$bestval, "\n")

cat(
  "\nSaved Table_S8_regional_parameters_and_uncertainty.csv\n"
)

