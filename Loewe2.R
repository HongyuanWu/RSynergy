library(data.table)

Loewe2 <- function (ChData, drug.row, drug.col, quiet = TRUE)
{
  response.mat <- synergyfinder::ReshapeData(ChData)$dose.response.mats[[1]]
  if (quiet) {
    options(warn = -1)
  }
  #drug.row <- ExtractSingleDrug(response.mat, dim = "row")
  #drug.col <- ExtractSingleDrug(response.mat, dim = "col")
  # con <- vapply(list(drug.col.model, drug.row.model), is.null, 
  #               logical(1))
  # if (!all(!con)) {
    drug.row.model <- synergyfinder::FitDoseResponse(drug.row)
    drug.col.model <- synergyfinder::FitDoseResponse(drug.col)
  # }
  drug.row.par <- stats::coef(drug.row.model)
  drug.row.type <- synergyfinder::FindModelType(drug.row.model)
  drug.col.par <- stats::coef(drug.col.model)
  drug.col.type <- synergyfinder::FindModelType(drug.col.model)
  drug.row$dose[drug.row$dose == 0] = 10^-10
  drug.col$dose[drug.col$dose == 0] = 10^-10
  loewe.mat <- response.mat
  eq <- switch(paste(drug.col.type, drug.row.type), `LL.4 LL.4` = synergyfinder:::eq.LL4.LL4, 
               `L.4 L.4` = synergyfinder:::eq.L4.L4, `LL.4 L.4` = synergyfinder:::eq.LL4.L4,
               `L.4 LL.4` = synergyfinder:::eq.L4.LL4)
  x <- max(drug.col.par[2], drug.row.par[2]) + 1
  for (i in seq_len((ncol(response.mat)) )) {
    for (j in seq_len((nrow(response.mat)) )) {
      x1 <- as.numeric(colnames(response.mat))[i]
      x2 <- as.numeric(rownames(response.mat))[j]
      options(warn = -1)
      slv <- tryCatch({
        slv <- nleqslv::nleqslv(x, eq, method = "Newton", 
                                x1 = x1, x2 = x2, drug.col.par = drug.col.par, 
                                drug.row.par = drug.row.par)
      }, error = function(e) {
        slv <- list(termcd = 999)
      })
      if (slv$termcd < 3) {
        y.loewe <- slv$x
      }
      else {
        y.loewe1 <- synergyfinder:::fun(x1, x2, drug.par = drug.col.par, 
                        model = drug.col.type)
        y.loewe2 <- synergyfinder:::fun(x1, x2, drug.par = drug.row.par, 
                        model = drug.row.type)
        y.loewe <- max(y.loewe1, y.loewe2)
      }
      loewe.mat[j, i] <- ifelse(y.loewe > 100, 
                                        100, y.loewe)
    }
  }
  synergy.mat <- response.mat - loewe.mat
  return(mean(synergy.mat))
  options(warn = 0)
  gc()
}

# # A test to check that the above modification is indeed working
# 
# data("mathews_screening_data")
# msd <- data.table(mathews_screening_data)
# 
# data <- ReshapeData(mathews_screening_data)
# response.mat <- data$dose.response.mats[[1]]
# Loewe(response.mat)
# mean(Loewe(response.mat)[-1,-1])
# max(Loewe(response.mat)[-1,-1])
# min(Loewe(response.mat)[-1,-1])
# 
# synergy.score <- CalculateSynergy(data = data, method = "Loewe", adjusted = FALSE)
# PlotSynergy(synergy.score, type = "all", pair.index = 1)
# 
# Loewe2( msd[block_id==1&conc_r>0&conc_c>0,],
#         msd[block_id==1&conc_c==0,.(response=100-response,dose=conc_r)],
#         msd[block_id==1&conc_r==0,.(response=100-response,dose=conc_c)] )
