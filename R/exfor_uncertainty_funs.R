#  exforUncertainty - Uncertainties for EXFOR Entries
#  Copyright (C) 2019  Georg Schnabel
#  
#  exforUncertainty is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  exforUncertainty is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>

# DATA-ERR and ERR-T never present together
# assume DATA-ERR and ERR-S is statistical uncertainty

# recipe: introduce all ERR-x uncertainties as fully correlated uncertainties, except ERR-S
# if there is a remainder to ERR-T, add the difference as an additional uncorrelated error contribution

.datatable.aware = TRUE

#' Check if EXFOR subentry has valid uncertainty info
#' 
#' @param subentList a list with EXFOR subentries
#' @param quiet should error information be printed
#'
#' @return boolean vector indicating whether elements in \code{subentList}
#'         contain valid uncertainty information.
#'
#' @export
hasValidUncertainties <- function(subentList, quiet=FALSE) {
  
  res <- rep(FALSE, length(subentList))
  
  for (i in seq_along(subentList)) {
    
    subent <- subentList[[i]]
    
    if (all(c("DATA-ERR","ERR-S") %in% subent$DATA$DESCR)) {
      cat(paste0("In entry ", subent$ID, " : both DATA-ERR and ERR-S present\n"))
      next
    }
    # a measurement without an error estimation is no measurement ! 
    if (!any(c("DATA-ERR","ERR-S","ERR-T") %in% subent$DATA$DESCR)) {
      cat(paste0("In entry ", subent$ID, " : DATA-ERR, ERR-S, or ERR-T must be present\n"))
      next
    }
    
    # all error checks passed
    res[i] <- TRUE
  }
  res
}


#' Extract systematic uncertainties from EXFOR subentry
#'
#' @param subent an EXFOR subentry
#' @param dataref reference cross sections used for relative uncertainties.
#'                if NULL (default) equals the cross sections in the subentyr. 
#'
#' @return a datatable with the columns \code{ERRTYPE = (sys-rel, sys-abs)},
#'         \code{UNC} and \code{INFO}. 
#'
#' @export
#'
getSystematicUncertainty <- function(subent, dataref = NULL) {
  
  # extract uncertainties
  sysErrIdcs <- grep("ERR-[0-9]+", subent$DATA$DESCR)
  sysUncLabels <- subent$DATA$DESCR[sysErrIdcs]
  numSysErr <- length(sysErrIdcs)
  
  dataIdx <- which("DATA" == subent$DATA$DESCR)
  origDataVals <- subent$DATA$TABLE[[dataIdx]]
  dataVals <- if (!is.null(dataref)) dataref else origDataVals
  dataUnit <- subent$DATA$UNIT[dataIdx]
  if (dataUnit != "MB")
    stop(paste0("Unit of DATA column is not MB for entry ", subent$ID))
  
  parVals <- numeric(numSysErr)
  errType <- character(numSysErr)
  
  # first, systematic uncertainties
  for (i in seq_len(numSysErr)) {
    curErrIdx <- sysErrIdcs[i]
    curUncLabel <- sysUncLabels[i]
    curUncUnit <- subent$DATA$UNIT[curErrIdx]
    
    if (curUncUnit == "PER-CENT") {
      tmp <- subent$DATA$TABLE[[curErrIdx]]
      if (!diff(range(tmp)) < 1e-1) {
        stop("Some problem with percent scale ", subent$ID, "\n")
      }
      errType[i] <- "sys-rel"
      parVals[i] <- max(tmp) / 100
    }
    else if (curUncUnit == "MB") {
      tmpAbs <- subent$DATA$TABLE[[curErrIdx]]
      tmpRel <- subent$DATA$TABLE[[curErrIdx]] / dataVals
      if (all(tmpAbs == tmpAbs[1])) {
        errType[i] <- "sys-abs"
        parVals[i] <- tmpAbs[1]
      }
      else if (diff(range(tmpRel)) < 1e-1) { # careful with double comparison
        errType[i] <- "sys-rel"
        parVals[i] <- max(tmpRel) / 100
      }
      else {
        stop("Some problem with normalization uncertainty", subent$ID, "\n")
      }
    }
    else {
      stop(paste0("Cannot handle other units than MB. Error appeared in column ",
                  curUncLabel, " of entry ", subent$ID))
    }
  }
  
  if (length(parVals) > 0)
    data.table(DIDX = 1L,
               ERRTYPE = errType,
               BLOCKSIZE = 1L,
               BLOCKID = NA_integer_,
               DATA = 0,
               UNC = parVals,
               INFO = sysUncLabels)
  else NULL
}

#' Extract statistical uncertainty from EXFOR subentry
#'
#' @param subent the EXFOR subentry
#' @param dataref the reference cross section used for relativce uncertainties.
#'        If \code{NULL} (default), the cross sections in the subentry are used
#'        as the reference.
#' @param debias.statunc a correction factor that can be applied to a cross
#'                       section uncertainty given on a relative scale.
#'
#' @return a datatable with columns 
#'         \tabular{ll}{
#'             DIDX    \tab line in subentry table associated with the uncertainty \cr
#'             ERRTYPE \tab either sys-abs or sys-rel \cr
#'             UNC     \tab value of the uncertainty \cr
#'             INFO    \tab additional explanation regarding the uncertainty 
#'         }
#' @export
#'
getStatisticalUncertainty <- function(subent, dataref = NULL, debias.statunc = TRUE) {
  
  dataIdx <- which("DATA" == subent$DATA$DESCR)
  dataErrIdx <- which("DATA-ERR" == subent$DATA$DESCR)
  statErrIdx <- which("ERR-S" == subent$DATA$DESCR)
  
  stopifnot(length(dataErrIdx) <= 1)
  stopifnot(length(statErrIdx) <= 1)
  stopifnot(length(dataIdx) == 1)
  stopifnot(length(dataErrIdx) == 0 || length(statErrIdx) == 0)
  stopifnot(subent$DATA$UNIT[dataIdx] == "MB")
  
  numStatErr <- nrow(subent$DATA$TABLE)
  origDataVals <- subent$DATA$TABLE[[dataIdx]]
  dataVals <- if (!is.null(dataref)) dataref else origDataVals
  # if (length(dataVals) != length(origDataVals)) browser()  # debug
  statErrs <- rep(0, numStatErr)
  statInfoStr <- character(numStatErr)
  statErrType <- rep("stat-abs", numStatErr)
  
  if (length(dataErrIdx) == 1) statErrIdx <- dataErrIdx
  
  if (length(statErrIdx) == 1) { 
    statErrs <- subent$DATA$TABLE[[statErrIdx]]
    statUnit <- subent$DATA$UNIT[statErrIdx]
    if (statUnit == "PER-CENT") {
      statErrs <- dataVals * statErrs / 100
    }
    else if (statUnit == "MB") {
      debiasFactor <- if(debias.statunc) dataVals / origDataVals else 1
      statErrs <- debiasFactor * subent$DATA$TABLE[[statErrIdx]]
    }
    else 
      stop(paste0("Unrecognized unit '", statUnit, "' for statistical error in column ", 
                  subent$DATA$DESCR[statErrIdx]))
    
    data.table(DIDX = seq_len(numStatErr),
               ERRTYPE = "stat-abs",
               BLOCKSIZE = 1L,
               BLOCKID = NA_integer_,
               DATA = 0,
               UNC = statErrs,
               INFO = "")
  }
  else
    NULL
}

#' Automatically construct uncertainties for an EXFOR entry
#'
#' Retrieves the uncertainty information in an EXFOR entry 
#' and applies some rules to construct uncertainties based 
#' on the information in the EXFOR subentry.
#' For instance, inconsistent uncertainty information is
#  penalized by the introduction of additional uncertainty.
#'
#' @param subent  the EXFOR subentry
#' @param sysDt   a datatable as obtained by 
#'                function \code{\link{getSystematicUncertainty}}
#' @param statDt  a datatable as obtained by
#'                function \code{\link{getStatisticalUncertainty}}
#' @param dataref reference cross sections used for relative uncertainties.
#'                If \code{NULL} (default) the cross sections in the subentry
#'                serve as the reference.
#'
#' @return a datatable of the same format as 
#'         \code{\link{getSystematicUncertainty}}
#'
#' @export
#'
getAutoUncertainty <- function(subent, sysDt, statDt, dataref = NULL) {
  
  # sum up normalization uncertainties
  sysErrs <- NULL
  statErrs <- NULL
  totErrs <- NULL
  addSysErr <- NULL
  addStatErr <- NULL
  
  if (!is.null(sysDt)) {
    relIdx <- which(sysDt$ERRTYPE == "sys-rel")
    tmp <- tcrossprod(subent$DATA$TABLE$DATA, sysDt$UNC[relIdx])
    sysErrs1 <- sqrt(colSums(tmp^2))
    absIdx <- which(sysDt$ERRTYPE == "sys-abs")
    sysErrs2 <- sqrt(sum(sysDt$UNC[absIdx]^2))
    sysErrs <- sqrt(sysErrs1^2 + sysErrs2^2)
  }
  
  if (!is.null(statDt))
  {
    stopifnot(all(statDt$ERRTYPE == "stat-abs"))
    statErrs <- statDt$UNC
  }
  
  dataIdx <- which(subent$DATA$DESCR == "DATA")
  if (length(dataIdx) == 0) browser()
  stopifnot(length(dataIdx) == 1)
  dataVals <- if (!is.null(dataref)) dataref else subent$DATA$TABLE[[dataIdx]]
  absDataVals <- abs(dataVals)
  
  totErrIdx <- which(subent$DATA$TABLE$DESCR == "ERR-T")
  stopifnot(length(totErrIdx) <= 1)
  if (length(totErrIdx) == 1)
    totErrs <- subent$DATA$TABLE[[totErrIdx]]
  
  # rules
  
  if (is.null(statErrs) && !is.null(totErrs) && !is.null(sysErrs)) {
    
    addStatErr <- sqrt(pmax(totErrs^2 - sysErrs^2, 0))
    if (any(addStatErr == 0))
      addStatErr <- absDataVals * 0.1
  }
  
  if (is.null(sysErrs) && !is.null(totErrs) && !is.null(statErrs)) {
    tmp <- sqrt(pmax(totErrs^2 - statErrs^2, 0))
    if (any(tmp) == 0)
      addSysErr <- 0.1
    else
      addSysErr <- max(tmp / absDataVals)
  }
  
  if (is.null(statErrs) && is.null(sysErrs) && !is.null(totErrs)) {
    addSysErr <- 0.1
    addStatErr <- absDataVals * 0.1
  }
  
  if (!is.null(statErrs) && !is.null(sysErrs) && !is.null(totErrs)) {
    diffUnc <- totErrs^2 - sysErrs^2 - statErrs^2
    tmp <- sqrt(abs(diffUnc))
    addSysErr <- max(abs / absDataVals)
    addStatErr <- abs(diffUnc)
  }
  
  if (is.null(totErrs) && is.null(statErrs)) {
    addStatErr <- absDataVals * 0.1
  }
  
  if (is.null(totErrs) && is.null(sysErrs)) {
    addSysErr <- 0.1
  }
  
  # making the datatable
  addStatDt <- NULL
  addSysDt <- NULL
  
  if (!is.null(addStatErr))
    addStatDt <- data.table(DIDX = seq_along(absDataVals),
                            ERRTYPE = "stat-abs",
                            BLOCKSIZE = 1L,
                            BLOCKID = NA_integer_,
                            DATA = 0,
                            UNC = addStatErr,
                            INFO = "auto-added")
  
  if (!is.null(addSysErr))
    addSysDt <- data.table(DIDX = 0L,
                           ERRTYPE = "sys-rel",
                           BLOCKSIZE = 1L,
                           BLOCKID = NA_integer_,
                           DATA = 0,
                           UNC = addSysErr,
                           INFO = "auto-added")
  
  rbind(addStatDt, addSysDt)
}


#' Constructs Uncertainties for a list of EXFOR subentries
#'
#' This function takes a datatable \code{expDt} which specifies
#' a selection of data points from different subentries
#' and returns a datataable with associated uncertainties.
#'
#' @param expDt       datatable which specifies a selection of data points in
#'                    the EXFOR database
#' @param subentList  a list of subentries that must contain the subentries 
#'                    indexed by rows in \code{expDt}
#' @param dataref.col name of the column in \code{expDt} containing the 
#'                    reference cross sections used as the basis for
#'                    relative uncertainties
#'
#' @note The returned uncertainties do not necessarily reflect
#'       one-to-one those in the EXFOR subentry. If uncertainties
#'       in a subentry are inconsistent, e.g., ERR-T does not match
#'       with individual contributions, extra uncertainties are
#'       introduced.
#'         
#' @return A datatable with uncertainty specifications.
#'         Relevant columns are \code{DIDX, ERRTYPE, UNC, INFO}.
#'
#' @export
#'
getAvailableUncertainties <- function(expDt, subentList, dataref.col = NULL) {
  
  # sanity checks
  subentIds <- sapply(subentList, function(x) x$ID)
  subentAvailable <- expDt$EXPID %in% subentIds
  if (!all(subentAvailable))
    stop(paste0("SubentList misses the needed subentries with IDs: ",
                paste0(expDt$EXPID[!subentAvailable], collapse=",")))
  
  neededSubent <- subentIds %in% expDt$EXPID
  validUncFlag <- hasValidUncertainties(subentList[neededSubent]) 
  if (!all(validUncFlag)) {
    stop(paste0("Subentries with invalid or missing uncertainty information, IDs: ",
                paste0(neededSubent[validUncFlag], collapse=",")))
  }
  
  dupRowIdx <- anyDuplicated(expDt[,list(EXPID,DIDX)])
  if (dupRowIdx > 0) {
    stop(paste0("expDt contains duplicated datasets, e.g., row ", dupRowIdx,
                " is a duplicate"))
  }
  
  # correlations within datasets
  resDt <- expDt[,{
    
    # retrieve relevant subentry
    subentIdx <- which(subentIds == EXPID)
    stopifnot(length(subentIdx) == 1)
    subent <- subentList[[subentIdx]]
    
    cat(paste0("PROCESSING ", subent$ID, "\n"))
    
    dataref <- rep(NA_real_, nrow(subent$DATA$TABLE))
    dataref[DIDX] <- if (!is.null(dataref.col))
      get(dataref.col, inherits = FALSE)
    else NULL
      
    curSysDt <- getSystematicUncertainty(subent, dataref = dataref)
    curStatDt <- getStatisticalUncertainty(subent, dataref = dataref)
    curAutoDt <- getAutoUncertainty(subent, curSysDt, curStatDt, dataref = dataref)
    
    rbind(curSysDt, curStatDt, curAutoDt)
  }, by="EXPID"]
  
  selIdx <- which(resDt[, ERRTYPE == "sys-rel" | ERRTYPE == "sys-abs"])
  resDt[selIdx, BLOCKID := seq_along(selIdx)]
  resDt[]
}

#' Merge uncertainty components 
#'
#' A datatable with uncertainties as returned by
#' \code{\link{getAvailableUncertainties}} may contain
#' several systematic uncertainty components for the 
#' same experiment. This function merges such components
#' to only one for more efficient statistical inference.
#'
#' @param uncDt a datatable with uncertainty specifications as obtained by 
#'        function \code{\link{getAvailableUncertainties}}
#'
#' @return a datatable of the same structure as \code{uncDt} 
#'         with a reduced number of systematic components
#'
#' @export
#'
compactifyUncDt <- function(uncDt) {
  
  stopifnot(all(uncDt[, ERRTYPE == "sys-rel" | ERRTYPE == "sys-abs" | ERRTYPE == "stat-abs"]))
  uncDt[, {
    list(BLOCKSIZE = 1L,
         BLOCKID = BLOCKID[1],
         DATA = 0,
         UNC = sqrt(sum(UNC^2)),
         INFO = "")
  }, by = c("EXPID", "DIDX", "ERRTYPE")]
}


#' Split uncertainty datatable
#'
#' Splits the uncertainty information in a datatable
#' as obtained by \code{\link{getAvailableUncertainties}}
#' into two datatables: one with systematic components
#' and another one with statistical components.
#'
#' @param uncDt a datatable with uncertainty information
#'
#' @return a list with two datatables for systematic and
#'         statistical uncertainties, respectively.
#'
#' @export
#'
splitUncDt <- function(uncDt) {
  
  statUncDt <- uncDt[ERRTYPE == "stat-abs",]
  sysUncDt <- uncDt[ERRTYPE %in% c("sys-rel", "sys-abs")]
  list(statUncDt = statUncDt,
       sysUncDt = sysUncDt)
}


#' Add statistical uncertainty to expDt
#'
#' @param expDt      a datatable which specifies a selection of
#'                   data points in the EXFOR database
#' @param statUncDt  a datatable with statistical uncertainties
#'
#' @note IMPORTANT: The datatable \code{expDt} given as input
#'                  argument is modified by this function!
#'
#' @return The modifed \code{expDt}, which is in principle not
#'         necessary because \code{expDt} is changed by reference.
#'
#' @export
#'
addStatUncToExpDt <- function(expDt, statUncDt) {

  stopifnot(all(statUncDt$ERRTYPE == "stat-abs"))  
  setkey(statUncDt, EXPID, DIDX)
  x1 <- expDt[, EXPID]
  x2 <- expDt[, DIDX]
  expDt[, UNC := statUncDt[J(x1, x2), UNC]]
  expDt[]
}





