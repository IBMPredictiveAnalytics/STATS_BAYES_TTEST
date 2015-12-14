#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2015
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "IBM SPSS, JKP"
# version__ = "1.0.0"

# History
# 23-sep-2015 Original Version


gtxt <- function(...) {
    return(gettext(...,domain="STATS_BAYES_TTEST"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_BAYES_TTEST"))
}

### MAIN ROUTINE ###
doBayesttest = function(variables, group=NULL, paired=FALSE, testvalue=NULL, 
    index=TRUE, rscalecont="medium", iterations=1000,
    modelfile=NULL, workspaceaction="clear", modelfileout=NULL) {
    # Estimate Bayes t test
    
    # The modelsource and modelfile
    # parameters are not implemented, awaiting requests for that functionality

    setuplocalization("STATS_BAYES_TTEST")
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Bayesian T Test")
    warningsprocname = gtxt("Bayesian T TEST: Warnings")
    omsid="STATSBAYESTTEST"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    tryCatch(library(BayesFactor), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.", "BayesFactor"),dostop=TRUE)
        }
    )
    if (!is.null(spssdictionary.GetWeightVariable())) {
        warns$warn(
            gtxt("The dataset is weighted, but case weights are not used in this procedure except for screening out cases with a nonpositive weight"),
            dostop=FALSE)
    }
    if (!is.null(spssdata.GetSplitVariableNames())) {
        warns$warn(
            gtxt("Split variables are not honored by this procedure"),
            dostop=FALSE)
    }
    # For unpaired t test, the data is a single variable (VARIABLES) and a grouping variable (GROUP),
    # which must have exactly two values ignoring missing values
    # or
    # a single variable and a test value (TESTVALUE)
    
    # For paired t test (PAIRED=YES), the data is two variables
    if (!paired) {
        if (length(variables) > 1 || (is.null(group) && is.null(testvalue))) {
            warns$warn(gtxt("For the unpaired t test, one test variable and a group variable or test value must be specified"),
                dostop=TRUE)
        }
        if (!is.null(group) && !is.null(testvalue)) {
            warns$warn(gtxt("Cannot specify both a group variable and a test value"),
                dostop=TRUE)
        }
    } else {
        if (length(variables) != 2 || !is.null(group)) {
            warns$warn(gtxt("For the paired t test, two test variables must be specified and no group variable"),
                dostop=TRUE)
        }
        if (!is.null(testvalue)) {
            warns$warn(gtxt("For the paired t test, a test value cannot be specified"),
                dostop=TRUE)
        }
    }

    alldata = c(variables, group)
    rscalelist = scales(rscalecont, warns)
    allargs = as.list(environment())
    # with only one variable, value is not a data frame
    dta = spssdata.GetDataFromSPSS(alldata, missingValueToNA=TRUE,
        factorMode="levels")
    if (any(as.logical(lapply(dta[, unlist(variables)],is.factor)))) {
        warns$warn(gtxt("Categorical variables cannot be tested in this procedure"),
            dostop=TRUE)
    }

    # The procedure does not allow missing values
    allargs$ncases = nrow(dta)
    # If the data frame has only one column, the extraction of
    # complete cases loses the data frame class :-(
    if (length(dta) == 1) {
        dta[2] = 1:allargs$ncases
    }
    dta = dta[complete.cases(dta),]
    allargs$nvalid = nrow(dta)
    
    reslist = list()
    tryCatch(
        {for (rindex in 1:length(rscalelist)) {
            if (!paired) {
                if (!is.null(testvalue)) {   # one variable plus a test value
                    reslist[rindex] = ttestBF(dta[,variables], mu=testvalue, rscale=rscalelist[[rindex]])
                } else {  # one variable plus grouping variable with exactly two values
                    tblrows = row.names(table(dta[[2]]))
                    if (length(tblrows) != 2) {
                        warns$warn(gtxt("The grouping variable must have exactly two non-missing values"),
                            dostop=TRUE)
                    }
                    reslist[rindex] = ttestBF(x=dta[1][dta[2] == tblrows[[1]]],
                                  y=dta[1][dta[2] == tblrows[[2]]], rscale=rscalelist[[rindex]])
                }
            } else {
                reslist[rindex] = ttestBF(x=dta[[1]], y=dta[[2]], paired=TRUE, rscale=rscalelist[[rindex]])
            }
        
        }
        }, error=function(e) {
            warns$warn(e, dostop=TRUE)
        }
    )
    post = doposterior(allargs, reslist[[1]], warns)

    displayresults(allargs, reslist, post, warns)
    
    if (!is.null(modelfileout)) {
        save(allargs, reslist, post, file=modelfileout)
    }
    if (workspaceaction == "retain") {
        assign("allargs", allargs, envir=.GlobalEnv)
        assign("reslist", reslist, envir=.GlobalEnv)
        assign("post", post, envir=.GlobalEnv)
    }
    warns$display()
}

doposterior = function(allargs, res, warns) {
    # calculate posterior distribution if model index specified
    
    if (!allargs$index) {
        return(NULL)
    }

    post = posterior(model=res, index=1, iterations=allargs$iterations,
        progress=FALSE)
    return(post)
}

slist = list("medium" = .7071, 'wide'=1.0, "ultrawide"=1.4142)
scales = function(scale, warns) {
    # return a numeric list of scale items, resolving certain strings
    
    if (is.null(scale)) {
        return(.7071)
    }

    s = tolower(strsplit(scale, split=" +")[1][[1]])
    if (s[[1]] == "") {
        s = s[-1]
    }
    sss = slist[s]
    for (index in 1:length(s)) {
        if (is.null(sss[[index]])) {
            sss[index] = as.numeric(s[index])
        }
    }
    if (any(is.na(sss)) || any(as.logical(lapply(sss, function(x) x <= 0))) ) {
        warns$warn(gtxt("An invalid value was given for the prior scale"),
            dostop=TRUE)
    }
    return(sss)
}

    
scaletrans=list("medium"=gtxtf("medium (-%s, %s)", .7071, .7071), 
        "wide"=gtxtf("wide (-%s, %s)", 1., 1.), 
        "ultrawide"=gtxtf("ultrawide (%s, %s)", 1.414, 1.414)
)
waction=list("clear"="clear", "retain"="retain")
outputmap = list("delta"=gtxt("Effect Size"), 
                 "mu"=gtxt("Mean"),
                 "beta (x - y)"=gtxt("Difference"),
                 "sig2"=gtxt("Sigma2"),
                 "g" = gtxt("g")
)

displayresults = function(allargs, reslist, post, warns) {
    # display results
    # allargs is the parameter set

    StartProcedure(allargs[["procname"]], allargs[["omsid"]])
    
    # summary results
    # input specifications
    # although groups can be specified (cengroup), separate results are not
    # produced.
    lbls = c(gtxt("Test Variables"),
             gtxt("Group Variable"),
             gtxt("Test Value"),
             gtxt("Paired"),
             gtxt("Number of Cases"),
             gtxt("Number of Valid Cases"),
             gtxt("Prior Scale on Standardized Slopes"),
             gtxt("Posterior Iterations"),
             gtxt("Workspace Action"),
             gtxt("Output Model File")
    )

    vals = c(
            paste(allargs$variables, collapse=" "),
            ifelse(is.null(allargs$group), 
                gtxt("--NA--"), allargs$group),
            ifelse(is.null(allargs$testvalue), gtxt("--NA--"), allargs$testvalue),
            ifelse(allargs$paired, gtxt("Yes"), gtxt("No")),
            allargs$ncases,
            allargs$nvalid,
            paste(allargs$rscalelist, collapse=" "),
            ifelse(!allargs$index, gtxt("--NA--"), allargs$iterations),
            waction[allargs$workspaceaction],
            ifelse(is.null(allargs$modelfile), gtxt("--NA--"), allargs$modelfile)
    )

    spsspivottable.Display(data.frame(cbind(vals), row.names=lbls), title = gtxt("Summary"),
        collabels=c(gtxt("Summary")), templateName="BAYESTTESTSUMMARY", outline=gtxt("Bayes T TEST Summary"),
        caption = gtxtf("Computations done by R package BayesFactor, version: %s", packageVersion("BayesFactor"))
    )
    # iterate over all the individual rscale results and produce a merged pivot table
    iter = 1
    for (res in reslist) {
    ressum = extractBF(res)

    bf = data.frame(seq(1: length(res)),ressum[1:2])
    bf[3] = bf[3] * 100.
    bf = data.frame(bf, length(res) - rank(bf[2]) + 1)
    # add in posterior probabilities excluding Intercept only model
    ###postprob = data.frame(as.BFprobability(newPriorOdds(res) * res))[-(nrow(bf)+1), 1]
    
    # construct posterior probabilities and merge with Bayes factors
    # The order for probabilities may not be the same as for the Bayes factors
    # which requires a few extra steps to get things merged
    # the BF data frame may not have the intercept row, so that row may be discarded
    postprob = data.frame(as.BFprobability(newPriorOdds(res) * res))[1]

    bf = merge(bf, postprob, by="row.names")
    bf = bf[order(bf[[2]]),]
    row.names(bf) = bf[["Row.names"]]
    bf = bf[-1]
    
    bf = bf[c(2,3,5)]

    names(bf) = c(
        gtxt("Bayes Factor"), gtxt("Error (+-%)"),
        gtxt("Posterior Probabilities (Equal Prior)"))
    row.names(bf) = c(allargs$rscalelist[iter])

    spsspivottable.Display(bf,
        title=gtxt("Bayes Factors"),
        rowdim=gtxt("Prior Scale"), 
        hiderowdimtitle=FALSE,
        templateName="BAYESTTESTFACTORS",
        outline=gtxt("Bayes Factors")
    )
    iter = iter + 1
    }

    if (allargs$index) {
        postsum = summary(post)
        postsumstats = data.frame(postsum$statistics[,-4])  # omit time series SEs
        names(postsumstats) = c(gtxt("Mean"), gtxt("Std. Deviation"), gtxt("SE Mean"))
        row.names(postsumstats) = outputmap[row.names(postsumstats)]
        spsspivottable.Display(
            postsumstats, 
            title=gtxtf("Posterior Summary Statistics (Prior scale: %s)", allargs$rscalelist[[1]]),
            rowdim=gtxt("Statistic"),
            hiderowdimtitle=TRUE,
            templateName="BAYESTTESTPOSTSTATS",
            outline=gtxt("Posterior Summary Statistics")
        )
        
        postsumquant = postsum$quantiles
        row.names(postsumquant) = outputmap[row.names(postsumquant)]
        spsspivottable.Display(
            postsumquant,
            title=gtxtf("Posterior Quantiles (Prior scale: %s)", allargs$rscalelist[[1]]),
            rowdim=gtxt("Variables"),
            hiderowdimtitle=TRUE,
            templateName="BAYESTTESTPOSTQUANTILES",
            outline=gtxt("Posterior Quantiles")
        )
    }

    spsspkg.EndProcedure()
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = mylist2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spsspkg.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

mylist2env = function(alist) {
    env = new.env()
    lnames = names(alist)
    for (i in 1:length(alist)) {
        assign(lnames[[i]],value = alist[[i]], envir=env)
    }
    return(env)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 
# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}



Run = function(args) {
    #Execute the STATS COMPRISK command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("VARIABLES", subc="", ktype="existingvarlist", var="variables", islist=TRUE),
        spsspkg.Template("GROUP", subc="", ktype="existingvarlist", var="group"),
        spsspkg.Template("PAIRED", subc="", ktype="bool", var="paired"),
        spsspkg.Template("TESTVALUE", subc="", ktype="float", var="testvalue"),
#         spsspkg.Template("MODELSOURCE", subc="", ktype="str", var="modelsource",
#             vallist=list("none", "workspace", "modelfile")),
#         spsspkg.Template("MODELFILE", subc="", ktype="literal", var="modelfile"),
        
        ###spsspkg.Template("COMPARISON", subc="OPTIONS", ktype="int", var="comparison"),
        ###spsspkg.Template("PLOTMODELS", subc="OPTIONS", ktype="bool", var="plotbayesf"),
        spsspkg.Template("POSTERIOR", subc="OPTIONS", ktype="bool", var="index"),
        spsspkg.Template("ITERATIONS", subc="OPTIONS", ktype="int", var="iterations",
            vallist=list(2)),
        spsspkg.Template("PRIORSCALE", subc="OPTIONS", ktype="literal", var="rscalecont"),
        spsspkg.Template("WORKSPACE", subc="SAVE", ktype="str", var="workspaceaction",
            vallist=list("retain", "clear")),
        spsspkg.Template("MODELFILE", subc="SAVE", ktype="literal", var="modelfileout")
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "doBayesttest")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
