rm(list=ls())
graphics.off()

library(reshape)
library(nlme)
## library(dplyr)
library(multcomp)
library(effects)

source("scoreMasc.r")

########################################################################################################################################################################################################
### START OF FUNCTIONS #################################################################################################################################################################################
########################################################################################################################################################################################################


## http://wiki.stdout.org/rcookbook/Graphs/Plotting%20means%20and%20error%20bars%20(ggplot2)/
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    require(plyr)

    ## New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }
    
    ## This is does the summary; it's not easy to understand...
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun= function(xx, col, na.rm) {
                       c(N    = length2(xx[,col], na.rm=na.rm),
                         mean = mean (xx[,col], na.rm=na.rm),
                         median = median (xx[,col], na.rm=na.rm),
                         IQR = IQR (xx[,col], na.rm=na.rm),                         
                         mad = mad (xx[,col], na.rm=na.rm),                       
                         sd   = sd   (xx[,col], na.rm=na.rm),
                         min  = min  (xx[,col], na.rm=na.rm),
                         max  = max  (xx[,col], na.rm=na.rm),                       
                         nacount  = sum  (is.na((xx[,col])))
                         )
                   },
                   measurevar,
                   na.rm
                   )
    
    ## Rename the "mean" column    
    datac <- rename(datac, c("mean"=measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  ## Calculate standard error of the mean
    
    ## Confidence interval multiplier for standard error
    ## Calculate t-statistic for confidence interval: 
    ## e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
}


readCsvFile <- function (inFilename, inSubjectColumnName="ID") {

    cat("*** Reading", inFilename, "\n")
    rCsv=read.csv(inFilename, header=T, na.strings = c(
                                            "NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "#NAME?", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", ".", "",
                                            "Test_Subject_", "Test_Subject", "NO SUBJECT", " NO SUBJECT", "NO SUBJECT ", " NO SUBJECT ", "Cancelled",
                                            "Head too big", "IGNORE THIS SUBJ", "NOT ELIGIBLE", "Permanent retainer", "PNS, fibromyalgia", "Psychosis/MDD"
    ))

    cat(sprintf("*** Read data for %s unique subjects\n",  length(unique(rCsv[, inSubjectColumnName]))))

    return(rCsv)
}

splitSubjectOrderIntoIdAndTimepoint <- function(inSubjectOrderTable) {

    new.subject.order.table=cbind(inSubjectOrderTable, as.data.frame(t(as.data.frame(strsplit(as.character(inSubjectOrderTable$subject), "_", fixed=TRUE)))))
    rownames(new.subject.order.table)=NULL
    colnames(new.subject.order.table)=c("id", "subject", "timepoint")

    ## new.subject.order.table <- subjectList %>% separate(subject, into=c("ID", "timepoint"), sep="_", remove=FALSE)
    ## colnames(new.subject.order.table)=c("id", "subject", "timepoint")
    
    return(new.subject.order.table)
}

fixSubjectIds <- function (inDataFrame, inColumnName=NULL) {
    if ( is.null(inColumnName) ) {
        stop("Got a NULL column name in fixSubjectIds. Stopping\n")
    }

    inDataFrame[, inColumnName]=as.character(inDataFrame[, inColumnName])
    subjectIds=                 as.character(inDataFrame[, inColumnName])
    if (length(pmatch(inDataFrame[, inColumnName], "300")) > 0 ) {
        inDataFrame[, inColumnName]=sub("^300", "169/300", inDataFrame[, inColumnName], fixed=FALSE)
        ## inDataFrame[, inColumnName]=sub("^0", "", inDataFrame[, inColumnName], fixed=FALSE)
    }
    inDataFrame[, inColumnName]=as.factor(inDataFrame[, inColumnName])
    return (inDataFrame)
}

fixDates <- function (inData) {
    ## this complicated looking regexp stuff cleans up the years in DOB
    ## and MRI with 4 digits to be just 2 digits
    ## month day year
    inData$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", inData$DOB)
    inData$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", inData$MRI)
    
    ## now convert to year/month/day
    inData$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/([0-9]{2})", "\\3/\\1/\\2", inData$DOB)
    inData$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/([0-9]{2})", "\\3/\\1/\\2", inData$MRI)
    
    inData$DOB=as.Date(inData$DOB, "%y/%m/%d")
    inData$MRI=as.Date(inData$MRI, "%y/%m/%d")

    return(inData)
}

computeAge <- function(inData) {

    age.in.weeks=difftime(inData$MRI, inData$DOB, units="weeks")
    age.in.weeks=as.numeric(age.in.weeks)

    inData$age.in.years=age.in.weeks/52

    return(inData)
}

find.usable.subjects <- function (in.data, in.column.name) {

    print (colnames(in.data))
    ss=subset(in.data, select=c("ID", "timepoint", "Grp", "age.in.years", in.column.name))
    rownames(ss)=NULL
    
    ss$is.in.subjectList=ss$ID %in% subjectList$subject

    ss.wide=reshape(ss,
        timevar="timepoint",
        idvar=c("ID", "Grp", "is.in.subjectList", "age.in.years"),
        direction="wide")
    ss.wide$at.both.timepoints=! ( is.na(ss.wide[ , paste(in.column.name, "A", sep=".") ]) | is.na(ss.wide[ , paste(in.column.name, "C", sep=".") ]) )
    rownames(ss.wide)=NULL
    
    usable.ss.wide=subset(ss.wide, is.in.subjectList & at.both.timepoints)
    rownames(usable.ss.wide)=NULL

    usable.ss.long=subset(ss, ss$ID %in% usable.ss.wide$ID)
    
    return (list("long"=ss, "usable"=usable.ss.wide, "usable.long"=usable.ss.long))
}

print.usable.subjects <- function (in.usable.subjects.list) {

    for (gg in levels(in.usable.subjects.list$Grp)) {
        cat (paste (as.vector ( in.usable.subjects.list[in.usable.subjects.list$Grp==gg, "ID"]), collapse=" "), "\n")
        cat (paste (as.vector ( in.usable.subjects.list[in.usable.subjects.list$Grp==gg, "Grp"    ]), collapse=" "), "\n")        
    }
}


save.data.table <- function (in.data.table, in.variable, in.seed.name, in.timepoint=NULL, in.file.type="csv") {

    ## ncl.list.filename=file.path(admin.data.dir, sprintf("new.ncl.%s.change.score.csv", in.variable))
    ## mdd.list.filename=file.path(admin.data.dir, sprintf("new.mdd.%s.change.score.csv", in.variable))

    file.types=c("csv", "tab")
    separators=c(",", "\t")
    file.extension=pmatch(in.file.type, file.types)    
    if (! is.na(file.extension)) {
        separator=separators[file.extension]
        file.extension=file.types[file.extension]
    } else {
        stop(sprintf("*** Unknown file extension %s. valid extensions are csv and tab\n", in.file.type)) 
    }
    
    ncl.list.filename=file.path(group.results.dir, sprintf("data.table.ncl.%s.%s.%s", in.variable, in.seed.name, file.extension))
    mdd.list.filename=file.path(group.results.dir, sprintf("data.table.mdd.%s.%s.%s", in.variable, in.seed.name, file.extension))

    ncl.list=in.data.table[in.data.table$Grp=="NCL", ]
    mdd.list=in.data.table[in.data.table$Grp=="MDD", ]
    rownames(ncl.list)=NULL
    rownames(mdd.list)=NULL    
    
    ncl.list$Grp=NULL
    mdd.list$Grp=NULL
    
    if (! is.null(in.timepoint) ) {
        ncl.list$Subj = paste(ncl.list$Subj, in.timepoint, sep="_")
        mdd.list$Subj = paste(mdd.list$Subj, in.timepoint, sep="_")        
    }
   
    cat("*** Writing NCL list to" , ncl.list.filename, "\n")
    ## print(ncl.list)
    write.table(ncl.list, ncl.list.filename, quote=FALSE, col.names=TRUE, row.names=FALSE, sep=separator)
    
    cat("*** Writing MDD list to" , mdd.list.filename, "\n")
    ## print(mdd.list)
    write.table(mdd.list, mdd.list.filename, quote=FALSE, col.names=TRUE, row.names=FALSE, sep=separator)
}


change.score.analyse <- function (in.subjects.lists, in.sl, y.axis, na.action="na.omit") {

    cat("####################################################################################################\n")
    cat("### Analysis of", y.axis, "\n")

    cat("### Summary of", y.axis, "\n")
    print(summarySE(in.subjects.lists[["usable.long"]], measure=y.axis, groupvars=c("Grp")))
    print(summarySE(in.subjects.lists[["usable.long"]], measure=y.axis, groupvars=c("timepoint")))
    smry=summarySE(in.subjects.lists[["usable.long"]], measure=y.axis, groupvars=c("Grp", "timepoint"))
    print(smry)
    cat(make.table.strings(smry, y.axis))
    cat("\n")

    y.axis=paste(y.axis, "diff", sep=".")
    ss=summarySE(in.sl, measure=y.axis, groupvars=c("Grp"))
    print(ss)
    gg="MDD"
    mdd.line.entry=sprintf("%s: \"%0.1f ± %0.1f (%0.1f, %0.1f)\"",
        gg, 
        ss[ss$Grp==gg, y.axis],
        ss[ss$Grp==gg, "sd"],
        ss[ss$Grp==gg, "min"],
        ss[ss$Grp==gg, "max"])
    gg="NCL"
    ncl.line.entry=sprintf("%s: \"%0.1f ± %0.1f (%0.1f, %0.1f)\"",
        gg, 
        ss[ss$Grp==gg, y.axis],
        ss[ss$Grp==gg, "sd"],
        ss[ss$Grp==gg, "min"],
        ss[ss$Grp==gg, "max"])
    cat(sprintf("%s, %s\n", mdd.line.entry, ncl.line.entry))
    
    cat("### Table of n\n")
    print(with(in.subjects.lists[["usable.long"]], addmargins(table(Grp, timepoint))))


    model.formula=as.formula(paste(y.axis, "~", "Grp"))
    random.formula=as.formula("~ 1 | ID")
    model.lme=lme(fixed=model.formula, data=in.sl, random=random.formula)

    model.anova=anova(model.lme)

    cat("\n### LME summary\n")
    print(summary(model.lme))
    cat("\n### ANOVA summary\n")
    print(model.anova)
    assign("model.formula", model.formula, envir=.GlobalEnv)
    assign("in.sl", in.sl, envir=.GlobalEnv)
    print(Effect("Grp", model.lme))
    remove("model.formula", envir=.GlobalEnv)
    remove("in.sl", envir=.GlobalEnv)          
}



make.graph <- function(in.summary) {

    my.base.size=12
    my.theme=
        theme_bw(base_size =  my.base.size) +
            theme(
                ## legend.position="none",
                legend.position="bottom",        
                ## panel.grid.major = element_blank(),
                ## panel.grid.minor = element_blank(),
                
                ##remove the panel border
                ## panel.border = element_blank(),
                
                ## add back the axis lines
                axis.line=element_line(colour = "grey50"),
                
                ##axis.title.x=element_blank(),
                axis.title.x = element_text(size=my.base.size, vjust=0),
                axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
                plot.title=element_text(size=my.base.size*1.2, vjust=1))


    my.dodge=position_dodge(.2)
    levels(in.summary$Grp)=c("MDD", "HCL")
    levels(in.summary$timepoint)=c("Baseline", "Follow-up")    
    ## graph=ggplot(data=in.summary, aes(x=Grp, y=CDRS.t.score, group=timepoint, shape=timepoint, color=timepoint)) +
    ##     geom_errorbar(aes(ymin=CDRS.t.score-sd, ymax=CDRS.t.score+sd), width=0.5, position=my.dodge) +
    ##         geom_line(position=my.dodge) +
    ##             geom_point(position=my.dodge, size=3, shape=21, fill="white") +
    ##                 scale_shape_discrete(name="Timepoint:",
    ##                                      breaks=c("A", "C"),
    ##                                      labels=c("Baseline", "3-month follow-up")) +
    ##                                          scale_color_brewer(name="Timepoint:", palette="Set1",
    ##                                                             breaks=c("A", "C"),
    ##                                                             labels=c("Baseline", "3-month follow-up")) +
    ##                                                                 labs(title="Change in CDRS-R Score Over Time",
    ##                                                                      x="Group",
    ##                                                                      y="CDRS-R Score") +
    ##                                                                          my.theme

    graph=ggplot(data=in.summary, aes(x=timepoint, y=CDRS.t.score, group=Grp, shape=Grp, color=Grp)) +
        geom_errorbar(aes(ymin=CDRS.t.score-sd, ymax=CDRS.t.score+sd), width=0.5, position=my.dodge) +
            geom_line(position=my.dodge) +
                geom_point(position=my.dodge, size=3, shape=21, fill="white") +
                    scale_shape_discrete(name="Group:") +
                        scale_color_brewer(name="Group:", palette="Set1") +
                            labs(title="Mean CDRS-R Score",
                                 x="Group",
                                 y="CDRS-R Score") +
                                     my.theme

    
    return(graph)
}

save.graph <- function(in.graph, in.variable) {

    imageFilename=file.path(group.results.dir, sprintf("%s.change.over.time.pdf", in.variable))
    cat(paste("*** Creating", imageFilename, "\n"))

    ggsave(imageFilename, in.graph, width=4, height=3.5, units="in")
    ## ggsave(imageFilename, in.graph)    
}

make.line.entry <- function (in.data, in.variable, in.group, in.timepoint) {

    nacount=in.data[in.data$Grp==in.group & in.data$timepoint==in.timepoint, "nacount"]
    nacount.string=""
    if(nacount > 0) {
        nacount.string=sprintf(" [%d]", nacount)
    }
    
    line.entry=sprintf("%s, %s: \"%0.1f ± %0.1f (%0.1f, %0.1f)%s\"",
        in.group, in.timepoint, 
        in.data[in.data$Grp==in.group & in.data$timepoint==in.timepoint, in.variable],
        in.data[in.data$Grp==in.group & in.data$timepoint==in.timepoint, "sd"],
        in.data[in.data$Grp==in.group & in.data$timepoint==in.timepoint, "min"],
        in.data[in.data$Grp==in.group & in.data$timepoint==in.timepoint, "max"],
        nacount.string)
    
    return(line.entry)
}

make.table.strings <- function (in.data, in.variable) {

    group.by.timepoint=expand.grid(unique(in.data$Grp), unique(in.data$timepoint))

    lines=apply(group.by.timepoint, 1,
        function(xx) {
            make.line.entry(in.data, in.variable, xx[1], xx[2])
        })

    return(lines)
}

compute.change.scores <- function(in.data, in.variable) {

    a.column.name=paste(in.variable, "A", sep=".")
    c.column.name=paste(in.variable, "C", sep=".")

    diff.column.name=paste(in.variable, "diff", sep=".")
    scaled.diff.column.name=paste(in.variable, "scaled.diff", sep=".")    
    percent.diff.column.name=paste(in.variable, "percent.diff", sep=".")
    
    ## ss.change$ss.diff       = ss.change[, paste(col, "C", sep=".")] - ss.change[, paste(col, "A", sep=".")] ## with(ss.change, ss.tscore.C - ss.tscore.A)
    ## ss.change$ss.percentDiff= (ss.change$ss.diff     / ss.change[, paste(col, "A", sep=".")]) * 100    

    ## vectorized code below that doesn't account for possible NaN and
    ## Inf introduced by unfortunate zeros
    ## in.data[, diff.column.name] = ( in.data[, c.column.name] -  in.data[, a.column.name] )
    ## in.data[, scaled.diff.column.name] = in.data[, diff.column.name] / in.data[, a.column.name]
    ## in.data[, percent.diff.column.name] = in.data[, scaled.diff.column.name] * 100

    in.data[ , diff.column.name]        = vector(mode="numeric", length=dim(in.data)[1])
    in.data[ , scaled.diff.column.name] = vector(mode="numeric", length=dim(in.data)[1])
    for (ii in 1:dim(in.data)[1]) {
        ## f = followup (time C) b = baseline (time A)
        f = in.data[ii, c.column.name]
        b = in.data[ii, a.column.name]

        if (f==0 && b==0) {
            diff=0
            scaled.diff=0
        } else {
            diff=f-b
            scaled.diff=diff/b
        }
        in.data[ii, diff.column.name]        = diff
        in.data[ii, scaled.diff.column.name] = scaled.diff        
    }

    return (in.data)
}

compute.residualized.scores <- function(in.data, in.variable) {

    a.column.name=paste(in.variable, "A", sep=".")
    c.column.name=paste(in.variable, "C", sep=".")

    form=paste(c.column.name, "~", a.column.name)
    std.residuals=rstandard(lm(as.formula(form), in.data))

    resid.column.name=paste(in.variable, "rstandard", sep=".")
    in.data[, resid.column.name]=std.residuals

    return (in.data)
}

find.usable.timepoint.a.subjects <- function (in.subjects.lists, in.group) {

    aa=in.subjects.lists$long[subjects.lists$long$Grp==in.group & subjects.lists$long$timepoint=="A" & subjects.lists$long$is.in.subjectList==TRUE, ]
    rownames(aa)=NULL
    return(aa)
}

scale.variable <- function(in.subjects.lists, in.variable) {

    variable.scaled=paste(in.variable, "scaled", sep=".")
    in.subjects.lists[, variable.scaled] = scale(in.subjects.lists[, in.variable], scale=FALSE)

    return(in.subjects.lists)
}    



new.analyse <- function (in.data, y.axis, na.action="na.omit") {

    cat("####################################################################################################\n")
    cat("### Analysis of", y.axis, "\n")

    cat("### Summary of", y.axis, "\n")
    print(summarySE(in.data, measure=y.axis, groupvars=c("Grp")))
    print(summarySE(in.data, measure=y.axis, groupvars=c("timepoint")))
    smry=summarySE(in.data, measure=y.axis, groupvars=c("Grp", "timepoint"))
    print(smry)
    cat(make.table.strings(smry, y.axis))
    cat("\n")

    cat("### Table of n\n")
    print(with(in.data, addmargins(table(Grp, timepoint))))

    ## setup the formulae used in the models
    model.formula.nointeraction=as.formula(paste(y.axis, "~", "Grp + timepoint"))
    model.formula.interaction=as.formula(paste(y.axis, "~", "Grp * timepoint"))    
    random.formula=as.formula("~ 1 | ID")

    in.data=in.data[complete.cases(in.data), ]

    ss.in.data = subset(in.data, Grp=="MDD")
    print(t.test(ss.in.data[ss.in.data$timepoint=="A", y.axis], ss.in.data[ss.in.data$timepoint=="C", y.axis], paired=TRUE))

    ## setup the LME models and ANOVAs derived therefrom
    model.lme.interaction=lme(data=in.data, fixed=model.formula.interaction, random=random.formula, correlation=corAR1(form=random.formula))
    ## print(model.lme.interaction)
    model.lme.nointeraction=lme(data=in.data, fixed=model.formula.nointeraction, random=random.formula, correlation=corAR1(form=random.formula))
    ## print(model.lme.nointeraction)
    model.anova.interaction=anova(model.lme.interaction)
    model.anova.nointeraction=anova(model.lme.nointeraction)

    cat("\n### Interaction LME summary\n")
    print(summary(model.lme.interaction))
    cat("\n### Nointeraction LME summary\n")
    print(summary(model.lme.nointeraction))    

    cat("\n### Interaction ANOVA summary\n")
    print(model.anova.interaction)
    cat("\n### Nointeraction ANOVA summary\n")
    print(model.anova.nointeraction)

    use.interaction.model=TRUE
    ## if (model.anova.interaction[4, "p-value"] < 0.05) {
    ##     use.interaction.model=TRUE
    ##     cat("\n")
    ##     cat("**************************************************************\n")
    ##     cat("*** The Grp * timepoint is significant: using interation model\n")
    ##     cat("**************************************************************\n")
    ## }

    
    ## ## original from
    ## ## http://cran.r-project.org/web/packages/multcomp/vignettes/multcomp-examples.pdf
    ## ## mod <- lm(breaks ~ wool * tension, data = warpbreaks)
    ## ## tmp <- expand.grid(tension = unique(warpbreaks$tension), wool = unique(warpbreaks$wool))
    ## ## X <- model.matrix(~ wool * tension, data = tmp)
    ## ## glht(mod, linfct = X)
    ## ## Tukey <- contrMat(table(warpbreaks$tension), "Tukey")
    ## ## K1 <- cbind(Tukey, matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)))
    ## ## rownames(K1) <- paste(levels(warpbreaks$wool)[1], rownames(K1), sep = ":")
    ## ## K2 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)), Tukey)
    ## ## rownames(K2) <- paste(levels(warpbreaks$wool)[2], rownames(K2), sep = ":")
    ## ## K <- rbind(K1, K2)
    ## ## colnames(K) <- c(colnames(Tukey), colnames(Tukey))
    ## ## ## and perform the tests via
    ## ## summary(glht(mod, linfct = K %*% X))
    
    if (use.interaction.model==FALSE) {
        cat("### Main effects in nointeraction model\n")
        assign("model.formula.nointeraction", model.formula.nointeraction, envir=.GlobalEnv)
        K1 <- glht(model.lme.nointeraction, mcp(Grp = "Tukey"))$linfct
        K2 <- glht(model.lme.nointeraction, mcp(timepoint = "Tukey"))$linfct
        print(summary(glht(model.lme.nointeraction, linfct = rbind(K1, K2))))
        remove("model.formula.nointeraction", envir=.GlobalEnv)        
    }

    if (use.interaction.model==TRUE) {
        cat("### Interaction effects\n")
        tmp <- expand.grid(timepoint = unique(in.data$timepoint), Grp = unique(in.data$Grp))
        X.a <- model.matrix(~ Grp * timepoint, data = tmp)
        ## print(glht(model.lme, linfct = X.a))
        Tukey.a <- contrMat(table(in.data$timepoint), "Tukey")
        Ka1 <- cbind(Tukey.a, matrix(0, nrow = nrow(Tukey.a), ncol = ncol(Tukey.a)))
        rownames(Ka1) <- paste(levels(in.data$Grp)[1], rownames(Ka1), sep = ":")
        Ka2 <- cbind(matrix(0, nrow = nrow(Tukey.a), ncol = ncol(Tukey.a)), Tukey.a)
        rownames(Ka2) <- paste(levels(in.data$Grp)[2], rownames(Ka2), sep = ":")
        K.a <- rbind(Ka1, Ka2)
        colnames(K.a) <- c(colnames(Tukey.a), colnames(Tukey.a))

        print(summary(glht(model.lme.interaction, linfct = K.a %*% X.a)))


        ##tmp <- expand.grid(Grp = unique(in.data$Grp), timepoint = unique(in.data$timepoint))
        X.b <- model.matrix(~ timepoint * Grp, data = tmp)
        ## print(glht(model.lme, linfct = X.b))
        Tukey.b <- contrMat(table(in.data$Grp), "Tukey")
        Kb1 <- cbind(Tukey.b, matrix(0, nrow = nrow(Tukey.b), ncol = ncol(Tukey.b)))
        rownames(Kb1) <- paste(levels(in.data$timepoint)[1], rownames(Kb1), sep = ":")
        Kb2 <- cbind(matrix(0, nrow = nrow(Tukey.b), ncol = ncol(Tukey.b)), Tukey.b)
        rownames(Kb2) <- paste(levels(in.data$timepoint)[2], rownames(Kb2), sep = ":")
        K.b <- rbind(Kb1, Kb2)
        colnames(K.b) <- c(colnames(Tukey.b), colnames(Tukey.b))
        ## and perform the tests via
        print(summary(glht(model.lme.interaction, linfct = K.b %*% X.b)))

        assign("in.data", in.data, envir=.GlobalEnv)
        assign("model.formula.interaction", model.formula.interaction, envir=.GlobalEnv)
        
        cat("*** Tables of means\n")
        print(allEffects(model.lme.interaction, formula=model.formula.interaction))

        remove("model.formula.interaction", envir=.GlobalEnv)
        remove("in.data", envir=.GlobalEnv)
    }

    return(smry)
}

## Reads the seed file and does the (crude) equivalent of BASH variable
## substitution
readSeedsFile <- function (inSeedsFile) {
    cat("*** Reading seed from", inSeedsFile, "\n")
    table=scan(inSeedsFile, what=character(), quiet=TRUE)
    table=gsub("$DATA", seeds.data.dir, table, fixed=TRUE)

    return (table)
}

## extracts the seed name from a file path name pointing to a NIfTI
## file containing the seed
getSeedName <- function(inSeedPath){
    name=basename(inSeedPath)
    if (grepl("\\.nii", name)) {
        return(gsub("\\.nii.*", "", name))
    } else if (grepl("\\+tlrc", name)) {
        return(gsub("\\+tlrc.*", "", name))
    } else {
        return (name)
    }
}

##########################################################################################################################################################################
### END OF FUNCTIONS #####################################################################################################################################################
##########################################################################################################################################################################

if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
    ncpus=as.integer(strsplit(system("sysctl hw.ncpu", intern=T), ' ')[[1]][2])
    cat(paste("Found" , ncpus, ifelse(ncpus == 1, "cpu", "cpus"), "\n"))
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
    ncpus=8
    cat(paste("Found" , ncpus, ifelse(ncpus == 1, "cpu", "cpus"), "\n"))    
} else {
    stop(paste("Sorry can't set data directories for this computer\n"))
}

scripts.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/scripts")
data.dir=file.path(root.dir, "sanDiego/rsfcGraphAnalysis/data/")
admin.data.dir=file.path(data.dir, "admin")
config.data.dir=file.path(data.dir, "config")
## group.data.dir=file.path(data.dir, "Group.data.withMedicated")
## group.results.dir=file.path(data.dir, "Group.results.withMedicated")
group.data.dir=file.path(data.dir, "Group.data")
group.results.dir=file.path(data.dir, "Group.results")
seeds.data.dir=file.path(data.dir, "seeds")

mgd.filename=file.path(admin.data.dir, "merged.demographics.and.neuropsych.csv")
mgd=readCsvFile(mgd.filename)

## only keep timepoints A and C
mgd=mgd[ mgd$timepoint %in% c("A", "C"), ]
mgd$timepoint=droplevels(mgd$timepoint)

mgd=mgd[ mgd$Grp %in% c("NCL", "MDD"), ]
mgd$Grp=droplevels(mgd$Grp)
mgd$DOB=as.Date(mgd$DOB)
mgd$MRI=as.Date(mgd$MRI)

mgd=computeAge(mgd)


subjectList=data.frame("subject"=dir(path="/Volumes/data/sanDiego/machLearnT1Analysis/data/", pattern="*_A*$"))
subjectList=splitSubjectOrderIntoIdAndTimepoint(subjectList)
subjectList=fixSubjectIds(subjectList, "id")

## we're missing the CGAS score for 136_C but for some reason it's
## entered as 0 rather than NA, so fix that
mgd[mgd$ID==136, "CGAS"]=NA

seeds.filename=file.path(config.data.dir, "juelich_whole_amygdala_seeds.txt")
seeds=readSeedsFile(seeds.filename)

## for (variable in c("CGAS", "CDRS.t.score", "MASC.tscore", "RADS.Total.tscore") ) {
## for (variable in c("CGAS") ) {
## for (variable in c("BDI.II.Total") ) {
## for (variable in c("MASC.tscore") ) {    


for (variable in c("CDRS.t.score") ) {

    cat("##############################################################################################################################\n")
    cat("*** Computing usable subjects for", variable, "\n")

    ss=subset(mgd, select=c("ID", "timepoint", "Grp", "age.in.years", variable))
    subjects.lists=find.usable.subjects(ss, variable)

    sl=subjects.lists[["usable"]]

    sl=compute.residualized.scores(compute.change.scores(sl, variable), variable)

    finite.values=is.finite(sl[ , paste(variable, "scaled.diff", sep=".") ])
    count.finite.values=sum(! finite.values)
    cat("*** Removing", count.finite.values, "infinite value(s)\n")
    sl=sl[ finite.values, ]
    
    cat("*** Number of usable subjects for", variable, "is:\n")
    print(addmargins(table(sl$Grp)))
    cat("*** The following subjects are usable for", variable, ":\n")
    print.usable.subjects(sl)
    sl = scale.variable(sl, c("age.in.years", paste(variable, c("A", "C"), sep=".")))

    sl$ExclusionFile=sapply(sl$ID, function (id) {
        ## rsfcPreprocessed/00_DO_NOT_ANALYSE_105_A_20percent.txt
        file.path(data.dir, paste(id, "_A", sep=""), "rsfcPreprocessed", paste("00_DO_NOT_ANALYSE", id, "A_20percent.txt", sep="_"))
    })
    sl$ExclusionFileExists=lapply(sl$ExclusionFile, function (ff) {file.exists(ff)})
    ## print(sl)

    ## exclude those subjects with too much motion
    sl=subset(sl, ExclusionFileExists==FALSE)
    
    for (seed in seeds) {
        seedName=getSeedName(seed)
 
        ## old AFNI+FSL analysis
        sl$InputFile=sapply(sl$ID, function (id) {
            file.path(data.dir, paste(id, "_A", sep=""), "rsfc", seedName, paste(seedName, "z-score+tlrc.HEAD", sep="."))
        })
       
        sl$InputFileExists=lapply(sl$InputFile, function (ff) {file.exists(ff)})
        sl=subset(sl, InputFileExists==TRUE)

        data.table=sl[, c("ID", "Grp", "age.in.years", paste(variable, "scaled.diff", sep="."), "InputFile")]
        names(data.table)[1]="Subj"

        save.data.table(data.table, variable, seedName, in.timepoint="A", in.file.type="tab")
        ## print(data.table)
    }
    
}
