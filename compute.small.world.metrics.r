if (! exists("g.attributes") ) {
    ## start with a clean slate if we've not already loaded the graph
    ## lists
    rm(list=ls())
    graphics.off()
}

devtools::load_all(file.path(Sys.getenv("HOME"), "src", "brainGraph"))

load.saved.data.structures <- function () {
    cat("*** Loading pregenerated graph and associated data structres from", saved.graph.data.structures.filename, "\n")
    pigz.load(saved.graph.data.structures.filename, envir = parent.frame(), verbose = TRUE)
    
    cat("*** Loaded saved graph structures\n")
    
    cat("*** The subject distribution is as follows:\n")
    subject.distribution=addmargins(xtabs(~ Group + Gender, data=characteristics.df))
    print(subject.distribution)
    
    cat("*** Loaded graphs for the following", length(graph.densities),
        ifelse(length(graph.densities) > 1, "densities:", "density:"),
        paste(graph.densities, collapse=" "), "\n")
    
    ## a bunch of sanity check to ensure the dimensions of various
    ## data strcutures match
    if (dim(characteristics.df)[1] != length(g)) {
        cat("*** The number of subjects (",
            length(dim(characteristics.df)[1]),
            ") and the number of subjects (",
            length(g), ") does not match\n", sep="")
    }

    if (length(graph.densities) != length(g[[1]])) {
        cat("*** The number of graph densities (",
            length(graph.densities),
            ") and the number of lists for those densities (",
            length(g[[1]]),
            ") does not match\n", sep="")
    }

    if (length(thresholded.matrices) != dim(characteristics.df)[1]) {
        cat("*** The number of thresholded matrices (",
            length(thresholded.matrices),
            ") and the number of subjects (",
            dim(characteristics.df)[1],
            ") does not match\n")
    }

    ## the data frames and arrays
    aa=cbind(as.character(characteristics.df$Study.ID), names(netcc[1, 1, ]), as.character(netcc.filenames.and.subjects.df$Study.ID))
    if (! all(apply(aa, 1, function (xx) { ! length(unique(xx)) > 1 } ))) {
        stop("*** The subject IDs in charactersitics.df, netcc[1, 1, ], and netcc.filenames.and.subjects.df$Study.ID do not match.\nCannot continue\n")
    } else {
        cat("*** Subject IDs in charactersitics.df, netcc, and netcc.filenames.and.subjects.df all match\n")
    }
    ## delete aa, it's no longer needed
    rm(aa)
}

simulate.random.networks <- function (in.g, N, compress.rds=TRUE, start.at.subject=1, end.at.subject=length(in.g), ...) {
    if ( ! all(unlist(
        lapply(in.g,
               function (aa) {
                   lapply(aa,
                          function (bb) {
                              all(c("Group", "name") %in% names(graph_attr(bb)))
                          })
               }))) ) {
        stop("*** Graphs do not have Group and name attributes. Cannot continue.\n")
    }

    if (is.null(names(in.g))) {
        stop("*** Elements of in.g must have names corresponding to the subject IDs\n")
    }

    if (start.at.subject < 1) {
        stop(paste("*** Subject enumeration starts at 1 not", start.at.subject, "\n"))
    }

    if (start.at.subject > length(in.g)) {
        stop(paste("*** Subject enumeration cannot start at", start.at.subject, "when there are only", length(in.g), "subjects\n"))
    }
 

    if (end.at.subject > length(in.g)) {
        stop(paste("*** Subject enumeration ends at", length(in.g), "not", end.at.subject, "\n"))
    }

    subject.sequence=seq.int(start.at.subject, end.at.subject)
    subject.count=length(subject.sequence)
    cat("*** Enumerating over", subject.count, "subjects\n")
    
    suppressMessages(library(chron))
 
    begin.at=Sys.time()
    cat(sprintf("*** Starting simulations of random networks at %s\n", begin.at))
   
    phi.norm <- foreach (ss = iter(subject.sequence)) %do% {
        start=Sys.time()
        cat(sprintf("*** Started subject %s (%03d of %03d) at %s. ",
                    names(in.g)[ss], ss-(start.at.subject-1), subject.count, start))
        ## set the subject name
        subject.name=names(in.g)[ss]
        
        foreach (dd = icount(length(in.g[[ss]]))) %dopar% {
            rand <- sim.rand.graph(in.g[[ss]][[dd]], N, use.parallel=FALSE, ...)
            saveRDS(rand, compress=compress.rds,
                    file=file.path(group.results.dir, '/',
                        sprintf('rand_%s_thr_%02i%s', subject.name, dd, '.rds')))
            ## return the rich club statistics data.table to be
            ## accumulated in phi.norm by foreach
            rich.club.norm(in.g[[ss]][[dd]], rand=rand)
            ## rm(rand)
            ## gc()
        }
        end=Sys.time()
        
        time.taken=format(as.chron(end) - as.chron(start))
        cat (sprintf("Ended at %s. Time taken %s. %0.2f%% completed\n", end, time.taken, ((ss-(start.at.subject-1))/subject.count)*100))
    }

    end.at=Sys.time()
    total.time.taken=format(as.chron(end.at) - as.chron(begin.at))
    cat (sprintf("*** Ended simulations of random networks at %s. Time taken %s\n", end, total.time.taken))

    return(phi.norm)
}    

compute.small.world.metrics <- function (in.g, N, compress.rds=TRUE) {

    if (is.null(names(in.g))) {
        stop("*** Elements of in.g must have names corresponding to the subject IDs\n")
    }

    ## Get all small-worldness and random measures, all thresholds
    ## ###########################################################
    ## groups <- covars[, levels(Group)]
    subjects = names(in.g)
    ## rand.dt <- small.dt <- vector('list', length=length(in.g))

    create.rand.dt <- function(rand, Study.ID, group, V1, kNumRand) {
        rand.mod <- sapply(rand, sapply, graph_attr, 'mod')
        rand.E.global <- sapply(rand, sapply, graph_attr, 'E.global')
        rand.Cp <- sapply(rand, sapply, graph_attr, 'Cp')
        rand.Lp <- sapply(rand, sapply, graph_attr, 'Lp')
        rand.dt <- data.table(V1=rep(V1, each=kNumRand),
                              Study.ID=rep(Study.ID, kNumRand),
                              Group=rep(group, kNumRand), mod=as.vector(rand.mod),
                              Cp=as.vector(rand.Cp), Lp=as.vector(rand.Lp),
                              E.global=as.vector(rand.E.global))
        return(rand.dt)
    }

    ## from
    ## http://stackoverflow.com/questions/19791609/saving-multiple-outputs-of-foreach-dopar-loop
    comb <- function(x, ...) {
        lapply(seq_along(x),
               function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
    }
    
    start=Sys.time()
    cat(sprintf("*** Starting to compute small world network metrics for %03d subjects at %s. ",
                length(in.g),  start))

    ##ss is the subject iterator
    ret.vals = foreach (ss = icount(length(in.g)), .combine='comb', .multicombine=TRUE,
        .init=list(list(), list())) %dopar% {
        subject.name=names(in.g)[ss]
        subjects.group=as.character(characteristics.df[characteristics.df$Study.ID == subject.name, "Group"])

        ## get the densities for each subject from thier graphs not
        ## from names on the list of densities (which are no longer
        ## ste by make.graphs.r), doing the latter can cause problems
        ## when thresholds are not rounded exactly to match the list
        ## of graph.densities originally used in make.graphs.r
        densities=vapply(in.g[[ss]], function(yy) { yy$density }, numeric(1), USE.NAMES=FALSE)
        
        ## the slicing at then end of the next statement ensures that
        ## we only take as many files as there are densities. We do
        ## this because the glob will match all files irrespective of
        ## how many densities per subject are actually in the in.g
        ## list. The situation may arrise where the length of the
        ## filenames list is longer than the number of densities per
        ## subject as might occur when testing the efficacy of the
        ## function with fewer densities than there filenames per
        ## subject. When this happens the small.dt array is populated
        ## by repeatedly recycling the in.g[[s]] elememnts with rand_all
        ## until there are as many elements in small.dt as there
        ## filenames (the length of rand_all being determined by the
        ## length of fnames)
        
        fnames <- list.files(group.results.dir,
                             sprintf('rand_%s_thr.*', subject.name), full.names=T)[1:length(densities)]
        rand_all <- lapply(fnames, readRDS)
        
        saveRDS(rand_all, compress=compress.rds,
                file=file.path(group.results.dir, paste('rand', subject.name, 'all.rds', sep="_")))
        small.dt <- small.world(in.g[[ss]], rand_all)
        small.dt[, Study.ID := rep(subject.name,   length(densities))]
        small.dt[, Group    := rep(subjects.group, length(densities))]
        rand.dt <- create.rand.dt(rand_all, subject.name, subjects.group, V1=densities, N)        
        setnames(rand.dt, "V1", "density")
        setkey(rand.dt, density)
        rm(rand_all)
        ## let's leave it to R to decide when to garbage collect
        ## rather than forcing it to happen
        ## gc()

        ## return the small.dt and rand.dt data tables ot be
        ## accumulated into ret.vals by foreach
        list(small.dt, rand.dt)
    }

    end=Sys.time()
    suppressMessages(library(chron))
    time.taken=format(as.chron(end) - as.chron(start))
    cat (sprintf("Ended at %s. Time taken %s.\n", end, time.taken))
    
    ## small.dt is a data table containing the small world properties
    ## for each subject
    ##
    ## rand.dt is a data table containing the mod, Cp, Lp, E.global
    ## metrics for each of the random networks (the number of which is
    ## dicted by N) used to computed the normalized Cp and Lp metrics
    ## stored in small.dt
    rand.dt <- rbindlist(ret.vals[[2]])
    small.dt <- rbindlist(ret.vals[[1]])
    setkey(small.dt, Study.ID, Group, density)

    setcolorder(small.dt, c("Study.ID", "Group", setdiff(colnames(small.dt), c("Study.ID", "Group"))))
    setcolorder(rand.dt,  c("Study.ID", "Group", setdiff(colnames(rand.dt),  c("Study.ID", "Group"))))    

    return(list("rand.dt"=rand.dt, "small.dt"=small.dt))

}

####################################################################################################
## END OF FUNCTIONS
####################################################################################################

### setup path variables
if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
    max.cpus=as.integer(strsplit(system("sysctl hw.ncpu", intern=T), ' ')[[1]][2])
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
    max.cpus=as.integer(system("getconf _NPROCESSORS_ONLN", intern=TRUE))
} else {
    cat(paste("Sorry can't set data directories for this computer\n"))
    max.cpus=NA
}

source("common.functions.r")

library(plyr)
library(data.table)
library(foreach)
library(iterators)

study.root.dir=file.path(root.dir, "sanDiego/machLearnT1Analysis")

scripts.dir=file.path(study.root.dir, "scripts")
data.dir=file.path(study.root.dir, "data")
group.results.dir=file.path(data.dir, "Group.results")


## use the partial correlation (PC) coefficients created by 3dNetCorr
matrix.type="PC"

saved.graph.data.structures.filename=file.path(group.results.dir, paste("graphs", matrix.type, "Rdata", sep="."))
if ( ! exists("g.attributes") ) {
    cat("*** g.attributes doesn't exist. Trying to load from saved data file\n")
    if (file.exists(saved.graph.data.structures.filename)) {
        load.saved.data.structures()
    } else {
        stop(paste("*** Cannot continue. No such file", saved.graph.data.structures.filename, "\n"))
    }
} else {
    cat("*** Found a pre-existing g.attributes in the environment. Skipping loading from file\n")
    cat("*** Execute rm(list=ls() to force reloading\n")
}

if (is.null(names(g.attributes))) {
    cat("*** Setting names on g.attributes\n")
    names(g.attributes) = names(g)
    g.attributes = lapply(g.attributes, function(xx) { names(xx) = names(g[[1]]) ; return(xx) })
}

parallel.executation=TRUE
if (parallel.executation) {
    cat("*** Enabling parallel processing\n")
    library(doMC)
    ## why 51? because there are 51 densities run per subject. It's
    ## only choosen to get as many densities processed in parallel as
    ## possible
    registerDoMC(cores = min(51, max.cpus))
    cat("*** Using", getDoParWorkers(), "CPU cores\n")
    cat("*** Plyr progress bars are disabled when parallel computation is enabled\n")
    progress.bar.type='none'        
} else {
    cat("*** Plyr progress bars are set to text\n")        
    progress.bar.type='text'
}

## pick only the first 4 densities for the first 4 subjects to have
## something quick to work with for testing purposes
## aa=lapply(g.attributes[1:3], function (xx) xx[1:4])

## pick only the first 4 subjects and all of their densities
## aa=g.attributes[1:4]

## small.world.properties(aa, 100)
## rand <- sim.rand.graph.par(aa[[1]][[1]], 100)
## tt=small.world.properties(aa, 100, clustering=FALSE, compress.rds=FALSE)
phi.norm=simulate.random.networks(g.attributes, 100, clustering=FALSE, compress.rds=FALSE, start.at.subject=1)
small.world.metrics=compute.small.world.metrics(g.attributes, 100, compress.rds=FALSE)
## add in the rich club statistics for compatability with
## analysis_random_graphs from brainGraph
if (exists("phi.norm"))
    small.world.metrics$rich = phi.norm

cat("*** Table of density threshold per subject\n")
cat("*** Columns should add to", length(g.attributes), "\n")
cat("*** Rows should add to", length(graph.densities[[1]]), "\n")
print(addmargins(xtabs(~ Study.ID + density, data=small.world.metrics$small.dt)))

small.world.metrics.filename=file.path(group.results.dir, paste("small.world.metrics", "Rdata", sep="."))
cat("*** Saving small.world.metrics to", small.world.metrics.filename, "\n")
pigz.save(small.world.metrics, file=small.world.metrics.filename)

