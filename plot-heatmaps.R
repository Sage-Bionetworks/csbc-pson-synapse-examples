suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("plyr"))

## This is a simple script intended to demonstrate:
## (1) Retrieving data from Synapse
## (2) Storing a result generated from that data back to Synapse
## (3) Setting the "provenance" of those results, which indicates
##     how the file was created (e.g., that this script was used to do so,
##     the inputs required of the script, and who executed the script)

## The script simply plots a heatmap for each expression matrix in the
## user-specified directory--which is passed as a SynapseID.
## i.e., the directory is a folder in Synapse and has expression matrices.

## e.g., to create heatmaps of the scRNA-seq glioblastoma data
## generated at Columbia stored at syn11678365 and to store those
## heatmaps in the dummy folder syn13364025
## Rscript ./plot-heatmaps.R --input-folder-synapse-id=syn11678365 --output-folder-synapse-id=syn13364025


## The Synapse R library.
## I needed to be root to install this, which I did from an R session as:
## install.packages("synapser", repos=c("https://sage-bionetworks.github.io/ran", "http://cran.fhcrc.org"))
## More information is here:
## http://docs.synapse.org/articles/getting_started.html#installing-synapse-clients
suppressPackageStartupMessages(library(synapser))

option_list <- list(
    make_option(c("--input-folder-synapse-id"), action="store",
                default=NULL,
                help="The Synapse ID of the folder holding expression matrices."),
    make_option(c("--output-folder-synapse-id"), action="store",
                default=NULL,
                help="The Synapse ID of the folder in which to store the heatmaps.")
)

descr <- "\
Create heatmaps of any expression matrices stored in the specified Synapse folder.  Output the heatmaps to output Synapse folder, setting their provenance to indicate that this script created them and the input parameters based to this script.

This script is intended as a demo to accessing/storing files in Synapse and to using provenance.
"

parser <- OptionParser(usage = "%prog [options]", option_list=option_list, description=descr)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if ( length(arguments$args) != 0 ) {
  print_help(parser)
  q(status=1)
}

input.folder.synapse.id <- opt$`input-folder-synapse-id`
output.folder.synapse.id <- opt$`output-folder-synapse-id`

if(is.null(input.folder.synapse.id) || is.null(output.folder.synapse.id)) {
    print_help(parser)
    q(status=1)
}

## Log in to Synapse.  You don't need to supply your password if you have
## set up a credential file.  See how to do so here:
## http://docs.synapse.org/articles/getting_started.html#using-a-config-file
synLogin()

## Iterate over the contents of the folder.
children <- synGetChildren(input.folder.synapse.id)
children <- as.list(children)

l_ply(children,
      .fun = function(entity) {
          ## Retrieve the annotations associated with this object
          annotations <- synGetAnnotations(entity$id)

          ## Only process this file if it is annotated as a gene expression
          ## file.  'dataType' is a required keyword.
          if(is.null(annotations$dataType) ||
             (annotations$dataType != "geneExpression")) {
              return()
          }

          ## Also require that it be a dataMatrix
          if(is.null(annotations$dataSubtype) ||
             (annotations$dataSubtype != "dataMatrix")) {
              return()
          }

          ## This is a gene expression matrix.  Let's retrieve it.
          obj <- synGet(entity$id)
          name <- obj$properties$name

          mat <- read.table(obj$path, sep=",", header=TRUE)

          ## Create a heatmap of a subset of the rows, just for
          ## efficiency sake.
          output.name <- gsub(name, pattern="^(.*)\\.[^\\.]+$", replacement="\\1.pdf")
          pdf(output.name)
          heatmap(as.matrix(mat[1:min(nrow(mat), 100), ]))
          d <- dev.off()

          ## Create a Synapse File object, linking it to the output folder
          output.file <- File(path = output.name, parent = output.folder.synapse.id)

          ## Provenance (i.e., how results were generated from data) will
          ## be described by an Activity.
          ## Indicate the script that was executed, i.e., the github
          ## url of this script.
          ## Also specify the input to this script in 'used'.
          act <- Activity(
              name = "clustering",
              description = "Heatmap with dendrogram clustering",
              used = c(input.folder.synapse.id, output.folder.synapse.id),
              executed = "https://github.com/Sage-Bionetworks/csbc-pson-synapse-examples/plot-heatmaps.R")
          
          ## Store the File in Synapse
          output.file <- synStore(output.file, activity = act)
      })

