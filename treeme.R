#!/usr/bin/env Rscript
options(warn=-1)
#packages = c("optparse", "crayon", "ggtree", "ggplot2", "treeio", "glue", "cowplot")

# if taxa labels are bunching up, use a bigger output

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(crayon))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(glue))
library(grid)

usage = "Create pretty tree using ggtree.
  Usage:
    treeme.R  -t [TREEFILE] \
              -o [OUTPUT.PDF] \
              -c [CLADES.TXT] \
              -m [META.TSV] \
              -l [VARIABLE] \
              -p [VARIABLE] \
              -g [TITLE]

  -t,  --tree         path and name of nhx tree to ingest.
  -o,  --output       output name and path of graphics output (will be a pdf).
  -c   --clades       tsv file containing clade_name\\tnode_number.
  -m   --metaFile     tsv meta-file containing meta data.
                      First column must be taxa names and must match exactly to tree.

  -l   --colorTaxa    Variable name (in meta file) to color the taxa labels.
  -p   --tipPoint     Variable name (in meta file) to color and shape the tip points.
  -s   --paperSize    Size of output pdf. A3p, A3l, A4l, [A4p]. p = portrait. l = landscape
  -g   --title        Title of plot - text will be added to to left of tree

  Output:
    Phylogenetic tree visualised using ggtree (pdf file).
  "

option_list = list(
  make_option(c("-t","--tree"),
              help = "input nhx tree",
              action = "store",
              type = "character", default = NA),
  make_option(c("-o","--output"),
              help = "name of output graphic, will be pdf",
              action = "store",
              type = "character", default = NA),
  make_option(c("-c","--cladesFile"),
              help = "clade(str) and node number in tsv format",
              action = "store",
              type = "character", default = NA),
  make_option(c("-m","--metaFile"),
              help = "tsv file containing meta data",
              action = "store",
              type = "character", default = NA),
  make_option(c("-l", "--colorTaxa"),
              help = "which variable (in meta file) to color taxa labels",
              action = "store",
              type = "character", default = NA),
  make_option(c("-p","--tipPoint"),
              help = "which variable (in meta file) to color/shape tip points",
              action = "store",
              type = "character", default = NA),
  make_option(c("-g","--title"),
              help = "The title of the final output",
              action = "store",
              type = "character", default = NA),
  make_option(c("-s", "--paperSize"),
              help = "output size - options: A3p, A3l, A4l, [A4p] - p/l is the orientation",
              action = "store",
              type = "character", default = 'A4p')
)

parser = OptionParser(
  usage = paste("%prog -t [TREEFILE] \
                -o [OUTPUT.PDF] \
                -c [CLADES.TXT] \
                -m [META.TSV] \
                -l [VARIABLE] \
                -p [VARIABLE]
                -s [A4l]
                -g [TITLE]",
                "Generate grahical tree with ggtree",
                sep="\n"),
  epilogue = "All options are required (maybe)",
  option_list = option_list)
##########################################
########## Parameter Validation ##########
##########################################

# custom function to stop quietly
stop_quietly = function(message, args) {

  opt = options(show.error.messages = FALSE)
  on.exit(options(opt))
  cat(message, sep = "\n")
  message(missing_message(args))
  quit()
}

# custom function to print missing arguments
missing_message = function(arguments) {
  missing_options = names(which(is.na(arguments)))
  m = glue_collapse(missing_options, sep = ", ")
  out = glue_col("{red Missing arguments:} {yellow {m}}")

  return(out)
}

tryCatch(
  expr = {

    arguments = parse_args(object = parser, positional_arguments = TRUE)$options
  },
  error = {
    function(e){
      print(e)
      stop(missing_message(arguments))}
  },
  finally = {
    if(any(is.na(arguments))) {
      opts_missing = which(is.na(arguments))
      stop_quietly(message = usage, args = arguments)
    }
  }
)

#########################################
############# Input checks ##############
#########################################

#check tree file
tryCatch(
  expr = {
    tree = read.beast(arguments$tree)
    message(crayon::green("Success: Tree file was read."))
  },
  error = function(e){
    print(e)
    stop(crayon::bgCyan("Reading tree file. Is it nhx?"))
  }
)

# read in meta file
tryCatch(
  expr = {
    meta = read.delim(arguments$metaFile, sep = '\t', na.strings = '', stringsAsFactors = F)
    cladesFile = read.delim(arguments$cladesFile, sep = "\t", na.strings = '', stringsAsFactors = F)
    message(crayon::green("Success: Meta file read."))
  },
  error = function(e){
    print(e)
    stop(crayon::bgCyan("Unable to read in meta or clade file.
                     Do they exist and is the path correct?"))
  }
)
#
# meta file variables match cli input
if (!(arguments$colorTaxa %in% names(meta))) {
  stop(crayon::bgCyan("CHECK YOUR META FILE. Check if this header is in the meta file:", arguments$colorTaxa))
}

if (!(arguments$tipPoint %in% names(meta))) {
  stop(crayon::bgCyan("CHECK YOUR META FILE. Check if this header is in the meta file:", arguments$tipPoint))
}

###############################################
########### Functions for Tree ################
###############################################

tree_title = function(title) {
  today = format(Sys.time(),  format = "%d %b %Y")
  out = paste0(as.character(title), '\n', today)
}

calc_text_size = function(height, lines, buffer) {
  max_line_height = (height)/lines #mm
  return(max_line_height)
}

check_taxa_names = function(tree_labs, meta_desig) {
  # tree_labs  :  tree@phylo$tip.label
  # meta_desig :  meta[, 'designation']
  mismatch <- tree_labs[!tree_labs %in% meta_desig]
  if (!length(mismatch)) {
    message(crayon::green("Good: All tree tip labels present in meta file"))
  } else {
    message(
      crayon::yellow("\n"),
      crayon::yellow("Warning - the following tip labels not present in meta file:\n"),
      crayon::yellow(paste("\t", mismatch, collapse = "\n")),
      crayon::yellow("\n")
    )
  }
}

# clade generator
add_clades = function(cladesFile, tree_data){
  # clade file must be a dataframe with vars:
  #   node_number
  #   clade_name
  # tree_data  :  tree@data from read.beast

  create_clade = function(node_number, clade_name, offset) {
    geom_cladelabel(node_number,
                    label = clade_name,
                    align = FALSE,
                    angle = 270,
                    hjust = 'center',
                    offset = offset,
                    offset.text = plot_dim_x * 0.01,
                    barsize = 1,
                    fontsize = 5,
                    extend = 0.2)
  }

  get_clade_node = function(muts, tree_data) {
    # returns node number given a mutation

    ind = which(tree_data$aa_muts %in% muts)
    if(length(ind) == 0) {
      return(NA_integer_)
    }
    node = tree_data[[ind, "node"]]
    return(node)
  }

  cladesFile$node_number = apply(cladesFile %>% dplyr::select(mutations), 1, get_clade_node, tree_data)
  cladesFile = cladesFile[complete.cases(cladesFile),]

  num_clade = nrow(cladesFile)
  min_offset = plot_dim_x * 0.02
  from = plot_dim_x - (plot_dim_x * 0.80)
  to = plot_dim_x + (plot_dim_x * 0.30)
  offset_values = seq(from, to, by = min_offset)[1:3]
  offset = rep(offset_values, ceiling(num_clade/3))[1:num_clade]

  return(purrr::pmap(list(
    as.integer(cladesFile$node_number),
    cladesFile$clade_name,
    offset),
    create_clade)
  )
}

get_paper_size = function(size_argument) {
  paperSizes = list(
    "A2p" = c(297*2, 420*2),
    "A2l" = c(420*2, 297*2),
    "A3p" = c(297, 420),
    "A3l" = c(420, 297),
    "A4p" = c(210, 297),
    "A4l" = c(297, 210)
  )
  if (!any(size_argument %in% names(paperSizes)) ){
    stop(crayon::red("Paper size unknown. Options:  A3p A3l A4p A4l"))
  }
  return(paperSizes[[size_argument]])
}

cat('Input arguments', '\n',
     '\tTree: ', arguments$tree, '\n',
     '\tOutput: ',{arguments$output}, '\n',
     '\tClade File: ',{arguments$cladesFile},'\n',
     '\tMeta File: ',{arguments$metaFile},'\n',
     '\tColor labels by: ',{arguments$colorTaxa},'\n',
     '\tColor node tips by: ',{arguments$tipPoint},'\n',
     '\tTitle: ', {arguments$title},'\n',
     '\tOutput Size: ', {arguments$paperSize}, '\n'
     )
message(crayon::green("~~~All Checks Okay - Plotting tree~~~"))

##############################################
################ Tree plotting ###############
##############################################

cols_file = read.delim(file = "scripts/colors.tsv", sep = "\t")
shapes_file = read.delim(file = "scripts/shapes.tsv", sep = "\t")

categories = cols_file$categories
color_pal = cols_file$color_pal
shapes_type = shapes_file$shapes_type
shape_colors = shapes_file$shape_colors
shape_cats = shapes_file$shape_cats
names(shapes_type) = shape_cats

check_taxa_names(tree@phylo$tip.label, meta[, "strain"])
output_size = get_paper_size(arguments$paperSize)

colorTaxa = arguments$colorTaxa
tipPoint = arguments$tipPoint
tipLabSize = calc_text_size(height = output_size[[1]],
                            lines = tree@phylo$Nnode,
                            buffer = 15)

treeplot = ggtree(tree)  %<+% meta +
  geom_tiplab(geom = "text",
              size = tipLabSize,
              aes(color = !! sym(colorTaxa)),
              key_glyph = cowplot::rectangle_key_glyph(fill = color,
                                                       padding = margin(0, 0, 0, 0),
                                                       color = 'black',
                                                       linetype = 3),
              offset = 0.00009,
              family = "Arial") +

  geom_tippoint(size = tipLabSize, aes(fill = !! sym(tipPoint),
                                shape = !! sym(tipPoint))) +

  scale_color_manual(colorTaxa,
                     limits = categories,
                     values = color_pal,
                     na.value = "#000000") +

  scale_fill_manual(tipPoint,
                    values = shape_colors,
                    limits =  shape_cats,
                    na.value = "#000000") +

  scale_shape_manual("Legend",
                     values = shapes_type,
                     breaks = shapes_type) +

  guides(fill = guide_legend(override.aes = list(size = tipLabSize * 1.5,
                                                 label = "",
                                                 shape = shapes_type))) +
  guides(color = guide_legend(override.aes = list(size = tipLabSize * 5,
                                                 label = "\u25A0",
                                                 linetype = 3))) +

  ggtitle(tree_title(arguments$title)) +

  theme(legend.position = c(0.1, 0.65),
        legend.key.size = unit(tipLabSize*2, "mm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"),
        legend.spacing.y = unit(0, "mm"),
        plot.title = element_text(hjust = 0.06, vjust = -15, size = 20),
        plot.subtitle = element_text(hjust = 0.02, vjust = -12, size = 20))

# mutations on branches
plot_dim_y <- ggplot_build(treeplot)$layout$panel_scales_y[[1]]$range$range[1]
nudge = plot_dim_y/3

treeplot = treeplot +
  geom_text(aes(x = branch,  label = aa_muts),
            size = tipLabSize - 0.5,
            nudge_y = nudge)

# Fix tip label clipping
plot_dim_x <- ggplot_build(treeplot)$layout$panel_scales_x[[1]]$range$range[2]


treeplot = treeplot +
  coord_cartesian(clip = 'off', expand = FALSE) +
  xlim(NA, ((0.40 * plot_dim_x) + plot_dim_x))

# add clades

treeplot = treeplot + add_clades(cladesFile, tree@data)
#
ggsave(arguments$output,
       plot = treeplot,
       device = cairo_pdf,
       width = output_size[1],
       height = output_size[2],
       units = "mm")

# output to svg. this is done to perserve text objects that seem to fail in inkscape
# current issue is that the font is not registered with svglite. this needs to be done with
# register_font(). see https://www.tidyverse.org/blog/2021/02/svglite-2-0-0/
# svglite::svglite(gsub(".pdf", ".svg", arguments$output),
#        width = 8.3,
#        height = 11.7)
# treeplot
# dev.off()

message(crayon::green("Done! Plotting was a success"))
