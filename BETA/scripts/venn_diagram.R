# To plot upregulated gene sets across UV timepoints
library(VennDiagram)

path <- '/work/daylab/ella/rna-seq-analysis/beta-plus-0.05/filtered-target-predictions/'
df_01hr <- read.csv(paste0(path,'filtered_01hr_uptarget.csv'))
df_24hr <- read.csv(paste0(path,'filtered_24hr_uptarget.csv'))
df_48hr <- read.csv(paste0(path,'filtered_48hr_uptarget.csv'))

# Get gene sets for each timepoint
uptarget_01hr <- as.list(df_01hr$GeneSymbol)
uptarget_24hr <- as.list(df_24hr$GeneSymbol)
uptarget_48hr <- as.list(df_48hr$GeneSymbol)

# Generate venn diagram
venn.plot <- venn.diagram(
  x = list(one_hr = uptarget_01hr, twenty_four_hr = uptarget_24hr, forty_eight_hr = uptarget_48hr),
  category.names = c('01 hr', '24 hr', '48 hr'),
  filename = NULL,
  output = TRUE,
  print.mode = c('raw'),
  fill = c('pink', 'lightblue', 'lightgreen'),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  label = TRUE,
  main = 'Uptargeted genes across UV timepoints'
)

# Display the Venn diagram in the plot pane
grid.draw(venn.plot)