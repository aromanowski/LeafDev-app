library(shiny)
library(shinythemes)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)

###############
## Functions ##
###############

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd+
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Plots 2 graphs of a specific gene:
# logCPM and logFC TP
# requires Gene ID, each of the datasets

plotComposite <- function(geneID, gene.logCPM.data, tp.de.data, md, dataset, plot_types){
  
  if (length(plot_types) == 0) {plot_types <- "logCPM"}
  if (length(dataset) == 0) {dataset <- "white light"}
  
  if ("logCPM" %in% plot_types) {
    # Select gene to plot logCPM data
    gene_to_plot <- gene.logCPM.data %>%
      filter(Gene_ID == geneID,
             treatment %in% dataset)
    
    # summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
    tgc <- summarySE(gene_to_plot, measurevar="logCPM", groupvars=c("treatment", "day"))
    
    geneplot <- ggplot(data = tgc, aes(x = day, y = logCPM, color = treatment, group = treatment)) +
      geom_errorbar(aes(ymin=logCPM-se, ymax=logCPM+se), width=.1) +
      geom_line() +
      geom_point() +
      theme_bw() +
      scale_color_manual(values=c("grey", "#A00000", "#F08080")) +
      theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(colour = "grey20", size = 12),
            text = element_text(size = 14))+
      xlab("D.A.S.") + ylab("Normalised expression (logCPM)") +
      scale_x_discrete(breaks = unique(gene_to_plot$day), labels = unique(gene_to_plot$day)) +
      ggtitle(paste0(gene_to_plot$Gene_ID, " Symbol: ", md$Symbol[which(md$Gene_ID == unique(gene_to_plot$Gene_ID))]),
              subtitle = "Leaf blade development - expression values")
  }
  
  if ("logFC" %in% plot_types) {
    # Select gene to plot TP DE data
    if (length(dataset) == 1){
      if (dataset == "white light") {
        dataset <- "early EoD FR"  
      }
    } else if (length(dataset) == 0){
      dataset <- "early EoD FR"  
    }
    
    gene_to_plot <- tp.de.data %>%
      filter(Gene_ID == geneID,
             treatment %in% dataset)
    # summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
    tgc <- summarySE(gene_to_plot, measurevar="logFC", groupvars=c("treatment", "day"))
    geneplot2 <- ggplot(data = tgc, aes(x = day, y = logFC, color = treatment, group = treatment)) +
      # geom_errorbar(aes(ymin=logCPM-se, ymax=logCPM+se), width=.1) +
      geom_col() +
      # geom_point() +
      # ylim(0,(max(tgc$logCPM)+max(tgc$se)+1)) +
      theme_bw() +
      scale_color_manual(values=c("#A00000", "#F08080")) +
      theme(axis.text.x = element_text(colour = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(colour = "grey20", size = 12),
            text = element_text(size = 14))+
      xlab("D.A.S.") + ylab("logFC") +
      scale_x_discrete(breaks = unique(gene_to_plot$day), labels = levels(day)) +
      ggtitle(paste0(unique(gene_to_plot$Gene_ID), " Symbol: ", md$Symbol[which(md$Gene_ID == unique(gene_to_plot$Gene_ID))]),
              subtitle = "Leaf blade development - logFC vs WL") +
      facet_wrap(~treatment)
  }
  
  if (length(plot_types) > 1) {
    composite <- ggarrange(geneplot, geneplot2,
                           labels = c("A", "B"),
                           ncol = 2, nrow = 1)
    return(composite)
  } else if (plot_types == "logCPM") {
    return(geneplot)
  } else if (plot_types == "logFC") {
      return(geneplot2)
    }
  
  
}


####################################
## Load the Leaf Development data ##
####################################

# Load metadata
md <- read_tsv(file = "data/araport11_metadata.txt")

# Create descriptive factors
treatment <- factor(c("white light", "early EoD FR", "late EoD FR"), levels = c("white light", "early EoD FR", "late EoD FR"))
day = factor(c("d13", "d16", "d20"), levels = c("d13", "d16", "d20"))

# Load DE data
gene.de.data.d13 <- read_tsv("data/00-primordiaDE_full.tab")
gene.de.data.d16 <- read_tsv("data/10-blade_d16_FR06_vs_WLDE_full.tab")
gene.de.data.d20 <- read_tsv("data/11-blade_d20_FR06_vs_WLDE_full.tab")
gene.de.data.d20.lt <- read_tsv("data/12-blade_d20_FR18_vs_WLDE_full.tab")

# Add treatment and day information to the data
gene.de.data.d13 <- gene.de.data.d13[,c(1,12:15)] %>% mutate(treatment = treatment[2], day = day[1])
gene.de.data.d16 <- gene.de.data.d16[,c(1,12:15)] %>% mutate(treatment = treatment[2], day = day[2])
gene.de.data.d20 <- gene.de.data.d20[,c(1,12:15)] %>% mutate(treatment = treatment[2], day = day[3])
gene.de.data.d20.lt <- gene.de.data.d20.lt[,c(1,12:15)] %>% mutate(treatment = treatment[3], day = day[3])

# Generate combined data for logFC plots
tp.de.data <- rbind(gene.de.data.d13, gene.de.data.d16, gene.de.data.d20, gene.de.data.d20.lt)

# Load logCPM expression data
gene.logCPM.data <- read_tsv("data/blade.genes.logCPM.normalised.values.tab")

# Reconstruct day by day data and add treatment and day information to the data
gene.logCPM.WL.d13.data <- rbind(gene.logCPM.data[,c(1,2)] %>% dplyr::rename(logCPM = Pr_D13_WL_01) %>% mutate(treatment = treatment[1], day = day[1])
                                 ,gene.logCPM.data[,c(1,3)] %>% dplyr::rename(logCPM = Pr_D13_WL_02) %>% mutate(treatment = treatment[1], day = day[1])
)

gene.logCPM.WL.d16.data <- rbind(gene.logCPM.data[,c(1,6)] %>% dplyr::rename(logCPM = B_D16_WL_01) %>% mutate(treatment = treatment[1], day = day[2])
                                 ,gene.logCPM.data[,c(1,7)] %>% dplyr::rename(logCPM = B_D16_WL_02) %>% mutate(treatment = treatment[1], day = day[2])
)
gene.logCPM.WL.d20.data <- rbind(gene.logCPM.data[,c(1,10)] %>% dplyr::rename(logCPM = B_D20_WL_01) %>% mutate(treatment = treatment[1], day = day[3])
                                 ,gene.logCPM.data[,c(1,11)] %>% dplyr::rename(logCPM = B_D20_WL_02) %>% mutate(treatment = treatment[1], day = day[3])
)
gene.logCPM.et.d13.data <- rbind(gene.logCPM.data[,c(1,4)] %>% dplyr::rename(logCPM = Pr_D13_FR06_01) %>% mutate(treatment = treatment[2], day = day[1])
                                 ,gene.logCPM.data[,c(1,5)] %>% dplyr::rename(logCPM = Pr_D13_FR06_02) %>% mutate(treatment = treatment[2], day = day[1])
)
gene.logCPM.et.d16.data <- rbind(gene.logCPM.data[,c(1,8)] %>% dplyr::rename(logCPM = B_D16_FR06_01) %>% mutate(treatment = treatment[2], day = day[2])
                                 ,gene.logCPM.data[,c(1,9)] %>% dplyr::rename(logCPM = B_D16_FR06_02) %>% mutate(treatment = treatment[2], day = day[2])
)
gene.logCPM.et.d20.data <- rbind(gene.logCPM.data[,c(1,12)] %>% dplyr::rename(logCPM = B_D20_FR06_01) %>% mutate(treatment = treatment[2], day = day[3])
                                 ,gene.logCPM.data[,c(1,13)] %>% dplyr::rename(logCPM = B_D20_FR06_02) %>% mutate(treatment = treatment[2], day = day[3])
)
gene.logCPM.lt.d20.data <- rbind(gene.logCPM.data[,c(1,14)] %>% dplyr::rename(logCPM = B_D20_FR18_01) %>% mutate(treatment = treatment[3], day = day[3])
                                 ,gene.logCPM.data[,c(1,15)] %>% dplyr::rename(logCPM = B_D20_FR18_02) %>% mutate(treatment = treatment[3], day = day[3])
)

# Generate combined ts data for logCPM plots
gene.logCPM.data <- rbind(gene.logCPM.WL.d13.data,
                          gene.logCPM.WL.d16.data,
                          gene.logCPM.WL.d20.data,
                          gene.logCPM.et.d13.data,
                          gene.logCPM.et.d16.data,
                          gene.logCPM.et.d20.data,
                          gene.logCPM.lt.d20.data)


# Create Shiny Interface

ui <- fluidPage(
  
  headerPanel('Arabidopsis leaf 3 gene expression dataset'),
  
  # Input() functions
  sidebarPanel(
    checkboxGroupInput("dataset", label = 'Select light conditions', choices = c("White light" = "white light",
                                                                               "EoD FR since day 06" = "early EoD FR",
                                                                               "EoD FR since day 18" = "late EoD FR"),
                       selected = c("white light", "early EoD FR")),
    checkboxGroupInput("plot_types", label = 'Select graph types', choices = c("LogCPM" = "logCPM",
                                                                               "LogFC" = "logFC"),
                       selected = "logCPM"),
    textInput("gene", label = "Gene Id (TAIR)", placeholder = "AT1G01060", value = "AT1G01060"),
    actionButton("Draw", "Plot this gene!")
  ),
  
  # Output() functions
  mainPanel(
    plotOutput('leaf_plot'),
    textOutput("text_out")
  )
)

server <- function(input, output, session) {
  observeEvent(input$Draw, {
      
    #  output$text_out <- renderText({
    #   paste("You chose: ",
    #         paste("Dataset: ", input$dataset,"\n", collapse = ", "),
    #         paste("Plot Type: ", input$plot_types, "\n", collapse = ", "),
    #         paste("Gene Id: ", input$gene),
    #         paste("Class:", class(input$dataset), "length: ", length(input$dataset)))
    # })
    
    output$leaf_plot <- renderPlot({
      par(mar = c(5.1, 4.1, 0, 1))
      plotComposite(input$gene, gene.logCPM.data, tp.de.data, md, input$dataset, input$plot_types)
    }, height = function() {
      session$clientData$output_leaf_plot_width/2
    })
    
  })
}

shinyApp(ui = ui, server = server)
