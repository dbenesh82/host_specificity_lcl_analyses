library(shiny)
library(ggplot2)
library(dplyr)
library(RColorBrewer)


data <- read.csv(file = "data_for_shiny.csv")
data <- mutate(data, maxLCL.c = as.character(maxLCL))

mypalette <- brewer.pal(5, "Set2")
theme.o <- theme_update(axis.text = element_text(colour="black", size = 15),
                        axis.title = element_text(colour="black", size = 18, face = "bold", lineheight=0.25),
                        axis.ticks = element_line(colour="black"),
                        panel.border = element_rect(colour = "black",fill=NA),
                        panel.grid.minor=element_blank(),
                        panel.grid.major=element_line(color="gray",linetype = "dotted"),
                        panel.background= element_rect(fill = NA))

tax.ranks <- c('genus', 'family', 'order', 'class', 'phylum') # for axis label



# define UI for slopegraph application
ui <- fluidPage(
   
   # Application title
   titlePanel("My First Shiny App"),
   
   # Input - LCL and host spec measure
   selectInput("trt", "Life cycle length", 
               choices = sort(unique(data$maxLCL.c))
               ),
   selectInput("hsi", "Host specificity measure",
               choices = list("Host range" = "hosts",
                           "Host specificity index" = "hs")
               ),
      
   # Show a plot
  mainPanel(
    plotOutput("slopesPlot")
      )
  )



server <- function(input, output) {

   output$slopesPlot <- renderPlot({
     
     if(input$hsi == "hosts") {
       plot <- ggplot(data = filter(data, hosts > 1, maxLCL.c != input$trt),
                      aes(x=factor(Host.no), y = hosts, group = Parasite.species)) +
         geom_line(alpha=0.25, color = "lightgray") +
         geom_point(size=3, alpha=0.2, color = "lightgray") +
         geom_point(data = filter(data, hosts > 1, maxLCL.c == input$trt),
                    size = 3, alpha = 0.4,
                    color = mypalette[ as.integer(input$trt) ]) +
         geom_line(data = filter(data, hosts > 1, maxLCL.c == input$trt),
                   alpha = 0.4, color = mypalette[ as.integer(input$trt) ]) +
         labs( x = "\nHost", y = "Host range\n") +
         scale_y_log10() +
         annotate('text', x = 4.5, y = 74.5,
                  label = paste(input$trt, '-host cycle', sep = ""),
                  color = mypalette[ as.integer(input$trt) ], size = 5)
     } else {
       plot <- ggplot(data = filter(data, hosts > 1, maxLCL.c != input$trt),
                      aes(x=factor(Host.no), y = hs, group = Parasite.species)) +
         geom_line(alpha=0.25, color = "lightgray") +
         geom_point(size=3, alpha = 0.2, color = "lightgray") +
         geom_point(data = filter(data, hosts > 1, maxLCL.c == input$trt),
                    size = 3, alpha = 0.4,
                    color = mypalette[ as.integer(input$trt) ]) +
         geom_line(data = filter(data, hosts > 1, maxLCL.c == input$trt),
                   alpha = 0.4, color = mypalette[ as.integer(input$trt) ]) +
         labs( x = "\nHost", y = "Host specificity index\n") +
         scale_y_continuous(limits = c(1,5), breaks = c(1:5), labels = tax.ranks) +
         annotate('text', x = 4.5, y = 4.5,
                  label = paste(input$trt, '-host cycle', sep = ""),
                  color = mypalette[ as.integer(input$trt) ], size = 5)
     }
     
     plot + 
       scale_x_discrete(expand=c(0.05,0.05))
       
   })
}

# Run the application 
shinyApp(ui = ui, server = server)
