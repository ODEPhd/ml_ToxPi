library(ggplot2)
library(gridExtra)


# Function for plotting the variable importance scores
mplot_importance <- function(var, imp, colours = NA, limit = 15, model_name = NA, subtitle = NA,
                             save = FALSE, file_name = "viz_importance.png", subdir = NA) {
  
  require(ggplot2)
  require(gridExtra)
  options(warn=-1)
  
  if (length(var) != length(imp)) {
    message("The variables and importance values vectors should be the same length.")
    stop(message(paste("Currently, there are",length(var),"variables and",length(imp),"importance values!")))
  }
  if (is.na(colours)) {
    colours <- "#0274BD" 
  }
  out <- data.frame(var = var, imp = imp, Type = colours)
  if (length(var) < limit) {
    limit <- length(var)
  }
  
  output <- out[1:limit,]
  
  p <- ggplot(output, 
              aes(x = reorder(var, imp), y = imp, 
                  label = round(imp, 1))) + 
    geom_col(aes(fill = Type), width = 0.3) +
    geom_point(aes(colour = Type), size = 10) + 
    coord_flip() + xlab('') +
    # Change the background theme
    theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    ylab('RF model VI scores') +
    xlab ('ToxPi variables') +
    theme(axis.title.x = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))+
    theme(axis.title.y = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))+
    # Asjust the size, face, and color of the plot points  
   geom_text(hjust = 0.5, size = 4, inherit.aes = TRUE, colour = "white")+
    geom_text(hjust = 0.5, size = 4.07, inherit.aes = TRUE, colour = "white") +
    geom_text(hjust = 0.5, size = 4.08, inherit.aes = TRUE, colour = "white") +
    # Change the x-axis text
    theme(axis.text.x = element_text(size=14, color = "Black", face="bold", angle=0, vjust = 0.5 , hjust = 0.5)) +
    # Change the y-axis text
    theme(axis.text.y = element_text(size=14, color = "Black", face="bold", angle=0, vjust = 0.5 , hjust = 1)) +
    labs(title ="")
    # labs(title = paste0("Variables Importances. (", limit, " / ", length(var), " plotted)"))
  
  if (length(unique(output$Type)) == 1) {
    p <- p + geom_col(fill = colours, width = 0.3) +
      geom_point(colour = colours, size = 10) + 
      guides(fill = FALSE, colour = FALSE) + 
      geom_text(hjust = 0.5, size = 4, inherit.aes = TRUE, colour = "white")+
      geom_text(hjust = 0.5, size = 4.07, inherit.aes = TRUE, colour = "white")+
      geom_text(hjust = 0.5, size = 4.08, inherit.aes = TRUE, colour = "white")
     }
  if(!is.na(model_name)) {
    p <- p + labs(caption = model_name)
  }
  if(!is.na(subtitle)) {
    p <- p + labs(subtitle = subtitle)
  }  
  if(save == TRUE) {
    if (!is.na(subdir)) {
      dir.create(file.path(getwd(), subdir))
      file_name <- paste(subdir, file_name, sep="/")
    }
    p <- p + ggsave(file_name, width=7, height=6)
  }
  
  return(p)
  
}


# Example dataset with new variables and importances
ToxPi_variables <- c("Detection frequency", "Water solubility", "Average peak areas",
                   "Half-lives", "Log Kow", "Exposure category", "Koc",
                   "Bioactivity ratio", "IARC category")
rf_importances <- c(100.00, 79.19, 60.97, 49.92, 37.25, 35.59, 28.68, 22.80, 1.00)

# Call the mplot_importance function with the new dataset
plot <- mplot_importance(var = ToxPi_variables, 
                         imp = rf_importances,
                         colours = "#0274BD", 
                         limit = 9,
                         model_name = "",
                         subtitle = "",
                         save = FALSE)
plot

#Save plot
ggsave("varImp2.png", width = 7, height = 6, dpi = 600)



############################################################
# Plot the ML-optimized ToxPi weighting factors
############################################################

# Load necessary libraries
library(ggplot2)

# Create a dataframe from the provided data
data <- data.frame(
  Variable = c("Average peak area", "Detection frequency", "Bioactivity ratio", 
               "Exposure category", "IARC category", "Log Kow", 
               "Water solubility", "Biodegradation half-life", "Koc"),
  Weight = c(14.9, 23.8, 5, 8.9, 1, 8.9, 18.8, 11.9, 6.9)
)

# Calculate percentages
data$Percentage <- (data$Weight / sum(data$Weight)) * 100

# Reorder the dataframe based on the percentages
data$Variable <- factor(data$Variable, levels = data$Variable[order(data$Percentage)])

# Define a custom color palette
custom_colors <- c("#FFA07A","#ffbf00", "#e77000",  "#32CD32", "#FF69B4",
                   "#9370DB", "#1E90FF" ,"#40E0D0", "#1f77b4")

# Create the plot
plot <- ggplot(data, aes(x = Variable, y = Percentage, fill = Variable)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "",
       x = "ToxPi variables",
       y = "Weighting factors (%)") +
  theme(axis.title.x = element_text(size = 14, face = "bold", color = "darkblue"),
        axis.title.y = element_text(size = 14, face = "bold", color = "darkblue"),
        axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.position = "none") +
  geom_text(aes(label = round(Percentage, 1)), 
            position = position_stack(vjust = 0.5), size = 5, color = "black", fontface = "bold") +
  ylim(0, 25)

# Display the plot
print(plot)

#Save plot
ggsave("ML-ptimized ToxPi weighting factors.png", width = 7, height = 6, dpi = 600)
