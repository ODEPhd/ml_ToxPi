library(ggplot2)
library(tidyr)

#########################################################################################

#Load desired datafile
LF <- read.csv(file.choose(),header = T)

# Manipulate the data (i.e., combine the areas of multiple sample region columns
#into one column and code them by each region)
LF1 <- gather(LF,ToxPi.slices, ToxPi.scores,2:10)

# Change the short names of the compounds to their full names
LF1$ToxPi.slices[LF1$ToxPi.slices=="Peak.area"]<-"Peak area"
LF1$ToxPi.slices[LF1$ToxPi.slices=="Detection.frequency"]<-"Detection frequency"
LF1$ToxPi.slices[LF1$ToxPi.slices=="Bioactivity.ratio"]<-"Bioactivity ratio"
LF1$ToxPi.slices[LF1$ToxPi.slices=="Exposure.category"]<-"Exposure category"
LF1$ToxPi.slices[LF1$ToxPi.slices=="IARC.catergory"]<-"IARC category"
LF1$ToxPi.slices[LF1$ToxPi.slices=="logkow"]<-"Log Kow"
LF1$ToxPi.slices[LF1$ToxPi.slices=="Water.solubility"]<-"Water solubility"
LF1$ToxPi.slices[LF1$ToxPi.slices=="Biodeg.half_life"]<-"Biodeg. half-life"
LF1$ToxPi.slices[LF1$ToxPi.slices=="Koc"]<-"Koc"

# Plot the basic stacked bar graph and decide the order of stacks
LF2 <- ggplot(LF1[order(LF1$ToxPi.slices, decreasing = T),],
               aes(fill = factor(ToxPi.slices, levels = c("Peak area", "Detection frequency",
                                                          "Bioactivity ratio", "Exposure category",
                                                          "IARC category", "Water solubility", "Log Kow",
                                                          "Biodeg. half-life", "Koc")), y = ToxPi.scores,
                   x = factor (Compound.names, levels = unique(Compound.names))))  +
  geom_bar(position = "fill", stat = "identity",width = 0.8) + scale_y_continuous(labels = scales::percent_format())+
  
  # Use a named vector with the group labels to manually specify a colour for each compound class
  #https://www.color-hex.com/color-palette
  scale_fill_manual(values = c("Peak area" = "#794044", "Detection frequency" = "#b167db", "Bioactivity ratio" = "#ffbf00",
                               "Exposure category" = "#19756a", "IARC category" = 	"#0000FF", "Water solubility" = "#35bfd7",
                               "Log Kow" = "#e77000", "Biodeg. half-life" = "#583660", "Koc" = "#38761d")) +
  
  # Change the background theme
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  
  # Change the x-axis text
  theme(axis.text.x = element_text(size=12, color = "Black", face="bold", angle=45, vjust = 1 , hjust = 1)) +
  
  # Change the y-axis text
  theme(axis.text.y = element_text(size=12, color = "Black", face="bold", angle=0, vjust = 0.5 , hjust = 1)) +    
  
  # Change the Legend position, title, and text
  theme(legend.position = "none", legend.direction = "horizontal" , legend.title = element_text(colour = "darkblue", size = 10, face = "bold")) +
  theme(legend.text = element_text(colour = "black", size = 10, face = "bold"))+
  
  # Group the compounds into separate classes and arrange the various classes in a particular order
  facet_grid(~factor(Config.level, levels = c('Level 1 (LF region)', 'Level 2 (LF region)')),scales = "free_x", space = "free_x") +
  
  # Change facet text font. Possible values for the font style:
  #'plain', 'italic', 'bold', 'bold.italic'.
  theme(strip.text.x = element_text(size=12, color="black",
                                    face="bold.italic", angle = 0)) +
  
  # Change the appearance of the rectangle around facet label
  theme(strip.background = element_rect(color="black", fill="white", 
                                        linewidth=1.5, linetype="solid")) +
  
  # # Change the Legend position, title, and text
  # theme(legend.position = "top", legend.title = element_text(color = "blue", size = 14, face = "bold")) +
  # theme(legend.text = element_text(color = "black", size = 14, face = "bold"))+
  # guides(fill = guide_legend(title = "ToxPi slices")) +
  
  # Add titles and axis labels to plot ### (\n is used to break words into different lines)
  # # ggtitle("Contribution of individual ToxPi slices to the total ToxPi score") +
  xlab("") +
  ylab("ml_ToxPi scores (Percent)") +
  
  # # Adjust style of plot titles ###
  # theme(plot.title = element_text(family ="serif", face = "bold", hjust = 0.5, vjust = 0.5, size = 20, color = "darkblue")) +
  
  ### Change font size of axis titles ###
  theme(axis.title.x = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))+
  theme(axis.title.y = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))


###############################################################################

#Load desired datafile
OC1 <- read.csv(file.choose(),header = T)

# Manipulate the data (i.e., combine the areas of multiple sample region columns into one column and code them by each region)
OC1_1 <- gather(OC1,ToxPi.slices, ToxPi.scores,2:10)

# Change the short names of the compounds to their full names
OC1_1$ToxPi.slices[OC1_1$ToxPi.slices=="Peak.area"]<-"Peak area"
OC1_1$ToxPi.slices[OC1_1$ToxPi.slices=="Detection.frequency"]<-"Detection frequency"
OC1_1$ToxPi.slices[OC1_1$ToxPi.slices=="Bioactivity.ratio"]<-"Bioactivity ratio"
OC1_1$ToxPi.slices[OC1_1$ToxPi.slices=="Exposure.category"]<-"Exposure category"
OC1_1$ToxPi.slices[OC1_1$ToxPi.slices=="IARC.catergory"]<-"IARC category"
OC1_1$ToxPi.slices[OC1_1$ToxPi.slices=="logkow"]<-"Log Kow"
OC1_1$ToxPi.slices[OC1_1$ToxPi.slices=="Water.solubility"]<-"Water solubility"
OC1_1$ToxPi.slices[OC1_1$ToxPi.slices=="Biodeg.half_life"]<-"Biodeg. half-life"
OC1_1$ToxPi.slices[OC1_1$ToxPi.slices=="Koc"]<-"Koc"


# Plot the basic stacked bar graph and decide the order of stacks
OC1_2 <- ggplot(OC1_1[order(OC1_1$ToxPi.slices, decreasing = T),],
              aes(fill = factor(ToxPi.slices, levels = c("Peak area", "Detection frequency",
                                                         "Bioactivity ratio", "Exposure category",
                                                         "IARC category", "Water solubility", "Log Kow",
                                                           "Biodeg. half-life", "Koc")), y = ToxPi.scores,
                  x = factor (Compound.names, levels = unique(Compound.names))))  +
      geom_bar(position = "fill", stat = "identity",width = 0.8) + scale_y_continuous(labels = scales::percent_format())+
 
# Use a named vector with the group labels to manually specify a colour for each compound class
#https://www.color-hex.com/color-palette
scale_fill_manual(values = c("Peak area" = "#794044", "Detection frequency" = "#b167db", "Bioactivity ratio" = "#ffbf00",
                             "Exposure category" = "#19756a", "IARC category" = 	"#0000FF", "Water solubility" = "#35bfd7",
                             "Log Kow" = "#e77000", "Biodeg. half-life" = "#583660", "Koc" = "#38761d")) +
  
# Change the background theme
  theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  
# Change the x-axis text
  theme(axis.text.x = element_text(size=12, color = "Black", face="bold", angle=45, vjust = 1 , hjust = 1)) +
  
# Change the y-axis text
  theme(axis.text.y = element_text(size=12, color = "Black", face="bold", angle=0, vjust = 0.5 , hjust = 1)) +    
  
# Change the Legend position, title, and text
  theme(legend.position = "right", legend.direction = "vertical" , legend.title = element_text(colour = "darkblue", size = 10, face = "bold")) +
  theme(legend.text = element_text(colour = "black", size = 10, face = "bold"))+
  
# Group the compounds into separate classes and arrange the various classes in a particular order
  facet_grid(~factor(Config.level, levels = c('Level 1 (OC1 region)', 'Level 2 (OC1 region)')),scales = "free_x", space = "free_x") +

# Change facet text font. Possible values for the font style:
#'plain', 'italic', 'bold', 'bold.italic'.
  theme(strip.text.x = element_text(size=12, color="black",
                                    face="bold.italic", angle = 0)) +

# Change the appearance of the rectangle around facet label
  theme(strip.background = element_rect(color="black", fill="white", 
                                        linewidth=1.5, linetype="solid")) +

# Change the Legend position, title, and text
  theme(legend.position = "right", legend.title = element_text(color = "blue", size = 14, face = "bold")) +
  theme(legend.text = element_text(color = "black", size = 14, face = "bold"))+
  guides(fill = guide_legend(title = "ToxPi slices")) +
  
# Add titles and axis labels to plot ### (\n is used to break words into different lines)
  # ggtitle("Contribution of individual ToxPi slices to the total ToxPi score") +
  xlab("") +
  ylab("") +
  
# Adjust style of plot titles ###
  theme(plot.title = element_text(family ="serif", face = "bold", hjust = 0.5, vjust = 0.5, size = 20, color = "darkblue")) +
  
### Change font size of axis titles ###
  theme(axis.title.x = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))+
  theme(axis.title.y = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))

 
# ### Save plot as a high quality png file ###
#    ggsave("Contribution of individual ToxPi slices_SP.png", width = 18, height = 10, dpi = 600)
  
   
    
   
 ############################################################################################
   
   #Load desired datafile
   OC2 <- read.csv(file.choose(),header = T)
   
   # Manipulate the data (i.e., combine the areas of multiple sample region columns into one column and code them by each region)
   OC2_1 <- gather(CK,ToxPi.slices, ToxPi.scores,2:10)
   
   # Change the short names of the compounds to their full names
   OC2_1$ToxPi.slices[OC2_1$ToxPi.slices=="Peak.area"]<-"Peak area"
   OC2_1$ToxPi.slices[OC2_1$ToxPi.slices=="Detection.frequency"]<-"Detection frequency"
   OC2_1$ToxPi.slices[OC2_1$ToxPi.slices=="Bioactivity.ratio"]<-"Bioactivity ratio"
   CK1OC2_1$ToxPi.slices[OC2_1$ToxPi.slices=="Exposure.category"]<-"Exposure category"
   OC2_1$ToxPi.slices[OC2_1$ToxPi.slices=="IARC.catergory"]<-"IARC category"
   OC2_1$ToxPi.slices[OC2_1$ToxPi.slices=="logkow"]<-"Log Kow"
   OC2_1$ToxPi.slices[OC2_1$ToxPi.slices=="Water.solubility"]<-"Water solubility"
   OC2_1$ToxPi.slices[OC2_1$ToxPi.slices=="Biodeg.half_life"]<-"Biodeg. half-life"
   OC2_1$ToxPi.slices[OC2_1$ToxPi.slices=="Koc"]<-"Koc"
   
   # Plot the basic stacked bar graph and decide the order of stacks
   OC2_2 <- ggplot(OC2_1[order(OC2_1$ToxPi.slices, decreasing = T),],
                  aes(fill = factor(ToxPi.slices, levels = c("Peak area", "Detection frequency",
                                                             "Bioactivity ratio", "Exposure category",
                                                             "IARC category", "Water solubility", "Log Kow",
                                                             "Biodeg. half-life", "Koc")), y = ToxPi.scores,
                      x = factor (Compound.names, levels = unique(Compound.names))))  +
     geom_bar(position = "fill", stat = "identity",width = 0.8) + scale_y_continuous(labels = scales::percent_format())+
     
     # Use a named vector with the group labels to manually specify a colour for each compound class
     #https://www.color-hex.com/color-palette
     scale_fill_manual(values = c("Peak area" = "#794044", "Detection frequency" = "#b167db", "Bioactivity ratio" = "#ffbf00",
                                  "Exposure category" = "#19756a", "IARC category" = 	"#0000FF", "Water solubility" = "#35bfd7",
                                  "Log Kow" = "#e77000", "Biodeg. half-life" = "#583660", "Koc" = "#38761d")) +
     
     # Change the background theme
     theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
     
     # Change the x-axis text
     theme(axis.text.x = element_text(size=12, color = "Black", face="bold", angle=45, vjust = 1 , hjust = 1)) +
     
     # Change the y-axis text
     theme(axis.text.y = element_text(size=12, color = "Black", face="bold", angle=0, vjust = 0.5 , hjust = 1)) +    
     
     # Change the Legend position, title, and text
     theme(legend.position = "none", legend.direction = "horizontal" , legend.title = element_text(colour = "darkblue", size = 10, face = "bold")) +
     theme(legend.text = element_text(colour = "black", size = 10, face = "bold"))+
     
     # Group the compounds into separate classes and arrange the various classes in a particular order
     facet_grid(~factor(Config.level, levels = c('Level 1 (OC2 region)', 'Level 2 (OC2 region)')),scales = "free_x", space = "free_x") +
     
     # Change facet text font. Possible values for the font style:
     #'plain', 'italic', 'bold', 'bold.italic'.
     theme(strip.text.x = element_text(size=12, color="black",
                                       face="bold.italic", angle = 0)) +
     
     # Change the appearance of the rectangle around facet label
     theme(strip.background = element_rect(color="black", fill="white", 
                                           linewidth=1.5, linetype="solid")) +
     
     # # Change the Legend position, title, and text
     # theme(legend.position = "top", legend.title = element_text(color = "blue", size = 14, face = "bold")) +
     # theme(legend.text = element_text(color = "black", size = 14, face = "bold"))+
     # guides(fill = guide_legend(title = "ToxPi slices")) +
     # 
     # # Add titles and axis labels to plot ### (\n is used to break words into different lines)
     # ggtitle("Contribution of individual ToxPi slices to the total ToxPi score") +
     xlab("") +
     ylab("ml_ToxPi scores (Percent)") +
     
     # # Adjust style of plot titles ###
     # theme(plot.title = element_text(family ="serif", face = "bold", hjust = 0.5, vjust = 0.5, size = 20, color = "darkblue")) +
      
     ### Change font size of axis titles ###
     theme(axis.title.x = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))+
     theme(axis.title.y = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))
   
   
   
   #########################################################################################
   
   #Load desired datafile
   AGR <- read.csv(file.choose(),header = T)
   
   # Manipulate the data (i.e., combine the areas of multiple sample region columns into one column and code them by each region)
   AGR1 <- gather(AGR,ToxPi.slices, ToxPi.scores,2:10)
   
   # Change the short names of the compounds to their full names
   AGR1$ToxPi.slices[AGR1$ToxPi.slices=="Peak.area"]<-"Peak area"
   AGR1$ToxPi.slices[AGR1$ToxPi.slices=="Detection.frequency"]<-"Detection frequency"
   AGR1$ToxPi.slices[AGR1$ToxPi.slices=="Bioactivity.ratio"]<-"Bioactivity ratio"
   AGR1$ToxPi.slices[AGR1$ToxPi.slices=="Exposure.category"]<-"Exposure category"
   AGR1$ToxPi.slices[AGR1$ToxPi.slices=="IARC.catergory"]<-"IARC category"
   AGR1$ToxPi.slices[AGR1$ToxPi.slices=="logkow"]<-"Log Kow"
   AGR1$ToxPi.slices[AGR1$ToxPi.slices=="Water.solubility"]<-"Water solubility"
   AGR1$ToxPi.slices[AGR1$ToxPi.slices=="Biodeg.half_life"]<-"Biodeg. half-life"
   AGR1$ToxPi.slices[AGR1$ToxPi.slices=="Koc"]<-"Koc"
   
   # Plot the basic stacked bar graph and decide the order of stacks
   AGR2 <- ggplot(AGR1[order(AGR1$ToxPi.slices, decreasing = T),],
                 aes(fill = factor(ToxPi.slices, levels = c("Peak area", "Detection frequency",
                                                            "Bioactivity ratio", "Exposure category",
                                                            "IARC category", "Water solubility", "Log Kow",
                                                            "Biodeg. half-life", "Koc")), y = ToxPi.scores,
                     x = factor (Compound.names, levels = unique(Compound.names))))  +
     geom_bar(position = "fill", stat = "identity",width = 0.8) + scale_y_continuous(labels = scales::percent_format())+
     
     # Use a named vector with the group labels to manually specify a colour for each compound class
     #https://www.color-hex.com/color-palette
     scale_fill_manual(values = c("Peak area" = "#794044", "Detection frequency" = "#b167db", "Bioactivity ratio" = "#ffbf00",
                                  "Exposure category" = "#19756a", "IARC category" = 	"#0000FF", "Water solubility" = "#35bfd7",
                                  "Log Kow" = "#e77000", "Biodeg. half-life" = "#583660", "Koc" = "#38761d")) +
     
     # Change the background theme
     theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
     
     # Change the x-axis text
     theme(axis.text.x = element_text(size=12, color = "Black", face="bold", angle=45, vjust = 1 , hjust = 1)) +
     
     # Change the y-axis text
     theme(axis.text.y = element_text(size=12, color = "Black", face="bold", angle=0, vjust = 0.5 , hjust = 1)) +    
     
     # Change the Legend position, title, and text
     theme(legend.position = "none", legend.direction = "horizontal" , legend.title = element_text(colour = "darkblue", size = 10, face = "bold")) +
     theme(legend.text = element_text(colour = "black", size = 10, face = "bold"))+
     
     # Group the compounds into separate classes and arrange the various classes in a particular order
     facet_grid(~factor(Config.level, levels = c('Level 1 (AGR region)', 'Level 2 (AGR region)')),scales = "free_x", space = "free_x") +
     
     # Change facet text font. Possible values for the font style:
     #'plain', 'italic', 'bold', 'bold.italic'.
     theme(strip.text.x = element_text(size=12, color="black",
                                       face="bold.italic", angle = 0)) +
     
     # Change the appearance of the rectangle around facet label
     theme(strip.background = element_rect(color="black", fill="white", 
                                           linewidth=1.5, linetype="solid")) +
     
     # Change the Legend position, title, and text
     theme(legend.position = "none", legend.title = element_text(color = "blue", size = 14, face = "bold")) +
     theme(legend.text = element_text(color = "black", size = 14, face = "bold"))+
     guides(fill = guide_legend(title = "ToxPi slices")) +

     # # Add titles and axis labels to plot ### (\n is used to break words into different lines)
     # # ggtitle("Contribution of individual ToxPi slices to the total ToxPi score") +
     xlab("") +
     ylab("") +
     
     # # Adjust style of plot titles ###
     # theme(plot.title = element_text(family ="serif", face = "bold", hjust = 0.5, vjust = 0.5, size = 20, color = "darkblue")) +
     
     ### Change font size of axis titles ###
     theme(axis.title.x = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))+
     theme(axis.title.y = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))
   
   
   #########################################################################################
   
   #Load desired datafile
   IND <- read.csv(file.choose(),header = T)
   
   # Manipulate the data (i.e., combine the areas of multiple sample region columns into one column and code them by each region)
   IND1 <- gather(IND,ToxPi.slices, ToxPi.scores,2:10)
   
   # Change the short names of the compounds to their full names
   IND1$ToxPi.slices[IND1$ToxPi.slices=="Peak.area"]<-"Peak area"
   IND1$ToxPi.slices[IND1$ToxPi.slices=="Detection.frequency"]<-"Detection frequency"
   IND1$ToxPi.slices[IND1$ToxPi.slices=="Bioactivity.ratio"]<-"Bioactivity ratio"
   IND1$ToxPi.slices[IND1$ToxPi.slices=="Exposure.category"]<-"Exposure category"
   IND1$ToxPi.slices[IND1$ToxPi.slices=="IARC.catergory"]<-"IARC category"
   IND1$ToxPi.slices[IND1$ToxPi.slices=="logkow"]<-"Log Kow"
   IND1$ToxPi.slices[IND1$ToxPi.slices=="Water.solubility"]<-"Water solubility"
   IND1$ToxPi.slices[IND1$ToxPi.slices=="Biodeg.half_life"]<-"Biodeg. half-life"
   IND1$ToxPi.slices[IND1$ToxPi.slices=="Koc"]<-"Koc"
   
   # Plot the basic stacked bar graph and decide the order of stacks
   IND2 <- ggplot(IND1[order(IND1$ToxPi.slices, decreasing = T),],
                 aes(fill = factor(ToxPi.slices, levels = c("Peak area", "Detection frequency",
                                                            "Bioactivity ratio", "Exposure category",
                                                            "IARC category", "Water solubility", "Log Kow",
                                                            "Biodeg. half-life", "Koc")), y = ToxPi.scores,
                     x = factor (Compound.names, levels = unique(Compound.names))))  +
     geom_bar(position = "fill", stat = "identity",width = 0.8) + scale_y_continuous(labels = scales::percent_format())+
     
     # Use a named vector with the group labels to manually specify a colour for each compound class
     #https://www.color-hex.com/color-palette
     scale_fill_manual(values = c("Peak area" = "#794044", "Detection frequency" = "#b167db", "Bioactivity ratio" = "#ffbf00",
                                  "Exposure category" = "#19756a", "IARC category" = 	"#0000FF", "Water solubility" = "#35bfd7",
                                  "Log Kow" = "#e77000", "Biodeg. half-life" = "#583660", "Koc" = "#38761d")) +
     
     # Change the background theme
     theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
     
     # Change the x-axis text
     theme(axis.text.x = element_text(size=12, color = "Black", face="bold", angle=45, vjust = 1 , hjust = 1)) +
     
     # Change the y-axis text
     theme(axis.text.y = element_text(size=12, color = "Black", face="bold", angle=0, vjust = 0.5 , hjust = 1)) +    
     
     # Change the Legend position, title, and text
     theme(legend.position = "none", legend.direction = "horizontal" , legend.title = element_text(colour = "darkblue", size = 10, face = "bold")) +
     theme(legend.text = element_text(colour = "black", size = 10, face = "bold"))+
     
     # Group the compounds into separate classes and arrange the various classes in a particular order
     facet_grid(~factor(Config.level, levels = c('Level 1 (IND region)', 'Level 2 (IND region)')),scales = "free_x", space = "free_x") +
     
     # Change facet text font. Possible values for the font style:
     #'plain', 'italic', 'bold', 'bold.italic'.
     theme(strip.text.x = element_text(size=12, color="black",
                                       face="bold.italic", angle = 0)) +
     
     # Change the appearance of the rectangle around facet label
     theme(strip.background = element_rect(color="black", fill="white", 
                                           linewidth=1.5, linetype="solid")) +
     
     # # Change the Legend position, title, and text
     # theme(legend.position = "top", legend.title = element_text(color = "blue", size = 14, face = "bold")) +
     # theme(legend.text = element_text(color = "black", size = 14, face = "bold"))+
     # guides(fill = guide_legend(title = "ToxPi slices")) +
     # 
     # # Add titles and axis labels to plot ### (\n is used to break words into different lines)
     # # ggtitle("Contribution of individual ToxPi slices to the total ToxPi score") +
     xlab("") +
     ylab("") +
     
     # # Adjust style of plot titles ###
     # theme(plot.title = element_text(family ="serif", face = "bold", hjust = 0.5, vjust = 0.5, size = 20, color = "darkblue")) +
     
     ### Change font size of axis titles ###
     theme(axis.title.x = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))+
     theme(axis.title.y = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))
   
 
   ##Combine all independedent plots into one 
   library(patchwork)
      
   (LF2 | OC1_2) / (OC2_2 | AGR2 | IND2) 
   
   # ### Save plot as a high quality png file ###
   ##(This will save to your selected work directory)
      ggsave("Contribution of individual ToxPi slices.png", width = 18, height = 10, dpi = 600)
      
      
      
################################################################################
      library(ggplot2)
      library(tidyr)
################################################################################
      
      
      #Load desired datafile
      GW_ToxPi <- read.csv(file.choose(),header = T)
      
      # Manipulate the data (i.e., combine the areas of multiple sample region columns into one column and code them by each region)
      GW_ToxPi1 <- gather(GW_ToxPi, Compound.class, ToxPi.rank.numbers,2:19)
      
      # Change the short names of the compounds to their full names
      GW_ToxPi1$Compound.class[GW_ToxPi1$Compound.class=="Alkyl_PAHs"]<-"Alkyl-PAHs"
      GW_ToxPi1$Compound.class[GW_ToxPi1$Compound.class=="Hetero_PAHs"]<-"Hetero-PAHs"
      GW_ToxPi1$Compound.class[GW_ToxPi1$Compound.class=="n_Alkanes"]<-"n-Alkanes"
      GW_ToxPi1$Compound.class[GW_ToxPi1$Compound.class=="Synthetic_musks"]<-"Synthetic musks"
      GW_ToxPi1$Compound.class[GW_ToxPi1$Compound.class=="Parent_PAHs"]<-"Parent PAHs"
      
      # Plot the basic stacked bar graph and decide the order of stacks
      GW_ToxPi2 <- ggplot(GW_ToxPi1[order(GW_ToxPi1$Compound.class, decreasing = T),],
                          aes(fill = factor(Compound.class, levels = c('Parent PAHs', 'Alkyl-PAHs', 'n-Alkanes', 'Hetero-PAHs', 'Hopanes',
                                                                       'Diamondoids', 'Sesquiterpenes', 'Synthetic musks', 'Phenols',
                                                                       'Plasticizers', 'BUVs', 'VMSs', 'OPFRs', 'Pesticides', 'Pharms',
                                                                       'PCPs', 'TPs', 'Others')), y = ToxPi.rank.numbers,
                                                                       x = factor (Groups, levels = unique(Groups))))  +
        geom_bar(position = "stack", stat = "identity",width = 0.6) +
        
        # Use a named vector with the group labels to manually specify a colour for each compound class
        #https://www.color-hex.com/color-palette
        scale_fill_manual(values = c('Parent PAHs'= "#1f77b4", 'Alkyl-PAHs'="#5c5c5c", 'n-Alkanes'="#CD0BBC",
                                     'Hetero-PAHs'='#61D04F', 'Hopanes'="#000000", 'Diamondoids'="#0000ff", 
                                     'Sesquiterpenes'="#ff0000", 'Synthetic musks'="#c6e2ff", 'Phenols'="#93AA00",
                                     'Plasticizers'="#794044", 'BUVs'="#40e0d0", 'VMSs'="#ffbf00", 'OPFRs'="#b167db",
                                     'Pesticides'="#e77000", 'Pharms'="#ffc3a0", 'PCPs'="#003366",
                                     'TPs'="#19756a", 'Others'="#800000")) +
        
        
        # Change the background theme
        theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
        
        # Change the x-axis text
        theme(axis.text.x = element_text(size=12, color = "Black", face="bold", angle=0, vjust = 0.5 , hjust = 0.5)) +

        # Change the y-axis text
        theme(axis.text.y = element_text(size=12, color = "Black", face="bold", angle=0, vjust = 0.5 , hjust = 0.5)) +    
        
        # Change the Legend position, title, and text
        theme(legend.position = "none", legend.direction = "horizontal" , legend.title = element_text(colour = "darkblue", size = 10, face = "bold")) +
        theme(legend.text = element_text(colour = "black", size = 10, face = "bold"))+
        
        # Group the compounds into separate classes and arrange the various classes in a particular order
        facet_grid(~factor(Region, levels = c('Landfill (LF) region', 'Oil-contaminated 1 (OC1) region',
                                              'Oil-contaminated 2 (OC2) region', 'Agricultural (AGR) region',
                                              'Industrial (IND) region')), scales = "free_x", space = "free_x") +
        
        # Change facet text font. Possible values for the font style:
        #'plain', 'italic', 'bold', 'bold.italic'.
        theme(strip.text.x = element_text(size=12, color="black",
                                          face="bold.italic", angle = 0)) +
        
        # Change the appearance of the rectangle around facet label
        theme(strip.background = element_rect(color="black", fill="white", 
                                              linewidth=1.5, linetype="solid")) +
        
        # Change the Legend position, title, and text
        theme(legend.position = "top", legend.title = element_text(color = "blue", size = 14, face = "bold")) +
        theme(legend.text = element_text(color = "black", size = 14, face = "bold"))+
        guides(fill = guide_legend(title = "Compounds category")) +
        
        # Add titles and axis labels to plot ### (\n is used to break words into different lines)
        # # ggtitle("Contribution of individual ToxPi slices to the total ToxPi score") +
        xlab("") +
        ylab("Number") +
        
        # # Adjust style of plot titles ###
        # theme(plot.title = element_text(family ="serif", face = "bold", hjust = 0.5, vjust = 0.5, size = 20, color = "darkblue")) +
        
        ### Change font size of axis titles ###
        theme(axis.title.x = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))+
        theme(axis.title.y = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))


      
      ################################################################################################
      
      # #Load desired datafile
      # GW_ToxPi <- read.csv(file.choose(),header = T)
       
      # # Manipulate the data (i.e., combine the areas of multiple sample region columns into one column and code them by each region)
      # GW_ToxPi1 <- gather(GW_ToxPi, Compound.class, ToxPi.rank.numbers,2:19)
      
      # Change the short names of the compounds to their full names
      # GW_ToxPi1$Compound.class[GW_ToxPi1$Compound.class=="Alkyl_PAHs"]<-"Alkyl-PAHs"
      # GW_ToxPi1$Compound.class[GW_ToxPi1$Compound.class=="Hetero_PAHs"]<-"Hetero-PAHs"
      # GW_ToxPi1$Compound.class[GW_ToxPi1$Compound.class=="n_Alkanes"]<-"n-Alkanes"
      # GW_ToxPi1$Compound.class[GW_ToxPi1$Compound.class=="Synthetic_musks"]<-"Synthetic musks"
      # GW_ToxPi1$Compound.class[GW_ToxPi1$Compound.class=="Parent_PAHs"]<-"Parent PAHs"
      #  
      # Plot the basic stacked bar graph and decide the order of stacks
      GW_ToxPi2_1 <- ggplot(GW_ToxPi1[order(GW_ToxPi1$Compound.class, decreasing = T),],
                    aes(fill = factor(Compound.class, levels = c('Parent PAHs', 'Alkyl-PAHs', 'n-Alkanes', 'Hetero-PAHs', 'Hopanes',
                                                                 'Diamondoids', 'Sesquiterpenes', 'Synthetic musks', 'Phenols',
                                                                 'Plasticizers', 'BUVs', 'VMSs', 'OPFRs', 'Pesticides', 'Pharms',
                                                                 'PCPs', 'TPs', 'Others')), y = ToxPi.rank.numbers,
                        x = factor (Groups, levels = unique(Groups))))  +
        geom_bar(position = "fill", stat = "identity",width = 0.6) + scale_y_continuous(labels = scales::percent_format())+
        
        # Use a named vector with the group labels to manually specify a colour for each compound class
        #https://www.color-hex.com/color-palette
        scale_fill_manual(values = c('Parent PAHs'= "#1f77b4", 'Alkyl-PAHs'="#5c5c5c", 'n-Alkanes'="#CD0BBC",
                                     'Hetero-PAHs'='#61D04F', 'Hopanes'="#000000", 'Diamondoids'="#0000ff", 
                                     'Sesquiterpenes'="#ff0000", 'Synthetic musks'="#c6e2ff", 'Phenols'="#93AA00",
                                     'Plasticizers'="#794044", 'BUVs'="#40e0d0", 'VMSs'="#ffbf00", 'OPFRs'="#b167db",
                                     'Pesticides'="#e77000", 'Pharms'="#ffc3a0", 'PCPs'="#003366",
                                     'TPs'="#19756a", 'Others'="#800000")) +
        
        
        # Change the background theme
        theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
        
        # Change the x-axis text
        theme(axis.text.x = element_text(size=12, color = "Black", face="bold", angle=0, vjust = 0.5 , hjust = 0.5)) +
        
        # Change the y-axis text
        theme(axis.text.y = element_text(size=12, color = "Black", face="bold", angle=0, vjust = 0.5 , hjust = 0.5)) +    
        
        # Change the Legend position, title, and text
        theme(legend.position = "none", legend.direction = "horizontal" , legend.title = element_text(colour = "darkblue", size = 10, face = "bold")) +
        theme(legend.text = element_text(colour = "black", size = 10, face = "bold"))+
        
        # Group the compounds into separate classes and arrange the various classes in a particular order
        facet_grid(~factor(Region, levels = c('Landfill (LF) region', 'Oil-contaminated 1 (OC1) region',
                                              'Oil-contaminated 2 (OC2) region', 'Agricultural (AGR) region',
                                              'Industrial (IND) region')), scales = "free_x", space = "free_x") +
        
        # Change facet text font. Possible values for the font style:
        #'plain', 'italic', 'bold', 'bold.italic'.
        theme(strip.text.x = element_text(size=12, color="black",
                                          face="bold.italic", angle = 0)) +
        
        # Change the appearance of the rectangle around facet label
        theme(strip.background = element_rect(color="black", fill="white", 
                                              linewidth=1.5, linetype="solid")) +
        
        # # Change the Legend position, title, and text
        # theme(legend.position = "top", legend.title = element_text(color = "blue", size = 14, face = "bold")) +
        # theme(legend.text = element_text(color = "black", size = 14, face = "bold"))+
        # guides(fill = guide_legend(title = "")) +
        # 
        # Add titles and axis labels to plot ### (\n is used to break words into different lines)
        # # ggtitle("Contribution of individual ToxPi slices to the total ToxPi score") +
        xlab("Priority ranking group") +
        ylab("Percentage (%)") +
        
        # # Adjust style of plot titles ###
        # theme(plot.title = element_text(family ="serif", face = "bold", hjust = 0.5, vjust = 0.5, size = 20, color = "darkblue")) +
        
        ### Change font size of axis titles ###
        theme(axis.title.x = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))+
        theme(axis.title.y = element_text(family ="arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue"))
      
  
      ##Combine all independedent plots into one    
      library(patchwork)
      
      GW_ToxPi2 / GW_ToxPi2_1 
      
      # ### Save plot as a high quality png file ###
      #(This will save to your selected work directory)
      ggsave("ToxPi pollutant distribution by groups.png", width = 16, height = 10, dpi = 600)
      