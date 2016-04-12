## Plotting VEGAN multivariate ordinations with ggplot
## D. Levy-Booth 04/12/2016

############################################################

### Suite of functions to plot four major multivariate analysis techniques. 
###
### More information on using the "capscape" and "adonis" functions for 
### multivariate analysis by Jari Oksanen: 
### CAPSCALE: http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/capscale.html
### ADONIS: http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html
###
### 1. PCoA - Multidimentional Scaling of a vegan "capscale" object.
### Unconstrained multivariate ordination of dissimilarity matricies
### Capscale notation:
### pcoa <- capscale(otu_table ~ 1, dist="bray")
### pcoa.plot <- plot_pcoa(pcoa, metadata$color, metadata$shape, "Plot Title")
###
### 
### 2. db-RDA - distance-based Redundancy Analysis (vegan).
### Constrained multivariate ordination with arrows indicating regressor loadings
### rda <- capscale(otu_table ~ ., continuous_metadata, dist="bray")
### rda.plot <- plot_rda(rda, metadata$color, metadata$shape, "Plot Title")
###
### 3. CAP - Constrained Analysis of Principal Coordinates (vegan). 
### Constrained ordination of catagorical variables
### cap <- capscale(otu_table ~ ., catagorical_metadata, dist="bray")
### cap.plot <- plot_cap(cap, metadata$color, metadata$shape, "Plot Title")
###
### 4. PERMANOVA - plot barplot of adonis function PERMANOVA(vegan).
### Permutational, non-parametric multivariate analysis of variance
### Statistically test influence of catagoical and continuous variables on dissimilarity matricies
### adonis <- adonis(otu_table ~ ., metadata, dist="bray")
### adonis.plot <- adonis_plot(adonis, "Plot Title")
### 
### Note that color and shape vectors can be "NULL" 
###
####################################################################################

#PCOA Plot 
plot_pcoa <-function(pcoa, metacolor, metashape, titleby) {

eig<- eigenvals(pcoa)
prop<-eig/sum(eig)
PCOA1 <- paste("PCoA1 ",100*round(prop[1],3),"%")
PCOA2 <- paste("PCoA2 ",100*round(prop[2],3),"%")

pcoa.sum <- summary(pcoa)
pcoa.sum.sites  <- data.frame(pcoa.sum$sites[,1:2])       # PC1 and PC2
pcoa.sum.species  <- data.frame(pcoa.sum$species[,1:2])     # loadings for PC1 and PC2

if(!is.null(metacolor)) {
pcoa.sum.sites["colorby"] <- metacolor
}
else { pcoa.sum.sites["colorby"] = rep("1", nrow(pcoa.sum.sites)) }

if(!is.null(metashape)) {
pcoa.sum.sites["shapeby"] <- metashape
}
else { pcoa.sum.sites["shapeby"] = rep("1", nrow(pcoa.sum.sites)) }

pcoa_plot <- ggplot(pcoa.sum.sites, aes(x=MDS1, y=MDS2)) + 
  geom_point(aes(color = colorby, shape = shapeby), size=4) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme_bw() +
  theme( panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank(),
         text = element_text(size=14), 
         legend.background = element_blank(), 
         legend.title = element_blank(), 
         legend.key = element_blank()) +
  ggtitle(titleby) + 
  xlab (PCOA1) + 
  ylab (PCOA2) + 
  coord_fixed()

return(pcoa_plot)
}


#####################################################################
#db-RDA Plot 


plot_rda <-function(rda, metacolor, metashape, titleby) {

  eig<- eigenvals(rda)
  prop<-eig/sum(eig)
  RDA1 <- paste("db-RDA1 ",100*round(prop[1],3),"%")
  RDA2 <- paste("db-RDA2 ",100*round(prop[2],3),"%")

    
  rda.sum <- summary(rda)
  rda.sum.sites  <- data.frame(rda.sum$sites[,1:2])       # PC1 and PC2
  rda.sum.biplot  <- data.frame(rda.sum$biplot[,1:2])     # loadings for PC1 and PC2

  rda_labels <- rownames(rda.sum.biplot)
  rda.sum.biplot["labels"] <- rda_labels

  
  if(!is.null(metacolor)) {
    rda.sum.sites["colorby"] <- metacolor
  }
  else { rda.sum.sites["colorby"] = rep("1", nrow(rda.sum.sites)) }
  
  if(!is.null(metashape)) {
    rda.sum.sites["shapeby"] <- metashape
  }
  else { rda.sum.sites["shapeby"] = rep("1", nrow(rda.sum.sites)) }
  
  if(length(unique(metashape)) == 4) { 
    
    shapes <- c(15, 16, 17, 18)
  
  rda_plot1 <- ggplot(rda.sum.sites, aes(x=CAP1, y=CAP2)) + 
    geom_point(aes(color = colorby, shape = shapeby), size=4) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    theme_bw() +
    scale_shape_manual(values = shapes) + 
    theme( panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank(),
           text = element_text(size=14), 
           legend.background = element_blank(), 
           legend.title = element_blank(), 
           legend.key = element_blank()) +
    ggtitle(titleby) + 
    xlab (RDA1) + 
    ylab (RDA2) + 
    coord_fixed() +
    geom_segment(data=rda.sum.biplot, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
                 color="grey20", arrow=arrow(length=unit(0.015,"npc")))
  
  rda_plot2 <- rda_plot1 +
    geom_text(data=rda.sum.biplot, 
              aes(x=CAP1,y=CAP2, label=labels,
                  hjust=0.6*(1-sign(CAP1)),vjust=0.6*(1-sign(CAP2))), 
              color="grey20", size=5)
  
  }
  
  else if(length(unique(metashape)) == 5) { 
    
    shapes <- c(21, 22, 23, 24, 25)
    
    rda_plot1 <- ggplot(rda.sum.sites, aes(x=CAP1, y=CAP2)) + 
      geom_point(aes(color = colorby, fill= colorby, shape = shapeby), size=4) +
      geom_hline(yintercept=0, linetype="dotted") +
      geom_vline(xintercept=0, linetype="dotted") +
      theme_bw() +
      scale_shape_manual(values = shapes) + 
      theme( panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             text = element_text(size=14), 
             legend.background = element_blank(), 
             legend.title = element_blank(), 
             legend.key = element_blank()) +
      ggtitle(titleby) + 
      xlab ("db-RDA1") + 
      ylab ("db-RDA2") + 
      coord_fixed() +
      geom_segment(data=rda.sum.biplot, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
                   color="grey20", arrow=arrow(length=unit(0.015,"npc")))
    
    rda_plot2 <- rda_plot1 +
      geom_text(data=rda.sum.biplot, 
                aes(x=CAP1,y=CAP2, label=labels,
                    hjust=0.6*(1-sign(CAP1)),vjust=0.6*(1-sign(CAP2))), 
                color="grey20", size=5)
    
  }
  
  else {
  rda_plot1 <- ggplot(rda.sum.sites, aes(x=CAP1, y=CAP2)) + 
    geom_point(aes(color = colorby, shape = shapeby), size=4) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    theme_bw() +
    theme( panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank(),
           text = element_text(size=14), 
           legend.background = element_blank(), 
           legend.title = element_blank(), 
           legend.key = element_blank()) +
    ggtitle(titleby) + 
    xlab (RDA1) + 
    ylab (RDA2) + 
    coord_fixed() +
    geom_segment(data=rda.sum.biplot, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
                 color="grey20", arrow=arrow(length=unit(0.015,"npc")))

rda_plot2 <- rda_plot1 +
    geom_text(data=rda.sum.biplot, 
            aes(x=CAP1,y=CAP2, label=labels,
                hjust=0.6*(1-sign(CAP1)),vjust=0.6*(1-sign(CAP2))), 
            color="grey20", size=5)

  }
  
  return(rda_plot2)
}



##########################################################################
#PERMANOVA (adonis) Plot 


adonis_plot <- function(adonis, adonis_title) {
  len <- length(rownames(adonis$aov.tab)) -2
  
  p.sym <- function(Pr){ if (Pr <= 0.001) { x<-"***" } else if (Pr <= 0.01) {x<-"**" } else if (Pr <=0.05) { x<-"*" } else x<-" "
  return(as.character(x))}
  
  t_list <- rownames(adonis$aov.tab)[1:len]
  r_list <- adonis$aov.tab$R2[1:len]
  v_list <- r_list*100
  p_list <- adonis$aov.tab$Pr[1:len]
  ps_list <- sapply(p_list, p.sym)
  output <- data.frame(t_list, r_list, v_list, p_list, ps_list)
  colnames(output) <- c("Term", "R2", "variation", "P", "P_")
  output <- output[order(-output$variation),]
  ## set the levels in order we want
  output$Term_ordered <- factor(output$Term, levels=output$Term)
  
  residuals <- paste("(Res: ", round(adonis$aov.tab$R2[len+1], 3), ")")
  
  # graph time!
  adonis_bar <- ggplot(output, aes(y=variation, x=Term_ordered)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_bw() +
    theme( panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           text = element_text(size=14),
           panel.margin = unit(1.2, "lines"),
           legend.background =element_blank(),
           legend.position = c(0, 1), 
           legend.justification = c(0.1, 1), 
           axis.title.y = element_text(vjust=1),
           axis.title.x = element_blank(),
           axis.text.y = element_text(size = rel(1.4)),
           axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size = rel(1.4)),
           axis.title = element_text(size = rel(1.4)))+ 
    ylab("Explained Variation (%)") +
    ggtitle(paste(adonis_title, residuals, sep=" ")) +
    geom_text(aes(label=P_), vjust=1, color="white", size = 10)

  return(adonis_bar)
} 

#####################################################################
# CAP Plot

plot_cap <-function(rda, metacolor, metashape, titleby) {
  
  eig<- eigenvals(rda)
  prop<-eig/sum(eig)
  CAP1 <- paste("CAP1 ",100*round(prop[1],3),"%")
  CAP2 <- paste("CAP2 ",100*round(prop[2],3),"%")
  
  rda.sum <- summary(rda)
  rda.sum.sites  <- data.frame(rda.sum$sites[,1:2])       # PC1 and PC2
  rda.sum.biplot  <- data.frame(rda.sum$biplot[,1:2])     # loadings for PC1 and PC2
  
  rda_labels <- rownames(rda.sum.biplot)
  rda.sum.biplot["labels"] <- rda_labels
  
  if(!is.null(metacolor)) {
    rda.sum.sites["colorby"] <- metacolor
  }
  else { rda.sum.sites["colorby"] = rep("1", nrow(rda.sum.sites)) }
  
  if(!is.null(metashape)) {
    rda.sum.sites["shapeby"] <- metashape
  }
  else { rda.sum.sites["shapeby"] = rep("1", nrow(rda.sum.sites)) }
  
  if(length(unique(metashape)) == 4) { 
    
    shapes <- c(15, 16, 17, 18)
    
    rda_plot1 <- ggplot(rda.sum.sites, aes(x=CAP1, y=CAP2)) + 
      geom_point(aes(color = colorby, shape = shapeby), size=4) +
      geom_hline(yintercept=0, linetype="dotted") +
      geom_vline(xintercept=0, linetype="dotted") +
      theme_bw() +
      scale_shape_manual(values = shapes) + 
      theme( panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             text = element_text(size=14), 
             legend.background = element_blank(), 
             legend.title = element_blank(), 
             legend.key = element_blank()) +
      ggtitle(titleby) + 
      xlab (CAP1) + 
      ylab (CAP2) + 
      coord_fixed() +
      geom_segment(data=rda.sum.biplot, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
                   color="grey20", arrow=arrow(length=unit(0.015,"npc")))
    
    rda_plot2 <- rda_plot1 +
      geom_text(data=rda.sum.biplot, 
                aes(x=CAP1,y=CAP2, label=labels,
                    hjust=0.6*(1-sign(CAP1)),vjust=0.6*(1-sign(CAP2))), 
                color="grey20", size=5)
    
  }
  else {
    rda_plot1 <- ggplot(rda.sum.sites, aes(x=CAP1, y=CAP2)) + 
      geom_point(aes(color = colorby, shape = shapeby), size=4) +
      geom_hline(yintercept=0, linetype="dotted") +
      geom_vline(xintercept=0, linetype="dotted") +
      theme_bw() +
      theme( panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             text = element_text(size=14), 
             legend.background = element_blank(), 
             legend.title = element_blank(), 
             legend.key = element_blank()) +
      ggtitle(titleby) + 
      xlab (CAP1) + 
      ylab (CAP2) + 
      coord_fixed()
    
  }
  
  return(rda_plot1)
}


#########################################################################
