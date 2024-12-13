
###################################
# R version 4.4.1                 #
# Operating System: Windows 10    #
# Code for data Visualization     #
# Supplement to: Li, C., Dallmeyer, A. & Herzschuh, U. Magnitude of biome shifts over the last 21,000 years reveals the impact of seasonality and intense warming on future vegetation (2024) #
# Contact: Chenzhi Li (chenzhi.li@awi.de)[Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research, Potsdam, Germany 2024] #
###################################

# Figure list:
# 1. Fig. 1: Spatial distribution of the fossil pollen records analyzed showing the proportion of biome shifts per site over the last 21,000 years
# 2. Fig. 2: Temporal patterns of biome and bioclimate change in the Northern Hemisphere over the last 21,000 years
# 3. Fig. 3: Gaussian kernel correlations between changes in biomes and composite bioclimatic variables
# 4. Fig. 4: Potential spatial patterns of bioclimates and biome changes in the MPI-ESM simulation over the next 300 years in the Northern Hemisphere
# 5. Supplementary Fig. 1: Record density for the last 21,000 years
# 6. Supplementary Fig. 2: Spatial patterns of biome distributions at 0 cal. ka BP and their agreement with modern potential natural biomes
# 7. Supplementary Fig. 3: Temporal patterns of simulated bioclimatic variables over the last 21,000 years in the Northern Hemisphere
# 8. Supplementary Fig. 4: Gaussian kernel correlations between compositional turnover and proportion of shifts in biomes over the last 21,000 years in the Northern Hemisphere
# 9. Supplementary Fig. 5: Gaussian kernel correlation coefficients between compositional turnover and proportion of shifts in biomes with simulated bioclimatic variables and their rates of change over the last 21,000 years in the Northern Hemisphere 

# Note: Biomes and their abbreviations and codes 
## 1 - Tropical forest (TRFO); 2 - Subtropical forest (WTFO); 3 - Temperate forest (TEFO); 4 - Boreal forest (BOFO); 
## 5 - (Warm) savanna and dry woodland (SAVA); 6 - Grassland and dry shrubland (STEP); 7 - (Warm) desert (DESE); 8 - Tundra and polar desert (TUND)

# Note: Bioclimatic variables
## BIO1 = Annual Mean Temperature; BIO4 = Temperature Seasonality (standard deviation ×100); BIO8 = Mean Temperature of Wettest Quarter; BIO9 = Mean Temperature of Driest Quarter
## BIO10 = Mean Temperature of Warmest Quarter; BIO11 = Mean Temperature of Coldest Quarter; BIO12 = Annual Precipitation; BIO13 = Precipitation of Wettest Month
## BIO14 = Precipitation of Driest Month; BIO15 = Precipitation Seasonality (Coefficient of Variation); BIO16 = Precipitation of Wettest Quarter; BIO17 = Precipitation of Driest Quarter
## BIO18 = Precipitation of Warmest Quarter; BIO19 = Precipitation of Coldest Quarter

# Note: Composite bioclimatic variables
## Bioclimate seasonality (seasonality) = BIO4 + BIO15; Non-growing season temperature (Tnon-growing) = BIO1 + BIO9 + BIO11; Growing season temperature (Tgrowing) = BIO8 + BIO10
## Growing season precipitation (Pgrowing) = BIO12 + BIO13 + BIO16 + BIO18; Non-growing season precipitation (Pnon-growing) = BIO14 + BIO17 + BIO19

# choose directory:
#setwd("~/Supplementary code/Visualization")  # select the folder "Supplementary code" from your directory

# install packages if not installed
install.packages(c("ggplot2", "rnaturalearth", "rnaturalearthdata", "cowplot", "ggpubr", "dplyr", "tidyr", "dtwclust", "ggh4x"))

# loading packages
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)
library(ggpubr)
library(dplyr)
library(tidyr)
library(dtwclust)
library(ggh4x)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 1. Fig. 1: Spatial distribution of the fossil pollen records analyzed showing the proportion of biome shifts per site over the last 21,000 years  ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
Biome_shift_proportion_spatial_0_21ka      <- read.csv2("data/Fig.1 data.1-The proportion of biome shifts per site over the last 21,000 year.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_potential_natural_biomes_5minutes   <- read.csv2("data/Fig.1 data.2-Modern potential natural biomes (spatial resolution 5 arc minutes).csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Biome_shift_proportion_spatial_0_21ka     <- type.convert(Biome_shift_proportion_spatial_0_21ka, as.is = TRUE) 
Modern_potential_natural_biomes_5minutes  <- type.convert(Modern_potential_natural_biomes_5minutes, as.is = TRUE) 


# Plot
{
  # Convert biome data into factors with specified levels
  Modern_potential_natural_biomes_5minutes$Biome<- factor(Modern_potential_natural_biomes_5minutes$Biome,
                                                          levels = c("1", "2", "3", "4", "5", "6", "7", "8")) 
  
  # Limit reconstruction values to 95th percentile for better visualization
  Biome_shift_proportion_spatial_0_21ka$Reconstruction[Biome_shift_proportion_spatial_0_21ka$Reconstruction >= quantile(Biome_shift_proportion_spatial_0_21ka$Reconstruction, probs = 0.95)] <- quantile(Biome_shift_proportion_spatial_0_21ka$Reconstruction, probs = 0.95)
  
  # Scale simulation data and limit to 95th percentile
  Biome_shift_proportion_spatial_0_21ka$Simulation <- Biome_shift_proportion_spatial_0_21ka$Simulation*2
  Biome_shift_proportion_spatial_0_21ka$Simulation[Biome_shift_proportion_spatial_0_21ka$Simulation >= quantile(Biome_shift_proportion_spatial_0_21ka$Simulation, probs = 0.95)] <- quantile(Biome_shift_proportion_spatial_0_21ka$Simulation, probs = 0.95)
  
  # Map settings and customization
  {
    # Define x and y axis labels for the map
    xlabs = seq(-150,150,50)
    ylabs = seq(-60,80,20)
    
    #define the cols:
    cols <- c("1"  = "#e31a1c",  "2"  = "#f781bf",  "3" = "#33a02c",   "4" = "#1f78b4",
              "5"  = "#ff7f00",  "6"  = "#fdbf6f",  "7"  = "#ffed6f",  "8" = "#6a3d9a")
    
    # Define colors for biome categories
    # !!! IMPORTANT: USE THE SAME ORDER LIKE THE COLORS !!!
    brks=c("1" ,  "2" ,  "3",  "4", "5" ,  "6" ,  "7",  "8")
    
    # Define biome categories, breaks, and labels
    # !!! IMPORTANT: USE THE SAME ORDER LIKE THE COLORS !!!
    labs=c( "Tropical forest",                  "Subtropical forest",           "Temperate forest",    "Boreal forest",                   
            "(Warm) savanna and dry woodland",  "Grassland and dry shrubland",  "(Warm) desert",       "Tundra and polar desert")
    
    # Load coastline data
    coastline <- ne_coastline(scale = "medium", returnclass = "sf")
    
  }
  
  ## (1) Pollen-based reconstruction
  {
    Biome_shift_proportion_spatial_pollen_0_21ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= Modern_potential_natural_biomes_5minutes, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols, guide = "none") +
      labs(title= 'Pollen-based reconstruction', x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = xlabs, labels = paste0(xlabs, "°E"), expand = c(0, 0)) +
      geom_point(data= Biome_shift_proportion_spatial_0_21ka, aes(x=Longitude, y=Latitude, color= Reconstruction), size = 0.8) +
      scale_colour_gradient(name = "Proportion", low = "yellow", high = "red", limits = c(0, max(Biome_shift_proportion_spatial_0_21ka$Reconstruction, Biome_shift_proportion_spatial_0_21ka$Simulation)), breaks = seq(0,0.6, by =0.2), labels = seq(0,0.6, by =0.2))+  
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt")),
            axis.text.x=element_text(size=14)) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
            axis.text.y=element_text(size=14)) +
      theme(legend.position="none",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
      guides(color =guide_colorbar(direction = "vertical", barwidth = 1, barheight = 20, title.position = "top", title.hjust = 0.5, label.hjust = 0.5, label.theme = element_text(size = 15)))
    
    Biome_shift_proportion_spatial_pollen_0_21ka_map
  }
  
  # (2) MPI-ESM_V3_GLAC1D
  {
    Biome_shift_proportion_spatial_MPI_0_21ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= Modern_potential_natural_biomes_5minutes, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols, guide = "none") +
      labs(title= 'ESM-based simulation', x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = xlabs, labels = paste0(xlabs, "°E"), expand = c(0, 0)) +
      geom_point(data= Biome_shift_proportion_spatial_0_21ka, aes(x=Longitude, y=Latitude, color= Simulation), size = 0.8, show.legend = FALSE) +
      scale_colour_gradient(name = "Proportion", low = "yellow", high = "red", limits = c(0, max(Biome_shift_proportion_spatial_0_21ka$Reconstruction, Biome_shift_proportion_spatial_0_21ka$Simulation)), breaks = seq(0,0.6, by =0.2), labels = seq(0,0.6, by =0.2))+  
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt")),
            axis.text.x=element_text(size=14)) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
            axis.text.y=element_text(size=14)) +
      theme(legend.position="bottom",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
      guides(fill=guide_legend(override.aes=list(), nrow=2, byrow=TRUE))
    
    Biome_shift_proportion_spatial_MPI_0_21ka_map
  }
  
  # (3) share legend
  {
    # Create a dummy dataset just to generate the legend
    dummy_data <- data.frame(
      Longitude = c(0),
      Latitude = c(0),
      Proportion = c(Biome_shift_proportion_spatial_0_21ka$Reconstruction, Biome_shift_proportion_spatial_0_21ka$Simulation)
    )
    
    # Generate a plot with only the color legend
    color_legend_plot <- ggplot(dummy_data, aes(x = Longitude, y = Latitude, color = Proportion)) +
      geom_point(size = 0) +  # Generates the legend but doesn't actually draw points
      scale_colour_gradient(name = "Proportion", low = "yellow", high = "red", limits = c(0, max(dummy_data$Proportion)), breaks = seq(0,0.6, by =0.2), labels = seq(0,0.6, by =0.2))+  
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
      theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt")),
            axis.text.x=element_text(size=14)) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
            axis.text.y=element_text(size=14)) +
      theme(legend.position="right",
            legend.text=element_text(size=14),
            legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
      guides(color =guide_colorbar(direction = "vertical", barwidth = 1, barheight = 25, title.position = "top", title.hjust = 0.5, label.hjust = 0.5, label.theme = element_text(size = 15)))
    
    # Extract only the legend
    color_legend <- get_legend(color_legend_plot)
    
    # Arrange the legend as a standalone plot
    legend_only_plot <- plot_grid(NULL, color_legend, ncol = 1, rel_heights = c(0, 1))
    legend_only_plot
    
  }
  
  
  # -------------------------------------------------------------------------------------------------
  
  # Arranging plots
  blank_plot <- ggplot() + theme_void()
  
  Biome_shift_proportion_spatial_pollen_MPI_0_21ka_map  <- ggarrange(ggarrange(Biome_shift_proportion_spatial_pollen_0_21ka_map, Biome_shift_proportion_spatial_MPI_0_21ka_map, 
                                                                               nrow = 2, ncol = 1, common.legend = F, heights = c(0.4675, 0.5325)),
                                                                     legend_only_plot, nrow = 1, ncol = 2, widths = c(0.9, 0.1), common.legend = F)
  
  # Save the final figure as a PNG file
  ggsave("result/Fig. 1-Spatial pattern of proportion of biome shifts over the last 21,000 years.png", width = 15, height = 13, units = "in", dpi = 300)
  
}



# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 2. Fig. 2: Temporal patterns of biome and bioclimate change in the Northern Hemisphere over the last 21,000 years  ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
Biome_compositional_turnover_NH_0_21ka  <- read.csv2("data/Fig.2 data.1-Biome compositional turnover over the last 21,000 year in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Biome_shift_proportion_NH_0_21ka        <- read.csv2("data/Fig.2 data.2-Proportion of biome shifts over the last 21,000 year in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Bioclimate_state_NH_0_21ka              <- read.csv2("data/Fig.2 data.3-Simulted bioclimate state over the last 21,000 year in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Bioclimate_change_rate_NH_0_21ka        <- read.csv2("data/Fig.2 data.4-Simulted bioclimate change rate over the last 21,000 year in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Biome_compositional_turnover_NH_0_21ka  <- type.convert(Biome_compositional_turnover_NH_0_21ka, as.is = TRUE) 
Biome_shift_proportion_NH_0_21ka        <- type.convert(Biome_shift_proportion_NH_0_21ka , as.is = TRUE) 
Bioclimate_state_NH_0_21ka              <- type.convert(Bioclimate_state_NH_0_21ka, as.is = TRUE) 
Bioclimate_change_rate_NH_0_21ka        <- type.convert(Bioclimate_change_rate_NH_0_21ka, as.is = TRUE) 

# Plot
{
  # (1) Biome compositional turnover
  {
    # Convert 'Biome' into a factor with specified order for plotting purposes
    Biome_compositional_turnover_NH_0_21ka$Biome <- factor(Biome_compositional_turnover_NH_0_21ka$Biome,
                                                           levels = c("NH", "TRFO", "WTFO", "TEFO", "BOFO", "SAVA", "STEP", "DESE", "TUND"))
    # Customized Y-axis scaling for each biome
    {
      # Create a data frame to hold Y-axis limits and properties for each biome
      y_scale_df <- data.frame(matrix(NA, nrow=9, ncol=4, byrow=TRUE))
      colnames(y_scale_df) <- c("Panel", "ymin", "ymax", "n")
      
      # Assign biome names and their respective Y-axis limits and break counts
      y_scale_df$Panel <-  c("NH", "TRFO", "WTFO", "TEFO", "BOFO", "SAVA", "STEP", "DESE", "TUND")
      y_scale_df$ymin  <-  c(0.15, rep(0.1, time = 8)) # Set minimum Y-axis values
      y_scale_df$ymax  <-  c(0.25, rep(0.45, time = 8)) # Set maximum Y-axis values
      y_scale_df$n     <-  c(3,    rep(4, time = 8)) # Number of breaks on the Y-axis
      
      # Create a list for individual scale settings for each biome
      y_scale_list <- list(NH   = y_scale_df[1,],
                           TRFO = y_scale_df[2,],
                           WTFO = y_scale_df[3,],
                           TEFO = y_scale_df[4,],
                           BOFO = y_scale_df[5,],
                           SAVA = y_scale_df[6,],
                           STEP = y_scale_df[7,],
                           DESE = y_scale_df[8,],
                           TUND = y_scale_df[9,])
      
      # Apply custom scale settings to each biome using a loop
      y_scale_list_final <- lapply(y_scale_list, function(x) {
        scale_y_continuous(limits = c(x$ymin, x$ymax), n.breaks = x$n, expand = c(0,0))
      })
      
    }
    
    # Plotting the biome compositional turnover
    {
      Biome_compositional_turnover_NH_0_21ka_fig <- ggplot(data=Biome_compositional_turnover_NH_0_21ka, aes(x= Timeslice_older, y = EMD)) +
        geom_line(aes(color = Biome, group = Biome), linetype="solid", linewidth=0.8) +
        facet_grid(Biome ~ ., scales = "free_y", switch = "y",
                   labeller = as_labeller(c("NH" = "Hemisphere", "TRFO" = "Tropical forest", "WTFO" = "Subtropical forest", "TEFO" = "Temperate forest", 
                                            "BOFO" = "Boreal forest", "SAVA" = "Savanna", "STEP" = "Steppe", "DESE" = "Desert", "TUND" = "Tundra"))) +
        ggh4x::facetted_pos_scales(y = y_scale_list_final) +
        scale_x_continuous(limits = c(0.5, 21), breaks = seq(0.5,20.5, by=2), expand = c(0,0),
                           labels= c("0-0.5", "2-2.5", "4-4.5", "6-6.5", "8-8.5", "10-10.5", "12-12.5", "14-14.5", "16-16.5", "18-18.5", "20-20.5")) +
        scale_color_manual(name="", values=c("NH"="black", "TRFO"="#e31a1c", "WTFO"="#f781bf", "TEFO"="#33a02c", "BOFO"="#1f78b4", "SAVA"="#ff7f00", "STEP"="#fdbf6f", "DESE"="#ffed6f", "TUND"="#6a3d9a"),
                           labels = c("NH", "TRFO", "WTFO", "TEFO", "BOFO", "SAVA", "STEP", "DESE", "TUND")) +
        labs(title= "(a) Biome compositional turnover", x="Age (ka cal BP)", y="") +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=13, hjust = 0, face = 'bold')) +
        theme(axis.title=element_text(size=12)) +
        theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt")),
              axis.text.x=element_text(size=10, angle = 45, vjust = 1, hjust = 1)) +
        theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 5, 0, 0), "pt")),
              axis.text.y=element_text(size=10)) +
        theme(legend.position="none",
              legend.text=element_text(size=6),
              legend.title = element_text(size=8, vjust = 0.5, hjust = 0.5),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white"),
              panel.spacing = unit(0.5, "cm")) +
        theme(strip.background=element_blank(),
              strip.text.y=element_text(size =12, angle=0, margin=unit(c(0, 3, 0, 0), "pt"), face = 'bold'),
              strip.placement="outside") +
        theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
      
      Biome_compositional_turnover_NH_0_21ka_fig
      
    }
    
  }
  
  
  # (2)  Proportion of biome shifts in the pollen-based reconstruction
  {
    # Convert 'Biome' into a factor with specified order for plotting purposes
    Biome_shift_proportion_NH_0_21ka$Biome <- factor(Biome_shift_proportion_NH_0_21ka$Biome,
                                                     levels = c("NH", "TRFO", "WTFO", "TEFO", "BOFO", "SAVA", "STEP", "DESE", "TUND"))
    
    # Customized Y-axis scaling for each biome
    {
      # Create a data frame to hold Y-axis limits and properties for each biome
      y_scale_df <- data.frame(matrix(NA, nrow=9, ncol=4, byrow=TRUE))
      colnames(y_scale_df) <- c("Panel", "ymin", "ymax", "n")
      
      # Assign biome names and their respective Y-axis limits and break counts
      y_scale_df$Panel <-  c("NH", "TRFO", "WTFO", "TEFO", "BOFO", "SAVA", "STEP", "DESE", "TUND")
      y_scale_df$ymin  <-  c(0.15, rep(0, time = 8)) # Set minimum Y-axis values
      y_scale_df$ymax  <-  c(0.40, rep(0.7, time = 8)) # Set maximum Y-axis values
      y_scale_df$n     <-  c(3,    rep(5, time = 8)) # Number of breaks on the Y-axis
      
      # Create a list for individual scale settings for each biome
      y_scale_list <- list(NH   = y_scale_df[1,],
                           TRFO = y_scale_df[2,],
                           WTFO = y_scale_df[3,],
                           TEFO = y_scale_df[4,],
                           BOFO = y_scale_df[5,],
                           SAVA = y_scale_df[6,],
                           STEP = y_scale_df[7,],
                           DESE = y_scale_df[8,],
                           TUND = y_scale_df[9,])
      
      # Apply custom scale settings to each biome using a loop
      y_scale_list_final <- lapply(y_scale_list, function(x) {
        scale_y_continuous(limits = c(x$ymin, x$ymax), n.breaks = x$n, expand = c(0,0))
      })
      
    }
    
    {
      Biome_shift_proportion_pollen_NH_0_21ka_fig <- ggplot(data=Biome_shift_proportion_NH_0_21ka, aes(x= Timeslice_older, y = Reconstruction)) +
        geom_line(aes(color = Biome, group = Biome), linetype="solid", linewidth=0.8) +
        facet_grid(Biome ~ ., scales = "free_y", switch = "y") +
        ggh4x::facetted_pos_scales(y = y_scale_list_final) +
        scale_x_continuous(limits = c(0.5, 21), breaks = seq(0.5,20.5, by=2), expand = c(0,0),
                           labels= c("0-0.5", "2-2.5", "4-4.5", "6-6.5", "8-8.5", "10-10.5", "12-12.5", "14-14.5", "16-16.5", "18-18.5", "20-20.5")) +
        scale_color_manual(name="", values=c("NH"="black", "TRFO"="#e31a1c", "WTFO"="#f781bf", "TEFO"="#33a02c", "BOFO"="#1f78b4", "SAVA"="#ff7f00", "STEP"="#fdbf6f", "DESE"="#ffed6f", "TUND"="#6a3d9a"),
                           labels = c("NH", "TRFO", "WTFO", "TEFO", "BOFO", "SAVA", "STEP", "DESE", "TUND")) +
        labs(title= "(b) Proportion of biome shifts in reconstruction", x="Age (ka cal BP)", y="") +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=13, hjust = 0, face = 'bold')) +
        theme(axis.title=element_text(size=12)) +
        theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt")),
              axis.text.x=element_text(size=10, angle = 45, vjust = 1, hjust = 1)) +
        theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 5, 0, 0), "pt")),
              axis.text.y=element_text(size=10)) + 
        theme(legend.position="none",
              legend.text=element_text(size=6),
              legend.title = element_text(size=8, vjust = 0.5, hjust = 0.5),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white"),
              panel.spacing = unit(0.5, "cm")) +
        theme(strip.background=element_blank(),
              strip.text.y=element_blank(),
              strip.placement="outside") +
        theme(plot.margin=unit(c(0.1,0.25,0.1,0.5),"cm")) 
      
      Biome_shift_proportion_pollen_NH_0_21ka_fig 
      
    }
    
  }
  
  
  # (3)  Proportion of biome shifts in the ESM-based simulation
  {
    # Customized Y-axis scaling for each biome
    {
      ## set the result matrix:
      y_scale_df <- data.frame(matrix(NA, nrow=9, ncol=4, byrow=TRUE))
      colnames(y_scale_df) <- c("Panel", "ymin", "ymax", "n")
      
      # Assign biome names and their respective Y-axis limits and break counts
      y_scale_df$Panel <-  c("NH", "TRFO", "WTFO", "TEFO", "BOFO", "SAVA", "STEP", "DESE", "TUND")
      y_scale_df$ymin  <-  c(0, rep(0, time = 8)) # Set minimum Y-axis values
      y_scale_df$ymax  <-  c(0.40, rep(0.605, time = 8)) # Set maximum Y-axis values
      y_scale_df$n     <-  c(3,    rep(4, time = 8)) # Number of breaks on the Y-axis
      
      # Create a list for individual scale settings for each biome
      y_scale_list <- list(NH   = y_scale_df[1,],
                           TRFO = y_scale_df[2,],
                           WTFO = y_scale_df[3,],
                           TEFO = y_scale_df[4,],
                           BOFO = y_scale_df[5,],
                           SAVA = y_scale_df[6,],
                           STEP = y_scale_df[7,],
                           DESE = y_scale_df[8,],
                           TUND = y_scale_df[9,])
      
      # Apply custom scale settings to each biome using a loop
      y_scale_list_final <- lapply(y_scale_list, function(x) {
        scale_y_continuous(limits = c(x$ymin, x$ymax), n.breaks = x$n, expand = c(0,0))
      })
      
    }
    
    {
      Biome_shift_proportion_MPI_NH_0_21ka_fig <- ggplot(data=Biome_shift_proportion_NH_0_21ka, aes(x= Timeslice_older, y = Simulation)) +
        geom_line(aes(color = Biome, group = Biome), linetype="solid", linewidth=0.8) +
        facet_grid(Biome ~ ., scales = "free_y", switch = "y") +
        ggh4x::facetted_pos_scales(y = y_scale_list_final) +
        scale_x_continuous(limits = c(0.5, 21), breaks = seq(0.5,20.5, by=2), expand = c(0,0),
                           labels= c("0-0.5", "2-2.5", "4-4.5", "6-6.5", "8-8.5", "10-10.5", "12-12.5", "14-14.5", "16-16.5", "18-18.5", "20-20.5")) +
        scale_color_manual(name="", values=c("NH"="black", "TRFO"="#e31a1c", "WTFO"="#f781bf", "TEFO"="#33a02c", "BOFO"="#1f78b4", "SAVA"="#ff7f00", "STEP"="#fdbf6f", "DESE"="#ffed6f", "TUND"="#6a3d9a"),
                           labels = c("NH", "TRFO", "WTFO", "TEFO", "BOFO", "SAVA", "STEP", "DESE", "TUND")) +
        labs(title= "(c) Proportion of biome shifts in simulation", x="Age (ka cal BP)", y="") +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=13, hjust = 0, face = 'bold')) +
        theme(axis.title=element_text(size=12)) +
        theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt")),
              axis.text.x=element_text(size=10, angle = 45, vjust = 1, hjust = 1)) +
        theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 5, 0, 0), "pt")),
              axis.text.y=element_text(size=10)) + 
        theme(legend.position="none",
              legend.text=element_text(size=6),
              legend.title = element_text(size=8, vjust = 0.5, hjust = 0.5),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white"),
              panel.spacing = unit(0.5, "cm")) +
        theme(strip.background=element_blank(),
              strip.text.y=element_blank(),
              strip.placement="outside") +
        theme(plot.margin=unit(c(0.1,0.25,0.1,0.25),"cm")) 
      
      Biome_shift_proportion_MPI_NH_0_21ka_fig
      
    }
    
  }
  
  
  # (4)  State of composite bioclimatic variables
  {
    # Convert 'Source' into a factor with specific levels and custom mathematical expression-based labels
    Bioclimate_state_NH_0_21ka$Source <- factor(Bioclimate_state_NH_0_21ka$Source,
                                                levels = c("Seasonality state", "Tnon-growing state", "Tgrowing state", "Pgrowing state", "Pnon-growing state"),
                                                labels = c(expression(bold(Seasonality)), expression(bold(T)[non-growing]), expression(bold(T)[growing]), expression(bold(P)[non-growing]), expression(bold(P)[growing])))
    
    # Convert 'Timeslice_younger' into a categorical factor with discrete levels for plotting purposes
    Bioclimate_state_NH_0_21ka$Timeslice_younger <- factor(Bioclimate_state_NH_0_21ka$Timeslice_younger,
                                                           levels = c("0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "6.5", "7", "7.5", "8", "8.5", "9", "9.5", "10", "10.5", "11", "11.5", 
                                                                      "12", "12.5", "13", "13.5", "14", "14.5", "15", "15.5", "16", "16.5", "17", "17.5", "18", "18.5", "19", "19.5", "20", "20.5"))
    
    Bioclimate_state_NH_0_21ka_fig <- ggplot(data = Bioclimate_state_NH_0_21ka, aes(x = Timeslice_younger, y = Value)) +
      geom_boxplot(aes(fill = Source), width = 0.8, position = position_dodge(width = 0.9), outlier.shape = NA, color = "grey70", alpha = 0.6, fatten = NULL) +
      stat_summary(fun = median, geom = "line", aes(color = Source, group = Source), position = position_dodge(width = 0.9), color = "black", linewidth = 0.8) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5), expand = c(0, 0)) +
      scale_x_discrete(breaks = seq(0, 20.5, by = 2), labels = c("0-0.5", "2-2.5", "4-4.5", "6-6.5", "8-8.5", "10-10.5", "12-12.5", "14-14.5", "16-16.5", "18-18.5", "20-20.5")) +
      scale_fill_manual(name = "", values = c("#a65628", "#fb9a99", "#e31a1c", "#1f78b4", "#a6cee3")) +
      facet_grid(rows = vars(Source), scales = "free_y", switch = "y", labeller = label_parsed) +  # Stack panels vertically
      labs(x="Age (ka cal BP)", y="", title = "(d) State of composite bioclimatic variables") +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=13, hjust = 0, face = 'bold')) +
      theme(axis.title=element_text(size=12)) +
      theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt")),
            axis.text.x=element_text(size=10, angle = 45, vjust = 1, hjust = 1)) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 5, 0, 0), "pt")),
            axis.text.y=element_text(size=10)) +  
      theme(legend.position="none",
            legend.text=element_text(size=6),
            legend.title = element_text(size=8, vjust = 0.5, hjust = 0.5),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white"),
            panel.spacing = unit(0.5, "cm")) +
      theme(strip.background=element_blank(),
            strip.text.y=element_text(size =12, angle=0, margin=unit(c(0, 3, 0, 0), "pt"), face = 'bold'),
            strip.placement="outside") +
      theme(plot.margin=unit(c(0.1,0.5,0.25,0),"cm")) 
    
    Bioclimate_state_NH_0_21ka_fig
    
  }
  
  
  # (5)  Rate of change of composite bioclimatic variables
  {
    # Convert 'Source' into a factor with specific levels and custom mathematical expression-based labels
    Bioclimate_change_rate_NH_0_21ka$Source <- factor(Bioclimate_change_rate_NH_0_21ka$Source,
                                                      levels = c("Seasonality change", "Tnon-growing change", "Tgrowing change", "Pgrowing change", "Pnon-growing change"),
                                                      labels = c(expression(bold(Seasonality)), expression(bold(T)[non-growing]), expression(bold(T)[growing]), expression(bold(P)[non-growing]), expression(bold(P)[growing])))
    
    # Convert 'Timeslice_younger' into a categorical factor with discrete levels for plotting purposes
    Bioclimate_change_rate_NH_0_21ka$Timeslice_younger <- factor(Bioclimate_change_rate_NH_0_21ka$Timeslice_younger,
                                                                 levels = c("0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "6.5", "7", "7.5", "8", "8.5", "9", "9.5", "10", "10.5", "11", "11.5", 
                                                                            "12", "12.5", "13", "13.5", "14", "14.5", "15", "15.5", "16", "16.5", "17", "17.5", "18", "18.5", "19", "19.5", "20", "20.5"))
    
    Bioclimate_change_rate_NH_0_21ka_fig <- ggplot(data = Bioclimate_change_rate_NH_0_21ka, aes(x = Timeslice_younger, y = Value)) +
      geom_boxplot(aes(fill = Source), width = 0.8, position = position_dodge(width = 0.9), outlier.shape = NA, color = "grey70", alpha = 0.6, fatten = NULL) +
      stat_summary(fun = median, geom = "line", aes(color = Source, group = Source), position = position_dodge(width = 0.9), color = "black", linewidth = 0.8) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5), expand = c(0, 0)) +
      scale_x_discrete(breaks = seq(0, 20.5, by = 2), labels = c("0-0.5", "2-2.5", "4-4.5", "6-6.5", "8-8.5", "10-10.5", "12-12.5", "14-14.5", "16-16.5", "18-18.5", "20-20.5")) +
      scale_fill_manual(name = "", values = c("#a65628", "#fb9a99", "#e31a1c", "#1f78b4", "#a6cee3")) +
      facet_grid(rows = vars(Source), scales = "free_y", switch = "y", labeller = label_parsed) +  # Stack panels vertically
      labs(x="Age (ka cal BP)", y="", title = "(e) Change rate of composite bioclimatic variables") +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=13, hjust = 0, face = 'bold')) +
      theme(axis.title=element_text(size=12)) +
      theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt")),
            axis.text.x=element_text(size=10, angle = 45, vjust = 1, hjust = 1)) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 5, 0, 0), "pt")),
            axis.text.y=element_text(size=10)) +  
      theme(legend.position="none",
            legend.text=element_text(size=6),
            legend.title = element_text(size=8, vjust = 0.5, hjust = 0.5),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white"),
            panel.spacing = unit(0.5, "cm")) +
      theme(strip.background=element_blank(),
            strip.text.y=element_text(size =12, angle=0, margin=unit(c(0, 3, 0, 0), "pt"), face = 'bold'),
            strip.placement="outside") +
      theme(plot.margin=unit(c(0.25,0.5,0.1,0),"cm")) 
    
    Bioclimate_change_rate_NH_0_21ka_fig
    
  }
  
  
  # -------------------------------------------------------------------------------------------------
  # Arranging plots
  Biome_Bioclimate_NH_0_21ka_fig  <- ggarrange(Biome_compositional_turnover_NH_0_21ka_fig, Biome_shift_proportion_pollen_NH_0_21ka_fig, Biome_shift_proportion_MPI_NH_0_21ka_fig,
                                               ggarrange(Bioclimate_state_NH_0_21ka_fig, Bioclimate_change_rate_NH_0_21ka_fig,
                                                         nrow = 2, ncol = 1, common.legend = F, legend = "none"),
                                               nrow = 1, ncol = 4, common.legend = F, legend = "none")
  
  # Save the final figure as a PNG file
  ggsave("result/Fig. 2-Temporal patterns of biome and bioclimate change in the Northern Hemisphere over the last 21,000 years.png", width = 52, height = 40, units = "cm", dpi = 300)
  
}



# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 3. Fig. 3: Gaussian kernel correlations between changes in biomes and composite bioclimatic variables  ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
GCK_biome_climate_state         <- read.csv2("data/Fig.3 data.1-Gaussian kernel correlations between changes in biomes and state of composite bioclimatic variables.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
GCK_biome_climate_change_rate   <- read.csv2("data/Fig.3 data.2-Gaussian kernel correlations between changes in biomes and change rate of composite bioclimatic variables.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
GCK_biome_climate_state         <- type.convert(GCK_biome_climate_state, as.is = TRUE) 
GCK_biome_climate_change_rate   <- type.convert(GCK_biome_climate_change_rate, as.is = TRUE) 

# Plot
{
  ## Gaussian kernel correlations between changes in biomes and state of composite bioclimatic variables
  {
    # Group by 'Source' and calculate the maximum absolute R-value for each group
    GCK_biome_climate_state_new <- GCK_biome_climate_state %>%
      group_by(Source) %>%
      mutate(R_value_Max = max(abs(R_value))) %>%
      ungroup()
    
    # Convert 'Source' to a factor with custom levels for plot ordering
    GCK_biome_climate_state_new$Source <- factor(GCK_biome_climate_state_new$Source,
                                                 levels = rev(c("Biome compositional turnover", "Proportion of biome shifts in reconstruction", "Proportion of biome shifts in simulation")))
    
    # Convert 'Variable' to a factor with custom levels for plot ordering
    GCK_biome_climate_state_new$Variable <- factor(GCK_biome_climate_state_new$Variable,
                                                   levels = c("Seasonality_mean", "T_growing_mean", "T_non_growing_mean", "P_growing_mean", "P_non_growing_mean"))
    
    {
      GCK_biome_climate_state_fig <- ggplot(data = GCK_biome_climate_state_new, aes(x = Variable, y = Source)) +
        geom_tile(aes(fill = R_value), alpha = 0.6) +  
        geom_text(aes(label = ifelse(P_value == "<0.001", paste0(sprintf("%.2f", R_value), " ***"),
                                     ifelse(P_value == "0.001-0.01", paste0(sprintf("%.2f", R_value), " **"),
                                            ifelse(P_value == "0.01-0.05", paste0(sprintf("%.2f", R_value), " *"), sprintf("%.2f", R_value)))),
                      vjust = 1, fontface = ifelse(abs(R_value) == R_value_Max, "bold", "plain")), size = 6)  +
        scale_x_discrete(labels = c("Seasonality_mean" = "Seasonality", "T_growing_mean" = expression(T[growth]), "T_non_growing_mean" = expression(T[non-growth]), 
                                    "P_growing_mean" = expression(P[growth]), "P_non_growing_mean" = expression(P[non-growth]))) +
        scale_fill_gradient2(name = "R-value", low = "blue", high = "red",
                             limits = c(range(GCK_biome_climate_state_new$R_value)[1], range(GCK_biome_climate_state_new$R_value)[2]),
                             breaks = seq(-0.6, 0.8, 0.2), labels = c("-0.6", "-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8")) +
        labs(title= '(a) Composite bioclimatic variable states', x="", y="") +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=18, hjust = 0, face = 'bold')) +
        theme(axis.title=element_text(size=12)) +
        theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt"), size=12),
              axis.text.x=element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold', 
                                       color=c("#a65628", "#e31a1c", "#fb9a99", "#1f78b4", "#a6cee3"))) +
        theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
              axis.text.y=element_text(size=16, vjust = 0.5, hjust = 0))  +  
        theme(legend.position="none",
              legend.text=element_text(size=12),
              legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(strip.background=element_blank(),
              strip.text=element_text(size=10, face="bold"),
              strip.text.y=element_text(angle=0, margin=unit(c(0, 10, 0, 0), "pt")),
              strip.placement="outside") +
        theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))  +
        guides(fill = guide_colorbar(barwidth = unit(0.6, "cm"), barheight = unit(20, "cm")))
      
      GCK_biome_climate_state_fig
      
    }
    
  }
  
  
  
  ## Gaussian kernel correlations between changes in biomes and change rate of composite bioclimatic variables
  {
    # Group by 'Source' and calculate the maximum absolute R-value for each group
    GCK_biome_climate_change_rate_new <- GCK_biome_climate_change_rate %>%
      group_by(Source) %>%
      mutate(R_value_Max = max(abs(R_value))) %>%
      ungroup()
    
    # Convert 'Source' to a factor with custom levels for plot ordering
    GCK_biome_climate_change_rate_new$Source <- factor(GCK_biome_climate_change_rate_new$Source,
                                                       levels = rev(c("Biome compositional turnover", "Proportion of biome shifts in reconstruction", "Proportion of biome shifts in simulation")))
    
    # Convert 'Variable' to a factor with custom levels for plot ordering
    GCK_biome_climate_change_rate_new$Variable <- factor(GCK_biome_climate_change_rate_new$Variable,
                                                         levels = c("Seasonality_mean", "T_growing_mean", "T_non_growing_mean", "P_growing_mean", "P_non_growing_mean"))
    
    {
      GCK_biome_climate_change_rate_fig <- ggplot(data = GCK_biome_climate_change_rate_new, aes(x = Variable, y = Source)) +
        geom_tile(aes(fill = R_value), alpha = 0.6) +  
        geom_text(aes(label = ifelse(P_value == "<0.001", paste0(sprintf("%.2f", R_value), " ***"),
                                     ifelse(P_value == "0.001-0.01", paste0(sprintf("%.2f", R_value), " **"),
                                            ifelse(P_value == "0.01-0.05", paste0(sprintf("%.2f", R_value), " *"), sprintf("%.2f", R_value)))),
                      vjust = 1, fontface = ifelse(abs(R_value) == R_value_Max, "bold", "plain")), size = 6)  +
        scale_x_discrete(labels = c("Seasonality_mean" = "Seasonality", "T_growing_mean" = expression(T[growth]), "T_non_growing_mean" = expression(T[non-growth]), 
                                    "P_growing_mean" = expression(P[growth]), "P_non_growing_mean" = expression(P[non-growth]))) +
        scale_fill_gradient2(name = "R-value", low = "blue", high = "red",
                             limits = c(range(GCK_biome_climate_change_rate_new$R_value)[1], range(GCK_biome_climate_change_rate_new$R_value)[2]),
                             breaks = seq(-0.6, 0.8, 0.2), labels = c("-0.6", "-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8")) +
        labs(title= '(b) Composite bioclimatic variable rates of change', x="", y="") +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=18, hjust = 0, face = 'bold')) +
        theme(axis.title=element_text(size=12)) +
        theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt"), size=12),
              axis.text.x=element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold', 
                                       color=c("#a65628", "#e31a1c", "#fb9a99", "#1f78b4", "#a6cee3"))) +
        theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
              axis.text.y=element_text(size=16, vjust = 0.5, hjust = 0))  +  
        theme(legend.position="none",
              legend.text=element_text(size=12),
              legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(strip.background=element_blank(),
              strip.text=element_text(size=10, face="bold"),
              strip.text.y=element_text(angle=0, margin=unit(c(0, 10, 0, 0), "pt")),
              strip.placement="outside") +
        theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))  +
        guides(fill = guide_colorbar(barwidth = unit(0.6, "cm"), barheight = unit(20, "cm")))
      
      GCK_biome_climate_change_rate_fig
      
    }
    
  }
  
  # -------------------------------------------------------------------------------------------------
  # Arranging plots
  GCK_biome_climate_fig  <- ggarrange(GCK_biome_climate_state_fig, GCK_biome_climate_change_rate_fig,
                                      nrow = 2, ncol = 1, common.legend = F)
  
  # Save the final figure as a PNG file
  ggsave("result/Fig. 3-Gaussian kernel correlations between changes in biomes and composite bioclimatic variables.png", width = 15, height = 13, units = "in", dpi = 300)
  
}



# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 4. Fig. 4: Potential spatial patterns of bioclimates and biome changes in the MPI-ESM simulation over the next 300 years in the Northern Hemisphere  ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
Bioclimatic_variable_state         <- read.csv2("data/Fig.4 data.1-Bioclimatic variable states over the next 300 years in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Bioclimatic_variable_change_rate   <- read.csv2("data/Fig.4 data.2-Bioclimatic variable rates of change over the next 300 years in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Biome_shift_magnitude_2250CE       <- read.csv2("data/Fig.4 data.3-Magnitude of biome shifts at timeslice 2250 C.E. in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Bioclimatic_variable_state        <- type.convert(Bioclimatic_variable_state, as.is = TRUE) 
Bioclimatic_variable_change_rate  <- type.convert(Bioclimatic_variable_change_rate, as.is = TRUE) 
Biome_shift_magnitude_2250CE      <- type.convert(Biome_shift_magnitude_2250CE, as.is = TRUE) 

# Plot
{
  # (1) The magnitude of both seasonality and growing season temperature rate of change
  {
    # Standardization
    {
      # Custom min-max normalization function
      min_max_normalize <- function(x) {
        return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
      }
      
      # Apply min-max normalization to specified columns
      Bioclimatic_variable_state_standardization <- Bioclimatic_variable_state %>%
        mutate(across(4:17, ~ min_max_normalize(.)))
      
      Bioclimatic_variable_change_rate_standardization <- Bioclimatic_variable_change_rate %>%
        mutate(across(5:18, ~ min_max_normalize(.)))
    }
    
    # Creat datasets
    {
      # Calculate datasets
      Bioclimatic_variable_state_standardization <- Bioclimatic_variable_state_standardization[ ,c(1:3,5,13)] %>%
        mutate(Seasonality_state = rowMeans(dplyr::select(., 4:5)))
      
      Bioclimatic_variable_change_rate_standardization <- Bioclimatic_variable_change_rate_standardization[ ,c(1:3,7,9)] %>%
        mutate(Tgrowing_rate = rowMeans(dplyr::select(., 4:5)))
      names(Bioclimatic_variable_change_rate_standardization)[3] <- "Timeslice"
      
      # combine datasets
      Bioclimatic_variable_state_change_rate_standardization <- full_join(Bioclimatic_variable_state_standardization[ ,c(1:3,6)], Bioclimatic_variable_change_rate_standardization[ ,c(1:3,6)], 
                                                                          by = c("Longitude", "Latitude", "Timeslice"))
    }
    
    # Create class
    {
      # Calculate quantiles for Seasonality state and Tgrowing_rate
      Seasonality_state_quantiles <- quantile(Bioclimatic_variable_state_change_rate_standardization$Seasonality_state, probs = seq(0, 1, by = 0.25))
      Tgrowing_rate_quantiles     <- quantile(Bioclimatic_variable_state_change_rate_standardization$Tgrowing_rate, probs = seq(0, 1, by = 0.25))
      
      # Categorize data based on quantiles
      Bioclimatic_variable_state_change_rate_standardization <- Bioclimatic_variable_state_change_rate_standardization %>%
        mutate(
          Seasonality_state_level = cut(Seasonality_state, breaks = Seasonality_state_quantiles, labels = FALSE, include.lowest = TRUE),
          Tgrowing_rate_level     = cut(Tgrowing_rate, breaks = Tgrowing_rate_quantiles, labels = FALSE, include.lowest = TRUE)
        ) %>%
        mutate(
          Class = (Seasonality_state_level - 1) * 4 + Tgrowing_rate_level
        ) %>%
        dplyr::select(everything(), Class)
      
    }
    
    # -------------------------------------------------------------------------------------------------
    # Plot
    {
      # Define gradient colors
      {
        color_gradient <- c("#ffffd4", "#d9f0a3", "#addd8e", "#008000", 
                            "#F7E5B0", "#FAD6A5", "#F2E1C1", "#B8D8D8",  
                            "#D7B272", "#FFD39B", "#FFB28C", "#FF6F6F", 
                            "#8B4513", "#E4B97F", "#FF7F7F", "#FF0000")
        
        # Map colors to a data frame
        color_mapping <- setNames(color_gradient, as.character(1:length(color_gradient)))
        
        # Retrieve coastline data
        coastline <- ne_coastline(scale = "medium", returnclass = "sf")
      }
      
      ## Next 100 years (-0.1 ka BP)
      {
        Bioclimatic_variable_state_change_rate_next100yr_NH_map <- ggplot() +
          geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
          geom_tile(data= subset(Bioclimatic_variable_state_change_rate_standardization, Timeslice == -0.1), aes(x=Longitude, y=Latitude, fill = factor(Class)), alpha = 0.6) +
          scale_fill_manual(name = "",values=color_mapping ) +
          labs(title= '(a) Magnitude of both seasonality and growing season temperature rate of change at 2050 C.E.', x="", y="") +
          coord_sf(ylim = c(0, 90), expand = FALSE) +
          scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
          theme_bw() +  # Reset the plot background
          theme(plot.title=element_text(size=14, hjust = 0, face = 'bold')) +
          theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt")),
                axis.text.x=element_text(size=12)) +
          theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
                axis.text.y=element_text(size=12)) +
          theme(legend.position="none",
                legend.text=element_text(size=14),
                legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
                legend.key.size = unit(1.2, "lines"),
                legend.key.height = unit(1, "lines"),
                legend.key.width = unit(2, "lines"),
                legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
                legend.box = "vertical",
                legend.justification = c(0.5, 0.5)) +
          theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.border=element_rect(linetype="solid", fill=NA),
                panel.background=element_rect(fill="white")) +
          theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))  
        
        Bioclimatic_variable_state_change_rate_next100yr_NH_map 
      }
      
      ## Next 200 years (-0.2 ka BP)
      {
        Bioclimatic_variable_state_change_rate_next200yr_NH_map <- ggplot() +
          geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
          geom_tile(data= subset(Bioclimatic_variable_state_change_rate_standardization, Timeslice == -0.2), aes(x=Longitude, y=Latitude, fill = factor(Class)), alpha = 0.6) +
          scale_fill_manual(name = "",values=color_mapping ) +
          labs(title= '(b) Magnitude of both seasonality and growing season temperature rate of change at 2150 C.E.', x="", y="") +
          coord_sf(ylim = c(0, 90), expand = FALSE) +
          scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
          theme_bw() +  # Reset the plot background
          theme(plot.title=element_text(size=14, hjust = 0, face = 'bold')) +
          theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt")),
                axis.text.x=element_text(size=12)) +
          theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
                axis.text.y=element_text(size=12)) +
          theme(legend.position="none",
                legend.text=element_text(size=14),
                legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
                legend.key.size = unit(1.2, "lines"),
                legend.key.height = unit(1, "lines"),
                legend.key.width = unit(2, "lines"),
                legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
                legend.box = "vertical",
                legend.justification = c(0.5, 0.5)) +
          theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.border=element_rect(linetype="solid", fill=NA),
                panel.background=element_rect(fill="white")) +
          theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))  
        
        Bioclimatic_variable_state_change_rate_next200yr_NH_map 
      }
      
      ## Next 300 years (-0.3 ka BP)
      {
        Bioclimatic_variable_state_change_rate_next300yr_NH_map <- ggplot() +
          geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
          geom_tile(data= subset(Bioclimatic_variable_state_change_rate_standardization, Timeslice == -0.3), aes(x=Longitude, y=Latitude, fill = factor(Class)), alpha = 0.6) +
          scale_fill_manual(name = "",values=color_mapping ) +
          labs(title= '(c) Magnitude of both seasonality and growing season temperature rate of change at 2250 C.E.', x="", y="") +
          coord_sf(ylim = c(0, 90), expand = FALSE) +
          scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
          theme_bw() +  # Reset the plot background
          theme(plot.title=element_text(size=14, hjust = 0, face = 'bold')) +
          theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt")),
                axis.text.x=element_text(size=12)) +
          theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
                axis.text.y=element_text(size=12)) +
          theme(legend.position="none",
                legend.text=element_text(size=14),
                legend.title = element_text(size=16, vjust = 0.5, hjust = 0.5, face = 'bold'),
                legend.key.size = unit(1.2, "lines"),
                legend.key.height = unit(1, "lines"),
                legend.key.width = unit(2, "lines"),
                legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
                legend.box = "vertical",
                legend.justification = c(0.5, 0.5)) +
          theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.border=element_rect(linetype="solid", fill=NA),
                panel.background=element_rect(fill="white")) +
          theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm"))  
        
        Bioclimatic_variable_state_change_rate_next300yr_NH_map 
      }
      
      ## Legend
      {
        # Define fill colors
        color_gradient <- c("#ffffd4", "#d9f0a3", "#addd8e", "#008000", 
                            "#F7E5B0", "#FAD6A5", "#F2E1C1", "#B8D8D8",  
                            "#D7B272", "#FFD39B", "#FFB28C", "#FF6F6F", 
                            "#8B4513", "#E4B97F", "#FF7F7F", "#FF0000")
        
        color_gradient_df <- data.frame(
          fill = color_gradient,
          id = 1:length(color_gradient))
        
        # Create a data frame with coordinates of four vertices for 16 small squares and their fill colors
        data <- expand.grid(
          x_start = seq(0, 0.75, by = 0.25),  # Starting x-coordinates
          y_start = seq(0, 0.75, by = 0.25)   # Starting y-coordinates
        ) %>%
          rowwise() %>%
          mutate(
            # Coordinates of the four vertices of each square
            x = list(c(x_start, x_start + 0.25, x_start + 0.25, x_start)),
            y = list(c(y_start, y_start, y_start + 0.25, y_start + 0.25))
          ) %>%
          unnest(c(x, y)) %>%
          group_by(x_start, y_start) %>%
          mutate(id = cur_group_id()) %>%
          ungroup()
        
        # Merge color data into the data frame
        data <- data %>%
          left_join(color_gradient_df, by = "id")
        
        # Calculate the center coordinates of each square
        centers <- data %>%
          group_by(id, fill) %>%
          summarize(
            x_center = mean(x),
            y_center = mean(y),
            .groups = 'drop'
          )
        
        # Plot squares and add text labels
        color_legend <- ggplot(data, aes(x = x, y = y, group = id, fill = fill)) +
          geom_polygon(alpha = 0.6) +
          scale_fill_identity() + 
          labs(title= '', x = "Seasonality state", y = expression(bold(T[growing] ~ "change rate"))) +
          geom_text(data = centers, aes(x = x_center, y = y_center, label = id), 
                    color = "black", size = 5, vjust = 0.5, hjust = 0.5) +  # Add IDs at the center
          scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("Weak          ", "50%", "Strong"), expand = c(0, 0)) +
          scale_y_continuous(limits = c(0, 1), breaks = c(0.5, 1), labels = c("50%", "Strong"), expand = c(0, 0)) +
          theme_bw() +  # Reset the plot background
          theme(plot.title = element_text(size = 14)) +
          theme(axis.title.x = element_text(margin = unit(c(10, 0, 0, 0), "pt"), face = 'bold', size = 11),
                axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 0.5, face = 'bold')) +
          theme(axis.title.y = element_text(angle = 90, margin = unit(c(0, 0, 10, 0), "pt"), face = 'bold', size = 11),
                axis.text.y = element_text(size = 10, vjust = 0.5, hjust = 0.5, face = 'bold')) +
          theme(axis.ticks.length = unit(0, "pt")) +
          theme(legend.position = "none",
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 16, vjust = 0.5, hjust = 0.5, face = 'bold'),
                legend.key.size = unit(1.2, "lines"),
                legend.key.height = unit(1, "lines"),
                legend.key.width = unit(2, "lines"),
                legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid", colour = "transparent"),
                legend.box = "vertical",
                legend.justification = c(0.5, 0.5)) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_rect(linetype = "solid", fill = NA),
                panel.background = element_rect(fill = "white")) +
          theme(plot.margin = unit(c(0.25, 0.75, 0.25, 0.25), "cm")) 
        
        color_legend  
        
      }
      
    }
    
  }
  
  
  # (2) Magnitude of biome shifts at timeslice 2250 C.E.
  {
    # Convert the 'Magnitude' column to a factor with specified levels
    Biome_shift_magnitude_2250CE$Magnitude <- factor(Biome_shift_magnitude_2250CE$Magnitude,
                                                     levels = c("0", "1", "2", "3")) 
    
    {
      Biome_shift_magnitude_2250CE_NH_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey90", size = 0.2) +
        geom_tile(data= Biome_shift_magnitude_2250CE, aes(x=Longitude, y=Latitude, fill = Magnitude), alpha = 0.6) +
        scale_fill_manual(name   = "Magnitude", values=c("0" = "grey70", "1" = "#ffed6f", "2" = "#D95F0E", "3" = "#FF0000")) +
        labs(title= '(d) Magnitude of biome shifts at 2250 C.E.', x="", y="") +
        coord_sf(ylim = c(0, 90), expand = FALSE) +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=14, hjust = 0, face = 'bold')) +
        theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt")),
              axis.text.x=element_text(size=12)) +
        theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
              axis.text.y=element_text(size=12)) +
        theme(legend.position="right",
              legend.text=element_text(size=12),
              legend.title = element_text(size=12, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
        guides(fill=guide_legend(override.aes=list(), nrow=4, byrow=TRUE))
      
      Biome_shift_magnitude_2250CE_NH_map 
    }
    
  }
  
  
  # (3) Magnitude of biome shifts at timeslice 2250 C.E.
  {
    # Merge datasets
    Bioclimatic_variable_state_change_rate_Biome_shift_magnitude_2250CE <- full_join(
      subset(Bioclimatic_variable_state_change_rate_standardization, Timeslice == -0.3), Biome_shift_magnitude_2250CE, by = c("Longitude", "Latitude", "Timeslice")) %>%
      na.omit() %>% # Remove any rows with missing values
      mutate(
        # Create a new 'Level' column based on the 'Class' values
        Level = case_when(
          Class %in% c(1, 2, 5, 6)     ~ 1,
          Class %in% c(3, 4, 7, 8)     ~ 2,
          Class %in% c(9, 10, 13, 14)  ~ 3,
          Class %in% c(11, 12, 15, 16) ~ 4))
    
    # Calculation by class and level
    Bioclimatic_variable_state_change_rate_Biome_shift_magnitude_2250CE_summary <- Bioclimatic_variable_state_change_rate_Biome_shift_magnitude_2250CE %>%
      group_by(Class) %>%
      summarise(
        Biome_shift_proportion = mean(Magnitude != "0", na.rm = TRUE),
        Magnitude_mean =  mean(as.numeric(as.character(Magnitude)), na.rm = TRUE),
        Level = first(Level)
      ) %>%
      # Further summarization by 'Level
      group_by(Level) %>%
      summarise(
        Biome_shift_proportion_mean  = mean(Biome_shift_proportion, na.rm = TRUE),
        Magnitude_mean_mean          = mean(Magnitude_mean, na.rm = TRUE)
      )
    
    # Pivot the data into a long format for plotting
    Bioclimatic_variable_state_change_rate_Biome_shift_magnitude_2250CE_summary_long <- Bioclimatic_variable_state_change_rate_Biome_shift_magnitude_2250CE_summary  %>%
      pivot_longer(cols = c(Biome_shift_proportion_mean, Magnitude_mean_mean), 
                   names_to = "Variable", 
                   values_to = "Value")
    
    # Convert 'Level' and 'Variable' to factors with custom levels
    Bioclimatic_variable_state_change_rate_Biome_shift_magnitude_2250CE_summary_long$Level <- factor(Bioclimatic_variable_state_change_rate_Biome_shift_magnitude_2250CE_summary_long$Level,
                                                                                                     levels = c("1", "2", "3", "4"))
    Bioclimatic_variable_state_change_rate_Biome_shift_magnitude_2250CE_summary_long$Variable <- factor(Bioclimatic_variable_state_change_rate_Biome_shift_magnitude_2250CE_summary_long$Variable,
                                                                                                        levels = c("Biome_shift_proportion_mean", "Magnitude_mean_mean"))
    
    Bioclimatic_variable_state_change_rate_Biome_shift_magnitude_2250CE_level_barplot <- ggplot(data=Bioclimatic_variable_state_change_rate_Biome_shift_magnitude_2250CE_summary_long, aes(x = Level, y = Value, fill = Variable)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.6) +
      scale_fill_manual(name = "", values = c("Biome_shift_proportion_mean" = "blue", "Magnitude_mean_mean" = "#FF0000"), 
                        labels = c("Biome_shift_proportion_mean" = "Proportion of biome shifts", "Magnitude_mean_mean" = "Magnitude of biome shifts")) +
      labs(title= '(e) Proportion of biome shifts and magnitude at 2250 C.E.', x="", y = "Value") +
      scale_y_continuous(limits = c(0,2.5), breaks = seq(0,2,by=0.5), labels = c("0.0", "0.5", "1.0", "1.5", "2.0"), expand = c(0,0)) +
      scale_x_discrete(labels = c(paste0("Weak seasonality state", "\n", "+", "\n", "Weak", " Tgrowing ", "change rate"),
                                  paste0("Weak seasonality state", "\n", "+", "\n", "Strong", " Tgrowing ", "change rate"),
                                  paste0("Strong seasonality state", "\n", "+", "\n", "Weak", " Tgrowing ", "change rate"),
                                  paste0("Strong seasonality state", "\n", "+", "\n", "Strong", " Tgrowing ", "change rate"))) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=14, hjust = 0, face = 'bold')) +
      theme(axis.title.x=element_text(size=12, margin=unit(c(5, 0, 0, 0), "pt"), face = 'bold'),
            axis.text.x=element_text(size=12)) +
      theme(axis.title.y=element_text(size=12, angle=90, margin=unit(c(0, 5, 0, 0), "pt"), face = 'bold'),
            axis.text.y=element_text(size=12)) +
      theme(legend.position="bottom",
            legend.text=element_text(size=12),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
      guides(fill=guide_legend(override.aes=list(), nrow=1, byrow=TRUE))
    
    Bioclimatic_variable_state_change_rate_Biome_shift_magnitude_2250CE_level_barplot
    
  }
  
  
  # -------------------------------------------------------------------------------------------------
  
  # Arrange plots
  blank_plot <- ggplot() + theme_void()
  
  Simulated_climate_biome_2250CE_NH_map <- ggarrange(ggarrange(ggarrange(ggarrange(Bioclimatic_variable_state_change_rate_next100yr_NH_map, Bioclimatic_variable_state_change_rate_next200yr_NH_map, 
                                                                                   Bioclimatic_variable_state_change_rate_next300yr_NH_map, ncol = 1, nrow = 3, common.legend = FALSE),
                                                                         ggarrange(blank_plot, blank_plot, color_legend, ncol = 1, nrow = 3, heights = c(0.425,0.425,0.15), common.legend = FALSE),
                                                                         ncol = 2, nrow = 1, common.legend = FALSE, widths = c(0.825,0.175)), 
                                                               ggarrange(Biome_shift_magnitude_2250CE_NH_map, blank_plot, ncol = 2, nrow = 1, widths = c(0.9,0.1), common.legend = FALSE),
                                                               ncol = 1, nrow = 2, heights = c(0.75,0.25)), 
                                                     ggarrange(blank_plot, Bioclimatic_variable_state_change_rate_Biome_shift_magnitude_2250CE_level_barplot, blank_plot, ncol = 3, nrow = 1, widths = c(0.015,0.8,0.185), common.legend = FALSE),
                                                     ncol = 1, nrow = 2, heights = c(0.8,0.2))
  
  # Save the final figure as a PNG file
  ggsave("result/Fig. 4-Potential spatial patterns of bioclimates and biome changes in the MPI-ESM simulation over the next 300 years in the Northern Hemisphere.png", width = 13, height = 21, units = "in", dpi = 300)
  
}


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 5. Supplementary Fig. 1: Record density for the last 21,000 years  ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
Timeslice_proportion_site_0_21ka           <- read.csv2("data/Supplementary Fig. 1 data.1-The proportion of 500-year timeslices available in the last 21,000 years for each site.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Site_number_biome_0_21ka                   <- read.csv2("data/Supplementary Fig. 1 data.2-Number of available sites at timeslices over the last 21,000 years in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_potential_natural_biomes_5minutes   <- read.csv2("data/Supplementary Fig. 1 data.3-Modern potential natural biomes (spatial resolution 5 arc minutes).csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Timeslice_proportion_site_0_21ka           <- type.convert(Timeslice_proportion_site_0_21ka, as.is = TRUE) 
Site_number_biome_0_21ka                   <- type.convert(Site_number_biome_0_21ka, as.is = TRUE) 
Modern_potential_natural_biomes_5minutes   <- type.convert(Modern_potential_natural_biomes_5minutes, as.is = TRUE) 

# Plot
{
  # (1) Spatial distribution of the proportion of 500-year timeslices available in the last 21,000 years
  {
    # Set Biome as a factor for consistent plotting order
    Modern_potential_natural_biomes_5minutes$Biome  <- factor(Modern_potential_natural_biomes_5minutes$Biome,
                                                              levels = c("1", "2", "3", "4", "5", "6", "7", "8")) 
    
    # Map settings and customization
    {
      # Define x and y axis labels for the map
      xlabs = seq(-150,150,50)
      ylabs = seq(-60,80,20)
      
      #define the cols:
      cols <- c("1"  = "#e31a1c",  "2"  = "#f781bf",  "3" = "#33a02c",   "4" = "#1f78b4",
                "5"  = "#ff7f00",  "6"  = "#fdbf6f",  "7"  = "#ffed6f",  "8" = "#6a3d9a")
      
      # Define colors for biome categories
      # !!! IMPORTANT: USE THE SAME ORDER LIKE THE COLORS !!!
      brks=c("1" ,  "2" ,  "3",  "4", "5" ,  "6" ,  "7",  "8")
      
      # Define biome categories, breaks, and labels
      # !!! IMPORTANT: USE THE SAME ORDER LIKE THE COLORS !!!
      labs=c( "Tropical forest",                  "Subtropical forest",           "Temperate forest",    "Boreal forest",                   
              "(Warm) savanna and dry woodland",  "Grassland and dry shrubland",  "(Warm) desert",       "Tundra and polar desert")
      
      # Load coastline data
      coastline <- ne_coastline(scale = "medium", returnclass = "sf")
      
    }
    
    
    Timeslice_proportion_site_0_21ka_map <- ggplot() +
      geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
      geom_tile(data= Modern_potential_natural_biomes_5minutes, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
      scale_fill_manual(name="", labels=labs, breaks=brks, values=cols) +
      labs(title= '(a) Spatial distribution of the proportion of 500-year timeslices available in the last 21,000 years', x="", y="") +
      scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
      geom_point(data= Timeslice_proportion_site_0_21ka, aes(x=Longitude, y=Latitude, color= Percentage), size = 0.8) +
      scale_colour_gradient2(name = "Percentage", low = "blue", mid = "yellow", high = "red", midpoint = median(Timeslice_proportion_site_0_21ka$Percentage),
                             limits = c(range(Timeslice_proportion_site_0_21ka$Percentage)[1], range(Timeslice_proportion_site_0_21ka$Percentage)[2]),
                             breaks = seq(0.2,1,by=0.2), labels = seq(0.2,1,by=0.2)) +
      coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=14, hjust = 0, face = 'bold')) +
      theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt")),
            axis.text.x=element_text(size=12)) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
            axis.text.y=element_text(size=12)) +
      theme(legend.position="right",
            legend.text=element_text(size=12),
            legend.title = element_text(size=12, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
      guides(color =guide_colorbar(direction = "vertical", barwidth = 1, barheight = 20, title.position = "top", title.hjust = 0.5, label.hjust = 0.5, label.theme = element_text(size = 14)),
             fill=guide_legend(override.aes=list(), nrow=2, byrow=TRUE))
    
    guide_color <- get_legend(Timeslice_proportion_site_0_21ka_map + guides(fill = "none"))
    
    Timeslice_proportion_site_0_21ka_map_final <- plot_grid(Timeslice_proportion_site_0_21ka_map +
                                                              guides(color = "none") +
                                                              theme(legend.position = "bottom"),
                                                            guide_color,
                                                            rel_widths = c(1, 0.075)) 
    
  }
  
  
  # (2) Number of available sites at timeslices over the last 21,000 years in the Northern Hemisphere
  {
    # Set Biome as a factor for consistent plotting order
    Site_number_biome_0_21ka$Biome <- factor(Site_number_biome_0_21ka$Biome,
                                             levels = c("1", "2", "3", "4", "5", "6", "7", "8"),
                                             labels = c("TRFO", "WTFO", "TEFO", "BOFO", "SAVA", "STEP", "DESE", "TUND"))
    
    Site_number_biome_0_21ka_fig <- ggplot(data=Site_number_biome_0_21ka, aes(x= Timeslice, y = Site_number)) +
      geom_line(aes(color = Biome, group = Biome), linetype="solid", size=0.8) +
      scale_y_log10() +
      scale_color_manual(name="", values=c("TRFO"="#e31a1c", "WTFO"="#f781bf", "TEFO"="#33a02c", "BOFO"="#1f78b4", "SAVA"="#ff7f00", "STEP"="#fdbf6f", "DESE"="#ffed6f", "TUND"="#6a3d9a"),
                         labels = c( "Tropical forest",                  "Subtropical forest",           "Temperate forest",    "Boreal forest",                   
                                     "(Warm) savanna and dry woodland",  "Grassland and dry shrubland",  "(Warm) desert",       "Tundra and polar desert")) +
      scale_x_continuous(limits = c(0, 21), breaks = seq(0,21, by=3), labels = seq(0,21, by=3), expand = c(0,0)) +
      labs(x="Age (ka cal BP)", y = "Number of sites (log)", title = "(b) Number of available sites at timeslices over the last 21,000 years") +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=14, hjust = 0, face = 'bold')) +
      theme(axis.title.x=element_text(size=14, margin=unit(c(5, 0, 0, 0), "pt")),
            axis.text.x=element_text(size=12)) +
      theme(axis.title.y=element_text(size=14, angle=90, margin=unit(c(0, 5, 0, 0), "pt")),
            axis.text.y=element_text(size=12)) +
      theme(legend.position="bottom",
            legend.text=element_text(size=12),
            legend.title = element_text(size=14, vjust = 1, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm")) +
      guides(color=guide_legend(override.aes=list(size=5), nrow=2, byrow=TRUE))
    
    Site_number_biome_0_21ka_fig
  }
  
  # -------------------------------------------------------------------------------------------------
  # Arranging plots
  blank_plot <- ggplot() + theme_void()
  
  Record_density_0_21ka_fig  <- ggarrange(Timeslice_proportion_site_0_21ka_map_final, 
                                          ggarrange(blank_plot, Site_number_biome_0_21ka_fig, blank_plot,
                                                    nrow = 1, ncol = 3, widths = c(0.2,0.6, 0.2)),
                                          nrow = 2, ncol = 1, common.legend = F)
  
  # Save the final figure as a PNG file
  ggsave("result/Supplementary Fig. 1-Record density for the last 21,000 years.png", width = 13, height = 13, units = "in", dpi = 300)
  
}


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 6. Supplementary Fig. 2: Spatial patterns of biome distributions at 0 cal. ka BP and their agreement with modern potential natural biomes  ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
Modern_potential_natural_biomes_5minutes      <- read.csv2("data/Supplementary Fig. 1 data.3-Modern potential natural biomes (spatial resolution 5 arc minutes).csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_potential_natural_biomes_3.75degrees   <- read.csv2("data/Supplementary Fig. 2 data.2-Modern potential natural biomes (spatial resolution 3.75 degrees).csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Biomes_reconstruction_0ka                     <- read.csv2("data/Supplementary Fig. 2 data.3-Reconstructed biome at 0 ka vs. Modern potential natural biomes (spatial resolution 5 arc minutes).csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Biomes_simulation_0ka                         <- read.csv2("data/Supplementary Fig. 2 data.4-Simulated biome at 0 ka vs. Modern potential natural biomes (spatial resolution 3.75 degrees).csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Modern_potential_natural_biomes_5minutes      <- type.convert(Modern_potential_natural_biomes_5minutes, as.is = TRUE) 
Modern_potential_natural_biomes_3.75degrees   <- type.convert(Modern_potential_natural_biomes_3.75degrees, as.is = TRUE) 
Biomes_reconstruction_0ka                     <- type.convert(Biomes_reconstruction_0ka, as.is = TRUE) 
Biomes_simulation_0ka                         <- type.convert(Biomes_simulation_0ka, as.is = TRUE) 

# Plot
{
  # Convert biome categories to factors for consistent mapping and visualization
  Modern_potential_natural_biomes_5minutes$Biome <- factor(Modern_potential_natural_biomes_5minutes$Biome,
                                                           levels = c("1", "2", "3", "4", "5", "6", "7", "8")) 
  Modern_potential_natural_biomes_3.75degrees$Biome  <- factor(Modern_potential_natural_biomes_3.75degrees$Biome,
                                                               levels = c("1", "2", "3", "4", "5", "6", "7", "8")) 
  
  # Map settings and customization
  {
    # Define x and y axis labels for the map
    xlabs = seq(-150,150,50)
    ylabs = seq(-60,80,20)
    
    #define the cols:
    cols <- c("1"  = "#e31a1c",  "2"  = "#f781bf",  "3" = "#33a02c",   "4" = "#1f78b4",
              "5"  = "#ff7f00",  "6"  = "#fdbf6f",  "7"  = "#ffed6f",  "8" = "#6a3d9a")
    
    # Define colors for biome categories
    # !!! IMPORTANT: USE THE SAME ORDER LIKE THE COLORS !!!
    brks=c("1" ,  "2" ,  "3",  "4", "5" ,  "6" ,  "7",  "8")
    
    # Define biome categories, breaks, and labels
    # !!! IMPORTANT: USE THE SAME ORDER LIKE THE COLORS !!!
    labs=c( "Tropical forest",                  "Subtropical forest",           "Temperate forest",    "Boreal forest",                   
            "(Warm) savanna and dry woodland",  "Grassland and dry shrubland",  "(Warm) desert",       "Tundra and polar desert")
    
    # Load coastline data
    coastline <- ne_coastline(scale = "medium", returnclass = "sf")
    
  }
  
  
  # (1) Spatial patterns of biome distributions at 0 cal. ka BP
  {
    # Convert reconstruction and simulation categories to factors for visualization
    Biomes_reconstruction_0ka$Reconstruction <- factor(Biomes_reconstruction_0ka$Reconstruction,
                                                       levels = c("1", "2", "3", "4", "5", "6", "7", "8")) 
    Biomes_simulation_0ka$Simulation    <- factor(Biomes_simulation_0ka$Simulation,
                                                  levels = c("1", "2", "3", "4", "5", "6", "7", "8"))
    {
      Biome_0ka_pollen_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_biomes_5minutes, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Biome", labels=labs, breaks=brks, values=cols) +
        labs(title= '(a) Biome distributions at 0 cal. ka BP from pollen-based reconstruction', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Biomes_reconstruction_0ka, aes(x=Longitude, y=Latitude, color= Reconstruction), size = 0.5) +
        scale_colour_manual(name="Biome", labels=labs, breaks=brks, values=cols) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=13, hjust = 0, face = 'bold')) +
        theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt")),
              axis.text.x=element_text(size=12)) +
        theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
              axis.text.y=element_text(size=12)) +
        theme(legend.position="bottom",
              legend.text=element_text(size=12),
              legend.title = element_text(size=12, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) +
        guides(color=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE),
               fill=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE)) 
      
      Biome_0ka_pollen_site_map
      
    }
    
    {
      Biome_0ka_MPI_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_biomes_3.75degrees, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Biome", labels=labs, breaks=brks, values=cols) +
        labs(title= '(b) Biome distributions at 0 cal. ka BP from ESM-based simulation', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Biomes_simulation_0ka, aes(x=Longitude, y=Latitude, color= Simulation), size = 0.5) +
        scale_colour_manual(name="Biome", labels=labs, breaks=brks, values=cols) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=13, hjust = 0, face = 'bold')) +
        theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt")),
              axis.text.x=element_text(size=12)) +
        theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
              axis.text.y=element_text(size=12)) +
        theme(legend.position="bottom",
              legend.text=element_text(size=12),
              legend.title = element_text(size=12, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) +
        guides(color=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE),
               fill=guide_legend(override.aes=list(size=6), nrow=2, byrow=TRUE)) 
      
      Biome_0ka_MPI_site_map
      
    }
    
  }
  
  
  # (2) Spatial patterns of biome agreement with modern potential natural biomes
  {
    # Convert reconstruction and simulation categories to factors for visualization
    Biomes_reconstruction_0ka$Agreement <- factor(Biomes_reconstruction_0ka$Agreement,
                                                  levels = c("yes", "no"))
    Biomes_simulation_0ka$Agreement    <- factor(Biomes_simulation_0ka$Agreement,
                                                 levels = c("yes", "no"))
  
    {
      Biome_agreement_0ka_pollen_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_biomes_5minutes, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Biome", labels=labs, breaks=brks, values=cols, guide = "none") +
        labs(title= '(c) Biome agreement of pollen-based reconstruction', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Biomes_reconstruction_0ka, aes(x=Longitude, y=Latitude, color= Agreement), size = 0.5) +
        scale_colour_manual(name="Status", values=c("yes" = "yellow", "no" = "red"), labels= c("yes" = "Agree", "no" = "Disagree")) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=13, hjust = 0, face = 'bold')) +
        theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt")),
              axis.text.x=element_text(size=12)) +
        theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
              axis.text.y=element_text(size=12)) +
        theme(legend.position="bottom",
              legend.text=element_text(size=12),
              legend.title = element_text(size=12, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) +
        guides(color=guide_legend(override.aes=list(size=6), nrow=1, byrow=TRUE)) 
      
      Biome_agreement_0ka_pollen_site_map
      
    }
  
    
    {
      Biome_agreement_0ka_MPI_site_map <- ggplot() +
        geom_sf(data = coastline, fill = NA, color = "grey80", size = 0.2) +
        geom_tile(data= Modern_potential_natural_biomes_3.75degrees, aes(x=Longitude, y=Latitude, fill = Biome), alpha = 0.3) +
        scale_fill_manual(name="Biome", labels=labs, breaks=brks, values=cols, guide = "none") +
        labs(title= '(d) Biome agreement of ESM-based simulation', x="", y="") +
        scale_x_continuous(limits = c(-180, 180), breaks = seq(-150,150,50), labels = paste0(seq(-150,150,50), "°E"), expand = c(0, 0)) +
        geom_point(data= Biomes_simulation_0ka, aes(x=Longitude, y=Latitude, color= Agreement), size = 0.5) +
        scale_colour_manual(name="Status", values=c("yes" = "yellow", "no" = "red"), labels= c("yes" = "Agree", "no" = "Disagree")) +
        coord_sf(ylim = c(-60, 83.5), expand = FALSE) +
        theme_bw() +  # Reset the plot background
        theme(plot.title=element_text(size=13, hjust = 0, face = 'bold')) +
        theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt")),
              axis.text.x=element_text(size=12)) +
        theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
              axis.text.y=element_text(size=12)) +
        theme(legend.position="bottom",
              legend.text=element_text(size=12),
              legend.title = element_text(size=12, vjust = 0.5, hjust = 0.5, face = 'bold'),
              legend.key.size = unit(1.2, "lines"),
              legend.key.height = unit(1, "lines"),
              legend.key.width = unit(2, "lines"),
              legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
              legend.box = "vertical",
              legend.justification = c(0.5, 0.5)) +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(linetype="solid", fill=NA),
              panel.background=element_rect(fill="white")) +
        theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) +
        guides(color=guide_legend(override.aes=list(size=6), nrow=1, byrow=TRUE)) 
      
      Biome_agreement_0ka_MPI_site_map
      
    }
    
  }
  
  #Arranging plots
  PNV_Biome_agreement_pollen_MPI_site_global_0ka_map   <- ggarrange(ggarrange(Biome_0ka_pollen_site_map, Biome_0ka_MPI_site_map, 
                                                                              nrow = 1, ncol = 2, common.legend = T, legend = "bottom"),
                                                                    ggarrange(Biome_agreement_0ka_pollen_site_map, Biome_agreement_0ka_MPI_site_map,
                                                                              nrow = 1, ncol = 2, common.legend = T, legend = "bottom"),
                                                                    nrow = 2, ncol = 1, heights = c(0.525, 0.475),common.legend = F)
  
  # Save the final figure as a PNG file
  ggsave("result/Supplementary Fig. 2-Spatial patterns of biome distributions at 0 cal. ka BP and their agreement with modern potential natural biomes.png", width = 14, height = 8.5, units = "in", dpi = 300)
  
}


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 7. Supplementary Fig. 3: Temporal patterns of simulated bioclimatic variables over the last 21,000 years in the Northern Hemisphere  ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
Bioclimatic_variables_state_standardization_NH_0_21ka        <- read.csv2("data/Supplementary Fig. 3 data.1-Standardization of Simulated bioclimatic variable state over the last 21,000 years in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Bioclimatic_variables_change_rate_standardization_NH_0_21ka  <- read.csv2("data/Supplementary Fig. 3 data.2-Standardization of Simulated bioclimatic variable change rate over the last 21,000 years in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Bioclimatic_variables_state_standardization_NH_0_21ka        <- type.convert(Bioclimatic_variables_state_standardization_NH_0_21ka, as.is = TRUE) 
Bioclimatic_variables_change_rate_standardization_NH_0_21ka  <- type.convert(Bioclimatic_variables_change_rate_standardization_NH_0_21ka, as.is = TRUE) 

# Plot
{
  # (1) Hierarchical clustering on the time series of the 14 hemispheric bioclimatic variables
  {
    # Covert dataset
    Bioclimatic_variables_state_standardization_NH_0_21ka_t <- Bioclimatic_variables_state_standardization_NH_0_21ka[ ,2:15] %>%
      t() %>%
      as.data.frame() %>%
      setNames(Bioclimatic_variables_state_standardization_NH_0_21ka$Timeslice)
    
    # Hierarchical clustering
    {
      {
        # Perform hierarchical clustering for temperature-related variables
        Temperature_state_clustering <- dtwclust::tsclust(series = Bioclimatic_variables_state_standardization_NH_0_21ka_t[c(1,3:6), ], 
                                                          type = "hierarchical", 
                                                          k = 2,  preproc = zscore,
                                                          distance = "dtw", centroid = shape_extraction,
                                                          control = hierarchical_control(method = "complete"))
        
        # Extract cluster groups and dendrogram data
        Temperature_state_clustering_hclus <- stats::cutree(Temperature_state_clustering, k = 2) %>% 
          as.data.frame(.) %>%
          dplyr::rename(.,cluster_group = .) %>%
          tibble::rownames_to_column("type_col") 
        
        Temperature_state_clustering_hcdata <- ggdendro::dendro_data(Temperature_state_clustering)
        names_order <- Temperature_state_clustering_hcdata$labels$label
        
        # Plot temperature clustering dendrogram
        Temperature_state_clustering_fig <- Temperature_state_clustering_hcdata %>%
          ggdendro::ggdendrogram(., rotate=TRUE, leaf_labels=F, labels = F) +
          theme(axis.text.x = element_blank()) +
          theme(axis.text.y = element_blank())
      }
      
      {
        # Perform hierarchical clustering for precipitation-related variables
        Precipitation_state_clustering <- dtwclust::tsclust(series = Bioclimatic_variables_state_standardization_NH_0_21ka_t[c(7:9,11:14), ], 
                                                            type = "hierarchical", 
                                                            k = 2,  preproc = zscore,
                                                            distance = "dtw", centroid = shape_extraction,
                                                            control = hierarchical_control(method = "complete"))
        
        # Extract cluster groups and dendrogram data
        Precipitation_state_clustering_hclus <- stats::cutree(Precipitation_state_clustering, k = 2) %>% 
          as.data.frame(.) %>%
          dplyr::rename(.,cluster_group = .) %>%
          tibble::rownames_to_column("type_col") 
        
        Precipitation_state_clustering_hcdata <- ggdendro::dendro_data(Precipitation_state_clustering)
        names_order <- Precipitation_state_clustering_hcdata$labels$label
        
        # Plot precipitation clustering dendrogram
        Precipitation_state_clustering_fig <- Precipitation_state_clustering_hcdata %>%
          ggdendro::ggdendrogram(., rotate=TRUE, leaf_labels=F, labels = F) +
          theme(axis.text.x = element_blank()) +
          theme(axis.text.y = element_blank())
        
      }
      
     }
    
  }
  
  
  # (2) Bioclimatic variables
  {
    ## Temperature-related
    {
      # Reshape dataset for visualization
      Temperature_state_standardization_NH_0_21ka_long <- pivot_longer(
        Bioclimatic_variables_state_standardization_NH_0_21ka[ ,c(1:2,4:7)],
        cols = starts_with("BIO"),         
        names_to = "BIO",                  
        values_to = "Value")
      
      # Factorize BIO variable with descriptive labels
      Temperature_state_standardization_NH_0_21ka_long$BIO <- factor(Temperature_state_standardization_NH_0_21ka_long$BIO,
                                                          levels = c("BIO11", "BIO1", "BIO9", "BIO10", "BIO8"),
                                                          labels = c("BIO 11", "BIO 1", "BIO 9", "BIO 10", "BIO 8"))
      
      {
        Temperature_state_standardization_NH_0_21ka_fig <- ggplot(data=Temperature_state_standardization_NH_0_21ka_long, aes(x= Timeslice, y = Value)) +
          geom_line(aes(color = BIO, group = BIO), linetype="solid", size=1) +
          facet_grid(BIO ~ ., scales = "free_y", switch = "y",
                     labeller = as_labeller(c("BIO 8"  = "Mean temperature of wettest quarter", 
                                              "BIO 10" = "Mean temperature of warmest quarter", 
                                              "BIO 9"  = "Mean temperature of driest quarter", 
                                              "BIO 11" = "Mean temperature of coldest quarter", 
                                              "BIO 1"  = "Annual mean temperature"))) +
          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by=0.5), labels = seq(0, 1, by=0.5), expand = c(0,0)) +
          scale_x_continuous(limits = c(0, 20.5), breaks = seq(0,20.5, by=2), expand = c(0,0),
                             labels= c("0-0.5", "2-2.5", "4-4.5", "6-6.5", "8-8.5", "10-10.5", "12-12.5", "14-14.5", "16-16.5", "18-18.5", "20-20.5")) +
          scale_color_manual(name="", values=c("BIO 8"="#e31a1c", "BIO 10"="#e31a1c", "BIO 9"="#fb9a99", "BIO 11"="#fb9a99", "BIO 1"="#fb9a99")) +
          labs(title= '', x="", y="") +
          theme_bw() +  # Reset the plot background
          theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
          theme(axis.title=element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold')) +
          theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt")),
                axis.text.x=element_blank()) +
          theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 5, 0, 0), "pt")),
                axis.text.y=element_text(size=13)) + 
          theme(legend.position="none",
                legend.text=element_text(size=6),
                legend.title = element_text(size=8, vjust = 0.5, hjust = 0.5),
                legend.key.size = unit(1.2, "lines"),
                legend.key.height = unit(1, "lines"),
                legend.key.width = unit(2, "lines"),
                legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
                legend.box = "vertical",
                legend.justification = c(0.5, 0.5)) +
          theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.border=element_rect(linetype="solid", fill=NA),
                panel.background=element_rect(fill="white"),
                panel.spacing = unit(0.5, "cm")) +
          theme(strip.background=element_blank(),
                strip.text=element_text(size=14, face="bold"),
                strip.text.y.left=element_text(size =14, vjust = 0.5, hjust = 0, angle=0, margin=unit(c(0, 5, 0, 0), "pt")),
                strip.placement="outside") +
          theme(plot.margin=unit(c(0,0.5,0,0.25),"cm")) 
        
        Temperature_state_standardization_NH_0_21ka_fig
        
      }
      
    }
    
    
    ## Precipitation-related
    {
      # Reshape dataset for visualization
      Precipitation_state_standardization_NH_0_21ka_long <- pivot_longer(
        Bioclimatic_variables_state_standardization_NH_0_21ka[ ,c(1,8:10,12:15)],
        cols = starts_with("BIO"),         
        names_to = "BIO",                  
        values_to = "Value")
      
      # Factorize BIO variable with descriptive labels
      Precipitation_state_standardization_NH_0_21ka_long$BIO <- factor(Precipitation_state_standardization_NH_0_21ka_long$BIO,
                                                            levels = c("BIO16", "BIO13", "BIO12", "BIO18", "BIO17", "BIO14", "BIO19"),
                                                            labels = c("BIO 16", "BIO 13", "BIO 12", "BIO 18", "BIO 17", "BIO 14", "BIO 19"))
      
      {
        Precipitation_state_standardization_NH_0_21ka_fig <- ggplot(data=Precipitation_state_standardization_NH_0_21ka_long, aes(x= Timeslice, y = Value)) +
          geom_line(aes(color = BIO, group = BIO), linetype="solid", size=1) +
          facet_grid(BIO ~ ., scales = "free_y", switch = "y",
                     labeller = as_labeller(c("BIO 16" = "Precipitation of wettest quarter", 
                                              "BIO 13" = "Precipitation of wettest month", 
                                              "BIO 12" = "Annual precipitation", 
                                              "BIO 18" = "Precipitation of warmest quarter", 
                                              "BIO 17" = "Precipitation of driest quarter", 
                                              "BIO 14" = "Precipitation of driest month", 
                                              "BIO 19" = "Precipitation of coldest quarter"))) +
          #ggh4x::facetted_pos_scales(y = y_scale_list_precipitation) +
          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by=0.5), labels = seq(0, 1, by=0.5), expand = c(0,0)) +
          scale_x_continuous(limits = c(0, 20.5), breaks = seq(0,20.5, by=2), expand = c(0,0),
                             labels= c("0-0.5", "2-2.5", "4-4.5", "6-6.5", "8-8.5", "10-10.5", "12-12.5", "14-14.5", "16-16.5", "18-18.5", "20-20.5")) +
          scale_color_manual(name="", values=c("BIO 16"="#1f78b4", "BIO 13"="#1f78b4", "BIO 12"="#1f78b4", "BIO 18"="#1f78b4", 
                                               "BIO 17"="#a6cee3", "BIO 14"="#a6cee3", "BIO 19"="#a6cee3")) +
          labs(title= '', x="Age (ka cal BP)", y="") +
          theme_bw() +  # Reset the plot background
          theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
          theme(axis.title=element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold')) +
          theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt")),
                axis.text.x=element_text(size=13, angle = 45, vjust = 1, hjust = 1)) +
          theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 5, 0, 0), "pt")),
                axis.text.y=element_text(size=13)) + 
          theme(legend.position="none",
                legend.text=element_text(size=6),
                legend.title = element_text(size=8, vjust = 0.5, hjust = 0.5),
                legend.key.size = unit(1.2, "lines"),
                legend.key.height = unit(1, "lines"),
                legend.key.width = unit(2, "lines"),
                legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
                legend.box = "vertical",
                legend.justification = c(0.5, 0.5)) +
          theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.border=element_rect(linetype="solid", fill=NA),
                panel.background=element_rect(fill="white"),
                panel.spacing = unit(0.5, "cm")) +
          theme(strip.background=element_blank(),
                strip.text=element_text(size=14, face="bold"),
                strip.text.y.left=element_text(size =14, vjust = 0.5, hjust = 0, angle=0, margin=unit(c(0, 12*3+5, 0, 0), "pt")),
                strip.placement="outside") +
          theme(plot.margin=unit(c(0,0.5,0,0.25),"cm")) 
        
        Precipitation_state_standardization_NH_0_21ka_fig
        
      }
      
    }
    
    ## Seasonality-related
    {
      # Reshape dataset for visualization
      Seasonality_state_standardization_NH_0_21ka_long <- pivot_longer(
        Bioclimatic_variables_state_standardization_NH_0_21ka[ ,c(1,3,11)],
        cols = starts_with("BIO"),         
        names_to = "BIO",                  
        values_to = "Value")
      
      # Factorize BIO variable with descriptive labels
      Seasonality_state_standardization_NH_0_21ka_long$BIO <- factor(Seasonality_state_standardization_NH_0_21ka_long$BIO,
                                                          levels = c("BIO4", "BIO15"),
                                                          labels = c("BIO 4", "BIO 15"))
      
      {
        Seasonality_state_standardization_NH_0_21ka_fig <- ggplot(data=Seasonality_state_standardization_NH_0_21ka_long, aes(x= Timeslice, y = Value)) +
          geom_line(aes(color = BIO, group = BIO), linetype="solid", size=1) +
          facet_grid(BIO ~ ., scales = "free_y", switch = "y",
                     labeller = as_labeller(c("BIO 4"  = "Temperature seasonality", 
                                              "BIO 15" = "Precipitation seasonality"))) +
          #ggh4x::facetted_pos_scales(y = y_scale_list_seasonality) +
          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by=0.5), labels = seq(0, 1, by=0.5), expand = c(0,0)) +
          scale_x_continuous(limits = c(0, 20.5), breaks = seq(0,20.5, by=2), expand = c(0,0),
                             labels= c("0-0.5", "2-2.5", "4-4.5", "6-6.5", "8-8.5", "10-10.5", "12-12.5", "14-14.5", "16-16.5", "18-18.5", "20-20.5")) +
          scale_color_manual(name="", values=c("BIO 4"="#a65628", "BIO 15"="#a65628")) +
          labs(title= '\n(a) Bioclimatic variables\n ', x="", y="") +
          theme_bw() +  # Reset the plot background
          theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
          theme(axis.title=element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold')) +
          theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt")),
                axis.text.x=element_blank()) +
          theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 5, 0, 0), "pt")),
                axis.text.y=element_text(size=13)) + 
          theme(legend.position="none",
                legend.text=element_text(size=6),
                legend.title = element_text(size=8, vjust = 0.5, hjust = 0.5),
                legend.key.size = unit(1.2, "lines"),
                legend.key.height = unit(1, "lines"),
                legend.key.width = unit(2, "lines"),
                legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
                legend.box = "vertical",
                legend.justification = c(0.5, 0.5)) +
          theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.border=element_rect(linetype="solid", fill=NA),
                panel.background=element_rect(fill="white"),
                panel.spacing = unit(0.5, "cm")) +
          theme(strip.background=element_blank(),
                strip.text=element_text(size=14, face="bold"),
                strip.text.y.left=element_text(size =14, vjust = 0.5, hjust = 0, angle=0, margin=unit(c(0, 12*7+5+3, 0, 0), "pt")),
                strip.placement="outside") +
          theme(plot.margin=unit(c(0,0.5,0,0.25),"cm")) 
        
        Seasonality_state_standardization_NH_0_21ka_fig
        
      }
      
    }
    
  }
  
  
  # (3) Rate of change of bioclimatic variables
  {
    ## Temperature-related
    {
      # Reshape dataset for visualization
      Temperature_change_rate_standardization_NH_0_21ka_long <- pivot_longer(
        Bioclimatic_variables_change_rate_standardization_NH_0_21ka[,c(1:2,4:7)],
        cols = starts_with("BIO"),         
        names_to = "BIO",                  
        values_to = "Value")
      
      # Factorize BIO variable with descriptive labels
      Temperature_change_rate_standardization_NH_0_21ka_long$BIO <- factor(Temperature_change_rate_standardization_NH_0_21ka_long$BIO,
                                                                           levels = c("BIO11", "BIO1", "BIO9", "BIO10", "BIO8"),
                                                                           labels = c("BIO 11", "BIO 1", "BIO 9", "BIO 10", "BIO 8"))
      
      {
        Temperature_change_rate_standardization_NH_0_21ka_fig <- ggplot(data=Temperature_change_rate_standardization_NH_0_21ka_long, aes(x= Timeslice, y = Value)) +
          geom_line(aes(color = BIO, group = BIO), linetype="solid", size=1) +
          facet_grid(BIO ~ ., scales = "free_y", switch = "y") +
          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by=0.5), labels = seq(0, 1, by=0.5), expand = c(0,0)) +
          scale_x_continuous(limits = c(0, 20.5), breaks = seq(0,20.5, by=2), expand = c(0,0),
                             labels= c("0-0.5", "2-2.5", "4-4.5", "6-6.5", "8-8.5", "10-10.5", "12-12.5", "14-14.5", "16-16.5", "18-18.5", "20-20.5")) +
          scale_color_manual(name="", values=c("BIO 8"="#e31a1c", "BIO 10"="#e31a1c", "BIO 9"="#fb9a99", "BIO 11"="#fb9a99", "BIO 1"="#fb9a99")) +
          labs(title= '', x="", y="") +
          theme_bw() +  # Reset the plot background
          theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
          theme(axis.title=element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold')) +
          theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt")),
                axis.text.x=element_blank()) +
          theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 5, 0, 0), "pt")),
                axis.text.y=element_text(size=13)) + 
          theme(legend.position="none",
                legend.text=element_text(size=6),
                legend.title = element_text(size=8, vjust = 0.5, hjust = 0.5),
                legend.key.size = unit(1.2, "lines"),
                legend.key.height = unit(1, "lines"),
                legend.key.width = unit(2, "lines"),
                legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
                legend.box = "vertical",
                legend.justification = c(0.5, 0.5)) +
          theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.border=element_rect(linetype="solid", fill=NA),
                panel.background=element_rect(fill="white"),
                panel.spacing = unit(0.5, "cm")) +
          theme(strip.background=element_blank(),
                strip.text=element_text(size=14, face="bold"),
                strip.text.y.left=element_blank(),
                strip.placement="outside") +
          theme(plot.margin=unit(c(0,0.5,0,0.25),"cm")) 
        
        Temperature_change_rate_standardization_NH_0_21ka_fig
        
      }
      
    }
    
    
    ## Precipitation-related
    {
      # Reshape dataset for visualization
      Precipitation_change_rate_standardization_NH_0_21ka_long <- pivot_longer(
        Bioclimatic_variables_change_rate_standardization_NH_0_21ka[,c(1,8:10,12:15)],
        cols = starts_with("BIO"),         
        names_to = "BIO",                  
        values_to = "Value")
      
      # Factorize BIO variable with descriptive labels
      Precipitation_change_rate_standardization_NH_0_21ka_long$BIO <- factor(Precipitation_change_rate_standardization_NH_0_21ka_long$BIO,
                                                            levels = c("BIO16", "BIO13", "BIO12", "BIO18", "BIO17", "BIO14", "BIO19"),
                                                            labels = c("BIO 16", "BIO 13", "BIO 12", "BIO 18", "BIO 17", "BIO 14", "BIO 19"))
      
      {
        Precipitation_change_rate_standardization_NH_0_21ka_fig <- ggplot(data=Precipitation_change_rate_standardization_NH_0_21ka_long, aes(x= Timeslice, y = Value)) +
          geom_line(aes(color = BIO, group = BIO), linetype="solid", size=1) +
          facet_grid(BIO ~ ., scales = "free_y", switch = "y") +
          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by=0.5), labels = seq(0, 1, by=0.5), expand = c(0,0)) +
          scale_x_continuous(limits = c(0, 20.5), breaks = seq(0,20.5, by=2), expand = c(0,0),
                             labels= c("0-0.5", "2-2.5", "4-4.5", "6-6.5", "8-8.5", "10-10.5", "12-12.5", "14-14.5", "16-16.5", "18-18.5", "20-20.5")) +
          scale_color_manual(name="", values=c("BIO 16"="#1f78b4", "BIO 13"="#1f78b4", "BIO 12"="#1f78b4", "BIO 18"="#1f78b4", 
                                               "BIO 17"="#a6cee3", "BIO 14"="#a6cee3", "BIO 19"="#a6cee3")) +
          labs(title= '', x="Age (ka cal BP)", y="") +
          theme_bw() +  # Reset the plot background
          theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
          theme(axis.title=element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold')) +
          theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt")),
                axis.text.x=element_text(size=13, angle = 45, vjust = 1, hjust = 1)) +
          theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 5, 0, 0), "pt")),
                axis.text.y=element_text(size=13)) + 
          theme(legend.position="none",
                legend.text=element_text(size=6),
                legend.title = element_text(size=8, vjust = 0.5, hjust = 0.5),
                legend.key.size = unit(1.2, "lines"),
                legend.key.height = unit(1, "lines"),
                legend.key.width = unit(2, "lines"),
                legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
                legend.box = "vertical",
                legend.justification = c(0.5, 0.5)) +
          theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.border=element_rect(linetype="solid", fill=NA),
                panel.background=element_rect(fill="white"),
                panel.spacing = unit(0.5, "cm")) +
          theme(strip.background=element_blank(),
                strip.text=element_text(size=14, face="bold"),
                strip.text.y.left=element_blank(),
                strip.placement="outside") +
          theme(plot.margin=unit(c(0,0.5,0,0.25),"cm")) 
        
        Precipitation_change_rate_standardization_NH_0_21ka_fig
        
      }
      
    }
    
    ## Seasonality-related
    {
      # Reshape dataset for visualization
      Seasonality_change_rate_standardization_NH_0_21ka_long <- pivot_longer(
        Bioclimatic_variables_change_rate_standardization_NH_0_21ka[,c(1,3,11)],
        cols = starts_with("BIO"),         
        names_to = "BIO",                  
        values_to = "Value")
      
      # Factorize BIO variable with descriptive labels
      Seasonality_change_rate_standardization_NH_0_21ka_long$BIO <- factor(Seasonality_change_rate_standardization_NH_0_21ka_long$BIO,
                                                          levels = c("BIO4", "BIO15"),
                                                          labels = c("BIO 4", "BIO 15"))
      
      {
        Seasonality_change_rate_standardization_NH_0_21ka_fig <- ggplot(data=Seasonality_change_rate_standardization_NH_0_21ka_long, aes(x= Timeslice, y = Value)) +
          geom_line(aes(color = BIO, group = BIO), linetype="solid", size=1) +
          facet_grid(BIO ~ ., scales = "free_y", switch = "y") +
          scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by=0.5), labels = seq(0, 1, by=0.5), expand = c(0,0)) +
          scale_x_continuous(limits = c(0, 20.5), breaks = seq(0,20.5, by=2), expand = c(0,0),
                             labels= c("0-0.5", "2-2.5", "4-4.5", "6-6.5", "8-8.5", "10-10.5", "12-12.5", "14-14.5", "16-16.5", "18-18.5", "20-20.5")) +
          scale_color_manual(name="", values=c("BIO 4"="#a65628", "BIO 15"="#a65628")) +
          labs(title= '\n(b) Change rate of bioclimatic variables\n ', x="", y="") +
          theme_bw() +  # Reset the plot background
          theme(plot.title=element_text(size=16, hjust = 0, face = 'bold')) +
          theme(axis.title=element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold')) +
          theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt")),
                axis.text.x=element_blank()) +
          theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 5, 0, 0), "pt")),
                axis.text.y=element_text(size=13)) + 
          theme(legend.position="none",
                legend.text=element_text(size=6),
                legend.title = element_text(size=8, vjust = 0.5, hjust = 0.5),
                legend.key.size = unit(1.2, "lines"),
                legend.key.height = unit(1, "lines"),
                legend.key.width = unit(2, "lines"),
                legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
                legend.box = "vertical",
                legend.justification = c(0.5, 0.5)) +
          theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.border=element_rect(linetype="solid", fill=NA),
                panel.background=element_rect(fill="white"),
                panel.spacing = unit(0.5, "cm")) +
          theme(strip.background=element_blank(),
                strip.text=element_text(size=14, face="bold"),
                strip.text.y.left=element_blank(),
                strip.placement="outside") +
          theme(plot.margin=unit(c(0,0.5,0,0.25),"cm")) 
        
        Seasonality_change_rate_standardization_NH_0_21ka_fig
        
      }
      
    }
    
  }
  
  
  # -------------------------------------------------------------------------------------------------
  
  # Arrange plots
  blank_plot <- ggplot() + theme_void()
  
  simulated_climate_state_shift_site_NH_0_21ka <- ggarrange(ggarrange(
    ggarrange(Seasonality_state_standardization_NH_0_21ka_fig, blank_plot,
              ncol = 2, nrow = 1, common.legend = FALSE, widths = c(0.9,0.1)),
    ggarrange(Temperature_state_standardization_NH_0_21ka_fig, 
              ggarrange(blank_plot, Temperature_state_clustering_fig, blank_plot,
                        ncol = 1, nrow = 3, common.legend = FALSE, heights = c(0.125, 0.8, 0.075)),
              ncol = 2, nrow = 1, common.legend = FALSE, widths = c(0.9,0.1)),
    ggarrange(Precipitation_state_standardization_NH_0_21ka_fig, 
              ggarrange(blank_plot, Precipitation_state_clustering_fig, blank_plot,
                        ncol = 1, nrow = 3, common.legend = FALSE, heights = c(0.05, 0.85, 0.1)),
              ncol = 2, nrow = 1, common.legend = FALSE, widths = c(0.9,0.1)),
    ncol = 1, nrow = 3, heights = c(2/14, 5/14, 7/14), common.legend = FALSE), 
    ggarrange(Seasonality_change_rate_standardization_NH_0_21ka_fig, Temperature_change_rate_standardization_NH_0_21ka_fig, Precipitation_change_rate_standardization_NH_0_21ka_fig,
              ncol = 1, nrow = 3, heights = c(2/14, 5/14, 7/14), common.legend = FALSE),
    ncol = 2, nrow = 1, widths = c(0.65,0.35), common.legend = FALSE)   
  
  # Save the final figure as a PNG file
  ggsave("result/Supplementary Fig. 3-Temporal patterns of simulated bioclimatic variables over the last 21,000 years in the Northern Hemisphere.png", width = 15, height = 21, units = "in", dpi = 300)
  
}


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 8. Supplementary Fig. 4: Gaussian kernel correlations between compositional turnover and proportion of shifts in biomes over the last 21,000 years in the Northern Hemisphere  ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
GKC_biome_NH_0_21ka        <- read.csv2("data/Supplementary Fig. 4 data.1-Gaussian kernel correlations between compositional turnover and proportion of shifts in biomes over the last 21,000 years in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
GKC_biome_NH_0_21ka        <- type.convert(GKC_biome_NH_0_21ka, as.is = TRUE) 

# Plot
{
  # Set custom order for Biome levels for plotting
  GKC_biome_NH_0_21ka$Biome <- factor(GKC_biome_NH_0_21ka$Biome,
                                      levels = c("Hemisphere", "Tropical forest", "Subtropical forest", "Temperate forest",
                                                 "Boreal forest", "Savanna", "Steppe", "Desert", "Tundra"))
 
  # Set custom order and labels for Variable levels for plotting
  GKC_biome_NH_0_21ka$Variable <- factor(GKC_biome_NH_0_21ka$Variable,
                                         levels = rev(c("Biome compositional turnover\nvs.\nReconstructed biome shift rate", "Biome compositional turnover\nvs.\nSimulated biome shift rate",
                                                        "Reconstructed biome shift rate\nvs.\nSimulated biome shift rate")),
                                         labels = rev(c("Biome compositional turnover\nvs.\nProportion of biome shifts in reconstruction", "Biome compositional turnover\nvs.\nProportion of biome shifts in simulation",
                                                        "Proportion of biome shifts in reconstruction\nvs.\nProportion of biome shifts in simulation")))
  
  GKC_biome_NH_0_21ka_fig <- ggplot(data = GKC_biome_NH_0_21ka, aes(x = Biome, y = Variable)) +
    geom_tile(aes(fill = R_value), alpha = 0.6) +  
    geom_text(aes(label = ifelse(P_value < 0.001, paste0(sprintf("%.2f", R_value), " ***"),
                                 ifelse(P_value < 0.01, paste0(sprintf("%.2f", R_value), " **"),
                                        ifelse(P_value < 0.05, paste0(sprintf("%.2f", R_value), " *"), sprintf("%.2f", R_value)))),
                  vjust = 1, fontface = ifelse(abs(R_value) == R_value_MaX, "bold", "plain")), size = 7)  +
    scale_fill_gradient2(name = "R-value", low = "blue", high = "red",
                         limits = c(range(GKC_biome_NH_0_21ka$R_value)[1], range(GKC_biome_NH_0_21ka$R_value)[2]),
                         breaks = seq(-0.2, 0.6, 0.2), labels = c("-0.2", "0.0", "0.2", "0.4", "0.6")) +
    labs(title= '', x="", y="") +
    theme_bw() +  # Reset the plot background
    theme(plot.title=element_text(size=12, hjust = 0, face = 'bold')) +
    theme(axis.title=element_text(size=12)) +
    theme(axis.title.x=element_text(margin=unit(c(5, 0, 0, 0), "pt"), size=12),
          axis.text.x=element_text(size=16, angle=45, vjust = 1, hjust = 1, face = 'bold')) +
    theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
          axis.text.y=element_text(size=18, vjust = 0.5, hjust = 0.5))  +  
    theme(legend.position="none",
          legend.text=element_text(size=12),
          legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
          legend.key.size = unit(1.2, "lines"),
          legend.key.height = unit(1, "lines"),
          legend.key.width = unit(2, "lines"),
          legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
          legend.box = "vertical",
          legend.justification = c(0.5, 0.5)) +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_rect(linetype="solid", fill=NA),
          panel.background=element_rect(fill="white")) +
    theme(strip.background=element_blank(),
          strip.text=element_text(size=10, face="bold"),
          strip.text.y=element_text(angle=0, margin=unit(c(0, 5, 0, 0), "pt")),
          strip.placement="outside") +
    theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))  +
    guides(fill = guide_colorbar(barwidth = unit(0.6, "cm"), barheight = unit(20, "cm")))
  
  GKC_biome_NH_0_21ka_fig
  
  # Save the final figure as a PNG file
  ggsave("result/Supplementary Fig. 4-Gaussian kernel correlations between compositional turnover and proportion of shifts in biomes over the last 21,000 years in the Northern Hemisphere.png", width = 21, height = 13, units = "in", dpi = 300)
  
}


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 9. Supplementary Fig. 5: Gaussian kernel correlation coefficients between compositional turnover and proportion of shifts in biomes with simulated bioclimatic variables and their rates of change over the last 21,000 years in the Northern Hemisphere  ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Load datasets
GKC_biome_climate_state_NH_0_21ka          <- read.csv2("data/Supplementary Fig. 5 data.1-GKC between compositional turnover and proportion of shifts in biomes with simulated bioclimatic variables over the last 21,000 years in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
GKC_biome_climate_change_rate_NH_0_21ka    <- read.csv2("data/Supplementary Fig. 5 data.2-GKC between compositional turnover and proportion of shifts in biomes with bioclimatic variable change rates over the last 21,000 years in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
GKC_biome_climate_state_NH_0_21ka          <- type.convert(GKC_biome_climate_state_NH_0_21ka, as.is = TRUE) 
GKC_biome_climate_change_rate_NH_0_21ka    <- type.convert(GKC_biome_climate_change_rate_NH_0_21ka, as.is = TRUE) 

# Plot
{
  # (1) Bioclimatic variables
  {
    # Convert variables and labels to factors for custom ordering in the plot
    GKC_biome_climate_state_NH_0_21ka$Variable <- factor(GKC_biome_climate_state_NH_0_21ka$Variable,
                                                         levels =c("Biome compositional turnover",  "Reconstructed megabiome shift rate", "Simulated megabiome shift rate"),
                                                         labels = c("Biome compositional turnover", "Proportion of biome shifts in reconstruction", "Proportion of biome shifts in simulation"))
    
    GKC_biome_climate_state_NH_0_21ka$BIO <- factor(GKC_biome_climate_state_NH_0_21ka$BIO,
                                                         levels =rev(c("Temperature seasonality", "Precipitation seasonality",
                                                                   "Mean temperature of coldest quarter", "Annual mean temperature", "Mean temperature of driest quarter",
                                                                   "Mean temperature of warmest quarter", "Mean temperature of wettest quarter",
                                                                   "Precipitation of wettest quarter", "Precipitation of wettest month", "Annual precipitation", "Precipitation of warmest quarter",
                                                                   "Precipitation of driest quarter", "Precipitation of driest month", "Precipitation of coldest quarter")))
    
    GKC_biome_climate_state_NH_0_21ka_fig <- ggplot(data = GKC_biome_climate_state_NH_0_21ka, aes(x = Variable, y = BIO)) +
      geom_tile(aes(fill = R_value), alpha = 0.6) +  
      geom_text(aes(label = ifelse(P_value < 0.001, paste0(sprintf("%.2f", R_value), " ***"),
                                   ifelse(P_value < 0.01, paste0(sprintf("%.2f", R_value), " **"),
                                          ifelse(P_value < 0.05, paste0(sprintf("%.2f", R_value), " *"), sprintf("%.2f", R_value)))),
                    vjust = 1, fontface = ifelse(abs(R_value) == R_value_Max, "bold", "plain")), size = 7) +
      scale_fill_gradient2(name = "R-value", low = "blue", high = "red",
                           limits = c(range(GKC_biome_climate_state_NH_0_21ka$R_value)[1], range(GKC_biome_climate_state_NH_0_21ka$R_value)[2]),
                           breaks = seq(-0.6, 0.8, 0.2), labels = c("-0.6", "-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8")) +
      labs(title= '(a) Bioclimatic variables', x="", y="") +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=20, hjust = 0, face = 'bold')) +
      theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt"), size=12),
            axis.text.x=element_text(size=16, angle=45, vjust = 1, hjust = 1)) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
            axis.text.y=element_text(size=16, vjust = 0.5, hjust = 0, 
                                     color= rev(c("Temperature seasonality" = "#a65628", "Precipitation seasonality" = "#a65628",
                                              "Mean temperature of coldest quarter" = "#fb9a99", "Annual mean temperature" = "#fb9a99", "Mean temperature of driest quarter"= "#fb9a99",
                                              "Mean temperature of warmest quarter" = "#e31a1c", "Mean temperature of wettest quarter" = "#e31a1c",
                                              "Precipitation of wettest quarter" = "#1f78b4", "Precipitation of wettest month" = "#1f78b4", "Annual precipitation" = "#1f78b4", "Precipitation of warmest quarter" = "#1f78b4",
                                              "Precipitation of driest quarter" = "#a6cee3", "Precipitation of driest month" = "#a6cee3", "Precipitation of coldest quarter" = "#a6cee3"))))  +  
      theme(legend.position="none",
            legend.text=element_text(size=12),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(strip.background=element_blank(),
            strip.text=element_text(size=10, face="bold"),
            strip.text.y=element_text(angle=0, margin=unit(c(0, 5, 0, 0), "pt")),
            strip.placement="outside") +
      theme(plot.margin=unit(c(0.25,0.5,0.25,0.25),"cm"))  +
      guides(fill = guide_colorbar(barwidth = unit(0.6, "cm"), barheight = unit(20, "cm")))
    
    GKC_biome_climate_state_NH_0_21ka_fig 
  }
  
  
  # (2) Change rate of bioclimatic variables
  {
    # Convert variables and labels to factors for custom ordering in the plot
    GKC_biome_climate_change_rate_NH_0_21ka$Variable <- factor(GKC_biome_climate_change_rate_NH_0_21ka$Variable,
                                                         levels =c("Biome compositional turnover",  "Reconstructed megabiome shift rate", "Simulated megabiome shift rate"),
                                                         labels = c("Biome compositional turnover", "Proportion of biome shifts in reconstruction", "Proportion of biome shifts in simulation"))
    
    GKC_biome_climate_change_rate_NH_0_21ka$BIO <- factor(GKC_biome_climate_change_rate_NH_0_21ka$BIO,
                                                    levels =rev(c("Temperature seasonality", "Precipitation seasonality",
                                                                  "Mean temperature of coldest quarter", "Annual mean temperature", "Mean temperature of driest quarter",
                                                                  "Mean temperature of warmest quarter", "Mean temperature of wettest quarter",
                                                                  "Precipitation of wettest quarter", "Precipitation of wettest month", "Annual precipitation", "Precipitation of warmest quarter",
                                                                  "Precipitation of driest quarter", "Precipitation of driest month", "Precipitation of coldest quarter")))
    
    GKC_biome_climate_change_rate_NH_0_21ka_fig <- ggplot(data = GKC_biome_climate_change_rate_NH_0_21ka, aes(x = Variable, y = BIO)) +
      geom_tile(aes(fill = R_value), alpha = 0.6) +  
      geom_text(aes(label = ifelse(P_value < 0.001, paste0(sprintf("%.2f", R_value), " ***"),
                                   ifelse(P_value < 0.01, paste0(sprintf("%.2f", R_value), " **"),
                                          ifelse(P_value < 0.05, paste0(sprintf("%.2f", R_value), " *"), sprintf("%.2f", R_value)))),
                    vjust = 1, fontface = ifelse(abs(R_value) == R_value_Max, "bold", "plain")), size = 7) +
      scale_fill_gradient2(name = "R-value", low = "blue", high = "red",
                           limits = c(range(GKC_biome_climate_change_rate_NH_0_21ka$R_value)[1], range(GKC_biome_climate_change_rate_NH_0_21ka$R_value)[2]),
                           breaks = seq(-0.6, 0.8, 0.2), labels = c("-0.6", "-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8")) +
      labs(title= '(b) Change rate of bioclimatic variables', x="", y="") +
      theme_bw() +  # Reset the plot background
      theme(plot.title=element_text(size=20, hjust = 0, face = 'bold')) +
      theme(axis.title.x=element_text(margin=unit(c(0, 0, 0, 0), "pt"), size=12),
            axis.text.x=element_text(size=16, angle=45, vjust = 1, hjust = 1)) +
      theme(axis.title.y=element_text(angle=90, margin=unit(c(0, 0, 0, 0), "pt")),
            axis.text.y=element_blank())  +  
      theme(legend.position="none",
            legend.text=element_text(size=12),
            legend.title = element_text(size=14, vjust = 0.5, hjust = 0.5, face = 'bold'),
            legend.key.size = unit(1.2, "lines"),
            legend.key.height = unit(1, "lines"),
            legend.key.width = unit(2, "lines"),
            legend.background=element_rect(fill="white", size=0.5, linetype="solid", colour="transparent"),
            legend.box = "vertical",
            legend.justification = c(0.5, 0.5)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(linetype="solid", fill=NA),
            panel.background=element_rect(fill="white")) +
      theme(strip.background=element_blank(),
            strip.text=element_text(size=10, face="bold"),
            strip.text.y=element_text(angle=0, margin=unit(c(0, 5, 0, 0), "pt")),
            strip.placement="outside") +
      theme(plot.margin=unit(c(0.25,0.5,0.25,0.25),"cm"))  +
      guides(fill = guide_colorbar(barwidth = unit(0.6, "cm"), barheight = unit(20, "cm")))
    
    GKC_biome_climate_change_rate_NH_0_21ka_fig 
  }
  
    
  # -------------------------------------------------------------------------------------------------
  
  # Arrange plots
  GKC_biome_climate_state_change_rate_NH_0_21ka_fig  <- ggarrange(GKC_biome_climate_state_NH_0_21ka_fig, GKC_biome_climate_change_rate_NH_0_21ka_fig,
                                                                  nrow = 1, ncol = 2, widths = c(0.625,0.375), common.legend = F)

  # Save the final figure as a PNG file
  ggsave("result/Supplementary Fig. 5-Gaussian kernel correlation coefficients between compositional turnover and proportion of shifts in biomes with simulated bioclimatic variables and their rates of change.png", width = 15, height = 18, units = "in", dpi = 300)
  
  
}
  
