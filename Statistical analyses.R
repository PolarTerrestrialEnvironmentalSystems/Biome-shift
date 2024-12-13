
###################################
# R version 4.4.1                 #
# Operating System: Windows 10    #
# Code for statistical analysis   #
# Supplement to: Li, C., Dallmeyer, A. & Herzschuh, U. Magnitude of biome shifts over the last 21,000 years reveals the impact of seasonality and intense warming on future vegetation (2024) #
# Contact: Chenzhi Li (chenzhi.li@awi.de)[Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research, Potsdam, Germany 2024] #
###################################

# Statistical analysis list:
## 1. Statistical analysis 1: Biome compositional turnover (EMD)_21-0ka
## 2. Statistical analysis 2: Biome shifts_21-0ka
## 3. Statistical analysis 3: Simulated bioclimatic variables_21-0ka (Asia as example)
## 4. Statistical analysis 4: Simulated climate and biome over the next 300 years in the Northern Hemisphere
## 5. Statistical analysis 5: Gaussian kernel correlation

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
#setwd("~/Supplementary code/Statistical analyses")  # select the folder "Supplementary code" from your directory

# install packages if not installed
install.packages(c("paleotools", "dplyr", "zoo", "PaleoSpec", "nest"))

# loading packages
library(paleotools)
library(dplyr)
library(zoo)
library(PaleoSpec)
library(nest)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 1. Statistical analysis 1: Biome compositional turnover (EMD)_21-0ka ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# 1.1 Calculation: Earth mover's distance (EMD) between consecutive timeslices (Asia as example)
# -------------------------------------------------------------------------------------------------
# Load datasets
Pollen_percentages_timeslice_Asia_0_21ka    <- read.csv2("data/Statistical analysis 1_data.1-Pollen percentages at timeslices in Aisa_normalized.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Pollen_percentages_timeslice_Asia_0_21ka    <- type.convert(Pollen_percentages_timeslice_Asia_0_21ka, as.is = TRUE) 

# Calculation
{
  # Get the unique Dataset_IDs from the dataset
  list_ID <- unique(Pollen_percentages_timeslice_Asia_0_21ka$Dataset_ID)
  
  # Initialize an empty result container
  Pollen_EMD_sample_site_0_21ka_Asia <- NULL
  
  # Loop through each unique Dataset_ID
  for (ID in list_ID) {
    
    # Print progress message for the current Dataset_ID
    print(paste0("+++++ Dataset_ID ",ID," | (",which(list_ID == ID),"/",length(list_ID),") +++++"))
    
    # Subset the data for the current Dataset_ID
    subset_site <- subset(Pollen_percentages_timeslice_Asia_0_21ka, Dataset_ID == ID) 
    
    # Initialize a data frame to store intermediate results for the current site
    Pollen_EMD_sample_site_0_21ka <- data.frame(matrix(NA, nrow=nrow(subset_site) -1 , ncol=8, byrow=TRUE))
    colnames(Pollen_EMD_sample_site_0_21ka) <- c(names(subset_site)[1:4], "Timeslice_younger","Timeslice_older","Timeslice_span", "EMD")
    
    # Populate the first 4 columns with the repeated metadata for each time slice
    Pollen_EMD_sample_site_0_21ka[ ,1:4] <- unique(subset_site[ ,1:4])[rep(seq_len(nrow(unique(subset_site[ ,1:4]))), each=nrow(subset_site) -1), ]
    
    # Loop through the rows of the subset and calculate the Earth Mover's Distance (EMD)
    for(i in 1: (nrow(subset_site)-1)) {
      
      # Assign the younger and older timeslices and calculate the span
      Pollen_EMD_sample_site_0_21ka[i,5]  <- subset_site$Timeslice[i]
      Pollen_EMD_sample_site_0_21ka[i,6]  <- subset_site$Timeslice[i+1]
      Pollen_EMD_sample_site_0_21ka[i,7]  <- subset_site$Timeslice[i+1] - subset_site$Timeslice[i]
      
      # Compute the EMD for the current and next row's data
      Pollen_EMD_sample_site_0_21ka[i, 8]  <- paleotools::EMD(subset_site[i, 6:ncol(subset_site)], subset_site[i+1, 6:ncol(subset_site)])
      
    }
    
    # Filter results for time slices with a span of 0.5 and ensure proper data types
    Pollen_EMD_sample_site_0_21ka <- type.convert(Pollen_EMD_sample_site_0_21ka[Pollen_EMD_sample_site_0_21ka$Timeslice_span == 0.5, ], as.is = TRUE) # select two adjacent time slices
    
    # Append the filtered results to the main result container
    Pollen_EMD_sample_site_0_21ka_Asia <- rbind(Pollen_EMD_sample_site_0_21ka_Asia, Pollen_EMD_sample_site_0_21ka)
    
    # Remove rows with missing values from the combined result
    Pollen_EMD_sample_site_0_21ka_Asia <- Pollen_EMD_sample_site_0_21ka_Asia[complete.cases(Pollen_EMD_sample_site_0_21ka_Asia), ]
  }
  
  # Save the final results as a CSV file
  write.csv(Pollen_EMD_sample_site_0_21ka_Asia, file="result/Statistical analysis 1_result.1-Biome compositional turnover (EMD) between consecutive timeslices in Asia.csv", row.names=FALSE)
  
}


# 1.2 Calculation: Temporal changes of biome compositional turnover (EMD) in the Northern Hemisphere over the last 21,000 years by modern potential natural biome distributions
# -------------------------------------------------------------------------------------------------

# Load datasets
Biome_compositional_turnover_site_timeslice_NH   <- read.csv2("data/Statistical analysis 1_data.2-Biome compositional turnover (EMD) between consecutive timeslices in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_potential_natural_biomes_site_NH          <- read.csv2("data/Statistical analysis 1_data.3-Modern potential natural biomes per site in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_potential_natural_biomes_grid_NH          <- read.csv2("data/Statistical analysis 1_data.4-Modern potential natural biomes per grid-cell (spatial resolution 5 arc minutes) in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Pollen_EMD_timeslice_site_NH                     <- read.csv2("data/Statistical analysis 1_data.5-Biome compositional turnover (EMD) between consecutive timeslices per site in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Biome_compositional_turnover_site_timeslice_NH   <- type.convert(Biome_compositional_turnover_site_timeslice_NH, as.is = TRUE) 
Modern_potential_natural_biomes_site_NH          <- type.convert(Modern_potential_natural_biomes_site_NH, as.is = TRUE) 
Modern_potential_natural_biomes_grid_NH          <- type.convert(Modern_potential_natural_biomes_grid_NH, as.is = TRUE) 
Pollen_EMD_timeslice_site_NH                     <- type.convert(Pollen_EMD_timeslice_site_NH, as.is = TRUE) 

# Calculation
{
  # Calculate the proportion of modern potential natural biomes
  {
    # Set up a results data frame
    Modern_potential_natural_biomes_proportion_NH <- data.frame(matrix(NA, nrow=8, ncol=3, byrow=TRUE))
    colnames(Modern_potential_natural_biomes_proportion_NH) <- c("Biome", "Grid_number", "Proportion")
    
    # Populate the Biome column with biome IDs (1 to 8)
    Modern_potential_natural_biomes_proportion_NH$Biome        <- 1:8
    # Calculate the number of grid cells for each biome
    Modern_potential_natural_biomes_proportion_NH$Grid_number  <- table(Modern_potential_natural_biomes_grid_NH$Biome)
    # Calculate the proportion of grid cells for each biome
    Modern_potential_natural_biomes_proportion_NH$Proportion   <- prop.table(table(Modern_potential_natural_biomes_grid_NH$Biome))
  }
  
  # Extract Dataset ID lists for each biome
  {
    DatasetID_list_NH_TRFO <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 1]
    DatasetID_list_NH_WTFO <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 2]
    DatasetID_list_NH_TEFO <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 3]
    DatasetID_list_NH_BOFO <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 4]
    DatasetID_list_NH_SAVA <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 5]
    DatasetID_list_NH_STEP <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 6]
    DatasetID_list_NH_DESE <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 7]
    DatasetID_list_NH_TUND <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 8]
  }
  
  # Calculation
  {
    # TRFO
    {
      # Define a sequence of timeslices from 0 to 20.5 with a step of 0.5
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store results
      Pollen_EMD_timeslice_site_NH_TRFO_0_21ka           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Pollen_EMD_timeslice_site_NH_TRFO_0_21ka) <- c("Source", "Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "EMD")
      
      # Fill metadata for this dataset
      Pollen_EMD_timeslice_site_NH_TRFO_0_21ka$Source <- "Biome compositional turnover (EMD)"
      Pollen_EMD_timeslice_site_NH_TRFO_0_21ka$Biome <- 'TRFO'
      
      # Iterate through each timeslice to calculate EMD statistics
      for (i in 1:length(list_timeslices)) {
        
        # Print progress information
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset data for the current timeslice
        Pollen_EMD_timeslice_site_NH_TRFO_0_21ka_subset <- subset(Pollen_EMD_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_TRFO & Timeslice_younger == list_timeslices[i])
        
        # Populate the results data frame based on the subset
        if(nrow(Pollen_EMD_timeslice_site_NH_TRFO_0_21ka_subset) > 0) {
          
          Pollen_EMD_timeslice_site_NH_TRFO_0_21ka[i,3]    <- unique(Pollen_EMD_timeslice_site_NH_TRFO_0_21ka_subset$Timeslice_younger)
          Pollen_EMD_timeslice_site_NH_TRFO_0_21ka[i,4]    <- unique(Pollen_EMD_timeslice_site_NH_TRFO_0_21ka_subset$Timeslice_older)
          Pollen_EMD_timeslice_site_NH_TRFO_0_21ka[i,5]    <- length(unique(Pollen_EMD_timeslice_site_NH_TRFO_0_21ka_subset$Dataset_ID))
          Pollen_EMD_timeslice_site_NH_TRFO_0_21ka[i,6]    <- median(Pollen_EMD_timeslice_site_NH_TRFO_0_21ka_subset$EMD)
          
        } else{
          # Handle cases where there is no data for the current timeslice
          Pollen_EMD_timeslice_site_NH_TRFO_0_21ka[i,3]    <- list_timeslices[i]
          Pollen_EMD_timeslice_site_NH_TRFO_0_21ka[i,4]    <- list_timeslices[i] + 0.5
          Pollen_EMD_timeslice_site_NH_TRFO_0_21ka[i,5]    <- 0
          Pollen_EMD_timeslice_site_NH_TRFO_0_21ka[i,6]    <- NA
        }
        
      }
      
    } 
    
    # WTFO
    {
      # Define a sequence of timeslices from 0 to 20.5 with a step of 0.5
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store results
      Pollen_EMD_timeslice_site_NH_WTFO_0_21ka           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Pollen_EMD_timeslice_site_NH_WTFO_0_21ka) <- c("Source", "Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "EMD")
      
      # Fill metadata for this dataset
      Pollen_EMD_timeslice_site_NH_WTFO_0_21ka$Source <- "Biome compositional turnover (EMD)"
      Pollen_EMD_timeslice_site_NH_WTFO_0_21ka$Biome <- 'WTFO'
      
      # Iterate through each timeslice to calculate EMD statistics
      for (i in 1:length(list_timeslices)) {
        
        # Print progress information
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset data for the current timeslice
        Pollen_EMD_timeslice_site_NH_WTFO_0_21ka_subset <- subset(Pollen_EMD_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_WTFO & Timeslice_younger == list_timeslices[i])
        
        # Populate the results data frame based on the subset
        if(nrow(Pollen_EMD_timeslice_site_NH_WTFO_0_21ka_subset) > 0) {
          
          Pollen_EMD_timeslice_site_NH_WTFO_0_21ka[i,3]    <- unique(Pollen_EMD_timeslice_site_NH_WTFO_0_21ka_subset$Timeslice_younger)
          Pollen_EMD_timeslice_site_NH_WTFO_0_21ka[i,4]    <- unique(Pollen_EMD_timeslice_site_NH_WTFO_0_21ka_subset$Timeslice_older)
          Pollen_EMD_timeslice_site_NH_WTFO_0_21ka[i,5]    <- length(unique(Pollen_EMD_timeslice_site_NH_WTFO_0_21ka_subset$Dataset_ID))
          Pollen_EMD_timeslice_site_NH_WTFO_0_21ka[i,6]    <- median(Pollen_EMD_timeslice_site_NH_WTFO_0_21ka_subset$EMD)
          
        } else{
          # Handle cases where there is no data for the current timeslice
          Pollen_EMD_timeslice_site_NH_WTFO_0_21ka[i,3]    <- list_timeslices[i]
          Pollen_EMD_timeslice_site_NH_WTFO_0_21ka[i,4]    <- list_timeslices[i] + 0.5
          Pollen_EMD_timeslice_site_NH_WTFO_0_21ka[i,5]    <- 0
          Pollen_EMD_timeslice_site_NH_WTFO_0_21ka[i,6]    <- NA
        }
        
      }
      
    } 
    
    # TEFO
    {
      # Define a sequence of timeslices from 0 to 20.5 with a step of 0.5
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store results
      Pollen_EMD_timeslice_site_NH_TEFO_0_21ka           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Pollen_EMD_timeslice_site_NH_TEFO_0_21ka) <- c("Source", "Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "EMD")
      
      # Fill metadata for this dataset
      Pollen_EMD_timeslice_site_NH_TEFO_0_21ka$Source <- "Biome compositional turnover (EMD)"
      Pollen_EMD_timeslice_site_NH_TEFO_0_21ka$Biome <- 'TEFO'
      
      # Iterate through each timeslice to calculate EMD statistics
      for (i in 1:length(list_timeslices)) {
        
        # Print progress information
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset data for the current timeslice
        Pollen_EMD_timeslice_site_NH_TEFO_0_21ka_subset <- subset(Pollen_EMD_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_TEFO & Timeslice_younger == list_timeslices[i])
        
        # Populate the results data frame based on the subset
        if(nrow(Pollen_EMD_timeslice_site_NH_TEFO_0_21ka_subset) > 0) {
          
          Pollen_EMD_timeslice_site_NH_TEFO_0_21ka[i,3]    <- unique(Pollen_EMD_timeslice_site_NH_TEFO_0_21ka_subset$Timeslice_younger)
          Pollen_EMD_timeslice_site_NH_TEFO_0_21ka[i,4]    <- unique(Pollen_EMD_timeslice_site_NH_TEFO_0_21ka_subset$Timeslice_older)
          Pollen_EMD_timeslice_site_NH_TEFO_0_21ka[i,5]    <- length(unique(Pollen_EMD_timeslice_site_NH_TEFO_0_21ka_subset$Dataset_ID))
          Pollen_EMD_timeslice_site_NH_TEFO_0_21ka[i,6]    <- median(Pollen_EMD_timeslice_site_NH_TEFO_0_21ka_subset$EMD)
          
        } else{
          # Handle cases where there is no data for the current timeslice
          Pollen_EMD_timeslice_site_NH_TEFO_0_21ka[i,3]    <- list_timeslices[i]
          Pollen_EMD_timeslice_site_NH_TEFO_0_21ka[i,4]    <- list_timeslices[i] + 0.5
          Pollen_EMD_timeslice_site_NH_TEFO_0_21ka[i,5]    <- 0
          Pollen_EMD_timeslice_site_NH_TEFO_0_21ka[i,6]    <- NA
        }
        
      }
      
    } 
    
    # BOFO
    {
      # Define a sequence of timeslices from 0 to 20.5 with a step of 0.5
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store results
      Pollen_EMD_timeslice_site_NH_BOFO_0_21ka           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Pollen_EMD_timeslice_site_NH_BOFO_0_21ka) <- c("Source", "Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "EMD")
      
      # Fill metadata for this dataset
      Pollen_EMD_timeslice_site_NH_BOFO_0_21ka$Source <- "Biome compositional turnover (EMD)"
      Pollen_EMD_timeslice_site_NH_BOFO_0_21ka$Biome <- 'BOFO'
      
      # Iterate through each timeslice to calculate EMD statistics
      for (i in 1:length(list_timeslices)) {
        
        # Print progress information
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset data for the current timeslice
        Pollen_EMD_timeslice_site_NH_BOFO_0_21ka_subset <- subset(Pollen_EMD_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_BOFO & Timeslice_younger == list_timeslices[i])
        
        # Populate the results data frame based on the subset
        if(nrow(Pollen_EMD_timeslice_site_NH_BOFO_0_21ka_subset) > 0) {
          
          Pollen_EMD_timeslice_site_NH_BOFO_0_21ka[i,3]    <- unique(Pollen_EMD_timeslice_site_NH_BOFO_0_21ka_subset$Timeslice_younger)
          Pollen_EMD_timeslice_site_NH_BOFO_0_21ka[i,4]    <- unique(Pollen_EMD_timeslice_site_NH_BOFO_0_21ka_subset$Timeslice_older)
          Pollen_EMD_timeslice_site_NH_BOFO_0_21ka[i,5]    <- length(unique(Pollen_EMD_timeslice_site_NH_BOFO_0_21ka_subset$Dataset_ID))
          Pollen_EMD_timeslice_site_NH_BOFO_0_21ka[i,6]    <- median(Pollen_EMD_timeslice_site_NH_BOFO_0_21ka_subset$EMD)
          
        } else{
          # Handle cases where there is no data for the current timeslice
          Pollen_EMD_timeslice_site_NH_BOFO_0_21ka[i,3]    <- list_timeslices[i]
          Pollen_EMD_timeslice_site_NH_BOFO_0_21ka[i,4]    <- list_timeslices[i] + 0.5
          Pollen_EMD_timeslice_site_NH_BOFO_0_21ka[i,5]    <- 0
          Pollen_EMD_timeslice_site_NH_BOFO_0_21ka[i,6]    <- NA
        }
        
      }
      
    } 
    
    # SAVA
    {
      # Define a sequence of timeslices from 0 to 20.5 with a step of 0.5
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store results
      Pollen_EMD_timeslice_site_NH_SAVA_0_21ka           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Pollen_EMD_timeslice_site_NH_SAVA_0_21ka) <- c("Source", "Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "EMD")
      
      # Fill metadata for this dataset
      Pollen_EMD_timeslice_site_NH_SAVA_0_21ka$Source <- "Biome compositional turnover (EMD)"
      Pollen_EMD_timeslice_site_NH_SAVA_0_21ka$Biome <- 'SAVA'
      
      # Iterate through each timeslice to calculate EMD statistics
      for (i in 1:length(list_timeslices)) {
        
        # Print progress information
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset data for the current timeslice
        Pollen_EMD_timeslice_site_NH_SAVA_0_21ka_subset <- subset(Pollen_EMD_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_SAVA & Timeslice_younger == list_timeslices[i])
        
        # Populate the results data frame based on the subset
        if(nrow(Pollen_EMD_timeslice_site_NH_SAVA_0_21ka_subset) > 0) {
          
          Pollen_EMD_timeslice_site_NH_SAVA_0_21ka[i,3]    <- unique(Pollen_EMD_timeslice_site_NH_SAVA_0_21ka_subset$Timeslice_younger)
          Pollen_EMD_timeslice_site_NH_SAVA_0_21ka[i,4]    <- unique(Pollen_EMD_timeslice_site_NH_SAVA_0_21ka_subset$Timeslice_older)
          Pollen_EMD_timeslice_site_NH_SAVA_0_21ka[i,5]    <- length(unique(Pollen_EMD_timeslice_site_NH_SAVA_0_21ka_subset$Dataset_ID))
          Pollen_EMD_timeslice_site_NH_SAVA_0_21ka[i,6]    <- median(Pollen_EMD_timeslice_site_NH_SAVA_0_21ka_subset$EMD)
          
        } else{
          # Handle cases where there is no data for the current timeslice
          Pollen_EMD_timeslice_site_NH_SAVA_0_21ka[i,3]    <- list_timeslices[i]
          Pollen_EMD_timeslice_site_NH_SAVA_0_21ka[i,4]    <- list_timeslices[i] + 0.5
          Pollen_EMD_timeslice_site_NH_SAVA_0_21ka[i,5]    <- 0
          Pollen_EMD_timeslice_site_NH_SAVA_0_21ka[i,6]    <- NA
        }
        
      }
      
    } 
    
    # STEP
    {
      # Define a sequence of timeslices from 0 to 20.5 with a step of 0.5
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store results
      Pollen_EMD_timeslice_site_NH_STEP_0_21ka           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Pollen_EMD_timeslice_site_NH_STEP_0_21ka) <- c("Source", "Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "EMD")
      
      # Fill metadata for this dataset
      Pollen_EMD_timeslice_site_NH_STEP_0_21ka$Source <- "Biome compositional turnover (EMD)"
      Pollen_EMD_timeslice_site_NH_STEP_0_21ka$Biome <- 'STEP'
      
      # Iterate through each timeslice to calculate EMD statistics
      for (i in 1:length(list_timeslices)) {
        
        # Print progress information
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset data for the current timeslice
        Pollen_EMD_timeslice_site_NH_STEP_0_21ka_subset <- subset(Pollen_EMD_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_STEP & Timeslice_younger == list_timeslices[i])
        
        # Populate the results data frame based on the subset
        if(nrow(Pollen_EMD_timeslice_site_NH_STEP_0_21ka_subset) > 0) {
          
          Pollen_EMD_timeslice_site_NH_STEP_0_21ka[i,3]    <- unique(Pollen_EMD_timeslice_site_NH_STEP_0_21ka_subset$Timeslice_younger)
          Pollen_EMD_timeslice_site_NH_STEP_0_21ka[i,4]    <- unique(Pollen_EMD_timeslice_site_NH_STEP_0_21ka_subset$Timeslice_older)
          Pollen_EMD_timeslice_site_NH_STEP_0_21ka[i,5]    <- length(unique(Pollen_EMD_timeslice_site_NH_STEP_0_21ka_subset$Dataset_ID))
          Pollen_EMD_timeslice_site_NH_STEP_0_21ka[i,6]    <- median(Pollen_EMD_timeslice_site_NH_STEP_0_21ka_subset$EMD)
          
        } else{
          # Handle cases where there is no data for the current timeslice
          Pollen_EMD_timeslice_site_NH_STEP_0_21ka[i,3]    <- list_timeslices[i]
          Pollen_EMD_timeslice_site_NH_STEP_0_21ka[i,4]    <- list_timeslices[i] + 0.5
          Pollen_EMD_timeslice_site_NH_STEP_0_21ka[i,5]    <- 0
          Pollen_EMD_timeslice_site_NH_STEP_0_21ka[i,6]    <- NA
        }
        
      }
      
    } 
    
    # DESE
    {
      # Define a sequence of timeslices from 0 to 20.5 with a step of 0.5
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store results
      Pollen_EMD_timeslice_site_NH_DESE_0_21ka           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Pollen_EMD_timeslice_site_NH_DESE_0_21ka) <- c("Source", "Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "EMD")
      
      # Fill metadata for this dataset
      Pollen_EMD_timeslice_site_NH_DESE_0_21ka$Source <- "Biome compositional turnover (EMD)"
      Pollen_EMD_timeslice_site_NH_DESE_0_21ka$Biome <- 'DESE'
      
      # Iterate through each timeslice to calculate EMD statistics
      for (i in 1:length(list_timeslices)) {
        
        # Print progress information
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset data for the current timeslice
        Pollen_EMD_timeslice_site_NH_DESE_0_21ka_subset <- subset(Pollen_EMD_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_DESE & Timeslice_younger == list_timeslices[i])
        
        # Populate the results data frame based on the subset
        if(nrow(Pollen_EMD_timeslice_site_NH_DESE_0_21ka_subset) > 0) {
          
          Pollen_EMD_timeslice_site_NH_DESE_0_21ka[i,3]    <- unique(Pollen_EMD_timeslice_site_NH_DESE_0_21ka_subset$Timeslice_younger)
          Pollen_EMD_timeslice_site_NH_DESE_0_21ka[i,4]    <- unique(Pollen_EMD_timeslice_site_NH_DESE_0_21ka_subset$Timeslice_older)
          Pollen_EMD_timeslice_site_NH_DESE_0_21ka[i,5]    <- length(unique(Pollen_EMD_timeslice_site_NH_DESE_0_21ka_subset$Dataset_ID))
          Pollen_EMD_timeslice_site_NH_DESE_0_21ka[i,6]    <- median(Pollen_EMD_timeslice_site_NH_DESE_0_21ka_subset$EMD)
          
        } else{
          # Handle cases where there is no data for the current timeslice
          Pollen_EMD_timeslice_site_NH_DESE_0_21ka[i,3]    <- list_timeslices[i]
          Pollen_EMD_timeslice_site_NH_DESE_0_21ka[i,4]    <- list_timeslices[i] + 0.5
          Pollen_EMD_timeslice_site_NH_DESE_0_21ka[i,5]    <- 0
          Pollen_EMD_timeslice_site_NH_DESE_0_21ka[i,6]    <- NA
        }
        
      }
      
    } 
    
    # TUND
    {
      # Define a sequence of timeslices from 0 to 20.5 with a step of 0.5
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store results
      Pollen_EMD_timeslice_site_NH_TUND_0_21ka           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Pollen_EMD_timeslice_site_NH_TUND_0_21ka) <- c("Source", "Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "EMD")
      
      # Fill metadata for this dataset
      Pollen_EMD_timeslice_site_NH_TUND_0_21ka$Source <- "Biome compositional turnover (EMD)"
      Pollen_EMD_timeslice_site_NH_TUND_0_21ka$Biome <- 'TUND'
      
      # Iterate through each timeslice to calculate EMD statistics
      for (i in 1:length(list_timeslices)) {
        
        # Print progress information
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset data for the current timeslice
        Pollen_EMD_timeslice_site_NH_TUND_0_21ka_subset <- subset(Pollen_EMD_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_TUND & Timeslice_younger == list_timeslices[i])
        
        # Populate the results data frame based on the subset
        if(nrow(Pollen_EMD_timeslice_site_NH_TUND_0_21ka_subset) > 0) {
          
          Pollen_EMD_timeslice_site_NH_TUND_0_21ka[i,3]    <- unique(Pollen_EMD_timeslice_site_NH_TUND_0_21ka_subset$Timeslice_younger)
          Pollen_EMD_timeslice_site_NH_TUND_0_21ka[i,4]    <- unique(Pollen_EMD_timeslice_site_NH_TUND_0_21ka_subset$Timeslice_older)
          Pollen_EMD_timeslice_site_NH_TUND_0_21ka[i,5]    <- length(unique(Pollen_EMD_timeslice_site_NH_TUND_0_21ka_subset$Dataset_ID))
          Pollen_EMD_timeslice_site_NH_TUND_0_21ka[i,6]    <- median(Pollen_EMD_timeslice_site_NH_TUND_0_21ka_subset$EMD)
          
        } else{
          # Handle cases where there is no data for the current timeslice
          Pollen_EMD_timeslice_site_NH_TUND_0_21ka[i,3]    <- list_timeslices[i]
          Pollen_EMD_timeslice_site_NH_TUND_0_21ka[i,4]    <- list_timeslices[i] + 0.5
          Pollen_EMD_timeslice_site_NH_TUND_0_21ka[i,5]    <- 0
          Pollen_EMD_timeslice_site_NH_TUND_0_21ka[i,6]    <- NA
        }
        
      }
      
    } 
    
    # NH
    {
      # Combine all biome-specific data frames into one for the Northern Hemisphere
      Pollen_EMD_timeslice_site_NH_all_0_21ka  <- rbind(Pollen_EMD_timeslice_site_NH_TRFO_0_21ka, Pollen_EMD_timeslice_site_NH_WTFO_0_21ka, 
                                                        Pollen_EMD_timeslice_site_NH_TEFO_0_21ka, Pollen_EMD_timeslice_site_NH_BOFO_0_21ka, 
                                                        Pollen_EMD_timeslice_site_NH_SAVA_0_21ka, Pollen_EMD_timeslice_site_NH_STEP_0_21ka, 
                                                        Pollen_EMD_timeslice_site_NH_DESE_0_21ka, Pollen_EMD_timeslice_site_NH_TUND_0_21ka)
      
      # Define a sequence of timeslices from 0 to 20.5 with a step of 0.5
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store results
      Pollen_EMD_timeslice_site_NH_NH_0_21ka           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Pollen_EMD_timeslice_site_NH_NH_0_21ka) <- c("Source", "Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "EMD")
      
      # Fill metadata for this dataset
      Pollen_EMD_timeslice_site_NH_NH_0_21ka$Source <- "Biome compositional turnover (EMD)"
      Pollen_EMD_timeslice_site_NH_NH_0_21ka$Biome <- 'Northern Hemisphere'
      
      # Iterate through each timeslice to calculate EMD statistics
      for (i in 1:length(list_timeslices)) {
        
        # Print progress information
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset data for the current timeslice
        Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset <- subset(Pollen_EMD_timeslice_site_NH_all_0_21ka, Timeslice_younger == list_timeslices[i])
        # Replace NA values with 0 in the subset to handle missing data
        Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset[is.na(Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset)] <- 0
        
        # Assign the timeslice bounds to the result data frame
        Pollen_EMD_timeslice_site_NH_NH_0_21ka[i,3]    <- unique(Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$Timeslice_younger)
        Pollen_EMD_timeslice_site_NH_NH_0_21ka[i,4]    <- unique(Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$Timeslice_older)
        # Sum the number of sites across all biomes for the current timeslice
        Pollen_EMD_timeslice_site_NH_NH_0_21ka[i,5]    <- sum(Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$Site_number)
        
        # Compute the EMD (Earth Mover's Distance) weighted by the proportion of each biome
        Pollen_EMD_timeslice_site_NH_NH_0_21ka[i,6]    <- (Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$EMD[Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$Biome == 'TRFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==1] + 
                                                             Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$EMD[Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$Biome == 'WTFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==2] +
                                                             Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$EMD[Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$Biome == 'TEFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==3] + 
                                                             Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$EMD[Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$Biome == 'BOFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==4] +
                                                             Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$EMD[Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$Biome == 'SAVA']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==5] +
                                                             Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$EMD[Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$Biome == 'STEP']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==6] +
                                                             Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$EMD[Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$Biome == 'DESE']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==7] +
                                                             Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$EMD[Pollen_EMD_timeslice_site_NH_NH_0_21ka_subset$Biome == 'TUND']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==8] )
        
      }
      
    }
    
  }
  
  # Combine result dataset
  Pollen_EMD_timeslice_site_NH_0_21ka_final <- rbind(Pollen_EMD_timeslice_site_NH_NH_0_21ka, Pollen_EMD_timeslice_site_NH_all_0_21ka)
  
  # Save the final results as a CSV file
  write.csv(Pollen_EMD_timeslice_site_NH_0_21ka_final, file="result/Statistical analysis 1_result.2-Temporal changes of biome compositional turnover (EMD) in the Northern Hemisphere over the last 21,000 years by modern biome distributions.csv", row.names=FALSE)
  
}


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 2. Statistical analysis 2: Biome shifts_21-0ka ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# 2.1 Calculation: Proportion and state of biome shift per site_Global
# -------------------------------------------------------------------------------------------------
# Load datasets
Biome_timeslice_pollen_MPI_site_global    <- read.csv2("data/Statistical analysis 2_data.1-Global reconstructed and simulated biome at timeslices per site.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Biome_timeslice_pollen_MPI_site_global    <- type.convert(Biome_timeslice_pollen_MPI_site_global, as.is = TRUE) 

# Calculation
{
  # (1) Biome shift state at timeslice per site
  
  {
    # Pollen-based reconstruction
    {
      # Extract a list of unique dataset IDs from the global dataset
      list_ID <- unique(Biome_timeslice_pollen_MPI_site_global$Dataset_ID)
      
      # Initialize an empty data frame to store results for all datasets
      Biome_shift_state_pollen <- NULL
      
      # Loop through each dataset ID
      for (ID in list_ID) {
        # Print progress information
        print(paste0("+++++ Dataset_ID ",ID," | (",which(list_ID == ID),"/",length(list_ID),") +++++"))
        
        # Subset the global dataset for the current Dataset_ID
        subset_site <- subset(Biome_timeslice_pollen_MPI_site_global, Dataset_ID == ID) 
        
        # Initialize a data frame to compare biomes between time slices for the current dataset
        biome_compared <- data.frame(matrix(NA, nrow=nrow(subset_site)-1, ncol=10, byrow=TRUE))
        colnames(biome_compared) <- c("Dataset_ID", "Longitude", "Latitude", "Timeslice_younger", "Biome_younger", "Timeslice_older", "Biome_older", "Timeslice_span", "Shifting", "Source")
        
        # Populate constant metadata fields
        biome_compared$Dataset_ID <- unique(subset_site$Dataset_ID)
        biome_compared$Longitude  <- unique(subset_site$Longitude)
        biome_compared$Latitude   <- unique(subset_site$Latitude)
        biome_compared$Source     <- "Pollen-based reconstruction"
        
        # Compare biomes for each pair of consecutive time slices
        for(i in 1:c(nrow(subset_site)-1)){ 
          # Assign younger time slice and biome
          biome_compared[i,4] <- subset_site[i,7] # Timeslice_younger
          biome_compared[i,5] <- subset_site[i,8] # Biome_younger
          
          # Assign older time slice and biome
          biome_compared[i,6] <- subset_site[i+1,7] # Timeslice_older
          biome_compared[i,7] <- subset_site[i+1,8] # Biome_older
          
          # Calculate the time span between slices
          biome_compared[i,8] <- (subset_site[i+1,7] - subset_site[i,7]) # Timeslice_span
          
          # Determine if there is a biome shift
          biome_compared[i,9] <- ifelse(subset_site[i,8] == subset_site[i+1,8], "no", "yes")
        }
        
        # Retain only rows where the timeslice span equals 0.5 (adjacent time slices)
        biome_compared <- type.convert(biome_compared[biome_compared$Timeslice_span == 0.5, ], as.is = TRUE) # select two adjacent time slices
        
        # Append the results for the current dataset to the final results data frame
        Biome_shift_state_pollen <- rbind(Biome_shift_state_pollen, biome_compared)
        
      }
      
    }
    
    
    # ESM-based simulation
    {
      # Extract a list of unique dataset IDs from the global dataset
      list_ID <- unique(Biome_timeslice_pollen_MPI_site_global$Dataset_ID)
      
      # Initialize an empty data frame to store results for all datasets
      Biome_shift_state_MPI <- NULL
      
      # Loop through each dataset ID
      for (ID in list_ID) {
        # Print progress information
        print(paste0("+++++ Dataset_ID ",ID," | (",which(list_ID == ID),"/",length(list_ID),") +++++"))
        
        # Subset the global dataset for the current Dataset_ID
        subset_site <- subset(Biome_timeslice_pollen_MPI_site_global, Dataset_ID == ID) 
        
        # Initialize a data frame to compare biomes between time slices for the current dataset
        biome_compared <- data.frame(matrix(NA, nrow=nrow(subset_site)-1, ncol=10, byrow=TRUE))
        colnames(biome_compared) <- c("Dataset_ID", "Longitude", "Latitude", "Timeslice_younger", "Biome_younger", "Timeslice_older", "Biome_older", "Timeslice_span", "Shifting", "Source")
        
        # Populate constant metadata fields
        biome_compared$Dataset_ID <- unique(subset_site$Dataset_ID)
        biome_compared$Longitude  <- unique(subset_site$Longitude)
        biome_compared$Latitude   <- unique(subset_site$Latitude)
        biome_compared$Source     <- "ESM-based simulation"
        
        # Compare biomes for each pair of consecutive time slices
        for(i in 1:c(nrow(subset_site)-1)){ 
          # Assign younger time slice and biome
          biome_compared[i,4] <- subset_site[i,7] # Timeslice_younger
          biome_compared[i,5] <- subset_site[i,9] # Biome_younger
          
          # Assign older time slice and biome
          biome_compared[i,6] <- subset_site[i+1,7] # Timeslice_older
          biome_compared[i,7] <- subset_site[i+1,9] # Biome_older
          
          # Calculate the time span between slices
          biome_compared[i,8] <- (subset_site[i+1,7] - subset_site[i,7]) # Timeslice_span
          
          # Determine if there is a biome shift
          biome_compared[i,9] <- ifelse(subset_site[i,9] == subset_site[i+1,9], "no", "yes")
        }
        
        # Retain only rows where the timeslice span equals 0.5 (adjacent time slices)
        biome_compared <- type.convert(biome_compared[biome_compared$Timeslice_span == 0.5, ], as.is = TRUE) # select two adjacent time slices
        
        # Append the results for the current dataset to the final results data frame
        Biome_shift_state_MPI <- rbind(Biome_shift_state_MPI, biome_compared)
        
      }
      
    }
    
    
    # Combine result dataset
    Biome_shift_state_pollen_MPI_global <- rbind(Biome_shift_state_pollen, Biome_shift_state_MPI)
    
    # Save the final results as a CSV file
    write.csv(Biome_shift_state_pollen_MPI_global, file="result/Statistical analysis 2_result.1-Biome shift state at timeslices per site in reconstruction and simulation globally.csv", row.names=FALSE)
    
  }
  
  
  # -------------------------------------------------------------------------------------------------
  
  # (2) Proportion of biome shifts per site over the last 21,000 years
  {
    # Pollen-based reconstruction
    {
      # Extract a list of unique dataset IDs from dataset
      list_ID <- unique(Biome_shift_state_pollen$Dataset_ID)
      
      # Initialize a data frame to store biome shift proportions for each dataset
      Biome_shift_proportion_pollen           <- data.frame(matrix(NA, nrow=length(list_ID), ncol=5, byrow=TRUE))
      colnames(Biome_shift_proportion_pollen) <- c("Dataset_ID", "Longitude", "Latitude", "Proportion", "Source")
      
      # Populate the Dataset_ID column with the unique IDs and set the source description
      Biome_shift_proportion_pollen$Dataset_ID <- list_ID
      Biome_shift_proportion_pollen$Source     <- "Pollen-based reconstruction"
      
      # Loop through each dataset ID to calculate and store the biome shift proportion
      for (ID in list_ID) {
        
        # Print progress message for the current dataset being processed
        print(paste0("+++++ SITE ",ID," | (",which(list_ID == ID),"/",length(list_ID),") +++++"))
        
        # Subset the data to include only records belonging to the current dataset ID
        subset_site <- subset(Biome_shift_state_pollen, Dataset_ID == ID) 
        
        # Retrieve and store the longitude and latitude for the current site
        Biome_shift_proportion_pollen[which(Biome_shift_proportion_pollen$Dataset_ID == ID), 2] <- unique(subset_site$Longitude)
        Biome_shift_proportion_pollen[which(Biome_shift_proportion_pollen$Dataset_ID == ID), 3] <- unique(subset_site$Latitude)
        
        # Calculate the proportion of shifting biomes for the current dataset
        Biome_shift_proportion_pollen[which(Biome_shift_proportion_pollen$Dataset_ID == ID), 4] <- sum(subset_site$Shifting == "yes")/nrow(subset_site)
        
      }
      
      # Filter the resulting data frame to include only rows with complete data
      Biome_shift_proportion_pollen <- Biome_shift_proportion_pollen[complete.cases(Biome_shift_proportion_pollen), ]
      
    }
    
    # ESM-based simulation
    {
      # Extract a list of unique dataset IDs from dataset
      list_ID <- unique(Biome_shift_state_MPI$Dataset_ID)
      
      # Initialize a data frame to store biome shift proportions for each dataset
      Biome_shift_proportion_MPI           <- data.frame(matrix(NA, nrow=length(list_ID), ncol=5, byrow=TRUE))
      colnames(Biome_shift_proportion_MPI) <- c("Dataset_ID", "Longitude", "Latitude", "Proportion", "Source")
      
      # Populate the Dataset_ID column with the unique IDs and set the source description
      Biome_shift_proportion_MPI$Dataset_ID <- list_ID
      Biome_shift_proportion_MPI$Source     <- "MPI-based reconstruction"
      
      # Loop through each dataset ID to calculate and store the biome shift proportion
      for (ID in list_ID) {
        
        # Print progress message for the current dataset being processed
        print(paste0("+++++ SITE ",ID," | (",which(list_ID == ID),"/",length(list_ID),") +++++"))
        
        # Subset the data to include only records belonging to the current dataset ID
        subset_site <- subset(Biome_shift_state_MPI, Dataset_ID == ID) 
        
        # Retrieve and store the longitude and latitude for the current site
        Biome_shift_proportion_MPI[which(Biome_shift_proportion_MPI$Dataset_ID == ID), 2] <- unique(subset_site$Longitude)
        Biome_shift_proportion_MPI[which(Biome_shift_proportion_MPI$Dataset_ID == ID), 3] <- unique(subset_site$Latitude)
        
        # Calculate the proportion of shifting biomes for the current dataset
        Biome_shift_proportion_MPI[which(Biome_shift_proportion_MPI$Dataset_ID == ID), 4] <- sum(subset_site$Shifting == "yes")/nrow(subset_site)
        
      }
      
      # Filter the resulting data frame to include only rows with complete data
      Biome_shift_proportion_MPI <- Biome_shift_proportion_MPI[complete.cases(Biome_shift_proportion_MPI), ]
      
    }
    
    
    # Combine result dataset
    Biome_shift_proportion_pollen_MPI_global <- rbind(Biome_shift_proportion_pollen, Biome_shift_proportion_MPI)
    
    # Save the final results as a CSV file
    write.csv(Biome_shift_proportion_pollen_MPI_global, file="result/Statistical analysis 2_result.2-Proportion of biome shifts per site in reconstruction and simulation over the last 21,000 years globally.csv", row.names=FALSE)
    
  }
  
}


# 2.2 Calculation: Temporal changes of proportion of biome shifts in the Northern Hemisphere over the last 21,000 years by modern potential natural biome distributions
# -------------------------------------------------------------------------------------------------

# Load datasets
Biome_shift_state_pollen_MPI_global        <- read.csv2("data/Statistical analysis 2_data.2-Biome shift state at timeslices per site in reconstruction and simulation globally.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_potential_natural_biomes_site_NH    <- read.csv2("data/Statistical analysis 2_data.3-Modern potential natural biomes per site in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_potential_natural_biomes_grid_NH    <- read.csv2("data/Statistical analysis 2_data.4-Modern potential natural biomes per grid-cell (spatial resolution 5 arc minutes) in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Biome_shift_state_pollen_MPI_global        <- type.convert(Biome_shift_state_pollen_MPI_global, as.is = TRUE) 
Modern_potential_natural_biomes_site_NH    <- type.convert(Modern_potential_natural_biomes_site_NH, as.is = TRUE) 
Modern_potential_natural_biomes_grid_NH    <- type.convert(Modern_potential_natural_biomes_grid_NH, as.is = TRUE) 

# Calculation
{
  # Calculate the proportion of modern potential natural biomes
  {
    # Set up a results data frame
    Modern_potential_natural_biomes_proportion_NH <- data.frame(matrix(NA, nrow=8, ncol=3, byrow=TRUE))
    colnames(Modern_potential_natural_biomes_proportion_NH) <- c("Biome", "Grid_number", "Proportion")
    
    # Populate the Biome column with biome IDs (1 to 8)
    Modern_potential_natural_biomes_proportion_NH$Biome        <- 1:8
    # Calculate the number of grid cells for each biome
    Modern_potential_natural_biomes_proportion_NH$Grid_number  <- table(Modern_potential_natural_biomes_grid_NH$Biome)
    # Calculate the proportion of grid cells for each biome
    Modern_potential_natural_biomes_proportion_NH$Proportion   <- prop.table(table(Modern_potential_natural_biomes_grid_NH$Biome))
  }
  
  # Extract Dataset ID lists for each biome
  {
    DatasetID_list_NH_TRFO <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 1]
    DatasetID_list_NH_WTFO <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 2]
    DatasetID_list_NH_TEFO <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 3]
    DatasetID_list_NH_BOFO <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 4]
    DatasetID_list_NH_SAVA <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 5]
    DatasetID_list_NH_STEP <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 6]
    DatasetID_list_NH_DESE <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 7]
    DatasetID_list_NH_TUND <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 8]
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # (1) Pollen-based reconstruction
  {
    # TRFO
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_pollen_NH_TRFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_pollen_NH_TRFO) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_pollen_NH_TRFO[ ,1] <- 'TRFO'
      Biome_shift_proportion_temporal_pollen_NH_TRFO$Source <- "Pollen-based reconstruction"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_pollen_NH_TRFO_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "Pollen-based reconstruction" & Dataset_ID %in% DatasetID_list_NH_TRFO & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_pollen_NH_TRFO_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_pollen_NH_TRFO[i,2]    <- unique(Biome_shift_state_pollen_NH_TRFO_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_pollen_NH_TRFO[i,3]    <- unique(Biome_shift_state_pollen_NH_TRFO_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_pollen_NH_TRFO[i,4]    <- length(unique(Biome_shift_state_pollen_NH_TRFO_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_pollen_NH_TRFO[i,5]    <- sum(Biome_shift_state_pollen_NH_TRFO_subset$Shifting == "yes")/length(unique(Biome_shift_state_pollen_NH_TRFO_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_pollen_NH_TRFO[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_pollen_NH_TRFO[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_pollen_NH_TRFO[i,4]    <- 0
          Biome_shift_proportion_temporal_pollen_NH_TRFO[i,5]    <- NA
        }
        
      }
      
    }
    
    # WTFO
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_pollen_NH_WTFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_pollen_NH_WTFO) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_pollen_NH_WTFO[ ,1] <- 'WTFO'
      Biome_shift_proportion_temporal_pollen_NH_WTFO$Source <- "Pollen-based reconstruction"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_pollen_NH_WTFO_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "Pollen-based reconstruction" & Dataset_ID %in% DatasetID_list_NH_WTFO & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_pollen_NH_WTFO_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_pollen_NH_WTFO[i,2]    <- unique(Biome_shift_state_pollen_NH_WTFO_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_pollen_NH_WTFO[i,3]    <- unique(Biome_shift_state_pollen_NH_WTFO_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_pollen_NH_WTFO[i,4]    <- length(unique(Biome_shift_state_pollen_NH_WTFO_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_pollen_NH_WTFO[i,5]    <- sum(Biome_shift_state_pollen_NH_WTFO_subset$Shifting == "yes")/length(unique(Biome_shift_state_pollen_NH_WTFO_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_pollen_NH_WTFO[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_pollen_NH_WTFO[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_pollen_NH_WTFO[i,4]    <- 0
          Biome_shift_proportion_temporal_pollen_NH_WTFO[i,5]    <- NA
        }
        
      }
      
    }
    
    # TEFO
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_pollen_NH_TEFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_pollen_NH_TEFO) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_pollen_NH_TEFO[ ,1] <- 'TEFO'
      Biome_shift_proportion_temporal_pollen_NH_TEFO$Source <- "Pollen-based reconstruction"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_pollen_NH_TEFO_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "Pollen-based reconstruction" & Dataset_ID %in% DatasetID_list_NH_TEFO & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_pollen_NH_TEFO_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_pollen_NH_TEFO[i,2]    <- unique(Biome_shift_state_pollen_NH_TEFO_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_pollen_NH_TEFO[i,3]    <- unique(Biome_shift_state_pollen_NH_TEFO_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_pollen_NH_TEFO[i,4]    <- length(unique(Biome_shift_state_pollen_NH_TEFO_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_pollen_NH_TEFO[i,5]    <- sum(Biome_shift_state_pollen_NH_TEFO_subset$Shifting == "yes")/length(unique(Biome_shift_state_pollen_NH_TEFO_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_pollen_NH_TEFO[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_pollen_NH_TEFO[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_pollen_NH_TEFO[i,4]    <- 0
          Biome_shift_proportion_temporal_pollen_NH_TEFO[i,5]    <- NA
        }
        
      }
      
    }
    
    
    # BOFO
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_pollen_NH_BOFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_pollen_NH_BOFO) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_pollen_NH_BOFO[ ,1] <- 'BOFO'
      Biome_shift_proportion_temporal_pollen_NH_BOFO$Source <- "Pollen-based reconstruction"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_pollen_NH_BOFO_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "Pollen-based reconstruction" & Dataset_ID %in% DatasetID_list_NH_BOFO & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_pollen_NH_BOFO_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_pollen_NH_BOFO[i,2]    <- unique(Biome_shift_state_pollen_NH_BOFO_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_pollen_NH_BOFO[i,3]    <- unique(Biome_shift_state_pollen_NH_BOFO_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_pollen_NH_BOFO[i,4]    <- length(unique(Biome_shift_state_pollen_NH_BOFO_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_pollen_NH_BOFO[i,5]    <- sum(Biome_shift_state_pollen_NH_BOFO_subset$Shifting == "yes")/length(unique(Biome_shift_state_pollen_NH_BOFO_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_pollen_NH_BOFO[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_pollen_NH_BOFO[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_pollen_NH_BOFO[i,4]    <- 0
          Biome_shift_proportion_temporal_pollen_NH_BOFO[i,5]    <- NA
        }
        
      }
      
    }
    
    # SAVA
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_pollen_NH_SAVA           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_pollen_NH_SAVA) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_pollen_NH_SAVA[ ,1] <- 'SAVA'
      Biome_shift_proportion_temporal_pollen_NH_SAVA$Source <- "Pollen-based reconstruction"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_pollen_NH_SAVA_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "Pollen-based reconstruction" & Dataset_ID %in% DatasetID_list_NH_SAVA & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_pollen_NH_SAVA_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_pollen_NH_SAVA[i,2]    <- unique(Biome_shift_state_pollen_NH_SAVA_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_pollen_NH_SAVA[i,3]    <- unique(Biome_shift_state_pollen_NH_SAVA_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_pollen_NH_SAVA[i,4]    <- length(unique(Biome_shift_state_pollen_NH_SAVA_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_pollen_NH_SAVA[i,5]    <- sum(Biome_shift_state_pollen_NH_SAVA_subset$Shifting == "yes")/length(unique(Biome_shift_state_pollen_NH_SAVA_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_pollen_NH_SAVA[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_pollen_NH_SAVA[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_pollen_NH_SAVA[i,4]    <- 0
          Biome_shift_proportion_temporal_pollen_NH_SAVA[i,5]    <- NA
        }
        
      }
      
    }
    
    
    # STEP
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_pollen_NH_STEP           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_pollen_NH_STEP) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_pollen_NH_STEP[ ,1] <- 'STEP'
      Biome_shift_proportion_temporal_pollen_NH_STEP$Source <- "Pollen-based reconstruction"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_pollen_NH_STEP_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "Pollen-based reconstruction" & Dataset_ID %in% DatasetID_list_NH_STEP & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_pollen_NH_STEP_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_pollen_NH_STEP[i,2]    <- unique(Biome_shift_state_pollen_NH_STEP_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_pollen_NH_STEP[i,3]    <- unique(Biome_shift_state_pollen_NH_STEP_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_pollen_NH_STEP[i,4]    <- length(unique(Biome_shift_state_pollen_NH_STEP_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_pollen_NH_STEP[i,5]    <- sum(Biome_shift_state_pollen_NH_STEP_subset$Shifting == "yes")/length(unique(Biome_shift_state_pollen_NH_STEP_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_pollen_NH_STEP[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_pollen_NH_STEP[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_pollen_NH_STEP[i,4]    <- 0
          Biome_shift_proportion_temporal_pollen_NH_STEP[i,5]    <- NA
        }
        
      }
      
    }
    
    
    # DESE
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_pollen_NH_DESE           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_pollen_NH_DESE) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_pollen_NH_DESE[ ,1] <- 'DESE'
      Biome_shift_proportion_temporal_pollen_NH_DESE$Source <- "Pollen-based reconstruction"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_pollen_NH_DESE_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "Pollen-based reconstruction" & Dataset_ID %in% DatasetID_list_NH_DESE & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_pollen_NH_DESE_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_pollen_NH_DESE[i,2]    <- unique(Biome_shift_state_pollen_NH_DESE_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_pollen_NH_DESE[i,3]    <- unique(Biome_shift_state_pollen_NH_DESE_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_pollen_NH_DESE[i,4]    <- length(unique(Biome_shift_state_pollen_NH_DESE_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_pollen_NH_DESE[i,5]    <- sum(Biome_shift_state_pollen_NH_DESE_subset$Shifting == "yes")/length(unique(Biome_shift_state_pollen_NH_DESE_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_pollen_NH_DESE[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_pollen_NH_DESE[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_pollen_NH_DESE[i,4]    <- 0
          Biome_shift_proportion_temporal_pollen_NH_DESE[i,5]    <- NA
        }
        
      }
      
    }
    
    
    # TUND
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_pollen_NH_TUND           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_pollen_NH_TUND) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_pollen_NH_TUND[ ,1] <- 'TUND'
      Biome_shift_proportion_temporal_pollen_NH_TUND$Source <- "Pollen-based reconstruction"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_pollen_NH_TUND_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "Pollen-based reconstruction" & Dataset_ID %in% DatasetID_list_NH_TUND & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_pollen_NH_TUND_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_pollen_NH_TUND[i,2]    <- unique(Biome_shift_state_pollen_NH_TUND_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_pollen_NH_TUND[i,3]    <- unique(Biome_shift_state_pollen_NH_TUND_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_pollen_NH_TUND[i,4]    <- length(unique(Biome_shift_state_pollen_NH_TUND_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_pollen_NH_TUND[i,5]    <- sum(Biome_shift_state_pollen_NH_TUND_subset$Shifting == "yes")/length(unique(Biome_shift_state_pollen_NH_TUND_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_pollen_NH_TUND[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_pollen_NH_TUND[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_pollen_NH_TUND[i,4]    <- 0
          Biome_shift_proportion_temporal_pollen_NH_TUND[i,5]    <- NA
        }
        
      }
      
    }
    
    # NH
    {
      # Combine all biome-specific data frames into one for the Northern Hemisphere
      Biome_shift_proportion_temporal_pollen_NH_all  <- rbind(Biome_shift_proportion_temporal_pollen_NH_TRFO, Biome_shift_proportion_temporal_pollen_NH_WTFO, 
                                                              Biome_shift_proportion_temporal_pollen_NH_TEFO, Biome_shift_proportion_temporal_pollen_NH_BOFO, 
                                                              Biome_shift_proportion_temporal_pollen_NH_SAVA, Biome_shift_proportion_temporal_pollen_NH_STEP, 
                                                              Biome_shift_proportion_temporal_pollen_NH_DESE, Biome_shift_proportion_temporal_pollen_NH_TUND)
      
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_pollen_NH_NH           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_pollen_NH_NH ) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_pollen_NH_NH [ ,1] <- 'Northern Hemisphere'
      Biome_shift_proportion_temporal_pollen_NH_NH $Source <- "Pollen-based reconstruction"
      
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_proportion_temporal_pollen_NH_NH_subset <- subset(Biome_shift_proportion_temporal_pollen_NH_all, Timeslice_younger == list_timeslices[i])
        Biome_shift_proportion_temporal_pollen_NH_NH_subset[is.na(Biome_shift_proportion_temporal_pollen_NH_NH_subset)] <- 0
        
        # Store timeslice information
        Biome_shift_proportion_temporal_pollen_NH_NH [i,2]    <- unique(Biome_shift_proportion_temporal_pollen_NH_NH_subset$Timeslice_younger)
        Biome_shift_proportion_temporal_pollen_NH_NH [i,3]    <- unique(Biome_shift_proportion_temporal_pollen_NH_NH_subset$Timeslice_older)
        
        # Store the number of unique datasets for the current timeslice
        Biome_shift_proportion_temporal_pollen_NH_NH [i,4]    <- sum(Biome_shift_proportion_temporal_pollen_NH_NH_subset$Site_number)
        
        # Calculate the proportion of shifting biomes and store the value
        Biome_shift_proportion_temporal_pollen_NH_NH [i,5]    <- (Biome_shift_proportion_temporal_pollen_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_pollen_NH_NH_subset$Biome == 'TRFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==1] + 
                                                                    Biome_shift_proportion_temporal_pollen_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_pollen_NH_NH_subset$Biome == 'WTFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==2] +
                                                                    Biome_shift_proportion_temporal_pollen_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_pollen_NH_NH_subset$Biome == 'TEFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==3] + 
                                                                    Biome_shift_proportion_temporal_pollen_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_pollen_NH_NH_subset$Biome == 'BOFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==4] +
                                                                    Biome_shift_proportion_temporal_pollen_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_pollen_NH_NH_subset$Biome == 'SAVA']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==5] +
                                                                    Biome_shift_proportion_temporal_pollen_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_pollen_NH_NH_subset$Biome == 'STEP']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==6] +
                                                                    Biome_shift_proportion_temporal_pollen_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_pollen_NH_NH_subset$Biome == 'DESE']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==7] +
                                                                    Biome_shift_proportion_temporal_pollen_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_pollen_NH_NH_subset$Biome == 'TUND']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==8] )
        
      }
      
    }
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # (2) ESM-based simulation
  {
    # TRFO
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_MPI_NH_TRFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_MPI_NH_TRFO) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_MPI_NH_TRFO[ ,1] <- 'TRFO'
      Biome_shift_proportion_temporal_MPI_NH_TRFO$Source <- "ESM-based simulation"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_MPI_NH_TRFO_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "ESM-based simulation" & Dataset_ID %in% DatasetID_list_NH_TRFO & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_MPI_NH_TRFO_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_MPI_NH_TRFO[i,2]    <- unique(Biome_shift_state_MPI_NH_TRFO_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_MPI_NH_TRFO[i,3]    <- unique(Biome_shift_state_MPI_NH_TRFO_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_MPI_NH_TRFO[i,4]    <- length(unique(Biome_shift_state_MPI_NH_TRFO_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_MPI_NH_TRFO[i,5]    <- sum(Biome_shift_state_MPI_NH_TRFO_subset$Shifting == "yes")/length(unique(Biome_shift_state_MPI_NH_TRFO_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_MPI_NH_TRFO[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_MPI_NH_TRFO[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_MPI_NH_TRFO[i,4]    <- 0
          Biome_shift_proportion_temporal_MPI_NH_TRFO[i,5]    <- NA
        }
        
      }
      
    }
    
    # WTFO
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_MPI_NH_WTFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_MPI_NH_WTFO) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_MPI_NH_WTFO[ ,1] <- 'WTFO'
      Biome_shift_proportion_temporal_MPI_NH_WTFO$Source <- "ESM-based simulation"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_MPI_NH_WTFO_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "ESM-based simulation" & Dataset_ID %in% DatasetID_list_NH_WTFO & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_MPI_NH_WTFO_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_MPI_NH_WTFO[i,2]    <- unique(Biome_shift_state_MPI_NH_WTFO_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_MPI_NH_WTFO[i,3]    <- unique(Biome_shift_state_MPI_NH_WTFO_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_MPI_NH_WTFO[i,4]    <- length(unique(Biome_shift_state_MPI_NH_WTFO_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_MPI_NH_WTFO[i,5]    <- sum(Biome_shift_state_MPI_NH_WTFO_subset$Shifting == "yes")/length(unique(Biome_shift_state_MPI_NH_WTFO_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_MPI_NH_WTFO[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_MPI_NH_WTFO[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_MPI_NH_WTFO[i,4]    <- 0
          Biome_shift_proportion_temporal_MPI_NH_WTFO[i,5]    <- NA
        }
        
      }
      
    }
    
    # TEFO
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_MPI_NH_TEFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_MPI_NH_TEFO) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_MPI_NH_TEFO[ ,1] <- 'TEFO'
      Biome_shift_proportion_temporal_MPI_NH_TEFO$Source <- "ESM-based simulation"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_MPI_NH_TEFO_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "ESM-based simulation" & Dataset_ID %in% DatasetID_list_NH_TEFO & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_MPI_NH_TEFO_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_MPI_NH_TEFO[i,2]    <- unique(Biome_shift_state_MPI_NH_TEFO_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_MPI_NH_TEFO[i,3]    <- unique(Biome_shift_state_MPI_NH_TEFO_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_MPI_NH_TEFO[i,4]    <- length(unique(Biome_shift_state_MPI_NH_TEFO_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_MPI_NH_TEFO[i,5]    <- sum(Biome_shift_state_MPI_NH_TEFO_subset$Shifting == "yes")/length(unique(Biome_shift_state_MPI_NH_TEFO_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_MPI_NH_TEFO[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_MPI_NH_TEFO[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_MPI_NH_TEFO[i,4]    <- 0
          Biome_shift_proportion_temporal_MPI_NH_TEFO[i,5]    <- NA
        }
        
      }
      
    }
    
    
    # BOFO
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_MPI_NH_BOFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_MPI_NH_BOFO) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_MPI_NH_BOFO[ ,1] <- 'BOFO'
      Biome_shift_proportion_temporal_MPI_NH_BOFO$Source <- "ESM-based simulation"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_MPI_NH_BOFO_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "ESM-based simulation" & Dataset_ID %in% DatasetID_list_NH_BOFO & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_MPI_NH_BOFO_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_MPI_NH_BOFO[i,2]    <- unique(Biome_shift_state_MPI_NH_BOFO_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_MPI_NH_BOFO[i,3]    <- unique(Biome_shift_state_MPI_NH_BOFO_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_MPI_NH_BOFO[i,4]    <- length(unique(Biome_shift_state_MPI_NH_BOFO_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_MPI_NH_BOFO[i,5]    <- sum(Biome_shift_state_MPI_NH_BOFO_subset$Shifting == "yes")/length(unique(Biome_shift_state_MPI_NH_BOFO_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_MPI_NH_BOFO[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_MPI_NH_BOFO[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_MPI_NH_BOFO[i,4]    <- 0
          Biome_shift_proportion_temporal_MPI_NH_BOFO[i,5]    <- NA
        }
        
      }
      
    }
    
    # SAVA
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_MPI_NH_SAVA           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_MPI_NH_SAVA) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_MPI_NH_SAVA[ ,1] <- 'SAVA'
      Biome_shift_proportion_temporal_MPI_NH_SAVA$Source <- "ESM-based simulation"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_MPI_NH_SAVA_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "ESM-based simulation" & Dataset_ID %in% DatasetID_list_NH_SAVA & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_MPI_NH_SAVA_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_MPI_NH_SAVA[i,2]    <- unique(Biome_shift_state_MPI_NH_SAVA_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_MPI_NH_SAVA[i,3]    <- unique(Biome_shift_state_MPI_NH_SAVA_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_MPI_NH_SAVA[i,4]    <- length(unique(Biome_shift_state_MPI_NH_SAVA_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_MPI_NH_SAVA[i,5]    <- sum(Biome_shift_state_MPI_NH_SAVA_subset$Shifting == "yes")/length(unique(Biome_shift_state_MPI_NH_SAVA_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_MPI_NH_SAVA[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_MPI_NH_SAVA[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_MPI_NH_SAVA[i,4]    <- 0
          Biome_shift_proportion_temporal_MPI_NH_SAVA[i,5]    <- NA
        }
        
      }
      
    }
    
    
    # STEP
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_MPI_NH_STEP           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_MPI_NH_STEP) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_MPI_NH_STEP[ ,1] <- 'STEP'
      Biome_shift_proportion_temporal_MPI_NH_STEP$Source <- "ESM-based simulation"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_MPI_NH_STEP_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "ESM-based simulation" & Dataset_ID %in% DatasetID_list_NH_STEP & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_MPI_NH_STEP_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_MPI_NH_STEP[i,2]    <- unique(Biome_shift_state_MPI_NH_STEP_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_MPI_NH_STEP[i,3]    <- unique(Biome_shift_state_MPI_NH_STEP_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_MPI_NH_STEP[i,4]    <- length(unique(Biome_shift_state_MPI_NH_STEP_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_MPI_NH_STEP[i,5]    <- sum(Biome_shift_state_MPI_NH_STEP_subset$Shifting == "yes")/length(unique(Biome_shift_state_MPI_NH_STEP_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_MPI_NH_STEP[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_MPI_NH_STEP[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_MPI_NH_STEP[i,4]    <- 0
          Biome_shift_proportion_temporal_MPI_NH_STEP[i,5]    <- NA
        }
        
      }
      
    }
    
    
    # DESE
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_MPI_NH_DESE           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_MPI_NH_DESE) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_MPI_NH_DESE[ ,1] <- 'DESE'
      Biome_shift_proportion_temporal_MPI_NH_DESE$Source <- "ESM-based simulation"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_MPI_NH_DESE_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "ESM-based simulation" & Dataset_ID %in% DatasetID_list_NH_DESE & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_MPI_NH_DESE_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_MPI_NH_DESE[i,2]    <- unique(Biome_shift_state_MPI_NH_DESE_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_MPI_NH_DESE[i,3]    <- unique(Biome_shift_state_MPI_NH_DESE_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_MPI_NH_DESE[i,4]    <- length(unique(Biome_shift_state_MPI_NH_DESE_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_MPI_NH_DESE[i,5]    <- sum(Biome_shift_state_MPI_NH_DESE_subset$Shifting == "yes")/length(unique(Biome_shift_state_MPI_NH_DESE_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_MPI_NH_DESE[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_MPI_NH_DESE[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_MPI_NH_DESE[i,4]    <- 0
          Biome_shift_proportion_temporal_MPI_NH_DESE[i,5]    <- NA
        }
        
      }
      
    }
    
    
    # TUND
    {
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_MPI_NH_TUND           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_MPI_NH_TUND) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_MPI_NH_TUND[ ,1] <- 'TUND'
      Biome_shift_proportion_temporal_MPI_NH_TUND$Source <- "ESM-based simulation"
      
      # Loop through each timeslice to calculate biome shift proportions
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_state_MPI_NH_TUND_subset <- subset(Biome_shift_state_pollen_MPI_global, Source == "ESM-based simulation" & Dataset_ID %in% DatasetID_list_NH_TUND & Timeslice_younger == list_timeslices[i])
        
        # Check if there are records for the current timeslice
        if(nrow(Biome_shift_state_MPI_NH_TUND_subset) > 0) {
          
          # Store timeslice information
          Biome_shift_proportion_temporal_MPI_NH_TUND[i,2]    <- unique(Biome_shift_state_MPI_NH_TUND_subset$Timeslice_younger)
          Biome_shift_proportion_temporal_MPI_NH_TUND[i,3]    <- unique(Biome_shift_state_MPI_NH_TUND_subset$Timeslice_older)
          
          # Store the number of unique datasets for the current timeslice
          Biome_shift_proportion_temporal_MPI_NH_TUND[i,4]    <- length(unique(Biome_shift_state_MPI_NH_TUND_subset$Dataset_ID))
          
          # Calculate the proportion of shifting biomes and store the value
          Biome_shift_proportion_temporal_MPI_NH_TUND[i,5]    <- sum(Biome_shift_state_MPI_NH_TUND_subset$Shifting == "yes")/length(unique(Biome_shift_state_MPI_NH_TUND_subset$Dataset_ID))
          
        } else{
          
          # If no records exist, fill with default values
          Biome_shift_proportion_temporal_MPI_NH_TUND[i,2]    <- list_timeslices[i]
          Biome_shift_proportion_temporal_MPI_NH_TUND[i,3]    <- list_timeslices[i] + 0.5
          Biome_shift_proportion_temporal_MPI_NH_TUND[i,4]    <- 0
          Biome_shift_proportion_temporal_MPI_NH_TUND[i,5]    <- NA
        }
        
      }
      
    }
    
    # NH
    {
      # Combine all biome-specific data frames into one for the Northern Hemisphere
      Biome_shift_proportion_temporal_MPI_NH_all  <- rbind(Biome_shift_proportion_temporal_MPI_NH_TRFO, Biome_shift_proportion_temporal_MPI_NH_WTFO, 
                                                           Biome_shift_proportion_temporal_MPI_NH_TEFO, Biome_shift_proportion_temporal_MPI_NH_BOFO, 
                                                           Biome_shift_proportion_temporal_MPI_NH_SAVA, Biome_shift_proportion_temporal_MPI_NH_STEP, 
                                                           Biome_shift_proportion_temporal_MPI_NH_DESE, Biome_shift_proportion_temporal_MPI_NH_TUND)
      
      # Define a sequence of timeslices
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize a data frame to store biome shift proportions over time
      Biome_shift_proportion_temporal_MPI_NH_NH           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 6, byrow=TRUE))
      colnames(Biome_shift_proportion_temporal_MPI_NH_NH ) <- c("Biome", "Timeslice_younger", "Timeslice_older",  "Site_number", "Proportion", "Source")
      
      # Set the biome name and specify the source
      Biome_shift_proportion_temporal_MPI_NH_NH [ ,1] <- 'Northern Hemisphere'
      Biome_shift_proportion_temporal_MPI_NH_NH $Source <- "ESM-based simulation"
      
      for (i in 1:length(list_timeslices)) {
        
        # Print progress message for the current timeslice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Subset the dataset for the current timeslice, biome, and matching conditions
        Biome_shift_proportion_temporal_MPI_NH_NH_subset <- subset(Biome_shift_proportion_temporal_MPI_NH_all, Timeslice_younger == list_timeslices[i])
        Biome_shift_proportion_temporal_MPI_NH_NH_subset[is.na(Biome_shift_proportion_temporal_MPI_NH_NH_subset)] <- 0
        
        # Store timeslice information
        Biome_shift_proportion_temporal_MPI_NH_NH [i,2]    <- unique(Biome_shift_proportion_temporal_MPI_NH_NH_subset$Timeslice_younger)
        Biome_shift_proportion_temporal_MPI_NH_NH [i,3]    <- unique(Biome_shift_proportion_temporal_MPI_NH_NH_subset$Timeslice_older)
        
        # Store the number of unique datasets for the current timeslice
        Biome_shift_proportion_temporal_MPI_NH_NH [i,4]    <- sum(Biome_shift_proportion_temporal_MPI_NH_NH_subset$Site_number)
        
        # Calculate the proportion of shifting biomes and store the value
        Biome_shift_proportion_temporal_MPI_NH_NH [i,5]    <- (Biome_shift_proportion_temporal_MPI_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_MPI_NH_NH_subset$Biome == 'TRFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==1] + 
                                                                 Biome_shift_proportion_temporal_MPI_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_MPI_NH_NH_subset$Biome == 'WTFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==2] +
                                                                 Biome_shift_proportion_temporal_MPI_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_MPI_NH_NH_subset$Biome == 'TEFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==3] + 
                                                                 Biome_shift_proportion_temporal_MPI_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_MPI_NH_NH_subset$Biome == 'BOFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==4] +
                                                                 Biome_shift_proportion_temporal_MPI_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_MPI_NH_NH_subset$Biome == 'SAVA']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==5] +
                                                                 Biome_shift_proportion_temporal_MPI_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_MPI_NH_NH_subset$Biome == 'STEP']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==6] +
                                                                 Biome_shift_proportion_temporal_MPI_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_MPI_NH_NH_subset$Biome == 'DESE']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==7] +
                                                                 Biome_shift_proportion_temporal_MPI_NH_NH_subset$Proportion[Biome_shift_proportion_temporal_MPI_NH_NH_subset$Biome == 'TUND']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==8] )
        
      }
      
    }
    
  }
  
  
  # Combine result dataset
  Biome_shift_proportion_temporal_pollen_MPI_NH_final <- rbind(Biome_shift_proportion_temporal_pollen_NH_NH, Biome_shift_proportion_temporal_pollen_NH_all, Biome_shift_proportion_temporal_MPI_NH_all, Biome_shift_proportion_temporal_MPI_NH_NH)
  
  # Save the final results as a CSV file
  write.csv(Biome_shift_proportion_temporal_pollen_MPI_NH_final, file="result/Statistical analysis 2_result.3-Temporal changes of proportion of biome shifts in the Northern Hemisphere over the last 21,000 years.csv", row.names=FALSE)
  
}


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 3. Statistical analysis 3: Simulated bioclimatic variables_21-0ka (Asia as example) ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# 3.1 Calculation: state and rate of change of bioclimatic variables at timeslice per site in Asia
# -------------------------------------------------------------------------------------------------

# Load datasets
Simulated_climate_timeslice_site_Asia    <- read.csv2("data/Statistical analysis 3_data.1-Simulated monthly temperature and precipitation at timeslice per site in Asia.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Simulated_climate_timeslice_site_Asia    <- type.convert(Simulated_climate_timeslice_site_Asia, as.is = TRUE) 

# Calculation
{
  # (1) Simulated bioclimatic variables
  {
    # Extract the unique Dataset_ID values from the data
    list_ID <- unique(Simulated_climate_timeslice_site_Asia$Dataset_ID)
    
    # Initialize an empty variable to store the results later
    Simulated_bioclimatic_variable_timeslice_site_Asia <- NULL
    
    # Create a sequence of time slices starting from 0 to 21 with a step size of 0.5
    list_timeslices <- seq(0,21,by=0.5) # seq (starting age, ending age, step)
    
    # Loop through each unique Dataset_ID
    for (ID in list_ID) {
      
      # Print progress for the current Dataset_ID
      print(paste0("+++++ Dataset_ID ",ID," | (",which(list_ID == ID),"/",length(list_ID),") +++++"))
      
      # Subset the data for the current Dataset_ID
      subset_site <- subset(Simulated_climate_timeslice_site_Asia, Dataset_ID == ID) 
      
      # Initialize a data frame for the time slice data for the current site
      bioclimatic_variables_timeslice_site_Asia           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 19, byrow=TRUE))
      colnames(bioclimatic_variables_timeslice_site_Asia) <- c("Source", "Dataset_ID",  "Longitude", "Latitude", "Timeslice", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Fill in metadata columns
      bioclimatic_variables_timeslice_site_Asia[ ,1] <- "MPI-ESM"
      bioclimatic_variables_timeslice_site_Asia[ ,2] <- ID
      bioclimatic_variables_timeslice_site_Asia[ ,3] <- unique(subset_site$Longitude)
      bioclimatic_variables_timeslice_site_Asia[ ,4] <- unique(subset_site$Latitude)
      
      # Loop through each unique timeslice
      for (i in 1:length(list_timeslices)) {
        
        # Subset data for the current timeslice
        subset_site_subset <- subset(subset_site, Timeslice == list_timeslices[i])
        
        # Assign the current timeslice identifier
        bioclimatic_variables_timeslice_site_Asia[i, 5] <- list_timeslices[i]
        
        
        # Separate Quarter_mean
        ## precipitation
        subset_site_subset$precipitation_1_2_3    <- mean(subset_site_subset$Precipitation[subset_site_subset$Month == 1] + subset_site_subset$Precipitation[subset_site_subset$Month == 2] + subset_site_subset$Precipitation[subset_site_subset$Month == 3])
        subset_site_subset$precipitation_2_3_4    <- mean(subset_site_subset$Precipitation[subset_site_subset$Month == 2] + subset_site_subset$Precipitation[subset_site_subset$Month == 3] + subset_site_subset$Precipitation[subset_site_subset$Month == 4])
        subset_site_subset$precipitation_3_4_5    <- mean(subset_site_subset$Precipitation[subset_site_subset$Month == 3] + subset_site_subset$Precipitation[subset_site_subset$Month == 4] + subset_site_subset$Precipitation[subset_site_subset$Month == 5])
        subset_site_subset$precipitation_4_5_6    <- mean(subset_site_subset$Precipitation[subset_site_subset$Month == 4] + subset_site_subset$Precipitation[subset_site_subset$Month == 5] + subset_site_subset$Precipitation[subset_site_subset$Month == 6])
        subset_site_subset$precipitation_5_6_7    <- mean(subset_site_subset$Precipitation[subset_site_subset$Month == 5] + subset_site_subset$Precipitation[subset_site_subset$Month == 6] + subset_site_subset$Precipitation[subset_site_subset$Month == 7])
        subset_site_subset$precipitation_6_7_8    <- mean(subset_site_subset$Precipitation[subset_site_subset$Month == 6] + subset_site_subset$Precipitation[subset_site_subset$Month == 7] + subset_site_subset$Precipitation[subset_site_subset$Month == 8])
        subset_site_subset$precipitation_7_8_9    <- mean(subset_site_subset$Precipitation[subset_site_subset$Month == 7] + subset_site_subset$Precipitation[subset_site_subset$Month == 8] + subset_site_subset$Precipitation[subset_site_subset$Month == 9])
        subset_site_subset$precipitation_8_9_10   <- mean(subset_site_subset$Precipitation[subset_site_subset$Month == 8] + subset_site_subset$Precipitation[subset_site_subset$Month == 9] + subset_site_subset$Precipitation[subset_site_subset$Month == 10])
        subset_site_subset$precipitation_9_10_11  <- mean(subset_site_subset$Precipitation[subset_site_subset$Month == 9] + subset_site_subset$Precipitation[subset_site_subset$Month == 10] + subset_site_subset$Precipitation[subset_site_subset$Month == 11])
        subset_site_subset$precipitation_10_11_12 <- mean(subset_site_subset$Precipitation[subset_site_subset$Month == 10] + subset_site_subset$Precipitation[subset_site_subset$Month == 11] + subset_site_subset$Precipitation[subset_site_subset$Month == 12])
        subset_site_subset$precipitation_11_12_1  <- mean(subset_site_subset$Precipitation[subset_site_subset$Month == 11] + subset_site_subset$Precipitation[subset_site_subset$Month == 12] + subset_site_subset$Precipitation[subset_site_subset$Month == 1])
        subset_site_subset$precipitation_12_1_2   <- mean(subset_site_subset$Precipitation[subset_site_subset$Month == 12] + subset_site_subset$Precipitation[subset_site_subset$Month == 1] + subset_site_subset$Precipitation[subset_site_subset$Month == 2])
        
        subset_site_subset$Wettest_Month    <- max(subset_site_subset$Precipitation)
        subset_site_subset$Driest_Month     <- min(subset_site_subset$Precipitation)
        
        subset_site_subset$Wettest_Quarter  <- apply(subset_site_subset[ ,which(names(subset_site_subset) %in% c("precipitation_1_2_3", "precipitation_2_3_4", "precipitation_3_4_5", "precipitation_4_5_6", "precipitation_5_6_7", "precipitation_6_7_8", "precipitation_7_8_9", "precipitation_8_9_10", "precipitation_9_10_11", "precipitation_10_11_12", "precipitation_11_12_1", "precipitation_12_1_2"))], 1, function(x) { names(x)[which.max(x)] })
        subset_site_subset$Driest_Quarter   <- apply(subset_site_subset[ ,which(names(subset_site_subset) %in% c("precipitation_1_2_3", "precipitation_2_3_4", "precipitation_3_4_5", "precipitation_4_5_6", "precipitation_5_6_7", "precipitation_6_7_8", "precipitation_7_8_9", "precipitation_8_9_10", "precipitation_9_10_11", "precipitation_10_11_12", "precipitation_11_12_1", "precipitation_12_1_2"))], 1, function(x) { names(x)[which.min(x)] })
        
        ## temperature
        subset_site_subset$temperature_1_2_3    <- mean(subset_site_subset$Temperature[subset_site_subset$Month == 1] + subset_site_subset$Temperature[subset_site_subset$Month == 2] + subset_site_subset$Temperature[subset_site_subset$Month == 3])
        subset_site_subset$temperature_2_3_4    <- mean(subset_site_subset$Temperature[subset_site_subset$Month == 2] + subset_site_subset$Temperature[subset_site_subset$Month == 3] + subset_site_subset$Temperature[subset_site_subset$Month == 4])
        subset_site_subset$temperature_3_4_5    <- mean(subset_site_subset$Temperature[subset_site_subset$Month == 3] + subset_site_subset$Temperature[subset_site_subset$Month == 4] + subset_site_subset$Temperature[subset_site_subset$Month == 5])
        subset_site_subset$temperature_4_5_6    <- mean(subset_site_subset$Temperature[subset_site_subset$Month == 4] + subset_site_subset$Temperature[subset_site_subset$Month == 5] + subset_site_subset$Temperature[subset_site_subset$Month == 6])
        subset_site_subset$temperature_5_6_7    <- mean(subset_site_subset$Temperature[subset_site_subset$Month == 5] + subset_site_subset$Temperature[subset_site_subset$Month == 6] + subset_site_subset$Temperature[subset_site_subset$Month == 7])
        subset_site_subset$temperature_6_7_8    <- mean(subset_site_subset$Temperature[subset_site_subset$Month == 6] + subset_site_subset$Temperature[subset_site_subset$Month == 7] + subset_site_subset$Temperature[subset_site_subset$Month == 8])
        subset_site_subset$temperature_7_8_9    <- mean(subset_site_subset$Temperature[subset_site_subset$Month == 7] + subset_site_subset$Temperature[subset_site_subset$Month == 8] + subset_site_subset$Temperature[subset_site_subset$Month == 9])
        subset_site_subset$temperature_8_9_10   <- mean(subset_site_subset$Temperature[subset_site_subset$Month == 8] + subset_site_subset$Temperature[subset_site_subset$Month == 9] + subset_site_subset$Temperature[subset_site_subset$Month == 10])
        subset_site_subset$temperature_9_10_11  <- mean(subset_site_subset$Temperature[subset_site_subset$Month == 9] + subset_site_subset$Temperature[subset_site_subset$Month == 10] + subset_site_subset$Temperature[subset_site_subset$Month == 11])
        subset_site_subset$temperature_10_11_12 <- mean(subset_site_subset$Temperature[subset_site_subset$Month == 10] + subset_site_subset$Temperature[subset_site_subset$Month == 11] + subset_site_subset$Temperature[subset_site_subset$Month == 12])
        subset_site_subset$temperature_11_12_1  <- mean(subset_site_subset$Temperature[subset_site_subset$Month == 11] + subset_site_subset$Temperature[subset_site_subset$Month == 12] + subset_site_subset$Temperature[subset_site_subset$Month == 1])
        subset_site_subset$temperature_12_1_2   <- mean(subset_site_subset$Temperature[subset_site_subset$Month == 12] + subset_site_subset$Temperature[subset_site_subset$Month == 1] + subset_site_subset$Temperature[subset_site_subset$Month == 2])
        
        subset_site_subset$Warmest_Month   <- max(subset_site_subset$Temperature)
        subset_site_subset$Coldest_Month   <- min(subset_site_subset$Temperature)
        
        subset_site_subset$Warmest_Quarter   <- apply(subset_site_subset[ ,which(names(subset_site_subset) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.max(x)] })
        subset_site_subset$Coldest_Quarter   <- apply(subset_site_subset[ ,which(names(subset_site_subset) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.min(x)] })
        
        # Calculate BIO
        ## BIO 1 = Annual Mean Temperature
        bioclimatic_variables_timeslice_site_Asia[i, 6] <- mean(subset_site_subset$Temperature)
        
        ## BIO 4 = Temperature Seasonality (standard deviation ×100)
        bioclimatic_variables_timeslice_site_Asia[i, 7] <- sd(subset_site_subset$Temperature)*100
        
        ## BIO 8 = Mean Temperature of Wettest Quarter
        bioclimatic_variables_timeslice_site_Asia[i, 8] <- (subset_site_subset$Temperature[subset_site_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_site_subset$Wettest_Quarter)))] +
                                                              subset_site_subset$Temperature[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_site_subset$Wettest_Quarter)))] +  
                                                              subset_site_subset$Temperature[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_site_subset$Wettest_Quarter)))])/3
        
        ## BIO 9 = Mean Temperature of Driest Quarter
        bioclimatic_variables_timeslice_site_Asia[i, 9] <- (subset_site_subset$Temperature[subset_site_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_site_subset$Driest_Quarter)))] +
                                                              subset_site_subset$Temperature[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_site_subset$Driest_Quarter)))] +  
                                                              subset_site_subset$Temperature[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_site_subset$Driest_Quarter)))])/3
        
        ## BIO10 = Mean Temperature of Warmest Quarter
        bioclimatic_variables_timeslice_site_Asia[i, 10] <- (subset_site_subset$Temperature[subset_site_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_site_subset$Warmest_Quarter)))] +
                                                               subset_site_subset$Temperature[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_site_subset$Warmest_Quarter)))] +  
                                                               subset_site_subset$Temperature[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_site_subset$Warmest_Quarter)))])/3
        
        ## BIO11 = Mean Temperature of Coldest Quarter
        bioclimatic_variables_timeslice_site_Asia[i, 11] <- (subset_site_subset$Temperature[subset_site_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_site_subset$Coldest_Quarter)))] +
                                                               subset_site_subset$Temperature[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_site_subset$Coldest_Quarter)))] +  
                                                               subset_site_subset$Temperature[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_site_subset$Coldest_Quarter)))])/3
        
        ## BIO12 = Annual Precipitation
        bioclimatic_variables_timeslice_site_Asia[i, 12] <- sum(subset_site_subset$Precipitation)
        
        ## BIO13 = Precipitation of Wettest Month
        bioclimatic_variables_timeslice_site_Asia[i, 13] <- max(subset_site_subset$Precipitation)
        
        ## BIO14 = Precipitation of Driest Month
        bioclimatic_variables_timeslice_site_Asia[i, 14] <- min(subset_site_subset$Precipitation)
        
        ## BIO15 = Precipitation Seasonality (Coefficient of Variation)
        bioclimatic_variables_timeslice_site_Asia[i, 15] <- sd(subset_site_subset$Precipitation)/(1+ mean(subset_site_subset$Precipitation))*100
        
        ## BIO16 = Precipitation of Wettest Quarter
        bioclimatic_variables_timeslice_site_Asia[i, 16] <- (subset_site_subset$Precipitation[subset_site_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_site_subset$Wettest_Quarter)))] +
                                                               subset_site_subset$Precipitation[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_site_subset$Wettest_Quarter)))] +  
                                                               subset_site_subset$Precipitation[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_site_subset$Wettest_Quarter)))])
        
        ## BIO17 = Precipitation of Driest Quarter 
        bioclimatic_variables_timeslice_site_Asia[i, 17] <- (subset_site_subset$Precipitation[subset_site_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_site_subset$Driest_Quarter)))] +
                                                               subset_site_subset$Precipitation[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_site_subset$Driest_Quarter)))] +  
                                                               subset_site_subset$Precipitation[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_site_subset$Driest_Quarter)))])
        
        ## BIO18 = Precipitation of Warmest Quarter
        bioclimatic_variables_timeslice_site_Asia[i, 18] <- (subset_site_subset$Precipitation[subset_site_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_site_subset$Warmest_Quarter)))] +
                                                               subset_site_subset$Precipitation[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_site_subset$Warmest_Quarter)))] +  
                                                               subset_site_subset$Precipitation[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_site_subset$Warmest_Quarter)))])
        
        ## BIO19 = Precipitation of Coldest Quarter
        bioclimatic_variables_timeslice_site_Asia[i, 19] <- (subset_site_subset$Precipitation[subset_site_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(subset_site_subset$Coldest_Quarter)))] +
                                                               subset_site_subset$Precipitation[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(subset_site_subset$Coldest_Quarter)))] +  
                                                               subset_site_subset$Precipitation[subset_site_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(subset_site_subset$Coldest_Quarter)))])
        
        
      }
      
      # Append the result to the main result container
      Simulated_bioclimatic_variable_timeslice_site_Asia <- rbind(Simulated_bioclimatic_variable_timeslice_site_Asia, bioclimatic_variables_timeslice_site_Asia)
      
      
    }
    
    # save data frame to csv:
    write.csv(Simulated_bioclimatic_variable_timeslice_site_Asia, file="result/Statistical analysis 3_result.1-Simulated bioclimatic variables at timeslice per site in Asia.csv", row.names=FALSE)
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # (2) Change rate of simulated bioclimatic variables
  {
    # Extract the unique Dataset_ID values from the data
    list_ID <- unique(Simulated_bioclimatic_variable_timeslice_site_Asia$Dataset_ID)
    
    # Initialize an empty variable to store the results later
    Simulated_bioclimatic_variable_change_rate_timeslice_site_Asia <- NULL
    
    # Loop through each unique Dataset_ID
    for (ID in list_ID) {
      
      # Print progress for the current Dataset_ID
      print(paste0("+++++ Dataset_ID ",ID," | (",which(list_ID == ID),"/",length(list_ID),") +++++"))
      
      # Subset the data for the current Dataset_ID
      subset_site <- subset(Simulated_bioclimatic_variable_timeslice_site_Asia, Dataset_ID == ID) 
      
      # Initialize a data frame for the time slice data for the current site
      bioclimatic_variables_change_rate_site_Asia <- data.frame(matrix(NA, nrow=nrow(subset_site)-1, ncol=21, byrow=TRUE))
      colnames(bioclimatic_variables_change_rate_site_Asia) <- c("Source", "Dataset_ID", "Longitude", "Latitude", "Timeslice_younger", "Timeslice_older", "Timeslice_span", 
                                                                 "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Fill in metadata columns
      bioclimatic_variables_change_rate_site_Asia$Source <- "MPI-ESM"
      bioclimatic_variables_change_rate_site_Asia$Dataset_ID <- unique(subset_site$Dataset_ID)
      bioclimatic_variables_change_rate_site_Asia$Longitude  <- unique(subset_site$Longitude)
      bioclimatic_variables_change_rate_site_Asia$Latitude   <- unique(subset_site$Latitude)
      
      # Loop through time slices for the current site
      for(i in 1:c(nrow(subset_site)-1)){ 
        
        # Assign the younger and older timeslices and compute the time span
        bioclimatic_variables_change_rate_site_Asia[i,5] <- subset_site[i,5]
        bioclimatic_variables_change_rate_site_Asia[i,6] <- subset_site[i+1,5]
        bioclimatic_variables_change_rate_site_Asia[i,7] <- subset_site[i+1,5] - subset_site[i,5]
        
        # Calculate the change rate for bioclimatic variables (columns 6:19 in the input dataset)
        for(k in 6:19){ 
          
          bioclimatic_variables_change_rate_site_Asia[i,k+2] <- abs((subset_site[i,k] - subset_site[i+1,k])/subset_site[i+1,k])*100 
          
        }
        
      }
      
      # Select rows where the timeslice span is exactly 0.5
      bioclimatic_variables_change_rate_site_Asia <- type.convert(bioclimatic_variables_change_rate_site_Asia[bioclimatic_variables_change_rate_site_Asia$Timeslice_span == 0.5, ], as.is = TRUE) 
      
      # Combine results for all sites
      Simulated_bioclimatic_variable_change_rate_timeslice_site_Asia <- rbind(Simulated_bioclimatic_variable_change_rate_timeslice_site_Asia, bioclimatic_variables_change_rate_site_Asia)
      # Remove rows with missing values
      Simulated_bioclimatic_variable_change_rate_timeslice_site_Asia  <- Simulated_bioclimatic_variable_change_rate_timeslice_site_Asia[complete.cases(Simulated_bioclimatic_variable_change_rate_timeslice_site_Asia), ]
      
    }
    
    
    # save data frame to csv:
    write.csv(Simulated_bioclimatic_variable_change_rate_timeslice_site_Asia, file="result/Statistical analysis 3_result.2-Change rate of simulated bioclimatic variables at timeslice per site in Asia.csv", row.names=FALSE)
    
  }
  
}


# 3.2 Calculation: Temporal changes of bioclimatic variables and their change rates in the Northern Hemisphere over the last 21,000 years by modern potential natural biome distributions
# -------------------------------------------------------------------------------------------------

# Load datasets
Simulated_bioclimatic_variable_timeslice_site_NH                <- read.csv2("data/Statistical analysis 3_data.2-Simulated bioclimatic variables at timeslice per site in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Simulated_bioclimatic_variable_change_rate_timeslice_site_NH    <- read.csv2("data/Statistical analysis 3_data.3-Change rate of simulated bioclimatic variables at timeslice per site in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_potential_natural_biomes_site_NH                         <- read.csv2("data/Statistical analysis 3_data.4-Modern potential natural biomes per site in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Modern_potential_natural_biomes_grid_NH                         <- read.csv2("data/Statistical analysis 3_data.5-Modern potential natural biomes per grid-cell (spatial resolution 5 arc minutes) in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Simulated_bioclimatic_variable_timeslice_site_NH                <- type.convert(Simulated_bioclimatic_variable_timeslice_site_NH, as.is = TRUE) 
Simulated_bioclimatic_variable_change_rate_timeslice_site_NH    <- type.convert(Simulated_bioclimatic_variable_change_rate_timeslice_site_NH, as.is = TRUE) 
Modern_potential_natural_biomes_site_NH                         <- type.convert(Modern_potential_natural_biomes_site_NH, as.is = TRUE) 
Modern_potential_natural_biomes_grid_NH                         <- type.convert(Modern_potential_natural_biomes_grid_NH, as.is = TRUE) 

# Calculation
{
  # Calculate the proportion of modern potential natural biomes
  {
    # Set up a results data frame
    Modern_potential_natural_biomes_proportion_NH <- data.frame(matrix(NA, nrow=8, ncol=3, byrow=TRUE))
    colnames(Modern_potential_natural_biomes_proportion_NH) <- c("Biome", "Grid_number", "Proportion")
    
    # Populate the Biome column with biome IDs (1 to 8)
    Modern_potential_natural_biomes_proportion_NH$Biome        <- 1:8
    # Calculate the number of grid cells for each biome
    Modern_potential_natural_biomes_proportion_NH$Grid_number  <- table(Modern_potential_natural_biomes_grid_NH$Biome)
    # Calculate the proportion of grid cells for each biome
    Modern_potential_natural_biomes_proportion_NH$Proportion   <- prop.table(table(Modern_potential_natural_biomes_grid_NH$Biome))
  }
  
  # Extract Dataset ID lists for each biome
  {
    DatasetID_list_NH_TRFO <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 1]
    DatasetID_list_NH_WTFO <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 2]
    DatasetID_list_NH_TEFO <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 3]
    DatasetID_list_NH_BOFO <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 4]
    DatasetID_list_NH_SAVA <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 5]
    DatasetID_list_NH_STEP <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 6]
    DatasetID_list_NH_DESE <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 7]
    DatasetID_list_NH_TUND <- Modern_potential_natural_biomes_site_NH$Dataset_ID[Modern_potential_natural_biomes_site_NH$Biome == 8]
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # (1) Simulated bioclimatic variables
  {
    # TRFO
    {
      # Define the range of time slices to analyze (from 0 to 21, in increments of 0.5)
      list_timeslices <- seq(0,21,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_state_timeslice_site_NH_TRFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 17, byrow=TRUE))
      colnames(Bioclimatic_variables_state_timeslice_site_NH_TRFO) <- c("Biome", "Timeslice", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_state_timeslice_site_NH_TRFO$Biome <- 'TRFO'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_state_timeslice_site_NH_TRFO_subset <- subset(Simulated_bioclimatic_variable_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_TRFO & Timeslice == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_state_timeslice_site_NH_TRFO_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_state_timeslice_site_NH_TRFO[i,2]    <- unique(Bioclimatic_variables_state_timeslice_site_NH_TRFO_subset$Timeslice)
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_state_timeslice_site_NH_TRFO[i,3]    <- length(unique(Bioclimatic_variables_state_timeslice_site_NH_TRFO_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 8:21) {
            
            Bioclimatic_variables_state_timeslice_site_NH_TRFO[i,j-4] <- median(Bioclimatic_variables_state_timeslice_site_NH_TRFO_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_state_timeslice_site_NH_TRFO[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_state_timeslice_site_NH_TRFO[i,3]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_state_timeslice_site_NH_TRFO[i,4:17] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # WTFO
    {
      # Define the range of time slices to analyze (from 0 to 21, in increments of 0.5)
      list_timeslices <- seq(0,21,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_state_timeslice_site_NH_WTFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 17, byrow=TRUE))
      colnames(Bioclimatic_variables_state_timeslice_site_NH_WTFO) <- c("Biome", "Timeslice", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_state_timeslice_site_NH_WTFO$Biome <- 'WTFO'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_state_timeslice_site_NH_WTFO_subset <- subset(Simulated_bioclimatic_variable_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_WTFO & Timeslice == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_state_timeslice_site_NH_WTFO_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_state_timeslice_site_NH_WTFO[i,2]    <- unique(Bioclimatic_variables_state_timeslice_site_NH_WTFO_subset$Timeslice)
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_state_timeslice_site_NH_WTFO[i,3]    <- length(unique(Bioclimatic_variables_state_timeslice_site_NH_WTFO_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 8:21) {
            
            Bioclimatic_variables_state_timeslice_site_NH_WTFO[i,j-4] <- median(Bioclimatic_variables_state_timeslice_site_NH_WTFO_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_state_timeslice_site_NH_WTFO[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_state_timeslice_site_NH_WTFO[i,3]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_state_timeslice_site_NH_WTFO[i,4:17] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # TEFO
    {
      # Define the range of time slices to analyze (from 0 to 21, in increments of 0.5)
      list_timeslices <- seq(0,21,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_state_timeslice_site_NH_TEFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 17, byrow=TRUE))
      colnames(Bioclimatic_variables_state_timeslice_site_NH_TEFO) <- c("Biome", "Timeslice", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_state_timeslice_site_NH_TEFO$Biome <- 'TEFO'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_state_timeslice_site_NH_TEFO_subset <- subset(Simulated_bioclimatic_variable_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_TEFO & Timeslice == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_state_timeslice_site_NH_TEFO_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_state_timeslice_site_NH_TEFO[i,2]    <- unique(Bioclimatic_variables_state_timeslice_site_NH_TEFO_subset$Timeslice)
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_state_timeslice_site_NH_TEFO[i,3]    <- length(unique(Bioclimatic_variables_state_timeslice_site_NH_TEFO_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 8:21) {
            
            Bioclimatic_variables_state_timeslice_site_NH_TEFO[i,j-4] <- median(Bioclimatic_variables_state_timeslice_site_NH_TEFO_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_state_timeslice_site_NH_TEFO[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_state_timeslice_site_NH_TEFO[i,3]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_state_timeslice_site_NH_TEFO[i,4:17] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # BOFO
    {
      # Define the range of time slices to analyze (from 0 to 21, in increments of 0.5)
      list_timeslices <- seq(0,21,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_state_timeslice_site_NH_BOFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 17, byrow=TRUE))
      colnames(Bioclimatic_variables_state_timeslice_site_NH_BOFO) <- c("Biome", "Timeslice", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_state_timeslice_site_NH_BOFO$Biome <- 'BOFO'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_state_timeslice_site_NH_BOFO_subset <- subset(Simulated_bioclimatic_variable_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_BOFO & Timeslice == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_state_timeslice_site_NH_BOFO_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_state_timeslice_site_NH_BOFO[i,2]    <- unique(Bioclimatic_variables_state_timeslice_site_NH_BOFO_subset$Timeslice)
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_state_timeslice_site_NH_BOFO[i,3]    <- length(unique(Bioclimatic_variables_state_timeslice_site_NH_BOFO_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 8:21) {
            
            Bioclimatic_variables_state_timeslice_site_NH_BOFO[i,j-4] <- median(Bioclimatic_variables_state_timeslice_site_NH_BOFO_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_state_timeslice_site_NH_BOFO[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_state_timeslice_site_NH_BOFO[i,3]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_state_timeslice_site_NH_BOFO[i,4:17] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # SAVA
    {
      # Define the range of time slices to analyze (from 0 to 21, in increments of 0.5)
      list_timeslices <- seq(0,21,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_state_timeslice_site_NH_SAVA           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 17, byrow=TRUE))
      colnames(Bioclimatic_variables_state_timeslice_site_NH_SAVA) <- c("Biome", "Timeslice", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_state_timeslice_site_NH_SAVA$Biome <- 'SAVA'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_state_timeslice_site_NH_SAVA_subset <- subset(Simulated_bioclimatic_variable_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_SAVA & Timeslice == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_state_timeslice_site_NH_SAVA_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_state_timeslice_site_NH_SAVA[i,2]    <- unique(Bioclimatic_variables_state_timeslice_site_NH_SAVA_subset$Timeslice)
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_state_timeslice_site_NH_SAVA[i,3]    <- length(unique(Bioclimatic_variables_state_timeslice_site_NH_SAVA_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 8:21) {
            
            Bioclimatic_variables_state_timeslice_site_NH_SAVA[i,j-4] <- median(Bioclimatic_variables_state_timeslice_site_NH_SAVA_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_state_timeslice_site_NH_SAVA[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_state_timeslice_site_NH_SAVA[i,3]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_state_timeslice_site_NH_SAVA[i,4:17] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # STEP
    {
      # Define the range of time slices to analyze (from 0 to 21, in increments of 0.5)
      list_timeslices <- seq(0,21,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_state_timeslice_site_NH_STEP           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 17, byrow=TRUE))
      colnames(Bioclimatic_variables_state_timeslice_site_NH_STEP) <- c("Biome", "Timeslice", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_state_timeslice_site_NH_STEP$Biome <- 'STEP'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_state_timeslice_site_NH_STEP_subset <- subset(Simulated_bioclimatic_variable_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_STEP & Timeslice == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_state_timeslice_site_NH_STEP_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_state_timeslice_site_NH_STEP[i,2]    <- unique(Bioclimatic_variables_state_timeslice_site_NH_STEP_subset$Timeslice)
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_state_timeslice_site_NH_STEP[i,3]    <- length(unique(Bioclimatic_variables_state_timeslice_site_NH_STEP_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 8:21) {
            
            Bioclimatic_variables_state_timeslice_site_NH_STEP[i,j-4] <- median(Bioclimatic_variables_state_timeslice_site_NH_STEP_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_state_timeslice_site_NH_STEP[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_state_timeslice_site_NH_STEP[i,3]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_state_timeslice_site_NH_STEP[i,4:17] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # DESE
    {
      # Define the range of time slices to analyze (from 0 to 21, in increments of 0.5)
      list_timeslices <- seq(0,21,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_state_timeslice_site_NH_DESE           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 17, byrow=TRUE))
      colnames(Bioclimatic_variables_state_timeslice_site_NH_DESE) <- c("Biome", "Timeslice", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_state_timeslice_site_NH_DESE$Biome <- 'DESE'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_state_timeslice_site_NH_DESE_subset <- subset(Simulated_bioclimatic_variable_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_DESE & Timeslice == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_state_timeslice_site_NH_DESE_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_state_timeslice_site_NH_DESE[i,2]    <- unique(Bioclimatic_variables_state_timeslice_site_NH_DESE_subset$Timeslice)
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_state_timeslice_site_NH_DESE[i,3]    <- length(unique(Bioclimatic_variables_state_timeslice_site_NH_DESE_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 8:21) {
            
            Bioclimatic_variables_state_timeslice_site_NH_DESE[i,j-4] <- median(Bioclimatic_variables_state_timeslice_site_NH_DESE_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_state_timeslice_site_NH_DESE[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_state_timeslice_site_NH_DESE[i,3]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_state_timeslice_site_NH_DESE[i,4:17] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # TUND
    {
      # Define the range of time slices to analyze (from 0 to 21, in increments of 0.5)
      list_timeslices <- seq(0,21,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_state_timeslice_site_NH_TUND           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 17, byrow=TRUE))
      colnames(Bioclimatic_variables_state_timeslice_site_NH_TUND) <- c("Biome", "Timeslice", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_state_timeslice_site_NH_TUND$Biome <- 'TUND'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_state_timeslice_site_NH_TUND_subset <- subset(Simulated_bioclimatic_variable_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_TUND & Timeslice == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_state_timeslice_site_NH_TUND_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_state_timeslice_site_NH_TUND[i,2]    <- unique(Bioclimatic_variables_state_timeslice_site_NH_TUND_subset$Timeslice)
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_state_timeslice_site_NH_TUND[i,3]    <- length(unique(Bioclimatic_variables_state_timeslice_site_NH_TUND_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 8:21) {
            
            Bioclimatic_variables_state_timeslice_site_NH_TUND[i,j-4] <- median(Bioclimatic_variables_state_timeslice_site_NH_TUND_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_state_timeslice_site_NH_TUND[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_state_timeslice_site_NH_TUND[i,3]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_state_timeslice_site_NH_TUND[i,4:17] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # NH
    {
      # Combine all biome-specific data frames into one for the Northern Hemisphere
      Bioclimatic_variables_state_timeslice_site_NH_all  <- rbind(Bioclimatic_variables_state_timeslice_site_NH_TRFO, Bioclimatic_variables_state_timeslice_site_NH_WTFO, 
                                                                  Bioclimatic_variables_state_timeslice_site_NH_TEFO, Bioclimatic_variables_state_timeslice_site_NH_BOFO, 
                                                                  Bioclimatic_variables_state_timeslice_site_NH_SAVA, Bioclimatic_variables_state_timeslice_site_NH_STEP, 
                                                                  Bioclimatic_variables_state_timeslice_site_NH_DESE, Bioclimatic_variables_state_timeslice_site_NH_TUND)
      
      # Define the range of time slices to analyze (from 0 to 21, in increments of 0.5)
      list_timeslices <- seq(0,21,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_state_timeslice_site_NH_NH           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 17, byrow=TRUE))
      colnames(Bioclimatic_variables_state_timeslice_site_NH_NH) <- c("Biome", "Timeslice", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_state_timeslice_site_NH_NH$Biome <- 'Northern Hemisphere'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_state_timeslice_site_NH_NH_subset <- subset(Bioclimatic_variables_state_timeslice_site_NH_all, Timeslice == list_timeslices[i])
        
        # Replace NA values with 0 in the subset to handle missing data
        Bioclimatic_variables_state_timeslice_site_NH_NH_subset[is.na(Bioclimatic_variables_state_timeslice_site_NH_NH_subset)] <- 0
        
        # Assign the current time slice to the "Timeslice" column
        Bioclimatic_variables_state_timeslice_site_NH_NH[i,2]    <- unique(Bioclimatic_variables_state_timeslice_site_NH_NH_subset$Timeslice)
        
        # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
        Bioclimatic_variables_state_timeslice_site_NH_NH[i,3]    <- sum(Bioclimatic_variables_state_timeslice_site_NH_NH_subset$Site_number)
        
        # Compute the bioclimatic variable weighted by the proportion of each biome
        for (j in 4:17) {
          
          Bioclimatic_variables_state_timeslice_site_NH_NH[i,j] <- (Bioclimatic_variables_state_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_state_timeslice_site_NH_NH_subset$Biome == 'TRFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==1] + 
                                                                      Bioclimatic_variables_state_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_state_timeslice_site_NH_NH_subset$Biome == 'WTFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==2] +
                                                                      Bioclimatic_variables_state_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_state_timeslice_site_NH_NH_subset$Biome == 'TEFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==3] + 
                                                                      Bioclimatic_variables_state_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_state_timeslice_site_NH_NH_subset$Biome == 'BOFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==4] +
                                                                      Bioclimatic_variables_state_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_state_timeslice_site_NH_NH_subset$Biome == 'SAVA']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==5] +
                                                                      Bioclimatic_variables_state_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_state_timeslice_site_NH_NH_subset$Biome == 'STEP']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==6] +
                                                                      Bioclimatic_variables_state_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_state_timeslice_site_NH_NH_subset$Biome == 'DESE']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==7] +
                                                                      Bioclimatic_variables_state_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_state_timeslice_site_NH_NH_subset$Biome == 'TUND']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==8] )
          
        }
        
      }
      
    }
    
    
    # Combine result dataset
    Bioclimatic_variables_state_timeslice_site_NH_final <- rbind(Bioclimatic_variables_state_timeslice_site_NH_NH, Bioclimatic_variables_state_timeslice_site_NH_all)
    
    # Save the final results as a CSV file
    write.csv(Bioclimatic_variables_state_timeslice_site_NH_final, file="result/Statistical analysis 3_result.3-Temporal changes of simulated bioclimatic variables in the Northern Hemisphere over the last 21,000 years by modern biome distributions.csv", row.names=FALSE)
    
  }
  
  
  # -------------------------------------------------------------------------------------------------
  
  # (2) Change rate of simulated bioclimatic variables
  {
    # TRFO
    {
      # Define the range of time slices to analyze
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 18, byrow=TRUE))
      colnames(Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO) <- c("Biome", "Timeslice_younger", "Timeslice_older", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO$Biome <- 'TRFO'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO_subset <- subset(Simulated_bioclimatic_variable_change_rate_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_TRFO & Timeslice_younger == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO[i,2]    <- list_timeslices[i]
          Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO[i,3]    <- list_timeslices[i]+0.5
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO[i,4]    <- length(unique(Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 9:22) {
            
            Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO[i,j-4] <- median(Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO[i,3]    <- list_timeslices[i]+0.5
          Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO[i,4]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO[i,5:18] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # WTFO
    {
      # Define the range of time slices to analyze
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 18, byrow=TRUE))
      colnames(Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO) <- c("Biome", "Timeslice_younger", "Timeslice_older", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO$Biome <- 'WTFO'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO_subset <- subset(Simulated_bioclimatic_variable_change_rate_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_WTFO & Timeslice_younger == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO[i,2]    <- list_timeslices[i]
          Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO[i,3]    <- list_timeslices[i]+0.5
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO[i,4]    <- length(unique(Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 9:22) {
            
            Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO[i,j-4] <- median(Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO[i,3]    <- list_timeslices[i]+0.5
          Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO[i,4]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO[i,5:18] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # TEFO
    {
      # Define the range of time slices to analyze
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 18, byrow=TRUE))
      colnames(Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO) <- c("Biome", "Timeslice_younger", "Timeslice_older", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO$Biome <- 'TEFO'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO_subset <- subset(Simulated_bioclimatic_variable_change_rate_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_TEFO & Timeslice_younger == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO[i,2]    <- list_timeslices[i]
          Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO[i,3]    <- list_timeslices[i]+0.5
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO[i,4]    <- length(unique(Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 9:22) {
            
            Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO[i,j-4] <- median(Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO[i,3]    <- list_timeslices[i]+0.5
          Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO[i,4]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO[i,5:18] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # BOFO
    {
      # Define the range of time slices to analyze
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 18, byrow=TRUE))
      colnames(Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO) <- c("Biome", "Timeslice_younger", "Timeslice_older", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO$Biome <- 'BOFO'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO_subset <- subset(Simulated_bioclimatic_variable_change_rate_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_BOFO & Timeslice_younger == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO[i,2]    <- list_timeslices[i]
          Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO[i,3]    <- list_timeslices[i]+0.5
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO[i,4]    <- length(unique(Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 9:22) {
            
            Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO[i,j-4] <- median(Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO[i,3]    <- list_timeslices[i]+0.5
          Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO[i,4]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO[i,5:18] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # SAVA
    {
      # Define the range of time slices to analyze
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 18, byrow=TRUE))
      colnames(Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA) <- c("Biome", "Timeslice_younger", "Timeslice_older", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA$Biome <- 'SAVA'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA_subset <- subset(Simulated_bioclimatic_variable_change_rate_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_SAVA & Timeslice_younger == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA[i,2]    <- list_timeslices[i]
          Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA[i,3]    <- list_timeslices[i]+0.5
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA[i,4]    <- length(unique(Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 9:22) {
            
            Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA[i,j-4] <- median(Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA[i,3]    <- list_timeslices[i]+0.5
          Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA[i,4]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA[i,5:18] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # STEP
    {
      # Define the range of time slices to analyze
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_change_rate_timeslice_site_NH_STEP           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 18, byrow=TRUE))
      colnames(Bioclimatic_variables_change_rate_timeslice_site_NH_STEP) <- c("Biome", "Timeslice_younger", "Timeslice_older", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_change_rate_timeslice_site_NH_STEP$Biome <- 'STEP'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_change_rate_timeslice_site_NH_STEP_subset <- subset(Simulated_bioclimatic_variable_change_rate_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_STEP & Timeslice_younger == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_change_rate_timeslice_site_NH_STEP_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_STEP[i,2]    <- list_timeslices[i]
          Bioclimatic_variables_change_rate_timeslice_site_NH_STEP[i,3]    <- list_timeslices[i]+0.5
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_STEP[i,4]    <- length(unique(Bioclimatic_variables_change_rate_timeslice_site_NH_STEP_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 9:22) {
            
            Bioclimatic_variables_change_rate_timeslice_site_NH_STEP[i,j-4] <- median(Bioclimatic_variables_change_rate_timeslice_site_NH_STEP_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_change_rate_timeslice_site_NH_STEP[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_change_rate_timeslice_site_NH_STEP[i,3]    <- list_timeslices[i]+0.5
          Bioclimatic_variables_change_rate_timeslice_site_NH_STEP[i,4]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_change_rate_timeslice_site_NH_STEP[i,5:18] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # DESE
    {
      # Define the range of time slices to analyze
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_change_rate_timeslice_site_NH_DESE           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 18, byrow=TRUE))
      colnames(Bioclimatic_variables_change_rate_timeslice_site_NH_DESE) <- c("Biome", "Timeslice_younger", "Timeslice_older", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_change_rate_timeslice_site_NH_DESE$Biome <- 'DESE'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_change_rate_timeslice_site_NH_DESE_subset <- subset(Simulated_bioclimatic_variable_change_rate_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_DESE & Timeslice_younger == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_change_rate_timeslice_site_NH_DESE_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_DESE[i,2]    <- list_timeslices[i]
          Bioclimatic_variables_change_rate_timeslice_site_NH_DESE[i,3]    <- list_timeslices[i]+0.5
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_DESE[i,4]    <- length(unique(Bioclimatic_variables_change_rate_timeslice_site_NH_DESE_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 9:22) {
            
            Bioclimatic_variables_change_rate_timeslice_site_NH_DESE[i,j-4] <- median(Bioclimatic_variables_change_rate_timeslice_site_NH_DESE_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_change_rate_timeslice_site_NH_DESE[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_change_rate_timeslice_site_NH_DESE[i,3]    <- list_timeslices[i]+0.5
          Bioclimatic_variables_change_rate_timeslice_site_NH_DESE[i,4]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_change_rate_timeslice_site_NH_DESE[i,5:18] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # TUND
    {
      # Define the range of time slices to analyze
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_change_rate_timeslice_site_NH_TUND           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 18, byrow=TRUE))
      colnames(Bioclimatic_variables_change_rate_timeslice_site_NH_TUND) <- c("Biome", "Timeslice_younger", "Timeslice_older", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_change_rate_timeslice_site_NH_TUND$Biome <- 'TUND'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_change_rate_timeslice_site_NH_TUND_subset <- subset(Simulated_bioclimatic_variable_change_rate_timeslice_site_NH, Dataset_ID %in% DatasetID_list_NH_TUND & Timeslice_younger == list_timeslices[i])
        
        # Check if the subset contains any data
        if(nrow(Bioclimatic_variables_change_rate_timeslice_site_NH_TUND_subset) > 0) {
          
          # Assign the current time slice to the "Timeslice" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_TUND[i,2]    <- list_timeslices[i]
          Bioclimatic_variables_change_rate_timeslice_site_NH_TUND[i,3]    <- list_timeslices[i]+0.5
          
          # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
          Bioclimatic_variables_change_rate_timeslice_site_NH_TUND[i,4]    <- length(unique(Bioclimatic_variables_change_rate_timeslice_site_NH_TUND_subset$Dataset_ID))
          
          # Calculate the median for each bioclimatic variable (BIO1 to BIO19)
          for (j in 9:22) {
            
            Bioclimatic_variables_change_rate_timeslice_site_NH_TUND[i,j-4] <- median(Bioclimatic_variables_change_rate_timeslice_site_NH_TUND_subset[ ,j])
            
          }
          
        } else{
          
          # If no data is available for the current time slice, assign default values
          Bioclimatic_variables_change_rate_timeslice_site_NH_TUND[i,2]    <- list_timeslices[i] # Set the Timeslice
          Bioclimatic_variables_change_rate_timeslice_site_NH_TUND[i,3]    <- list_timeslices[i]+0.5
          Bioclimatic_variables_change_rate_timeslice_site_NH_TUND[i,4]    <- 0 # Set Site_number to 0
          Bioclimatic_variables_change_rate_timeslice_site_NH_TUND[i,5:18] <- NA # Set bioclimatic variables to NA
        }
        
      }
      
    }  
    
    # NH
    {
      # Combine all biome-specific data frames into one for the Northern Hemisphere
      Bioclimatic_variables_change_rate_timeslice_site_NH_all  <- rbind(Bioclimatic_variables_change_rate_timeslice_site_NH_TRFO, Bioclimatic_variables_change_rate_timeslice_site_NH_WTFO, 
                                                                        Bioclimatic_variables_change_rate_timeslice_site_NH_TEFO, Bioclimatic_variables_change_rate_timeslice_site_NH_BOFO, 
                                                                        Bioclimatic_variables_change_rate_timeslice_site_NH_SAVA, Bioclimatic_variables_change_rate_timeslice_site_NH_STEP, 
                                                                        Bioclimatic_variables_change_rate_timeslice_site_NH_DESE, Bioclimatic_variables_change_rate_timeslice_site_NH_TUND)
      
      # Define the range of time slices to analyze (from 0 to 21, in increments of 0.5)
      list_timeslices <- seq(0,20.5,by=0.5) # seq (starting age, ending age, step)
      
      # Initialize an empty data frame to store the results
      Bioclimatic_variables_change_rate_timeslice_site_NH_NH           <- data.frame(matrix(NA, nrow=length(list_timeslices), ncol= 18, byrow=TRUE))
      colnames(Bioclimatic_variables_change_rate_timeslice_site_NH_NH) <- c("Biome", "Timeslice_younger", "Timeslice_older", "Site_number", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the biome name 
      Bioclimatic_variables_change_rate_timeslice_site_NH_NH$Biome <- 'Northern Hemisphere'
      
      # Iterate through each time slice
      for (i in 1:length(list_timeslices)) {
        
        # Print progress for the current time slice
        print(paste0("+++++ Timeslice ", list_timeslices[i], " ka", " | (",i,"/",length(list_timeslices),") +++++"))
        
        # Filter the dataset for the current time slice and Dataset_IDs
        Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset <- subset(Bioclimatic_variables_change_rate_timeslice_site_NH_all, Timeslice_younger == list_timeslices[i])
        
        # Replace NA values with 0 in the subset to handle missing data
        Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset[is.na(Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset)] <- 0
        
        # Assign the current time slice to the "Timeslice" column
        Bioclimatic_variables_change_rate_timeslice_site_NH_NH[i,2]    <- list_timeslices[i]
        Bioclimatic_variables_change_rate_timeslice_site_NH_NH[i,3]    <- list_timeslices[i]+0.5
        
        # Count the number of unique sites (Dataset_IDs) and assign it to the "Site_number" column
        Bioclimatic_variables_change_rate_timeslice_site_NH_NH[i,4]    <- sum(Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset$Site_number)
        
        # Compute the bioclimatic variable weighted by the proportion of each biome
        for (j in 5:18) {
          
          Bioclimatic_variables_change_rate_timeslice_site_NH_NH[i,j] <- (Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset$Biome == 'TRFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==1] + 
                                                                            Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset$Biome == 'WTFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==2] +
                                                                            Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset$Biome == 'TEFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==3] + 
                                                                            Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset$Biome == 'BOFO']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==4] +
                                                                            Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset$Biome == 'SAVA']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==5] +
                                                                            Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset$Biome == 'STEP']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==6] +
                                                                            Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset$Biome == 'DESE']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==7] +
                                                                            Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset[ ,j][Bioclimatic_variables_change_rate_timeslice_site_NH_NH_subset$Biome == 'TUND']*Modern_potential_natural_biomes_proportion_NH$Proportion[Modern_potential_natural_biomes_proportion_NH$Biome ==8] )
          
        }
        
      }
      
    }
    
    
    # Combine result dataset
    Bioclimatic_variables_change_rate_timeslice_site_NH_final <- rbind(Bioclimatic_variables_change_rate_timeslice_site_NH_NH, Bioclimatic_variables_change_rate_timeslice_site_NH_all)
    
    # Save the final results as a CSV file
    write.csv(Bioclimatic_variables_change_rate_timeslice_site_NH_final, file="result/Statistical analysis 3_result.3-Temporal changes of change rates of simulated bioclimatic variables in the Northern Hemisphere over the last 21,000 years by modern biome distributions.csv", row.names=FALSE)
    
  }
  
}


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 4. Statistical analysis 4: Simulated climate and biome over the next 300 years in the Northern Hemisphere ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# 4.1 Calculation: biome shift per grid-cell
# -------------------------------------------------------------------------------------------------
# Load datasets
Biome_timeslice_MPI_grid_NH_PI_next300yr    <- read.csv2("data/Statistical analysis 4_data.1-Simluated biome since PI per grid-cell in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Biome_timeslice_MPI_grid_NH_PI_next300yr    <- type.convert(Biome_timeslice_MPI_grid_NH_PI_next300yr, as.is = TRUE) 


# calculation
{
  # (1) State of biome shift
  {
    # Extract the unique grid coordinates (Longitude and Latitude)
    coordinate_list <- unique(Biome_timeslice_MPI_grid_NH_PI_next300yr[ ,1:2])
    
    # Extract the unique time slices, selecting the third to fifth entries for analysis
    timeslice_list <- unique(Biome_timeslice_MPI_grid_NH_PI_next300yr$Timeslice)[3:5]
    
    # Initialize an empty data frame to store the results
    Biome_shift_state_MPI_PI_next300yr <- NULL
    
    # Loop through each grid cell (based on unique coordinates)
    for (i in 1:nrow(coordinate_list)) {
      
      # Print progress for the current grid cell
      print(paste0(" +++++ grid-cell ", i, "/", nrow(coordinate_list), " +++++"))
      
      # Subset the data for the current grid cell
      grid_subset <- subset(Biome_timeslice_MPI_grid_NH_PI_next300yr, Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2]) 
      
      # Initialize a result data frame for comparing biome states for the current grid cell
      biome_compared <- data.frame(matrix(NA, nrow=length(timeslice_list), ncol=7, byrow=TRUE))
      colnames(biome_compared) <- c("Longitude", "Latitude", "Timeslice_older", "Biome_older", "Timeslice_younger", "Biome_younger",  "Shifting")
      
      # Assign the longitude and latitude of the current grid cell to all rows
      biome_compared[ ,1:2] <- coordinate_list[i, ]
      # Assign the reference time slice (0.1) to the "Timeslice_older" column
      biome_compared$Timeslice_older  <- 0.1
      # Assign the biome state at the reference time slice (0.1) to the "Biome_older" column
      biome_compared$Biome_older  <- grid_subset$Biome[grid_subset$Timeslice == 0.1]
      
      # Loop through each time slice to compare biome states
      for (t in 1:length(timeslice_list)) {
        
        # Assign the current time slice to the "Timeslice_younger" column
        biome_compared[t, 5] <- timeslice_list[t]
        # Assign the biome state at the current time slice to the "Biome_younger" column
        biome_compared[t, 6] <- grid_subset$Biome[grid_subset$Timeslice == timeslice_list[t]]
        
        # Determine if the biome state has shifted compared to the reference time slice (0.1)
        biome_compared[t, 7] <- ifelse(grid_subset$Biome[grid_subset$Timeslice == timeslice_list[t]] == grid_subset$Biome[grid_subset$Timeslice == 0.1], "no", "yes")
        
      }
      
      # Append the results for the current grid cell to the final result data frame
      Biome_shift_state_MPI_PI_next300yr <- rbind(Biome_shift_state_MPI_PI_next300yr, biome_compared)
      
    }
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # (2) Magnitude of biome shifts
  {
    # Define the penalty matrix (Magnitude_matrix: defines the magnitude of biome shift between different biome types)
    Magnitude_matrix <- matrix(c(0,	1,	2,	3,	1,	2,	3,	4,
                                 1,	0,	1,	2,	2,	2,	3,	3,
                                 2,	1,	0,	1,	3,	2,	3,	2,
                                 3,	2,	1,	0,	4,	3,	2,	1,
                                 1,	2,	3,	4,	0,	1,	2,	4,
                                 2,	2,	2,	3,	1,	0,	1,	1,
                                 3,	3,	3,	2,	2,	1,	0,	1,
                                 4,	3,	2,	1,	4,	1,	2,	0), 
                               ncol = 8, nrow = 8, byrow = TRUE) 
    
    {
      # Add a new column "Magnitude" to the data frame to store shift magnitudes
      Biome_shift_state_MPI_PI_next300yr$Magnitude <- NA
      
      # Loop through each row of the data frame to calculate the shift magnitude
      for (i in 1:nrow(Biome_shift_state_MPI_PI_next300yr)) {
        
        # Print progress for the current row
        print(paste0(" +++++ row ", i, "/", nrow(Biome_shift_state_MPI_PI_next300yr), " +++++"))
        
        # Determine the biome shift magnitude using the penalty matrix
        Biome_shift_state_MPI_PI_next300yr[i,8] <- Magnitude_matrix[Biome_shift_state_MPI_PI_next300yr[i,4], Biome_shift_state_MPI_PI_next300yr[i,6]]
        
      }
      
    }
    
  }
  
  # Save the final results as a CSV file
  write.csv(Biome_shift_state_MPI_PI_next300yr, file="result/Statistical analysis 4_result.1-State and magnitude of biome shift over next 300 years relative to pre-industrial per grid in the Northern Hemisphere.csv", row.names=FALSE)
  
}  


# 4.2 Calculation: state and rate of change of bioclimatic variables at timeslice per grid-cell
# -------------------------------------------------------------------------------------------------
# Load datasets
Simulated_climate_timeslice_grid_NH_next300yr    <- read.csv2("data/Statistical analysis 4_data.2-Simluated monthly temperature and precipitation over next 300 years per grid-cell in the Northern Hemisphere.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Simulated_climate_timeslice_grid_NH_next300yr    <- type.convert(Simulated_climate_timeslice_grid_NH_next300yr, as.is = TRUE) 

# Calculation
{
  # (1) Simulated bioclimatic variables
  {
    # Extract unique grid-cell coordinates (Longitude and Latitude)
    coordinate_list <- unique(Simulated_climate_timeslice_grid_NH_next300yr[ ,1:2])
    
    # Extract unique timeslices from the data
    timeslice_list <- unique(Simulated_climate_timeslice_grid_NH_next300yr$Timeslice)
    
    # Initialize an empty data frame to store results
    bioclimatic_variables_state_timeslice_grid_NH_next300yr <- NULL
    
    # Loop through each timeslice to process data
    for(t in 1:length(timeslice_list)){
      
      # Create a data frame to store the results for the current timeslice
      bioclimatic_variables_state_timeslice_grid <- data.frame(matrix(NA, nrow=nrow(coordinate_list), ncol=17, byrow=TRUE))
      colnames(bioclimatic_variables_state_timeslice_grid) <- c("Longitude","Latitude", "Timeslice", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the coordinates to the first two columns
      bioclimatic_variables_state_timeslice_grid[ ,1:2] <- coordinate_list
      
      # Loop through each grid-cell in the coordinate list
      for (i in 1:nrow(coordinate_list)) {
        
        # Print progress to the console
        print(paste0("+++++ timeslice ", timeslice_list[t], " ka  ", t, "/", length(timeslice_list), " +++++ grid-cell ", i, "/", nrow(coordinate_list), " +++++"))
        
        # Subset data for the current timeslice and grid-cell
        timeslice_subset <- subset(Simulated_climate_timeslice_grid_NH_next300yr, Timeslice == timeslice_list[t] & Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2])
        
        # Assign the current timeslice to the Timeslice column
        bioclimatic_variables_state_timeslice_grid[i, 3] <- timeslice_list[t]
        
        # Separate Quarter_mean
        ## precipitation
        timeslice_subset$precipitation_1_2_3    <- mean(timeslice_subset$Precipitation[timeslice_subset$Month == 1] + timeslice_subset$Precipitation[timeslice_subset$Month == 2] + timeslice_subset$Precipitation[timeslice_subset$Month == 3])
        timeslice_subset$precipitation_2_3_4    <- mean(timeslice_subset$Precipitation[timeslice_subset$Month == 2] + timeslice_subset$Precipitation[timeslice_subset$Month == 3] + timeslice_subset$Precipitation[timeslice_subset$Month == 4])
        timeslice_subset$precipitation_3_4_5    <- mean(timeslice_subset$Precipitation[timeslice_subset$Month == 3] + timeslice_subset$Precipitation[timeslice_subset$Month == 4] + timeslice_subset$Precipitation[timeslice_subset$Month == 5])
        timeslice_subset$precipitation_4_5_6    <- mean(timeslice_subset$Precipitation[timeslice_subset$Month == 4] + timeslice_subset$Precipitation[timeslice_subset$Month == 5] + timeslice_subset$Precipitation[timeslice_subset$Month == 6])
        timeslice_subset$precipitation_5_6_7    <- mean(timeslice_subset$Precipitation[timeslice_subset$Month == 5] + timeslice_subset$Precipitation[timeslice_subset$Month == 6] + timeslice_subset$Precipitation[timeslice_subset$Month == 7])
        timeslice_subset$precipitation_6_7_8    <- mean(timeslice_subset$Precipitation[timeslice_subset$Month == 6] + timeslice_subset$Precipitation[timeslice_subset$Month == 7] + timeslice_subset$Precipitation[timeslice_subset$Month == 8])
        timeslice_subset$precipitation_7_8_9    <- mean(timeslice_subset$Precipitation[timeslice_subset$Month == 7] + timeslice_subset$Precipitation[timeslice_subset$Month == 8] + timeslice_subset$Precipitation[timeslice_subset$Month == 9])
        timeslice_subset$precipitation_8_9_10   <- mean(timeslice_subset$Precipitation[timeslice_subset$Month == 8] + timeslice_subset$Precipitation[timeslice_subset$Month == 9] + timeslice_subset$Precipitation[timeslice_subset$Month == 10])
        timeslice_subset$precipitation_9_10_11  <- mean(timeslice_subset$Precipitation[timeslice_subset$Month == 9] + timeslice_subset$Precipitation[timeslice_subset$Month == 10] + timeslice_subset$Precipitation[timeslice_subset$Month == 11])
        timeslice_subset$precipitation_10_11_12 <- mean(timeslice_subset$Precipitation[timeslice_subset$Month == 10] + timeslice_subset$Precipitation[timeslice_subset$Month == 11] + timeslice_subset$Precipitation[timeslice_subset$Month == 12])
        timeslice_subset$precipitation_11_12_1  <- mean(timeslice_subset$Precipitation[timeslice_subset$Month == 11] + timeslice_subset$Precipitation[timeslice_subset$Month == 12] + timeslice_subset$Precipitation[timeslice_subset$Month == 1])
        timeslice_subset$precipitation_12_1_2   <- mean(timeslice_subset$Precipitation[timeslice_subset$Month == 12] + timeslice_subset$Precipitation[timeslice_subset$Month == 1] + timeslice_subset$Precipitation[timeslice_subset$Month == 2])
        
        timeslice_subset$Wettest_Month    <- max(timeslice_subset$Precipitation)
        timeslice_subset$Driest_Month     <- min(timeslice_subset$Precipitation)
        
        timeslice_subset$Wettest_Quarter  <- apply(timeslice_subset[ ,which(names(timeslice_subset) %in% c("precipitation_1_2_3", "precipitation_2_3_4", "precipitation_3_4_5", "precipitation_4_5_6", "precipitation_5_6_7", "precipitation_6_7_8", "precipitation_7_8_9", "precipitation_8_9_10", "precipitation_9_10_11", "precipitation_10_11_12", "precipitation_11_12_1", "precipitation_12_1_2"))], 1, function(x) { names(x)[which.max(x)] })
        timeslice_subset$Driest_Quarter   <- apply(timeslice_subset[ ,which(names(timeslice_subset) %in% c("precipitation_1_2_3", "precipitation_2_3_4", "precipitation_3_4_5", "precipitation_4_5_6", "precipitation_5_6_7", "precipitation_6_7_8", "precipitation_7_8_9", "precipitation_8_9_10", "precipitation_9_10_11", "precipitation_10_11_12", "precipitation_11_12_1", "precipitation_12_1_2"))], 1, function(x) { names(x)[which.min(x)] })
        
        ## temperature
        timeslice_subset$temperature_1_2_3    <- mean(timeslice_subset$Temperature[timeslice_subset$Month == 1] + timeslice_subset$Temperature[timeslice_subset$Month == 2] + timeslice_subset$Temperature[timeslice_subset$Month == 3])
        timeslice_subset$temperature_2_3_4    <- mean(timeslice_subset$Temperature[timeslice_subset$Month == 2] + timeslice_subset$Temperature[timeslice_subset$Month == 3] + timeslice_subset$Temperature[timeslice_subset$Month == 4])
        timeslice_subset$temperature_3_4_5    <- mean(timeslice_subset$Temperature[timeslice_subset$Month == 3] + timeslice_subset$Temperature[timeslice_subset$Month == 4] + timeslice_subset$Temperature[timeslice_subset$Month == 5])
        timeslice_subset$temperature_4_5_6    <- mean(timeslice_subset$Temperature[timeslice_subset$Month == 4] + timeslice_subset$Temperature[timeslice_subset$Month == 5] + timeslice_subset$Temperature[timeslice_subset$Month == 6])
        timeslice_subset$temperature_5_6_7    <- mean(timeslice_subset$Temperature[timeslice_subset$Month == 5] + timeslice_subset$Temperature[timeslice_subset$Month == 6] + timeslice_subset$Temperature[timeslice_subset$Month == 7])
        timeslice_subset$temperature_6_7_8    <- mean(timeslice_subset$Temperature[timeslice_subset$Month == 6] + timeslice_subset$Temperature[timeslice_subset$Month == 7] + timeslice_subset$Temperature[timeslice_subset$Month == 8])
        timeslice_subset$temperature_7_8_9    <- mean(timeslice_subset$Temperature[timeslice_subset$Month == 7] + timeslice_subset$Temperature[timeslice_subset$Month == 8] + timeslice_subset$Temperature[timeslice_subset$Month == 9])
        timeslice_subset$temperature_8_9_10   <- mean(timeslice_subset$Temperature[timeslice_subset$Month == 8] + timeslice_subset$Temperature[timeslice_subset$Month == 9] + timeslice_subset$Temperature[timeslice_subset$Month == 10])
        timeslice_subset$temperature_9_10_11  <- mean(timeslice_subset$Temperature[timeslice_subset$Month == 9] + timeslice_subset$Temperature[timeslice_subset$Month == 10] + timeslice_subset$Temperature[timeslice_subset$Month == 11])
        timeslice_subset$temperature_10_11_12 <- mean(timeslice_subset$Temperature[timeslice_subset$Month == 10] + timeslice_subset$Temperature[timeslice_subset$Month == 11] + timeslice_subset$Temperature[timeslice_subset$Month == 12])
        timeslice_subset$temperature_11_12_1  <- mean(timeslice_subset$Temperature[timeslice_subset$Month == 11] + timeslice_subset$Temperature[timeslice_subset$Month == 12] + timeslice_subset$Temperature[timeslice_subset$Month == 1])
        timeslice_subset$temperature_12_1_2   <- mean(timeslice_subset$Temperature[timeslice_subset$Month == 12] + timeslice_subset$Temperature[timeslice_subset$Month == 1] + timeslice_subset$Temperature[timeslice_subset$Month == 2])
        
        timeslice_subset$Warmest_Month   <- max(timeslice_subset$Temperature)
        timeslice_subset$Coldest_Month   <- min(timeslice_subset$Temperature)
        
        timeslice_subset$Warmest_Quarter   <- apply(timeslice_subset[ ,which(names(timeslice_subset) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.max(x)] })
        timeslice_subset$Coldest_Quarter   <- apply(timeslice_subset[ ,which(names(timeslice_subset) %in% c("temperature_1_2_3", "temperature_2_3_4", "temperature_3_4_5", "temperature_4_5_6", "temperature_5_6_7", "temperature_6_7_8", "temperature_7_8_9", "temperature_8_9_10", "temperature_9_10_11", "temperature_10_11_12", "temperature_11_12_1", "temperature_12_1_2"))], 1, function(x) { names(x)[which.min(x)] })
        
        # Calculate BIO
        ## BIO 1 = Annual Mean Temperature
        bioclimatic_variables_state_timeslice_grid[i, 4] <- mean(timeslice_subset$Temperature)
        
        ## BIO 4 = Temperature Seasonality (standard deviation ×100)
        bioclimatic_variables_state_timeslice_grid[i, 5] <- sd(timeslice_subset$Temperature)*100
        
        ## BIO 8 = Mean Temperature of Wettest Quarter
        bioclimatic_variables_state_timeslice_grid[i, 6] <- (timeslice_subset$Temperature[timeslice_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(timeslice_subset$Wettest_Quarter)))] +
                                                               timeslice_subset$Temperature[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(timeslice_subset$Wettest_Quarter)))] +  
                                                               timeslice_subset$Temperature[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(timeslice_subset$Wettest_Quarter)))])/3
        
        ## BIO 9 = Mean Temperature of Driest Quarter
        bioclimatic_variables_state_timeslice_grid[i, 7] <- (timeslice_subset$Temperature[timeslice_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(timeslice_subset$Driest_Quarter)))] +
                                                               timeslice_subset$Temperature[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(timeslice_subset$Driest_Quarter)))] +  
                                                               timeslice_subset$Temperature[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(timeslice_subset$Driest_Quarter)))])/3
        
        ## BIO10 = Mean Temperature of Warmest Quarter
        bioclimatic_variables_state_timeslice_grid[i, 8] <- (timeslice_subset$Temperature[timeslice_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(timeslice_subset$Warmest_Quarter)))] +
                                                               timeslice_subset$Temperature[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(timeslice_subset$Warmest_Quarter)))] +  
                                                               timeslice_subset$Temperature[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(timeslice_subset$Warmest_Quarter)))])/3
        
        ## BIO11 = Mean Temperature of Coldest Quarter
        bioclimatic_variables_state_timeslice_grid[i, 9] <- (timeslice_subset$Temperature[timeslice_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(timeslice_subset$Coldest_Quarter)))] +
                                                               timeslice_subset$Temperature[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(timeslice_subset$Coldest_Quarter)))] +  
                                                               timeslice_subset$Temperature[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(timeslice_subset$Coldest_Quarter)))])/3
        
        ## BIO12 = Annual Precipitation
        bioclimatic_variables_state_timeslice_grid[i, 10] <- sum(timeslice_subset$Precipitation)
        
        ## BIO13 = Precipitation of Wettest Month
        bioclimatic_variables_state_timeslice_grid[i, 11] <- max(timeslice_subset$Precipitation)
        
        ## BIO14 = Precipitation of Driest Month
        bioclimatic_variables_state_timeslice_grid[i, 12] <- min(timeslice_subset$Precipitation)
        
        ## BIO15 = Precipitation Seasonality (Coefficient of Variation)
        bioclimatic_variables_state_timeslice_grid[i, 13] <- sd(timeslice_subset$Precipitation)/(1+ mean(timeslice_subset$Precipitation))*100
        
        ## BIO16 = Precipitation of Wettest Quarter
        bioclimatic_variables_state_timeslice_grid[i, 14] <- (timeslice_subset$Precipitation[timeslice_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(timeslice_subset$Wettest_Quarter)))] +
                                                                timeslice_subset$Precipitation[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(timeslice_subset$Wettest_Quarter)))] +  
                                                                timeslice_subset$Precipitation[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(timeslice_subset$Wettest_Quarter)))])
        
        ## BIO17 = Precipitation of Driest Quarter 
        bioclimatic_variables_state_timeslice_grid[i, 15] <- (timeslice_subset$Precipitation[timeslice_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(timeslice_subset$Driest_Quarter)))] +
                                                                timeslice_subset$Precipitation[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(timeslice_subset$Driest_Quarter)))] +  
                                                                timeslice_subset$Precipitation[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(timeslice_subset$Driest_Quarter)))])
        
        ## BIO18 = Precipitation of Warmest Quarter
        bioclimatic_variables_state_timeslice_grid[i, 16] <- (timeslice_subset$Precipitation[timeslice_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(timeslice_subset$Warmest_Quarter)))] +
                                                                timeslice_subset$Precipitation[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(timeslice_subset$Warmest_Quarter)))] +  
                                                                timeslice_subset$Precipitation[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(timeslice_subset$Warmest_Quarter)))])
        
        ## BIO19 = Precipitation of Coldest Quarter
        bioclimatic_variables_state_timeslice_grid[i, 17] <- (timeslice_subset$Precipitation[timeslice_subset$Month == as.numeric(gsub("\\D*(\\d).*", "\\1", unique(timeslice_subset$Coldest_Quarter)))] +
                                                                timeslice_subset$Precipitation[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)_.*", "\\1", unique(timeslice_subset$Coldest_Quarter)))] +  
                                                                timeslice_subset$Precipitation[timeslice_subset$Month == as.numeric(gsub(".*_(\\d+)$", "\\1", unique(timeslice_subset$Coldest_Quarter)))])
        
        
        
      }
      
      # Append the result to the main result container
      bioclimatic_variables_state_timeslice_grid_NH_next300yr <- rbind(bioclimatic_variables_state_timeslice_grid_NH_next300yr, bioclimatic_variables_state_timeslice_grid)
      
    }
    
    # save data frame to csv:
    write.csv(bioclimatic_variables_state_timeslice_grid_NH_next300yr, file="result/Statistical analysis 4_result.3-Simulated bioclimatic variables at timeslice over next 300 years per grid-cell in the Northern Hemisphere.csv", row.names=FALSE)
    
  }
  
  
  # -------------------------------------------------------------------------------------------------
  
  # (2) Change rate of simulated bioclimatic variables
  {
    # Extract unique grid-cell coordinates (Longitude and Latitude)
    coordinate_list <- unique(bioclimatic_variables_state_timeslice_grid_NH_next300yr[ ,1:2])
    
    # Initialize an empty data frame to store results
    bioclimatic_variables_change_rate_timeslice_grid_NH_PI_next300yr <- NULL
    
    # Loop through each grid-cell in the coordinate list
    for (i in 1:nrow(coordinate_list)) {
      
      # Print progress for the current grid-cell
      print(paste0(" +++++ grid-cell ", i, "/", nrow(coordinate_list), " +++++"))
      
      # Subset data for the current grid-cell
      grid_subset <- subset(bioclimatic_variables_state_timeslice_grid_NH_next300yr, Longitude == coordinate_list[i, 1] & Latitude == coordinate_list[i, 2])
      
      # Set up a result matrix for storing change rates for the current grid-cell
      bioclimatic_variables_change_rate_timeslice_grid <- data.frame(matrix(NA, nrow=nrow(grid_subset)-1, ncol=18, byrow=TRUE))
      colnames(bioclimatic_variables_change_rate_timeslice_grid) <- c("Longitude", "Latitude", "Timeslice_younger", "Timeslice_older", "BIO1", "BIO4", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19")
      
      # Assign the coordinates to the result matrix
      bioclimatic_variables_change_rate_timeslice_grid[ ,1:2] <- coordinate_list[i,]
      bioclimatic_variables_change_rate_timeslice_grid$Timeslice_older <- 0.1
      
      # Loop through the rows of the grid-cell subset to calculate change rates 
      for(j in 2:nrow(grid_subset)){ 
        
        # Assign the younger timeslice to the result matrix
        bioclimatic_variables_change_rate_timeslice_grid[j-1,3] <- grid_subset[j,3]
        
        # Calculate the change rate for each bioclimatic variable
        for(k in 4:17){ 
          
          # Calculate the percentage change relative to the oldest timeslice (0.1)
          bioclimatic_variables_change_rate_timeslice_grid[j-1,k+1] <- abs((grid_subset[j,k] - grid_subset[grid_subset$Timeslice == 0.1,k])/grid_subset[grid_subset$Timeslice == 0.1,k])*100
          
        }
        
      }
      
      # Append the result for the current grid-cell to the final data frame
      bioclimatic_variables_change_rate_timeslice_grid_NH_PI_next300yr <- rbind(bioclimatic_variables_change_rate_timeslice_grid_NH_PI_next300yr, bioclimatic_variables_change_rate_timeslice_grid)
      
    }
    
    # save data frame to csv:
    write.csv(bioclimatic_variables_change_rate_timeslice_grid_NH_PI_next300yr, file="result/Statistical analysis 4_result.4-Change rate of simulated bioclimatic variables at timeslice over next 300 years per grid-cell in the Northern Hemisphere.csv", row.names=FALSE)
    
  }
  
  
}


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
## 5. Statistical analysis 5: Gaussian kernel correlation ##
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# Example in Gaussian kernel correlations between proportion of reconstructed biome shifts and bioclimatic variables over the last 21,000 years in the Northern Hemisphere

# Load datasets
Biome_shift_proportion_temporal_pollen_NH             <- read.csv2("data/Statistical analysis 5_data.1-Temporal changes of proportion of biome shifts in the Northern Hemisphere over the last 21,000 years.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')
Bioclimatic_variables_change_rate_temporal_pollen_NH  <- read.csv2("data/Statistical analysis 5_data.2-Temporal changes of change rates of simulated bioclimatic variables in the Northern Hemisphere over the last 21,000 years.csv", stringsAsFactors=FALSE, header=TRUE, sep = ',')

# Data Type Conversion
Biome_shift_proportion_temporal_pollen_NH             <- type.convert(Biome_shift_proportion_temporal_pollen_NH, as.is = TRUE) 
Bioclimatic_variables_change_rate_temporal_pollen_NH  <- type.convert(Bioclimatic_variables_change_rate_temporal_pollen_NH, as.is = TRUE) 

# Combine datasets
Biome_shift_Bioclimatic_variables_change_temporal_pollen_NH <- full_join(Biome_shift_proportion_temporal_pollen_NH[ ,c(2:3,5)], Bioclimatic_variables_change_rate_temporal_pollen_NH[ ,c(2:3,5:18)], by = c("Timeslice_younger", "Timeslice_older"))

# calculation
{
  # Extract the list of bioclimatic variables
  BIO_list <- names(Biome_shift_Bioclimatic_variables_change_temporal_pollen_NH)[4:17]
  
  # Initialize a data frame to store correlation results
  Gaussian_kernel_correlations_pollen_NH <- data.frame(matrix(NA, nrow=length(BIO_list), ncol=4, byrow=TRUE))
  
  colnames(Gaussian_kernel_correlations_pollen_NH) <- c("Correlation", "BIO", "P-value", "R-value")
  
  # Define the correlation description
  Gaussian_kernel_correlations_pollen_NH$Correlation <- "Proportion of biome shifts in reconstruction vs. Bioclimatic variables in MPI-ESM"
  Gaussian_kernel_correlations_pollen_NH$BIO         <- BIO_list
  
  # Loop through each bioclimatic variable
  for (i in 1:length(BIO_list)) {
      
    # Prepare time-series data for the biome shifts (x) and the bioclimatic variable (y)
      x <- zoo(Biome_shift_Bioclimatic_variables_change_temporal_pollen_NH[ ,3], order.by= Biome_shift_Bioclimatic_variables_change_temporal_pollen_NH$Timeslice_younger)
      y <- zoo(Biome_shift_Bioclimatic_variables_change_temporal_pollen_NH[ ,i+3], order.by= Biome_shift_Bioclimatic_variables_change_temporal_pollen_NH$Timeslice_younger)
      
      # Compute the correlation coefficient (r-value) and p-value using the Gaussian kernel method
      Gaussian_kernel_correlations_pollen_NH[i, 3] <- nexcf_ci(x,y)$rxy  ## Correlation coefficient
      Gaussian_kernel_correlations_pollen_NH[i, 4] <- nexcf_ci(x,y)$pval ## P-value
      
  }
  
  # save data frame to csv:
  write.csv(Gaussian_kernel_correlations_pollen_NH, file="result/Statistical analysis 5_result.1-Gaussian kernel correlations between proportion of reconstructed biome shifts and bioclimatic variables in the Northern Hemisphere.csv", row.names=FALSE)
  
} 


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# END IN HERE #
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

