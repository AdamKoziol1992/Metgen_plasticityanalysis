###Auxillary functions used for the hill_numbers.Rmd file
###Author: Adam Koziol

#############Takes the dataframe and joins the metadata to create a useful csv from it
make_overall_csv <- function(first_df){
  first_df %>%
    mutate(individual_id = as.character(individual_id)) %>%
    left_join(metadata, by = c('individual_id')) %>%
    unique %>%
    mutate(beta = 1 - as.numeric(beta.dissimilarity)) %>%
    mutate(across(c(individual_id, species, cage), as.character),
           across(c(beta, beta.dissimilarity), as.numeric))}

############Filters individuals that do not have all data points
which_ind_2keep <- function(metadata, species){metadata %>%
    group_by(individual_id) %>%
    tally %>%
    filter(n == 5) %>%
    filter(str_detect(individual_id, species)) %>%
    left_join(metadata, by = 'individual_id')}

############Plotting function to visualise the overall dissimilarities
plot_overall_betas <- function(input_df){
  Overall_beta <- ggplot(input_df, aes(x=species, y=beta.dissimilarity, fill = species)) +
    geom_boxplot() +
    geom_jitter(alpha=0.2, position=position_jitter(0.2), size=2) +
    scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1")) +
    ylab('Overall dissimilarity') +
    xlab('Species') +
    labs(fill = 'Species') +
    theme(legend.title = element_text(size = 14),
          legend.text = element_text(size=13),
          panel.background = element_blank()) +
    ylim(0,1) +
    ggtitle(deparse(substitute(input_df))) +
    stat_compare_means(method = 'wilcox.test')
}

#############Creates the Overall dissimilarity dataframe for later downstream use
Create_dissims <- function(input_df, input_metadata, input_list, stored_data, hillR_func, extra_input, extra_file){
  input_df <- as.matrix(input_df)
  stored_data <- data.frame()
  if(extra_input == TRUE){
    for (m in input_list){
      print(m)
      PRE <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Day1', 'sample_id'])
      POST <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Acc', 'sample_id'])
      HEAT <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Heat', 'sample_id'])
      COLD <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Cold', 'sample_id'])
      DIET <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Diet', 'sample_id'])
      my_matrix <- input_df[,c(PRE,POST,HEAT,COLD,DIET)]
      my_beta <- hillR_func(t(my_matrix), extra_file, q = 1) %>%
        mutate(dissimilarity = 1-local_similarity) %>%
        select(dissimilarity)
      row <- c(individual_id=m, beta=my_beta)
      stored_data <- rbind(stored_data,row)
    }
    return(stored_data)
  }
  else{
    for (m in input_list){
      print(m)
      PRE <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Day1', 'sample_id'])
      POST <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Acc', 'sample_id'])
      HEAT <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Heat', 'sample_id'])
      COLD <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Cold', 'sample_id'])
      DIET <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Diet', 'sample_id'])
      my_matrix <- input_df[,c(PRE,POST,HEAT,COLD,DIET)]
      my_beta <- hillR_func(t(my_matrix), q = 1) %>%
        mutate(dissimilarity = 1-local_similarity) %>%
        select(dissimilarity)
      row <- c(individual_id=m, beta=my_beta)
      stored_data <- rbind(stored_data,row)
    }
    return(stored_data)
  }
}

###################Creates the dissimilarity table between time-point pairs
beta_div_func <- function(input_df, input_metadata, input_list, input_function, stored_data, extra_input, extra_file, phylo){
  stored_data <- data.frame()
  if(extra_input == TRUE & phylo == TRUE){
    for (m in input_list){
      print(m)
      Day1 <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Day1', 'sample_id'])
      Acc <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Acc', 'sample_id'])
      Heat <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Heat', 'sample_id'])
      Cold <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Cold', 'sample_id'])
      Diet <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Diet', 'sample_id'])
      
      pair1 <- input_df[,c(Day1, Acc)] %>%
        filter(rowSums(.) != 0)
      T1 <- extra_file %>%
        keep.tip(rownames(pair1))
      
      pair2 <- input_df[,c(Acc, Heat)] %>%
        filter(rowSums(.) != 0)
      T2 <- extra_file %>%
        keep.tip(rownames(pair2))
      
      pair3 <- input_df[,c(Heat, Cold)] %>%
        filter(rowSums(.) != 0)
      T3 <- extra_file %>%
        keep.tip(rownames(pair3))
      
      pair4 <- input_df[,c(Cold, Diet)] %>%
        filter(rowSums(.) != 0)
      T4 <- extra_file %>%
        keep.tip(rownames(pair4))
      
      Day1_Acc <- input_function(t(pair1), T1, q=1)
      Acc_Heat <- input_function(t(pair2), T2, q=1)
      Heat_Cold <- input_function(t(pair3), T3, q=1)
      Cold_Diet <- input_function(t(pair4), T4, q=1)
      stored_data <- rbind(stored_data,Day1_Acc=Day1_Acc, Acc_Heat=Acc_Heat, Heat_Cold=Heat_Cold, Cold_Diet=Cold_Diet)
    }
    stored_data <- stored_data %>%
      rownames_to_column('Treatment_pair') %>%
      mutate(individual_id = rep(input_list, each=4))
    return(stored_data)
  }
  else if(extra_input == TRUE & phylo == FALSE){
    for (m in input_list){
      print(m)
      Day1 <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Day1', 'sample_id'])
      Acc <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Acc', 'sample_id'])
      Heat <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Heat', 'sample_id'])
      Cold <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Cold', 'sample_id'])
      Diet <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Diet', 'sample_id'])
      
      pair1 <- input_df[,c(Day1, Acc)] %>%
        filter(rowSums(.) != 0)
      T1 <- extra_file %>%
        .[rownames(pair1),]
      
      pair2 <- input_df[,c(Acc, Heat)] %>%
        filter(rowSums(.) != 0)
      T2 <- extra_file %>%
        .[rownames(pair2),]
      
      pair3 <- input_df[,c(Heat, Cold)] %>%
        filter(rowSums(.) != 0)
      T3 <- extra_file %>%
        .[rownames(pair3),]
      
      pair4 <- input_df[,c(Cold, Diet)] %>%
        filter(rowSums(.) != 0)
      T4 <- extra_file %>%
        .[rownames(pair4),]
      
      Day1_Acc <- input_function(t(pair1), T1, q=1)
      Acc_Heat <- input_function(t(pair2), T2, q=1)
      Heat_Cold <- input_function(t(pair3), T3, q=1)
      Cold_Diet <- input_function(t(pair4), T4, q=1)
      stored_data <- rbind(stored_data,Day1_Acc=Day1_Acc, Acc_Heat=Acc_Heat, Heat_Cold=Heat_Cold, Cold_Diet=Cold_Diet)
    }
    stored_data <- stored_data %>%
      rownames_to_column('Treatment_pair') %>%
      mutate(individual_id = rep(input_list, each=4))
    return(stored_data) 
  }
  else{
    for (m in input_list){
      print(m)
      Day1 <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Day1', 'sample_id'])
      Acc <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Acc', 'sample_id'])
      Heat <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Heat', 'sample_id'])
      Cold <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Cold', 'sample_id'])
      Diet <- as.character(input_metadata[input_metadata$individual_id == m & input_metadata$treatment == 'Diet', 'sample_id'])
      pair1 <- input_df[,c(Day1, Acc)] %>%
        filter(rowSums(.) != 0)
      pair2 <- input_df[,c(Acc, Heat)] %>%
        filter(rowSums(.) != 0)
      pair3 <- input_df[,c(Heat, Cold)] %>%
        filter(rowSums(.) != 0)
      pair4 <- input_df[,c(Cold, Diet)] %>%
        filter(rowSums(.) != 0)
      Day1_Acc <- input_function(pair1, q=1)
      Acc_Heat <- input_function(pair2, q=1)
      Heat_Cold <- input_function(pair3, q=1)
      Cold_Diet <- input_function(pair4, q=1)
      stored_data <- rbind(stored_data,Day1_Acc=Day1_Acc, Acc_Heat=Acc_Heat, Heat_Cold=Heat_Cold, Cold_Diet=Cold_Diet)
    }
    stored_data <- stored_data %>%
      rownames_to_column('Treatment_pair') %>%
      mutate(individual_id = rep(input_list, each=4))
    return(stored_data)
  }
}
