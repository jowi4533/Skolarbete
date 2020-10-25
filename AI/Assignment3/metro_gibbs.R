#made by Jonas Wikstrom, Jonas Jons, Axel Johansson



metro_gibbs <-function(network,case){
  #possible diseases: Pneumonia, Tuberculosis, Lung Cancer and Bronchitis
  
  #set of variables:Pn, VTB, TB, XR, Te, Dy
  #unknown variables: TB, Dy, BR, LC , pn
  #known variables: Te, VTB, Sm, XR, Dy
  
  #set known variables values 
  val_Te = case[["Te"]]
  val_VTB = case[["VTB"]]
  val_SM = case[["Sm"]]
  val_XR = case[["XR"]]
  val_Dy = case[["Dy"]]
  
  #unknown variables
  old_Pn = rbinom(1, 1,.5)
  old_TB = rbinom(1, 1,.5)
  old_LC  = rbinom(1, 1,.5)
  old_BR = rbinom(1, 1,.5)
  
  current_sample <- c(val_Te,val_VTB,val_SM,val_XR,val_Dy,old_TB,old_BR,old_LC,old_Pn)
  p_old = j_dist(network,current_sample)
  samples_df = setNames(data.frame(matrix(ncol = 9, nrow = 0)), c("Te", "Vtb", "Sm","Xr","Dy","Tb", "Br", "Lc","Pn"))
  rand_list = runif(4000)
  counter = 1
  for (i in 1:1000){
    #new pn variable?  
    new_Pn = 1- old_Pn
    rand_no =  rand_list[counter]
    current_sample <- c(val_Te,val_VTB,val_SM,val_XR,val_Dy,old_TB,old_BR,old_LC,new_Pn)
    p_new = j_dist(network,current_sample)
    counter = counter+ 1
    if (p_new >= p_old){
      old_Pn <- new_Pn
      p_old <- p_new
    }
    else if (rand_no  <= (p_new/p_old)){
      old_Pn <- new_Pn
      p_old <- p_new
    }
    #new tb variable?
    new_TB = 1- old_TB
    rand_no =  rand_list[counter]
    counter = counter+ 1
    current_sample <- c(val_Te,val_VTB,val_SM,val_XR,val_Dy,new_TB,old_BR,old_LC,old_Pn)
    p_new =  j_dist(network,current_sample)
    if (p_new  > p_old){
      old_TB <- new_TB
      p_old <- p_new
      
    }
    else if(rand_no  <= (p_new/p_old)){
      old_TB <- new_TB
      p_old <- p_new
    }
    #new lc variable?
    new_LC = 1- old_LC
    rand_no =  rand_list[counter]
    counter = counter+ 1
    current_sample <- c(val_Te,val_VTB,val_SM,val_XR,val_Dy,old_TB,old_BR,new_LC,old_Pn)
    p_new =  j_dist(network,current_sample)
    if (p_new >= p_old){
      old_LC <- new_LC
      p_old <- p_new
    }
    else if(rand_no  <= (p_new/p_old)){
      old_LC <- new_LC
      p_old <- p_new
    }
    #new br?
    new_BR = 1- old_BR
    rand_no =  rand_list[counter]
    counter = counter+ 1
    current_sample <- c(val_Te,val_VTB,val_SM,val_XR,val_Dy,old_TB,new_BR,old_LC,old_Pn)
    p_new =  j_dist(network,current_sample)
    if(p_new >= p_old){
      old_BR <- new_BR
      p_old <- p_new
    }
    else if(rand_no  <= (p_new/p_old)){
      old_BR <- new_BR
      p_old <- p_new
    }
    
    sample <- c(val_Te,val_VTB,val_SM,val_XR,val_Dy,old_TB,old_BR,old_LC,old_Pn)
    samples_df[nrow(samples_df) + 1,] = sample
  }
  #remove burn in: 100
  
  samples_df = samples_df[-c(1:100),]
  
  #probablility of diseases
  total_cases = nrow(samples_df)
  
  #pnemonia
  pneumonia_cases = samples_df[samples_df$Pn == TRUE, "Pn"]
  pn_true = length(pneumonia_cases)/total_cases
  pn_false = 1 - pn_true
  
  #tb
  Tuberculosis_cases = samples_df[samples_df$Tb == TRUE, "Tb"]
  tb_true = length(Tuberculosis_cases)/total_cases
  tb_false = 1 - tb_true
  
  #lc
  Lung_Cancer_cases = samples_df[samples_df$Lc == TRUE, "Lc"]
  lc_true = length(Lung_Cancer_cases)/total_cases
  lc_false = 1 - lc_true
  
  #br
  Bronchitis_cases = samples_df[samples_df$Br == TRUE, "Br"]
  br_true = length(Bronchitis_cases)/total_cases
  br_false = 1 - br_true
  
  
  return(c(pn_true,tb_true,lc_true,br_true))
}

j_dist <- function(network,sample){
  # TE:1 , VTB:2 , SM:3, XR:4, DY:5, TB: 6, BR: 7, LC: 8, PN: 9 
  p_pn = network[[1]]["1",toString(sample[[9]])]
  p_sm = network[[5]]["1", toString(sample[[3]])]
  p_vtb = network[[3]]["1", toString(sample[[2]])]
  p_lc = network[[6]][toString(sample[[3]]),toString(sample[[8]])]
  p_br = network[[7]][toString(sample[[3]]),toString(sample[[7]])]
  p_dy = network[[9]][i_s(c(sample[[8]],sample[[7]])),toString(sample[[5]])]
  p_tb = network[[4]][toString(sample[[2]]),toString(sample[[6]])]
  p_xr = network[[8]][i_s(c(sample[[9]],sample[[6]],sample[[8]])),toString(sample[[4]])]
  p_te = getDensity(sample[[9]],sample[[1]])
  p = p_pn * p_sm * p_vtb * p_lc * p_br  * p_dy * p_tb * p_xr * p_te
  return(p)
}

#converts binary representation to string
i_s <- function(val){
  return(paste(val, collapse='' ))
}
#returns density of temperature given pn
getDensity <- function(pn, te) {
  #Calc mean and sd of temp given pn = 0 and pn = 1
  mean_te_pn0 = mean(hist$Te[hist$Pn == 0])
  mean_te_pn1 = mean(hist$Te[hist$Pn == 1])
  sd_te_pn0 = sd(hist$Te[hist$Pn == 0])
  sd_te_pn1 = sd(hist$Te[hist$Pn == 1])
  #Gets densities for temp given pn = 0 and pn = 1, this is a vector of 5 densities
  if(pn == 0) {
    te_pn = dnorm(te, mean_te_pn0, sd_te_pn0)
  }
  else if(pn == 1) {
    te_pn = dnorm(te, mean_te_pn1, sd_te_pn1)
  }
  return(te_pn)
}

diagnose <- function(network,cases){
  disease_mat <- matrix(nrow = 10, ncol = 4)
  for(i in 1:nrow(cases)){
    disease_mat[i,] <- metro_gibbs(network,cases[i,])
    
  }
  return(disease_mat)
}
learn = function(historical_data){
  
  #true_row <- c(1,1,1,1,1,1,1,1,1)
  #false_row <- c(1,1,1,1,1,1,1,1,1)
  #historical_data[nrow(historical_data) + 1,] = true_row
  #historical_data[nrow(historical_data) + 1,] = false_row
  
  # historical data is a R dataframe, starts at row 1, column 1
  # From historical data get the cases where Pneumonia is present
  # Second param says get only the values of that, ignore rest of dataframe
  #-------------------------------
  # Pn = Pneumonia, TB = Tuberculosis, LC = Lung Cancer, Br = Bronchitis
  # Sm = Smoker, XR = Xray result, Dy = Dyspnea, Te = Temperature, VTB = Visited a TB location
  # History of 10000 patients
  
  # Leaf nodes
  Pn_parents = NULL
  VTB_parents = NULL
  Sm_parents = NULL
  # Has parents, order is important
  Te_parents = "Pn"
  TB_parents = "VTB"
  LC_parents = "Sm"
  Br_parents = "Sm"
  Dy_parents = c("LC", "Br")
  XR_parents = c("Pn", "TB", "LC")
  parents_list = list(Pn_parents, Te_parents, VTB_parents, TB_parents, Sm_parents, LC_parents, Br_parents, XR_parents, Dy_parents)
  amount_of_parents = c()
  for(i in 1:length(parents_list)){
    amount_of_parents = c(amount_of_parents, length(parents_list[[i]]))
  }
  node_network = createNetwork(historical_data, amount_of_parents)
  node_network = trainNetwork(node_network, historical_data, parents_list)
  return(node_network)
}
createNetwork = function(historical_data, amount_of_parents){
  matrix_vector = list()
  for(node in 1:ncol(historical_data)){
    node_matrix = createNodeMatrix(amount_of_parents[node])
    matrix_vector[[node]] = node_matrix
  }
  return(matrix_vector)
}
trainNetwork = function(nodeNetwork, historical_data, parents_list){
  total_cases = nrow(historical_data)
  for (i in 1:total_cases){
    for (node in 1:length(nodeNetwork)){
      nodeNetwork[[node]] = calcProbabilities(nodeNetwork[[node]], historical_data[i, node], parents_list[[node]], historical_data, i)
    }
  }
  nodeNetwork = normalizeValues(nodeNetwork, total_cases)
  return(nodeNetwork)
}
calcProbabilities = function(node, case, parents, historical_data, caseIndex){
  # Case has values 0 or 1, if 0 put in column 1, 1 put in column 2
  
  
  if(case == 1 || case == 0){
    col = case + 1
    if(is.null(parents)){
      # No parents, simply update from data
      node[1, col] = node[1, col] + 1
    } else {
      parent_values = c()
      for(parent in parents){
        column = which(colnames(historical_data) == parent)
        row = caseIndex
        parent_values = c(parent_values, historical_data[row,column])
      }
      binaryEncodedRow = paste(parent_values, collapse="")
      
      
      node[binaryEncodedRow, col] = node[binaryEncodedRow, col] + 1
    }
  }
  return(node)
}
normalizeValues = function(nodeNetwork, total_cases){
  #for(node in 1:length(nodeNetwork)){
  #nodeNetwork[[node]] = nodeNetwork[[node]][]/10000 
  #}
  for(index in 1:length(nodeNetwork)){
    if(index != 2){
      rows = rownames(nodeNetwork[[index]])
      cols = colnames(nodeNetwork[[index]])
      for(row in rows){
        row_sum = nodeNetwork[[index]][row, cols[1]] + nodeNetwork[[index]][row, cols[2]]
        for(col in cols){
          nodeNetwork[[index]][row, col] = nodeNetwork[[index]][row, col] / row_sum
        }
      }
    }
  }
  
  return(nodeNetwork)
}
createNodeMatrix = function(parents){
  matrix_dimension = 2^(parents + 1)
  rows = 2^parents
  matrix_columns = 2 # Lockzed due to possible observations being 0 or 1
  prob_matrix = matrix(1:matrix_dimension, ncol = matrix_columns)
  if(parents > 0){
    binary_value_matrix = expand.grid(replicate(parents, 0:1, simplify = FALSE))
    bin_values = c()
    for (i in 1:nrow(binary_value_matrix)){
      binCombination = paste(binary_value_matrix[i,], collapse="")
      bin_values = c(bin_values, binCombination)
    }
  } else {
    bin_values = "1"
  }
  colnames(prob_matrix) = c("0", "1")
  rownames(prob_matrix) = bin_values
  prob_matrix[] = 0
  return(prob_matrix)
}
