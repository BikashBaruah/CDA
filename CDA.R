getwd()
Dataset = read.csv("SRP300091.csv")

n_cmnty = 0
repeat
{
  n_cmnty = n_cmnty + 1
  ###### Calculate correlation matrix ##########
  library(Biobase)
  
  geneData = Dataset
  dim(geneData)
  
  rownames(geneData) = geneData$Gene.Symbol
  geneData = geneData[c(-1)]
  
  geneData = t(geneData)
  geneData = as.data.frame(geneData)
  
  # correlation of samples
  my_cor_matrix <- cor(geneData)
  my_cor_matrix = abs(my_cor_matrix)
  
  ### Convert Correlation matrix into Binary correlation Matrix ######
  bin_cor_matrix = my_cor_matrix
  n_r = nrow(my_cor_matrix)
  n_c = nrow(my_cor_matrix)
  
  for(i in 1:n_r)
  {
    for(j in 1:n_c)
    {
      if (is.na(my_cor_matrix[i,j]) == TRUE)
      {
        bin_cor_matrix[i,j] = 0
      }
      else if(my_cor_matrix[i,j] >= 0.5)
      {
        bin_cor_matrix[i,j] = 1
      }
      else
      {
        bin_cor_matrix[i,j] = 0
      }
    }
    print(paste("i = ", i))
  }
  
  ######### Sort the Degree #######
  Dataset_Degree = data.frame()
  nr = nrow(bin_cor_matrix)
  nc = ncol(bin_cor_matrix)
  for(i in 1:nr)
  {
    Dataset_Degree[i,1] = Dataset[i,1]
    temp = 0
    for(j in 1:nc)
    {
      if(bin_cor_matrix[i,j] == 1)
      {
        temp = temp + 1
      }
    }
    Dataset_Degree[i,2] = temp
    print(paste("Dataset Degree = ", i))
  }
  colnames(Dataset_Degree)[1] = "Gene Symbol"
  colnames(Dataset_Degree)[2] = "Degree"  
  df = Dataset_Degree[order(decreasing = TRUE, Dataset_Degree$Degree), ]  
  Ordered_Dataset_Degree = df
  
  ############## Define a community 
  CI = 0
  Community = data.frame()
  P1_P2_CI = data.frame()
  abc = 1
  P1_P2_CI[abc,1] = Ordered_Dataset_Degree[1,1]
  P1_P2_CI[abc,2] = 0
  P1_P2_CI[abc,3] = 0
  P1_P2_CI[abc,4] = 0
  colnames(P1_P2_CI)[1] = "Gene Symbol"
  colnames(P1_P2_CI)[2] = "P1"
  colnames(P1_P2_CI)[3] = "P2"
  colnames(P1_P2_CI)[4] = "CI"
  
  ###### Add the first node in the community ######
  Count = 1
  m = as.numeric(row.names(Ordered_Dataset_Degree)[1])
  Community[Count, c(1:ncol(Dataset))] = Dataset[m, c(1:ncol(Dataset))]
  
  ######### Add the next node in the community
  for(a in 1:3)
    #for(a in 1:nrow(Dataset))
  {
    if(a <= nrow(Community))
    {
      Temp = data.frame()
      Count = 0
      Ordered_Temp = data.frame()
      for(j in 1:nrow(bin_cor_matrix))
      {
        if(Community[a,1] == row.names(bin_cor_matrix)[j])
        {
          for(k in 1:ncol(bin_cor_matrix))
          {
            if(bin_cor_matrix[j,k] == 1)
            {
              node = colnames(bin_cor_matrix)[k]
              Entry = 1
              for(m in 1:nrow(Community))
              {
                if(node == Community[m,1])
                {
                  Entry = 0
                  break
                }
              }
              if(Entry == 1)
              {
                Count = Count + 1
                Temp[Count,1] = Dataset_Degree[k,1]
                Temp[Count,2] = Dataset_Degree[k,2]
                row.names(Temp)[Count] = row.names(Dataset_Degree)[k]
              }
            }
          }
          break
        }
      }
      if(nrow(Temp) > 0)
      {
        ####Sort the Temp data frame in non increasing order and store in Ordered_Temp.
        Ordered_Temp = Temp[order(decreasing = TRUE, Temp$V2), ]
        for (b in 1:nrow(Ordered_Temp))
        {
          ##### Calculate the CI Score.
          ######## Check the CI Score of the community with the new node ###
          next_node = data.frame()
          m = as.numeric(row.names(Ordered_Temp)[b])
          next_node[1, c(1:ncol(Dataset))] = Dataset[m, c(1:ncol(Dataset))]
          Community_Temp = Community
          Community_Temp[nrow(Community)+1, c(1:ncol(Dataset))] = Dataset[m, c(1:ncol(Dataset))]
          
          #### Edge count multiplied by 2 is the degree ######
          library(Biobase)
          geneData = Community_Temp
          dim(geneData)
          
          rownames(geneData) = geneData$Gene.Symbol
          geneData = geneData[c(-1)]
          
          geneData = t(geneData)
          geneData = as.data.frame(geneData)
          
          # correlation of samples
          my_cor_matrix <- cor(geneData)
          my_cor_matrix = abs(my_cor_matrix)
          #my_cor_matrix <- abs(my_cor_matrix)
          
          degree_count = 0
          n_r = nrow(my_cor_matrix)
          n_c = ncol(my_cor_matrix)
          
          for(i in 1:n_r)
          {
            for(j in 1:n_c)
            {
              if (is.na(my_cor_matrix[i,j]) == TRUE)
              {
                my_cor_matrix[i,j] = 0
              }
              if(my_cor_matrix[i,j] >= 0.5)
              {
                degree_count = degree_count + 1
                #        print(degree_count)
              }
            }
          }
          
          ##### Add the degrees of the nodes of the community w.r.t. the dataset #####
          n_r_c = nrow(Community_Temp)
          n_r_d = nrow(Dataset_Degree)
          add_degree = 0
          
          for(i in 1:n_r_c)
          {
            for(j in 1:n_r_d)
            {
              if(Community_Temp[i,1] == Dataset_Degree[j,1])
              {
                add_degree = add_degree + Dataset_Degree[j,2]
              }
            }
          }
          
          ##### Caculate First Part #######
          P_1 = degree_count / add_degree 
          
          
          ########## 2nd part started #######
          ############## Calculate BaryCentre #############
          Bary = data.frame()
          Bary[1,1] = "BaryCentre"
          
          for(j in 2:ncol(Community))
          {
            Bary[1,j] = mean(Community[ ,j]) 
          }
          
          ###### Correlation between BaryCentre and new gene i.i. the P_2 #####
          Bary_Temp = Bary[ ,-c(1)]
          next_node_Temp = next_node[ ,-c(1)]
          Bary_Temp = as.matrix(Bary_Temp)
          next_node_Temp = as.matrix(next_node_Temp)
          P_2 = cor(next_node_Temp[1, ], Bary_Temp[1, ])
          P_2 = abs(P_2)
          
          ########### Calculate the new Community Index (CI_New) #######
          CI_New = P_1 * P_2
          print(paste("Indentifying Community = ", n_cmnty))
          print(paste("a = ", a))
          print(paste("b = ", b))
          print(CI_New)
          print(paste("abc = ", abc))
          
          CI = CI - (CI * 0.06)
          if(CI_New >= CI)
          {
            ######Insert the new node in the community
            Community = Community_Temp
            CI = CI_New
            abc = abc + 1
            P1_P2_CI[abc,1] = Ordered_Temp[b,1]
            P1_P2_CI[abc,2] = P_1
            P1_P2_CI[abc,3] = P_2
            P1_P2_CI[abc,4] = CI
            print(paste("New node Inserted = ", abc))
          }
        }
      }
    }
    #  print(paste("a = ",a))
  }
  
  if(nrow(Community) <= 30)
  {
    break
  }
  write.csv(P1_P2_CI, paste("P1_P2_CI_DEG_SRP300091_N_06%_",n_cmnty,".csv"), row.names = FALSE)
  write.csv(Community, paste("Community_DEG_SRP300091_N_06%_",n_cmnty,".csv"), row.names = FALSE)
  
  #### Removes the nodes of identified communities from the dataset ##
  Dataset_new = data.frame()
  n_r_d = nrow(Dataset)
  n_r_c = nrow(Community)
  Count = 0
  
  for(i in 1:n_r_d)
  {
    Entry = 1
    for(j in 1:n_r_c)
    {
      if(Dataset[i,1] == Community[j,1])
      {
        Entry = 0
        break
      }
    }
    if(Entry == 1)
    {
      Count = Count + 1
      Dataset_new[Count,c(1:ncol(Dataset))] = Dataset[i,c(1:ncol(Dataset))] 
      row.names(Dataset_new)[Count] = Count
    }
    print(paste("Count = ", Count))
  }
  Dataset = Dataset_new
  write.csv(Dataset, paste("DEG_SRP300091_N_06%_for_",n_cmnty+1,"th_Community.csv"), row.names = FALSE)
}



