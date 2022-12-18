library(testit)
library(tidyverse)
library(ggplot2)
### parameters
ppl_Arriv = 3000
lambda_extract = 20
lambda_detect = 60
lambda_pool = 300


# ------- pooled sample testing  -----------------------

pooled_sample_testing = function(T = 100, lab = 2, assistants = 2, n_samples = 25, positivity_const = 0.03){
  #initialization 
  t = 0
  # server times in lab_i = t_s[(assistant*(i-1) + 1):(assistant*i)]
  t_s = rep(Inf, lab*assistants)
  arrival_counter = 0
  test_counter = 0
  in_server_ID = vector(mode = 'list', length = lab*assistants)
  arrival_time_record = c()
  positivity_record = c()
  departure_time_record = c()
  arrival_queue = c()
  internal_queue = vector(mode = 'list', length = lab)
  for(i in 1:lab){
    internal_queue[[i]] = list()
  }
  t_arrival = t + rexp(1, rate = ppl_Arriv)
  positivity = runif(1) < positivity_const
  
  while(1){
    if(t_arrival < min(t_s) & t_arrival <= T){
      arrival_counter = arrival_counter + 1
      arrival_queue = c(arrival_queue,arrival_counter)
      ID = c(ID, arrival_counter)
      arrival_time_record = c(arrival_time_record, t_arrival)
      positivity_record = c(positivity_record, positivity)
      departure_time_record = c(departure_time_record, Inf)
      t = t_arrival
      # generate new arrival
      t_arrival = t + rexp(1, rate = ppl_Arriv)
      positivity = runif(1) < positivity_const
      
      if(max(t_s) == Inf){ # there is an server available
        lab_ind = (which.max(t_s)-1) %/% assistants + 1
        server_ind = (which.max(t_s)-1) %% assistants + 1
        # add to internal queue
        if (length(arrival_queue) >= n_samples){
          #create a pooled sample  
          pool = arrival_queue[1:n_samples]
          arrival_queue = arrival_queue[-(1:n_samples)]
          #internal_queue[[lab_ind]] = c(internal_queue[[lab_ind]], list(pool) )
          if (n_samples == 1){
            Z = 0
          } else {
            Z = rexp(1,rate = lambda_pool/n_samples)
          }
        } else {
          pool = arrival_queue[1:length(arrival_queue)]
          #internal_queue[[lab_ind]] = c(internal_queue[[lab_ind]], list(pool) )
          if (length(arrival_queue) == 1){ 
            Z = 0 
          } else {
            Z = rexp(1,rate = lambda_pool/length(arrival_queue))
          }
          arrival_queue = c()
        }
        in_server_ID[[which.max(t_s)]] = pool
        X = rexp(1,rate = lambda_extract)
        Y = rexp(1,rate = lambda_detect)
        test_counter = test_counter + 1
        t_s[which.max(t_s)] = t + X + Y + Z
      }
      
    } # end case 1
    else if (min(t_s) < t_arrival | (min(c(t_s,t_arrival)) > T) ){ 
      ind  = which.min(t_s)
      if(min(t_s) == Inf){
        return(list( 'wait_time' = departure_time_record - arrival_time_record,'test_counter' = test_counter ))
      }
      lab_ind = (which.min(t_s)-1) %/% assistants + 1
      # check results:
      if (sum(positivity_record[in_server_ID[[ind]]]) == 0 | length(in_server_ID[[ind]]) == 1){#negative or pool has 1 sample
        departure_time_record[in_server_ID[[ind]]] = t_s[ind]
        t = t_s[ind]
        # add the front of internal queue to server
        if (length(internal_queue[[lab_ind]]) > 0){
          in_server_ID[[ind]] = internal_queue[[lab_ind]][[1]]
          internal_queue[[lab_ind]] = internal_queue[[lab_ind]][-1]
          if (length(in_server_ID[[ind]]) == 1){ 
            Z = 0 
          } else {
            Z = rexp(1,rate = lambda_pool/length(in_server_ID[[ind]]))
          }
          X = rexp(1,rate = lambda_extract)
          Y = rexp(1,rate = lambda_detect)
          test_counter  = test_counter + 1
          t_s[ind] = t + X + Y + Z
        } else{ # internal queue empty, add from the input queue
          if (length(arrival_queue) >= n_samples){
            #create a pooled sample  
            pool = arrival_queue[1:n_samples]
            arrival_queue = arrival_queue[-(1:n_samples)]
            #internal_queue[[lab_ind]] = c(internal_queue[[lab_ind]], list(pool) )
            if (n_samples == 1){
              Z = 0
            } else {
              Z = rexp(1,rate = lambda_pool/n_samples)
            }
            X = rexp(1,rate = lambda_extract)
            Y = rexp(1,rate = lambda_detect)
            in_server_ID[[ind]] = pool
            test_counter= test_counter + 1
            t_s[ind] = t + X + Y + Z
          } else if (length(arrival_queue) > 0 & length(arrival_queue) < n_samples){
            pool = arrival_queue[1:length(arrival_queue)]
            #internal_queue[[lab_ind]] = c(internal_queue[[lab_ind]], list(pool) )
            if (length(arrival_queue) == 1){
              Z = 0
            } else {
              Z = rexp(1,rate = lambda_pool/length(arrival_queue))
            }
            arrival_queue = c()
            in_server_ID[[ind]] = pool
            X = rexp(1,rate = lambda_extract)
            Y = rexp(1,rate = lambda_detect)
            test_counter = test_counter + 1
            t_s[which.min((t_s))] = t + X + Y + Z
          } else { # input queue also empty
            t_s[which.min((t_s))] = Inf
          }
        }
      } else{ #split into 2 partitions
        l = length(in_server_ID[[ind]])
        pat1 = in_server_ID[[ind]][1:round(l/2)]
        pat2 = in_server_ID[[ind]][-(1:round(l/2))]
        #append to the front
        internal_queue[[lab_ind]] = c(list(pat2),internal_queue[[lab_ind]])
        # add pat1 to server
        in_server_ID[[ind]] = pat1
        X = rexp(1,rate = lambda_extract)
        Y = rexp(1,rate = lambda_detect)
        if (length(pat1) == 1){
          Z = 0
        } else {
          Z = rexp(1,rate = lambda_pool/length(pat1))
        }
        test_counter = test_counter + 1
        t_s[which.min((t_s))] = t + X + Y + Z
      }
    } # end case 2
  }
}


bisection_pooled_sample_testing = function(T = 100, lab = 2, assistants = 2, n_samples = 5, positivity_const = 0.03){
  #initialization 
  t = 0
  # server times in lab_i = t_s[(assistant*(i-1) + 1):(assistant*i)]
  t_s = rep(Inf, lab*assistants)
  arrival_counter = 0
  test_counter = 0
  in_server_ID = vector(mode = 'list', length = lab*assistants)
  arrival_time_record = c()
  positivity_record = c()
  departure_time_record = c()
  arrival_queue = c()
  internal_queue = vector(mode = 'list', length = lab)
  for(i in 1:lab){
    internal_queue[[i]] = list()
  }
  t_arrival = t + rexp(1, rate = ppl_Arriv)
  positivity = runif(1) < positivity_const
  
  while(1){
    if(t_arrival < min(t_s) & t_arrival <= T){
      arrival_counter = arrival_counter + 1
      arrival_queue = c(arrival_queue,arrival_counter)
      ID = c(ID, arrival_counter)
      arrival_time_record = c(arrival_time_record, t_arrival)
      positivity_record = c(positivity_record, positivity)
      departure_time_record = c(departure_time_record, Inf)
      t = t_arrival
      # generate new arrival
      t_arrival = t + rexp(1, rate = ppl_Arriv)
      positivity = runif(1) < positivity_const
      if(max(t_s) == Inf){ # there is an server available
        lab_ind = (which.max(t_s)-1) %/% assistants + 1
        server_ind = (which.max(t_s)-1) %% assistants + 1
        # add to server
        if (length(arrival_queue) >= n_samples**2){
          #create a pooled sample 
          pool = arrival_queue[1:n_samples**2]
          arrival_queue = arrival_queue[-(1:n_samples**2)]
          #internal_queue[[lab_ind]] = c(internal_queue[[lab_ind]], list(pool) )
          Z = rexp(1,rate = lambda_pool/(2*n_samples^2))
          
          in_server_ID[[which.max(t_s)]] = pool
          X = rexp(1,rate = lambda_extract/(2*n_samples) )
          Y = rexp(1,rate = lambda_detect/(2*n_samples))
          t_s[which.max((t_s))] = t + X + Y + Z
        } else { # conventional direct testing
          sample = arrival_queue[1]
          #internal_queue[[lab_ind]] = c(internal_queue[[lab_ind]], list(pool) )
          Z = 0
          in_server_ID[[which.max(t_s)]] = sample
          X = rexp(1,rate = lambda_extract)
          Y = rexp(1,rate = lambda_detect)
          t_s[which.max((t_s))] = t + X + Y + Z
        }
      }
      
    } # end case 1
    else if (min(t_s) < t_arrival | (min(c(t_s,t_arrival)) > T) ){ 
      ind  = which.min(t_s)
      lab_ind = (which.min(t_s)-1) %/% assistants + 1
      if(min(t_s) == Inf){
        #return(departure_time_record - arrival_time_record)
        return(list( 'wait_time' = departure_time_record - arrival_time_record,'test_counter' = test_counter ))
      }
      lab_ind = (which.min(t_s)-1) %/% assistants + 1
      # check results:
      if(length(in_server_ID[[ind]]) == n_samples**2){
        test_counter = test_counter + 2 * n_samples
        # create the matrix
        A = matrix(data = positivity_record[in_server_ID[[ind]]], nrow = n_samples , ncol = n_samples)
        # find the intersections
        rows = rep(0,n_samples)
        cols = rep(0,n_samples)
        for(i in 1:n_samples){
          rows[i] = sum(A[i,]) 
          cols[i] = sum(A[,i])
        }
        # list of potential positives:
        potential = c()
        if (sum(rows) > 1 & sum(cols) > 1 ){
          for(i in 1:n_samples){
            if (rows[i] != 0){ #positive
              for(j in 1:n_samples){
                if (cols[j] != 0){
                  t = (i-1)*n_samples + j
                  id = in_server_ID[[ind]][t]
                  assert(id > 0)
                  potential = c(potential, id)
                  # add to internal queue
                  internal_queue[[lab_ind]] = c(list(id), internal_queue[[lab_ind]])
                }#positive
              }
            }
          }
        }
        # dismiss the remaining samples
        dis = in_server_ID[[ind]][!(in_server_ID[[ind]] %in% potential)]
        departure_time_record[dis] = t_s[ind]
        #
      } else {
        test_counter = test_counter + 1
        # conventional direct testing
        assert(length(in_server_ID[[ind]]) == 1)
        index = in_server_ID[[ind]]
        departure_time_record[index] = t_s[ind]
      }
      #update time
      t = t_s[ind]
      # now this server is free update t_s[ind]
      if (length(internal_queue[[lab_ind]]) > 0){# pop from internal queue
        # print(internal_queue[[lab_ind]])
        in_server_ID[[ind]] = internal_queue[[lab_ind]][[1]]
        internal_queue[[lab_ind]] = internal_queue[[lab_ind]][-1]
        # condition on number of the samples
        if (length(in_server_ID[[ind]]) == n_samples**2){
          
          Z = rexp(1,rate = lambda_pool/(2*n_samples^2) )
          X = rexp(1,rate = lambda_extract/(2*n_samples))
          Y = rexp(1,rate = lambda_detect/(2*n_samples))
          t_s[ind] = t + X + Y + Z
        } else { # 1 sample 
          assert(length(in_server_ID[[ind]]) == 1)
          Z = 0
          X = rexp(1,rate = lambda_extract)
          Y = rexp(1,rate = lambda_detect)
          t_s[ind] = t + X + Y + Z
        }
      } else if (length(arrival_queue) > 0 ){ # pop from arrival queue
        if (length(arrival_queue) >= n_samples**2){
          #create a pooled sample  
          pool = arrival_queue[1:n_samples**2]
          arrival_queue = arrival_queue[-(1:n_samples**2)]
          #internal_queue[[lab_ind]] = c(internal_queue[[lab_ind]], list(pool) )
          Z = rexp(1,rate = lambda_pool/(2*n_samples^2))
          in_server_ID[[ind]] = pool
          X = rexp(1,rate = lambda_extract/(2*n_samples))
          Y = rexp(1,rate = lambda_detect/(2*n_samples))
          t_s[ind] = t + X + Y + Z
        } else {
          sample = arrival_queue[1]
          arrival_queue = arrival_queue[-1]
          #internal_queue[[lab_ind]] = c(internal_queue[[lab_ind]], list(pool) )
          Z = 0
          in_server_ID[[ind]] = sample
          X = rexp(1,rate = lambda_extract)
          Y = rexp(1,rate = lambda_detect)
          t_s[ind] = t + X + Y + Z
        }
      } else {# set to Inf
        t_s[ind] = Inf
      }
      
    } # end case 2
  }
}



p = seq(0,0.3,0.01)
avg_wait_time_conv = rep(0,length(p))
avg_wait_time_pool_25 = rep(0,length(p))
avg_wait_time_pool_10 = rep(0,length(p))
avg_wait_time_pool_5 = rep(0,length(p))
avg_wait_time_bisec55 = rep(0,length(p))

test_conv = rep(0,length(p))
test_pool_25 = rep(0,length(p))
test_pool_10 = rep(0,length(p))
test_pool_5 = rep(0,length(p))
test_bisec = rep(0,length(p))
for (i in seq_along(p)){
  ll = pooled_sample_testing(T = 40, lab= 2, assistants = 2,n_samples = 1)
  avg_wait_time_conv[i] = mean(ll$wait_time)
  test_conv[i] = ll$test_counter
  ll = pooled_sample_testing(T = 40, lab= 2, assistants = 2,n_samples = 25, positivity_const = p[i])
  avg_wait_time_pool_25[i] = mean(ll$wait_time)
  test_pool_25[i] = ll$test_counter
  ll = pooled_sample_testing(T = 40, lab= 2, assistants = 2,n_samples = 10, positivity_const = p[i])
  avg_wait_time_pool_10[i] = mean(ll$wait_time)
  test_pool_10[i] = ll$test_counter
  ll = pooled_sample_testing(T = 40, lab= 2, assistants = 2,n_samples = 5, positivity_const = p[i])
  avg_wait_time_pool_5[i] = mean(ll$wait_time)
  test_pool_5[i] = ll$test_counter
  
  ll = bisection_pooled_sample_testing(T = 40, lab= 2, assistants = 2,n_samples = 5, positivity_const = p[i])
  avg_wait_time_bisec55[i] = mean(ll$wait_time)
  test_bisec[i] = ll$test_counter
  #avg_wait_time_bisec[i] = mean(bisection_pooled_sample_testing(T = 40, lab= 3, assistants = 5,n_samples = 5,  positivity_const = p[i]))
  print(i)
}

df = tibble(avg_wait_time_conv, avg_wait_time_pool_25, avg_wait_time_pool_10, avg_wait_time_pool_5,'avg_wait_time_bisec5by5' = avg_wait_time_bisec55,'sample_positivity' = p)
df %>%
  pivot_longer(c(avg_wait_time_bisec5by5,avg_wait_time_pool_25,avg_wait_time_conv, avg_wait_time_pool_10, avg_wait_time_pool_5), names_to = "Diagnostic testing strategies", values_to = 'average wait time')%>%
  ggplot(aes(x = sample_positivity, y = `average wait time`, color = `Diagnostic testing strategies`)) + geom_point() + geom_line()
write.csv(df,'wait_time.csv')


df = tibble('conventional_direct' = test_conv,'pool25' = test_pool_25,'pool10' =test_pool_10,'pool5' =test_pool_5,'bidirectional5by5' = test_bisec,'sample_positivity' = p)
df %>%
  pivot_longer(c(conventional_direct,pool10,pool25,pool5,bidirectional5by5), names_to = "Diagnostic testing strategies", values_to = 'Amount of extract/detect kit used')%>%
  ggplot(aes(x = sample_positivity, y = `Amount of extract/detect kit used`, color = `Diagnostic testing strategies`)) + geom_point() + geom_line()
write.csv(df,'kit_consumption.csv')
