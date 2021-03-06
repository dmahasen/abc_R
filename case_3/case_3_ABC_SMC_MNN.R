library(tmvtnorm)
 
source('case_3_preamble.R')
data <- read.csv(file.path("..", "data", "citrus_tristeza_data.csv"))
# Status of data classified as: 0-  susceptable, 1 - infected in 1981; 2 - infected in 1982
#  Indexes of trees infected in 1981
inf0<-data$id[data$status==1] 
#  Indexes of trees infected in 1982
inf1<-data$id[data$status==2] 

# Number of neighbours for covariance matrix calculations
M <- 50

# Number of particles
N <- 1000

# Number of simulations for each parameter set
n <- 1

# Thresholds
epsilon<- c(20, 15, 10, 7.5, 5, 4.5)

# Number of generations
G <- length(epsilon)

#  Lower and upper boundaries for priors
lm.low<-c(0.01, 0)
lm.upp<-c(5, 10)

# Empty matrices to store results (population plus 5 model parameters)
res.old<-matrix(ncol=2,nrow=N)
res.new<-matrix(ncol=2,nrow=N)

# Empty vectors to store weights
w.old<-rep(NA, N)
w.new<-rep(NA, N)

# Calculate the distances between trees
xy <- cbind(data$x, data$y)
d<- as.matrix(dist(xy, method = 'euclidean', diag = TRUE, upper = TRUE))
#  Find a set of all distances between trees
uniq.dist<-sort(unique(c(d)))
# Calculate a number of observed minimal distances
near.dist.obs<-sapply(1:length(inf1), function(a) min(as.numeric(d[inf1[a], inf0])))
n.obs<-sapply(1:length(uniq.dist), function(a) length(which(near.dist.obs==uniq.dist[a])))

for(g in 1:G){  

	#Initiate counter
	i<-1	
	while(i <= N){ # While the number of accepted particles is less than N_particles
   		if(g==1){
    			# Sample from prior distributions 
 			alpha<-runif(1,min=lm.low[1], max=lm.upp[1])
			beta<-runif(1, min=lm.low[2], max=lm.upp[2])
		} else {
			#  Select particle from previous generation
			p<-sample(seq(1,N),1,prob=w.old)		
			sigma<-Sigma[[p]]
			par<- rK(as.numeric(res.old[p,]),sigma, lm.low, lm.upp)
			alpha<-par[1]
			beta<-par[2]
		}
      	#  Test if prior non zero
      	if(prior.non.zero(c(alpha,beta))) {
    			# Set number of accepted simulations to zero
    			m<-0
    			distance <-rep(NA,n)
    			for(j in 1:n){
    				D_star<-run_model(alpha, beta)     
    				# Calculate distances 
    				calc.dist<-calc_distance(D_star)
    				distance[j] <-calc.dist    
    				if(calc.dist <= epsilon[g]){ # If distance is less than tolerance
    					m<-m+1
    				}
    			}	
    			if (m>0){
    				# Store results
    				res.new[i,]<-c(alpha, beta)  
      			# Calculate weights
      			w1<-prod(sapply(1:2, function(b) dunif(res.new[i,b], min=lm.low[b], max=lm.upp[b])))
				if(g==1){
					w2<-1
				} else {
					w2<-sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
				}
      			w.new[i] <- (m/n)*w1/w2
      			# Update counter
      			i <- i+1
      			print(paste0('Generation: ', g, ", particle: ", i))
      			}
    		} 
    	}
    	Sigma <- list(NA, N)
    for(p in 1:N){
    	Sigma[[p]]<- getSigmaNeighbours(M, res.new[p,], res.new) 
    }
    	res.old<-res.new
	w.old<-w.new/sum(w.new)

 	write.csv(res.new, file = paste("results_case_3_ABC_SMC_MNN_gen_",g,".csv",sep=""), row.names=FALSE)
}
