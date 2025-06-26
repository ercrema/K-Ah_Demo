dDoubleFlat  <- nimbleFunction(run = function(x = integer(0), a = double(0), b = double(0), mu = double(0), eta = double(0), log = integer(0)) 
			       {
				       returnType(double(0))
				       mu = round(mu)
				       t1 = 1:(abs(mu-a))
				       t2 = 1:(abs(b-mu)+1)
				       t1final = abs(mu-a)
				       t2final = abs(b-mu)+1
				       n1 = numeric(abs(mu-a))
				       n2 = numeric(abs(b-mu)+1)
				       eta1 = 1/(1+exp(eta))
				       eta2 = 1-eta1
				       for (i in 1:t1final)
				       {
					       n1[i] = eta1
				       }
				       for (i in 1:t2final)
				       {
					       n2[i] = eta2
				       }
				       n = c(n1,n2)
				       p = n/sum(n)
				       # This last bit would be the same for any model
				       logProb = dcat(a-x+1,prob=p,log=TRUE)
				       if(log) {
					       return(logProb)
				       } else {
					       return(exp(logProb))
				       }
			       })   


rDoubleFlat  <- nimbleFunction(run = function(n = integer(0), a = double(0), b = double(0), mu = double(0), eta = double(0)) 
			       {
				       returnType(double(0))
				       mu = round(mu)
				       t1 = 1:(abs(mu-a))
				       t2 = 1:(abs(b-mu)+1)
				       t1final = abs(mu-a)
				       t2final = abs(b-mu)+1
				       pop1 = numeric(abs(mu-a))
				       pop2 = numeric(abs(b-mu)+1)
				       eta1 = 1/(1+exp(eta))
				       eta2 = 1 - eta1
				       for (i in 1:t1final)
				       {
					       pop1[i] = eta1
				       }
				       for (i in 1:t2final)
				       {
					       pop2[i] = eta2
				       }
				       pop = c(pop1,pop2)
				       p = pop/sum(pop)
				       res=a-rcat(n=1,prob=p)+1
				       return(res)
			       })   

registerDistributions(list(
      dDoubleFlat = list(
        BUGSdist = "dDoubleFlat(a,b,mu,eta)",
        Rdist = "dDoubleFlat(a,b,mu,eta)",
        pqAvail = FALSE
      )))
