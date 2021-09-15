
#################################################################
#
#################################################################
generate.param.membership = function(n, k){
	nk = n%/%k
	clu.sizes = rep(nk, k)
	nb.remaining = n%%k
	remaining = rep(0, k)
	if(nb.remaining > 0)
		remaining[1:nb.remaining] = 1
	clu.sizes = clu.sizes + remaining
	membership <- rep(1:k,clu.sizes)
	return(membership)
}



# we suppose that the graph is complete. Because, we have not applied yet density here.
compute.prop.neg = function(n, d, k, prop.mispl){

	prop.mispl = as.numeric(prop.mispl)
	n = as.numeric(n)
	k = as.numeric(k)
	
	# membership <- rep(1:k,each=n%/%k)
	membership <- generate.param.membership(n, k)
	pext <- sum(apply(t(combn(x=max(membership),m=2,simplify=TRUE)), 1, function(r)
					{	n1 <- length(which(membership==r[1]))
						n2 <- length(which(membership==r[2]))
						n1 * n2
					})) / (n*(n-1)/2)
	
			
	if(d == 1){
		prop.neg = pext
	} else {
		prop.neg = ((1-prop.mispl) - (1-pext))/(1-2*prop.mispl) # proportion of negative links
	}
	
	return(prop.neg)
}