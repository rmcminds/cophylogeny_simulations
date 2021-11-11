n_sims <- 10
tree_size_range <- c(50,100)
rate_range <- c(0.1,10)

host <- lapply(sample(tree_size_range[[1]]:tree_size_range[[2]],n_sims),ape::rcoal)
symbiont <- lapply(sample(tree_size_range[[1]]:tree_size_range[[2]],n_sims),ape::rcoal)

host_trait <- lapply(host, function(h) ape::rTraitDisc(h, k=2, rate=exp(runif(1,log(rate_range[[1]]),log(rate_range[[2]])))))
symbiont_trait <- lapply(symbiont, function(h) ape::rTraitDisc(h, k=2, rate=exp(runif(1,log(rate_range[[1]]),log(rate_range[[2]])))))

data_matrices <- lapply(1:n_sims, function(x) {
  nhosts <- length(host[[x]]$tip.label)
  nsymbionts <- length(symbiont[[x]]$tip.label)
  mat <- matrix(0, nsymbionts, nhosts, dimnames=list(symbiont[[x]]$tip.label,host[[x]]$tip.label))
  for(h in 1:nhosts) {
    for(s in 1:nsymbionts) {
      if(host_trait[[x]][[h]] == 'B' & symbiont_trait[[x]][[s]] == 'B') {
        mat[s,h] <- 1
      }
    }
  }
  return(mat)
})


randfacts <- c('host', 'host.phy', 'symbiont', 'symbiont.phy', 'host.symbiont', 'host.phy.symbiont', 'host.symbiont.phy', 'host.symbiont.cophy')
rand <- as.formula(paste0('~ ',paste(randfacts, collapse=' + ')))
priorC <- list(B=list(mu=0, V=1e+8), R=list(V=1, nu=0))
priorC$G <-lapply(1:length(randfacts), function(x) {list(V=1, nu=1, alpha.mu=0, alpha.V=1000)})
names(priorC$G) <- paste("G", 1:length(randfacts), sep="")

hadfield_res <- lapply(1:n_sims, function(x) {
  
  hostA <- MCMCglmm::inverseA(host[[x]])$Ainv
  symbiontA <- MCMCglmm::inverseA(symbiont[[x]])$Ainv
  
  host.symbiontA<-as(kronecker(hostA, symbiontA), "dgCMatrix")                   # coevolutionary effect
  host.symbiontAS<-as(kronecker(hostA, Matrix::Diagonal(nrow(symbiontA))), "dgCMatrix")  # host evolutionary effect
  host.symbiontSA<-as(kronecker(Matrix::Diagonal(nrow(hostA)), symbiontA), "dgCMatrix")  # parasite evolutionary effect
  
  rownames(host.symbiontA)<-apply(expand.grid(rownames(symbiontA), rownames(hostA)), 1, function(x){paste(x[2],x[1], sep=".")})
  rownames(host.symbiontAS)<-rownames(host.symbiontSA)<-rownames(host.symbiontA)
  
  input <- reshape2::melt(data_matrices[[x]],as.is=T)
  colnames(input) <- c('symbiont','host','present')
  input$present <- factor(input$present, levels=c(0,1))
  input$host.phy <- input$host
  input$symbiont.phy <- input$symbiont
  input$host.symbiont <- paste(input$host,input$symbiont,sep='.')
  input$host.phy.symbiont <- input$host.symbiont
  input$host.symbiont.phy <- input$host.symbiont
  input$host.symbiont.cophy <- input$host.symbiont
  
  return(MCMCglmm::MCMCglmm(present ~ 1,
                            random = rand,
                            family="categorical",
                            data=input,
                            ginverse=list(host.phy=hostA, symbiont.phy=symbiontA, host.phy.symbiont=host.symbiontAS, host.symbiont.phy=host.symbiontSA, host.symbiont.cophy=host.symbiontA),
                            prior=priorC,
                            nitt=125000,
                            thin=5,
                            burnin=25000,
                            slice=T,
                            pr=T))
  
})

save.image('~/outputs/hadfield_simulations.RData')


n <- 1
heatmap(data_matrices[[n]], Rowv=as.dendrogram(as.hclust(ape::ladderize(symbiont[[n]]))), Colv=as.dendrogram(as.hclust(ape::ladderize(host[[n]],right=FALSE))), scale='none')
