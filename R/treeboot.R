#.make_faux <- function(file) {
#  n <- 1000
#
#  set.seed(1)
#  Sigma <- matrix(c(1,0.8,0.2,0.8,1,-0.4,0.2,-0.4,1),3,3)
#  traits <- mvrnorm(n, rep(0,3), Sigma)
#  traits[,1] <- ifelse(traits[,1] < qnorm(0.60), 0, 1)
#  traits[,2] <- ifelse(traits[,2] < qnorm(0.95), 0, 1)
#  traits[,3] <- ifelse(traits[,3] < qnorm(0.75), 0, 1)
#
#  samp <- list()
#  graph <- erdos.renyi.game(n, 5*n, "gnm")
#  samp$nodes <- 1:n
#  samp$edges <- data.frame(node1=as_edgelist(graph)[,1], node2=as_edgelist(graph)[,2])
#
#  traits <- cbind(1:n, traits)
#  colnames(traits) <- c("id","X","Y","Z")
#  adj.mat <- matrix(0, n, n)
#  adj.mat[cbind(samp$edges$node1, samp$edges$node2)] <- 1
#  adj.mat[cbind(samp$edges$node2, samp$edges$node1)] <- 1
#
#  faux.network <- list(traits=data.frame(traits), adj.mat=adj.mat)
#
#  save(faux.network, file=file)
#}

.prop.vh <- function(samp, trait, pi) sum(trait[samp]/pi[samp], na.rm=T)/sum((!is.na(trait[samp]))/pi[samp])
.weight.vh <- function(samp, trait, pi) sum((!is.na(trait[samp]))/pi[samp])

.TBS <- function(samp, B) {
  n <- length(samp$nodes)

  samp.adj.mat <- matrix(F, n, n)
  samp.adj.mat[cbind(samp$edges$node1, samp$edges$node2)] <- T
  samp.adj.list <- list()
  for(i in 1:n) samp.adj.list[[i]] <- which(samp.adj.mat[i,])

  seeds <- which(apply(samp.adj.mat, 2, sum) == 0)
  num.seeds <- length(seeds)
  leaves <- apply(samp.adj.mat, 1, sum) == 0

  resamp <- list()
  for(b in 1:B) {
    nodes <- rep(NA, 2*n)
    nodes[1:num.seeds] <- seeds[sample.int(num.seeds, replace=T)]

    curr <- 1
    end <- length(seeds) + 1
    while(curr < end) {
      adj <- samp.adj.list[[nodes[curr]]]
      num.adj <- length(adj)
      if(end+num.adj-1 > length(nodes)) nodes <- c(nodes, rep(NA, n))
      if(num.adj > 0) nodes[end:(end+num.adj-1)] <- adj[sample.int(num.adj, replace=T)]
      curr <- curr + 1
      end <- end + num.adj
    }
    resamp[[b]] <- nodes[1:(end-1)]
  }

  return(resamp)
}

sample.RDS <- function(traits, adj.mat, n=100, num.seeds=1, num.samp=2, num.prob=c(rep(0,num.samp),1), replace=TRUE) {
  N <- dim(adj.mat)[1]
  pi <- apply(adj.mat, 1, sum)

  nodes <- rep(NA, n)
  edges <- data.frame(node1=rep(NA, n-num.seeds), node2=rep(NA, n-num.seeds))

  seeds <- sample(N, num.seeds, prob=pi, replace=replace)
  nodes[1:num.seeds] <- seeds
  if(!replace) adj.mat[,seeds] <- 0

  curr <- 1
  end <- num.seeds + 1
  while(curr < end) {
    adj <- which(adj.mat[nodes[curr],] == 1)
    if(replace) num <- min(sample.int(num.samp+1, 1, prob=num.prob)-1, n-end+1)
    else num <- min(sample.int(num.samp+1, 1, prob=num.prob)-1, length(adj), n-end+1)

    if(num > 0) {
      samp.adj <- adj[sample(length(adj), num, replace=replace)]
      nodes[end:(end+num-1)] <- samp.adj
      edges[(end-num.seeds):(end-num.seeds+num-1),] <- cbind(curr, end:(end+num-1))
      if(!replace) {
        adj.mat[nodes[curr],] <- 0
        adj.mat[,samp.adj] <- 0
      }
    }

    curr <- curr + 1
    end <- end + num
  }

  if(end <= n) {
    nodes <- nodes[-(end:n)]
    edges <- edges[-(end:n),]
  }

  return(list(nodes=traits[nodes,1], edges=edges, degree=pi[nodes], traits=traits[nodes,-1]))
}

treeboot.RDS <- function(samp, quant=c(0.025, 0.10, 0.90, 0.975), B=2000) {
  resamp <- .TBS(samp, B)

  results <- matrix(NA, dim(samp$traits)[2], length(quant))
  for(t in 1:dim(samp$traits)[2]) {
    p.TBS <- sapply(resamp, .prop.vh, samp$traits[,t], samp$degree)
    w.TBS <- sapply(resamp, .weight.vh, samp$traits[,t], samp$degree)
    for(q in 1:length(quant))
      results[t,q] <- sort(p.TBS)[max(which(cumsum(w.TBS[order(p.TBS)]/sum(w.TBS)) <= quant[q]))]
  }

  rownames(results) <- colnames(samp$traits)
  colnames(results) <- quant
  return(results)
}

