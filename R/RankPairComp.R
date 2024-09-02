### pair.comp_{ij} = 1{Y_i < Y_j} ###

#' Compute Pairwise Comparison Matrix for Full Ranking List of the Entities
#'
#' Compute the pairwise comparison matrix from the ranking lists of the ranked entities.
#' @param rank.vec A full ranking list containing the ranks of all the \eqn{N} entities. Note that here we follow the usual definition of rank in R, that is, the larger the evaluation score of an entity, the larger this entity's rank is. Specifically, for a full ranking list of \eqn{N} entities, the rank of an entity equals \eqn{N+1} minus its ranked position.
#' @return An \eqn{N} by \eqn{N} pairwise comparison for all \eqn{N} entities, where the (\eqn{i},\eqn{j}) element equals 1 if \eqn{i} is ranked higher than \eqn{j}, and 0 if \eqn{i} is ranker lower than \eqn{j}. Note that the diagonal elements (\eqn{i},\eqn{i})'s are set to NA.
#' @export
FullRankToPairComp <- function( rank.vec, n = length(rank.vec) ){
  pair.comp <- matrix(NA, n, n)
  for(i in 1:n){
    j = which(rank.vec == i)
    pair.comp[j,  rank.vec > i] = 1
    pair.comp[j,  rank.vec < i] = 0
  }
  return(pair.comp)
}

# rank.vec = c(2,1,5,3,4)
# FullRankToPairComp(rank.vec)
# FullRankToPairComp(rank.vec) + t(FullRankToPairComp(rank.vec))



### pair.comp_{ij} = 1{Y_i < Y_j} ###
### rank matrix: each column denotes ranks of a subset of items ###
### overlap means whether there is overlapping items in these subsets ###

#' Compute Pairwise Comparison Matrix for Partial Ranking Lists of the Entities
#'
#' Compute the pairwise comparison matrix from partial ranking lists of the ranked entities.
#' @param rank.matrix An \eqn{N} by \eqn{K} matrix containing \eqn{K} partial ranking lists for the \eqn{N} ranked entities. Here each column of \code{rank.matrix} contains a partial ranking list for the \eqn{N} entities, where the value of an entity is NA if its pairwise comparison information is missing in this partial list. Note that here we follow the usual definition of rank in R, that is, the larger the evaluation score of an entity, the larger this entity's rank is. Specifically, for a full ranking list of \eqn{N} entities, the rank of an entity equals \eqn{N+1} minus its ranked position.
#' @return An \eqn{N} by \eqn{N} pairwise comparison for all \eqn{N} entities, where the (\eqn{i},\eqn{j}) element equals 1 if \eqn{i} is ranked higher than \eqn{j}, 0 if \eqn{i} is ranker lower than \eqn{j}, and NA if the relation between \eqn{i} and \eqn{j} is missing. Note that the diagonal elements (\eqn{i},\eqn{i})'s are set to NA.
#' @export
PartRankToPairComp <- function( rank.matrix, n = nrow(rank.matrix), overlap = TRUE){
  pair.comp <- matrix(NA, n, n)
  nR = ncol(rank.matrix)
  for(i in 1:n){
    set1 = NULL
    set0 = NULL
    for(j in 1:nR){
      if( !is.na( rank.matrix[i, j]) ){
        set1 = union(set1, which( rank.matrix[, j] > rank.matrix[i, j] ) )
        set0 = union(set0, which( rank.matrix[, j] < rank.matrix[i, j] ) )
      }
    }
    pair.comp[i, set1] = 1
    pair.comp[i, set0] = 0
  }

  if(overlap){
    for(i in 1:n){
      set1 = which( pair.comp[i, ] == 1 )
      set1.old = c()
      while( !setequal(set1.old, set1) ){
        set1.old = set1
        set1.add = c()
        for(j in set1){
          set1.add = union( set1.add, which( pair.comp[j, ] == 1 ) )
        }
        set1 = union(set1, set1.add)
      }

      set0 = which( pair.comp[i, ] == 0 )
      set0.old = c()
      while( !setequal(set0.old, set0) ){
        set0.old = set0
        set0.add = c()
        for(j in set0){
          set0.add = union( set0.add, which( pair.comp[j, ] == 0 ) )
        }
        set0 = union(set0, set0.add)
      }

      pair.comp[i, set1] = 1
      pair.comp[i, set0] = 0
    }
  }


  check = pair.comp + t(pair.comp)
  if( sum( check[ !is.na(check) ] != 1 ) > 0 ){
    print("Inconsistency in the ranking lists!")
  }
  return(pair.comp)
}

# rank.matrix = matrix(NA, nrow = 5, ncol = 3)
# rank.matrix[,1] = c(2, 1, NA, NA, NA)
# rank.matrix[,2] = c(NA, NA, 3, 1, 2)
# rank.matrix[,3] = c(NA, 2, 1, NA, NA)
# PartRankToPairComp(rank.matrix)
# PartRankToPairComp(rank.matrix) + t(PartRankToPairComp(rank.matrix))
#
# rank.matrix = matrix(NA, nrow = 5, ncol = 2)
# rank.matrix[,1] = c(2, 1, NA, NA, NA)
# rank.matrix[,2] = c(NA, NA, 3, 1, 2)
# PartRankToPairComp(rank.matrix, overlap = FALSE)

# rank.vec = c(2,1,5,3,4)
# FullRankToPairComp(rank.vec)
# PartRankToPairComp( as.matrix(rank.vec, ncol=1) )



### RankToplist
#' Compute Ranks, Top List and Evaluation Scores from each other
#'
#' @param Rank A full ranking list containing the ranks of all the \eqn{N} ranked entities. Note that here we follow the usual definition of rank in R, that is, the larger the evaluation score of an entity, the larger this entity's rank is. Specifically, for a full ranking list of \eqn{N} entities, the rank of an entity equals \eqn{N+1} minus its ranked position.
#' @param Toplist A ordered vector listing the indices of the \eqn{N} entities with decreasing ranks or evaluation scores.
#' @param mu A vector of the evaluation scores for all the \eqn{N} ranked entities.
#' @return A list containing the ranks, top list and evaluation scores for all the \eqn{N} entities.
#' @export
RankToplist <- function(Rank = NULL, Toplist = NULL, mu = NULL){
  if(!is.null(Rank)){
    Toplist = order( Rank, decreasing = TRUE )
  }
  if(!is.null(Toplist)){
    Rank = length(Toplist) + 1 - sort(Toplist, decreasing = FALSE, index.return = TRUE)$ix
  }
  if(!is.null(mu)){
    Rank = rank(mu)
    Toplist = order(mu, decreasing = TRUE)
  }
  return(list(Rank=Rank, Toplist=Toplist, mu=mu))
}


# mu = rnorm(20)
# Rank = RankToplist(mu = mu)$Rank
# Toplist = RankToplist(mu = mu)$Toplist
#
# identical( RankToplist(Rank = Rank)$Toplist, Toplist )
# identical( RankToplist(Toplist = Toplist)$Rank, Rank )


### kendall tau distance
#' Compute the Normalized Kendall Tau Distance from two Full Ranking Lists
#'
#' @param Rank1 A full ranking list containing the ranks of all the \eqn{N} entities. Note that here we follow the usual definition of rank in R, that is, the larger the evaluation score of an entity, the larger this entity's rank is. Specifically, for a full ranking list of \eqn{N} entities, the rank of an entity equals \eqn{N+1} minus its ranked position.
#' @param Rank2 Another full ranking list containing the ranks of all the \eqn{N} entities.
#' @return The normalized Kendall tau distance from the two full ranking lists.
#' @export
RankDist <- function(Rank1, Rank2, method = "kendall"){
  return( ( 1 - cor(Rank1, Rank2, method = method) )/2 )
}


### kendall tau distance
#' Compute the Normalized Kendall Tau Distance from two Full Ranking Lists
#'
#' @param Rank1 A full ranking list containing the ranks of all the \eqn{N} entities. Note that here we follow the usual definition of rank in R, that is, the larger the evaluation score of an entity, the larger this entity's rank is. Specifically, for a full ranking list of \eqn{N} entities, the rank of an entity equals \eqn{N+1} minus its ranked position.
#' @param Rank2 Another full ranking list containing the ranks of all the \eqn{N} entities.
#' @return The normalized Kendall tau distance from the two full ranking lists.
#' @export
Partialtauhat <- function(Rank1, Rank2, method = "kendall"){
  # itemsinboth <- intersect(Rank1,Rank2)
  # Rank1subvec <- Rank1[Rank1 %in% Rank2]
  # Rank2subvec <- Rank2[Rank2 %in% Rank1]
  if(length(Rank1) != length(Rank2)){
    stop("Rank bectors must be of the same length")
  }

  # itemsinboth <- intersect(which(!is.na(Rank1)), which(!is.na(Rank2)))
  #
  # Rank1subvec <- Rank1[itemsinboth]
  # Rank2subvec <- Rank2[itemsinboth]

  # print("cor(Rank1subvec, Rank2subvec, method = method) = ")
  # print(cor(Rank1subvec, Rank2subvec, method = method) )
  #
  # print("(1 - cor(Rank1subvec, Rank2subvec, method = method))/2 = ")
  # print((1 - cor(Rank1subvec, Rank2subvec, method = method))/2 )
  #
  # print("Rank1subvec = ")
  # print(Rank1subvec)
  #
  # print("Rank2subvec = ")
  # print(Rank2subvec)

  # itemsinboth <- intersect(Rank1,Rank2)
  # Rank1subvec <- Rank1[Rank1 %in% Rank2]
  # Rank2subvec <- Rank2[Rank2 %in% Rank1]

  ntilde <- length(intersect(which(!is.na(Rank1)), which(!is.na(Rank2))))
  # ntilde <- length(unique(c(which(!is.na(Rank1)), which(!is.na(Rank2)))))

  # print("ntilde = ")
  # print(ntilde)

  tempsum <- 0
  for(i in 1:length(Rank1)){
    if(is.na(Rank1[i]) | is.na(Rank2[i])){
      next
    }
    for(j in 1:length(Rank1)){
      if((i==j)| (is.na(Rank1[j]) | is.na(Rank2[j]) )   ){
        next
      }

      if(  ((Rank1[i] <= Rank1[j]) & ( Rank2[i] <= Rank2[j] ))  ){
        tempsum = tempsum + 1
      }
      if(  ((Rank1[i] > Rank1[j]) & ( Rank2[i] > Rank2[j] ))  ){
        tempsum = tempsum + 1
      }

      if(  ((Rank1[i] > Rank1[j]) & ( Rank2[i] <= Rank2[j] ))  ){
        tempsum = tempsum - 1
      }
      if(  ((Rank1[i] <= Rank1[j]) & ( Rank2[i] > Rank2[j] ))  ){
        tempsum = tempsum - 1
      }

    }
  }

  # print("tempsum = ")
  # print(tempsum)
  #
  # # print("tempsum/(ntilde*(ntilde-1)/2) = ")
  # # print(tempsum/(ntilde*(ntilde-1)/2))
  #
  # print("tempsum/(ntilde*(ntilde-1)) = ")
  # print(tempsum/(ntilde*(ntilde-1)))
  #
  #
  # print(" (1 - tempsum/(ntilde*(ntilde-1)))  /2 = ")
  # print((1 - tempsum/(ntilde*(ntilde-1)))  /2 )

  tempres <- (1 - tempsum/(ntilde*(ntilde-1)))  /2

  # This works when there are no ties

  # return( ( 1 - cor(Rank1subvec, Rank2subvec, method = method) )/2 )
  return( tempres )

  }





### kendall tau distance
#' Compute the Normalized Kendall Tau Distance from two Full Ranking Lists
#'
#' @param Rank1 A full ranking list containing the ranks of all the \eqn{N} entities. Note that here we follow the usual definition of rank in R, that is, the larger the evaluation score of an entity, the larger this entity's rank is. Specifically, for a full ranking list of \eqn{N} entities, the rank of an entity equals \eqn{N+1} minus its ranked position.
#' @param Rank2 Another full ranking list containing the ranks of all the \eqn{N} entities.
#' @return The normalized Kendall tau distance from the two full ranking lists.
#' @export
ModKTau <- function(Rank1, Rank2, method = "kendall"){
  if(length(Rank1) != length(Rank2)){
    stop("Rank vectors must be of the same length")
  }
  # This works when there are no ties

  if(any(is.na(c(Rank1,Rank2)))){
    stop("This function is not designed for Ranks with NA values. Please input top K ranks and enter ranks of all items outside top K as K+1.")
  }

  n.item <- length(Rank1)


  tempsum <- 0
  for(i in 1:length(Rank1)){
    # if(is.na(Rank1[i]) | is.na(Rank2[i])){
    #   next
    # }
    for(j in 1:length(Rank1)){
      # if((i==j)| (is.na(Rank1[j]) | is.na(Rank2[j]) )   ){
      #   next
      # }

      # if(  ((Rank1[i] <= Rank1[j]) & ( Rank2[i] <= Rank2[j] ))  ){
      #   tempsum = tempsum + 1
      # }
      # if(  ((Rank1[i] > Rank1[j]) & ( Rank2[i] > Rank2[j] ))  ){
      #   tempsum = tempsum + 1
      # }

      if(  ((Rank1[i] > Rank1[j]) & ( Rank2[i] < Rank2[j] ))  ){
        tempsum = tempsum + 1
      }
      if(  ((Rank1[i] < Rank1[j]) & ( Rank2[i] > Rank2[j] ))  ){
        tempsum = tempsum + 1
      }

    }
  }

  # print("modified distance = ")
  # print(tempsum)
  #
  # print("modified correlation = ")
  # print( 1 - 4*tempsum/(n.item *( n.item - 1 ))   )

  modcorr <- 1 - 4*tempsum/(n.item *( n.item - 1 ))

  # print("modified normalized distance = ")
  # print( (1 - modcorr)/2)
  # print("Alt calc, modified normalized distance = ")
  # print( 2*tempsum/(n.item *( n.item - 1 )  ) )

  # return( ( 1 - cor(Rank1subvec, Rank2subvec, method = method) )/2 )
  return( (1 - modcorr)/2 )

}


