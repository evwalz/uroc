# This function computes the uROC curve as an weighted average over the ROC curves from the associated rocm

Exact = function(response, predictor){

  n = length(response)
  thresholds <- unique(response)[-1]
  N = length(thresholds)+1

  ncontrols <- (which(duplicated(response)==FALSE)-1)[-1]
  ncases <- (n-ncontrols)
  weights_all <- (ncases*ncontrols)
  weights = sum(weights_all)

  pre.order <- order(predictor, decreasing=TRUE)
  predictor <- predictor[pre.order]
  response <- response[pre.order]
  dups <- rev(duplicated(rev(predictor)))
  
  response.binary = lapply(thresholds, function(x, response){as.numeric(response>=x)}, response)
  Hitrate.all = lapply(response.binary, function(x, dups){c(0,(cumsum(x==1)*sum(x==0))[!dups])}, dups)
  Farate.all = lapply(response.binary, function(x, dups){c(0,(cumsum(x==0)/sum(x==0))[!dups])}, dups)
  farate.unique = sort(unique(c(unlist(Farate.all),0)))

  #Initialization
  final.hit.min.new = rep(0,length(farate.unique))
  final.hit.max.new = rep(0,length(farate.unique))
  hit.min.new = rep(0,length(farate.unique))
  hit.max.new = rep(0,length(farate.unique))
  
  for(i in 1:(N-1)){
    hitrate = Hitrate.all[[i]]
    farate = Farate.all[[i]]

    far.group.min <- which(as.numeric(duplicated(farate))==0)
    far.group.max = !duplicated(farate, fromLast = TRUE)

    #only select far-values with steps in graph and corresponding min and max value of hitrate
    far.max = farate[far.group.min]
    hit.min = hitrate[far.group.min]
    hit.max = hitrate[far.group.max]

    indx.ties = which(hit.min[-1] != hit.max[-length(hit.max)])

    cuts.min = c(-Inf,far.max)
    cuts.max = c(far.max,Inf)

    intervals.min = cut(farate.unique, breaks=cuts.min, labels=FALSE,right=TRUE)
    intervals.max = cut(farate.unique, breaks=cuts.max, labels=FALSE,right=FALSE)

    hit.max.new = hit.max[intervals.max]
    hit.min.new = hit.min[intervals.min]


    indx.ties.length = length(indx.ties)

    if(indx.ties.length>0){

      far.min.tied = far.max[indx.ties]
      far.max.tied = far.max[(indx.ties+1)]

      hit.max.tied = hit.min[(indx.ties+1)]
      hit.min.tied = hit.max[indx.ties]

      slope = (hit.max.tied -hit.min.tied)/(far.max.tied - far.min.tied)

      for(j in 1:indx.ties.length){
        indx = which(farate.unique > far.min.tied[j] & farate.unique < far.max.tied[j])
        hit.max.new[indx] = hit.min.new[indx] = slope[j]*(farate.unique[indx]-far.min.tied[j])+hit.min.tied[j]
      }
    }

    final.hit.min.new <- final.hit.min.new + hit.min.new
    final.hit.max.new <- final.hit.max.new + hit.max.new

  }

  return(list(Farate = sort(c(farate.unique,farate.unique)), Hitrate = sort(c(final.hit.min.new/weights,final.hit.max.new/weights))))
}


# This function computes an approximation of the uROC curve and takes into account ties in the predictor

Approx1 = function(response, predictor,space.size){

  n = length(response)
  N = length(unique(response))

  ncontrols <- (which(duplicated(response)==FALSE)-1)[-1]
  ncases <- (n-ncontrols)
  weights_all <- (ncases*ncontrols)
  weights = sum(weights_all)

  pre.order <- order(predictor, decreasing=TRUE)
  predictor <- predictor[pre.order]
  dups <- rev(duplicated(rev(predictor)))

  if(is.null(space.size)){
    space.size = 100
  }

  InterPoint = seq(0,1,1/space.size)
  Hit.weighted = rep(0,length(InterPoint))

  for(i in 1:(N-1)){
    controls = ncontrols[i]
    cases = ncases[i]
    response.binary = c(rep(0,controls),rep(1,cases))[pre.order]
    hitrate = c(0,cumsum(response.binary==1)[!dups]*controls)
    farate = c(0,cumsum(response.binary==0)[!dups]/controls)

    far.group.min = !duplicated(farate)
    far.group.max = !duplicated(farate, fromLast = TRUE)
    far.max = farate[far.group.max]

    hit.min = hitrate[far.group.min]
    hit.max = hitrate[far.group.max]

    # identify ties
    indx.ties = which(hit.min[-1] != hit.max[-length(hit.max)])

    if(length(indx.ties)>0){
      start.f = far.max[indx.ties]
      end.f = far.max[indx.ties+1]

      start.h = hit.max[indx.ties]
      end.h = hit.max[indx.ties+1]

      slope = (end.h-start.h)/(end.f-start.f)

      for(j in 1:length(start.f)){
        indx = which(InterPoint>start.f[j] & InterPoint<end.f[j])
        far.max = c(far.max, InterPoint[indx])
        hit.new = slope[j]*(InterPoint[indx]-start.f[j])+start.h[j]
        hit.max = c(hit.max,hit.new)
      }
      far.max = sort(far.max)
      hit.max = sort(hit.max)
    }
    intervals.max = cut(InterPoint, breaks=c(far.max,Inf), labels=FALSE,right=FALSE)
    Hit.weighted <- Hit.weighted + hit.max[intervals.max]
  }
  return(list(Farate = c(0,InterPoint), Hitrate = sort(c(0,Hit.weighted/weights))))
}


# This function computes an approximation of the uROC curve by ignoring ties in the predictor

Approx2 = function(response, predictor, space.size){

  n = length(response)
  N = length(unique(response))

  ncontrols <- (which(duplicated(response)==FALSE)-1)[-1]
  ncases <- (n-ncontrols)
  weights_all <- (ncases*ncontrols)
  weights = sum(weights_all)

  pre.order <- order(predictor, decreasing=TRUE)
  predictor <- predictor[pre.order]
  dups <- rev(duplicated(rev(predictor)))

  if(is.null(space.size)){
    space.size = 500
  }

  InterPoint = seq(0,1,1/space.size)
  Hit.weighted = rep(0,length(InterPoint))

  for (i in 1:(N-1)) {
    controls = ncontrols[i]
    cases = ncases[i]
    resBin = c(rep(0,controls),rep(1,cases))[pre.order]
    hitrate = c(0,cumsum(resBin==1)[!dups]*controls)
    farate = c(0,cumsum(resBin==0)[!dups]/controls)

    roc_curve = stepfun(x = farate[-1], y = hitrate)  
    Hit.weighted <- Hit.weighted + roc_curve(InterPoint)
  }

  return(list(Farate = c(0, InterPoint), Hitrate = sort(c(0, Hit.weighted / weights))))
}

