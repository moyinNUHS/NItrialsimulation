getoutcome<-function(vector.outcome1, vector.outcome0, intervention){
  outcome<-c()
  for (i in 1:length(vector.outcome0)){
    if (intervention[i]==1) {outcome[[i]]=vector.outcome1[i]} else {outcome[[i]]=vector.outcome0[i]}
  }
  return(unlist(outcome))
}