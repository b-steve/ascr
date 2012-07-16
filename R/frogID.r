make.frog.captures=function(mics,clicks,dt){
#---------------------------------------------------------------------------------------
# Assigns new IDs to clicks on the basis of how close they are in time. If two consecutive
# clicks on mics i and j are separated by less time than it takes sound to travel from
# i to j they are assumed to have come from the same source click. Consecutive clicks on
# the same microphone are always assumed to be from different source clicks.
#
# DLB 14-06-12
#---------------------------------------------------------------------------------------
  K=dim(mics)[1]
  captures=clicks
  ct=rep(-Inf,K) # times of clicks in current set
  ID=1 # counter for click
  ct[clicks$trap[1]]=clicks$tim[1] # store time of click by mic 1 in ct[mic]
  new=FALSE # indicater that is true when click can't be part of current set
  nclicks=length(clicks$tim)
  for(i in 2:nclicks){
#    cat("doing record ",i,"\n")
    if(ct[clicks$trap[i]]>-Inf){ # next click on a mic already in current set so end this call
      nd=length(which(ct>-Inf)) # number clicks in current
      captures$ID[(i-nd):(i-1)]=ID # make all but last of those in current set part of same capture history
      ct=rep(-Inf,K) # re-initialise
      ct[clicks$trap[i]]=clicks$tim[i] # store time of click by mic i in ct[mic]
      ID=ID+1
      if(i==nclicks) captures$ID[i]=ID # write last record with new ID      
    }
    else { # next click on a mic not in current set
      ct[clicks$trap[i]]=clicks$tim[i] # store time of click by mic i in ct[mic]
      ctset=which(ct>-Inf)
      dts=dt[ctset,clicks$trap[i]] # times between mics in current set and mic i
      cts=-(ct[ctset]-clicks$tim[i]) # times between clicks in current set and click i
      if(any((cts-dts)>0)) new=TRUE
      if(new) { 
        nd=length(which(ct>-Inf))-1 # number clicks in current set before new click added
        captures$ID[(i-nd):(i-1)]=ID # make all but last of those in current set part of same capture history
        ct=rep(-Inf,K) # re-initialise 
        ct[clicks$trap[i]]=clicks$tim[i] # store time of click by mic i in ct[mic]
        ID=ID+1
        new=FALSE
        if(i==nclicks) captures$ID[i]=ID # write last record with new ID
      } else if(i==nclicks){
        nd=length(which(ct>-Inf)) # number clicks in current set
        captures$ID[(i-nd+1):i]=ID # make all in current set part of same capture history      
      }
    }
  }
  return(captures)
}
