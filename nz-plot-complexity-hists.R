t1 = 1000
t2 = 5000
t3 = 10000
maxc = 45
res.factor = 2
probs = T

if(probs == T) {
  hist.freq = F
  trace.norm = T
} else {
  hist.freq = T
  trace.norm = F
}

for(expt in 1:3) {
  if(expt == 1) {
    png("nz-complexity-hists-directed.png", width=800*res.factor, height=800*res.factor, res=72*res.factor)
    df = read.csv("stats-out-blend-1-0.100-10-16-16-64-1-0-1-5000-1e+08.csv", header=T)
  } else if(expt == 2) {
    png("nz-complexity-hists-undirected.png", width=800*res.factor, height=800*res.factor, res=72*res.factor)
    df = read.csv("stats-out-blend-0-0.100-10-16-16-64-1-0-1-5000-1e+08.csv", header=T) 
  } else {
    png("nz-complexity-hists-sampled.png", width=700*res.factor, height=300*res.factor, res=72*res.factor)
    df = read.csv("stats-out-blend-1-0.100-10-16-16-64-0-0-1-0-1e+08.csv", header=T)
  }

if(expt == 3) {
    par(mfrow=c(1,2))
    hist(df$CMNNonzero, xlim=c(0,maxc), xlab="~                                          \nK (number of interface types)", ylab = "Frequency", main="All structures", breaks=0:45, freq=hist.freq)
    hist(df$CMNNonzero[df$Size == 16], xlim=c(0,maxc), xlab="~                                          \nK (number of interface types)", ylab = "Frequency", main="All 16-mers", breaks=0:45, freq=hist.freq)
    dev.off()
  }
  else {
    par(mfrow=c(4,4))
    
    hist(df$CMNNonzero[df$MinDiscoveryTime < t1], xlim=c(0,maxc), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="Without occurrence counts:\nAll structures\ndiscovered t<t1", breaks=0:45, freq=hist.freq)
    hist(df$CMNNonzero[df$MinDiscoveryTime < t2], xlim=c(0,maxc), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="All structures\ndiscovered t<t2", breaks=0:45, freq=hist.freq)
    hist(df$CMNNonzero[df$MinDiscoveryTime < t3], xlim=c(0,maxc), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="All structures\ndiscovered t<t3", breaks=0:45, freq=hist.freq)
    hist(df$CMNNonzero, xlim=c(0,maxc), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="All structures", breaks=0:45, freq=hist.freq)
    
    if(length(df$CMNNonzero[df$Size == 16 & df$MinDiscoveryTime < t1]) == 0) { plot(0,type='n',axes=FALSE,ann=FALSE) }
    else { hist(df$CMNNonzero[df$Size == 16 & df$MinDiscoveryTime < t1], xlim=c(0,maxc), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="16-mers\ndiscovered t<t1", breaks=0:45, freq=hist.freq) }
    hist(df$CMNNonzero[df$Size == 16 & df$MinDiscoveryTime < t2], xlim=c(0,maxc), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="16-mers\ndiscovered t<t2", breaks=0:45, freq=hist.freq)
    hist(df$CMNNonzero[df$Size == 16 & df$MinDiscoveryTime < t3], xlim=c(0,maxc), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="16-mers\ndiscovered t<t3", breaks=0:45, freq=hist.freq)
    hist(df$CMNNonzero[df$Size == 16], xlim=c(0,maxc), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="All 16-mers", breaks=0:45, freq=hist.freq)
    
    new.df = data.frame(Complexity=NULL, Time=NULL, Count=NULL)
    for(i in seq(from=0,to=maxc)) {
      for(j in 1:4) {
        new.df = rbind(new.df, data.frame(Complexity = i, Time = ifelse(j == 1, t1, ifelse(j == 2, t2, ifelse(j==3, t3, 20000))), Count=0))
      }
    }
    new.16.df = new.df
    
    for(i in 1:nrow(df)) {
      t = df$MinDiscoveryTime[i]
      if(t < t1) { ref = which(new.df$Complexity == df$CMNNonzero[i] & new.df$Time == t1); new.df$Count[ref] = new.df$Count[ref] + df$DiscoveryCount[i]  }
      if(t < t2) { ref = which(new.df$Complexity == df$CMNNonzero[i] & new.df$Time == t2); new.df$Count[ref] = new.df$Count[ref] + df$DiscoveryCount[i]  }
      if(t < t3) { ref = which(new.df$Complexity == df$CMNNonzero[i] & new.df$Time == t3); new.df$Count[ref] = new.df$Count[ref] + df$DiscoveryCount[i]  }
      ref = which(new.df$Complexity == df$CMNNonzero[i] & new.df$Time == 20000); new.df$Count[ref] = new.df$Count[ref] + df$DiscoveryCount[i]  
      if(df$Size[i] == 16) {
        if(t < t1) { ref = which(new.16.df$Complexity == df$CMNNonzero[i] & new.16.df$Time == t1); new.16.df$Count[ref] = new.16.df$Count[ref] + df$DiscoveryCount[i]  }
        if(t < t2) { ref = which(new.16.df$Complexity == df$CMNNonzero[i] & new.16.df$Time == t2); new.16.df$Count[ref] = new.16.df$Count[ref] + df$DiscoveryCount[i]  }
        if(t < t3) { ref = which(new.16.df$Complexity == df$CMNNonzero[i] & new.16.df$Time == t3); new.16.df$Count[ref] = new.16.df$Count[ref] + df$DiscoveryCount[i]  }
        ref = which(new.16.df$Complexity == df$CMNNonzero[i] & new.16.df$Time == 20000); new.16.df$Count[ref] = new.16.df$Count[ref] + df$DiscoveryCount[i]  
      }
    }
    
    plot(new.df$Complexity[new.df$Time == t1], new.df$Count[new.df$Time == t1]/ifelse(trace.norm, sum(new.df$Count[new.df$Time == t1]), 1), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="With occurrence counts:\nAll structures\ndiscovered t<t1", type="b")
    plot(new.df$Complexity[new.df$Time == t2], new.df$Count[new.df$Time == t2]/ifelse(trace.norm, sum(new.df$Count[new.df$Time == t2]), 1), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="All structures\ndiscovered t<t2", type="b")
    plot(new.df$Complexity[new.df$Time == t3], new.df$Count[new.df$Time == t3]/ifelse(trace.norm, sum(new.df$Count[new.df$Time == t3]), 1), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="All structures\ndiscovered t<t3", type="b")
    plot(new.df$Complexity[new.df$Time == 20000], new.df$Count[new.df$Time == 20000]/ifelse(trace.norm, sum(new.df$Count[new.df$Time == 20000]), 1), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="All structures", type="b")
    
    plot(new.16.df$Complexity[new.16.df$Time == t1], new.16.df$Count[new.16.df$Time == t1]/ifelse(trace.norm, sum(new.16.df$Count[new.16.df$Time == t1]), 1), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="16-mers\ndiscovered t<t1", type="b")
    plot(new.16.df$Complexity[new.16.df$Time == t2], new.16.df$Count[new.16.df$Time == t2]/ifelse(trace.norm, sum(new.16.df$Count[new.16.df$Time == t2]), 1), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="16-mers\ndiscovered t<t2", type="b")
    plot(new.16.df$Complexity[new.16.df$Time == t3], new.16.df$Count[new.16.df$Time == t3]/ifelse(trace.norm, sum(new.16.df$Count[new.16.df$Time == t3]), 1), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="16-mers\ndiscovered t<t3", type="b")
    plot(new.16.df$Complexity[new.16.df$Time == 20000], new.16.df$Count[new.16.df$Time == 20000]/ifelse(trace.norm, sum(new.16.df$Count[new.16.df$Time == 20000]), 1), xlab="~                                                \nK (number of interface types)", ylab="Frequency", main="All 16-mers", type="b")
    
    dev.off()
  }
}