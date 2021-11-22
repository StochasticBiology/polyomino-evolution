# sprintf(outstr, "out-%i-%.3f-%i-%i-%i-%i-%i-%i-%i-%li.csv", DIRECTED, MUT, NPAR, TARGETSIZE, NTILE, NCOL, OUTPUTALL, CONFOUND, RSEED, NUMR)

library(ggplot2)
library(gridExtra)
library(GGally)

### complexity histograms for many different cases

experiments = c("stats-out-blend-1-0.100-10-16-16-64-0-0-1-5000-1e+08.csv",
                "stats-out-blend-1-0.100-10-16-16-64-1-0-1-5000-1e+08.csv",
                "stats-out-blend-1-0.100-10-8-16-64-0-0-1-5000-1e+08.csv",
		"stats-out-blend-1-0.100-10-15-16-64-0-0-1-5000-1e+08.csv",
		"stats-out-blend-1-0.010-10-16-16-64-0-0-1-5000-1e+08.csv",
		"stats-out-blend-1-1.000-10-16-16-64-0-0-1-5000-1e+08.csv",
		"stats-out-blend-0-0.100-10-16-16-64-1-0-1-5000-1e+08.csv",
		"stats-out-blend-1-0.100-10-16-16-64-0-1-1-5000-1e+08.csv",
                "stats-out-blend-0-0.100-10-16-16-64-1-1-1-5000-1e+08.csv",
		"stats-out-blend-1-0.100-10-16-16-64-0-2-1-5000-1e+08.csv",
                "stats-out-blend-0-0.100-10-16-16-64-1-2-1-5000-1e+08.csv",
		"stats-out-blend-1-0.100-10-16-16-64-0-3-1-5000-1e+08.csv",
                "stats-out-blend-0-0.100-10-16-16-64-1-3-1-5000-1e+08.csv")

titles = c("Default: F=16, targets only, mu = 0.1",
           "all structures",
           "F=8",
	   "F=15",
	   "mu = 0.01",
	   "mu = 1.0",
	   "undirected, all structures",
	   "confound 1",
	   "undirected, confound 1, all structures",
	   "confound 2",
	   "undirected, confound 2, all structures",
	   "confound 3",
	   "undirected, confound 3, all structures")

plots = list()
for(i in 1:length(experiments)) {
  df = read.csv(experiments[i], header=T)
  plots[[length(plots)+1]] = ggplot(df, aes(x=CMNNonzero, y=log(DiscoveryCount))) +
      geom_jitter(width=0.1, height = 0.01) +
      labs(title = titles[i]) +
      theme(legend.position = "none") 
      
}

png("nz-plotblend-0.png", width=800, height=800)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]], plots[[11]], plots[[12]], plots[[13]], nrow=4)
dev.off()

### complexity histograms for many different cases

experiments = c("stats-out-blend-1-0.100-10-16-16-64-0-0-1-5000-1e+08.csv",
		"stats-out-blend-0-0.100-10-16-16-64-1-0-1-5000-1e+08.csv",
		"stats-out-blend-1-0.100-10-16-16-64-0-1-1-5000-1e+08.csv",
                "stats-out-blend-0-0.100-10-16-16-64-1-1-1-5000-1e+08.csv",
		"stats-out-blend-1-0.100-10-16-16-64-0-2-1-5000-1e+08.csv",
                "stats-out-blend-0-0.100-10-16-16-64-1-2-1-5000-1e+08.csv",
		"stats-out-blend-1-0.100-10-16-16-64-0-3-1-5000-1e+08.csv",
                "stats-out-blend-0-0.100-10-16-16-64-1-3-1-5000-1e+08.csv")

titles = c("Default: F=16, targets only, mu = 0.1",
	   "undirected, all structures",
	   "confound 1",
	   "undirected, confound 1, all structures",
	   "confound 2",
	   "undirected, confound 2, all structures",
	   "confound 3",
	   "undirected, confound 3, all structures")

plots = list()
for(i in 1:length(experiments)) {
  df = read.csv(experiments[i], header=T)
  plots[[length(plots)+1]] = ggplot(df, aes(x=CMNNonzero, y=log(DiscoveryCount))) +
      geom_jitter(width=0.1, height = 0.01) +
      labs(title = titles[i]) +
      theme(legend.position = "none") 
      
}

png("nz-plotblend-0b.png", width=800, height=800)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], nrow=4)
dev.off()

### sampling effects

experiments = c("stats-out-blend-1-0.100-10-16-16-64-0-0-1-0-1e+07.csv",
		"stats-out-blend-1-0.100-10-16-16-64-0-0-1-0-1e+08.csv",
		"stats-out-blend-1-0.100-10-16-16-64-0-0-1-0-1e+09.csv")
plot.df = data.frame(Complexity = NULL, Frequency = NULL, Samples = NULL)
sample.labels = c("1e7", "1e8", "1e9")
for(i in 1:length(experiments)) {
  df = read.csv(experiments[i], header=T)
  plot.df = rbind(plot.df, data.frame(Complexity = df$CMNNonzero, Frequency = df$SampleCount, Samples = sample.labels[i]))
}

png("nz-plotblend-1.png", width=800, height=800)
ggplot(plot.df, aes(x=Complexity, y=log(Frequency), colour=Samples)) + geom_point()
dev.off()

### symmetry

experiments = c("stats-out-blend-1-0.100-10-16-16-64-0-0-1-5000-1e+08.csv",
                "stats-out-blend-0-0.100-10-16-16-64-1-0-1-5000-1e+08.csv",
		"stats-out-blend-1-0.100-10-16-16-64-0-0-1-5000-1e+08.csv")

symmgroups = c("D4", "C4", "D2", "C2", "D1", "C1")

df = read.csv(experiments[1])
df$SymmGroup = symmgroups[df$Symmetry+1]
symm.df = data.frame(SymmGroup = rep(symmgroups, 2), Class=c(rep("Evolved", length(symmgroups)), rep("Structures", length(symmgroups))), Count=0)
for(i in 1:nrow(df)) {
  ref = which(symm.df$SymmGroup == df$SymmGroup[i] & symm.df$Class == "Evolved")
  symm.df$Count[ref] = symm.df$Count[ref] + df$AdaptCount[i]/sum(df$AdaptCount)
  ref = which(symm.df$SymmGroup == df$SymmGroup[i] & symm.df$Class == "Structures")
  symm.df$Count[ref] = symm.df$Count[ref] + 1/nrow(df)
}

symm2.df = data.frame(SymmGroup = rep(symmgroups, 2), Class=c(rep("Sampled 16-mers", length(symmgroups)), rep("Sampled under 16", length(symmgroups)), rep("Directed 16-mers", length(symmgroups)), rep("Undirected all", length(symmgroups))), Count=0)
df = read.csv(experiments[2])
df$SymmGroup = symmgroups[df$Symmetry+1]
for(i in 1:nrow(df)) {
  if(df$Size[i] == 16) {
    ref = which(symm2.df$SymmGroup == df$SymmGroup[i] & symm2.df$Class == "Sampled 16-mers")
    symm2.df$Count[ref] = symm2.df$Count[ref] + 1/nrow(df)
  }
  if(df$Size[i] < 16) {
    ref = which(symm2.df$SymmGroup == df$SymmGroup[i] & symm2.df$Class == "Sampled under 16")
    symm2.df$Count[ref] = symm2.df$Count[ref] + 1/nrow(df)
  }
  ref = which(symm2.df$SymmGroup == df$SymmGroup[i] & symm2.df$Class == "Undirected all")
  symm2.df$Count[ref] = symm2.df$Count[ref] + df$DiscoveryCount[i]/sum(df$DiscoveryCount)
}
df = read.csv(experiments[3])
df$SymmGroup = symmgroups[df$Symmetry+1]
for(i in 1:nrow(df)) {
    ref = which(symm2.df$SymmGroup == df$SymmGroup[i] & symm2.df$Class == "Directed 16-mers")
    symm2.df$Count[ref] = symm2.df$Count[ref] + df$DiscoveryCount[i]/sum(df$DiscoveryCount)
}

symm.baseline = floor(min(log(symm.df$Count)))
symm.df$transLogCount = log(symm.df$Count)-symm.baseline
symm.plot.1 = ggplot(symm.df, aes(x=factor(SymmGroup, levels=symmgroups), y=transLogCount, fill=Class, colour=Class)) +
  geom_col(position="dodge") +
  scale_y_continuous(labels = function(y) y + symm.baseline)

symm2.baseline = floor(min(log(symm2.df$Count[symm2.df$Count != 0])))
symm2.df$transLogCount = log(symm2.df$Count)-symm2.baseline
symm.plot.2 = ggplot(symm2.df, aes(x=factor(SymmGroup, levels=symmgroups), y=transLogCount, fill=factor(Class, levels = c("Sampled 16-mers", "Sampled under 16", "Directed 16-mers", "Undirected all")))) +
  geom_col(position="dodge") +
  scale_y_continuous(labels = function(y) y + symm2.baseline) +
  labs(fill = "Class")


png("nz-plotblend-2.png", width=800, height=800)
grid.arrange(symm.plot.1, symm.plot.2, nrow=2)
dev.off()

### complexity dynamics

experiments = c("stats-out-blend-1-0.100-10-16-16-64-1-0-1-5000-1e+08.csv",
                "stats-out-blend-0-0.100-10-16-16-64-1-0-1-5000-1e+08.csv")
		
plots = list()
for(expt in 1:2) {
  df = read.csv(experiments[expt], header=T)
  pop.df = read.csv(gsub("stats", "pop", experiments[expt]), header=F)
  pop.df = pop.df[pop.df$V1 < 100,]
  complexity.df = pop.df
  for(i in 3:ncol(complexity.df)) {
    complexity.df[,i] = unlist(lapply(pop.df[,i], function(x) ifelse(x == -1, 0, df$CMNNonzero[x+1])))
  }
  complexity.stats.df = data.frame(Run=complexity.df$V1, Time=complexity.df$V2, Min=0, Mean=0, Max=0, HistMax=0)
  tmp.hist.max = 0
  for(i in 1:nrow(complexity.df)) {
    if(i > 1) {
      if(complexity.stats.df$Run[i] != complexity.stats.df$Run[i-1]) {
        tmp.hist.max = 0
      }
    }
    complexity.stats.df$Min[i] = min(as.numeric(complexity.df[i,3:ncol(complexity.df)]))
    tmp.max = max(as.numeric(complexity.df[i,3:ncol(complexity.df)]))
    complexity.stats.df$Max[i] = tmp.max
    if(tmp.max > tmp.hist.max) {
      tmp.hist.max = tmp.max
    }
    complexity.stats.df$HistMax[i] = tmp.hist.max
    complexity.stats.df$Mean[i] = mean(as.numeric(complexity.df[i,3:ncol(complexity.df)]))
  }

  expt.title = ifelse(expt == 1, "Directed", "Undirected") 
  plots[[length(plots)+1]] = ggplot(complexity.stats.df, aes(x=Time, y=Mean, col=factor(Run))) +
    geom_line() +
    labs(title = paste(c(expt.title, "Mean"))) +
    theme(legend.position = "none") 
  plots[[length(plots)+1]] = ggplot(complexity.stats.df, aes(x=Time, y=Max, col=factor(Run))) +
    geom_line() +
    labs(title = paste(c(expt.title, "Max"))) +
    theme(legend.position = "none") 
  plots[[length(plots)+1]] = ggplot(complexity.stats.df, aes(x=Time, y=HistMax, col=factor(Run))) +
    geom_line() +
    labs(title = paste(c(expt.title, "Historical Max"))) +
    theme(legend.position = "none") 
}

png("nz-plotblend-3.png", width=800, height=800)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], nrow=2)
dev.off()

### fitness dynamics

experiments = c("stats-out-blend-1-0.100-10-16-16-64-1-0-1-5000-1e+08.csv")

plots = list()
expt = 1
df = read.csv(experiments[expt], header=T)
pop.df = read.csv(gsub("stats", "pop", experiments[expt]), header=F)
pop.df = pop.df[pop.df$V1 < 10,]
fitness.df = pop.df
for(i in 3:ncol(fitness.df)) {
  fitness.df[,i] = unlist(lapply(pop.df[,i], function(x) ifelse(x == -1, 0, 1/(abs(df$Size[x+1]-16)+1) ) ))
}
fitness.stats.df = data.frame(Run=fitness.df$V1, Time=fitness.df$V2, Min=0, Mean=0, Max=0, HistMax=0)
tmp.hist.max = 0
for(i in 1:nrow(fitness.df)) {
  if(i > 1) {
    if(fitness.stats.df$Run[i] != fitness.stats.df$Run[i-1]) {
      tmp.hist.max = 0
    }
  }
  fitness.stats.df$Min[i] = min(as.numeric(fitness.df[i,3:ncol(fitness.df)]))
  tmp.max = max(as.numeric(fitness.df[i,3:ncol(fitness.df)]))
  fitness.stats.df$Max[i] = tmp.max
  if(tmp.max > tmp.hist.max) {
    tmp.hist.max = tmp.max
  }
  fitness.stats.df$HistMax[i] = tmp.hist.max
  fitness.stats.df$Mean[i] = mean(as.numeric(fitness.df[i,3:ncol(fitness.df)]))
}

expt.title = "Fitness"
plots[[length(plots)+1]] = ggplot(fitness.stats.df, aes(x=Time, y=Mean, col=factor(Run))) +
  geom_line() +
  labs(title = paste(c(expt.title, "Mean"))) +
  theme(legend.position = "none") 
plots[[length(plots)+1]] = ggplot(fitness.stats.df, aes(x=Time, y=Max, col=factor(Run))) +
  geom_line() +
  labs(title = paste(c(expt.title, "Max"))) +
  theme(legend.position = "none") 
plots[[length(plots)+1]] = ggplot(fitness.stats.df, aes(x=Time, y=HistMax, col=factor(Run))) +
  geom_line() +
  labs(title = paste(c(expt.title, "Historical Max"))) +
  theme(legend.position = "none") 

png("nz-plotblend-4.png", width=800, height=800)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], nrow=2)
dev.off()

#### big trellis plot 1

experiments = c("stats-out-blend-1-0.100-10-16-16-64-0-0-1-5000-1e+08.csv",
                "stats-out-blend-0-0.100-10-16-16-64-1-0-1-5000-1e+08.csv",
		"stats-out-blend-1-0.100-10-16-16-64-0-0-1-0-1e+08.csv",
		"stats-out-blend--1-0.100-10-16-16-64-1-0-1-5000-1e+08.csv")
		
directed.df = read.csv(experiments[1], header=T)
undirected.df = read.csv(experiments[2], header=T)
sampling.df = read.csv(experiments[3], header=T)
random.df = read.csv(experiments[4], header=T)

plot.a = ggplot(sampling.df, aes(x = CMNNonzero, y = log(SampleCount))) + geom_point()
plot.b = ggplot(undirected.df, aes(x = CMNNonzero, y = log(DiscoveryCount))) + geom_point()
plot.c = ggplot(undirected.df, aes(x = log(SampleCount), y = log(DiscoveryCount))) + geom_point()
plot.d = ggplot(random.df, aes(x = CMNNonzero, y = log(DiscoveryCount))) + geom_point()
plot.e = ggplot(undirected.df[undirected.df$Size == 16,], aes(x = CMNNonzero, y = log(DiscoveryCount))) + geom_point()
plot.f = ggplot(directed.df, aes(x = log(SampleCount), y = log(DiscoveryCount))) + geom_point()

png("nz-plotblend-5.png", width=800, height=800)
grid.arrange(plot.a, plot.b, plot.c, plot.d, plot.e, plot.f, nrow=2)
dev.off()

#### search spaces

ss.experiments = c("stats-out-blend-1-0.100-10-16-2-8-0-0-1-0-1e+08.csv",
                   "stats-out-blend-1-0.100-10-16-4-16-0-0-1-0-1e+08.csv",
		   "stats-out-blend-1-0.100-10-16-8-32-0-0-1-0-1e+08.csv",
		   "stats-out-blend-1-0.100-10-16-16-64-0-0-1-0-1e+08.csv")

results.df = data.frame(nb = NULL, stat = NULL, val = NULL)
nbs = c(2, 4, 8, 16)
for(i in 1:4) {
  df.list = read.csv(ss.experiments[i], header=T)
  results.df = rbind(results.df, data.frame(nb = nbs[i], stat = "ST", val = df.list$SampleCount[which(df.list$Size == 1)]))
  results.df = rbind(results.df, data.frame(nb = nbs[i], stat = "SX", val = df.list$SampleCount[which(df.list$Size == 5 & df.list$Symmetry == 0)]))
  DX.count = df.list$SampleCount[which(df.list$Size == 9 & df.list$Symmetry == 0)]
  if(length(DX.count) == 0) {DX.count = 0 }
  results.df = rbind(results.df, data.frame(nb = nbs[i], stat = "DX", val = DX.count))
  results.df = rbind(results.df, data.frame(nb = nbs[i], stat = "D4", val = sum(df.list$SampleCount[df.list$Symmetry == 0])))
  results.df = rbind(results.df, data.frame(nb = nbs[i], stat = "C4", val = sum(df.list$SampleCount[df.list$Symmetry == 1])))
  results.df = rbind(results.df, data.frame(nb = nbs[i], stat = "D2", val = sum(df.list$SampleCount[df.list$Symmetry == 2])))
  results.df = rbind(results.df, data.frame(nb = nbs[i], stat = "C2", val = sum(df.list$SampleCount[df.list$Symmetry == 3])))
  results.df = rbind(results.df, data.frame(nb = nbs[i], stat = "D1", val = sum(df.list$SampleCount[df.list$Symmetry == 4])))
  results.df = rbind(results.df, data.frame(nb = nbs[i], stat = "C1", val = sum(df.list$SampleCount[df.list$Symmetry == 5])))
}

results.df$stat = factor(results.df$stat)
results.df$val[results.df$val == 0] = NA
structures.plot = ggplot(results.df[results.df$stat %in% c("ST", "SX", "DX"),], aes(x=nb, y=log(val), fill=stat, colour = stat)) + geom_point()
symms.plot = ggplot(results.df[results.df$stat %in% c("D4", "C4", "D2", "C2", "D1", "C1"),], aes(x=nb, y=log(val), fill=stat, colour=stat)) + geom_point()

png("nz-plotblend-6.png", width=800, height=800)
grid.arrange(structures.plot, symms.plot, nrow=2)
dev.off()

#### big trellis plot 2

experiments = c("stats-out-blend-1-0.100-10-16-16-64-0-0-1-5000-1e+08.csv",
		"stats-out-blend-1-0.100-10-16-16-64-0-0-1-0-1e+08.csv")

symmgroups = c("D4", "C4", "D2", "C2", "D1", "C1")

directed.df = read.csv(experiments[1], header=T)
sampling.df = read.csv(experiments[2], header=T)

sampling.df = sampling.df[sampling.df$Size == 16,]

symmgroups = c("D4", "C4", "D2", "C2", "D1", "C1")

sampling.df$SymmGroup = symmgroups[sampling.df$Symmetry+1]
directed.df$SymmGroup = symmgroups[directed.df$Symmetry+1]
symm.df = data.frame(SymmGroup = rep(symmgroups, 2), Class=c(rep("Directed", length(symmgroups)), rep("Sampled", length(symmgroups))), Count=0)
for(i in 1:nrow(sampling.df)) {
  ref = which(symm.df$SymmGroup == sampling.df$SymmGroup[i] & symm.df$Class == "Sampled")
  symm.df$Count[ref] = symm.df$Count[ref] + 1
}
for(i in 1:nrow(directed.df)) {
  ref = which(symm.df$SymmGroup == directed.df$SymmGroup[i] & symm.df$Class == "Directed")
  symm.df$Count[ref] = symm.df$Count[ref] + directed.df$DiscoveryCount[i]
}
symm.df$Count[symm.df$Class == "Directed"] = symm.df$Count[symm.df$Class == "Directed"] / sum(symm.df$Count[symm.df$Class == "Directed"])
symm.df$Count[symm.df$Class == "Sampled"] = symm.df$Count[symm.df$Class == "Sampled"] / sum(symm.df$Count[symm.df$Class == "Sampled"])

symm.plot = ggplot(symm.df ,aes(x=factor(SymmGroup, levels=symmgroups), y=log(Count), colour=Class)) + geom_point()

max.mod = max(max(sampling.df$Modularity), max(directed.df$Modularity))
mod.df = data.frame(Mod = seq(from=1,to=max.mod), Class = c(rep("Directed", max.mod), rep("Sampled", max.mod)), Count=0)
for(i in 1:nrow(sampling.df)) {
  ref = which(mod.df$Mod == sampling.df$Modularity[i] & mod.df$Class == "Sampled")
  mod.df$Count[ref] = mod.df$Count[ref] + 1
}
for(i in 1:nrow(directed.df)) {
  ref = which(mod.df$Mod == directed.df$Modularity[i] & mod.df$Class == "Directed")
  mod.df$Count[ref] = mod.df$Count[ref] + directed.df$DiscoveryCount[i]
}
mod.df$Count[mod.df$Class == "Directed"] = mod.df$Count[mod.df$Class == "Directed"] / sum(mod.df$Count[mod.df$Class == "Directed"])
mod.df$Count[mod.df$Class == "Sampled"] = mod.df$Count[mod.df$Class == "Sampled"] / sum(mod.df$Count[mod.df$Class == "Sampled"])
mod.df$Mod = 1/mod.df$Mod

mod.plot = ggplot(mod.df ,aes(x=log(Mod), y=Count, fill=Class)) + geom_col(position="dodge")

mod.stats.df = data.frame(Mod = NULL, Class = NULL, Complexity=NULL, SymmGroup=NULL)
for(i in 1:nrow(sampling.df)) {
  mod.stats.df = rbind(mod.stats.df, data.frame(Mod = sampling.df$Modularity[i], Class = "Sampled", Complexity = sampling.df$CMNNonzero[i], SymmGroup = symmgroups[sampling.df$Symmetry[i]+1]))
}
#for(i in 1:nrow(directed.df)) {
#  mod.stats.df = rbind(mod.stats.df, data.frame(Mod = directed.df$Modularity[i], Class = "Directed", Complexity = directed.df$CMNNonzero[i], SymmGroup = symmgroups[directed.df$Symmetry[i]+1]))
#}
mod.stats.df$Mod = 1/mod.stats.df$Mod

mod.dist = ggplot(mod.stats.df, aes(x=factor(log(Mod)), y = Complexity, colour=Class)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

symm.dist = ggplot(mod.stats.df, aes(x=factor(SymmGroup, levels=symmgroups), y = Complexity, fill=Class)) + geom_boxplot()

png("nz-plotblend-7.png", width=800, height=800)
grid.arrange(symm.plot, mod.plot, mod.dist, symm.dist, nrow=4)
dev.off()

#### pairs plot

experiments = c("stats-out-blend-1-0.100-10-16-16-64-0-0-1-0-1e+08.csv")

sampling.df = read.csv(experiments[1], header=T)
sampling.df$logP = log(sampling.df$SampleCount)
to.plot = subset(sampling.df, select=c(logP,CMNBlock,CMNNonzero,CMNNonzero,CMLZW,CMNecklace,Size))
colnames(to.plot) = gsub("CM", "", colnames(to.plot))

png("nz-plotblend-8.png", width=800, height=800)
ggpairs(to.plot, lower = list(continuous = "cor", combo = "box_no_facet", discrete = "count", na = "na"), upper = list(continuous = "points", combo = "facethist", discrete = "facetbar", na = "na"))
dev.off()

### complexity dynamics for bigger structures

experiments = c("stats-out-blend-1-0.100-10-32-16-64-1-0-1-10-1e+08.csv",
                "stats-out-blend-1-0.100-10-35-16-64-1-0-1-10-1e+08.csv",
                "stats-out-blend-1-0.100-10-40-16-64-1-0-1-10-1e+08.csv")
		
plots = list()
for(expt in 1:3) {
  df = read.csv(experiments[expt], header=T)
  pop.df = read.csv(gsub("stats", "pop", experiments[expt]), header=F)
  pop.df = pop.df[pop.df$V1 < 100,]
  complexity.df = pop.df
  for(i in 3:ncol(complexity.df)) {
    complexity.df[,i] = unlist(lapply(pop.df[,i], function(x) ifelse(x == -1, 0, df$CMNNonzero[x+1])))
  }
  complexity.stats.df = data.frame(Run=complexity.df$V1, Time=complexity.df$V2, Min=0, Mean=0, Max=0, HistMax=0)
  tmp.hist.max = 0
  for(i in 1:nrow(complexity.df)) {
    if(i > 1) {
      if(complexity.stats.df$Run[i] != complexity.stats.df$Run[i-1]) {
        tmp.hist.max = 0
      }
    }
    complexity.stats.df$Min[i] = min(as.numeric(complexity.df[i,3:ncol(complexity.df)]))
    tmp.max = max(as.numeric(complexity.df[i,3:ncol(complexity.df)]))
    complexity.stats.df$Max[i] = tmp.max
    if(tmp.max > tmp.hist.max) {
      tmp.hist.max = tmp.max
    }
    complexity.stats.df$HistMax[i] = tmp.hist.max
    complexity.stats.df$Mean[i] = mean(as.numeric(complexity.df[i,3:ncol(complexity.df)]))
  }

  expt.title = switch(expt == 1, "s* = 32", "s* = 35", "s* = 40")
  plots[[length(plots)+1]] = ggplot(complexity.stats.df, aes(x=Time, y=Mean, col=factor(Run))) +
    geom_line() +
    labs(title = paste(c(expt.title, "Mean"))) +
    theme(legend.position = "none") 
  plots[[length(plots)+1]] = ggplot(complexity.stats.df, aes(x=Time, y=Max, col=factor(Run))) +
    geom_line() +
    labs(title = paste(c(expt.title, "Max"))) +
    theme(legend.position = "none") 
  plots[[length(plots)+1]] = ggplot(complexity.stats.df, aes(x=Time, y=HistMax, col=factor(Run))) +
    geom_line() +
    labs(title = paste(c(expt.title, "Historical Max"))) +
    theme(legend.position = "none") 
}

png("nz-plotblend-9.png", width=800, height=800)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], nrow=2)
dev.off()
