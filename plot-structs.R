#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1) {
  stop("Need experiment!")
}

#args = c("out-blend-1-0.100-10-16-16-64-0-0-1-5000-1e+08")

library(ggplot2)

ARR = 16
shadow = 0.3
max.t = 1e4
n.plot = 64
x.plot = 8
space = 10

infile = paste(c("lib-", args[1], ".csv.txt"), collapse="")
t = scan(infile, what="character", sep="\n")
struct.df = data.frame(Structure=NULL,x=NULL,y=NULL)
this.struct = 0
this.y = 1
if(length(t) < max.t) { max.t = length(t) }
for(i in 1:max.t) {
  if(!grepl("Library", t[i])) {
    positives = setdiff(seq(1:ARR), gregexpr("[.]", t[i])[[1]])
    for(xs in positives) {
      struct.df = rbind(struct.df, data.frame(Structure=this.struct, x=xs, y=this.y))
    }
    this.y = this.y+1
  } else {
    this.struct = this.struct + 1
    this.y = 1
}
}


statsinfile = paste(c("stats-", args[1], ".csv"), collapse="")

stat.df = read.csv(statsinfile, header=T)
text.df = data.frame(Structure = NULL, Label=NULL)
for(i in 1:nrow(stat.df)) {
  text.df = rbind(text.df, data.frame(Structure = i, Label = paste(c(switch(stat.df$Symmetry[i]+1, "D₄", "C₄ ", "D₂", "C₂", "D₁", "C₁"), " ", sprintf("%.1f%%", 100*stat.df$AdaptCount[i]/sum(stat.df$AdaptCount))), collapse="")))
}

ranking = nrow(stat.df)-rank(stat.df$AdaptCount, ties="random")

xoffset = function(i) {
  return(i %% x.plot)
}
yoffset = function(i) {
  return(floor(i/x.plot))
}

p.struct.df = struct.df[ranking[struct.df$Structure] <= n.plot,]
p.text.df = text.df[ranking[text.df$Structure] <= n.plot,]

outfile = paste(c("structs-", args[1], ".png"), collapse="")

png(outfile, width=800, height=800)
ggplot() +
  geom_tile(data = p.struct.df,
    aes(x=x+space*xoffset(ranking[Structure])+shadow, y=y-space*yoffset(ranking[Structure])-shadow), fill="#AAAAAA") +
  geom_tile(data = p.struct.df,
    aes(x=x+space*xoffset(ranking[Structure]), y=y-space*yoffset(ranking[Structure])), fill="#EEEEEE",colour="#000000", lwd=0.25) +
  geom_text(data = p.text.df, aes(x = space*xoffset(ranking[Structure])+space/4, y = -space*yoffset(ranking[Structure])-1, label=Label)) +
  theme_void()

dev.off()