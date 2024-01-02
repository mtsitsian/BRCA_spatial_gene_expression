#Function to plot chromosomes
human.genome.plot<-function(chromosomes, data, dataname)
{
  #load data like this
  #chromosomes<-read.table("/home/christoforos/Dropbox/Data/hg18annot/hg18_R.genome")[,1:4]
  #data<-read.table("/home/christoforos/Dropbox/Data/hg18annot/hg18_genes.bed")[,1:3]
  chrom<-as.matrix(chromosomes)
  # adding a simple number in case of score-less datasets
  if(dim(data)[2]<4){data<-cbind(data,1)}
  #
  sites<-as.data.frame(data)
  title<-as.character(dataname)
  par(mar=c(1,1,2,1))
  margin<-250;
  width<-50;
  offset<-15;
  x<-seq(1,margin);
  y<-seq(1,width,by=width/margin);
  y<-c(y,1,1,1,1);
  scale<-1000000;
  plot(x, y, type="n", xlab="", ylab="", axes=F, main=title, cex.main=1.6); #plotting an empty frame
  #color coding
  values=sites[,4];
  quantile(-log(values)/log(10), probs=seq(0,1,by=0.05))->quants
  quantile(values[values<=0], probs=seq(0,1,by=0.05))->quants1 # separate quantiles to split color code around 0
  quantile(values[values>0], probs=seq(0,1,by=0.05))->quants2
  quants<-c(quants1,quants2)
  quants<-levels(as.factor(sites[,4]))
  cols<-rainbow(length(quants))
  cols<-colorRampPalette(c("steelblue4","firebrick4"))(length(quants))
  colors<-0;
  # logspace
  for(k in 1:length(quants)){which(-log(sites[,4])/log(10)>=quants[k] & -log(sites[,4])/log(10)<=quants[k+1])->ind; colors[ind]<-cols[k]}
  # normal space
  for(k in 1:length(quants)){which(sites[,4]>=quants[k] & sites[,4]<=quants[k+1])->ind; colors[ind]<-cols[k]}
  # extreme values
  for(k in 1:length(quants)){which(sites[,4]==quants[k])->ind; colors[ind]<-cols[k]}
  for(i in 1:length(chrom[,1])){
    name<-chrom[i,2];
    size<-as.numeric(chrom[i,3])/scale;
    cent<-as.numeric(chrom[i,4])/scale;
    coord1<-width-1-2*i;
    coord2<-width-2*i;
    rect(offset, coord1+1.3, size+offset, coord2-1.3, col="dark grey", border=F);
    segments(cent+offset,coord1+1.3-0.2,cent+offset,coord2-1.3+0.2, lwd=2, col="black");
    text(1, coord1+0.75, labels=name, cex=0.8);
    subset<-which(sites[,1]==name);
    positions<-(as.numeric(sites[subset,2])+as.numeric(sites[subset,3]))/2
    positions<-positions/scale
    pos1<-(as.numeric(sites[subset,2]));
    pos2<-(as.numeric(sites[subset,3]));
    pos3<-ceiling((pos1+pos2)/2);
    pos1<-pos1/scale;
    pos2<-pos2/scale;
    pos3<-pos3/scale;
    cols<-colors[subset]
    if(length(pos1)>0){
      for (j in 1:length(pos1)) {
        	 # rect(pos3[j]+offset, coord1+1, pos3[j]+offset, coord2-0.5, col="black", border="black", lwd=3)
        rect(pos1[j]+offset, coord1+1.3, pos2[j]+offset, coord2-1.3, col=cols[j], border=cols[j])
      }
    }
  }
}
