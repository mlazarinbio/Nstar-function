
# require vegan 
library(vegan)



Nstar<- function (datum,nset=0,rep=500,ExpAll=FALSE,perms=1000,method = "exact")
{

#Nstar Main
#___________________________________________________
#Subroutines
#___________________________________________________
#function to find the minimum N when I > S for the first time
#_____________________________________________________
NstarS<-function(n,sn,a)
{
for( i in n:2) 
   {     if (i*a/2 <= sn[i]) {
           break     }
   }
if (i==n) 
{NstarS<- -i} 
else 
{
b1<-sn[i+1]-sn[i]
a1<-sn[i+1]-b1*(i+1)
a2<-((a/2) -b1)
if (a2==0) {Nstar<- -1} else {NstarS<-a1/a2}

}
}

#Function to estimate Nstar of a random pattern
#___________________________________________________
Rnstar<-function(bw){
	# start from an approximate value
	f<-1/bw
	nrt<-0.5+1.6/f
	repeat{
	nrt1= (2/f) *(1-(1-f)^nrt)
	if (abs(nrt1-nrt)<0.001) {break} else {nrt<-nrt1}
		}
return (nrt)
}

#___________________________________________________

EstNstarS<-function(datum,nset ,dName,plot ,method ,perms)
{
#select a random subset of nset samples
if (nset>0) 
{Idx <- sample(nrow(datum), size=nset)
datum<-datum[Idx,]}
# Species accumulation curve
x<-specaccum(datum,  method,  perms,TRUE,  "jack1")

#get Richness and sd in separate vectors
sn<- x$richness
sd<-x$sd
N<-length(sn)
if (N<3){
   return ("N* is not estimable with fewer than 3 samples")
   break}
#set Non estimable SD to zero
sd[sd=="NaN"]=0
#alpha
a<-sn[1]
#estimate Whittaker and Harrisons beta, est N* of a random data set

bw<-sn[N]/a

bh<-(bw-1)/(N-1)
nsrand<-Rnstar(bw)
nc<-1
# lower and upper bounds for Sn
sl<-sn-sd*nc
sh<- sn+sd*nc
Nav<-NstarS(N,sn,a)
Nmin<-NstarS(N,sl,a)
Nmax<-NstarS(N,sh,a)
z<-floor(Nav)

Sns<-sn[z]+(sn[z+1]-sn[z])*Nav
#round the values
Nmin<- round(Nmin, 3)
Nmax<-round(Nmax, 3)
Nav<-round(Nav, 3)

#prepare output
cells <- c(N,round(sn[N],3),round(sn[1],3),Nmin,Nav,Nmax, round(nsrand,3),round(Sns,3),round(bw,3),round(bh,3))
rnames <- c("Samples","Species","mean_alpha","NL", "Naver","NU","N*_Rnd","SN*","Whittaker","Harrison")
cnames<-(c(dName)) 
if (Nav > 0) {
  y <- matrix(cells, nrow=1, ncol=10, byrow=FALSE,
  dimnames=list(cnames, rnames))
   } else {
	y <-"N* > N"}
# Should make a plot?
if (plot) {
#set plot limit
	if (N>30){lmt<-1.3*Nav}else{lmt<-N}
	lmt1<-sn[lmt]
#plot 
	plot(x,  ci = nc, ci.type = "bar", col = par("fg"), ci.col = 1, ci.lty = 1, xlab = "Number of sampling units (N)", ylab = "Number of Species (S)",xlim=range(1,lmt), ylim=range (1,lmt1))
	abline(0,a/2,col="red")
	abline(v=Nmin,col=4,lty=3)
	abline(v=Nav,col=4,lty=1)
	abline(v=Nmax,col=4,lty=3)
}
#remove temp objects
rm(x)
rm (datum)
rm(sn)
rm(sd)
rm(sl)
rm(sh)
return(y)
}

#___________________________________________________
#Main
#___________________________________________________
#get the name of the dataset
dName<-toString( match.call())
dName <- strsplit(dName,",")
dName<-sapply(dName,"[",2)

#Use subsets?
if (nset > 2) 
{Idx <- sample(nrow(datum), size=nset)
datumr<-datum[Idx,]}
else
{
# use all samples 
  y<- EstNstarS (datum,0,dName,TRUE,method,perms)
return (y)
   break}
#use subsets

if (rep<2) 
# witout repetition
{ y<- EstNstarS (datum,nset,dName,TRUE,method,perm)
  return (y)
  break} 
#repeat the estimation  rep times
Ntot<-0
Nrtot<-0
Nt2<-0
Nrt2<-0
sa<-0
rp<-0
 y <- matrix(rep, nrow=rep, ncol=1)
 yr <- matrix(rep, nrow=rep, ncol=1)
 ya<-matrix(rep, nrow=rep, ncol=1)
for (i in 1:rep)
{
Idx <- sample(nrow(datum), size=nset)
datumr<-datum[Idx,]
x<-specaccum(datumr, method , perms,TRUE,  "jack1")

sn<- x$richness
#N<-length(sn)
N<-nset
a<-sn[1]

bw<-sn[N]/a

Nav<-NstarS(N,sn,a)
Nrav<-Rnstar(bw)
if (Nav>0) {

rp<-rp+1
y[rp]<-Nav
yr[rp]<-Nrav
ya[rp]<-a
}
}
rm (datumr)
y<-y[1:rp]
yr<-yr[1:rp]
ya<-ya[1:rp]
Ntot<-Ntot/rp
Nrtot<-Nrtot/rp
#remove temp objects
rm(x)
rm (datum)
rm(sn)
#structure to export
if (ExpAll==TRUE)
{
return (y)}
else
{
# ratio N*/Nr
yr<-y/yr
#N* statistics
#med<-median(y)
ci<-sqrt(var(y))

Nav<-mean(y)
q<-quantile(y)
q2<-q[2]
q4<-q[4]
med<-q[3]

#Nr statistics

Nrav<-mean(yr)
cir<-sqrt(var(yr))
q<-quantile(yr)
q2r<-q[2]
q4r<-q[4]
medr<-q[3]

#alpha statistics
sa<-mean(ya)
cia<-sqrt(var(ya))
q<-quantile(ya)
q2a<-q[2]
q4a<-q[4]
meda<-q[3]


rp<-100*rp/rep
#c1<-c(med,Nav,ci95,ci99)
c1<-c(Nav,ci,med,q2,q4)
c2<-c(Nrav,cir,medr,q2r,q4r)
c3<-c(sa,cia,meda,q2a,q4a)
cells <- c(c1,c2,c3)
y <- matrix(cells, nrow=3, ncol=5, byrow=TRUE)
colnames (y) <- c("Average","std","Median","Q(25%)","Q(75%)")
rownames(y)<-c("N*","Nr","alpha")

ex<-list("data"=dName,"% valid estimations"=rp,"Results"=y)
return (ex)
}
}



