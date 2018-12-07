## surveillance guided rx (i.e. use specific abx according to probability of being effective * averted mortality for each disease)
#, or CRP testing without or with surveillance data
#smart surveillance - choice of abx based on expected value of mortality (prv*mr)

strt<-c("ET","CRP","SV","CRP+SV")
qrt<-c("Q1","Q2","Q3","Q4")
dis=c("Flu","Dng","JEV","ST","Lpt","Bact")

pop<-880000
vhw<-1500
#predicted prev of each dis in population of 100k (inc fever of 33 per 100 person-years from Capeding et al 2013).
pr_l <- read.table("qrt_inc5.tab",header=T)
pdx<-sum(pr_l[1:6,]) #proportion of total inc with confirmed dx
allinc <- 0.33 #inc fever pppa
p_att<-1 #probability febrile patient seeks care with HCW
inc_l<- ceiling(allinc*pop*p_att*pr_l)

rownames(inc_l)<-c(dis,'No Dx')
colnames(inc_l)<-qrt 

#treatment options
opt<-c("NoRx","Amox","Doxy")
crp_h<-c(0.20,0.12,0.42,0.70,0.81,0.84,0.36) #prob crp high with a threshold of 40mg/l based on data from Lubell et al, BMC ID 2015

p_abx<-0.49 #probability abx in those without known infection (check laos data)
etmt<-read.table("etlaos.tab",header=F)
etmt<-etmt*100
rownames(etmt)<-c(dis,'No Dx')
colnames(etmt)<-opt 
etall<-colSums((etmt))
p_bl<- etall[2]/sum(etall[2:3]) #probability of beta-lactam (currently grouping all non-doxy as bl) in those prescribed an abx in absence of surveillance


#cost of drug (mention c_AMR but leave out of main analysis?)
c_am <- 2
c_dx <- 2
cABX <-c(0,c_am,c_dx) #vector for costs of antibiotics in each treatment option (norx, bl, dx)

# 1 in smpn get tested for surveillance data
smpn<-50

c_crp <-2 #cost crp test
# Cost of surveillance system assuming 0.33 fevers ppa, 1/50 fevers tested, cost/test=179, and mobile phones and anciallry costs for VHWs of 200 annually)
c_sv <-ceiling(((179*(sum(inc_l))/smpn)+vhw*200)/sum(inc_l)) 


q<-4

sim=30 #Number of scenarios
#Matrices for costs and deaths by quarter within each simulated scenario
sim_qrt_cost_et<-matrix(rep(0,4*sim),sim,4)
sim_qrt_cost_crp<-matrix(rep(0,4*sim),sim,4)
sim_qrt_cost_sv<-matrix(rep(0,4*sim),sim,4)
sim_qrt_cost_sv_crp<-matrix(rep(0,4*sim),sim,4)

sim_qrt_mr_bl<-matrix(rep(0,4*sim),sim,4)
sim_qrt_mr_et<-matrix(rep(0,4*sim),sim,4)
sim_qrt_mr_crp<-matrix(rep(0,4*sim),sim,4)
sim_qrt_mr_sv<-matrix(rep(0,4*sim),sim,4)
sim_qrt_mr_sv_crp<-matrix(rep(0,4*sim),sim,4)

#Vectors for proportion treated
ptrt_et<- rep(0,4)
ptrt_crp<- rep(0,4)
ptrt_sv<- rep(0,4)
ptrt_crp_sv<- rep(0,4)

##vectors for total annual mortality at baseline and morality averted with each strategy, and costs of strategies
ann_mr_bl<- rep(0,sim)
ann_ma_et <- rep(0,sim)
ann_ma_crp <- rep(0,sim)
ann_ma_sv <- rep(0,sim)
ann_ma_crp_sv <- rep(0,sim)

ann_nc_et <- rep(0,sim)
ann_nc_crp <- rep(0,sim)
ann_nc_sv <- rep(0,sim)
ann_nc_crp_sv <-rep(0,sim)

# proportion of times SV was correct
psvcr <- rep(0,sim)
svtrt<-rep(0,4*sim) #optimal treatment by quarter

for (s in 1:sim){
inc_s<-inc_l
#Mortality rate for untreated infection (randomly sampled from beta distributions)
mr_untrx=c(rbeta(1,1,999),rbeta(1,1,999),rbeta(1,1,999),rbeta(1,6,94),rbeta(1,2.2,97.8),rbeta(1,15,85),rbeta(1,0.5,99.5))
#Matrix for mortality rates using above sample values
mrt<-matrix(c(mr_untrx[1],mr_untrx[1],mr_untrx[1],mr_untrx[2],mr_untrx[2],mr_untrx[2],mr_untrx[3],mr_untrx[3],mr_untrx[3], mr_untrx[4],mr_untrx[4],0, mr_untrx[5],0,0, mr_untrx[6],0,mr_untrx[6],mr_untrx[7],mr_untrx[7],mr_untrx[7]),7,3,byrow=TRUE) #matrix for corrext treatment.
colnames(mrt)<-opt
rownames(mrt)<-c(dis,'No Dx')
  
for (q in 1:4){
  inc_s[1,q]<-ceiling(rgamma(1,0.01*inc_l[1,q],0.01))
  inc_s[2,q]<-ceiling(rgamma(1,0.01*inc_l[2,q],0.01))
  inc_s[3,q]<-ceiling(rgamma(1,0.01*inc_l[3,q],0.01))
  inc_s[4,q]<-ceiling(rgamma(1,0.01*inc_l[4,q],0.01))
  inc_s[5,q]<-ceiling(rgamma(1,0.01*inc_l[5,q],0.01))
  inc_s[6,q]<-ceiling(rgamma(1,0.01*inc_l[6,q],0.01))
}

#simulated number of febrile cases per quarter 
poolq1<-rep(dis,inc_s[1:6,1]) 
poolq2<-rep(dis,inc_s[1:6,2])
poolq3<-rep(dis,inc_s[1:6,3])
poolq4<-rep(dis,inc_s[1:6,4])
inq<-c(length(poolq1),length(poolq2),length(poolq3),length(poolq4))

ttl_nc_et<-rep(0,4)
ttl_nc_sv<-rep(0,4)
ttl_nc_crp<-rep(0,4)
ttl_nc_crp_sv<-rep(0,4)

ttl_mr_bl<-rep(0,4)
ttl_mr_et<-rep(0,4)
ttl_mr_sv<-rep(0,4)
ttl_mr_crp<-rep(0,4)
ttl_mr_crp_sv<-rep(0,4)
svcorr<-rep(0,4)

for (qt in 1:4){
  
  inc <- inc_s[1:6,qt]
  incp<-inc/sum(inc)
  
  
  # random sample of smpn out of all fevers in that quarter
  svp<-sample(rep(c(dis,'No Dx'),inc_s[,qt]),inq[qt]/smpn, replace=TRUE) 
  # choice of drug given sample
  svc<- which.min(c(sum(mrt[svp,'Amox']),sum(mrt[svp,'Doxy'])))
  # actual optimal choice
  svac<- which.min(c(sum(mrt[rep(dis,inc_s[1:6,qt]),'Amox']),sum(mrt[rep(dis,inc_s[1:6,qt]),'Doxy'])))
  svcorr[qt]<-svac==svc
  
  
  simp<-inq[qt]
  s_pt<-rep(0,simp) #vector for number of patients in the qrt
  sim_crp<-rep(0,simp) #vector for simulated crp low/high
  s_et<-rep(0,simp) #vector for simulated trt option for each strategy
  s_crp<-rep(0,simp)
  s_sv<-rep(0,simp)
  s_crp_sv<-rep(0,simp)
  
  nc_et<-rep(0,simp)#vectors for net cost per simulated patient 
  nc_crp<-rep(0,simp)
  nc_sv<-rep(0,simp)
  nc_crp_sv<-rep(0,simp)
  
  #vectors for mortality outcome per stratgey
  mr_bl<-rep(0,simp)
  mr_et<-rep(0,simp)
  mr_sv<-rep(0,simp)
  mr_crp<-rep(0,simp)
  mr_crp_sv<-rep(0,simp)
  
  for (i in 1:inq[qt]){
        
    s_pt[i]<-sample(rep(dis,inc_s[1:6,qt]),1, replace=TRUE) #simulated patient
    mr_bl[i]<-mrt[s_pt[i],1] #mortality in absence of treatment
    sim_crp[i]<-runif(1,0,1)<crp_h[grep(s_pt[i],dis)] #CRP level in simulated patient, FALSE=Low
    
    #simulated treatment decision
    s_et[i]<-sample(rep(opt,etmt[which(s_pt[1]==dis),]),1)
    #outcome of case/decision (net cost) from matrix
   
    
    nc_et[i]<-cABX[which(opt==s_et[i])]
    mr_et[i]<-mrt[s_pt[i],s_et[i]]
    
    s_crp[i]<-ifelse(sim_crp[i]==FALSE, opt[1],ifelse(runif(1,0,1)<p_bl,opt[2],opt[3]))
    nc_crp[i]<-c_crp+cABX[which(opt==s_crp[i])]
    mr_crp[i]<-mrt[s_pt[i],s_crp[i]]
    
    #smart sv identifies treatment option with least anticipated MR (i.e. not necessarily most prevalent pathogen/abx)
    s_sv[i]<-ifelse(s_et[i]==opt[1],opt[1],opt[1+svc]) 
     
    nc_sv[i]<- c_sv+cABX[which(opt==s_sv[i])]
    mr_sv[i]<-mrt[s_pt[i],s_sv[i]]
    
    s_crp_sv[i]<-ifelse(sim_crp[i]==FALSE, opt[1],opt[1+svc]) 
    nc_crp_sv[i]<- (c_crp+c_sv)+cABX[which(opt==s_crp_sv[i])]
    mr_crp_sv[i]<-mrt[s_pt[i],s_crp_sv[i]]
   
  } #end patient loop
   #record treatment option under sv
   svtrt[(4*s-4)+qt] <-opt[1+which.max(c(incp[5]*mr_untrx[5]+incp[6]*mr_untrx[6],incp[4]*mr_untrx[4]+incp[5]*mr_untrx[5]))]

  ttl_nc_et[qt]<-sum(nc_et)+(inq[qt]/pdx-inq[qt])*p_abx*c_am
  ttl_nc_crp[qt]<-sum(nc_crp)+(inq[qt]/pdx-inq[qt])*(c_crp+crp_h[7]*c_am) #adding in costs for patients without one of the key pathogens
  ttl_nc_sv[qt]<-sum(nc_sv)+(inq[qt]/pdx-inq[qt])*c_sv+(inq[qt]/pdx-inq[qt])*p_abx*c_am
  ttl_nc_crp_sv[qt]<-sum(nc_crp_sv)+(inq[qt]/pdx-inq[qt])*((c_crp+crp_h[7]*c_am)+c_sv)
  
  ttl_mr_bl[qt]<-sum(mr_bl)
  ttl_mr_et[qt]<-sum(mr_et)
  ttl_mr_crp[qt]<-sum(mr_crp)
  ttl_mr_sv[qt]<-sum(mr_sv)
  ttl_mr_crp_sv[qt]<-sum(mr_crp_sv)

  ptrt_et[qt]<- sum(s_et!="NoRx")
  ptrt_crp[qt]<- sum(s_crp!="NoRx")
  ptrt_sv[qt]<- sum(s_sv!="NoRx")
  ptrt_crp_sv[qt]<- sum(s_crp_sv!="NoRx")
  
  } #end quarter loop  
  ann_mr_bl[s] <- sum(ttl_mr_bl)
  ann_ma_et[s] <- sum(ttl_mr_bl)-sum(ttl_mr_et)
  ann_ma_crp[s] <- sum(ttl_mr_bl)-sum(ttl_mr_crp)
  ann_ma_sv[s] <- sum(ttl_mr_bl)-sum(ttl_mr_sv)
  ann_ma_crp_sv[s] <- sum(ttl_mr_bl)-sum(ttl_mr_crp_sv)
  
  ann_nc_et[s] <- sum(ttl_nc_et)
  ann_nc_crp[s] <- sum(ttl_nc_crp)
  ann_nc_sv[s] <- sum(ttl_nc_sv)
  ann_nc_crp_sv[s] <- sum(ttl_nc_crp_sv)

sim_qrt_cost_et[s,]<-ttl_nc_et[1:4]
sim_qrt_cost_crp[s,]<-ttl_nc_crp[1:4]
sim_qrt_cost_sv[s,]<-ttl_nc_sv[1:4]
sim_qrt_cost_sv_crp[s,]<-ttl_nc_crp_sv[1:4]

sim_qrt_mr_bl[s,]<-ttl_mr_bl[1:4]
sim_qrt_mr_et[s,]<-ttl_mr_et[1:4]
sim_qrt_mr_crp[s,]<-ttl_mr_crp[1:4]
sim_qrt_mr_sv[s,]<-ttl_mr_sv[1:4]
sim_qrt_mr_sv_crp[s,]<-ttl_mr_crp_sv[1:4]

p_trtan_et<-sum(ptrt_et)/sum(inq)
p_trtan_crp<-sum(ptrt_crp)/sum(inq)
p_trtan_sv<-sum(ptrt_sv)/sum(inq)
p_trtan_crp_sv<-sum(ptrt_crp_sv)/sum(inq)

 # prop in which SV was correct
psvcr[s]<- mean(svcorr)

 }


###plot mortality
colst<- c("gray", "red", "orange", "green", "blue")

boxplot(sim_qrt_mr_bl+(1-p_att)*sim_qrt_mr_bl, boxwex = 0.1, at = 1:4 - 0.3, col = colst[1],
        main = "Predicted mortality in rural Savannakhet per strategy", xlab = "Quarter", ylab = "Mortality", ylim=c(0,max(sim_qrt_mr_bl)+max((1-p_att)*sim_qrt_mr_bl)),yaxs = "i")
boxplot(sim_qrt_mr_et+(1-p_att)*sim_qrt_mr_bl, add = TRUE, boxwex = 0.1, at = 1:4 - 0.2, col =  colst[2],xaxt="n")
boxplot(sim_qrt_mr_crp+(1-p_att)*sim_qrt_mr_bl, add = TRUE, boxwex = 0.1, at = 1:4 - 0.1, col = colst[3], xaxt="n")
boxplot(sim_qrt_mr_sv+(1-p_att)*sim_qrt_mr_bl, add = TRUE, boxwex = 0.1, at = 1:4, col = colst[4], xaxt="n")
boxplot(sim_qrt_mr_sv_crp+(1-p_att)*sim_qrt_mr_bl, add=TRUE, boxwex = 0.1, at = 1:4+.1, col = colst[5], xaxt="n")

legend("topleft",c("Baseline", "Current practice", "CRP testing", "Surveillance", "CRP+Surveillance"), fill = colst,cex=0.85, ncol=2, bty="n", x.intersp=0.4, y.intersp=0.7)
legend("topright",dis, text.col= rainbow(6), lty=c(1:6), lwd=2,col=rainbow(6),cex=1, y.intersp=0.5, bty="n")
par(new=TRUE)
matplot(t(inc_l[1:6,]),type="l", col=rainbow(6), lty=c(1:6), lwd=2,axes=F,xlab=NA,ylab=NA,xlim = c(0.5,4.5))
axis(4, ylim=c(0,max(1.5*inc_l)))
mtext("Quarterly incidence", side=4, line=-1)

sim_mr_mean<-c(mean(rowSums(sim_qrt_mr_bl)),mean(rowSums(sim_qrt_mr_et)), mean(rowSums(sim_qrt_mr_crp)), mean(rowSums(sim_qrt_mr_sv)), mean(rowSums(sim_qrt_mr_sv_crp)))
sim_nc_mean<-c(mean(rowSums(sim_qrt_cost_et)), mean(rowSums(sim_qrt_cost_crp)), mean(rowSums(sim_qrt_cost_sv)), mean(rowSums(sim_qrt_cost_sv_crp)))

trtmt<-cbind(s_pt,sim_crp, s_et,s_crp,s_sv,s_crp_sv) #for viewing trt selected in each strt 
# View(trtmt)

YLL<-50
WTP <- 1230 #0.5*GDP/capita
#CEA matrix
ceamt<-matrix(c(sim_nc_mean,YLL*c(sim_mr_mean[1]-sim_mr_mean[2],sim_mr_mean[1]-sim_mr_mean[3],sim_mr_mean[1]-sim_mr_mean[4],sim_mr_mean[1]-sim_mr_mean[5])),4,2)
colnames(ceamt)<-c("Cost", "DALYs averted")
rownames(ceamt)<-strt
ceamt<-ceamt[order(ceamt[,2]),] #order by effectiveness
dmt<-c(ceamt[1,1]>0,ceamt[2,1]<ceamt[3,1],ceamt[3,1]<ceamt[4,1],ceamt[4,1]>ceamt[3,1]) #mark dominated interventions (FALSE=dominated)
ceamt<-cbind(ceamt,dmt)
ceamt<-ceamt[ceamt[,3] == 1,] #remove dominated interventions
nic<-length(ceamt[,1]) #number of ICERs to calculate
icer<-rep(0,nic)

icer[1]<-ceamt[1,1]/(ceamt[1,2])
for (ip in 2:nic){
  icer[ip]<-(ceamt[ip,1]-ceamt[1,1])/((ceamt[ip,2]-ceamt[1,2]))
}
ceamt<-cbind(ceamt,icer)
exdmt<-c(1,rep(0,nic-1))
for (ex in 2:nic){
exdmt[ex]<-c(ceamt[ex,4]<=ceamt[nic,4])
}
ceamt<-cbind(ceamt,exdmt)
ceamt<-ceamt[ceamt[,5] == 1,] #remove dominated interventions

plot(YLL*ann_ma_crp,ann_nc_crp, col=colst[3],pch=18, cex=1.2, xlab="DALYs averted compared with baseline", ylab="Total cost (USD)", xlim=c(0,max(YLL*c(ann_ma_et,ann_ma_crp,ann_ma_sv,ann_ma_crp_sv))),ylim=c(0,max(c(ann_nc_et,ann_nc_crp,ann_nc_sv,ann_nc_crp_sv))))
points(YLL*ann_ma_et,ann_nc_et, col=colst[2], pch=18, cex=1.2)
points(YLL*ann_ma_sv,ann_nc_sv, col=colst[4],pch=18, cex=1.2)
points(YLL*ann_ma_crp_sv,ann_nc_crp_sv, col=colst[5],pch=18, cex=1.2)
abline(0,WTP, col="purple",lty=4)
legend("topleft",c("Current practice", "CRP testing", "Surveillance", "CRP+Surveillance"), fill = colst[2:5],cex=0.85, ncol=2, bty="n", x.intersp=0.4, y.intersp=0.7)
legend("bottomright", c(paste("WTP $",WTP, " per DALY averted", sep=""), "Cost-effectiveness frontier"), lty=4, col=c("purple","black"), cex=0.8)
lines(c(0,ceamt[1,2]),c(0,ceamt[1,1]), lty=3, col=colst[which(strt==rownames(ceamt)[1])+1],lwd=2)
text(x=1.5*ceamt[1,2], y=0.9*ceamt[1,1], pos=1, labels=paste("ICER ", rownames(ceamt)[1], " $",ceiling(ceamt[1,1]/(ceamt[1,2]))," per DALY",sep =""), col=c(colst[which(strt==rownames(ceamt)[1])+1]))

# ceamt<-rbind(c(0,0,0,0,0),ceamt) #add in row for baseline
nic<-length(ceamt[,1]) #number of ICERs to calculate
for (ip in 2:nic){
  icrlab <- paste("ICER ", rownames(ceamt)[ip], " $",ceiling(ceamt[ip,4])," per DALY",sep ="")
 lines(c(ceamt[ip-1,2],ceamt[ip,2]), c(ceamt[ip-1,1],ceamt[ip,1]), lty=3, col=colst[which(strt==rownames(ceamt)[ip])+1],lwd=2)
text(x=1.5*ceamt[ip,2], y=0.9*ceamt[ip,1], pos=1, labels=icrlab,col=c(colst[which(strt==rownames(ceamt)[ip])+1]))
}

table(svtrt) #table for frequency of doxy and amox under SV / SV+CRP
ver<-cbind(trtmt[1:20],mr_et[1:20],mr_crp[1:20],mr_sv[1:20],mr_crp_sv[1:20])


####CEAC for CRP+SV vs CRP alone
ceamrt<-cbind(YLL*ann_ma_et,YLL*ann_ma_crp,YLL*ann_ma_sv,YLL*ann_ma_crp_sv,ann_nc_et,ann_nc_crp,ann_nc_sv,ann_nc_crp_sv)
l_ceac<-100
nbmt<-matrix(rep(0,sim*l_ceac),sim,l_ceac)
wtp_ceac <- seq(0,WTP,length.out=l_ceac) # A vector of 100 elements of WTP thresholds ranging from $0-2000

for (l in 1:l_ceac){
nmb_et<-wtp_ceac[l]*ceamrt[,1]-ceamrt[,5]
nmb_crp<-wtp_ceac[l]*ceamrt[,2]-ceamrt[,6]
nmb_sv<-wtp_ceac[l]*ceamrt[,3]-ceamrt[,7]
nmb_crp_sv<-wtp_ceac[l]*ceamrt[,4]-ceamrt[,8]
  
for (nb in 1:sim){nbmt[nb,l]<-which.max(c(nmb_et[nb],nmb_crp[nb],nmb_sv[nb],nmb_crp_sv[nb]))}
}
#Vectors for ceaf
pet<-rep(0,l_ceac)
pcrp<-rep(0,l_ceac)
psv<-rep(0,l_ceac)
pcrp_sv<-rep(0,l_ceac)

for (cf in 1:l_ceac) {
  pet[cf]<-sum(nbmt[,cf]==1)/length(nbmt[,cf])
  pcrp[cf]<-sum(nbmt[,cf]==2)/length(nbmt[,cf])
  psv[cf]<-sum(nbmt[,cf]==3)/length(nbmt[,cf])
  pcrp_sv[cf]<-sum(nbmt[,cf]==4)/length(nbmt[,cf])
}

plot(wtp_ceac,pet, col=colst[2], type="l", ylab="Probability of being cost-effective", xlab="Willingness to pay threshold")
lines(wtp_ceac,pcrp, col=colst[3])
lines(wtp_ceac,psv, col=colst[4])
lines(wtp_ceac,pcrp_sv, col=colst[5])
legend("right",c("Current practice", "CRP testing", "Surveillance", "CRP+Surveillance"), fill = colst[2:5],cex=0.85, ncol=2, bty="n", x.intersp=0.4, y.intersp=0.7)


#summary for text
#In the absence of any treatment total mortality would be
sum(rowSums(inc_l[1:7,])*mr_untrx)/p_att 

#Broken down by pathogen
ceiling(rowSums(inc_l[1:7,])*mr_untrx/p_att)

#optimal treatment in each qrt&simulation
trtqt<-matrix(svtrt,sim,4,byrow=T)
# treatment choice in each qrt (returns error in low numvber of sims if 100% doxy/amox)
trt_sim<-matrix(c(prop.table(table(trtqt[,1])),prop.table(table(trtqt[,2])),prop.table(table(trtqt[,3])) ,prop.table(table(trtqt[,4]))),4,2,byrow=T)
colnames(trt_sim)<-c("Beta-lactam", "Tetreacycline")
rownames(trt_sim)<-qrt
trt_sim

#C/E of CRP and SV vs ET
(sim_nc_mean[2]-sim_nc_mean[1])/((sim_mr_mean[1]-sim_mr_mean[2])*YLL)
(sim_nc_mean[3]-sim_nc_mean[1])/((sim_mr_mean[2]-sim_mr_mean[4])*YLL)

mr_all<-data.frame(c(rowSums(sim_qrt_mr_bl)),rowSums(sim_qrt_mr_et), rowSums(sim_qrt_mr_crp), rowSums(sim_qrt_mr_sv), rowSums(sim_qrt_mr_sv_crp))
colnames(mr_all)<-c('BL',"ET","CRP","SV","CRPSV")

#means and 95% CI by bootstrapping
library(rcompanion)
groupwiseMean(BL~1, data= mr_all, conf= 0.95,digits = 3, R= 10000, boot= TRUE,traditional= FALSE,normal= FALSE,basic= FALSE, percentile  = FALSE,bca= TRUE)
groupwiseMean(ET~1, data= mr_all, conf= 0.95,digits = 3, R= 10000, boot= TRUE,traditional= FALSE,normal= FALSE,basic= FALSE, percentile  = FALSE,bca= TRUE)
groupwiseMean(CRP~1, data= mr_all, conf= 0.95,digits = 3, R= 10000, boot= TRUE,traditional= FALSE,normal= FALSE,basic= FALSE, percentile  = FALSE,bca= TRUE)
groupwiseMean(SV~1, data= mr_all, conf= 0.95,digits = 3, R= 10000, boot= TRUE,traditional= FALSE,normal= FALSE,basic= FALSE, percentile  = FALSE,bca= TRUE)
groupwiseMean((CRPSV)~1, data= mr_all, conf= 0.95,digits = 3, R= 10000, boot= TRUE,traditional= FALSE,normal= FALSE,basic= FALSE, percentile  = FALSE,bca= TRUE)

#For summary of prescribing in each strategy
table(trtmt[,3])
table(trtmt[,4])
table(trtmt[,5])
table(trtmt[,6])
