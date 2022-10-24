if homedir() == "/home/z840"    #for downward compatibility ;)
    localpath::String = "$(homedir())/2019_Lego_Evo/EnigmaEvo/src/";
elseif isfile("$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/loadfuncs.jl")
    localpath::String = "$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/src/";
else                            #othserwise use relative paths
    localpath::String = "src/";
end;

#using Revise    #helps with debugging in REPL (automatically tracks changes eg in files included with "includet" (include and track))
include(localpath*"loadfuncs.jl");

include(localpath*"set_up_params.jl");

initpoolnet::ENIgMaGraph = setuppool(S,lambda,SSprobs,SOprobs);

# EVOLUTIONARY VERSION
@time poolnet,colnet,sprich,rich,pool,mstrength,evolvedstrength,clock,CID,#=intm_evo,=#mutstep,freqe,freqn,events =
    assemblyevo(initpoolnet, rates0, maxits, cm,cn,ce,cpred, diverse, restrict_colonization, logging);


collapsetime = clock[maxits - findall(!iszero,reverse(diff(sprich)))[1]];

lineplot(clock,sprich,xlabel="clock time",ylabel="species richness",title="Species richness of new model")
lineplot(clock,freqe,xlabel="clock time",ylabel="freq eat",title="Frequency of eat interactions of new model")
lineplot(clock,freqn,xlabel="clock time",ylabel="freq need",title="Frequency of need interactions of new model")

#lineplot(clock,meansprichnew,xlabel="clock time",ylabel="species richness",title="Species richness of new model")
#lineplot(clock,meanfreqenew,xlabel="clock time",ylabel="freq eat",title="Frequency of eat interactions of new model")
#lineplot(clock,meanfreqnnew,xlabel="clock time",ylabel="freq need",title="Frequency of need interactions of new model")


R"""
plot($clock,$sprich,type='l')
points($(collapsetime),$(last(sprich)),cex=2,col='red')
"""

#------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------should work until here----------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#


R"""
barplot(cbind($(sum(eb)),$(sum(eb_evo)),$(sum(nb0)),$(sum(nb0_evo)),$(sum(mb)),$(sum(mb_evo))),names=c('eb','ebevo','nb0','nb0evo','mb','mb0evo'))
"""

sum(mutstep)/length(findall(x->x>=4,events))

R"""
plot($(mstrength),type='l',log='y')
"""



















@time sprich,obrich,clock,edgelist,cid,evID,mutstep,freqe,freqn,events = assemblyevo_diverse(edgelist_origin,sID,oID,e_t,n_t,maxits,probmut,cn,ce,cp,div_t);
R"plot($clock,$sprich,type='l')"
evopos = findall(x->floor(x)==4,events);
R"points($(clock[evopos]),$(sprich[evopos]),pch=16,col='red',cex=0.8)"
lineplot(log.(freqn))

reps = 30;
spr = SharedArray{Int64}(maxits,reps);
clk = SharedArray{Float64}(maxits,reps);
@sync @distributed for i=1:reps
    edgelist_origin,sID,oID = intmatrixv5(S,lambda,SSprobs,SOprobs,OOprobs);
    sprich,obrich,clock,edgelist,cid,evID,mutstep,freqe,freqn,events = assemblyevo_diverse(edgelist_origin,sID,oID,e_t,n_t,maxits,probmut,cn,ce,cp,div_t);
    spr[:,i] = sprich;
    clk[:,i] = clock;
end
R"""
plot($(clk[:,1]),$(spr[:,1]),type='l',ylim=c(0,200),xlim=c(1,400),log='x')
"""
for i = 2:reps
    R"lines($(clk[:,i]),$(spr[:,i]))"
end


intm, eb, nb, nb0 = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);






















# e_b,
# n_b,
# i_b,
# m_b,
# n_b0,
# sp_v,
# int_id = preamble_defs(intm);

# @time sprich,rich,clock,CID,intm_evo,mutstep,freqe,freqn,events = assemblyevo(
    # intm,e_b,n_b,i_b,m_b,n_b0,sp_v,int_id,lambda,
    # athresh,nthresh,maxits,probmut,cn,ce,cp);
intm, tp_m, tind_m, mp_m, mind_m = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);
intm_origin = copy(intm);
e_b_origin,
n_b_origin,
i_b_origin,
m_b_origin,
n_b0_origin,
sp_v,
int_id = preamble_defs(intm_origin);
# e_b_evo,
# n_b_evo,
# i_b_evo,
# m_b_evo,
# n_b0_evo,
# sp_v,
# int_id = preamble_defs(intm_evo);

#Rerun *evolved* system with eco-assembly
@time sprich_evoeco,rich_evoeco,clock_evoeco,CID_evoeco,events_evoeco = assemblyeco(
intm_evo,e_b_evo,n_b_evo,i_b_evo,m_b_evo,n_b0_evo,sp_v,int_id,lambda,athresh,nthresh,maxits,cn,ce,cp);

#Rerun from initial state with eco-assembly
@time sprich_eco,rich_eco,clock_eco,CID_eco,events_eco = assemblyeco(
intm_origin,e_b_origin,n_b_origin,i_b_origin,m_b_origin,n_b0_origin,sp_v,int_id,lambda,athresh,nthresh,maxits,cn,ce,cp);

R"""
plot($sprich / $sprich_eco,pch='.',ylim=c(0,1.5))
points($sprich_evoeco / $sprich_eco,pch='.',col='red')
lines(seq(1,$maxits),rep(1,$maxits),lty=3)
"""


degree_orig = sum(e_b_origin,dims=1);
degree_evo = sum(e_b_evo,dims=1);
degree_mut_orig = sum(n_b0_origin,dims=1);
degree_mut_evo = sum(n_b0_evo,dims=1);
R"""
par(mfrow=c(2,2))
hist($(degree_orig),breaks=20,col='gray')
hist($degree_evo,breaks=20,col='gray')
hist($(degree_mut_orig),breaks=20,col='gray')
hist($degree_mut_evo,breaks=20,col='gray')
"""


R"""
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
par(mfrow=c(2,1))
plot($freqe,pch='.',col=pal[1])
plot($freqn,pch='.',col=pal[2])
"""







namespace = "$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/figures/sprich_web.pdf"
R"""
#pdf($namespace,width=10,height=5)
par(mfrow=c(1,1))
plot($clock,$sprich,type='l',lwd=3,xlab='Time',ylab='Sp/Ob richness',ylim=c(0,max($([sprich;rich.-sprich]))),col='black')
lines($clock,$(rich .- sprich),col='gray')
points($clock,$mutstep * $sprich)
#dev.off()
"""


R"""
plot($freqe,$freqn,pch='.')
"""


md = Array{Float64}(undef,maxits);
for t=1:maxits
    md[t] = std(sum(e_b[CID[:,t],CID[:,t]],dims=2));
end



m_origin = copy(e_b_origin);
osort = sortperm(vec(sum(m_origin,dims=2)))
osort2 = sortperm(vec(sum(m_origin[osort,osort],dims=1)))
m_evo = copy(e_b_evo);
esort = sortperm(vec(sum(m_evo,dims=2)))
esort2 = sortperm(vec(sum(m_evo[esort,esort],dims=1)))
R"""
par(mfrow=c(1,3))
image($(m_origin[osort2,osort]))
image($(m_evo[esort2,esort]))
image($(m_evo[findall(isodd,vec(CID[:,maxits])),findall(isodd,vec(CID[:,maxits])) ]))
# image($m_origin - $m_evo)
"""


@time sprich,rich,clock,CID,intm_evo,mutstep,freqe,freqn = assembly(
    intm_evo,e_b_evo,n_b_evo,i_b_evo,m_b,n_b0_evo,sp_v,int_id,tp_m,tind_m,lambda,
    athresh,nthresh,maxits,probmut);







filename = "figures/intm_evo1.pdf"
namespace = smartpath(filename);
R"""
pdf($namespace,height=5,width=10)
"""
#Reorganize to clump objects
objects = deleteat!(findall(x->x=='i',diag(intm_origin)),1);
objectssort = objects[sortperm(vec(sum(m_b[objects,objects],dims=1)),rev=true)];
# objectssort2 = objectssort[sortperm(vec(sum(m_b[objectssort,objectssort].+e_b[objectssort,objectssort].+n_b[objectssort,objectssort],2)),rev=false)];
species = setdiff(collect(1:length(diag(intm_origin))),objects);
speciessort = species[sortperm(vec(sum(e_b[species,species].+n_b[species,species],dims=1)),rev=true)];
speciessort2 = speciessort[sortperm(vec(sum(e_b[speciessort,speciessort].+n_b[speciessort,speciessort],dims=2)),rev=false)];
int_msort = intm_origin[[speciessort2;objectssort],[speciessort2;objectssort]];
int_v = Array{Int64}(undef,length(int_msort[1,:]),length(int_msort[1,:]));
int_v[(LinearIndices(int_msort))[findall(x->x=='a',int_msort)]].=1;
int_v[(LinearIndices(int_msort))[findall(x->x=='n',int_msort)]].=2;
int_v[(LinearIndices(int_msort))[findall(x->x=='i',int_msort)]].=3;
int_v[(LinearIndices(int_msort))[findall(x->x=='m',int_msort)]].=4;
R"""
library(igraph)
library(plotrix)
library(RColorBrewer)
pal=brewer.pal(5,'Set1')
pal[4] = pal[3]
pal[3] = '#ffffff'
num_play = length($(int_v[1,:]))
xx=matrix(as.numeric(as.factor($int_v)),c(num_play,num_play))
xx2=xx;
xx2[which(xx==1)] = pal[1];
xx2[which(xx==2)] = pal[2];
xx2[which(xx==3)] = pal[3];
xx2[which(xx==4)] = pal[4];
#shade made objects
darken <- function(color, factor=1.2){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
}
objects = which(diag(xx)==3)[-1];
for (i in 1:length(objects)) {
    for (j in 1:num_play) {
        col = xx2[objects[i],j];
        xx2[objects[i],j] = darken(col);
        if (length(intersect(objects,j)) == 0) {
            col = xx2[j,objects[i]];
            xx2[j,objects[i]] = darken(col);
            }
        }
    }
par(mar=c(1,1,1,4))
int_types=c('e','n','i','m')
color2D.matplot(xx,extremes=c(1:length(int_types)), border='white', axes=FALSE, xlab='', ylab='',main='',cellcolors=xx2)
legend(x=num_play+1,y=num_play,legend=int_types,pch=22,pt.bg=pal,xpd=TRUE, bty='n')
#text(x=rep(-0.8,length(objects)),y=num_play-objects+0.5,labels='o', xpd=TRUE,cex=0.6)
#text(x=objects-0.5,y=rep(num_play+0.8,length(objects)),labels='o', xpd=TRUE,cex=0.6)
text(x=-0.8,y=num_play-0.8,labels=expression(paste('1',degree)), xpd=TRUE,cex=0.6)
text(x=0.8,y=num_play+0.8,labels=expression(paste('1',degree)), xpd=TRUE,cex=0.6)
dev.off()
"""

filename = "figures/intm_evo2.pdf"
namespace = smartpath(filename);
R"""
pdf($namespace,height=5,width=10)
"""
e_b,
n_b,
i_b,
m_b,
n_b0,
sp_v,
int_id = preamble_defs(intm_evo);

#Reorganize to clump objects
objects = deleteat!(findall(x->x=='i',diag(intm_evo)),1);
objectssort = objects[sortperm(vec(sum(m_b[objects,objects],dims=1)),rev=true)];
# objectssort2 = objectssort[sortperm(vec(sum(m_b[objectssort,objectssort].+e_b[objectssort,objectssort].+n_b[objectssort,objectssort],2)),rev=false)];
species = setdiff(collect(1:length(diag(intm_evo))),objects);
speciessort = species[sortperm(vec(sum(e_b[species,species].+n_b[species,species],dims=1)),rev=true)];
speciessort2 = speciessort[sortperm(vec(sum(e_b[speciessort,speciessort].+n_b[speciessort,speciessort],dims=2)),rev=false)];
int_msort = intm_evo[[speciessort2;objectssort],[speciessort2;objectssort]];
int_v = Array{Int64}(undef,length(int_msort[1,:]),length(int_msort[1,:]));
int_v[(LinearIndices(int_msort))[findall(x->x=='a',int_msort)]].=1;
int_v[(LinearIndices(int_msort))[findall(x->x=='n',int_msort)]].=2;
int_v[(LinearIndices(int_msort))[findall(x->x=='i',int_msort)]].=3;
int_v[(LinearIndices(int_msort))[findall(x->x=='m',int_msort)]].=4;
R"""
library(igraph)
library(plotrix)
library(RColorBrewer)
pal=brewer.pal(5,'Set1')
pal[4] = pal[3]
pal[3] = '#ffffff'
num_play = length($(int_v[1,:]))
xx=matrix(as.numeric(as.factor($int_v)),c(num_play,num_play))
xx2=xx;
xx2[which(xx==1)] = pal[1];
xx2[which(xx==2)] = pal[2];
xx2[which(xx==3)] = pal[3];
xx2[which(xx==4)] = pal[4];
#shade made objects
darken <- function(color, factor=1.2){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
}
objects = which(diag(xx)==3)[-1];
for (i in 1:length(objects)) {
    for (j in 1:num_play) {
        col = xx2[objects[i],j];
        xx2[objects[i],j] = darken(col);
        if (length(intersect(objects,j)) == 0) {
            col = xx2[j,objects[i]];
            xx2[j,objects[i]] = darken(col);
            }
        }
    }
par(mar=c(1,1,1,4))
int_types=c('e','n','i','m')
color2D.matplot(xx,extremes=c(1:length(int_types)), border='white', axes=FALSE, xlab='', ylab='',main='',cellcolors=xx2)
legend(x=num_play+1,y=num_play,legend=int_types,pch=22,pt.bg=pal,xpd=TRUE, bty='n')
#text(x=rep(-0.8,length(objects)),y=num_play-objects+0.5,labels='o', xpd=TRUE,cex=0.6)
#text(x=objects-0.5,y=rep(num_play+0.8,length(objects)),labels='o', xpd=TRUE,cex=0.6)
text(x=-0.8,y=num_play-0.8,labels=expression(paste('1',degree)), xpd=TRUE,cex=0.6)
text(x=0.8,y=num_play+0.8,labels=expression(paste('1',degree)), xpd=TRUE,cex=0.6)
dev.off()
"""







extinctions = findall(x->x==-1,diff(sprich));
dt = diff(clock)[extinctions];
extrate = 1 ./ dt;
R"hist(1/$dt,col='red',xlim=c(20,100),add=TRUE)"

persistance = sum(CID[2:S,:],dims=2) ./ maxits;
neats = sum(e_b[2:S,2:S],dims=2);
npreds = sum(e_b[2:S,2:S],dims=1);
nneeds = sum(n_b[2:S,2:S],dims=2);
comp = vec(neats) .- vec(nneeds) .- vec(npreds)
trophic = trophicalc(2:S,tind_m);

filename = "figures/persistance.pdf"
namespace = smartpath(filename);
R"""
pdf($namespace,height=8,width=8)
par(mfrow=c(2,2))
plot($neats,$persistance,pch=16,cex=1)
plot($npreds,$persistance,pch=16,cex=1)
plot($nneeds,$persistance,pch=16,cex=1)
plot($trophic,$persistance,pch=16,cex=1)
dev.off()
"""


tstep = maxits;
cid = findall(isodd,CID[:,tstep]);
deg,troph = structure(S,cid,sp_v,tind_m);
spcid = intersect(sp_v,cid);
spcid_ind = indexin(spcid,[1;sp_v]);
#Degree distribution
# degrees = vec(sum(tind_m[spcid_ind,spcid_ind],2));
adjmatrix = tind_m[[1;spcid_ind],[1;spcid_ind]];
indmatrix = adjmatrix .- tp_m[[1;spcid_ind],[1;spcid_ind]];
dirmatrix = tp_m[[1;spcid_ind],[1;spcid_ind]];


connectance = sum(adjmatrix)/(sprich[maxits]^2);
A = nichemodelweb(sprich[maxits],connectance);
adjmatrix_niche = A[1];

R"""
par(mfrow=c(1,2))
image($(adjmatrix_niche));
image($adjmatrix);
"""

outdegrees_niche = sum(adjmatrix_niche,dims=1)
outdegrees = sum(adjmatrix,dims=1)

indegrees_niche = sum(adjmatrix_niche,dims=2);
indegrees = sum(adjmatrix,dims=2);

R"""
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
par(mfrow=c(1,2))
hist($outdegrees_niche,xlim=c(0,20),col=paste(pal[2],50,sep=''),breaks=10)
hist($outdegrees,add=TRUE,col=paste(pal[1],50,sep=''),breaks=10)
hist($indegrees_niche,xlim=c(0,20),col=paste(pal[2],50,sep=''),breaks=10)
hist($indegrees,add=TRUE,col=paste(pal[1],50,sep=''),breaks=10)
"""

# g = DiGraph(adjmatrix');
# paths = gdistances(g,1);
# keepnodes = findall(x->x<N+1,paths);

#Visualize graph
#namespace = "$(homedir())/Dropbox/PostDoc/2019_Lego_Evo/EnigmaEvo/figures/webnet.pdf"
R"""
library(igraph)
library(RColorBrewer)
#pdf($namespace,width=6,height=5)
pal <- brewer.pal(3,"Set1")
fw_g <- graph.adjacency($(adjmatrix'));
basal_pos <- 1
trophic = as.numeric($([0;troph[1:size(adjmatrix)[1]-1]]));
#trophic = as.numeric($([0;paths[keepnodes[2:length(keepnodes)]]]));
keepnodes = c(1,which(trophic>0.9))"""; @rget keepnodes; keepnodes = Int64.(keepnodes);
R"""
#keepnodes = $keepnodes;
trophic2 = trophic[keepnodes];
coords <- cbind(runif(length(keepnodes)),trophic2);
coords[basal_pos,1] <- 0.5
fw_g = graph.adjacency($(adjmatrix[keepnodes,keepnodes]'))
plot(fw_g,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color='#6495ED',vertex.label=NA,vertex.frame.color=NA, vertex.color=c(pal[1],rep(pal[2],vcount(fw_g)-1)))
#main=ecount(fw_g)/$(size(adjmatrix)[1])^2,
fw_ind <- graph.adjacency($(indmatrix[keepnodes,keepnodes]'));
#plot(fw_ind,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color='red',vertex.label=NA,vertex.frame.color=NA, vertex.color=c(pal[1],rep(pal[2],vcount(fw_g)-1)),add=TRUE)
dev.off()
"""



R"""
library(bipartite)
nest = networklevel($e_b,index="NODF")
"""
@rget nest;
