#Random.seed!(83);


thorough = true;
seed = 23;
test(thorough,seed)


for seed in rand(1:10_000,10)
	println("Seed = $seed")
	test(thorough,seed)
end



#used to compare with old version of code
intm,eb,nb,nb0,mb,SSpwp,SOpwp = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);

poolnet = converttoENIgMaGraph(intm);

"cid = $cid"
"col = $col"
"primext = $primext"
"secext = $secext"
freqeold = freqe[it]
freqnold = freqn[it]


colonizers = Int[]; #later used to store ids of potential colonizers
secextspec = Int[]; #later used to store ids of spec that could go extinct secondarily
extobj = Int[];

basalres = ENIgMaVert();
addn!(basalres,1);

colnet = ENIgMaGraph(poolnet.estsize,poolnet.idmanager);
addmod!(colnet,1,basalres);

for id in cid
	if id <= 200
		colonize!(poolnet,colnet,id);
	end
end

getpotcolonizers!(poolnet,colnet,colonizers);
getsecext!(colnet,secextspec,extobj);
primextspec = getprimext(poolnet,colnet,ce,cn,cm,cpred);

freqenew = sum([length(colnet[id].eat) for id in colnet.spec])/length(colnet.vert);
freqnnew = sum([length(colnet[id].need) for id in colnet.spec])/length(colnet.vert);


setdiff(col,colonizers)
setdiff(colonizers,col)

setdiff(secext,secextspec)
setdiff(secextspec,secext)

setdiff(primext,primextspec)
setdiff(primextspec,primext)

freqenew - freqeold
freqnnew - freqnold

cid = [2, 4, 8, 12, 19, 32, 38, 42, 44, 49, 50, 57, 68, 69, 71, 73, 79, 80, 82, 96, 99, 114, 125, 135, 157, 162, 164, 166, 174, 177, 181, 183, 190, 196, 198, 202, 203];

col = [31, 40, 78, 85, 87, 89, 97, 102, 119, 127, 130, 132, 142, 147, 168, 169, 185, 188, 191, 192, 193, 199];

primext = [49, 68, 80, 164];

secext = [114];

freqeold = 0.7297297297297297;
freqnold = 0.05405405405405406;



cid = [2, 4, 8, 12, 31, 32, 38, 42, 44, 50, 54, 57, 59, 66, 69, 70, 71, 73, 77, 78, 79, 80, 82, 85, 91, 96, 99, 114, 121, 123, 125, 129, 132, 135, 142, 148, 151, 155, 157, 158, 161, 162, 166, 168, 169, 174, 177, 181, 182, 183, 188, 190, 192, 195, 196, 197, 198, 199, 201, 202, 203, 205, 206];

col = [7, 10, 11, 33, 37, 40, 46, 49, 53, 61, 68, 72, 74, 83, 84, 87, 88, 89, 94, 97, 106, 112, 119, 127, 130, 138, 147, 153, 160, 164, 170, 175, 180, 184, 185, 186, 191, 193];

primext = [12, 31, 71, 77, 80, 121, 142, 161, 188];

secext = [78, 114, 174];

freqeold = 0.9206349206349206;
freqnold = 0.1111111111111111;



cid = [2, 4, 8, 10, 32, 37, 38, 42, 44, 53, 54, 57, 59, 66, 68, 69, 70, 73, 79, 82, 85, 87, 88, 89, 91, 96, 97, 99, 112, 121, 123, 125, 129, 131, 132, 135, 142, 148, 151, 153, 157, 158, 161, 164, 166, 169, 175, 177, 181, 182, 183, 185, 186, 190, 191, 192, 193, 195, 196, 197, 198, 199, 201, 202, 203, 205, 206];

col = [7, 11, 12, 33, 40, 46, 49, 61, 62, 71, 72, 74, 77, 78, 80, 83, 98, 102, 104, 106, 114, 116, 119, 127, 138, 139, 147, 149, 160, 162, 168, 170, 174, 180, 184, 188];

primext = [38, 54, 89, 121, 125, 142, 148, 164, 169, 183, 185, 191];
secext = [123, 131, 157, 158, 161, 175];
freqeold = 0.8656716417910447;
freqnold = 0.16417910447761194;








#test setuppool (compare to intmatrixv4)

numit = 10000;
freqe = zeros(numit);
freqn = zeros(numit);
mperspec = zeros(numit);

for i = 1:numit
	poolnet = setuppool(S,lambda,SSprobs,SOprobs);
	#freqe[i] = ((sum(length(poolnet[id].eat) for id in poolnet.spec ) - length(poolnet[1].feed))/(length(poolnet) - 1))
	#freqn[i] = ((sum(length(setdiff(poolnet[id].need,1)) for id in poolnet.spec ))/(length(poolnet) - 1))
	#mperspec[i] = (sum(length(poolnet[id].make) for id in poolnet.spec ))/numspec(poolnet);
end

meanfreqenew = mean(freqe)
stdfreqenew = std(freqe)
meanfreqnnew = mean(freqn)
stdfreqnnew = std(freqn)
meanlambda = mean(mperspec)


for i in 1:numit
	intm,eb,nb,nb0,mb,SSpwp,SOpwp = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);
	freqe[i] = sum(eb[2:201,2:end])/(size(intm)[1]-1)
	freqn[i] = sum(nb0[2:201,2:end])/(size(intm)[1]-1)
end

meanfreqeold = mean(freqe)
stdfreqeold = std(freqe)
meanfreqnold = mean(freqn)
stdfreqnold = std(freqn)

#look at degree distribution
degreedists = zeros(30,10000);

for i = 1:10000
	degreedist = zeros(30);
	poolnet = setuppool(S,lambda,SSprobs,SOprobs);
	for id in poolnet.spec
		spec = poolnet[id];
		degree = length(spec.eat) + length(spec.need) + length(spec.make);
		degreedist[degree+1] += 1;
	end
	degreedists[:,i] = degreedist;
end

meandegreedist = mean(degreedists,dims = 2)
meandegreedist ./= sum(meandegreedist)

lineplot(0:15,meandegreedist[1:16],xlabel="degree")
println(meandegreedist[1:20])