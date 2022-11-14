function get_extinction_size_distrib_test()
    working = true;
    function isDist(dist,caseNum)
        if sum(dist) != 1
            println("Error: extinction_size_distrib does not return probability distribution in test case $(caseNum).")
            println("\tsum(dist) != 1")
            return false
        end
        return true
    end

    dist = get_extinction_size_distrib(zeros(10),true)
    working &= isDist(dist,1)
    if dist[0] != 1
        working = false;
        println("Error: extinction_size_distrib failed test case 1.")
    end
    dist = get_extinction_size_distrib(1:10,true)
    working &= isDist(dist,2)
    if dist[9] != 1
        working = false;
        println("Error: extinction_size_distrib failed test case 2.")
    end
    dist = get_extinction_size_distrib([7,7,6,5,4,3,3,4,5,4],true)
    working &= isDist(dist,3)
    if dist[[-4,-1,0,2]] != [1,1,2,1.]/5
        working = false;
        println("Error: extinction_size_distrib failed test case 3.")
    end
    return working
    dist = get_extinction_size_distrib([1,2,3,2,3,4,3,2,3,4,5,6,5,4,5,4,3],true)
    working &= isDist(dist,4)
    if dist[[-2,-1,1,2,4]] != [3,1,1,2,1]/8
        working = false;
        println("Error: extinction_size_distrib failed test case 4.")
    end
    return working
end

function testtrophlvls()
	intm = [0 0 0 0 0 0 0 0
			1 2 2 0 0 0 0 0		#lvl 1
			1 0 2 0 0 0 0 0		#lvl 1
			0 1 0 2 0 0 0 0		#lvl 2
			0 0 1 1 2 0 0 0		#lvl 2
			0 0 0 0 1 2 0 0		#lvl 3
			0 0 1 1 1 1 2 0		#lvl 2
			0 0 0 0 1 0 0 2];	#lvl 3
	net = converttoENIgMaGraph(intm);
	troph_lvls = gettrophiclevels(net);
	exact_lvls = [Set{Int}([2,3]),Set{Int}([4,5,7]),Set{Int}([6,8])];
	if troph_lvls != exact_lvls
		println("Test failed: gettrophiclevels is not working as expected.")
		println("\tgettrophiclevels(testnet) == $troph_lvls\n\t!= $exact_lvls, the precomputed exact soulution.")
		return false
	end
	return true;
end

function testtrophlvls(net)
	isworking = testtrophlvls()
	troph_lvls = gettrophiclevels(net)
	for i in eachindex(troph_lvls)
		for j in 1:(i-1)
			intersection = intersect(troph_lvls[i],troph_lvls[j])
			if !isempty(intersection)
				isworking = false
				println("Test failed: gettrophiclevels is not working as expected.")
				println("\tThe elements of $intersection are categorized as having trophic level $i and $j.")
			end 
		end
	end
	return isworking
end

function checkbasicconsistency(net::ENIgMaGraph)
	consistent = true;
	if !(net.estsize == length(net.hasspec) == length(net.hasv) >= net.idmanager.maxid)
		consistent = false;
		println("Inconsistency:\n\testsize == $(net.estsize), length(hasspec) == $(length(net.hasspec)), length(hasv) == $(length(net.hasv)), maxid == $(net.idmanager.maxid)");
	end
	maxid = maximum(keys(net.vert))
	if net.idmanager.maxid < maxid
		consistent = false;
		println("Inconsistency:\n\tmaxid == $(net.idmanager.maxid) < maximum(keys(net.vert)) == $maxid.")
	end
	spec = findall(net.hasspec)
	if !issetequal(spec,net.spec)
		consistent = false;
		println("Inconsistency:\n\t!issetequal(findall(net.hasspec),net.spec)");
		println("\tsetdiff(findall(net.hasspec),net.spec) = $(setdiff(spec,net.spec))")
		println("\tsetdiff(net.spec,findall(net.hasspec)) = $(setdiff(net.spec,spec))")
	end
	vert = findall(net.hasv)
	if !issetequal(vert,keys(net.vert))
		consistent = false;
		println("Inconsistency:\n\t!issetequal(findall(net.hasv),keys(net.vert))");
		println("\tsetdiff(findall(net.hasv),keys(net.vert)) = $(setdiff(vert,keys(net.vert)))")
		println("\tsetdiff(keys(net.vert),findall(net.hasv)) = $(setdiff(keys(net.vert),vert))")
	end
	for (id,vert) in net
		for fid in vert.feed
			if !net.hasv[fid]
				consistent = false;
				println("Inconsistency:\n\t$fid in net[$id].feed but hasv[$fid] == false.")
			elseif !(id in net[fid].eat)
				consistent = false;
				println("Inconsistency:\n\t$fid in net[$id].feed but $(id) not in net[$fid].eat.")
			end
		end
		for eid in vert.eat
			if !net.hasv[eid]
				consistent = false;
				println("Inconsistency:\n\t$eid in net[$id].eat but hasv[$eid] == false.")
			elseif !(id in net[eid].feed)
				consistent = false;
				println("Inconsistency:\n\t$eid in net[$id].eat but $(id) not in net[$eid].feed.")
			end
		end
		for mid in vert.make
			if !net.hasv[mid]
				consistent = false;
				println("Inconsistency:\n\t$mid in net[$id].make but hasv[$mid] == false.")
			elseif !(id in net[mid].need)
					consistent = false;
					println("Inconsistency:\n\t$mid in net[$id].make but $(id) not in net[$mid].need.")
			end
		end
		if !net.hasspec[id] && id != 1
			for nid in vert.need
				if !net.hasspec[nid]
					consistent = false;
					println("Inconsistency:\n\t$nid in net[$id].need but hasspec[$nid] == false.")
				end
				if !net.hasv[nid]
					consistent = false;
					println("Inconsistency:\n\t$nid in net[$id].need but hasv[$nid] == false.")
				elseif !(id in net[nid].make)
					consistent = false;
					println("Inconsistency:\n\t$nid in net[$id].need but $(id) not in net[$nid].make.")
				end
			end
		end

		#check that no double links exist
		if intersect(vert.eat, vert.need) |> !isempty
			consistent = false;
			println("Inconsistency: Double link:")
			println("\tintersect(net[$id].eat, net[$id].need) == $(intersect(vert.eat, vert.need))")
		end
		if intersect(vert.eat, vert.make) |> !isempty
			consistent = false;
			println("Inconsistency: Double link:")
			println("\tintersect(net[$id].eat, net[$id].make) == $(intersect(vert.eat, vert.make))")
		end
		if intersect(vert.need, vert.make) |> !isempty
			consistent = false;
			println("Inconsistency: Double link:")
			println("\tintersect(net[$id].need, net[$id].make) == $(intersect(vert.need, vert.make))")
		end
	end
	return consistent;
end

function checkinitpoolconsistency(poolnet)
	consistent = checkbasicconsistency(poolnet);
	for (id,vert) in poolnet
		for nid in vert.need
			if !poolnet.hasv[nid]
				consistent = false;
				println("Inconsistency:\n\t$nid in poolnet[$id].need but hasv[$nid] == false.")
			end
		end
	end
	return consistent;
end

function checkconsistency(poolnet,colnet)
	consistent = checkbasicconsistency(colnet);
	consistent &= checkbasicconsistency(poolnet)
	if !issubset(colnet.spec,poolnet.spec)
		consistent = false;
		println("Inconsistency:\n\t!issubset(colnet.spec,poolnet.spec)")
		println("\n\tsetdiff(colnet.spec,poolnet.spec) == $(setdiff(colnet.spec,poolnet.spec)).")
	end
	if !issubset(keys(colnet.vert),keys(poolnet.vert))
		consistent = false;
		println("Inconsistency:\n\t!issubset(keys(colnet.vert),keys(poolnet.vert))")
		println("\tsetdiff(keys(colnet.vert),keys(poolnet.vert)) ==\n\t$(setdiff(keys(colnet.vert),keys(poolnet.vert))).")
	end
	for (id,colvert) in colnet
		poolvert = poolnet[id]
		intersection = intersect(poolvert.eat,keys(colnet.vert));
		if !issetequal(colvert.eat, intersection)
			consistent = false;
			println("Inconsistency:\n\t!issetequal(colnet[$id].eat, intersect(poolnet[$id].eat,keys(colnet.vert)))")
			println("\n\tsetdiff(colnet[$id].eat, intersect(poolnet[$id].eat,keys(colnet.vert))) ==\n\t\t$(setdiff(colvert.eat, intersection)).")
			println("\n\tsetdiff(intersect(poolnet[$id].eat,keys(colnet.vert)), colnet[$id].eat) ==\n\t\t$(setdiff(intersection, colvert.eat)).")
		end
		intersection = intersect(poolvert.feed,keys(colnet.vert));
		if !issetequal(colvert.feed, intersection)
			consistent = false;
			println("Inconsistency:\n\t!issetequal(colnet[$id].feed, intersect(poolnet[$id].feed,keys(colnet.vert)))")
			println("\n\tsetdiff(colnet[$id].feed, intersect(poolnet[$id].feed,keys(colnet.vert))) ==\n\t\t$(setdiff(colvert.feed, intersection)).")
			println("\n\tsetdiff(intersect(poolnet[$id].feed,keys(colnet.vert)), colnet[$id].feed) ==\n\t\t$(setdiff(intersection, colvert.feed)).")
		end
		if colnet.hasspec[id]|| id == 1 		
			if !issetequal(colvert.need, poolvert.need)
				consistent = false;
				println("Inconsistency:\n\t!issetequal(colnet[$id].need, poolnet[$id].need)")
				println("\n\tsetdiff(colnet[$id].need, poolnet[$id].need) == $(setdiff(colvert.need, poolvert.need)).")
				println("\n\tsetdiff(poolnet[$id].need, colnet[$id].need) == $(setdiff(poolvert.need, colvert.need)).")
			end
		else
			intersection = intersect(poolvert.need,colnet.spec);
			if !issetequal(colvert.need, intersection)
				consistent = false;
				println("Inconsistency:\n\t!issetequal(colnet[$id].need, intersect(poolnet[$id].need,colnet.spec))")
				println("\n\tsetdiff(colnet[$id].need, intersect(poolnet[$id].need,colnet.spec)) ==\n\t\t$(setdiff(colvert.need, intersection)).")
				println("\n\tsetdiff(intersect(poolnet[$id].need,colnet.spec), colnet[$id].need) ==\n\t\t$(setdiff(intersection, colvert.need)).")
			end
		end
		if !issetequal(colvert.make, poolvert.make)
			consistent = false;
			println("Inconsistency:\n\t!issetequal(colnet[$id].make, poolnet[$id].make)")
			println("\n\tsetdiff(colnet[$id].make, poolnet[$id].make) == $(setdiff(colvert.make, poolvert.make)).")
			println("\n\tsetdiff(poolnet[$id].make, colnet[$id].make) == $(setdiff(poolvert.make, colvert.make)).")
		end			
	end
	return consistent;
end

function testmutation(ce,cn,cm,cf,thorough=false)
	getinteractiontype = ENIgMaGraphs.getinteractiontype;
	spints = [0,1,2];
    obints = [0,1,2,3];
	
    evolutiontable = [[0 0 0 1 1 1 2 2 2 3 3 3];[1 2 3 0 2 3 0 1 3 0 1 2]];
    tallytable = [4.1 4.2 5.1 4.3 4.4 5.2 4.5 4.6 5.3 6.1 6.2 6.3];

	
	#Setup toy network that has basal resource 1 species 2,3 for mutations and object 4 for mutation.
	#Further species 5,6, and object 7 to have all sorts of interactions and 8,9 not colonized species
	function setup_toy_net()
		
		intm = [0 0 0 0 0 0 0 0 0
				1 2 0 0 2 0 3 1 0
				0 0 2 3 0 2 1 0 1
				0 0 2 0 0 0 0 2 0
				0 1 2 1 2 0 3 0 2
				0 2 1 2 1 2 0 0 1
				0 2 0 0 2 0 0 0 0
				1 2 1 3 0 0 0 2 0
				0 1 2 0 0 1 0 0 2];

		pool_net = converttoENIgMaGraph(intm)
		col_net = ENIgMaGraph(pool_net.estsize, pool_net.idmanager)
		basalres = ENIgMaVert();
		addn!(basalres,1);
		addmod!(col_net,1,basalres);

		for id in [2,3,5,6]
			colonize!(pool_net,col_net,id)
		end
		return pool_net,col_net
	end

	if thorough #is setup_toy_net working?
		poolnet,colnet = setup_toy_net()
		isworking = checkconsistency(poolnet,colnet)
		if !isworking
			println("Error: setup_toy_net is creating an inconsistent network.")
		end
	else
		isworking = true;
	end

	#go over all kinds of outer circumstances and types of mutations and check if the result is as expected
	for diverse in [true,false]	

		#spec-spec interaction
		for old_out_int in spints
			for old_in_int in spints
				for change_in_int in [true,false]
					old_int = change_in_int ? old_in_int : old_out_int;
					for new_int in setdiff(spints,old_int)
						poolnet,colnet = setup_toy_net(); #too lazy to write a deepcopy function for ENIgMaGraphs
						mut_spec_pool = poolnet[2];
						mut_spec_col = colnet[2];
						int_spec_pool = poolnet[3];
						int_spec_col = colnet[3];

						if old_in_int == 1							#establish old interactions
							adde!(mut_spec_pool,2,int_spec_pool,3)
							adde!(mut_spec_col,2,int_spec_col,3)
						elseif old_in_int == 2
							addn!(mut_spec_pool,3)
							addn!(mut_spec_col,3)
						end

						if old_out_int == 1
							adde!(int_spec_pool,3,mut_spec_pool,2)
							adde!(int_spec_col,3,mut_spec_col,2)
						elseif old_out_int == 2
							addn!(int_spec_pool,2)
							addn!(int_spec_col,2)
						end

						if thorough	#did it work
							if !(getinteractiontype(poolnet,2,3) == getinteractiontype(colnet,2,3) == old_in_int) || !(getinteractiontype(poolnet,3,2) ==  getinteractiontype(colnet,3,2) == old_out_int)
								isworking = false;
								println("Error: Interactions weren't added properly in testmutation.")
							end
						end

						mutate!(poolnet,colnet,2,3, change_in_int, old_int, new_int, diverse, tallytable, evolutiontable, ce, cn, cm, cf);

						if !checkconsistency(poolnet,colnet)
							isworking = false;
							println("Error: The mutation created an inconsistent network.")
						end

						if thorough	#have any non related interactions been altered?
							secondpool, secondcol = setup_toy_net()
							for i in 1:7
								for j in 1:7
									if (3,2) != (i,j) != (2,3)
										if getinteractiontype(poolnet,i,j) != getinteractiontype(secondpool,i,j) || getinteractiontype(colnet,i,j) != getinteractiontype(secondcol,i,j)
											isworking = false;
											println("Error: An interaction that shouldn't have been touched was altered during mutation with spec.")
										end
									end
								end
								if diverse && ( 3 != i != 2 )
									if getinteractiontype(poolnet,i,10) != getinteractiontype(secondpool,i,2) != 2 || getinteractiontype(colnet,i,10) != getinteractiontype(secondcol,i,2) != 2 || getinteractiontype(poolnet,10,i) != getinteractiontype(secondpool,2,i) || getinteractiontype(colnet,10,i) != getinteractiontype(secondcol,2,i)
										isworking = false;
										println("Error: An interaction with mutated spec that shouldn't have been touched was altered during mutation with spec.")
									end
								end
							end
						end

						#check if result is as expected
						new_in_int, new_out_int = change_in_int ? (new_int,old_out_int == 2 && diverse ? 0 : old_out_int) : (old_in_int, new_int)	#what are the new interactions? (Remember that outgoing needs are dropped in diverse mode)
						if diverse
							if !(getinteractiontype(poolnet,2,3) == getinteractiontype(colnet,2,3) == old_in_int) || !(getinteractiontype(poolnet,3,2) ==  getinteractiontype(colnet,3,2)  == old_out_int)
								isworking = false;
								println("Error in testmutation: The original interaction with spec was changed in the diverse mode.")
							end
							if !(getinteractiontype(poolnet,10,3) == getinteractiontype(colnet,10,3) == new_in_int) || !(getinteractiontype(poolnet,3,10) ==  getinteractiontype(colnet,3,10)  == new_out_int)
								isworking = false;
								println("Error in testmutation: The new interaction with spec wasn't added properly in diverse mode.")
							end
						else
							if ENIgMaGraphs.getdeltastrength(old_int,new_int,change_in_int,ce,cn,cm,cf) >= 0
								if !(getinteractiontype(poolnet,2,3) == getinteractiontype(colnet,2,3) == new_in_int) || !(getinteractiontype(poolnet,3,2) ==  getinteractiontype(colnet,3,2)  == new_out_int)
									isworking = false;
									println("Error in testmutation: The new interaction with spec wasn't properly added in the non-diverse mode.")
								end
							else
								if !(getinteractiontype(poolnet,2,3) == getinteractiontype(colnet,2,3) == old_in_int) || !(getinteractiontype(poolnet,3,2) ==  getinteractiontype(colnet,3,2)  == old_out_int)
									isworking = false;
									println("Error in testmutation: The interaction with a spec was changed in the non-diverse mode although the competition strength decreased.")
								end
							end
						end
					end
				end
			end
		end

		#spec-object interaction
		for old_int in obints
			for new_int in setdiff(obints,old_int)
				poolnet,colnet = setup_toy_net(); #too lazy to write a deepcopy function for ENIgMaGraphs
				mut_spec_pool = poolnet[2];
				mut_spec_col = colnet[2];
				int_obj_pool = poolnet[4];
				int_obj_col = colnet[4];

				if old_int == 1	#setup old interaction
					adde!(mut_spec_pool,2,int_obj_pool,4)
					adde!(mut_spec_col,2,int_obj_col,4)
				elseif old_int == 2
					addn!(mut_spec_pool,4)
					addn!(mut_spec_col,4)
				elseif old_int == 3
					addm!(mut_spec_pool,2,int_obj_pool,4)
					addm!(mut_spec_col,2,int_obj_col,4)
				end
				
				if thorough #did it work
					if !(getinteractiontype(poolnet,2,4) == getinteractiontype(colnet,2,4) == old_int)
						isworking = false;
						println("Error: Interactions weren't added properly in testmutation with object.")
					end
				end

				mutate!(poolnet, colnet, 2, 4, true, old_int, new_int, diverse, tallytable, evolutiontable, ce, cn, cm, cf);
				
				if !checkconsistency(poolnet,colnet)
					isworking = false;
					println("Error: The mutation with an object as interactee created an inconsistent network.")
				end

				if thorough	#have any non related interactions been altered?
					secondpool, secondcol = setup_toy_net()
					for i in 1:7
						for j in 1:7
							if (4,2) != (i,j) != (2,4)
								if getinteractiontype(poolnet,i,j) != getinteractiontype(secondpool,i,j) || getinteractiontype(colnet,i,j) != getinteractiontype(secondcol,i,j)
									isworking = false;
									println("Error: An interaction that shouldn't have been touched was altered during mutation with obj.")
								end
							end
						end

						if diverse && (4 != i != 2)
							if getinteractiontype(poolnet,i,10) != getinteractiontype(secondpool,i,2) != 2 || getinteractiontype(colnet,i,10) != getinteractiontype(secondcol,i,2) != 2 || getinteractiontype(poolnet,10,i) != getinteractiontype(secondpool,2,i) || getinteractiontype(colnet,10,i) != getinteractiontype(secondcol,2,i)
								isworking = false;
								println("Error: An interaction with mutated spec that shouldn't have been touched was altered during mutation with obj.")
							end
						end
					end
				end

				if diverse
					if !(getinteractiontype(poolnet,2,4) == getinteractiontype(colnet,2,4) == old_int)
						isworking = false;
						println("Error in testmutation: The original interaction with object was changed in the diverse mode.")
					end
					if !(getinteractiontype(poolnet,10,4) == getinteractiontype(colnet,10,4) == new_int)
						isworking = false;
						println("Error in testmutation: The new interaction with object wasn't added properly in diverse mode.")
					end
				else
					if ENIgMaGraphs.getdeltastrength(old_int,new_int,true,ce,cn,cm,cf) >= 0
						if !(getinteractiontype(poolnet,2,4) == getinteractiontype(colnet,2,4) == new_int)
							isworking = false;
							println("Error in testmutation: The new interaction with object wasn't properly added in the non-diverse mode.")
						end
					else
						if !(getinteractiontype(poolnet,2,4) == getinteractiontype(colnet,2,4) == old_int)
							isworking = false;
							println("Error in testmutation: The interaction with an object was changed in the non-diverse mode although the competition strength decreased.")
						end
					end	
				end
			end
		end
	end
	return isworking;
end

function test(thorough,random_seed = 0)
	Random.seed!(random_seed);
	include("set_up_params.jl")

	iserrorfree = testmutation(ce,cn,cm,cpred,thorough);
    iserrorfree &= get_extinction_size_distrib_test();

	for logging in [true,false]
		for restrict_colonization in [true,false]
			for diverse in [0,1]
				initpoolnet::ENIgMaGraph = setuppool(S,lambda,SSprobs,SOprobs);

				thorough && (iserrorfree &= checkinitpoolconsistency(initpoolnet));

				poolnet,colnet,sprich,rich,pool,mstrength,evolvedstrength,clock,CID,maxids,globextspec,mutstep,freqe,freqn,freqe_pool,freqn_pool,events =
				assemblyevo(initpoolnet,rates0,maxits,cm,cn,ce,cpred,diverse,restrict_colonization,logging);

				if !checkconsistency(poolnet,colnet)
					iserrorfree = false;
					println("The above inconsistencies where found in the system after assembly with the following parameters:")
					println("\tdiverse == $diverse, restrict_colonization == $restrict_colonization, logging == $logging")
				end
			end
		end
	end

	return iserrorfree;
end

works = true
for seed in 234:250
	works &= test(true,seed)
end