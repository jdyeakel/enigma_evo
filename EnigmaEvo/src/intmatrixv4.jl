function intmatrixv4(S, lambda, SSprobs, SOprobs, OOprobs)

		#NOTE: In this version, different sets of probabilities are used to fill up interactions in each matrix quadrant (SS, SO, OO)

    #NOTE In this version, interactions are randomly assigned and there is no inputted distribution for trophic or need interactions

		#The SO quadrant has all e-n-i-m interactions, so do this one first

    SOp_n=SOprobs.p_n;
    SOp_a=SOprobs.p_a;

    #Draw the number of objects made per species
    #Done many times, the mean will be lambda*S
    #A species is an engineer if number of objects > 0
    pdist = Poisson(lambda);
    OpS = rand(pdist,S-1);
    OpS = [0;OpS]; # basal resource does not make anything
    engind = findall(!iszero,OpS);
    E = length(engind);

    #Some of these objects may not be unique
    obline = collect(1:sum(OpS));
    obindpS = zeros(Int64,E,sum(OpS));
    for i = 1:E
        o = sample(obline,OpS[engind][i],replace=false);
        obindpS[i,o] .= 1;
    end

    obindpS = obindpS[:,findall(!iszero,vec(sum(obindpS,dims=1)))];
    O = size(obindpS)[2];

    #Expected size of the system
    # O = convert(Int64,round(lambda*S,0));
    N = S + O;
    spind = collect(1:S);
    obind = collect(S+1:N);
    engind = collect(S-E+1:S);

    SOp_m = (sum(obindpS))/(N^2);
    SOp_i = 1 - (SOp_n + SOp_a + SOp_m);


    # exp_p_m = S*lambda/((S + S*lambda*(1-exp(-1)))^2)

		#SO interactions: E-N-I-M

    SOpwp = (
      na = SOp_n*(SOp_a/(SOp_a+SOp_n+SOp_i+SOp_m)) + SOp_a*(SOp_n/(SOp_a+SOp_i+SOp_n)),
      nn = SOp_n*(SOp_n/(SOp_a+SOp_n+SOp_i+SOp_m)),
      ni = SOp_n*(SOp_i/(SOp_a+SOp_n+SOp_i+SOp_m)) + SOp_i*(SOp_n/(SOp_a+SOp_n+SOp_i)),
      nm = SOp_n*(SOp_m/(SOp_a+SOp_n+SOp_i+SOp_m)) + SOp_m*1, #(SOp_n/SOp_n),
      ia = SOp_i*(SOp_a/(SOp_a+SOp_n+SOp_i)) + SOp_a*(SOp_i/(SOp_a+SOp_i+SOp_n)),
      ii = SOp_i*(SOp_i/(SOp_a+SOp_n+SOp_i)),
      aa = SOp_a*(SOp_a/(SOp_i+SOp_n+SOp_a))
    );

		#SS interactions: E-N-I
		SSp_n = SSprobs.p_n;
    SSp_a = SSprobs.p_a;
		SSp_i = 1 - SSp_a - SSp_n;

		SSpwp = (
      na = SSp_n*(SSp_a/(SSp_a+SSp_n+SSp_i)) + SSp_a*(SSp_n/(SSp_a+SSp_i+SSp_n)),
      nn = SSp_n*(SSp_n/(SSp_a+SSp_n+SSp_i)),
      ni = SSp_n*(SSp_i/(SSp_a+SSp_n+SSp_i)) + SSp_i*(SSp_n/(SSp_a+SSp_n+SSp_i)),
      ia = SSp_i*(SSp_a/(SSp_a+SSp_n+SSp_i)) + SSp_a*(SSp_i/(SSp_a+SSp_i+SSp_n)),
      ii = SSp_i*(SSp_i/(SSp_a+SSp_n+SSp_i)),
      aa = SSp_a*(SSp_a/(SSp_i+SSp_n+SSp_a))
    );

		#OO interactions: N-I
		OOp_n = OOprobs.p_n;
		OOp_i = 1 - OOp_n;

		OOpwp = (
      nn = OOp_n*(OOp_n/(OOp_n+OOp_i)),
      ni = OOp_n*(OOp_i/(OOp_n+OOp_i)) + OOp_i*(OOp_n/(OOp_n+OOp_i)),
      ii = OOp_i*(OOp_i/(OOp_n+OOp_i))
    );

    #Create an empty character array with dimensions equal to the number of players
    #Region 1) Upper left SxS are species-species interactions
    #Region 2) Upper right SxO are species-object interactions
    #Region 3) Lower left OxS are object-species interacions (n only)
    #Region 4) Lower right OxO are object-object interactions (i only)
    intm = Array{Char}(undef,N,N);
    intm .= Ref('0');

    #Assign the make-need interactions
    #Maybe faster way to do this?
    for i = 1:E
        for j=1:O
            if obindpS[i,j] == 1
                intm[engind[i],obind[j]] = 'm';
                intm[obind[j],engind[i]] = 'n';
            end
        end
    end

    #Calculate the distribution of trophic links per species

    #Matrix of trophic-only interactions
    #The number of species in the trophic and mutualistic matrices will all start out the sames size as intm, however, the objects will be trimmed later
    tp_m = zeros(Int64,S,S);
    #Matrix of mutualistic interactions
    mp_m = zeros(Int64,S,S);

    #The first true species (row/col 2) is always a primary producer
    intm[2,1] = 'e'
    tp_m[2,1] = 1;

    #Fill in diagonal
    #Index 1 is the basal resource ('i')
    #Indices 2:S are species ('i')
    #Indices S+1:N are objects ('i')
    diagindices = diagind(intm);

    intm[diagindices[2:S]] .= Ref('n');
		intm[diagindices[S+1:N]] .= Ref('i');

    #Deal with the basal resource
    intm[1,:] .= Ref('i');
    tp_m[1,:] .= 0;
    intm[findall(x->x=='0',intm[:,1]),1] .= Ref('i');


    #NOTE: Interactions BETWEEN SPECIES

    #Eliminate n-a, i-a, a-a, n-m interactions, which are already determined

    # 1) pr_na
    # 2) pr_nn **
    # 3) pr_ni **
    # 4) pr_nm
    # 5) pr_ia
    # 6) pr_ii **
    # 7) pr_aa
    # deleteat!(pw_prob_new,[1,4,5,7])
    # pw_prob_new = pw_prob_init[[1,2,3,5,6,7]];
    pw_prob_new = [SSpwp.na,SSpwp.nn,SSpwp.ni,SSpwp.ia,SSpwp.ii,SSpwp.aa];
    pw_prob_new = pw_prob_new/sum(pw_prob_new);


    prob_line = cumsum(pw_prob_new);



    for i = 2:S
        for j = 2:S
            #Only choose interactions for empty elements (intm = '0')
            #Only do this for the lower triangular part of the matrix
            if intm[i,j] == '0' && i > j

                r_draw = rand()

                #N:A - asymmetric mutualism
                if r_draw < prob_line[1]
                  rr_draw = rand();
                  if rr_draw < 0.5
                      intm[i,j] = 'e';
                      tp_m[i,j] = 1;
                      intm[j,i] = 'n';
                      mp_m[j,i] = 1;
                  else
                      intm[i,j] = 'n';
                      mp_m[i,j] = 1;
                      intm[j,i] = 'e';
                      tp_m[j,i] = 1;
                  end
                end

                #N:N - symmetric mutualism
                if r_draw > prob_line[1] && r_draw < prob_line[2]
                  intm[i,j] = 'n'
                  intm[j,i] = 'n'
                  #Update mutualistic network
                  mp_m[i,j] = 1;
                  mp_m[j,i] = 1;
                end

                #N:I - commensalism
                if r_draw > prob_line[2] && r_draw < prob_line[3]
                    rr_draw = rand();
                    if rr_draw < 0.5
                        intm[i,j] = 'i';
                        intm[j,i] = 'n';
                        mp_m[j,i] = 1;
                    else
                        intm[i,j] = 'n';
                        mp_m[i,j] = 1;
                        intm[j,i] = 'i';
                    end
                end

                #I:A - predation
                if r_draw > prob_line[3] && r_draw < prob_line[4]
                    rr_draw = rand();
                    if rr_draw < 0.5
                        intm[i,j] = 'e';
                        tp_m[i,j] = 1;
                        intm[j,i] = 'i';
                    else
                        intm[i,j] = 'i';
                        intm[j,i] = 'e';
                        tp_m[j,i] = 1;
                    end
                end

                #I:I - neutral interaction
                if r_draw > prob_line[4] && r_draw < prob_line[5]
                  intm[i,j] = 'i'
                  intm[j,i] = 'i'
                end

                #I:A - predation
                #N:N - symmetric mutualism
                if r_draw > prob_line[5] && r_draw < prob_line[6]
                  intm[i,j] = 'e'
                  intm[j,i] = 'e'
                  #Update mutualistic network
                  tp_m[i,j] = 1;
                  tp_m[j,i] = 1;
                end

            end #if
        end #j
    end #i

    #We could assume that any species without recorded trophic interactions is a primary producer
    total_trophic = vec(sum(tp_m,dims=2));
    prim_prod = deleteat!(findall(iszero,total_trophic),1); #eliminate row 1
    intm[prim_prod,1] .= Ref('e');
    tp_m[prim_prod,1] .= 1;

    #SPECIES-OBJECT INTERACTIONS

    # 1) pr_na
    # 2) pr_nn
    # 3) pr_ni **
    # 4) pr_nm
    # 5) pr_ia
    # 6) pr_ii **
    # 7) pr_aa

    # so_pw_prob = pw_prob_init[[3,5,6]];
    so_pw_prob = [SOpwp.ni,SOpwp.ia,SOpwp.ii];
    so_pw_prob = so_pw_prob/sum(so_pw_prob);
    so_prob_line = cumsum(so_pw_prob);

    for i=2:S
        for j=S+1:N
            if intm[i,j] == '0'
                r_draw = rand();

                #need-ignore
                if r_draw < so_prob_line[1]
                    intm[i,j] = 'n';
                    intm[j,i] = 'i';
                end

                #assimilate-ignore
                if r_draw > so_prob_line[1] && r_draw < so_prob_line[2]
                    intm[i,j] = 'e';
                    intm[j,i] = 'i';
                end

                #ignore-ignore
                if r_draw > so_prob_line[2] && r_draw < so_prob_line[3]
                    intm[i,j] = 'i';
                    intm[j,i] = 'i';
                end
            end
        end
    end

		    #Object-OBJECT INTERACTIONS

    oo_pw_prob = [OOpwp.ni,OOpwp.nn,OOpwp.ii];
    oo_pw_prob = oo_pw_prob/sum(oo_pw_prob);
    oo_prob_line = cumsum(oo_pw_prob);

    for i=(S+1):N
        for j=(S+1):N
            if intm[i,j] == '0'
                r_draw = rand();

                #need-ignore
                if r_draw < oo_prob_line[1]
                    rr_draw = rand();
										if rr_draw < 0.5
                    	intm[i,j] = 'n';
                    	intm[j,i] = 'i';
										else
											intm[i,j] = 'i';
											intm[j,i] = 'n';
										end
                end

                #need-need
                if r_draw > oo_prob_line[1] && r_draw < oo_prob_line[2]
                    intm[i,j] = 'n';
                    intm[j,i] = 'n';
                end

                #ignore-ignore
                if r_draw > oo_prob_line[2] && r_draw < oo_prob_line[3]
                    intm[i,j] = 'i';
                    intm[j,i] = 'i';
                end
            end
        end
    end

    #A matrix for Direct + Indirect trophic interactions
    tind_m = copy(tp_m);
    #A matrix for Direct + Indirect mutualistic interactions
    mind_m = copy(mp_m);

    #Force Basal primary producer (Row/Col 2) to not 'need' anything
    colonizer_n = deleteat!(findall(x->x=='n',intm[2,:]),1);

    intm[2,colonizer_n] .= Ref('i');
    mp_m[2,:] .= 0;

    #Document the indirect interactions
    for i=2:S
        for j=1:O
            if intm[spind[i],obind[j]] == 'e'
                makers = findall(x->x=='m',intm[:,obind[j]]);
                tind_m[spind[i],makers] .= 1;
            end
            if intm[spind[i],obind[j]] == 'n'
                makers = findall(x->x=='m',intm[:,obind[j]]);
                mind_m[spind[i],makers] .= 1;
            end
        end
    end


    return(intm, tp_m, tind_m, mp_m, mind_m)

end
