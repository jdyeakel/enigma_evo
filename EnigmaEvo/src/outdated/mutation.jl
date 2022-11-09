function mutation(spmut,intmut,spints,obints,intm,spv,tallytable,evolutiontable)

    oldint_in = intm[spmut,intmut];
    #how intmut interacts with spmut
    oldint_out = intm[intmut,spmut];

    intm_mut = copy(intm);
    
    if in(intmut,spv)
        newint_in = rand(setdiff(spints,oldint_in))[1];
        newint_out = rand(setdiff(spints,oldint_out))[1];
        # Select evolution of either in degree interaction or out-degree interaction
        evol_type_draw = rand();
        if evol_type_draw < 0.5
            newint = newint_in;
            oldint = oldint_in;
            intm_mut[spmut,intmut] = newint;
            evol_type = 1.;
        else
            newint = newint_out;
            oldint = oldint_out;
            intm_mut[intmut,spmut] = newint;
            evol_type = 2.;
        end
    #If the mutated interaction is with an object
    else
        #For objects, randomly choose ignore, eat, need, make
        newint_in = rand(setdiff(obints,oldint_in))[1];
        newint = newint_in;
        oldint = oldint_in;
        intm_mut[spmut,intmut] = newint;
        evol_type = 1.;

        #RECIPROCAL OBJECT NEEDS
        #If an object is newly made, it now needs the species
        if newint == 3;
            intm_mut[intmut,spmut] = 2;
        end
        #If an object is newly unmade, it now does not need the species
        if oldint == 3;
            intm_mut[intmut,spmut] = 0;
        end
    end

    # ACCEPT REGARDLESS
    ebmut,nbmut,nb0mut,mbmut = intbool(intm_mut);
    tally = tallytable[(evolutiontable[1,:] .== oldint) .& (evolutiontable[2,:] .== newint)][1];

    #Record evolution type (in degree vs. out degree) on backend of tally
    tally += evol_type*0.01;

    return(intm_mut, ebmut, nbmut, nb0mut, mbmut, tally)

end
