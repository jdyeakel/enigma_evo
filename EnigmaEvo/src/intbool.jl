function intbool(intm)

    #Boolian matrices
    e_b = (intm .== 1);
    n_b = (intm .== 2);
    # i_b = (intm .== 'i')*1;
    # m_b = (intm .== 'm')*1;

    #copy of need binary matrix with diag = 0
    n_b0 = copy(n_b);
    n_b0[diagind(n_b0)] .= false;

    # #Vector and length of species IDs (over all intm)
    # sp_v = findall(isodd,diag(n_b));
    # l_sp = length(sp_v);

    # int_id = collect(1:size(intm)[1]);

    return(
    e_b,
    n_b,
    n_b0
    )
end
