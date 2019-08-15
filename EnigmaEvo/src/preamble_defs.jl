function preamble_defs(intm)

    #Boolian matrices
    e_b = (intm .== 'e')*1;
    n_b = (intm .== 'n')*1;
    i_b = (intm .== 'i')*1;
    m_b = (intm .== 'm')*1;

    #copy of need binary matrix with diag = 0
    n_b0 = copy(n_b);
    n_b0[diagind(n_b0)] .= 0;

    #Vector and length of species IDs (over all intm)
    sp_v = findall(isodd,diag(n_b));
    l_sp = length(sp_v);

    int_id = collect(1:size(intm)[1]);

    return(
    e_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id
    )
end
