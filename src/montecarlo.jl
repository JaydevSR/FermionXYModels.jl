function metropolis_update!(model::FermionXYModel1D, h::Float64)
    site = rand(1:model.L)
    P_old = config_probability(model)
    model.sites[site] = -model.sites[site]  # do the flip
    P_new = config_probability(model)

    if rand() < P_new / P_old
        model.sites[site] = -model.sites[site]  # revert the flip
    end
    return model
end

function config_probabilty(model::FermionXYModel1D)
    G_corr = correlation_matrix(model)
    I_g = Diagonal(model.sites)*G_corr
    return (1/2)^model.L * det(I - I_g)
end