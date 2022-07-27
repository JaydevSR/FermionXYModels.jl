function metropolis_update!(model::FermionXYModel1D; P_old::Float64=config_probability(model))
    site = rand(1:model.L)
    model.sites[site] = -model.sites[site]  # do the flip
    P_new = config_probability(model)
    if rand() > (P_new / P_old)
        model.sites[site] = -model.sites[site]  # revert the flip
    end
    return model, P_new
end

function equilibrate!(model::FermionXYModel1D, steps::Int=1000)
    p = config_probability(model)
    for _ in 1:steps
        model, p = metropolis_update!(model, P_old=p)
    end
    return model
end

function config_probability(model::FermionXYModel1D)
    return det(probability_matrix(model))
end