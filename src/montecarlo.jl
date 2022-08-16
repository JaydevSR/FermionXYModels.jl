export metropolis_update!, config_probability

@inbounds function metropolis_update!(model::FermionXYChain; P_old::Float64=config_probability(model))
    site = rand(1:model.L)
    model.sites[site] *= -1  # do the flip
    P_new = config_probability(model)
    if rand() > (P_new / P_old)
        model.sites[site] *= -1  # revert the flip
    end
    return model, P_new
end

function equilibrate!(model::FermionXYChain, steps::Int=1000)
    p = config_probability(model)
    for _ in 1:steps
        model, p = metropolis_update!(model, P_old=p)
    end
    return model
end

config_probability(model::FermionXYChain) = det(probability_matrix(model))
