module DiffEqNoiseProcessOptimExt

using DiffEqNoiseProcess
import DiffEqNoiseProcess: constrained_optimization_problem, linear_interpolation_wedges
import Optim

function constrained_optimization_problem(densf, fij, fij2, fij3, fij4, ri, ai, Δr, Δa)
    function difference(x)
        return densf(x[1], x[2]) -
            linear_interpolation_wedges(
            fij, fij2, fij3, fij4, x[1], x[2], ri, ai, Δr,
            Δa
        )
    end
    ϵijmax = Optim.optimize(
        x -> -difference(x), [ri, ai], [ri + Δr, ai + Δa],
        [ri + Δr / 2, ai + Δa / 2], Optim.Fminbox(Optim.NelderMead())
    )
    ϵijmin = Optim.optimize(
        difference, [ri, ai], [ri + Δr, ai + Δa],
        [ri + Δr / 2, ai + Δa / 2], Optim.Fminbox(Optim.NelderMead())
    )
    return Optim.minimum(ϵijmin), Optim.minimum(ϵijmax)
end

end
