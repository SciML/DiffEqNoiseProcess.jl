module DiffEqNoiseProcessChainRulesCoreExt

using SciMLBase: AbstractNoiseProcess
using ChainRulesCore: ChainRulesCore, NoTangent, ProjectTo

ChainRulesCore.ProjectTo(::AbstractNoiseProcess) = ProjectTo{NoTangent}()

end
