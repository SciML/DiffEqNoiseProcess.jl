module DiffEqNoiseProcessRandom123Ext

using DiffEqNoiseProcess
import Random123

function DiffEqNoiseProcess._resolve_vbt_rng(::DiffEqNoiseProcess._DefaultVBTRNG)
    return Random123.Threefry4x()
end

function DiffEqNoiseProcess._set_counter!(rng::Random123.AbstractR123, seed)
    return Random123.set_counter!(rng, seed)
end

end
