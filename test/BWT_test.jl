#@testset "BoxWedgeTail tests" begin

using DiffEqNoiseProcess, Test, Random
import Distributions
using Cubature
using Statistics

# Test generation of boxes
W = BoxWedgeTail(0.0,zeros(2), box_grouping = :Columns)
# number of boxes with column-wise packing
@test length(W.boxes.boxes) == 2894
@test sum(W.boxes.probability) ≈ 0.911857 rtol=1e-5
@test DiffEqNoiseProcess.entropy(W.boxes.probability) ≈ 10.13 rtol=1e-5


W = BoxWedgeTail(0.0,zeros(2), box_grouping = :none)
# number of boxes of smallest size
@test length(W.boxes.boxes) == 119519
@test minimum(2*W.Δr*W.Δa*W.Δz .== W.boxes.probability)
@test sum(W.boxes.probability) ≈ 0.911857 rtol=1e-5
@test DiffEqNoiseProcess.entropy(W.boxes.probability) ≈ Distributions.entropy(W.boxes.dist,2) rtol=1e-5

W = BoxWedgeTail(0.0,zeros(2), box_grouping = :MinEntropy)
# number of boxes with min-entropy packing
@test length(W.boxes.boxes) == 2975
@test sum(W.boxes.probability) ≈ 0.911857 rtol=1e-5
@test DiffEqNoiseProcess.entropy(W.boxes.probability) ≈ 7.14 rtol=1e-3

# test generation of wedges

(val,err) = hcubature(x-> W.jpdf(x[1],x[2]), [zero(W.rM),zero(W.aM)],  [W.rM,W.aM])
@test 2*val ≈ 0.99732 rtol=1e-5

# number of columns
@test (W.aM/W.Δa) * (W.rM/W.Δr) == 4096
# probability for wedges
@test 2*val - sum(W.boxes.probability) ≈ 0.0855 rtol=1e-3
# on average less than half of the samples have to be rejected
@test sum(W.wedges.probability) <= 2*(2*val - sum(W.boxes.probability))
@test length(W.wedges.boxes) == 4096

# without sqeezing method
W = BoxWedgeTail(0.0,zeros(2), box_grouping = :MinEntropy, sqeezing=false)
# probability for wedges
@test 2*val - sum(W.boxes.probability) ≈ 0.0855 rtol=1e-3
# on average less than half of the samples have to be rejected
@test sum(W.wedges.probability) <= 2*(2*val - sum(W.boxes.probability))
@test length(W.wedges.boxes) == 4096

# check probabilities for tail approximation
# region 1
val1, _ = hcubature(x->W.jpdf(x[1],x[2]), [W.rM,zero(W.aM)], [12*one(W.rM),W.aM])
@test 2*val1 ≈ 2.98*1e-4 rtol=1e-3
@test W.tails.tail1.p ≈ 2*val1 rtol=1e-10
# region 2
val2, _ = hcubature(x->W.jpdf(x[1],x[2]), [W.rM,W.aM], [8*one(W.rM),8*one(W.aM)])
@test 2*val2 ≈ 3.65*1e-5 rtol=1e-3
@test W.tails.tail2.p ≈ 2*val2 rtol=1e-10
# region 3
val3, _ = hcubature(x->W.jpdf(x[1],x[2]), [2*one(W.rM),W.aM], [W.rM,8*one(W.aM)])
@test 2*val3 ≈ 1.80*1e-3 rtol=1e-2
@test W.tails.tail3.p ≈ 2*val3 rtol=1e-10
# region 4
val4, _ = hcubature(x->W.jpdf(x[1],x[2]), [zero(W.rM),W.aM], [one(W.rM)/2,6*one(W.aM)])
@test 2*val4 ≈ 1.61*1e-6 rtol=1e-2
@test W.tails.tail4.p ≈ 2*val4 rtol=1e-10
# region 5
val5, _ = hcubature(x->W.jpdf(x[1],x[2]), [one(W.rM)/2,W.aM], [one(W.rM),6*one(W.aM)])
@test 2*val5 ≈ 1.96*1e-5 rtol=1e-2
@test W.tails.tail5.p ≈ 2*val5 rtol=1e-10
# region 6
val6, _ = hcubature(x->W.jpdf(x[1],x[2]), [one(W.rM),W.aM], [3*one(W.rM)/2,6*one(W.aM)])
@test 2*val6 ≈ 1.20*1e-4 rtol=1e-2
@test W.tails.tail6.p ≈ 2*val6 rtol=1e-10
# region 7
val7, _ = hcubature(x->W.jpdf(x[1],x[2]), [3*one(W.rM)/2,W.aM], [2*one(W.rM),6*one(W.aM)])
@test 2*val7 ≈ 3.86*1e-4 rtol=1e-3
@test W.tails.tail7.p ≈ 2*val7 rtol=1e-10
# region 8
val8, _ = hcubature(x->W.jpdf(x[1],x[2]), [one(W.rM),6*one(W.aM)], [2*one(W.rM),8*one(W.aM)])
@test 2*val8 ≈ 6.57*1e-6 rtol=1e-3
@test W.tails.tail8.p ≈ 2*val8 rtol=1e-10
# region 9
val9, _ = hcubature(x->W.jpdf(x[1],x[2]), [2*one(W.rM),8*one(W.aM)], [5*one(W.rM),10*one(W.aM)])
@test 2*val9 ≈ 4.15*1e-6 rtol=1e-3
@test W.tails.tail9.p ≈ 2*val9 rtol=1e-10

# test remainder
@test 1.0-(2*(val+val1+val2+val3+val4+val5+val6+val7+val8+val9)) ≈ 3.81e-7 atol=1e-3

# tests if sampling works (samples from expected domain)

# boxes
samples = [DiffEqNoiseProcess.sample_box(W, W.boxes) for i=1:100_000]
@test minimum(getindex.(samples, 1)) > 0
@test minimum(getindex.(samples, 2)) > 0
@test maximum(getindex.(samples, 1)) < W.rM
@test maximum(getindex.(samples, 2)) < W.aM

# wedges
samples = [DiffEqNoiseProcess.sample_wedge(W, W.wedges) for i=1:100_000]
@test minimum(getindex.(samples, 1)) > 0
@test minimum(getindex.(samples, 2)) > 0
@test maximum(getindex.(samples, 1)) < W.rM
@test maximum(getindex.(samples, 2)) < W.aM

# tails
samples = [DiffEqNoiseProcess.sample_tail(W.rng, W.jpdf, W.tails.tail1) for i=1:100_000]
@test minimum(getindex.(samples, 1)) > 4
@test minimum(getindex.(samples, 2)) > 0
@test maximum(getindex.(samples, 1)) < 12
@test maximum(getindex.(samples, 2)) < 4

samples = [DiffEqNoiseProcess.sample_tail(W.rng, W.jpdf, W.tails.tail2) for i=1:100_000]
@test minimum(getindex.(samples, 1)) > 4
@test minimum(getindex.(samples, 2)) > 4
@test maximum(getindex.(samples, 1)) < 8
@test maximum(getindex.(samples, 2)) < 8

samples = [DiffEqNoiseProcess.sample_tail(W.rng, W.jpdf, W.tails.tail3) for i=1:100_000]
@test minimum(getindex.(samples, 1)) > 2
@test minimum(getindex.(samples, 2)) > 4
@test maximum(getindex.(samples, 1)) < 4
@test maximum(getindex.(samples, 2)) < 8

samples = [DiffEqNoiseProcess.sample_tail(W.rng, W.jpdf, W.tails.tail4) for i=1:100_000]
@test minimum(getindex.(samples, 1)) > 0
@test minimum(getindex.(samples, 2)) > 4
@test maximum(getindex.(samples, 1)) < 0.5
@test maximum(getindex.(samples, 2)) < 6

samples = [DiffEqNoiseProcess.sample_tail(W.rng, W.jpdf, W.tails.tail5) for i=1:100_000]
@test minimum(getindex.(samples, 1)) > 0.5
@test minimum(getindex.(samples, 2)) > 4
@test maximum(getindex.(samples, 1)) < 1
@test maximum(getindex.(samples, 2)) < 6

samples = [DiffEqNoiseProcess.sample_tail(W.rng, W.jpdf, W.tails.tail6) for i=1:100_000]
@test minimum(getindex.(samples, 1)) > 1
@test minimum(getindex.(samples, 2)) > 4
@test maximum(getindex.(samples, 1)) < 1.5
@test maximum(getindex.(samples, 2)) < 6

samples = [DiffEqNoiseProcess.sample_tail(W.rng, W.jpdf, W.tails.tail7) for i=1:100_000]
@test minimum(getindex.(samples, 1)) > 1.5
@test minimum(getindex.(samples, 2)) > 4
@test maximum(getindex.(samples, 1)) < 2
@test maximum(getindex.(samples, 2)) < 6

samples = [DiffEqNoiseProcess.sample_tail(W.rng, W.jpdf, W.tails.tail8) for i=1:100_000]
@test minimum(getindex.(samples, 1)) > 1
@test minimum(getindex.(samples, 2)) > 6
@test maximum(getindex.(samples, 1)) < 2
@test maximum(getindex.(samples, 2)) < 8

samples = [DiffEqNoiseProcess.sample_tail(W.rng, W.jpdf, W.tails.tail9) for i=1:100_000]
@test minimum(getindex.(samples, 1)) > 2
@test minimum(getindex.(samples, 2)) > 8
@test maximum(getindex.(samples, 1)) < 5
@test maximum(getindex.(samples, 2)) < 10


#end
