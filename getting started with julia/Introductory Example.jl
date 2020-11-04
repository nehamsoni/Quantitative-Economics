### A Pluto.jl notebook ###
# v0.12.6

using Markdown
using InteractiveUtils

# ╔═╡ c9f9318a-1e89-11eb-2939-f3a5f1289a12
begin
	import Pkg; Pkg.add("ForwardDiff")
	using ForwardDiff
end

# ╔═╡ 652e2558-1e7a-11eb-2483-cb1203507751
begin
	using LinearAlgebra
	using InstantiateFromURL
	using Plots
	gr(fmt=:png);
end

# ╔═╡ e30d5e5c-1e7b-11eb-0173-5dbf06d875bb
md"""
# Introductory Examples
### Solutions to the Exercise"""

# ╔═╡ 02f30026-1e7c-11eb-1bdf-713cc1950a22
md"""
##### `Exercise 1 - 4 are quite naive and will only help you to get familiarity with the platform` """

# ╔═╡ c0b1a1c0-1e8e-11eb-11e4-9187b5708dad
html"""<br><br>"""

# ╔═╡ 0297fd2c-1e7c-11eb-20de-ed44583cc810
md"""
##### Exercise 5
Simulate and plot the correlated time series

xt+1 = αxt + ϵt+1 wherex0 =0 and t= 0,…,n.  
The sequence of shocks {ϵt} is assumed to be iid and standard normal.
Set n=200 and α=0.9."""

# ╔═╡ c5a9b326-1e7c-11eb-043e-6b700836e7df
begin
	path=accumulate(1:200;init=0) do old,i
		old*0.9 +randn()
	end
	plot(path,legend=false)
end

# ╔═╡ 79411808-1e7e-11eb-1073-d9d8e2b1a662
html"""<hr><br><br>"""

# ╔═╡ 55e7aa80-1e7d-11eb-3447-810cb162aa9b
md"""
##### Exercise 6
Plot three simulated time series, one for each of the cases α=0, α=0.8 and α=0.98.

(The figure will illustrate how time series with the same one-step-ahead conditional volatilities, as these three processes have, can have very different unconditional volatilities)"""

# ╔═╡ 96c557f0-1e7d-11eb-0577-1bef484e817a
let p=plot()
	for α in [0,0.8,0.98]
		pa=accumulate(1:200;init=1) do old,i
			old*α +randn()
		end
		plot!(p,pa,label="α = $α")
	end
	p
end	

# ╔═╡ 9438bd3c-1e7e-11eb-038d-c759ea6cf550
html"""<hr><br><br>"""

# ╔═╡ b83ea3d6-1e7e-11eb-3c70-f5e2349bd425
md"""
##### Exercise 7
*This exercise is more challenging.*

Take a random walk, starting from x0=1
xt+1 = αxt + σϵt+1 where x0=1 and t= 0,…,tmax

* Furthermore, assume that the xtmax=0 (i.e. at tmax, the value drops to zero, regardless of its current state).
* The sequence of shocks {ϵt} is assumed to be iid and standard normal.

* For a given path {xt} define a first-passage time as Ta=min{t|xt≤a}, where by the assumption of the process Ta≤tmax.
Start with σ=0.2,α=1.0. 
calculate the first-passage time, T0, for 100 simulated random walks – to a tmax=200 and plot a histogram
plot the sample mean of T0 from the simulation for α∈{0.8,1.0,1.2}"""

# ╔═╡ ae782114-1e7f-11eb-31bb-9f1b009f6853
let left= plot(), right=plot(), σ=0.2
	bigg=[]
	for j in [0.8,1.0,1.2]
		result_=Vector{}()
		for i ∈ 1:100
			pathing=accumulate(1:200;init=1) do old, i
				j*old + σ*randn()
			end
			push!(result_,findfirst(x->x<0.0,pathing))
		end
		push!(bigg,replace(result_,nothing=>200))
	end
	left=bar(bigg[2],bins=:scott,legend=false)
	right=bar(sum.(bigg)./100,bins=:scott,label=[0.8,1.0,1.2])
	plot(left,right)
end

# ╔═╡ a56d8f26-1e89-11eb-2c21-791307ccc9c2
html"""<hr><br><br>"""

# ╔═╡ 6167139e-1e89-11eb-26f4-55cbdbc70de2
md"""
##### Exercise 8
*This exercise is more challenging.*

The root of a univariate function f(⋅) is an x such that f(x)=0.
One solution method to find local roots of smooth functions is called Newton’s method.
Starting with an x0 guess, a function f(⋅) and the first-derivative f′(⋅), the algorithm is to repeat

xn+1=xn−f(xn)f′(xn)

until |xn+1−xn| is below a tolerance

Use a variation of the fixedpointmap code to implement Newton’s method, where the function would accept arguments f, f_prime, x_0, tolerance, maxiter.
Test it with f(x)=(x−1)3 and another function of your choice where you can analytically find the derivative."""

# ╔═╡ 3928f3c2-1e8c-11eb-2244-3f9c7c2857d4
begin
	f(x) = (x-1)^3
	f_d(x) = ForwardDiff.derivative(f, x)
	
	f(1.7),f_d(1.7)
end

# ╔═╡ 4b431106-1e8b-11eb-00cd-814c8e98e0b7
let
	maxitr=200
	ini=rand(1:10)
	second=0
	for i in 1:maxitr
		second= ini - (f(ini)/f_d(ini))
		if abs(second-ini)< 1e-7
			break
		else
			ini=second
		end
	end
	second
end

# ╔═╡ Cell order:
# ╟─e30d5e5c-1e7b-11eb-0173-5dbf06d875bb
# ╟─02f30026-1e7c-11eb-1bdf-713cc1950a22
# ╟─c0b1a1c0-1e8e-11eb-11e4-9187b5708dad
# ╟─0297fd2c-1e7c-11eb-20de-ed44583cc810
# ╠═c5a9b326-1e7c-11eb-043e-6b700836e7df
# ╟─79411808-1e7e-11eb-1073-d9d8e2b1a662
# ╟─55e7aa80-1e7d-11eb-3447-810cb162aa9b
# ╠═96c557f0-1e7d-11eb-0577-1bef484e817a
# ╟─9438bd3c-1e7e-11eb-038d-c759ea6cf550
# ╟─b83ea3d6-1e7e-11eb-3c70-f5e2349bd425
# ╠═ae782114-1e7f-11eb-31bb-9f1b009f6853
# ╟─a56d8f26-1e89-11eb-2c21-791307ccc9c2
# ╟─6167139e-1e89-11eb-26f4-55cbdbc70de2
# ╠═3928f3c2-1e8c-11eb-2244-3f9c7c2857d4
# ╠═4b431106-1e8b-11eb-00cd-814c8e98e0b7
# ╟─c9f9318a-1e89-11eb-2939-f3a5f1289a12
# ╟─652e2558-1e7a-11eb-2483-cb1203507751
