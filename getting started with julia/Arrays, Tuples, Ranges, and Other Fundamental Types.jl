### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ 1c7e00f8-330d-11eb-36a3-fb25098f61f0
begin
	import Pkg; Pkg.add("RollingFunctions")
	using LinearAlgebra
	using Plots
	using QuantEcon
	using RollingFunctions
end

# ╔═╡ 9bb7e430-3309-11eb-045b-45984b90b650
md"""
# Arrays, Tuples, Ranges, and Other Fundamental Types
### Solutions to the Exercise
"""

# ╔═╡ 6b36bba0-330a-11eb-2f33-6d4b9f648302

html"""<a href="https://julia.quantecon.org/getting_started_julia/fundamental_types.html#Exercise-1">Exercise 1 (click)</a>"""



# ╔═╡ c9846d96-330b-11eb-05ec-7774cfa3e09e
function ComputedUnconditionVariables(A::AbstractArray,Σ::AbstractArray;max_err =1e-7,maxitr=1000)
	
	S=zeros(size(A))
	i=0
	while true
		i+=1
		newS = A*S*A' + Σ*Σ'
		if norm(newS-S) < max_err && i>maxitr
			break
		end
		S = newS
	end
	S
end

# ╔═╡ 60fe778c-330e-11eb-3efd-37fe60e5becf
let 
	A =[0.8 -0.2;-0.1 0.7]
	Σ = [0.5 0.4;0.4 0.6]
	norm(ComputedUnconditionVariables(A,Σ)-solve_discrete_lyapunov(A,Σ*Σ'))
end

# ╔═╡ 3ea21b48-330f-11eb-1fab-c396a34b1089
html"""
<a href="https://julia.quantecon.org/getting_started_julia/fundamental_types.html#Exercise-2">Exercise 2 (click)</a>"""

# ╔═╡ 382c3762-3312-11eb-3371-f71f6337cd13
Θ =[0.8,0.9,0.98]

# ╔═╡ 75bf8eb2-330f-11eb-2313-fb73eb133cbb
function CalculateY(T::Int64,θ::Float64;σ=1,γ=1,y=0)
	A=Array{Float64}(undef,0)
	for i in 1:T
		y = γ + θ*y +σ*randn()
		push!(A,y)
	end
	A
end

# ╔═╡ c8e9e352-3310-11eb-35dd-f72acb3b5373
function CalculateY(T::Int64;σ=1,γ=1,y=0)
	ans=[]
	for θ in Θ
		push!(ans,CalculateY(T,θ;σ=1,γ=1,y=0))
	end
	ans
end

# ╔═╡ f53d440c-3311-11eb-217a-2f5e52af604b
let
	res = CalculateY(150)
	p=plot(1:150,res)
	rollingres = rolling.(mean,res,10)
	plot!(p,1:150,cat.(Ref(zeros(9)),rollingres,dims=1),legend=false)
end
	

# ╔═╡ 75978c48-3317-11eb-09b8-d33d9c09cb4f
let
	answers = Dict()
	for θ in Θ
		answers[θ]=Array{Float64}(undef,0)
	end
	for N in 1:200
		p=CalculateY(150)
		for (i,v) in enumerate(Θ)
			push!(answers[v],p[i][end])
		end
	end
	all=cat(collect(values(answers))...,dims=1)
	t1=mean(all)# answer for last part doing it for values of Θ collectively
	t2=mean(x->(x^2),all)# answer of last part
	@show t1
	@show t2-t1^2
	histogram(answers;alpha=0.5,label = keys(answers))
end

# ╔═╡ 1365fe9c-3321-11eb-1b03-1bcb8ec30dc8
html"""
<a 
href="https://julia.quantecon.org/getting_started_julia/fundamental_types.html#Exercise-3">Exercise 3 (click)</a>"""

# ╔═╡ 2df0c788-3321-11eb-08c3-cf636bb7d128
begin
	parameters=[0.1,0.2,0.5,1.0,1.0]
	Y=[]
	X=[]
	for M in 1:20
		x1=randn(50)
		x2=randn(50)
		w=randn(50)
		push!(X,cat(x1,x1.^2,x2,ones(50),w,dims=2))
		push!(Y,X[end]*parameters)
	end
	X,Y
end

# ╔═╡ d44bb630-3324-11eb-1dd8-47e7ce22a550
function find_parameters(x::AbstractArray,y::AbstractArray)
	(x'*x)\(x'*y)
end

# ╔═╡ 5f790744-3325-11eb-2629-f52ace14e7bb
para=[ find_parameters(X[M],Y[M]) for M in 1:20]

# ╔═╡ 108c296c-3326-11eb-04a8-3bdb958ddeef
let
	para2=Matrix(reshape(cat(para...,dims=1),(5,20))')# for display
	histogram(para2,bins=0:0.1:1.1;alpha=.2,label=['a','b','c','d','σ'])
end

# ╔═╡ 518146e4-330a-11eb-1720-d35a5e832af2
bigbreak=html"""<br><br><br>"""

# ╔═╡ 697a48ea-330f-11eb-1650-1fb81c301c71
bigbreak

# ╔═╡ 2fe0b666-3321-11eb-0598-eb28cd43af60
bigbreak

# ╔═╡ Cell order:
# ╟─9bb7e430-3309-11eb-045b-45984b90b650
# ╟─6b36bba0-330a-11eb-2f33-6d4b9f648302
# ╠═c9846d96-330b-11eb-05ec-7774cfa3e09e
# ╠═60fe778c-330e-11eb-3efd-37fe60e5becf
# ╟─697a48ea-330f-11eb-1650-1fb81c301c71
# ╟─3ea21b48-330f-11eb-1fab-c396a34b1089
# ╠═382c3762-3312-11eb-3371-f71f6337cd13
# ╠═75bf8eb2-330f-11eb-2313-fb73eb133cbb
# ╠═c8e9e352-3310-11eb-35dd-f72acb3b5373
# ╟─f53d440c-3311-11eb-217a-2f5e52af604b
# ╠═75978c48-3317-11eb-09b8-d33d9c09cb4f
# ╟─2fe0b666-3321-11eb-0598-eb28cd43af60
# ╟─1365fe9c-3321-11eb-1b03-1bcb8ec30dc8
# ╠═2df0c788-3321-11eb-08c3-cf636bb7d128
# ╠═d44bb630-3324-11eb-1dd8-47e7ce22a550
# ╠═5f790744-3325-11eb-2629-f52ace14e7bb
# ╠═108c296c-3326-11eb-04a8-3bdb958ddeef
# ╟─1c7e00f8-330d-11eb-36a3-fb25098f61f0
# ╟─518146e4-330a-11eb-1720-d35a5e832af2
