### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ be75fc5e-3247-11eb-1c9f-6f30ea72331c
begin
	import Pkg; Pkg.add("ForwardDiff")
	Pkg.add("Plots")
	using Plots
end

# ╔═╡ 5f085b68-323f-11eb-3ee9-13b634011eac
md"""
# Julia Essentials
### Solutions to the Exercise
"""

# ╔═╡ e6599ddc-323f-11eb-3f27-a7f05eb26ae5
md"""
##### `Exercise 1 - 4 are quite naive and will only help you to get familiarity with the platform, you are advised to permorn on your own` """

# ╔═╡ e57c7f4c-323f-11eb-26dd-5b88be11bed5
html"""<br><br>"""

# ╔═╡ c994bf8e-3245-11eb-259f-a791ae0734dc
md"""
##### Exercise 5
The Julia libraries include functions for interpolation and approximation.

Nevertheless, let’s write our own function approximation routine as an exercise.

In particular, write a function linapprox that takes as arguments

A function f mapping some interval [a,b] into ℝ.
two scalars a and b providing the limits of this interval.
An integer n determining the number of grid points.
A number x satisfying a ≤ x ≤ b.
and returns the piecewise linear interpolation of f at x, based on n evenly spaced grid points a = point[1] < point[2] < ... < point[n] = b.

Aim for clarity, not efficiency.

Hint: use the function range to linearly space numbers."""

# ╔═╡ 20cd089a-3246-11eb-368f-fdbc5129f16c
function linapprox(f,a,b,n,x)
	Δ = (b-a)/(n-1)
	k = (x-a)÷Δ
	u,v = k*Δ +a,(k+1)*Δ +a
	δ= (f(v)-f(u))/(Δ)
	f(u) + (x-u)*δ
end
	

# ╔═╡ 4a812cde-3249-11eb-3d46-c1937d75ef1e
html"""<br><br>"""

# ╔═╡ 42fd5d3c-3248-11eb-0dac-43b98063620e
md"""
##### Exercise 6 
`should be easy, try to do it yourself`
"""

# ╔═╡ 58b6b208-3249-11eb-110e-1f2d5f82756c
html"""<br><br>"""

# ╔═╡ 46a8a40c-3249-11eb-345a-8143e3ac48cc
md"""
##### Exercise 7
Redo Exercise 5 except

Pass in a range instead of the a, b, and n. Test with a range such as nodes = -1.0:0.5:1.0.
Instead of the while used in the solution to Exercise 5, find a better way to efficiently bracket the x in the nodes.
Hints:

Rather than the signature as function linapprox(f, a, b, n, x), it should be called as function linapprox(f, nodes, x).
step(nodes), length(nodes), nodes[1], and nodes[end] may be useful.
Type ?÷ into jupyter to explore quotients from Euclidean division for more efficient bracketing."""

# ╔═╡ 78f845ea-3249-11eb-3ac1-fb51649acd04
function linapprox(f,nodes,x)
	a,b = nodes[1],nodes[end]
	n = length(nodes)
	linapprox(f,a,b,n,x)
end

# ╔═╡ 9dbef2fe-3247-11eb-2827-1b0878e2cd0c
let
	f_ex5(x) = x^2
	g_ex5(x) = linapprox(f_ex5, -1, 1, 3, x)
	x_grid = range(-1.0, 1.0, length = 100)
	y_vals = f_ex5.(x_grid)
	y = g_ex5.(x_grid)
	plot(x_grid, y_vals, label = "true")
	plot!(x_grid, y, label = "approximation")
end

# ╔═╡ 9e334f98-324a-11eb-024f-93e55e491e58
let
	f_ex5(x) = x^2
	g_ex5(x) = linapprox(f_ex5, -1:5, x)
	x_grid = range(-1.0, 5.0, length = 100)
	y_vals = f_ex5.(x_grid)
	y = g_ex5.(x_grid)
	plot(x_grid, y_vals, label = "true")
	plot!(x_grid, y, label = "approximation")
end

# ╔═╡ Cell order:
# ╟─5f085b68-323f-11eb-3ee9-13b634011eac
# ╟─e6599ddc-323f-11eb-3f27-a7f05eb26ae5
# ╟─e57c7f4c-323f-11eb-26dd-5b88be11bed5
# ╟─c994bf8e-3245-11eb-259f-a791ae0734dc
# ╠═20cd089a-3246-11eb-368f-fdbc5129f16c
# ╟─9dbef2fe-3247-11eb-2827-1b0878e2cd0c
# ╟─4a812cde-3249-11eb-3d46-c1937d75ef1e
# ╟─42fd5d3c-3248-11eb-0dac-43b98063620e
# ╟─58b6b208-3249-11eb-110e-1f2d5f82756c
# ╟─46a8a40c-3249-11eb-345a-8143e3ac48cc
# ╠═78f845ea-3249-11eb-3ac1-fb51649acd04
# ╟─9e334f98-324a-11eb-024f-93e55e491e58
# ╟─be75fc5e-3247-11eb-1c9f-6f30ea72331c
