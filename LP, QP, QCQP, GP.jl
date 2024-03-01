#
# Written by: Milad Dehghani Filabadi

# In this code, I have written the codes of 4 problems nameley:
#   Example 1: linear programming (LP)
#   Example 2: quadratic programming (QP)
#   Example 3: quadratically constrained quadratic programming (QCQP)
#   Example 4: geometric programming (GP)
# using JuMP and MoosekTools.
# The original formulations as well as conic formulations are available in this file.

using JuMP, MosekTools
using LinearAlgebra    # for decomposing Matrices

# Example 1
#        LP formulation                 Conic formulation
#        max c'x                        max c'x
#        s.t: Ax <= b (-Ax+b >= 0)      s.t: -Ax + b in K
#             x >= 0                    K is Nonnegatives (orthant)
#
function LP(solver, A :: Array{Float64,2}, c :: Vector{Float64}, b :: Vector{Float64})
    model = Model(optimizer_with_attributes(solver))

    n = length(c)
    @variable(model, x[1:n] >= 0)

    # Define the objective function
    @objective(model, Max, sum(c .* x))   # . in julia is used for element-wise operation

    # Define the constraints
    @constraint(model, A * x .<= b)

    # Solve the LP problem
    optimize!(model)

    # return the following expressions
    termination_status(model), objective_value(model), [value(item) for item in x]

end

# Example 2
# In the following, a quadratic programming (QP) is solved where Q is
# positive semi-definite and can be written as Q = F'F.
#
#        QP formulation                 Conic formulation
#        max c'x                        max c'x
#   (1)  s.t: (1/2)xQx + a'x <= b       s.t: t + a'x - b = 0
#             x >= 0 and Integer          (t,1,Fx) in 𝒬r^(m+2)
#                                         𝒬r^(m+2) is RotatedSecondOrderCone
#                                         x Int
#
#  constraint (1) is: 5x1^2 -2x1*x2*+4x2^2 - 4x1 +3x2 <= 2
#  =>  Q = [5   -1;  -1   4],  a' = [4  3], b = 2

function QP(solver, F :: Array{Float64,2},c :: Vector{Float64},
        a :: Vector{Float64}, b :: Float64)
    model = Model(optimizer_with_attributes(solver))

    n = length(c)  # or m = size(F,1)

    @variable(model, x[1:n] >= 0, Int)
    @variable(model, zQP)
    @variable(model, t)

    @constraint(model, t + a'*x-b == 0)
    @constraint(model, [t;1;F*x] in MOI.RotatedSecondOrderCone(n+2))
    @constraint(model, zQP == c'*x)

    @objective(model, Max, zQP)

    optimize!(model)

    termination_status(model), objective_value(model), [value(item) for item in x]

end


# Example 3
# In the following, a quadratically constrained quadratic programming (QCQP)
# is solved where Q and P are positive semi-definite and can be written
# as Q = F'F, P = V'V
#
#      QCQP formulation                Conic formulation
#      max (1/2)xPx + c'x              max s + c'x
#  (1) s.t: (1/2)xQx + a'x <= b        s.t: t + a'x - b = 0
#           x >= 0                      (t,1,Fx) in 𝒬r^(m+2)
#                                       (s,1,Vx) in 𝒬r^(m+2)
#                                       𝒬r^(m+2) is RotatedSecondOrderCone
#  Assume:
# constraint (1) is: 5x1^2 -2x1*x2*+4x2^2 - 4x1 +3x2 <= 2
#  =>  Q = [5   -1;  -1   4],  a' = [4  3], b = 2
# objective function is: 3x1^2 +2x1*x2*+2x2^2 + x1 +x2
#  =>  P = [3   1;  1   2],  c' = [1  1]

function QCQP(solver, F :: Array{Float64,2}, V :: Array{Float64,2},
    c :: Vector{Float64}, a :: Vector{Float64}, b :: Float64 )

    model = Model(optimizer_with_attributes(solver))

    n = length(c)  # or m = size(F,1)

    @variable(model, x[1:n] >= 0)
    @variable(model, zQCQP)
    @variable(model, t)
    @variable(model, 0 <= s <= 3)

    @constraint(model, t + a'*x-b == 0)
    @constraint(model, [t;1;F*x] in MOI.RotatedSecondOrderCone(n+2))
    @constraint(model, zQCQP == s + c'*x)
    @constraint(model, [s;1;V*x] in MOI.RotatedSecondOrderCone(n+2))

    @objective(model, Max, zQCQP)

    optimize!(model)

    termination_status(model), objective_value(model), [value(item) for item in x]

end

# Example 4
# In the following, I solve a Geometric Programming (GP) using exponential cones.
# The problem we deal with is as follows:
#
#   GP formulation
#  max   x + z*y^(2)                      (o)
# s.t: 0.1*x^(0.5) + 2*y^(-1) <= 1        (1)
#       z^(-1) + y*x^(-2)  <= 1           (2)
#
# To write the conic reformulation, define variables x = e^u, y = e^v, z = e^w.
# Simplyfing the formulation and using the definition of Exponential Cone as
# x1 >= x2 * e^(x3/x2)   <===> (x1, x2, x3) in ExponentialCone(3),
# we get the following conic reformulation.
#
# NOTE1: See section 5.3.1 of https://docs.mosek.com/modeling-cookbook/expo.html
# to see the details of this reformulation
#
#      Conic reformulation
#      Min  t
# (o):     p1 + q1 <= 1  or (  [-p1-q1+1;1] in Nonnegatives(2) )
#          (p1, 1, u-t) & (q1, 1, 2v+w−t) in ExponentialCone(3)
# (1):     p2 + q2 <= 1
#          (p2, 1, 0.5*u+log0.1) & (q2, 1, −v+log2) in ExponentialCone(3)
# (2):     p3 + q3 <= 1
#          (p3, 1, −w) & (q3, 1, v−2u) in ExponentialCone(3)
# the optimal solution for GP should be: [ x,y,z ] = [ 3.14 ,2.43,1.32]
#
#################################################################
###############  IMPORTANT NOTES FOR IMPLEMENTATION #############
#################################################################
# Baed on most of conic programming references, the form (x,y,z) in
# ExponentialCone() corresponds to: x >= y * e^(z/y)
# HOWEVER, TO WRITE THE CODE IN MathOptInterface, based on REF1, Mosek knows
# another standard form. For example, x >= y * e^(z/y) is modeled as follows
# in mosek:
# @constraint(model, [ z ; y ; x ] in MOI.ExponentialCone()  )
#
# REF1 (paper): MathOptInterface: a data structure for mathematical
# optimization problems. Legat et al (2020)

function GP(solver) # since this example is small, we directly use
                    # the example date, instead of passing the date
    model= Model(optimizer_with_attributes(solver))

    @variable(model, t)
    @variable(model, p[1:3])
    @variable(model, q[1:3])
    @variable(model, u )
    @variable(model, v)
    @variable(model, w)
# (o)
    @constraint(model, [ - p[1]-q[1]+1 ] in MOI.Nonnegatives(1))
    @constraint(model, [ u-t; 1; p[1]] in MOI.ExponentialCone())
    @constraint(model, [ 2*v+w-t; 1; q[1]] in MOI.ExponentialCone())
# Note: No need to define dimentions for Exponential Conic constraints

# (1)
    @constraint(model, [ -p[2]-q[2]+1 ] in MOI.Nonnegatives(1))
    @constraint(model, [ 0.5*u + log(0.1); 1;p[2] ] in MOI.ExponentialCone())
    @constraint(model, [ -v + log(2) ; 1; q[2] ] in MOI.ExponentialCone())
# Note: the base for log(.) function is ℯ by default

# (2)
    @constraint(model, [ -p[3]-q[3]+1 ] in MOI.Nonnegatives(1))
    @constraint(model, [ -w; 1; p[3] ] in MOI.ExponentialCone())
    @constraint(model, [ v-2*u; 1; q[3] ] in MOI.ExponentialCone())

    @objective(model, Min, t)

    optimize!(model)

    ### Recall to return x and y and z based on u, v, w so that:
    #       x = e^u,    y = e^v,    z = e^w
    termination_status(model),objective_value(model),ℯ^value(u),ℯ^value(v),ℯ^value(w)

end





# Call solver
solver = Mosek.Optimizer
# use Cholesky factorization method for matrix decomposition
LinearAlgebra.Cholesky


## Data:
# Data for example 1
 c = [8 , 10, 9.0]
 b = [14, 12.5]
 A = [1     3       2;
      1     5       3.0]

#Data for example 2
c2 = [5 , 1.]
a2 = [4 , 3.]
b2 = 22.
Q = [5   -1.;  -1   4]
# If we have matrix F such that Q= F*F', we can directly use F and pass it to
# the function. Otherwise, the following code uses the cholesky decomposition
# to find F.
w=size(Q,1)
Dec1 = cholesky(Q)
SaveF = Dec1.L
F = zeros(w,w)
for i=1:w
    for j=1:w
        if (SaveF != 0)
            F[i,j]=SaveF[i,j]
        end
    end
end
# Note: F*F' = Dec1.L *Dec1.U = Q

# Additional data for example 3
#P = [22. 10. -16.; 10. 37. -43.; -16. -43. 98.]
P = [3   1;  1   2]
Dec2 = cholesky(P)
SaveV = Dec2.L
w2= size(P,1)
V = zeros(w2,w2)
for i=1:w2
    for j=1:w2
        if (SaveV != 0)
            V[i,j]=SaveV[i,j]
        end
    end
end
# note: V*V' = Dec1.L *Dec1.U = P


# Data of example 4 were directly implemmented in the code

# Solve problems
s1,z1,x1 = LP(solver,A,c,b)
s2,z2,x2 = QP(solver, F, c2, a2, b2)
s3,z3,x3 = QCQP(solver, F, V, c2 , a2, b2)
s4,z4,xx,yy,zz = GP(solver)

# Output of LP
println("Solution of the LP")
println("Status: ", s1)
println("objective value = $z1")
println("  x = $x1")
println("")
println("")
# Output of QP
println("Solution of the QP")
println("Status: ", s2)
println("objective value = $z2")
println("  x =  ", x2 , " (integer)")
#println("  x = $x2")
println("")
println("")
# Output of QCQP
println("Solution of the QCQP")
println("Status: ",s3)
println("objective value = $z3")
println("  x = $x3")
println("")
println("")
# Output of GP
println("Solution of the GP")
println("Status: ",s4)
println("objective value = $z4")
println(" [ x,y,z ] = [ ", xx , " , ", yy, " , ",zz, "  ]" )
println("")
println("")
println("# Written by: Milad Dehghani Filabadi")
println("")
println("")
