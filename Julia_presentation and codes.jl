# Written by: Milad Dehghani Filabadi


# Example 1: Linear programming model

#𝑀𝑎𝑥     8𝑥_1+10𝑥_2+9𝑥_3
#𝑠.𝑡.        𝑥_1+3𝑥_2+2𝑥_3≤14
#           𝑥_1+5𝑥_2+3𝑥_3≤12.5
#             𝑥_1,𝑥_2,𝑥_3≥0


# Example 2: Second-order conic programming
    # 𝑀𝑖𝑛     𝑦_1+𝑦_2
    # 𝑠.𝑡.     𝑦_1^2 + 𝑦_2^2 ≤ 4
    #         (𝑦_1-3)^2 + 𝑦_2^2 ≤ 4
    #         𝑒^(𝑦_1−2) ≤ 𝑦_2
    #         𝑦_1,𝑦_2  free 

# Example 3: Conic programming: Second-order cone and exponential cone
    # 𝑀𝑖𝑛     𝑦_1+𝑦_2
    # 𝑠.𝑡.     𝑦_1^2 + 𝑦_2^2 ≤ 4
    #         (𝑦_1-3)^2 + 𝑦_2^2 ≤ 4
    #         𝑒^(𝑦_1−2) ≤ 𝑦_2
    #         𝑦_1,𝑦_2  free 


using JuMP
using Mosek
using MosekTools

#  
function solve_LP()
    # Define LP coefficients
    c = [8, 10, 9.0]
    b = [14, 12.5]
    A = [1 3 2; 1 5 3.0]

    # Create a JuMP model with solver
    model = Model(optimizer_with_attributes(Mosek.Optimizer))

    # Define variables
    n = length(c)
    @variable(model, x[1:n] >= 0)

    # Define the objective function
    @objective(model, Max, sum(c .* x))

    # Define the constraints
    @constraint(model, A * x .<= b)
    for i in 1:2
        @constraint(model, A[i,:]'*x <= b[i])
    end

    # Solve the LP problem
    optimize!(model)

    # Display the results
    println("SOLVE LP EXAMPLE 1")
    print("optimal x: ",  value.(x))
    print("optimal z: ", objective_value(model))
    print("status: ", termination_status(model))
end

function solve_so_conic_program()
    # Example 3:
    # 𝑀𝑖𝑛     𝑦_1+𝑦_2
    # 𝑠.𝑡.     𝑦_1^2 + 𝑦_2^2 ≤ 4
    #         (𝑦_1-3)^2 + 𝑦_2^2 ≤ 4
    #         𝑒^(𝑦_1−2) ≤ 𝑦_2
    #         𝑦_1,𝑦_2  free 
    model = Model(optimizer_with_attributes(Mosek.Optimizer))

    @variable(model, x[1:2] )

    @objective(model, Min, x[1]+x[2])

    @constraint(model, [2;x[1];x[2]] in MOI.SecondOrderCone(3))
    @constraint(model, [2;x[1]-3;x[2]] in MOI.SecondOrderCone(3))


    @constraint(model, [x[1]-2;1;x[2]] in MOI.ExponentialCone())
    
    optimize!(model)
    println("SOLVE Conic Program of EXAMPLE 3")
    print("optimal x: ",  value.(x))
    print("optimal z: ", objective_value(model))
    print("status: ", termination_status(model))
end


function solve_exp_conic_program()
    # Example 3:
    # 𝑀𝑖𝑛     𝑦_1+𝑦_2
    # 𝑠.𝑡.     𝑦_1^2 + 𝑦_2^2 ≤ 4
    #         (𝑦_1-3)^2 + 𝑦_2^2 ≤ 4
    #         𝑒^(𝑦_1−2) ≤ 𝑦_2
    #         𝑦_1,𝑦_2  free 
    model = Model(optimizer_with_attributes(Mosek.Optimizer))

    @variable(model, x[1:2] )

    @objective(model, Min, x[1]+x[2])

    @constraint(model, [2;x[1];x[2]] in MOI.SecondOrderCone(3))
    @constraint(model, [2;x[1]-3;x[2]] in MOI.SecondOrderCone(3))


    @constraint(model, [x[1]-2;1;x[2]] in MOI.ExponentialCone())
    
    optimize!(model)
    println("SOLVE Conic Program of EXAMPLE 3")
    print("optimal x: ",  value.(x))
    print("optimal z: ", objective_value(model))
    print("status: ", termination_status(model))
end

# call which problem you would like to solve
solve_LP()
solve_so_conic_program()
solve_exp_conic_program()
