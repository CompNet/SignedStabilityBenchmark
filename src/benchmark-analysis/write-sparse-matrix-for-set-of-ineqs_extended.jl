using DelimitedFiles;
using DataStructures;
using SparseArrays;
using JLD;


# ==========================================
# vec: float vector where the indexes represent
#       edge variables. The values correspond to the best linear relaxation
# optSol: binary vector where the indexes represent
#       edge variables. The values correspond to opt sol
# ==========================================
function is_constraint_equal(vec, value, optSol)
    res = vec' * optSol;
    #println(round(res[1],digits=2));
    if round(res[1],digits=2) == value
		return(true);
    end
	return(false);
end



# ==========================================
#  Prepares the names of the variables, which corresponds to the column names of a matrix.
#   This function assumes that the vertex pairs are used as variables in the ILP formulation
#   (contrary to edge variables)
# ==========================================
function prepare_colnames(n)
    nb_pair = Int32((n*(n-1))/2);
    colnames = Array{String}(undef, nb_pair);

    counter = 0;
    for i in 1:(n-1)
        for j in (i+1):n
            counter += 1;
            colnames[counter] = string("x_",j-1,"_",i-1) ## first index is always greater than the second index
        end
    end
    return(colnames);
end


# ==========================================
# Prepares a single row vector of the matrix in output. Note that there are as many rows as constraints.
# This function assumes that the vertex pairs are used as variables in the ILP formulation
#   (contrary to edge variables)
#   So, there must be (n*(n-1)*(n-2))/2 triangle ineqs in the model.
#   The model otput by ExCC already satisfies this requirement.
# ==========================================
function prepare_row_vec(n, vars, coefs, colnames)
    nb_pair = Int32((n*(n-1))/2);
    vec = Array{Int32}(zeros(nb_pair));
    row_vec = OrderedDict(zip(colnames,vec));
    for i in 1:size(vars,1)
        var = vars[i];
        coef = coefs[i];
        row_vec[var] = Int32(coef);
    end
    return(row_vec);
end



# ==========================================
# Handles integrality ineqs. There are (n*(n-1))/2 such ineqs. 
#
# Outputs mainly three vectors I, J and V, which is inline with the sparse matrix format.
#   Vector I: row index values for the matrix
#   Vector J: column index values for the matrix. Note that there are (n*(n-1))/2 columns, e.g. x_10, x_20, ..
#   Vector V: the coefficient of the ij.th entry, i.e. variable, in the matrix
# ==========================================

function handle_integrality_ineqs(n, colnames, varsOptSol)
    I = Int32.([]);
    J = Int32.([]);
    V = Int32.([]);
    bounds = Int32.([]);

    counter = 0;
    rhs_coef = 0;
    for i in 0:(n-2)
	    for j in (i+1):(n-1)
		    term_str = string("x_",j,"_",i);
            vars = [term_str];
            coefs = [-1];
		    new_row = prepare_row_vec(n, vars, coefs, colnames);
		    if is_constraint_equal(new_row.vals, rhs_coef, varsOptSol)
                counter += 1;
                sp = sparsevec(new_row.vals);
		        append!(I, Array{Int32}(fill(counter, size(coefs,1))));
                append!(J, Array{Int32}(sp.nzind));
                append!(V, Array{Int32}(sp.nzval));
                append!(bounds, 0);
		    end
	    end
    end

    return I, J, V, bounds, counter;
end



# ==========================================
# Handles triangle ineqs. This function assumes that
#   the redundant triangle inequalities are not filtered from the ILP model.
#   So, there must be (n*(n-1)*(n-2))/2 triangle ineqs in the model.
#   The model otput by ExCC already satisfies this requirement.
#
# Outputs mainly three vectors I, J and V, which is inline with the sparse matrix format.
#   Vector I: row index values for the matrix
#   Vector J: column index values for the matrix. Note that there are (n*(n-1))/2 columns, e.g. x_10, x_20, ..
#   Vector V: the coefficient of the ij.th entry, i.e. variable, in the matrix
# ==========================================
function handle_triangle_ineqs(n, colnames, varsOptSol, counter)
    I = Int32.([]);
    J = Int32.([]);
    V = Int32.([]);
    bounds = Int32.([]);

    rhs_coef = 1;
    for i in 0:(n-3)
	    for j in (i+1):(n-2)
		    for k in (j+1):(n-1)
			    term1_str = string("x_",j,"_",i)
			    term2_str = string("x_",k,"_",i)
			    term3_str = string("x_",k,"_",j)
                vars = [term1_str, term2_str, term3_str];

                # ------- 1 --------
                coefs = [1,1,-1];
                new_row = prepare_row_vec(n, vars, coefs, colnames);
		        if is_constraint_equal(new_row.vals, rhs_coef, varsOptSol)
                    counter += 1;
                    sp = sparsevec(new_row.vals);
		            append!(I, Array{Int32}(fill(counter, size(coefs,1))));
                    append!(J, Array{Int32}(sp.nzind));
                    append!(V, Array{Int32}(sp.nzval));
                    append!(bounds, 1);
		        end

                # ------- 2 --------
                coefs = [1,-1,1];
                new_row = prepare_row_vec(n, vars, coefs, colnames);
		        if is_constraint_equal(new_row.vals, rhs_coef, varsOptSol)
                    counter += 1;
                    sp = sparsevec(new_row.vals);
		            append!(I, Array{Int32}(fill(counter, size(coefs,1))));
                    append!(J, Array{Int32}(sp.nzind));
                    append!(V, Array{Int32}(sp.nzval));
                    append!(bounds, 1);
		        end

                # ------- 3 --------
                coefs = [-1,1,1];
                new_row = prepare_row_vec(n, vars, coefs, colnames);
		        if is_constraint_equal(new_row.vals, rhs_coef, varsOptSol)
                    counter += 1;
                    sp = sparsevec(new_row.vals);
		            append!(I, Array{Int32}(fill(counter, size(coefs,1))));
                    append!(J, Array{Int32}(sp.nzind));
                    append!(V, Array{Int32}(sp.nzval));
                    append!(bounds, 1);
		        end
                
            end
        end
    end

    return I, J, V, bounds, counter;
end





# ===============================================
# Handles other ineqs: those inserted into the model during the Cutting Plane method in ExCC
#   This function assumes that the redundant triangle inequalities are not filtered from the ILP model.
#   So, there must be (n*(n-1)*(n-2))/2 triangle ineqs in the model.
#   The model otput by ExCC already satisfies this requirement.
#
# Outputs mainly three vectors I, J and V, which is inline with the sparse matrix format.
#   Vector I: row index values for the matrix
#   Vector J: column index values for the matrix. Note that there are (n*(n-1))/2 columns, e.g. x_10, x_20, ..
#   Vector V: the value of the ij.th entry in the matrix
# ===============================================
function handle_other_ineqs(in_folder, n, colnames, varsOptSol, lp_filepath, counter)

    I = Int32.([]);
    J = Int32.([]);
    V = Int32.([]);
    bounds = Int32.([]);

    nb_triple = Int32((n*(n-1)*(n-2))/6)*3;

    lines = readlines(lp_filepath);

    process1 = false;
    process = false;
    is_whole_ineq = false;
    ineq = "";
    ineq_ids = [];
    inner_counter = 0;
    for line in lines
        if line == "Subject To"
            process1 = true;
        end

	    if process1
		    inner_counter = inner_counter + 1;
        end
        #println(inner_counter);

        if process && line != "Bounds"
            is_begining_ineq = any(contains.(line,":"));
            if(is_begining_ineq)
                ineq = line;
            else #if(!is.begining.ineq)
                ineq = string(ineq, line);
            end

            if( any(contains.(line,"<=")) )
                is_whole_ineq = true;
            end
            
            if is_whole_ineq
                ineq = replace(ineq, r" +" => "");
                ineq_id = split(ineq, ":")[1];
                push!(ineq_ids, ineq_id);

                main_part = split(ineq, ":")[2];
                rhs_coef = split(main_part, "<=")[2];
                vars = split(main_part, "<=")[1];
                vars = replace(vars, r"\+" => ";+");
                vars = replace(vars, r"\-" => ";-");
                var_list = split(vars, ";");
                if var_list[1] == ""
                   deleteat!(var_list, 1);
                end

                # prepare the coefficients
                coefs = Int32.([]);
                vars = [];
                for var in var_list
                    coef = split(var, "x_")[1];
                    if coef == "" || coef == "+"
                        push!(coefs, +1);
                    else
                        push!(coefs, -1);
                    end

                    var = string("x_",split(var, "x_")[2]);
                    push!(vars, var);
                end

                # prepare row vector
                new_row = prepare_row_vec(n, vars, coefs, colnames);

               #if is_constraint_equal(new_row.vals, rhs_coef, varsOptSol)
                    counter += 1;
                    sp = sparsevec(new_row.vals);
		            append!(I, Array{Int32}(fill(counter, size(coefs,1))));
                    append!(J, Array{Int32}(sp.nzind));
                    append!(V, Array{Int32}(sp.nzval));
                    append!(bounds, parse(Int,rhs_coef));
		        #end

                is_whole_ineq = false;
            end
            
        end


	    if process1 && inner_counter>nb_triple
		    process = true;
        end

        if line == "Bounds"
		    process1 = false;
            process = false;
	    end

    end

    return I, J, V, bounds, counter;
end




####################################################################################
# Main
####################################################################################


nb_args = size(ARGS,1);

if nb_args==5
    in_folder = ARGS[1]
    out_folder = ARGS[2]
    graph_filepath = ARGS[3]
    membership_filepath = ARGS[4]
    lp_filepath = ARGS[5]
    # frac_graph_filepath = ARGS[6]

    if isdir(in_folder) && isdir(out_folder)

        open(`head -n1 $graph_filepath`) do io
            global n = convert(Int,readdlm(io,'\t')[1,1]);
        end
        println(n)

        membership = Array{Int32}(readdlm(membership_filepath));

        varsOptSol = Int32.([]);
        for i in 1:(n-1)
            for j in (i+1):n
                val = 0;
                if membership[i] == membership[j]
                    val = 1;
                end
                push!(varsOptSol, val);
            end
        end

        open(string(out_folder,"/vars.txt"), "w") do io
           writedlm(io, varsOptSol, ',')
        end

        # ------------------------------------

        #varsOptSol = readdlm("vars.txt",',');
        colnames = prepare_colnames(n);
        nb_pair = Int32((n*(n-1))/2);
        nb_triple = Int32((n*(n-1)*(n-2))/6);

        I1, J1, V1, bounds1, counter = handle_integrality_ineqs(n, colnames, varsOptSol);
        I2, J2, V2, bounds2, counter = handle_triangle_ineqs(n, colnames, varsOptSol, counter);
        I3, J3, V3, bounds3, counter = handle_other_ineqs(in_folder, n, colnames, varsOptSol, lp_filepath, counter);

        I = Array{Int32}(vcat(I1, I2, I3));
        J = Array{Int32}(vcat(J1, J2, J3));
        V = Array{Int32}(vcat(V1, V2, V3));
        b = Array{Int32}(vcat(bounds1, bounds2, bounds3));

        nrow = counter;
        ncol = size(colnames,1);
        S = sparse(I,J,V,nrow,ncol);
        save(string(out_folder,"/S.jld"),"x",S)

        # ---------------------------------

        #S2 = load("/home/nejat/S.jld","x");
        #b = Array{Int32}(vcat(fill(0, nb_pair), fill(1, nb_triple), bounds));
        open(string(out_folder,"/b.txt"), "w") do io
           writedlm(io, b, ',')
        end

    else
        println("one of the input or output folders does not exits: ", out_folder);
    end

else
    println("you should correctly indicate the input parameters ...");
end

