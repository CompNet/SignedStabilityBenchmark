using DelimitedFiles;
using DataStructures;
using SparseArrays;
using JLD;
using JuMP, CPLEX;



# Our goal is to increase the amount of imbalance: so we want to invert the sign of the edges correctly classified, i.e. well-placed ones
# TODO: Since there can be many optimal solutions, we need to know all optimal solutions. Maybe, in one solution, the link seems to be well-placed, whereas it is not in another one


# ==========================================
# Vectorizes the weights. The length of this vector is (n*(n-1))/2.
# ==========================================
function vectorize_weights(n, graph_filepath, weights)
    G = readdlm(graph_filepath,'\t', skipstart=1);
    E = Int32.(G[:,1:2]);

    m = Int32.(zeros(n,n));
    for i in 1:size(weights,1)
        m[E[i,1]+1,E[i,2]+1] = weights[i];
        m[E[i,2]+1,E[i,1]+1] = weights[i];
    end

    new_weights = Int32.([]);
    for i in 1:(n-1)
        for j in (i+1):n
            append!(new_weights, m[i,j]);
        end
    end

    return new_weights;
end



# ==========================================
# Checks if the perturbed weight vector satisfy the optimality condition
#   i.e. whether it is in the stability range or not. 
#
# See this reference for more information:
#   Nowozin, S., Jegelka, S.: Solution stability in linear programming relaxations. 
#        In: Proceedings of the 26th Annual International Conference on Machine Learning - ICML. ACM Press (2009).
#        DOI 10.1145/1553374.1553473
# ==========================================
function check_optimality(A, w, cplex_threads)
    m, n = A.m, A.n;

    #model = Model(CPLEX.Optimizer);
    model = Model(with_optimizer(CPLEX.Optimizer));
    MOI.set(model, MOI.Silent(), true);
    with_optimizer(CPLEX.Optimizer, logLevel = 0)
    set_optimizer_attribute(model, "CPXPARAM_Threads", cplex_threads) # TODO I am not sure if it makes the execution faster

    @objective(model, Min, 0);
    @variable(model, x[1:m] >= 0);
    @constraint(model, A' * x .== w);

    status = optimize!(model); # solves the model
    status = Int(termination_status(model));

    if(status == 1)
        return true;
    end
    return false;
end


# ==========================================
# Writes the current valid, i.e. optimal, weight vector into a file.
#   This prevents from having information loss if the program is interrupted for some reasons.
# ==========================================
function write_perturbed_graph_file(graph_filepath, pert_graph_filepath, w_perturbed)
    first_line = open(graph_filepath) do f
        readline(f)
    end

    # 1) write into file without the first line
    graph_content = readdlm(graph_filepath,'\t', skipstart=1);
    graph_content[:,3] = w_perturbed;
    open(pert_graph_filepath, "w") do io 
        writedlm(io, Int32.(graph_content), '\t') 
    end 

    # 2) then append the first line
    s=open(pert_graph_filepath, "a+");
    seekstart(s)
    write(s,first_line,"\n");
    lines = readlines(pert_graph_filepath);
    for line in lines
        write(s,line,"\n");
    end
    close(s);
end


# ==========================================
# Returns all the indexes susceptible to be perturbed for internal edges.
# ==========================================
function get_remaining_indexes_for_internal_links(membership, source_ids, target_ids)
    remaining_indxs = [];
    for i in 1:size(source_ids,1)
       if membership[source_ids[i]+1] == membership[target_ids[i]+1] # inside cluster
            push!(remaining_indxs,i);
        end
    end
    return remaining_indxs;
end


# ==========================================
# Returns all the indexes susceptible to be perturbed for external edges.
# ==========================================
function get_remaining_indexes_for_external_links(membership, source_ids, target_ids)
    remaining_indxs = [];
    for i in 1:size(source_ids,1)
       if membership[source_ids[i]+1] != membership[target_ids[i]+1] # inside cluster
            push!(remaining_indxs,i);
       end
    end
    return remaining_indxs;
end

# ==========================================
# Returns all the indexes susceptible to be perturbed for internal and external edges.
# ==========================================
function get_remaining_indexes_for_all_links(membership, source_ids, target_ids)
    remaining_indxs = [];
    for i in 1:size(source_ids,1)
            push!(remaining_indxs,i);
    end
    return remaining_indxs;
end


# ==========================================
# # TODO take all opt solutions into account
# ==========================================
#function get_remaining_indexes_ver2(source_ids, target_ids, weights, membership)
#    remaining_indxs = [];
#    for i in 1:size(source_ids,1)
#        ok1 = membership[source_ids[i]+1] != membership[target_ids[i]+1] && weights[i]<0; # external well-placed links
#        ok2 = membership[source_ids[i]+1] == membership[target_ids[i]+1] && weights[i]>0; # internal well-placed links
#        if ok1 || ok2
#            push!(remaining_indxs,i);
#        end
#    end
#    return(remaining_indxs);
#end



# ==========================================
# This is the main method, which perturbs the weights of the initial perfectly balanced signed graph.
#   Note that we perturb only the weights of the existing edges, i.e. we dont add any new links if the graph is not complete.
#
# During the perturbation process, we do not try to perturb every edge weight. If there is no progress at some point, i.e. the last pertrubation makes the new weight vector suboptimal, we attempt to change this situation 'n' times, where n is the graph order.
#
# In this function, there are 3 ways to perturb the edge weights.
#   - "internal": perturbing only the weights of the internal edges. So, this strategy increases the proportion of negative edges.
#   - "external": perturbing only the weights of the external edges. So, this strategy increases the proportion of positive edges.
#   - "all": the hybrid strategy. It perturbs the weights of internal and external edges. It does so by alternating: 1) internal, 2) external, 3) internal, etc. 
#
# ==========================================
function perturb_signed_graph(out_folder, graph_filepath, membership_filepath, vars_filepath, constrs_filepath, bounds_filepath, nb_attemps_for_termination, target_link_type, n, cplex_threads)
    pert_graph_filepath = string(out_folder, "/perturbed_graph_",target_link_type,".G");

    #--------------------
    # LOAD DATA
    #--------------------
    b = readdlm(bounds_filepath,',');
    #A = convert(Array{Float64,2}, readdlm("A.txt",',', skipstart=1)[:,2:end]);
    A = load(constrs_filepath,"x"); # it is a sparse matrix
    #w = vectorize_weights(n, graph_filepath);
    varsOptSol = readdlm(vars_filepath,',');
    #varsOptSol = readdlm(string(in_folder,"/fractionalGraph.G"),'\t', skipstart=1)[:,3];
    membership = Int32.(readdlm(membership_filepath,','));

    println("data loaded !!!");

    w = readdlm(graph_filepath,'\t', skipstart=1)[:,3];
    w_vec = vectorize_weights(n, graph_filepath, w);
    start = check_optimality(A, w_vec, cplex_threads)
    println(start)

    if(start)
        global nb_attemps = nb_attemps_for_termination; 

        println("start to perturb the signed graph ...")
        source_ids = Int32.(readdlm(graph_filepath,'\t', skipstart=1)[:,1]);
        target_ids = Int32.(readdlm(graph_filepath,'\t', skipstart=1)[:,2]);

        global remaining_indxs = get_remaining_indexes_for_internal_links(membership, source_ids, target_ids);
        if target_link_type == "external"
            global remaining_indxs = get_remaining_indexes_for_external_links(membership, source_ids, target_ids);
        elseif target_link_type == "all"
             global remaining_indxs = get_remaining_indexes_for_all_links(membership, source_ids, target_ids);
        end

        w_perturbed = readdlm(graph_filepath,'\t', skipstart=1)[:,3];
        global delta = 0;
        indxs = []

        while size(remaining_indxs,1)>0 && nb_attemps>0
            #sleep(1);

            if mod(size(remaining_indxs,1),2) == 0
                println("writing into file")
                write_perturbed_graph_file(graph_filepath, pert_graph_filepath, w_perturbed)
            end

            global w_perturbed = readdlm(graph_filepath,'\t', skipstart=1)[:,3];
            if size(indxs,1) != 0
                global w_perturbed[indxs] = -w_perturbed[indxs];
                # remove the indexes that are already selected => preparation for the next rand()
                for i in indxs
                    global remaining_indxs = filter(x -> x!=i, remaining_indxs);
                end
            end

            # allow only links having negative weights located between clusters
            global i = 0;
            while size(remaining_indxs,1)>0
                global i = rand(remaining_indxs);
                #if(w[i]<0 && membership[source_ids[i]+1] != membership[target_ids[i]+1])
                #if(membership[source_ids[i]+1] != membership[target_ids[i]+1])
                println("indx:",i,", source_id:",source_ids[i],", target_id:",target_ids[i],", init weight:",w[i]," (nb current:",size(indxs,1),", nb total:",size(remaining_indxs,1),")")
                break
                #end
            end; 
            #i = rand(remaining_indxs);
            global w_perturbed[i] = -w_perturbed[i];
            
            
            new_w_perturbed = vectorize_weights(n, graph_filepath, w_perturbed);
            if check_optimality(A, new_w_perturbed, cplex_threads)
                global nb_attemps = nb_attemps_for_termination;
                println("ok -> indx:",i,", source_id:",source_ids[i],", target_id:",target_ids[i],", init weight:",w[i])
                push!(indxs,i);

                if membership[source_ids[i]+1] == membership[target_ids[i]+1] # note that vertex ids start from 0
                    global delta += w[i]; # if w[i]>0, then it will be negative. So, this will increase imbalance
                elseif membership[source_ids[i]+1] != membership[target_ids[i]+1]
                    global delta -= w[i]; # if w[i]>0, then it will be negative. So, this will decrease imbalance
                end
                #break;
            else
                global nb_attemps -= 1;
                println("nb attemps:",nb_attemps,",");
                global remaining_indxs = filter(x -> x!=i, remaining_indxs);
            end

        end
        println("end. Delta: ", delta,". Writing into file ...")

        global w_perturbed = readdlm(graph_filepath,'\t', skipstart=1)[:,3];
        if size(indxs,1) != 0
            global w_perturbed[indxs] = -w_perturbed[indxs];
        end
        write_perturbed_graph_file(graph_filepath, pert_graph_filepath, w_perturbed)

    else
        println("there is a problem, without perturbation. Verify the input vectors..");  
    end

end




# ============================================================================



# ==========================================
# Main
# ==========================================


nb_args = size(ARGS,1);


if nb_args==8
    out_folder = ARGS[1];
    graph_filepath = ARGS[2];
    membership_filepath = ARGS[3];
    vars_filepath = ARGS[4];
    constrs_filepath = ARGS[5];
    bounds_filepath = ARGS[6];
    cplex_threads = parse(Int32, ARGS[7])
    target_link_type = ARGS[8] # 3 OPTIONS AVAILABLE: "internal", "external", "all"


    open(`head -n1 $graph_filepath`) do io
        global n = convert(Int,readdlm(io,'\t')[1,1]);
    end
    println(n)

    nb_attemps_for_termination = n; #2*n;

    if isdir(out_folder)

        perturb_signed_graph(out_folder, graph_filepath, membership_filepath, vars_filepath, constrs_filepath, bounds_filepath, nb_attemps_for_termination, target_link_type, n, cplex_threads);

    else
        println("one of the input or output folders does not exits: ", out_folder);
    end

else
    println("you should correctly indicate the input parameters ...");
end
