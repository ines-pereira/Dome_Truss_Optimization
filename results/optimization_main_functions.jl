function load_results(filenames; base_folder=results_folder,
                      has_header=true, keep_header=false,
                      sep=files_sep, usecols=collect(map(x->x[1], names_mapping)),
                      max_lines=max_evals)
    """Loads the data from the specified `base_folder` using the `filenames`.
    It assumes the filenames
    """

    filepaths = [joinpath(base_folder, f) for f in filenames]

    limit_csv(csv) = csv[1:min(size(csv, 1), max_lines), usecols]

    [limit_csv(CSV.read(f)) for f in filepaths]

end


function get_run_indices(dfs, run, n_algorithms=n_algorithms)
    """Returns the dataframes that correspond to the specified `run`.
    This assumes that the dataframes are read per run, that means that
    if we run two algorithms (SMPSO and SPEA2) for 3 runs each, this
    method assumes that it was read in the following order:

    > SMPSO_run1
    > SPEA2_run1
    > SMPSO_run2
    > SPEA2_run2
    > SMPSO_run3
    > SPEA2_run3

    In that case, this invocation `get_run_indices(dfs, 1, n_algorithms=2)`
    will return:
        dfs[0:2]
    """
    run = 1
    dfs[run*n_algorithms:(run+1)*n_algorithms, :]
end

# **IMPORTANT NOTE**: This function assumes that your problem is a minimization problem for every objective dimension.
function weakly_dominates(v0, v1)
    """Computes whether v0 dominates v1, i.e., whether at least one objective
    is better (in this case, smaller) than some other)
    """
    all(t-> t[1]<=t[2], zip(v0, v1)) && any(t-> t[1]<t[2], zip(v0, v1))
#     all(v0 .<= v1) && any(v0 .< v1)
end

function get_non_dominated(V, dominance=weakly_dominates)
    """Computes the optimal and non-optimal solutions.
    Optimal solutions are called non-dominated and non-optimal
    solutions are called denominated.
    """
    nsols, nobjs = size(V)

    dominated = zeros(nsols)
    dominated_by = zeros(nsols)
    for i in 1:nsols
        for j in 1:nsols
            if i != j
                if dominance(V[j, :], V[i, :])
                    dominated[i] = 1
                    dominated_by[i] = j
                    break
                end
            end
        end
    end

   dominated, dominated_by
end

function add_isdominated_cols(d; cols=objs_cols)
    """Adds to the provided dataframe columns for Pareto optimal solution."""
    A = d[:, cols]
    B, C = get_non_dominated(A)
    d[!, :isDominated] = B
    d[!, :dominatedBy] = C
    println(by(d, :isDominated, nrow))
    d
end

function get_combined_PF(dfs; drop_cols=relevant_cols, objs_cols=objs_cols)
    """Computes the combined Pareto front based on a set of input dataframes"""
    all_data = vcat(dfs...)
    if !isempty(drop_cols)
        all_data = unique(all_data, drop_cols)
    end
    add_isdominated_cols(all_data, cols=objs_cols)
end

ensure_array(x) = x isa Union{Array, Tuple} ? x : [x]

function broadcasting_oper(df, cols, oper, value)
    cols = ensure_array(cols)
    for col in cols
        df[:, col] = oper(df[:, col], value)
    end
    df
end

broadcasting_multi(df, cols, value) = broadcasting_oper(df, cols, *, value)

broadcasting_div(df, cols, value) = broadcasting_oper(df, cols, /, value)

function get_symmetric(df, cols)
    """Computes the symmetric value of the provided `cols` and returns a
    copy of the original dataframe where the values of the specified `cols`
    are symmetric."""
    broadcasting_oper(df, cols, *, -1)
end

 """Computes the real value used for the variables during a discrete optimization process."""
unscale(df, cols, minimum, step) =
    broadcasting_oper(df, cols, (v, step)->minimum .+ step .* v, step)

# Default layout for the pareto fronts graphs
layout = Layout(
    template="plotly_white",
    autosize=true,
#     legend=Dict(
#         :orientation=>"h"
#     ),
    legend_orientation="h",
    # Define axis
    xaxis=Dict(
        :autorange=>true,
        :showgrid=>true,
        :zeroline=>false,
        :showline=>true,
        :ticks=>"",
        :showticklabels=>true,
        :tickformat=>"."
    ),
    yaxis=Dict(
        :autorange=>true,
        :showgrid=>true,
        :zeroline=>false,
        :showline=>true,
        :ticks=>"",
        :showticklabels=>true,
        :tickformat=>"."
    )
)

function create_pf(; pf, name, x, y, nd_color="rgb(0,0,255)", ln_width=1.5, marker_size=5.5, d_color=nothing)
    traces = []
    custom_data = []

    # Create the non dominated trace (in a different color, as specified in *nd_color*)
    non_dominated = filter(row->row[:isDominated] == 0, pf)
    dominated = filter(row->row[:isDominated] == 1, pf)
    x_pf = non_dominated[:, x]
    y_pf = non_dominated[:, y]
    x_pf, y_pf = zip(sort(collect(zip(x_pf, y_pf)), by=t->t[1])...)

    # Get information about the objectives and variables - relevant for an iteractive PF
    info = vcat(vars_cols, objs_cols)
    info_pf = non_dominated[:, info] # it returns a DF, in python returns an array
    info_npf = dominated[:, info]

    push!(traces,
          scatter(
            x=x_pf,
            y=y_pf,
            mode="lines+markers",
            name=name,
            # name=name * " NonDominated",
#           name = name,
            opacity=1,

            # Layout do marker
            marker=Dict(
                :color=>nd_color,
                :size=>marker_size
            ),
            line=Dict(
                :color=>nd_color,
                :width=>ln_width
            )))

    push!(custom_data, info_pf)

    if ! isnothing(d_color)
        x_npf = dominated[:, x]
        y_npf = dominated[:, y]

        # Create the dominated trace (in a different color, as specified in *d_color*)
        push!(traces,
              scatter(
                x=x_npf,
                y=y_npf,
                mode="markers",
                name=name * " Dominated",
                # showlegend=false,
                opacity=0.5,

                # Layout do Marker
                marker=Dict(
    #                       :color = d_color,
                        :color=>nd_color,
                        :size=>marker_size * 0.7
                )))

        push!(custom_data, info_npf)
    end

    traces, custom_data
end

function get_traces(pfs, x, y, draw_dominated=true,
               names=all_algorithms, colorscale="viridis", colors=nothing,
               tpf=nothing, tpf_name="Combined_PF", tpf_color="rgb(0,0,0)",
               layout=layout)

    pfs = ensure_array(pfs)
    names = ensure_array(names)
    n_pfs = length(pfs)
    colors = isnothing(colorscale) ? colors : get_colors(n_pfs, colorscale)

    traces = []
    custom_data = []

    if ! isnothing(tpf)
        sub_traces, sub_custom_data = create_pf(pf=tpf, name=tpf_name, x=x, y=y, ln_width=4, marker_size=10, nd_color=tpf_color)
        append!(traces, sub_traces)
        append!(custom_data, sub_custom_data)
    end

    for (i, pf) in enumerate(pfs)
        d_color = draw_dominated ? colors[i] : nothing
        sub_traces, sub_custom_data = create_pf(pf=pf, name=names[i], x=x, y=y, nd_color=colors[i], d_color=d_color)
        append!(traces, sub_traces)
        append!(custom_data, sub_custom_data)
    end

    traces, custom_data
end

"""
Right now this this function is adapated to be used for the waved façade case study (PLEA 2020).

If needed change/comment the variables after the on function.
"""
function create_pfs(pfs; x, y, draw_dominated=true,
               names=all_algorithms, colorscale="viridis", colors=nothing,
               tpf=nothing, tpf_name="Combined_PF", tpf_color="rgb(0,0,0)",
               layout=layout)

    traces, custom_data = get_traces(pfs, x, y, draw_dominated, names, colorscale, colors, tpf, tpf_name, tpf_color, layout)

    fig = PlotlyJS.plot([traces...], layout)

    # display(fig)
    on(fig.scope["click"]) do data
        let subdata = data["points"][1],
            point_idx = subdata["pointIndex"],
            curve_idx = subdata["curveNumber"],
            info = custom_data[curve_idx + 1][point_idx + 1, :]
    #         m_stripes = Int(info[:m_stripes]),
    #         rotation = info[:rotation],
    #         thickness = info[:thickness]
    #
          println(info)
    #       facade_pf(; m_stripes_south03=m_stripes, angle_rotation=rotation, f_thicken=thickness)
        end
    end
    fig
end
