"""
    average_time_series( time_series_name::String, simulation_name, param_name, param_vals, repetitions_per_param )

    Averages one or an array of time series with name 'time_series_name' from the simulation 'simulation_name'
    where 'param_name' is varied over 'param_vals' with 'repetitions_per_param' repetitions per parameter (surprise).

"""
function average_time_series( time_series_name::String, simulation_name, param_name, param_vals, repetitions_per_param )
    local maxits = load("data/$(simulation_name)/$(param_name)=$(param_vals[1])_repet=1.jld2","maxits")
    time_series = zeros(maxits,length(param_vals));
    for (param_val_it,param_val) in enumerate(param_vals)
        for repet in 1:repetitions_per_param
            file_name = "data/$(simulation_name)/$(param_name)=$(param_val)_repet=$(repet).jld2"
            time_series[:,param_val_it] += load(file_name,time_series_name)/repetitions_per_param;
        end
    end
    return time_series;
end

function average_time_series( time_series_names::Vector{String}, simulation_name, param_name, param_vals, repetitions_per_param )
    local maxits,clock = load("data/$(simulation_name)/$(param_name)=$(param_vals[1])_repet=1.jld2","maxits","clock")
    time_series = NTuple{length(time_series_names),Array{Float64}}(zeros(maxits,length(param_vals)) for i in 1:length(time_series_names) )
    for (param_val_it,param_val) in enumerate(param_vals)
        for repet in 1:repetitions_per_param
            file_name = "data/$(simulation_name)/$(param_name)=$(param_val)_repet=$(repet).jld2"
            temp_ts = load(file_name,time_series_names...)
            for (i,ts) in enumerate(temp_ts)
                time_series[i][:,param_val_it] += ts/repetitions_per_param;
            end
        end
    end
    return time_series;
end