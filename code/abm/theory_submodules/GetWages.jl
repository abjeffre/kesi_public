function GetWages(effort, wages::Vector{Float64}, agents::DataFrame = agents) wages[agents.gid].*effort[:,1] end 
function GetWages(effort, wages::Float64) wages.*effort[:,1] end
function GetWages(effort, wages::Float64, n, ngroups, labor_market) 
    (1-labor_market*(sum(effort[:,1])/(n*ngroups)))*wages.*effort[:,1] 
end
