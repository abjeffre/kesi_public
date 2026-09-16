function GetCongestion(effort, loc, ngroups, congestion)
X =zeros(ngroups)
for i in 1:ngroups
    X[i] =sum(effort[loc .== i])*congestion
end
return(X)
end
