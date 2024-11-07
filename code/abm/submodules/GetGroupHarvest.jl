  # Harvest
      function GetGroupHarvest(effort, loc, K, kmax, tech, labor, degrade, ngroups)
       b = zeros(ngroups)
       X =zeros(ngroups)
       for i in 1:ngroups
         X[i] =tech[i]*((sum(effort[loc .== i])^labor[i])*K[i]^degrade[i])
       end
       return(X)
      end


        # Harvest
        function GetGroupHarvest(effort, loc, K, kmax, tech::Float64, labor::Float64, degrade::Float64, ngroups)
         X =zeros(ngroups)
          for i in 1:ngroups
            X[i] =tech*((sum(effort[loc .== i])^labor)*K[i]^degrade)
          end
          return(X)
         end
   