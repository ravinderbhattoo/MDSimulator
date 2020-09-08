export write_trajectory, paraview_save

function write_trajectory(filename, res::DiffEqBase.ODESolution)
    if typeof(res) <: DiffEqBase.ODESolution
        for i in 1:length(res.t)
            t = res.t[i]
            data = Array(reshape(Array(res(t)),(3,:)))'
            N = Int(size(data, 1)/2)
            v = data[1:N,:]
            x = data[N+1:end,:]
            writedata("$(filename)_$i.data", (x, v))
        end
    else
        throw("Object is not of type DiffEqBase.ODESolution.")
    end
    print("Trajectory has been written.")
end

function save_trajectory(filename, res::DiffEqBase.ODESolution)
    jldopen("$filename", "w") do file
        g = g_create(file, "res") # create a group
        g["u"] = res.u              # create a scalar dataset inside the group
        g["t"] = res.t
        g["interatomic_potentials"] = res.prob.p[2].interatomic_potentials
        g["boundary_condition"] = res.prob.p[2].boundary_condition
        g["units"] = string(UNITS)
        g["global_vals"] = ustrip(res.prob.p[2].global_vals)
        g["acceleration"] = ustrip(res.prob.p[2].acc)
    end
    return nothing
end


function paraview_save(filename, res::DiffEqBase.ODESolution)
    if typeof(res) <: DiffEqBase.ODESolution
        pvd = paraview_collection(filename)
        mkpath("./files_$filename/")
        for i in 1:length(res.t)
            t = res.t[i]
            data = Array(reshape(Array(res(t)),(3,:)))'
            N = Int(size(data, 1)/2)
            v = data[1:N,:]
            x = data[N+1:end,:]
            a = res.prob.p[2].acc[:,:,i]'
            vtkfile = vtk_grid("./files_$filename/frame_$i", x[:,1], x[:,2], x[:,3], MeshCell[])
            vtkfile["Velocity"] = (v[:,1], v[:,2], v[:,3])
            vtkfile["Acceleration"] = (a[:,1], a[:,2], a[:,3])
            pvd[i] = vtkfile
        end
        vtk_save(pvd)
    else
        throw("Object is not of type DiffEqBase.ODESolution.")
    end
    print("Trajectory has been written.")
end


# outfiles = vtk_multiblock("my_vtm_file") do vtm
#     for i in 1:length(xs)
#         begin
#             vtkfile = vtk_grid(vtm, xs[i][:,1], xs[i][:,2], xs[i][:,3], MeshCell[])
#             vtkfile["Velocity"] = (xs[i][:,4], xs[i][:,5], xs[i][:,6])
#         end
#     end
# end



#
