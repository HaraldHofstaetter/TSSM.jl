using HDF5, JLD
using PyPlot

function get_wf(filename, dset_real="psi_real", dset_imag="psi_imag")
    fid = h5open(filename, "r")
    xmin = read(attrs(fid)["xmin"])
    xmax = read(attrs(fid)["xmax"])
    nx = read(attrs(fid)["nx"])
    ix0 = read(attrs(fid)["ix0"])
    sx = size(fid[dset_real])[1]    
    dx = (xmax-xmin)/nx
    x_grid = linspace(xmin+dx*ix0,xmin+dx*(ix0+sx-1),sx)
    if exists(fid,dset_imag)
        psi = read(fid[dset_real])+1im* read(fid[dset_imag])
    else    
        psi = read(fid[dset_real])
    end     
    close(fid)
    return x_grid, psi
end


function plot_wf(filename)
    x_grid, psi = get_wf(filename)
    plot(x_grid, abs(psi).^2)
end

function plot_potential(filename)
    x_grid, pot = get_wf(filename, "potential", "dummy_imag")
    plot(x_grid, pot)
end


