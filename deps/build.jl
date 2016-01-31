cd(dirname(@__FILE__))
run(`make`)
if (!ispath("usr"))
    run(`mkdir usr`)
end
if (!ispath("usr/lib"))
    run(`mkdir usr/lib`)
end

run(`mv TSSM/OPENMP/libtssm.$(Libdl.dlext) usr/lib`)

