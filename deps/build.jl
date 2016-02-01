cd(dirname(@__FILE__))
if (!ispath("TSSM"))
    run(`git clone git://github.com/HaraldHofstaetter/TSSM.git`)
else
    cd("TSSM")
    run(`git pull`)
    cd("..")
end

run(`make -C TSSM/OPENMP JULIA_BUILD=1 all`)

if (!ispath("usr"))
    run(`mkdir usr`)
end
if (!ispath("usr/lib"))
    run(`mkdir usr/lib`)
end

run(`mv TSSM/OPENMP/libtssm.$(Libdl.dlext) usr/lib`)

