cd(dirname(@__FILE__))

if (!ispath("usr"))
    run(`mkdir usr`)
end
if (!ispath("usr/lib"))
    run(`mkdir usr/lib`)
end

if (   !ispath("usr/lib/libfftw3q.so.3.4.4")
    || !ispath("usr/lib/libfftw3q_omp.so.3.4.4"))
    if (!ispath("fftw-3.3.4"))
        download("http://www.fftw.org/fftw-3.3.4.tar.gz","fftw3.tgz")
        run(`tar xzf fftw3.tgz`)
        rm("fftw3.tgz")
    end
    install_dir = pwd()
    cd("fftw-3.3.4")
    run(`./configure --prefix=$install_dir/usr --enable-quad-precision --enable-openmp --enable-threads --enable-shared`)
    run(`make install`)
    cd("..")
    rm("fftw-3.3.4", recursive=true)
end

if (!ispath("TSSM"))
    run(`git clone git://github.com/HaraldHofstaetter/TSSM.git`)
else
    cd("TSSM")
    run(`git pull`)
    cd("..")
end

if searchindex(readall(`uname -a`), "juliabox")>0
    run(`make -C TSSM/OPENMP JULIABOX_BUILD=1 all`)
    run(`mv TSSM/OPENMP/libtssm.$(Libdl.dlext) usr/lib`)
    run(`make -C TSSM/OPENMP clean`)
    run(`make -C TSSM/QUADPREC_OPENMP JULIABOX_BUILD=1 all`)
    run(`mv TSSM/QUADPREC_OPENMP/libtssmq.$(Libdl.dlext) usr/lib`)
    run(`make -C TSSM/QUADPREC_OPENMP clean`)
else    
    run(`make -C TSSM/OPENMP all`)
    run(`make -C TSSM/QUADPREC_OPENMP all`)
    run(`mv TSSM/OPENMP/libtssm.$(Libdl.dlext) usr/lib`)
    run(`mv TSSM/QUADPREC_OPENMP/libtssmq.$(Libdl.dlext) usr/lib`)
end


