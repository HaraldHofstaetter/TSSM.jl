cd(dirname(@__FILE__))
if (!ispath("fftw-3.3.4"))
    download("http://www.fftw.org/fftw-3.3.4.tar.gz","fftw3.tgz")
    run(`tar xzf fftw3.tgz`)
end
cd("fftw-3.3.4")
install_dir = pwd()
run(`./configure --prefix=$install_dir --enable-quad-precision --enable-openmp --enable-threads --enable-shared`)
run(`make install`)
cd("..")

if (!ispath("TSSM"))
    run(`git clone git://github.com/HaraldHofstaetter/TSSM.git`)
else
    cd("TSSM")
    run(`git pull`)
    cd("..")
end

if searchindex(readall(`uname -a`), "juliabox")>0
    run(`make -C TSSM/OPENMP JULIABOX_BUILD=1 all`)
else    
    run(`make -C TSSM/OPENMP all`)
end

if (!ispath("usr"))
    run(`mkdir usr`)
end
if (!ispath("usr/lib"))
    run(`mkdir usr/lib`)
end

run(`mv TSSM/OPENMP/libtssm.$(Libdl.dlext) usr/lib`)

