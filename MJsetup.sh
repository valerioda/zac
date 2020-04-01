export PERL_BASE="/lfs/l1/legend/users/dandrea/software/perl"
export PERL5LIB="$PERL_BASE/perl-5.24.1/lib/perl5:$PERL_BASE/perl-5.24.1/lib/perl5/x86_64-linux"
export PATH="$PERL_BASE/perl-5.24.1/bin:$PATH"

export PATH="/lfs/l1/legend/users/dandrea/software/swmod/bin:$PATH"
. swmod.sh init
export SWMOD_INST_BASE=/lfs/l1/gerda/kermaidy/Analysis/software
export SWMOD_MODPATH=/lfs/l1/gerda/kermaidy/Analysis/software
export SWMOD_HOSTSPEC=linux-scientific-7.2-x86_64
export GERDAINSTALL=$GERDA_SOFTWARE/gerda/linux-scientific-7.2-x86_64/master
swmod load gerda@master

