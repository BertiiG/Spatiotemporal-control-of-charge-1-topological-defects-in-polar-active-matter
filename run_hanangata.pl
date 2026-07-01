#!/usr/bin/perl -w

my $cores = 4; 
my $Gd = 50; #argv1
my $maxtime = 2160; #argv2
my $timesteps = 90000; #argv3
my $wr_at = 100; #argv4
my $dkK = 0; #argv5
my $psi0 = 0.016; #argv6
my $Verreps = 0.07; #argv7
my $timetol = 0.007; #argv8
my $timereltol = 0.007; #argv9
my $boxsize = 50; #argv10
my $zeta = -1; #argv11
my $ri = 0.32; #argv12
my $lambda = 15; #argv13
my $act = -60; 
my $Ks = 5;
my $nu = 2;
my $gamma = 0.5;

# calculate delmu from old non-dimensionalized activity and skip invalid parameters
my $denom = ($zeta + $nu * $gamma * $lambda);
if (abs($denom) < 1e-12) {
    logmsg("SKIPPING invalid parameters: zeta=$zeta lambda=$lambda (denominator ~0)");
    next;
}else{
	my $delmu = ($act * $Ks)/(($Gd/2.)*($Gd/2.) * $denom); #argv14

	my $setpoint = 14; #argv15

	system("mpirun -np ${cores} ./../non_homogeneous_activity_hanangata ${Gd} ${maxtime} ${timesteps} ${wr_at} ${dkK} ${psi0} ${Verreps} ${timetol} ${timereltol} ${boxsize} ${zeta} ${ri} ${lambda} ${delmu} ${setpoint} -- -pc_type gasm > outputfile");
}