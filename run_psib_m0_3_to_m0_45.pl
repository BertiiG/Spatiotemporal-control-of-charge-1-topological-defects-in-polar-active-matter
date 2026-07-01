#!/usr/bin/perl
use strict;
use warnings;
use Time::Piece;
use File::Path qw(make_path);

# ------------------------
# Fixed parameters
# ------------------------

my $cores = 4;
my $Gd = 50;
my $maxtime = 200;
my $timesteps = 30000;
my $wr_at = 100;
my $Verreps = 7e-2;
my $timetol = 7e-3;
my $timereltol = 7e-3;
my $boxsize = 50;

# ------------------------
# Parameter sweeps
# ------------------------

my @psi0_vals   = (-0.3, -0.35, -0.4, -0.45);
my @ri_vals     = (0.32);
my @zeta_vals   = (-1);
my @lambda_vals = (15);
my @dK_vals     = (0);

# ------------------------
# Logging
# ------------------------

# my $script_name = "run_psib_m0_3_to_m0_45";
(my $script_name = $0) =~ s!.*/|\.pl$!!g;
my $logfile = "simulations_log_${script_name}.txt";

open(my $LOG, ">>", $logfile) or die "Cannot open logfile";

sub logmsg {
    my ($msg) = @_;
    my $t = localtime->datetime;
    print "$t  $msg\n";
    print $LOG "$t  $msg\n";
}

# ------------------------
# Convert number to folder-safe format
# ------------------------

sub format_param {

    my ($val) = @_;

    # round nicely
    $val = sprintf("%.4g", $val);

    if ($val < 0) {
        $val = abs($val);
        return "m$val";
    }

    return $val;
}

# ------------------------
# Compute folder size
# ------------------------

sub folder_size {

    my ($dir) = @_;

    my $size = `du -sh $dir 2>/dev/null`;
    chomp($size);

    $size =~ s/\s.*//;

    return $size;
}

# ------------------------
# Total simulation count
# ------------------------

my $total =
  scalar(@psi0_vals)
* scalar(@ri_vals)
* scalar(@zeta_vals)
* scalar(@lambda_vals)
* scalar(@dK_vals);

my $count = 0;

logmsg("Starting simulation sweep: $total total simulations");

# ------------------------
# Main loop
# ------------------------

for my $psi0 (@psi0_vals) {
for my $ri (@ri_vals) {
for my $zeta (@zeta_vals) {
for my $lambda (@lambda_vals) {
for my $dK (@dK_vals) {

$count++;

my $psi_str = format_param($psi0);
my $ri_str = format_param($ri);
my $zeta_str = format_param($zeta);
my $lambda_str = format_param($lambda);
my $dK_str = format_param($dK);

my $dirname =
"psib_${psi_str}_ri_${ri_str}_zeta_${zeta_str}_lambda_${lambda_str}_dK_${dK_str}";

# ------------------------
# Skip finished simulations
# ------------------------

if (-e "$dirname/outputfile") {

logmsg("[$count/$total] SKIPPING completed simulation $dirname");
next;

}

logmsg("[$count/$total] Preparing simulation $dirname");

make_path($dirname) unless -d $dirname;

# ------------------------
# Write run.pl
# ------------------------

my $runfile = "$dirname/run.pl";

open(my $RUN, ">", $runfile) or die "Cannot write run.pl";

print $RUN <<"EOF";
#!/usr/bin/perl -w

my \$cores = $cores;
my \$Gd = $Gd;
my \$maxtime = $maxtime;
my \$timesteps = $timesteps;
my \$wr_at = $wr_at;
my \$dkK = $dK;
my \$psi0 = $psi0;
my \$Verreps = $Verreps;
my \$timetol = $timetol;
my \$timereltol = $timereltol;
my \$boxsize = $boxsize;
my \$zeta = $zeta;
my \$ri = $ri;
my \$lambda = $lambda;

system("mpirun -np \${cores} ./../non_homogeneous_activity \${Gd} \${maxtime} \${timesteps} \${wr_at} \${dkK} \${psi0} \${Verreps} \${timetol} \${timereltol} \${boxsize} \${zeta} \${ri} \${lambda} -- -pc_type gasm > outputfile");
EOF

close($RUN);

chmod 0755, $runfile;

# ------------------------
# Run simulation
# ------------------------

my $start = time();

logmsg("[$count/$total] START simulation $dirname");

system("cd $dirname && ./run.pl");

my $end = time();
my $runtime = $end - $start;

my $size = folder_size($dirname);

logmsg("[$count/$total] FINISHED simulation $dirname  runtime=${runtime}s  folder_size=$size");

}}}}}

logmsg("All simulations finished.");

close($LOG);
