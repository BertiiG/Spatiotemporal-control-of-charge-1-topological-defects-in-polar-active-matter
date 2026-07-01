#!/usr/bin/perl
use strict;
use warnings;
use Time::Piece;
use File::Path qw(make_path);
use File::Path qw(remove_tree);

# ------------------------
# Fixed parameters
# ------------------------

my $cores = 4;
my $Gd = 50;
my $maxtime = 600;
my $timesteps = 30000;
my $wr_at = 100;
my $Verreps = 7e-2;
my $timetol = 7e-3;
my $timereltol = 7e-3;
my $boxsize = 50;
my $Ks = 5.;
my $nu = 2.;
my $gamma = 0.5;
my $setpoint = 8;
my $kp_dist = 0.0005;
my $ki_dist = 0.005;

# ------------------------
# Parameter sweeps
# ------------------------

my @psi0_vals   = (0.25);
my @ri_vals     = (0.32);
my @zeta_vals   = (-1);
my @lambda_vals = (15);
my @act_vals    = (-100, -90, -80, -70, -60, -50, -40);
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
* scalar(@act_vals)
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
for my $act (@act_vals) {
for my $dK (@dK_vals) {

$count++;

my $psi_str = format_param($psi0);
my $ri_str = format_param($ri);
my $zeta_str = format_param($zeta);
my $lambda_str = format_param($lambda);
my $act_str = format_param($act);
my $dK_str = format_param($dK);
my $setpoint_str = format_param($setpoint);

my $dirname =
"psib_${psi_str}_ri_${ri_str}_zeta_${zeta_str}_lambda_${lambda_str}_act_${act_str}_dK_${dK_str}_setpoint_${setpoint_str}";

# calculate delmu from old non-dimensionalized activity and skip invalid parameters
my $denom = ($zeta + $nu * $gamma * $lambda);
if (abs($denom) < 1e-12) {
    logmsg("[$count/$total] SKIPPING invalid parameters: zeta=$zeta lambda=$lambda (denominator ~0)");
    next;
}
my $delmu = ($act * $Ks)/(($Gd/2.)*($Gd/2.) * $denom);

# ------------------------
# Skip finished simulations
# ------------------------

if (-e "$dirname/outputfile") {

    print "[$count/$total] Simulation $dirname already exists. Redo? (y/n) [auto-yes in 3 min]: ";

    my $answer;
    eval {
        local $SIG{ALRM} = sub { die "timeout\n" };
        alarm 180;  # 3 minutes

        chomp($answer = <STDIN>);
        alarm 0;    # cancel alarm if input received
    };

    if ($@) {
        if ($@ eq "timeout\n") {
            print "\nNo input after 5 minutes → proceeding with redo.\n";
            $answer = 'y';
        } else {
            die "Unexpected error: $@";
        }
    }

    if (!defined $answer || lc($answer) ne 'y') {
        logmsg("Skipping $dirname");
        next;
    }

    logmsg("Cleaning directory $dirname before rerun");
    remove_tree($dirname);   # deletes everything inside

    logmsg("[$count/$total] REDOING completed simulation $dirname (everything will be overwritten)");
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
my \$delmu = $delmu;
my \$setpoint = $setpoint;
my \$kp_dist = $kp_dist;
my \$ki_dist = $ki_dist;

system("mpirun -np \${cores} ./../non_homogeneous_activity_control_rd \${Gd} \${maxtime} \${timesteps} \${wr_at} \${dkK} \${psi0} \${Verreps} \${timetol} \${timereltol} \${boxsize} \${zeta} \${ri} \${lambda} \${delmu} \$setpoint \$kp_dist \$ki_dist -- -pc_type gasm > outputfile");
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

}}}}}}

logmsg("All simulations finished.");

close($LOG);
