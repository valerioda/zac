#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use YAML::Tiny 'LoadFile';
use YAML::Tiny;
use Data::Dumper;
use POSIX qw(ceil floor);

use constant false => 0;
use constant true  => 1;

my $localdir = $ENV{PWD};

my $calibDir = $ARGV[0]; # example: /nfs/gerda5/gerda-data/blind/stable/gen/tier1/ged/cal/
my $run = $ARGV[1];
my $date = $ARGV[2];
my $len = @ARGV;
#print "lenght- > $len\n";
my $extF = 0;
my $extFileList;
if($len > 3 ) {
    print "Length = $len\n";
    $extF = 1;
    $extFileList = $ARGV[3];
}
my $skip = 0;
if($len > 4 ) {
    print "Length = $len\n";
    $skip = 1;
}

#print "EXTF = $extF \n";
print "Skip = $skip\n";

#Open the configuration file
my $configfile = YAML::Tiny::LoadFile( $localdir . "/GERDA.CUSPopt.config.yml");

my $DAQconfigFile = $localdir . "/" . $configfile->{DAQconfig};

my $Emin = $configfile->{Energy}[0];
my $Emax = $configfile->{Energy}[1];

my $TP = $configfile->{Pulser};

my $binning = $configfile->{DAQBin};
my $suffix = $configfile->{DirSuffix};

my $filterLength = $configfile->{FilterLength};

my $startFlatTop = $configfile->{FlatTop}[0];
my $stopFlatTop = $configfile->{FlatTop}[1];
my $stepFlatTop = $configfile->{FlatTop}[2];

my $startSigma = $configfile->{Sigma}[0];
my $stopSigma = $configfile->{Sigma}[1];
my $stepSigma = $configfile->{Sigma}[2];

my $startTau = $configfile->{Tau}[0];
my $stopTau = $configfile->{Tau}[1];
my $stepTau = $configfile->{Tau}[2];


my @expr = ( "!filelist!", "!rootout!",
	     "!logfile!", "!emin!", "!emax!", "!filterfile!");

my @cmd = ( "!filtercommand!", "!mvcommand!", "!gelatiocommand!",  "!cmdanalysis!", "!chmodcmd!");

my $myUser = $ENV{'USER'};

my $tRun = "run";
if($run<100) {
    $tRun = $tRun . "00" . $run;
}
elsif($run<1000) {
    $tRun = $tRun . "0" . $run;
}
else {
    $tRun = $tRun . "" . $run;
}

# save current umask
my $old_umask = umask;
umask 0000;

my $thisCalibDir = $localdir . "/" . $tRun;
mkdir "$thisCalibDir", 0770 unless -d "$thisCalibDir";

$thisCalibDir = $thisCalibDir . "/calib" . $date . $suffix;
mkdir "$thisCalibDir", 0770 unless -d "$thisCalibDir";

$calibDir = $calibDir . "/" .$tRun;

my $filelist = $thisCalibDir . "/GERDA_calib_" . $tRun . "_". $date . ".list";
my $filelistTmp = $thisCalibDir . "/tmpGERDA_calib_" . $tRun . "_". $date . ".list";
#print $filelistTmp . "\n";
#print $filelist . "\n";
if($extF==0) {
    print "HERE\n";
    open OUT, ">$filelistTmp" or die "Can't write on file $filelist: $!\n";
    opendir(DIR, $calibDir) or die $!;


    while (my $file = readdir(DIR)) {
	if($file =~ /$date/ && $file =~ /\.root/) {
	    print OUT $calibDir . "/". $file . "\n";
	    print $file . "\n";
	}
    }
    close OUT;
    closedir(DIR);
}
else {
    $filelistTmp = $thisCalibDir . "/tmp" . $extFileList;
    $filelist = $thisCalibDir . "/" . $extFileList;

}
my $sortListCmd = "sort " . $filelistTmp . " > " . $filelist;
system($sortListCmd);
if($extF==0) {
    unlink $filelistTmp or warn "Could not unlink $filelistTmp: $!";
}

#### Create the results directory
my $resdir = $thisCalibDir . "/Results";
mkdir "$resdir", 0770 unless -d "$resdir";

#### Loop on filter parameters

my $pNumb = 0;
my $currFlatTop=$startFlatTop;
if($skip==0) {
    while( $currFlatTop <= $stopFlatTop)
    {
	
	my $FTdir = $resdir . "/FT_" . $currFlatTop . "mus";
	mkdir "$FTdir", 0770 unless -d "$FTdir";

	for(my $currSigma=$startSigma; $currSigma<=$stopSigma; $currSigma+=$stepSigma)
	{

	    my $Sigmadir = $FTdir . "/Sigma_" . $currSigma . "mus";
	    mkdir "$Sigmadir", 0770 unless -d "$Sigmadir";

	    for(my $currTau=$startTau; $currTau<=$stopTau; $currTau+=$stepTau)
	    {

		my $Taudir = $Sigmadir . "/Tau_" . $currTau . "mus";
		mkdir "$Taudir", 0770 unless -d "$Taudir";

		print "Submitting -> Flat Top: " . $currFlatTop . "  - Sigma: " . $currSigma . " - Tau: " . $currTau . "\n";

		#### Running ZACFilter
		my $filterName = "ZACfilter_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".txt";
		my $filterCommand = $localdir . "/ZACfilter " . $binning . " " . $filterLength . " " . $currSigma . " " . $currFlatTop . " " . $currTau . " " . $filterName;

		my $filterFile = $Taudir . "/" . $filterName;

		#############

		my $outputfile = $Taudir . "/CUSPopt_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".tier2.root";
		my $outputlogfile = $Taudir . "/CUSPopt_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".log";

		my $initobecopied = $localdir . "/GERDA.CUSPopt.template.ini";
		open IN, $initobecopied or die "Can't read source file $initobecopied: $!\n";

		my $currini = $Taudir . "/CUSPopt_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".ini";
		open OUT, ">$currini" or die "Can't write on file $currini: $!\n";

		while (<IN>) {
		    s/$expr[0]/$filelist/g;
		    s/$expr[1]/$outputfile/g;
		    s/$expr[2]/$outputlogfile/g;
		    s/$expr[3]/$Emin/g;
		    s/$expr[4]/$Emax/g;
		    s/$expr[5]/$filterFile/g;
		    print OUT $_;
		}

		close IN;
		close OUT;
		
		my $scriptlog = $Taudir . "/CUSPopt_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".out";
		my $scripterr = $Taudir . "/CUSPopt_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".err";
		my $jobName = $date . "_" . $tRun . "_" . $pNumb;
		my $command = "execModuleIni " . $currini . " >& " . $scriptlog;
		
		my $calibProgram = $localdir . "/energyCalibration";
		my $rootout = $Taudir . "/CUSPopt_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".analysis.out";
		
		my $currFileList = $Taudir . "/GERDA_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".list";
		open OUT, ">$currFileList" or die "Can't write on file $currFileList: $!\n";
		print OUT $outputfile . "\n";
		close OUT;
		
		my $cmd_analysis = $calibProgram . " -t " . $currFileList . " -o -f 0 -m 1 -c " . $DAQconfigFile . " -F 1 -r " . $run . " -T " . $TP . " -z " . $currSigma . "," . $currFlatTop . "," . $currTau . "," . $filterLength . " > " . $rootout;
		
		my $scripttobecopied = $localdir . "/script.template.sh";
		open IN, $scripttobecopied or die "Can't read source file $scripttobecopied: $!\n";
		
		my $currscript = $Taudir . "/script.CUSPopt_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".sh";
		open OUT, ">$currscript" or die "Can't write on file $currscript: $!\n";

		
		while (<IN>) {
		    s/$cmd[0]/$filterCommand/g;
		    s/$cmd[2]/$command/g;
		    s/$cmd[3]/$cmd_analysis/g;
		    print OUT $_;
		}

		close IN;
		close OUT;

		### submit job to queue

		my $QUEUEcmd = "qsub -N " . $jobName . " -q gerda -V -d " . $Taudir . " -m abe -e localhost:". $scripterr . " -o localhost:" . $scriptlog . " -l mem=4000mb " . $currscript;
		system($QUEUEcmd);
		print $QUEUEcmd . "\n";
		$pNumb++;
	    }
	}

	$currFlatTop += $stepFlatTop;
    }
}


my $plotType = 0;
my $savePlots = 1;
my $plot = 0;
my $chi2Limit = 100;
my $totalChannels = 41;

open my $file, '<', $filelist;
my $firstLine = <$file>;
close $file;
print $firstLine . "\n";
my $position = index($firstLine, $date);
print "position -> " . $position . "\n";
my $time = substr($firstLine, $position+9, 6);
print $time . "\n";


#my $countJob = "qstat -u " . $myUser . " |  grep " . $userName . " | wc -l ";
my $countJob = "qstat -u " . $myUser . " |  grep " . $date . "_" . $tRun . " | wc -l ";
print $countJob . "\n";
my $inQueue = 1;

my $scandir = $thisCalibDir;
while ($inQueue) {

    if($skip==0) {
	sleep 300;
    }
    my $actualJob = `$countJob`;
    print $actualJob;
    if($actualJob==0) {
	$inQueue=0;

	#### Create the results directory
	$scandir = $scandir . "/scanResults/";
	mkdir "$scandir", 0770 unless -d "$scandir";

	#### Loop on channels
	for(my $channel = 0; $channel < $totalChannels; $channel++)  {
	    #if ( !($run > 68 && $channel == 7) ){
	    my $currFile = $scandir . "/CUSPopt_analysis_channel" . $channel . "_length155.000000.dat";
	    my $cutCmd = "cat " . $resdir . "/FT_*/Sigma_*/Tau_*/CUSPopt_analysis_channel" . $channel . "_length*.dat > " . $currFile;
	    system($cutCmd);
	    #}
	}

	my $analysisCmd = $localdir . "/scanOptZACGraph " . $scandir . " " . $run . " " . $totalChannels . " " . $date  . " " . $time  . " " . " " . $plotType . " " . $savePlots  . " " . $plot . " " . $chi2Limit . " >& " .$scandir . "/scanOptZACGraph.out";

	system($analysisCmd);

    }
}

# Reset umask to old value
umask $old_umask;

my $from = 'optimizeZAC';
my $to = 'francesco.salamida\@aquila.infn.it, valerio.dandrea\@lngs.infn.it';
my $subject = $tRun . "-" . $date . "T" . $time . "Z";
my $jsonfile = $scandir . "gerda-" . $tRun . "-" . $date . "T" . $time . "Z-cal-ged-tier2-calib.json";
my $body = "Hi, \n this is an automatic message from ZAC filter optimization of calibration:\n". $tRun . " " . $date . " " . $time . " \n\n";
$body = $body . "In attachment you find the json file for the tier2 production. \n\n";
$body = $body . "If this file is corrupted you can find it at LNGS: \n" . $jsonfile . "\n\n";
my $sendMail = "echo \"". $body . "\" | mailx -s \"". $subject . "\" -r " . $from . " -a \"". $jsonfile . "\" " .  $to;
system($sendMail);
