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
my $configfile = YAML::Tiny::LoadFile( $localdir . "/config.yml");

my $DAQ = $configfile->{DAQ};

my $DAQconfigFile = $localdir . "/" . $configfile->{DAQconfig};

my $Emin = $configfile->{Energy}[0];
my $Emax = $configfile->{Energy}[1];

my $TP = $configfile->{Pulser};

my $binning = $configfile->{DAQBin};
#my $suffix = $configfile->{DirSuffix};

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
if($run==119) {$tRun = "new-head";}

# save current umask
my $old_umask = umask;
umask 0000;


my $thisCalibDir;
if ( $DAQ == 1 ) {$thisCalibDir = $localdir . "/FlashCam";}
if ( $DAQ == 0 ) {$thisCalibDir = $localdir . "/Struck";}
mkdir "$thisCalibDir", 0770 unless -d "$thisCalibDir";

$thisCalibDir = $thisCalibDir . "/" . $tRun;
mkdir "$thisCalibDir", 0770 unless -d "$thisCalibDir";

$thisCalibDir = $thisCalibDir . "/calib" . $date;# . $suffix;
mkdir "$thisCalibDir", 0770 unless -d "$thisCalibDir";

#$calibDir = $calibDir . "/" .$tRun;

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
    
    my $sortListCmd = "sort " . $filelistTmp . " > " . $filelist;
    system($sortListCmd);
    unlink $filelistTmp or warn "Could not unlink $filelistTmp: $!";
}
else {
    $filelist = $localdir . "/" . $extFileList;
}
print "filelist: " . $filelist . "\n";
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
		my $filterFile = $Taudir . "/ZACfilter_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".txt";
		my $filterout = $Taudir . "/ZACfilter_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".out";
		my $filterCommand = $localdir . "/ZACfilter " . $Taudir . " " . $binning . " " . $filterLength . " " . $currSigma . " " . $currFlatTop . " " . $currTau . " " . $DAQ . " " . $filterFile . " > " . $filterout;

		#############

		my $outputfile = $Taudir . "/CUSPopt_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".tier2.root";
		my $outputlogfile = $Taudir . "/CUSPopt_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".log";

		my $initobecopied;
		if ( $DAQ==1) {$initobecopied = $localdir . "/template.ini";}
		else {$initobecopied = $localdir . "/template_GERDA.ini";}
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
		#my $jobName = $date . "_" . $tRun . "_" . $pNumb;
		my $jobName = "s" . $date . $pNumb;
		if ($DAQ==1) {$jobName = "f" . $date . $pNumb; }
		my $command = "execModuleIni " . $currini . " >& " . $scriptlog;
		
		my $calibProgram = $localdir . "/energyCalibration";
		my $rootout = $Taudir . "/CUSPopt_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".analysis.out";
		
		my $currFileList = $Taudir . "/GERDA_L". $filterLength . "_sigma" . $currSigma . "_FT" . $currFlatTop . "_tau" . $currTau . ".list";
		open OUT, ">$currFileList" or die "Can't write on file $currFileList: $!\n";
		print OUT $outputfile . "\n";
		close OUT;
		
		my $cmd_analysis = $calibProgram . " -D " . $Taudir . " -t " . $currFileList . " -o -f 0 -m 1 -c " . $DAQconfigFile . " -F 1 -r " . $run . " -T " . $TP . " -z " . $currSigma . "," . $currFlatTop . "," . $currTau . "," . $filterLength . " > " . $rootout;
		
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

		my $QUEUEcmd = "qsub -P short -N " . $jobName . " -e localhost:". $scripterr . " -o localhost:" . $scriptlog . " " . $currscript;
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
my $chi2Limit = 1;
my $firstChn = 0;
if ($DAQ==1) { $firstChn = 24;}
my $totalChannels = 41 + $firstChn;
if ($run==119) {$totalChannels = 12 + $firstChn;}

open my $file, '<', $filelist;
my $firstLine = <$file>;
close $file;
print $firstLine . "\n";
my $position = index($firstLine, $date);
print "position -> " . $position . "\n";
my $time = substr($firstLine, $position+9, 6);
print $time . "\n";


#my $countJob = "qstat -u " . $myUser . " |  grep " . $date . "_" . $tRun . " | wc -l ";
my $countJob = "qstat -u " . $myUser . " |  grep s" . $date . " | wc -l ";
if ($DAQ==1) { $countJob = "qstat -u " . $myUser . " |  grep f" . $date . " | wc -l ";}
print $countJob . "\n";
my $inQueue = 1;

my $scandir = $thisCalibDir;
while ($inQueue) {

    if($skip==0) {
	sleep 30;
    }
    my $actualJob = `$countJob`;
    print $actualJob;
    if($actualJob==0) {
	$inQueue=0;

	#### Create the results directory
	$scandir = $scandir . "/scanResults/";
	mkdir "$scandir", 0770 unless -d "$scandir";

	#### Loop on channels
	
	for(my $channel = $firstChn; $channel < $totalChannels; $channel++)  {
	    #if ( !($run > 68 && $channel == 7) ){
	    my $isInverted = 0;
	    if ( $run>=95 && $channel>=$firstChn+36 ){ $isInverted = 1;}
	    if ( ($run==119) && ($channel==32 || $channel==33) ){ $isInverted = 1;}
	    my $currFile = $scandir . "/ZAC-FWHM_channel" . $channel . ".txt";
	    my $cutCmd = "cat " . $resdir . "/FT_1.0mus/Sigma_*/Tau_*/ZAC-FWHM_chn" . $channel . ".txt > " . $currFile;
	    if ($isInverted) {$cutCmd = "cat " . $resdir . "/FT_1.5mus/Sigma_*/Tau_*/ZAC-FWHM_chn" . $channel . ".txt > ".$currFile;}
	    system($cutCmd);
	    #}
	}

	my $analysisCmd = $localdir . "/scanOptZACGraph " . $scandir . " " . $run . " " . $totalChannels . " " . $firstChn . " " . $date  . " " . $time  . " " . " " . $plotType . " " . $savePlots  . " " . $plot . " " . $chi2Limit . " >& " .$scandir . "/scanOptZACGraph.out";
	print $analysisCmd . "\n";
	system($analysisCmd);

    }
}

# Reset umask to old value
umask $old_umask;

my $from = 'optimizeZAC';
my $to = 'valerio.dandrea\@lngs.infn.it';
my $subject = $tRun . "-" . $date . "T" . $time . "Z";
my $jsonfile = $scandir . "gerda-" . $tRun . "-" . $date . "T" . $time . "Z-cal-ged-tier2-calib.json";
my $resultsfile = $scandir . "/ZACanalysis.txt";
my $body = "Hi, \n this is an automatic message from ZAC filter optimization of calibration:\n". $tRun . " " . $date . " " . $time . " \n\n";
$body = $body . "In attachment you find the json file for the tier2 production. \n\n";
$body = $body . "If this file is corrupted you can find it on the mpi-hd.mpg.de cluster: \n" . $jsonfile . "\n\n";
my $sendMail = "echo \"". $body . "\" | mailx -s \"". $subject . "\" -r " . $from . " -a \"". $resultsfile . "\" -a \"". $jsonfile . "\" " .  $to;
system($sendMail);
