#!/usr/bin/env perl

#--------------------------------------------------------------
# Francesco Santanastasio  <francesco.santanastasio@cern.ch>
#--------------------------------------------------------------

print "Starting...\n";

use Time::Local;
use Getopt::Std;


## input info

my $storageDir;
my $tagname;
my $inputList;
my $templateCrab;
my $myCMSSWconfig;
my $OutputUserName;
my $publishName;
#my $OutputUserStorageDir;
my $runNumber = -999; #optional
my $jsonFile;         #optional

getopts('h:d:v:i:t:c:n:p:u:r:j:');

if(!$opt_d) {help();}
if(!$opt_v) {help();}
if(!$opt_i) {help();}
if(!$opt_t) {help();}
if(!$opt_c) {help();}
if(!$opt_n) {help();}
if(!$opt_p) {help();}
#if(!$opt_u) {help();}
#if(!$opt_r) {help();}
#if(!$opt_j) {help();}

if($opt_h) {help();}
if($opt_d) {$storageDir = $opt_d;}
if($opt_v) {$tagname = $opt_v;}
if($opt_i) {$inputList = $opt_i;}
if($opt_t) {$templateCrab = $opt_t;}
if($opt_c) {$myCMSSWconfig = $opt_c;}
if($opt_n) {$OutputUserName = $opt_n;}
if($opt_p) {$publishName = $opt_p;}
#if($opt_u) {$OutputUserStorageDir = $opt_u;}
if($opt_r) {$runNumber = $opt_r;}
if($opt_j) {$jsonFile = $opt_j;}

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
$year = $year + 1900;
$mon = $mon + 1;

if($year<10){ $year = "0$year";
	  #   print "change year --> $year\n";
	  }
if($mon<10){ $mon = "0$mon";
	  #   print "change mon --> $mon\n";
	  }
if($mday<10){ $mday = "0$mday";
	  #   print "change mday --> $mday\n";
	  }
if($hour<10){ $hour = "0$hour";
	  #   print "change hour --> $hour\n";
	  }
if($min<10){ $min = "0$min";
	  #   print "change min --> $min\n";
	  }
if($sec<10){ $sec = "0$sec";
	  #   print "change sec --> $sec\n";
	  }

my $date = "$year$mon$mday\_$hour$min$sec";

## create directories

#prepare OUTPUT directories

#-- first letter of username -- this is not needed for OUTPUT (HS)
#$FirstChar = $OutputUserName;
#$FirstChar = substr($FirstChar, 0, 1);
#$FirstChar =~ tr/[A-Z]/[a-z]/;
#print $FirstChar;

###################################################################

## THIS PART MIGHT CHANGE ACCORDINGLY WITH THE OUTPUTDIR ON SE ##

# $OutputPath = "/pnfs/roma1.infn.it/data/cms/store/";
# $OutputUserPath1 = "/user/".$OutputUserName;
# #$OutputUserPath1 = "/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/".$OutputUserName;
# #$OutputUserPath1 = "/cmst3/group/exovv/".$OutputUserName;

$OutputTagName = $tagname."\_".$date;
# $OutputUserPath2 = $OutputUserPath1."/".$OutputTagName;
# $fullOutputPath = $OutputPath.$OutputUserPath2;

###################################################################

##create OUTPUT dir
#print("mkdir -p $fullOutputPath \n");
#system("mkdir -p $fullOutputPath");

# #setting permission on OUTPUT dir
# $currentPath = $OutputPath.$OutputUserPath1;
# print("chmod -r 755 $currentPath \n");
# system("chmod -r 755 $currentPath");
# #my @directories = split(/\//, $OutputUserStorageDir);
# my @directories = split(/\//, "");
# foreach $dir(@directories)
# {
#     $currentPath = $currentPath."/".$dir;
#     print("chmod 755 $currentPath \n");
#     system("chmod 755 $currentPath");
# }
# print("chmod 755 $fullOutputPath \n");
# system("chmod 755 $fullOutputPath");


#other directories
my $productionDir = $storageDir."\/".$tagname."\_".$date;
my $cfgfilesDir = $productionDir."\/"."cfgfiles";
my $outputDir = $productionDir."\/"."output";
my $workDir = $productionDir."\/"."workdir";

print("mkdir -p $productionDir \n");
print("mkdir -p $cfgfilesDir \n");
print("mkdir -p $outputDir \n");
print("mkdir -p $workDir \n");

system("mkdir -p $productionDir");
system("mkdir -p $cfgfilesDir");
system("mkdir -p $outputDir");
system("mkdir -p $workDir");

system("cp $inputList $productionDir\/inputList.txt");

## read input list

open (INPUTLIST, "<$inputList") || die ("...error opening file $inputList $!");
@inputListFile = <INPUTLIST>;
#print @inputListFile;
close(INPUTLIST);


## loop over datasets in the list

foreach $inputListLine(@inputListFile)
{
    chomp($inputListLine);
    #print $inputListLine;

    ## split each line
    my ($dataset, $Nevents, $Njobs) = split(/\s+/, $inputListLine);
    my @datasetParts = split(/\//, $dataset);
    shift @datasetParts; #remove the first element of the list which is an empty-space

    #print "$dataset\n";
    #print "@datasetParts\n";
    #print "$datasetParts[1]\n";
    #print "$Nevents\n";
    #print "$Njobs\n";
    #print "\n";

    print "\n";
    print "processing dataset : $dataset ... \n";


    ## create datasetname

    my $datasetName="";
    my $datasetNameWithRunNumber="";
    my $counter=0;
    foreach $name(@datasetParts)
    {
	$counter++;
	if( $counter < scalar(@datasetParts) ) {$datasetName=$datasetName.$name."__";}
	else {$datasetName=$datasetName.$name;}

    }
    $counter=0;

    ## create workdir for this dataset
    my $thisWorkDir=$workDir."/".$datasetName;
    print "mkdir -p $thisWorkDir\n";
    system "mkdir -p $thisWorkDir";

    ## outputfile .root
    my $outputfile = $datasetName.""."\.root";
#     if($runNumber != -999 )
#     {
# 	$outputfile = $datasetName."\_\_\_\_run".$runNumber.""."\.root";
#     }
    #print "outputfilename : $outputfile ... \n";

    ## read template CMSSW config file
    open (TEMPLATECMSSW, "<$myCMSSWconfig") || die ("...error opening file $myCMSSWconfig $!");
    @templateCMSSWFile = <TEMPLATECMSSW>;
    #print @templateCMSSWFile;
    close(TEMPLATECMSSW);


    ## create new CMSSW config file
    my $newCMSSWconfig = $cfgfilesDir."/".$datasetName."\_"."cmssw"."\_cfg"."\.py";
    print "creating $newCMSSWconfig ... \n";

    open(NEWCMSSWCONFIG,">$newCMSSWconfig");

    foreach $templateCMSSWFileLine(@templateCMSSWFile)
    {

	chomp ($templateCMSSWFileLine);
	#print("$templateCMSSWFileLine\n");

	if($templateCMSSWFileLine=~/THISROOTFILE/)
	{
	    #$templateCMSSWFileLine = "    fileName = cms.string\(\'$outputfile\'\)"; #old version..

	    $substituteString = "\'$outputfile\'";
  	    $templateCMSSWFileLine =~ s/THISROOTFILE/$substituteString/g;

	    #print("$templateCMSSWFileLine\n");
	}
	print NEWCMSSWCONFIG "$templateCMSSWFileLine\n";

    }

    close(NEWCMSSWCONFIG);


    ## read template crab config file
    open (TEMPLATECRAB, "<$templateCrab") || die ("...error opening file $templateCrab $!");
    @templateCrabFile = <TEMPLATECRAB>;
    #print @templateCrabFile;
    close(TEMPLATECRAB);


    ## create new crab config file
    my $newcrabconfig = $cfgfilesDir."/".$datasetName."\_"."crab"."\.cfg";
    print "creating $newcrabconfig ... \n";

    open(NEWCRABCONFIG,">$newcrabconfig");

    foreach $templateCrabFileLine(@templateCrabFile)
    {

	chomp ($templateCrabFileLine);
	#print("$templateCrabFileLine\n");

	if($templateCrabFileLine=~/THISDATASET/)
	{
	    $templateCrabFileLine = "datasetpath = $dataset";
	    #print("$templateCrabFileLine\n");
	}

	if($templateCrabFileLine=~/THISCMSSWCONFIGFILE/)
	{
	    $templateCrabFileLine = "pset = $newCMSSWconfig";
	    #print("$templateCrabFileLine\n");
	}

	if($templateCrabFileLine=~/THISNEVENTS/)
	{
	    $templateCrabFileLine = "total_number_of_events = $Nevents";
	    #print("$templateCrabFileLine\n");
	}

	if($templateCrabFileLine=~/THISNLUMIS/)
	{
	    $templateCrabFileLine = "total_number_of_lumis = $Nevents";
	    #print("$templateCrabFileLine\n");
	}

	if($templateCrabFileLine=~/THISNJOBS/)
	{
	    $templateCrabFileLine = "number_of_jobs = $Njobs";
	    #print("$templateCrabFileLine\n");
	}

	if($templateCrabFileLine=~/THISOUTPUTFILE/)
	{
	    $templateCrabFileLine = "output_file = $outputfile";
	    #print("$templateCrabFileLine\n");
	}

	if($templateCrabFileLine=~/THISUIWORKINGDIR/)
	{
	    $templateCrabFileLine = "ui_working_dir = $thisWorkDir";
	    #print("$templateCrabFileLine\n");
	}

#	if($templateCrabFileLine=~/THISOUTPUTDIR/)
#	{
#	    $templateCrabFileLine = "outputdir = $outputDir";
#	    #print("$templateCrabFileLine\n");
#	}

	if($templateCrabFileLine=~/THISUSERREMOTEDIR/)
	{
	    #$templateCrabFileLine = "user_remote_dir = $OutputUserPath2";
	    $templateCrabFileLine = "user_remote_dir = $OutputTagName";
	    #print("$templateCrabFileLine\n");
	}

	if($templateCrabFileLine=~/THISPUBLISHNAME/)
	{

	    $templateCrabFileLine = "publish_data_name = $publishName";
	    #print("$templateCrabFileLine\n");
	}

	if($templateCrabFileLine=~/THISRUNNUMBER/)
	{
	    if($runNumber != -999 ){
		$templateCrabFileLine = "runselection = $runNumber";
		#print("$templateCrabFileLine\n");
	    }
	    else
	    {
		$templateCrabFileLine = "  "
	    }
	}

	if($templateCrabFileLine=~/THISJSONFILE/)
	{
	    $templateCrabFileLine = "lumi_mask = $jsonFile";
	    #print("$templateCrabFileLine\n");
	}

	print NEWCRABCONFIG "$templateCrabFileLine\n";

    }


    close(NEWCRABCONFIG);


    ## create crab jobs for this dataset
    print "creating jobs for dataset $dataset ... \n";

    print "crab -create -cfg $newcrabconfig\n";
    system "crab -create -cfg $newcrabconfig";

    ## submit crab jobs for this dataset
    #print "submitting jobs for dataset $dataset ... \n";
    #print "crab -submit -c $thisWorkDir\n";
    #system "crab -submit -c $thisWorkDir";

}


#---------------------------------------------------------#

sub help(){
    print "Usage: ./createJobsWithCrab_CopyToT2.pl -d <storageDir> -v <tagname> -i <inputList> -t <templateCrab> -c <myCMSSWconfig> -n <UserName> -p <publishName> [-r <runNumber> optional] [-j <jsonFile> mandatory in JSON MODE] [-h <help?>] \n";
    print "See https://twiki.cern.ch/twiki/bin/edit/CMS/ExoDijet13TeVSubmitJobsWithCrab2014 for more details \n";
    print "Example1 (NORMAL MODE): ./createJobsWithCrab_CopyToT2.pl -d `pwd`/RootNtuple -v V00-03-07-DATA-xxxxxx-yyyyyy -i inputList.txt -t template_crab_grid.cfg -c rootTupleMaker_CRAB_DATA_2012_53X_cfg.py -n <UserName> -p <publishName> \n";
    print "Example2 (NORMAL MODE + RUN SELECTION): ./createJobsWithCrab_CopyToT2.pl -d `pwd`/RootNtuple -v V00-03-07-DATA-xxxxxx-yyyyyy -i inputList.txt -t template_crab_grid.cfg -c rootTupleMaker_CRAB_DATA_2012_53X_cfg.py -n <UserName> -p <publishName> -r 132440 \n";
    print "Example3 (JSON MODE): ./createJobsWithCrab_CopyToT2.pl -d `pwd`/RootNtuple -v V00-03-07-DATA-xxxxxx-yyyyyy -i inputList.txt -t template_crab_grid_JSON.cfg -c rootTupleMaker_CRAB_DATA_2012_53X_cfg.py -n <UserName> -p <publishName> -j JSON.txt \n";

    print "Options:\n";
    print "-d <storageDir>:           choose the local storage directory\n";
    print "-v <tagname>:              choose the tagname of RootNtupleMaker\n";
    print "-i <inputList>:            choose the input list with the datasets\n";
    print "-t <templateCrab>:         choose the crab template\n";
    print "-c <myCMSSWconfig>:        choose the CMSSW config file\n";
    print "-n <UserName>:          choose the user name\n";
    print "-p <publishName>:          choose the name for publication on DBS/DAS (shorter than tagname above since crab will add anyway a long identifier based on the datasetname to the published dataset) \n";
    print "-r <runNumber>:            choose the run number to process (optional, use only in NORMAL MODE)\n";
    print "-j <jsonFile>:             choose the json file (*mandatory* in JSON MODE, don't use in NORMAL MODE)\n";
    print "-h <yes>:                  to print the help\n";
    die "please, try again...\n";
}
