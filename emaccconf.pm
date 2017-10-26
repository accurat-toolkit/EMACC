# EMACC configuration file. Change this before running!
#
# (C) ICIA 2011, Radu ION.
#
# ver 1.0, 30.09.2011, Radu ION: portable for Windows and Linux.
# ver 1.1, 10.10.2011, Radu ION: added EMACCMODE.

package emaccconf;

use strict;
use warnings;
use File::Spec;
use File::Path;
use Sys::Hostname;

sub checkDir( $$ );
sub checkFile( $$ );
sub checkClusterFile( $$ );
sub checkInitMode( $$ );
sub checkPLexType( $$ );
sub checkEmaccMode( $$ );
sub checkProbType( $$ );
sub checkInt( $$ );
sub checkProb( $$ );
sub checkReal( $$ );
sub checkBool( $$ );
sub checkLang( $$ );
sub new;
sub addValue( $$$$ );
sub genClusterFile();
sub portableCopyFileToDir( $$ );
sub portableRemoveFile( $ );
sub portableRemoveFileFromDir( $$ );
sub portableRemoveAllFilesFromDir( $ );
sub portableForkAndDetach( $ );

##################
#CONFIG FILE######
##################

#Only change values between 'BEGIN CONF' and 'END CONF'!
sub new {
	my( $classname ) = shift();
	my( $conf ) = shift();
	my( $this ) = {};
	
	############################## BEGIN CONF ##############################################################################
	#Source language
	addValue( $this, $conf, "SRCL", "en" );
	#Target language
	addValue( $this, $conf, "TRGL", "ro" );
	#How many EM loops should EMACC execute
	addValue( $this, $conf, "EMLOOPS", 5 );
	#EMACC mode: 'emacc-full' (the EM algorithm) or 'emacc-simple' (a simple greedy alignment in which all pairs above a certain prob are output).
	#Currently, if too many possible document pairs, the models are too big.
	addValue( $this, $conf, "EMACCMODE", 'emacc-full' );
	#How many document alignments we permit.
	#If MAXTARGETALIGNMENTS = 1, then we have only 1:1 alignments, if = 2, then we permit 1:2 alignments...
	addValue( $this, $conf, "MAXTARGETALIGNMENTS", 3 );
	#What is the translation equivalents probability above which we should re-estimate their probabilities in the EM process?
	addValue( $this, $conf, "TEQPUPDATETHR", 0.4 );
	#What percent from the alignments to output? (usually 1)
	#But for 'emacc-simple' this is usually lower (e.g. 0.5)
	addValue( $this, $conf, "OUTPUTPERCENT", 1 );
	#What percent from the top document alignments found to be used for translation equivalents re-estimation (see the paper from BUCC 2011)
	addValue( $this, $conf, "UPDATEPERCENT", $this->{"OUTPUTPERCENT"} );
	#How small is the smallest value? (best left unchanged)
	addValue( $this, $conf, "DEFSMALLP", 1e-30 );
	#Initial document alignment distribution:
	#- uniform: D1
	#- document similarity precomputed: D2
	#Best results: D2!
	addValue( $this, $conf, "INIDISTRIB", "D2" );
	#This variable holds the name of the corpus we are processing. It will appear in the output model names.
	addValue( $this, $conf, "CORPUSNAME", "emacc2-run" );
	#The name of the output file (the file that will store the document alignments)
	#It is advised to keep that name that contains the values of important parameters.
	addValue( $this, $conf, "DALOUTFILE", $this->{"SRCL"} . "-" . $this->{"TRGL"} . "-" . $this->{"CORPUSNAME"} . "-INID_" . $this->{"INIDISTRIB"} . "-TEQUPD_" . $this->{"TEQPUPDATETHR"} . "-OUTP_" . $this->{"OUTPUTPERCENT"} . ".al" );
	#What is the minimum probability of a translation equivalents pair in order to be used in the alignment process?
	#GIZA++ scores: 0..1
	addValue( $this, $conf, "LEXALSCORE", 0.4 );
	#What measure to use when computing the probability of the initial alignment of two documents:
	#'fair': the measure is given by the total number of translations over all possible pairs;
	#'biased': D2 measure (see the BUCC 2011 paper from ACCURAT publications page).
	addValue( $this, $conf, "PLEXTYPE", 'biased' );
	#If to try and stem words (1) or not (0).
	#Lemmas are tried with suffix matching ... time consuming.
	addValue( $this, $conf, "LEMMAS", 1 );
	#The 'clsuter.info' file which holds the number of available machines/cpus
	addValue( $this, $conf, "CLUSTERFILE", 'generate' );
	#The NFS mount point (the NFS directory which contains the documents do be aligned and the lingustic resources)
	#It can be a local directory if no clustering is involved.
	#This directory must exist (on each cluster node) and must be writable and readable by 'rion' user.
	addValue( $this, $conf, "NFSMOUNTPOINT", "." );
	#This directory must exist (on each cluster node) and must be writable and readable by 'rion' user.
	#An attempt is made to create it if does not exist.
	my( $tmpdir ) = File::Spec->catdir( File::Spec->tmpdir(), "tmpalign" );
        
	mkpath( $tmpdir );
	addValue( $this, $conf, "LOCALMOUNTPOINT", $tmpdir );
	#The pre-computed D(1|2) distribution file name (the '.pc' file)
	addValue( $this, $conf, "PRECOMPMODELFILE", $this->{"SRCL"} . "-" . $this->{"TRGL"} . "-" . $this->{"CORPUSNAME"} . "-" . $this->{"PLEXTYPE"} . "-LEM_" . $this->{"LEMMAS"} . "-TEQ_" . $this->{"LEXALSCORE"} . ".pc" );
	#The pre-computed translation equivalents pairs per document pair (the '.f' file)
	addValue( $this, $conf, "DPAIRTEQMODELFILE", $this->{"SRCL"} . "-" . $this->{"TRGL"} . "-" . $this->{"CORPUSNAME"} . "-LEM_" . $this->{"LEMMAS"} . "-TEQ_" . $this->{"LEXALSCORE"} . ".f" );
	#When using GIZA++ probabilities
	#'giza' - retain the original probability;
	#'comp' - recompute a different probabiliry.
	addValue( $this, $conf, "PROBTYPE", 'giza' );
	#Stop words files:
	addValue( $this, $conf, "STOPWSRCL", File::Spec->catfile( $this->{"NFSMOUNTPOINT"}, "res", "stopwords_" . $this->{"SRCL"} . ".txt" ) );
	addValue( $this, $conf, "STOPWTRGL", File::Spec->catfile( $this->{"NFSMOUNTPOINT"}, "res", "stopwords_" . $this->{"TRGL"} . ".txt" ) );
	#Inflectional endings files:
	addValue( $this, $conf, "INFLSRCL", File::Spec->catfile( $this->{"NFSMOUNTPOINT"}, "res", "endings_" . $this->{"SRCL"} . ".txt" ) );
	addValue( $this, $conf, "INFLTRGL", File::Spec->catfile( $this->{"NFSMOUNTPOINT"}, "res", "endings_" . $this->{"TRGL"} . ".txt" ) );
	#We have en_<lang> dictionaries.
	#If you want <lang>_en alignments, put 1 here (and be sure to set SRCL and TRGL properly!). Else 0.
	addValue( $this, $conf, "DICTINVERSE", 0 ); checkBool( "DICTINVERSE", $this );
	
	my( $DICTFILE );
	my( $DICTINVERSE ) = $this->{"DICTINVERSE"};

	if ( $DICTINVERSE ) {
		$DICTFILE = File::Spec->catfile( $this->{"NFSMOUNTPOINT"}, "dict", $this->{"TRGL"} . "_" . $this->{"SRCL"} );
	}
	else {
		$DICTFILE = File::Spec->catfile( $this->{"NFSMOUNTPOINT"}, "dict", $this->{"SRCL"} . "_" . $this->{"TRGL"} );
	}
	
	addValue( $this, $conf, "DICTFILE", $DICTFILE );
	############################## END CONF ################################################################################
	
	checkLang( "SRCL", $this );
	checkLang( "TRGL", $this );
	checkReal( "DEFSMALLP", $this );
	checkProb( "TEQPUPDATETHR", $this );
	checkProb( "OUTPUTPERCENT", $this );
	checkProb( "LEXALSCORE", $this );
	checkProb( "UPDATEPERCENT", $this );
	checkInt( "EMLOOPS", $this );
	checkInt( "MAXTARGETALIGNMENTS", $this );
	checkInitMode( "INIDISTRIB", $this );
	checkPLexType( "PLEXTYPE", $this );
	checkEmaccMode( "EMACCMODE", $this );
	checkProbType( "PROBTYPE", $this );
	checkBool( "LEMMAS", $this );
	checkClusterFile( "CLUSTERFILE", $this );
	checkDir( "NFSMOUNTPOINT", $this );
	checkFile( "STOPWSRCL", $this );
	checkFile( "STOPWTRGL", $this );
	checkFile( "INFLSRCL", $this );
	checkFile( "INFLTRGL", $this );
	checkFile( "DICTFILE", $this );
	checkDir( "NFSMOUNTPOINT", $this );
	checkDir( "LOCALMOUNTPOINT", $this );

	bless( $this, $classname );
	return $this;
}

#The rest of these functions are not to be called through the object interface.
sub addValue( $$$$ ) {
	my( $this, $conf, $varname, $vardefaultvalue ) = @_;

	if ( exists( $conf->{$varname} ) && $conf->{$varname} ne "" ) {
		$this->{$varname} = $conf->{$varname};
	}
	else {
		$this->{$varname} = $vardefaultvalue;
	}
}

sub checkDir( $$ ) {
	my( $varname, $this ) = @_;
	my( $dir ) = $this->{$varname};
	my( $testfile ) = File::Spec->catfile( $dir, "test2476blah2144" );  
	
	open( TF, ">", $testfile ) or die( "emaccconf::checkDir: '$varname' has issues: '$! (" . int( $! ) . ")'.\n" );
	close( TF );
	
	unlink( $testfile ) or warn( "emaccconf::checkDir: could not remove '$testfile' because '$!'.\n" );
}

sub checkFile( $$ ) {
	my( $varname, $this ) = @_;
	my( $file ) = $this->{$varname};
	
	if ( ! -f ( $file ) ) {
		die( "emaccconf::checkFile: '$varname' does not exist !\n" );
	} 
}

sub checkClusterFile( $$ ) {
	my( $varname, $this ) = @_;

	if ( $this->{$varname} eq "generate" ) {
		genClusterFile();
		$this->{$varname} = "cluster-autogen.info";
	}

	checkFile( $varname, $this );
}

sub checkInitMode( $$ ) {
	my( $varname, $this ) = @_;
	my( $smode ) = $this->{$varname};
	
	if ( $smode ne "D1" && $smode ne "D2" ) {
		die( "emaccconf::checkInitMode: invalid value for '$varname' (either 'D1' or 'D2') !\n" );
	}
}

sub checkPLexType( $$ ) {
	my( $varname, $this ) = @_;
	my( $smode ) = $this->{$varname};
	
	if ( $smode ne "fair" && $smode ne "biased" ) {
		die( "emaccconf::checkPLexType: invalid value for '$varname' (either 'fair' or 'biased') !\n" );
	}
}

sub checkEmaccMode( $$ ) {
	my( $varname, $this ) = @_;
	my( $smode ) = $this->{$varname};
	
	if ( $smode ne "emacc-full" && $smode ne "emacc-simple" ) {
		die( "emaccconf::checkEmaccMode: invalid value for '$varname' (either 'emacc-full' or 'emacc-one') !\n" );
	}
}

sub checkProbType( $$ ) {
	my( $varname, $this ) = @_;
	my( $smode ) = $this->{$varname};
	
	if ( $smode ne "giza" && $smode ne "comp" ) {
		die( "emaccconf::checkProbType: invalid value for '$varname' (either 'giza' or 'comp') !\n" );
	}
}

sub checkInt( $$ ) {
	my( $varname, $this ) = @_;
	my( $int ) = $this->{$varname};
	
	if ( $int !~ /^[0-9]+$/ ) {
		die( "emaccconf::checkInt: invalid value for '$varname' !\n" );
	}
}

sub checkProb( $$ ) {
	my( $varname, $this ) = @_;
	my( $prob ) = $this->{$varname};
	
	if ( $prob !~ /^[0-9]+(?:\.[0-9]+)?(?:[eE]-?[0-9]+)?$/ ) {
		die( "emaccconf::checkProb: invalid value for '$varname' (real number) !\n" );
	}
	
	if ( $prob < 0 || $prob > 1 ) {
		die( "emaccconf::checkProb: invalid value for '$varname' ([0..1]) !\n" );
	}
}

sub checkReal( $$ ) {
	my( $varname, $this ) = @_;
	my( $real ) = $this->{$varname};
	
	if ( $real !~ /^[0-9]+(?:\.[0-9]+)?(?:[eE]-?[0-9]+)?$/ ) {
		die( "emaccconf::checkReal: invalid value for '$varname' !\n" );
	}
}

sub checkBool( $$ ) {
	my( $varname, $this ) = @_;
	my( $bool ) = $this->{$varname};
	
	if ( $bool !~ /^[01]$/ ) {
		die( "emaccconf::checkBool: invalid value for '$varname' (either '0' or '1') !\n" );
	}
}

sub checkLang( $$ ) {
	my( $varname, $this ) = @_;
	my( $lang ) = $this->{$varname};
	
	if ( $lang !~ /^(?:en|ro|de|lt|lv|sl|el|hr|et)$/ ) {
		die( "emaccconf::checkLang: invalid value for '$varname' !\n" );
	}
}

sub genClusterFile() {
	my( $thishostname ) = hostname();

	open( CLF, ">", "cluster-autogen.info" ) or die( "pdataextractconf::genClusterFile: cannot open file 'cluster-autogen.info' !\n" );

	print( CLF "#This is a comment.\n" );
	print( CLF "#This autogenerated file will NOT work if a cluster run is desired!\n" );
	print( CLF "#Line format (tab separated fields):\n" );
	print( CLF "#- hostname of the machine in cluster (run 'hostname' command)\n" );
	print( CLF "#- IP of the machine\n" );
	print( CLF "#- ID (string) of one CPU core\n\n" );

	#Linux systems...
	if ( -f ( "/proc/cpuinfo" ) ) {
		open( CPU, "<", "/proc/cpuinfo" ) or die( "pdataextractconf::genClusterFile: cannot open file '/proc/cpuinfo' !\n" );

		while ( my $line = <CPU> ) {
			$line =~ s/^\s+//;
			$line =~ s/\s+$//;

			next if ( $line !~ /:/ );

			my( $variable, $value ) = split( /\s*:\s*/, $line );

			$variable =~ s/^\s+//;
			$variable =~ s/\s+$//;
			$value =~ s/^\s+//;
			$value =~ s/\s+$//;

			if ( $variable eq "processor" ) {
				print( CLF $thishostname . "\t" . "127.0.0.1" . "\t" . "cpu$value" . "\n" );
			}
		}

		close( CPU );
	}
	#Windows systems...
	else {
		#Don't know. 1 core :D
		print( CLF $thishostname . "\t" . "127.0.0.1" . "\t" . "cpu0" . "\n" );
	}
	
	close( CLF );
}

sub portableCopyFileToDir( $$ ) {
	my( $file, $dir ) = @_;
	
	#Windows run
	if ( $^O =~ /^MSWin(?:32|64)$/i ) {
		warn( "`copy \/Y ${file} ${dir}\\'\n" );
		qx/copy \/Y ${file} ${dir}\\/;
	}
	#Linux run
	elsif ( $^O =~ /^Linux$/i || $^O =~ /^Cygwin$/i || $^O =~ /^MSys$/i ) {
		qx/cp -fv ${file} ${dir}\/ 1>&2/;
	}
	else {
		die( "emaccconf::portableCopyFileToDir: unsupported operating system '$^O' !\n" );
	}
}

sub portableRemoveFile( $ ) {
	my( $file ) = $_[0];
	
	#Windows run
	if ( $^O =~ /^MSWin(?:32|64)$/i ) {
		warn( "`del \/F \/Q ${file}'\n" );
		qx/del \/F \/Q ${file}/;
	}
	#Linux run
	elsif ( $^O =~ /^Linux$/i || $^O =~ /^Cygwin$/i || $^O =~ /^MSys$/i ) {
		qx/rm -fv ${file} 1>&2/;
	}
	else {
		die( "emaccconf::portableRemoveFile: unsupported operating system '$^O' !\n" );
	}
}

sub portableRemoveFileFromDir( $$ ) {
	my( $dir, $file ) = @_;
	
	#Windows run
	if ( $^O =~ /^MSWin(?:32|64)$/i ) {
		warn( "`del \/F \/Q ${dir}\\${file}'\n" );
		qx/del \/F \/Q ${dir}\\${file}/;
	}
	#Linux run
	elsif ( $^O =~ /^Linux$/i || $^O =~ /^Cygwin$/i || $^O =~ /^MSys$/i ) {
		qx/rm -fv ${dir}\/${file} 1>&2/;
	}
	else {
		die( "emaccconf::portableRemoveFileFromDir: unsupported operating system '$^O' !\n" );
	}
}

sub portableRemoveAllFilesFromDir( $ ) {
	my( $dir ) = $_[0];
	
	#Windows run
	if ( $^O =~ /^MSWin(?:32|64)$/i ) {
		warn( "`/del \/F \/Q ${dir}\\'\n" );
		qx/del \/F \/Q ${dir}\\/;
	}
	#Linux run
	elsif ( $^O =~ /^Linux$/i || $^O =~ /^Cygwin$/i || $^O =~ /^MSys$/i ) {
		qx/rm -fv ${dir}\/* 1>&2/;
	}
	else {
		die( "emaccconf::portableRemoveAllFilesFromDir: unsupported operating system '$^O' !\n" );
	}
}

sub portableForkAndDetach( $ ) {
	my( $cmd ) = $_[0];
	
	#Windows run
	if ( $^O =~ /^MSWin(?:32|64)$/i ) {
		warn( "`start /B ${cmd}'\n" );
		system( "start /B ${cmd}" );
	}
	#Linux run
	elsif ( $^O =~ /^Linux$/i || $^O =~ /^Cygwin$/i || $^O =~ /^MSys$/i ) {
		warn( "`${cmd} &'\n" );
		system( "${cmd} &" );	
	}
	else {
		die( "emaccconf::portableForkAndDetach: unsupported operating system '$^O' !\n" );
	}
}

1;
