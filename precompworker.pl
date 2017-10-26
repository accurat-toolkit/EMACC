#!/usr/bin/perl -w

# ver 0.1, Radu ION, 10.01.2011, worker for the parallel version of the script. Created.
# ver 0.2, Radu ION, 10.01.2011, added time limits for precomputeLexAlProbsSimple(). 
# ver 0.3, Radu ION, 09.02.2011, improved time limits for precomputeLexAlProbsSimple().
# ver 0.4, Radu ION, 24.08.2011, added user comfort modifications.
# ver 0.5, Radu ION, 30.09.2011, portable Linux and Windows.
# ver 0.6, Radu ION, 07.10.2011, fixed some memory inefficient parts (like storing all the document pairs in memory).

use strict;
use warnings;
use IO::Handle;
use Time::HiRes qw( time alarm sleep );
use File::Spec;
use File::Path;

sub portableCopyFileToDir( $$ );
sub portableRemoveFile( $ );
sub probXTAlign( $$ );
sub precomputeLexAlProbsSimple( $$$$$$ );
sub readTextDocument( $$$ );
sub readGIZAPPDict( $$ );
sub readInputParams( $ );
sub readStopWordsList( $ );
sub readInflectionList( $ );
sub lemmatizeWord( $$ );
sub timerStart( $ );
sub timerStop( $ );

if ( scalar( @ARGV ) != 1 ) {
	die( "precompworker.pl <input file from emacc.pl>\n" );
}

my( $batchfile ) = $ARGV[0];
my( $outfilebasename, $emaccconf, $howmanydocpairs ) = readInputParams( $batchfile );

##### CONFIG SECTION #########################
#Please do not modify here! Modify 'emaccconf.pm' instead.
my( $SRCL ) = $emaccconf->{"SRCL"};
my( $TRGL ) = $emaccconf->{"TRGL"};
my( $LEXALSCORE ) = $emaccconf->{"LEXALSCORE"};
my( $PLEXTYPE ) = $emaccconf->{"PLEXTYPE"};
my( $PROBTYPE ) = $emaccconf->{"PROBTYPE"};
my( $LEMMAS ) = $emaccconf->{"LEMMAS"};
my( $NFSMOUNTPOINT ) = $emaccconf->{"NFSMOUNTPOINT"};
my( $LOCALMOUNTPOINT ) = $emaccconf->{"LOCALMOUNTPOINT"};

#Try to create the temp dir if it does not exist...
mkpath( $LOCALMOUNTPOINT );

#This flag says if we copy the results to the NFS location. Default yes (1) but on cygwin, this is 0 and the copying is done by hand.
#Don't touch this.
my( $COPYONNFSONFINISH ) = 1;

my( %STOPWEN ) = readStopWordsList( $emaccconf->{"STOPWSRCL"} );
my( %STOPWRO ) = readStopWordsList( $emaccconf->{"STOPWTRGL"} );
my( %INFLEN ) = readInflectionList( $emaccconf->{"INFLSRCL"} );
my( %INFLRO ) = readInflectionList( $emaccconf->{"INFLTRGL"} );

my( $DICTINVERSE ) = $emaccconf->{"DICTINVERSE"};
my( $ENRODICT );

if ( $DICTINVERSE ) {
	$ENRODICT = readGIZAPPDict( $emaccconf->{"DICTFILE"}, "inverse" );
}
else {
	$ENRODICT = readGIZAPPDict( $emaccconf->{"DICTFILE"}, "direct" );
}

my( $DEFSMALLP ) = $emaccconf->{"DEFSMALLP"};
##### END CONFIG SECTION #####################

#Reopen STDERR to see what's wrong.
my( $errfileLG ) = File::Spec->catfile( $LOCALMOUNTPOINT, $outfilebasename . ".err" );

open( FERROR, "> $errfileLG" ) or die( "precompworker::main: cannot open FERROR!\n" );
STDERR->fdopen( \*FERROR, 'w' ) or die( "precompworker::main: cannot reopen STDERR!\n" );
STDERR->autoflush( 1 );

my( $docfileLG ) = File::Spec->catfile( $LOCALMOUNTPOINT, $outfilebasename . "-dpp.out" );
my( $teqfileLG ) = File::Spec->catfile( $LOCALMOUNTPOINT, $outfilebasename . "-teq.out" );

precomputeLexAlProbsSimple( $batchfile, $SRCL, $TRGL, $docfileLG, $teqfileLG, $howmanydocpairs );

close( FERROR );

if ( $COPYONNFSONFINISH ) {
	#Copy stuff on NFS...
	portableCopyFileToDir( $docfileLG, $NFSMOUNTPOINT );
	portableCopyFileToDir( $teqfileLG, $NFSMOUNTPOINT );
	portableRemoveFile( $errfileLG );
	portableRemoveFile( $docfileLG );
	portableRemoveFile( $teqfileLG );
}

#Guard against reading incomplete files.
my( $readyfileLG ) = File::Spec->catfile( $LOCALMOUNTPOINT, $outfilebasename . ".ready" );

open( RDY, "> $readyfileLG" ) or die( "precompworker::main: cannot open file '$readyfileLG' !\n" );
close( RDY );

if ( $COPYONNFSONFINISH ) {
	#Copy READY onto NFS...
	portableCopyFileToDir( $readyfileLG, $NFSMOUNTPOINT );
	portableRemoveFile( $readyfileLG );
}

#################### from emaccconf.pm ############################################
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
#################### end from emaccconf.pm ########################################

sub probXTAlign( $$ ) {
	my( $docen, $docro ) = @_;
	my( $comm ) = 0;
	my( $N ) = $docen->{"_N"} + $docro->{"_N"};
	
	foreach my $ie ( keys( %{ $docen } ) ) {
		next if ( $ie eq "_N" );
		
		foreach my $ir ( keys( %{ $docro } ) ) {
			next if ( $ir eq "_N" );
			
			if ( $ie eq $ir ) {
				$comm += $docen->{$ie} + $docro->{$ir};
				last;
			}
		}
	}
	
	#0 div guard.
	return $DEFSMALLP if ( $N == 0 );
	return $comm / $N;
}

#Simplu, fara clustere
sub precomputeLexAlProbsSimple( $$$$$$ ) {
	#en means source and ro means target
	my( $jobfile, $langen, $langro, $docfileL, $teqfileL, $docpairsno ) = @_;
	my( $crtpaircnt ) = 0;
	
	open( DPF, "> $docfileL" ) or die( "precompworker::precomputeLexAlProbsSimple: cannot open file '$docfileL' !\n" );
	binmode( DPF, ":utf8" );
			
	open( TEQD, "> $teqfileL" ) or die( "precompworker::precomputeLexAlProbsSimple: cannot open file '$teqfileL' !\n" );
	binmode( TEQD, ":utf8" );
	
	open( BATCH, "<", $jobfile ) or die( "precompworker::precomputeLexAlProbsSimple: cannot open file '$jobfile' !\n" );
	binmode( BATCH, ":utf8" );
	
	#Basename of all the files (ignored).
	my $outfbn = <BATCH>;
	
	while ( my $line = <BATCH> ) {
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		next if ( $line =~ /^--param\s/ );
		
		my( $de, $dr ) = split( /\t+/, $line );
		$crtpaircnt++;
			
		warn( "precompworker::precomputeLexAlProbsSimple[$outfilebasename]: " . $crtpaircnt . "/" . $docpairsno . "\n" )
			if (  $crtpaircnt % 100 == 0 );

		print( TEQD $de . "#:#" . $dr . "\n" );
		
		#Config options !
		my( $xcesde ) = readTextDocument( $de, \%STOPWEN, \%INFLEN );
		my( $xcesdr ) = readTextDocument( $dr, \%STOPWRO, \%INFLRO );
		my( %writtenpairs ) = ();
		my( $plexsum ) = 0;
		my( $plexsums, $plexsumt ) = ( 0, 0 );
		my( %cen, %cro ) = ();
		
		#Skip _N !!
		foreach my $dek ( keys( %{ $xcesde } ) ) {
			next if ( $dek eq "_N" );
			
			#The English word has Romanian translations...
			if ( exists( $ENRODICT->{$dek} ) ) {
				#See if we can find them in our Romanian document...
				my( @drkeys ) = keys( %{ $ENRODICT->{$dek} } );
				
				foreach my $drk ( @drkeys ) {
					#Is the translation in the dictionary?
					if ( exists( $xcesdr->{$drk} ) && $ENRODICT->{$dek}->{$drk} >= $LEXALSCORE ) {
						#0 div guard.
						my( $teqp ) = 0;
						my( $totalpairs ) = $xcesde->{"_N"} * $xcesdr->{"_N"};
						
						if ( $totalpairs > 0 ) {
							#'giza' to retain the prob from the dict, 'comp' to recompute it according to doc pair.
							SWPTYPE: {
								$PROBTYPE eq "comp" and do {
									$teqp = ( $xcesde->{$dek} * $xcesdr->{$drk} ) / $totalpairs;
									last;
								};
							
								$PROBTYPE eq "giza" and do {
									$teqp = $ENRODICT->{$dek}->{$drk};
									last;
								};
							}
						}
					
						#They are too small ... maybe scale all by 10000 ?...
						print( TEQD $dek . "#" . $drk . "\t" . $teqp . "\n" ) if ( ! exists( $writtenpairs{$dek . "#" . $drk} ) );
						$writtenpairs{$dek . "#" . $drk} = 1;
						
						$plexsum += ( $xcesde->{$dek} * $xcesdr->{$drk} );
						
						$plexsums += $xcesde->{$dek} if ( ! exists( $cen{$dek} ) );
						$cen{$dek} = 1;
						$plexsumt += $xcesdr->{$drk} if ( ! exists( $cro{$drk} ) );
						$cro{$drk} = 1;
					} #end dek and drk are TEQ
				} #end all Romanian toks that English tok has as translations.
			} #end if English tok has translations.
		} #end English.
		
		#_N marks the total frequency
		#0 div guard.
		my( $totalpairs ) = $xcesde->{"_N"} * $xcesdr->{"_N"};
		my( $transprobfair ) = 0;
		my( $transprobbiased ) = 0;

		if ( $totalpairs > 0 ) {
			$transprobfair = $plexsum / $totalpairs;
			$transprobbiased = ( $plexsums * $plexsumt ) / $totalpairs;
		}
			
		#print( STDERR $transprobfair . "\n" );
		#print( STDERR $transprobbiased . "\n\n" );
			
		my( $transprob ) = 0;
			
		$transprob = $transprobfair if ( $PLEXTYPE eq "fair" );
		$transprob = $transprobbiased if ( $PLEXTYPE eq "biased" );
		
		print( DPF $de . "\t" . $dr . "\t" . $transprob . "\n" );
	} #end all pairs of docs
	
	close( BATCH );
	close( TEQD );
}

sub readTextDocument( $$$ ) {
	my( $corpus, $stopwlist, $infllist ) = @_;
	my( %document ) = ( "_N" => 0 );
	
	open( COR, "< $corpus" ) or die( "precompworker::readTextDocument: cannot open file '$corpus' !\n" );
	binmode( COR, ":utf8" );
	
	while ( my $line = <COR> ) {
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		my( @toks ) = split( /\s+/, $line );
		
		foreach my $t ( @toks ) {
			$t = lc( $t );
			
			#Remove punctuation from beginnig and end
			$t =~ s/^\W+//;
			$t =~ s/\W+$//;
			
			#Next if punctuation ...
			next if ( $t eq "" );
			#No stop words ...
			next if ( exists( $stopwlist->{$t} ) );
			#Other skip conditions go here ...
			
			#Lemmatize $t ...
			if ( $LEMMAS ) {
				$t = lemmatizeWord( $t, $infllist );
			}

			if ( ! exists( $document{$t} ) ) {
				$document{$t} = 1;
			}
			else {
				$document{$t}++;
			}
			
			$document{"_N"}++;
		} #end crt fragment
	} #end while
	
	close( COR );
	
	return \%document;
}

sub readGIZAPPDict( $$ ) {
	my( %gizapp ) = ();
	my( $gizafile, $dir ) = @_;

	open( DICT, "< $gizafile" ) or die( "precompworker::readGIZAPPDict : cannot open file \'$gizafile\' !\n" );
	binmode( DICT, ":utf8" );
	
	while ( my $line = <DICT> ) {
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		next if ( $line =~ /^$/ );

		my( @toks ) = split( /\s+/, $line );
		#en - SRC, ro - TRG
		my( $enw ) = lc( $toks[0] );
		my( $row ) = lc( $toks[1] );
		my( $score ) = $toks[2];
		
		if ( $LEMMAS ) {
			$enw = lemmatizeWord( $enw, \%INFLEN );
			$row = lemmatizeWord( $row, \%INFLRO );
		}		
		
		SWTD: {
			$dir eq "direct" and do {
				if ( ! exists( $gizapp{$enw} ) ) {
					$gizapp{$enw} = { $row => $score };
				}
				else {
					$gizapp{$enw}->{$row} = $score;
				}
				
				last;
			};

			$dir eq "inverse" and do {
				if ( ! exists( $gizapp{$row} ) ) {
					$gizapp{$row} = { $enw => $score };
				}
				else {
					$gizapp{$row}->{enw} = $score;
				}
				
				last;
			};

			die( "precompworker::readGIZAPPDict : $dir is not valid !!\n" );
		}
	}
	
	close( DICT );
	return \%gizapp;
}

sub readInputParams( $ ) {
	my( $outfbn );
	my( @docsp ) = ();
	my( %conf ) = ();
	my( $npairs ) = 0;

	open( IN, "< " . $_[0] ) or die( "precompworker::readInputParams: cannot open file " . $_[0] . " !\n" );

	$outfbn = <IN>;

	$outfbn =~ s/^\s+//;
	$outfbn =~ s/\s+$//;

	while ( my $line = <IN> ) {
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;

		if ( $line =~ /^--param\s/ ) {
			$line =~ s/^--param\s//;

			my( $param, $value ) = split( /\s*=\s*/, $line );
			
			$conf{$param} = $value;
			next;
		}
		else {
			$npairs++;
		}
	}

	close( IN );

	return ( $outfbn, \%conf, $npairs );
}

sub readStopWordsList( $ ) {
	my( %swl ) = ();
	
	open( SWL, "< $_[0]" ) or die( "precompworker::readDocList: cannot open file '$_[0]' !\n" );
	binmode( SWL, ":utf8" );
	
	while ( my $line = <SWL> ) {
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		$swl{lc( $line )} = 1;
	}
	
	close( SWL );
	return %swl;
}

sub readInflectionList( $ ) {
	my( %infl ) = ( "LONGEST" => 0 );
	
	open( INF, "< $_[0]" ) or die( "precompworker::readInflectionList: cannot open file '$_[0]' !\n" );
	binmode( INF, ":utf8" );
	
	while ( my $line = <INF> ) {
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		$infl{lc( $line )} = length( lc( $line ) );
		
		if ( $infl{"LONGEST"} < length( lc( $line ) ) ) {
			$infl{"LONGEST"} = length( lc( $line ) );
		}
	}
	
	close( INF );
	return %infl;
}

sub lemmatizeWord( $$ ) {
	my( $word, $infllist ) = @_;
	
	#Endings are in lowercase
	$word = lc( $word );
	
	if ( length( $word ) > $infllist->{"LONGEST"} ) {
		my( @wordlett ) = split( //, $word );
		my( @wordend ) = @wordlett[$#wordlett - $infllist->{"LONGEST"} + 1 .. $#wordlett];
		
		#Match the longest suffix ...
		for ( my $i = 0; $i < scalar( @wordend ); $i++ ) {
			my( $crtsfx ) = join( "", @wordend[$i .. $#wordend] );
			
			if ( exists( $infllist->{$crtsfx} ) ) {
				$word =~ s/${crtsfx}$//;
		
				return $word if ( $word ne "" );
			}
		}
	}
	
	return $word;
}

sub timerStart( $ ) {
	my( $timer ) = $_[0];
	my( $tm ) = time();
	
	$$timer = $tm;
}

sub timerStop( $ ) {
	my( $timer ) = $_[0];
	my( $tm ) = time();
	
	return $tm - $$timer;
}
