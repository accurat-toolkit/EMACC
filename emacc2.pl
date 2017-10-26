#!/usr/bin/perl -w

# EMACC: textual unit alignment with Expectation-Maximization.
# 
# (C) ICIA 2011, Radu ION.
#
# ver 0.1, Radu ION, 18.01.2011, am eliminat %dptable ca sa nu mai incarce tot in memorie.
# ver 0.2, Radu ION, 19.01.2011, am elimiant %teq2docpair urmand ca sa o calculez la timpul potrivit.
# ver 0.3, Radu ION, 21.01.2011, am adaugat update percent: numai alinierile de documente din prima parte 
#   (ca acest procent din toate alinierile generate) sunt date la iesire iar updateul se face numai pentru acestea la fiecare pas.
# ver 0.4, Radu ION, 21.01.2011, EM iese daca probabilitatea aliniamentului curent e mai proasta decat cea anterioara.
# ver 0.5, Radu ION, 4.02.2011, added createZTableD1 and D2
# ver 1.0, Radu ION, 24.08.2011, added emacc-precompute-p.pl into this file.
# ver 1.1, Radu ION, 24.08.2011, parallel version (-p) of that script. Created.
# ver 1.2, Radu ION, 24.08.2011, modified readDocList().
# ver 1.3, Radu ION, 24.08.2011, added distributeEvenly().
# ver 1.4, Radu ION, 30.09.2011, added M:N alignments (by permitting 1:N alignments).
#	Now, it does not matter which list of documents contains fewer documents.
# ver 1.5, Radu ION, 30.09.2011, fixed the random behaviour of EMACC where when ran with the same parameters, returned different results.
#	Thanks Nikos.
# ver 1.6, Radu ION, 06.10.2011, fixed distributeEvenly().
# ver 1.7, Radu ION, 08.10.2011, memory usage optimization.
# ver 1.8, Radu ION, 10.10.2011, added EMACCMODE and emaccOne().
# ver 2.0, Radu ION, 10.01.2011, heavy memory optimization. Working with file (on HDD) matrices now...
# ver 2.1, Radu ION, 25.10.2011, some more memory optimization.
# ver 2.2, Radu ION: 17.11.2011, added emaccSimple for EMACCMODE.

use strict;
use warnings;
use Sys::Hostname;
use File::Spec;
use Time::HiRes qw( time alarm sleep );
use emaccconf;
use hddmatrix;

sub emaccFull();
sub emaccSimple();
#Initial distribuion D2 (see the BUCC 2011 paper) - with the document alignments cache
sub createZTableD2( $$$ );
#Initial distribution D1 (see the BUCC 2011 paper) - uniform distribution
sub createZTableD1( $$ );
sub normZTable( $$ );
sub probPLexAlign( $$$ );
sub generateRandomAlignment( $$ );
sub readAlignTEQPs( $ );
sub getAllTeq4Dp( $$$ );
sub normalizePLexTable( $ );
sub updatePLexTable( $$ );
sub checkEqual( $$$$$ );
sub readDocList1( $$ );
sub timerStart( $ );
sub timerStop( $ );
## emacc-precomp-p.pl
sub emaccPrecompP( $$ );
sub runInParallel( $$$$ );
sub readDocList2( $ );
sub readDocLists2( $$ );
sub readCluster( $ );
sub distributeEvenly( $$$ );
## command line
sub normalizeLang( $ );
sub readCmdLineArguments( @ );
sub cleanTemp();

if ( scalar( @ARGV ) < 2 ) {
	die( "emacc.pl \\
	[--source en] [--target ro] \\
	[--param EMLOOPS=5] \\
	[--param EMACCMODE=emacc-full] \\
	[--param MAXTARGETALIGNMENTS=3] \\
	[--param INIDISTRIB=D2] \\
	[--param TEQPUPDATETHR=0.4] \\
	[--param LEXALSCORE=0.4] \\
	[--param PROBTYPE=giza] \\
	[--param DICTINVERSE=0] \\
	[--param CLUSTERFILE=generate] \\
	--input <source language document list file> \\
	--input <target language document list file> \\
	[--output <output file name>]\n" );
}

my( $cmdlineconfh ) = readCmdLineArguments( @ARGV );
my( $emaccconfh ) = emaccconf->new( $cmdlineconfh );

###### CONFIG1 ###########################################
#Please do not modify here! Modify 'emaccconf.pm' instead.
my( $SRCL ) = $emaccconfh->{"SRCL"};
my( $TRGL ) = $emaccconfh->{"TRGL"};
my( $EMLOOPS ) = $emaccconfh->{"EMLOOPS"};
my( $EMACCMODE ) = $emaccconfh->{"EMACCMODE"};
my( $DEFSMALLP ) = $emaccconfh->{"DEFSMALLP"};
my( $TEQPUPDATETHR ) = $emaccconfh->{"TEQPUPDATETHR"};
my( $OUTPUTPERCENT ) = $emaccconfh->{"OUTPUTPERCENT"};
my( $UPDATEPERCENT ) = $emaccconfh->{"UPDATEPERCENT"};
my( $INIDISTRIB ) = $emaccconfh->{"INIDISTRIB"};
my( $DALOUTFILE ) = $emaccconfh->{"DALOUTFILE"};
###### END CONFIG1 #######################################

##### CONFIG2 ############################################
#Please do not modify here! Modify 'emaccconf.pm' instead.
my( $MAXTARGETALIGNMENTS ) = $emaccconfh->{"MAXTARGETALIGNMENTS"};
my( $LEXALSCORE ) = $emaccconfh->{"LEXALSCORE"};
my( $PLEXTYPE ) = $emaccconfh->{"PLEXTYPE"};
my( $LEMMAS ) = $emaccconfh->{"LEMMAS"};
my( $CLUSTERFILE ) = $emaccconfh->{"CLUSTERFILE"};
my( $NFSMOUNTPOINT ) = $emaccconfh->{"NFSMOUNTPOINT"};
my( $LOCALMOUNTPOINT ) = $emaccconfh->{"LOCALMOUNTPOINT"};
my( $PCFILEOUT ) = $emaccconfh->{"PRECOMPMODELFILE"};
my( $TEQFILEOUT ) = $emaccconfh->{"DPAIRTEQMODELFILE"};
##### END CONFIG2 ########################################

cleanTemp();

my( %DOCALGS ) = readDocList1( $cmdlineconfh->{"INPUTFILE1"}, $cmdlineconfh->{"INPUTFILE2"} );
my( @DOCSEN ) = keys( %{ $DOCALGS{$SRCL}->{"u"} } );
my( @DOCSRO ) = keys( %{ $DOCALGS{$TRGL}->{"u"} } );

#1. Precompute distributions
emaccPrecompP( $cmdlineconfh->{"INPUTFILE1"}, $cmdlineconfh->{"INPUTFILE2"} );

#2. Run EMACC
SWEMMODE: {
	#If you have up to around 1M of possible document pairs, do full EM re-estimation.
	$EMACCMODE eq "emacc-full" and do { emaccFull(); last; };
	#else, do a simple, maximal probability alignment based on alignment probabilities (models too big to fit in memory).
	$EMACCMODE eq "emacc-simple" and do { emaccSimple(); last; };
}

#Remove intermediary files.
#Well, don't. We may need them.
#emaccconf::portableRemoveFile( $PCFILEOUT );
#emaccconf::portableRemoveFile( $TEQFILEOUT );

##### End main ###################

sub emaccSimple() {
	my( %probfreq ) = ();
	my( $DP ) = 0;
	my( $linecount ) = 0;
	
	#1. Construct probability frequencies
	open( PCF, "<", $PCFILEOUT ) or die( "emacc2::emaccSimple: cannot open file '$PCFILEOUT' !\n" );
	binmode( PCF, ":utf8" );
	
	while ( my $line = <PCF> ) {
		$linecount++;
		
		print( STDERR "emacc2::emaccSimple[first read]: read $linecount lines...\n" )
			if ( $linecount % 1000000 == 0 );
		
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		my( $srcdoc, $trgdoc, $prob ) = split( /\t+/, $line );
		my( $prob5dec ) = int( $prob * 100000 ) / 100000;
		
		if ( ! exists( $probfreq{$prob5dec} ) ) {
			$probfreq{$prob5dec} = 1;
		}
		else {
			$probfreq{$prob5dec}++;
		}
		
		$DP++;
	}
	
	close( PCF );
	
	#2. Sort the probabilities and determine the threshold for which we have OUTPUTPERCENT alignments
	my( @docalignprobs ) = sort { $b <=> $a } keys( %probfreq );
	my( $DPOUT ) = 0;
	my( $PTHR ) = 0;
	
	foreach my $p ( @docalignprobs ) {
		$DPOUT += $probfreq{$p};
		
		if ( $DPOUT / $DP >= $OUTPUTPERCENT ) {
			$PTHR = $p;
			last;
		}
	}
	
	$linecount = 0;
	
	#3. Select all document pairs for which p is larger than or equal to $PTHR
	open( PCF, "<", $PCFILEOUT ) or die( "emacc2::emaccSimple: cannot open file '$PCFILEOUT' !\n" );
	binmode( PCF, ":utf8" );

	open( DAL, ">", $DALOUTFILE ) or die( "emacc2::emaccSimple: cannot open alignment file !\n" );
	binmode( DAL, ":utf8" );
	
	while ( my $line = <PCF> ) {
		$linecount++;
		
		print( STDERR "emacc2::emaccSimple[second read]: read $linecount lines...\n" )
			if ( $linecount % 1000000 == 0 );
		
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		my( $srcdoc, $trgdoc, $prob ) = split( /\t+/, $line );
		my( $prob5dec ) = int( $prob * 100000 ) / 100000;
		
		if ( $prob5dec >= $PTHR ) {
			print( DAL $srcdoc . "\t" . $trgdoc . "\t" . $prob . "\n" );
		}
	}
	
	close( PCF );
	close( DAL );
}

sub emaccFull() {
	my( $ztable );

	SWINIDIST: {
		$INIDISTRIB eq "D1" and do {
			$ztable = createZTableD1( \@DOCSEN, \@DOCSRO );
			last;
		};
	
		$INIDISTRIB eq "D2" and do {
			$ztable = createZTableD2( \@DOCSEN, \@DOCSRO, $PCFILEOUT );
			last;
		};
	
		die( "emacc2::emaccFull: unknown distribution '$INIDISTRIB' !\n" );
	}

	my( $lexalcache );
	my( $plextable );
	my( @bestalignment ) = ();

	$plextable = readAlignTEQPs( $TEQFILEOUT );
	normalizePLexTable( $plextable );

	for ( my $i = 1; $i <= $EMLOOPS; $i++ ) {
		#1. Compute all document pairs alignment probabilities.
		#Update table for the document pair probabilities ...
		my( $ztablecnt ) = hddmatrix->newKey( \@DOCSEN, \@DOCSRO, "d" );
		my( $cnt ) = 1;
		my( $tmref ) = 0;

		timerStart( \$tmref );
		print( STDERR "emacc2::emaccFull: Iteration $i: step 1 ...\n" );
	
		my( $iterdpfile ) = File::Spec->catfile( $LOCALMOUNTPOINT, "$SRCL-$TRGL-docpairsfile-$i.dpout" );
	
		open( DPF, ">", $iterdpfile ) or die( "emacc2::emaccFull[$i]: cannot open file '$iterdpfile' !\n" );
		binmode( DPF, ":utf8" );
		
		EN:
		foreach my $den ( @DOCSEN ) {

			RO:
			foreach my $dro ( @DOCSRO ) {
				$cnt++;
				print( STDERR "Iteration $i: $cnt/" . ( scalar( @DOCSEN ) * scalar( @DOCSRO ) ) . "\n" )
					if ( $cnt % 10000 == 0 );
				
				#Compute E-step ...
				#Here we go through z's values ...
				#z=true only
				my( $sterm ) = 0;
				my( $uterm ) = 0;

				#New version
				$uterm = probPLexAlign( $den, $dro, $plextable ) * $ztable->read( $den, $dro );

				if ( $uterm > 0 ) {
					$sterm = log( $uterm );
				}
				else {
					$sterm = -1e+20;
				}

				$ztablecnt->write( $den, $dro, $uterm );

				my( $estep ) = $sterm;

				print( DPF $den . "\t" . $dro . "\t" . $estep . "\n" );
			} #end all ro docs
		} #end all en docs
		
		close( DPF );
	
		print( STDERR "\n" );
		print( STDERR "emacc2::emaccFull: step 1 elapsed time: " . timerStop( \$tmref ) . ".\n" );
	
		print( STDERR "\n" );
		print( STDERR "emacc2::emaccFull: Iteration $i: step 2 ...\n" );
	
		timerStart( \$tmref );

		#2. Choose the assignment with the greatest alignment probability.
		#ver 1.4 with M:N alignments.
		#The new assignment (greedy method)
		my( @theta ) = ();
		my( %alreadyaligneden ) = ();
		my( %alreadyalignedro ) = ();		

		open( DPF, "<", $iterdpfile ) or die( "emacc2::emaccFull[$i]: cannot open file '$iterdpfile' !\n" );
		binmode( DPF, ":utf8" );
		
		while ( my $line = <DPF> ) {
			$line =~ s/^\s+//;
			$line =~ s/\s+$//;
	
			my( $de, $dr, $p ) = split( /\t+/, $line );
			my( $subaddalignment ) = sub {
				my( $alreadyaligned, $end, $rod, $prob ) = @_;
			
				if ( ! exists( $alreadyaligned->{$end} ) ) {
					$alreadyaligned->{$end} = { $rod => $prob };
				}
				#We did not reach the required target alignments
				elsif ( scalar( keys( %{ $alreadyaligned->{$end} } ) ) < $MAXTARGETALIGNMENTS ) {
					$alreadyaligned->{$end}->{$rod} = $prob;
				}
				#Let's store the best MAXTARGETALIGNMENTS alignments.
				else {
					my( @raligns ) = map { [ $_, $alreadyaligned->{$end}->{$_} ] } keys( %{ $alreadyaligned->{$end} } );
			
					push( @raligns, [ $rod, $prob ] );
				
					my( @sraligns ) = sort { $b->[1] <=> $a->[1] } @raligns;
			
					pop( @sraligns );
			
					$alreadyaligned->{$end} = {};
			
					foreach my $drpp ( @sraligns ) {
						my( $ndr, $np ) = @{ $drpp };
				
						$alreadyaligned->{$end}->{$ndr} = $np;
					}
				}
			};
		
			$subaddalignment->( \%alreadyaligneden, $de, $dr, $p );
			$subaddalignment->( \%alreadyalignedro, $dr, $de, $p );
		} #end all doc pairs
		
		close( DPF );
		
		my( %alreadyalignedenro ) = ();
	
		#Intersect alignments from EN to RO with those from RO to EN keeping our direction.
		foreach my $de ( keys( %alreadyaligneden ) ) {
			foreach my $dr ( keys( %{ $alreadyaligneden{$de} } ) ) {
				if ( exists( $alreadyalignedro{$dr}->{$de} ) ) {

					if ( ! exists( $alreadyalignedenro{$de} ) ) {
						$alreadyalignedenro{$de} = { $dr => $alreadyaligneden{$de}->{$dr} };
					}
					else {
						$alreadyalignedenro{$de}->{$dr} = $alreadyaligneden{$de}->{$dr};
					}
				}
			}
			#If nothing in intersection... do nothing. Unreliable alignments.
		}
	
		#Compose @theta and @thetasc
		foreach my $de ( keys( %alreadyalignedenro ) ) {
			foreach my $dr ( keys( %{ $alreadyalignedenro{$de} } ) ) {
				my( $p ) = $alreadyalignedenro{$de}->{$dr};
				
				push( @theta, [ $de, $dr, $p ] );
			}
		}
	
		print( STDERR "emacc2::main: step 2 elapsed time: " . timerStop( \$tmref ) . ".\n" );

		#3. Update the tables and repeat the EM process ...
		#Performing update of the ztable ...
		normZTable( $ztablecnt, $ztable );

		#New version
		updatePLexTable( $plextable, \@theta );
		normalizePLexTable( $plextable );

		#Update parameters
		@bestalignment = @theta;
		#And go with the re-estimation !!
		
		#Remove intermediate calculations.
		emaccconf::portableRemoveFile( $iterdpfile );
	} #end EM loop

	my( @sortedbestalignment ) = sort { $b->[2] <=> $a->[2] } @bestalignment;

	open( DAL, ">", $DALOUTFILE ) or die( "emacc2::emaccFull: cannot open alignment file!\n" );

	my( $normsum ) = 0;

	for ( my $i = 0; $i < scalar( @bestalignment ); $i++ ) {
		#We have output exactly $OUTPUTPERCENT of the pairs...
		last if ( $i / scalar( @bestalignment ) > $OUTPUTPERCENT );
	
		$normsum += exp( $sortedbestalignment[$i]->[2] );
	}

	#Final output from EMACC
	for ( my $i = 0; $i < scalar( @bestalignment ); $i++ ) {
		#We have output exactly $OUTPUTPERCENT of the pairs...
		last if ( $i / scalar( @bestalignment ) > $OUTPUTPERCENT );
	
		my( $srcd ) = $sortedbestalignment[$i]->[0];
		my( $trgd ) = $sortedbestalignment[$i]->[1];
	
		#Get back to documents from ints
		print( DAL $srcd . "\t" . $trgd . "\t" . exp( $sortedbestalignment[$i]->[2] ) / $normsum . "\n" );
	}

	close( DAL );
} #end emaccFull

#Doing general clean-up from a previous run so as to start clean.
sub cleanTemp() {
	emaccconf::portableRemoveFileFromDir( $NFSMOUNTPOINT, "*.in" );
	emaccconf::portableRemoveFileFromDir( $NFSMOUNTPOINT, "*.ready" );
	emaccconf::portableRemoveFileFromDir( $NFSMOUNTPOINT, "*.out" );
	emaccconf::portableRemoveFileFromDir( $NFSMOUNTPOINT, "*.err" );
	emaccconf::portableRemoveAllFilesFromDir( $LOCALMOUNTPOINT );
}

####################################### Begin EM optimization #######################################################

#Uniform distribution ztable...
#hddmatrix ready
sub createZTableD1( $$ ) {
	my( $docsen, $docsro ) = @_;
	my( $ztbl ) = hddmatrix->newKey( $docsen, $docsro, "d" );
	my( $nsum ) = scalar( @{ $docsen } ) * scalar(  @{ $docsro } );

	#1. Create ztable...
	foreach my $de ( @{ $docsen } ) {
		foreach my $dr ( @{ $docsro } ) {
			$ztbl->write( $de, $dr, 1 / $nsum );
		} #end ro doc
	} #end en docs
	
	#Normalize ZTable ... already normalized.
	return $ztbl;
}

#ok
#This also receives the doc alignment probabilities as a seed.
#hddmatrix ready
sub createZTableD2( $$$ ) {
	my( $docsen, $docsro, $docalfile ) = @_;
	my( $ztbl ) = hddmatrix->newKey( $docsen, $docsro, "d" );
	my( $nsum ) = 0;
	my( $lcnt ) = 0;

	#1. Read document aligments probabilities...
	open( DAL, "< $docalfile" ) or die( "emacc2::createZTableD2: cannot open file '$docalfile' !\n" );
	
	while ( my $line = <DAL> ) {
		$lcnt++;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		print( STDERR "emacc2::createZTableD2: read '$lcnt' lines ...\n" ) if ( $lcnt % 1000000 == 0 );

		my( $srcd, $trgd, $p ) = split( /\t+/, $line );

		$ztbl->write( $srcd, $trgd, $p );
		$nsum += $p;
	}
	
	close( DAL );
	
	#Normalize ZTable
	foreach my $de ( @{ $docsen } ) {
		foreach my $dr ( @{ $docsro } ) {
			my( $p ) = $ztbl->read( $de, $dr );

			$ztbl->write( $de, $dr, $p / $nsum ) if ( $nsum > 0 );
		}
	}
	
	return $ztbl;
}

#hddmatrix ready
sub normZTable( $$ ) {
	my( $update, $ztable ) = @_;
	my( $ztsum ) = 0;
	my( $nodp ) = 0;
	
	foreach my $de ( @{ $ztable->getRowKeys() } ) {
		foreach my $dr ( @{ $ztable->getColumnKeys() } ) {
			$nodp++;
			print( STDERR "emacc2::normZTable: read '$nodp' pairs ...\n" ) if ( $nodp % 1000000 == 0 );
			
			my( $p ) = $ztable->read( $de, $dr );

			$p += $update->read( $de, $dr );
			$ztable->write( $de, $dr, $p );
			$ztsum += $p;
		}
	}
	
	$nodp = 0;
	
	foreach my $de ( @{ $ztable->getRowKeys() } ) {
		foreach my $dr ( @{ $ztable->getColumnKeys() } ) {
			$nodp++;
			print( STDERR "emacc2::normZTable: written '$nodp' pairs ...\n" ) if ( $nodp % 1000000 == 0 );

			my( $p ) = $ztable->read( $de, $dr );
			
			#0 div guard.
			$ztable->write( $de, $dr, $p / $ztsum ) if ( $ztsum > 0 );
		}
	}
}

#memory efficient ready (ver 2.1).
sub readAlignTEQPs( $ ) {
	#TEQ are stored as ints
	my( %teqpairs ) = ( "INT2PROB" => {}, "TEQ2INT" => { "_COUNT" => 0 } );
	my( %srcd2teq ) = ( "INT2TEQ" => {}, "DOC2INT" => { "_COUNT" => 0 }, "INT2DOC" => {} );
	my( %trgd2teq ) = ( "INT2TEQ" => {}, "DOC2INT" => { "_COUNT" => 0 }, "INT2DOC" => {} );
	my( %dpkey ) = ();
	my( $readdpairs ) = 0;
	
	open( TEQ, "< $_[0]" ) or die( "emacc2::readAlignTEQPs: cannot open file $_[0] !\n" );
	
	while ( my $line = <TEQ> ) {
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		next if ( $line eq "" );
		
		if ( $line =~ /#:#/ ) {
			$readdpairs++;
			print( STDERR "emacc2::readAlignTEQPs: read '$readdpairs' document pairs ...\n" ) if ( $readdpairs % 1000000 == 0 );
			
			my( $srcd, $trgd ) = split( /#:#/, $line );
			
			%dpkey = ( "SRCD" => $srcd, "TRGD" => $trgd );
			next;
		}
		else {
			my( $pair, $prob ) = split( /\s+/, $line );
			
			if ( ! exists( $teqpairs{"TEQ2INT"}->{$pair} ) ) {
				$teqpairs{"TEQ2INT"}->{$pair} = $teqpairs{"TEQ2INT"}->{"_COUNT"};
				$teqpairs{"INT2PROB"}->{$teqpairs{"TEQ2INT"}->{$pair}} = { "PROB" => $prob, "DPAIRS" => 1 };
				$teqpairs{"TEQ2INT"}->{"_COUNT"}++;
			}
			else {
				#In how many document pairs this TEQ appears.
				$teqpairs{"INT2PROB"}->{$teqpairs{"TEQ2INT"}->{$pair}}->{"DPAIRS"}++;
			}
			
			my( $intpair ) = $teqpairs{"TEQ2INT"}->{$pair};
			my( $srcdoc ) = $dpkey{"SRCD"};
			my( $trgdoc ) = $dpkey{"TRGD"};
			
			if ( ! exists( $srcd2teq{"DOC2INT"}->{$srcdoc} ) ) {
				$srcd2teq{"DOC2INT"}->{$srcdoc} = $srcd2teq{"DOC2INT"}->{"_COUNT"};
				$srcd2teq{"INT2DOC"}->{$srcd2teq{"DOC2INT"}->{$srcdoc}} = $srcdoc;
				$srcd2teq{"DOC2INT"}->{"_COUNT"}++;
			}
			
			if ( ! exists( $trgd2teq{"DOC2INT"}->{$trgdoc} ) ) {
				$trgd2teq{"DOC2INT"}->{$trgdoc} = $trgd2teq{"DOC2INT"}->{"_COUNT"};
				$trgd2teq{"INT2DOC"}->{$trgd2teq{"DOC2INT"}->{$trgdoc}} = $trgdoc;
				$trgd2teq{"DOC2INT"}->{"_COUNT"}++;
			}
			
			my( $intsrcdoc ) = $srcd2teq{"DOC2INT"}->{$srcdoc};
			my( $inttrgdoc ) = $trgd2teq{"DOC2INT"}->{$trgdoc};
			
			if ( ! exists( $srcd2teq{"INT2TEQ"}->{$intsrcdoc} ) ) {
				$srcd2teq{"INT2TEQ"}->{$intsrcdoc} = { $intpair => 1 };
			}
			else {
				$srcd2teq{"INT2TEQ"}->{$intsrcdoc}->{$intpair} = 1;
			}
			
			if ( ! exists( $trgd2teq{"INT2TEQ"}->{$inttrgdoc} ) ) {
				$trgd2teq{"INT2TEQ"}->{$inttrgdoc} = { $intpair => 1 };
			}
			else {
				$trgd2teq{"INT2TEQ"}->{$inttrgdoc}->{$intpair} = 1;
			}			
		}
	}
	
	close( TEQ );
	return [ \%srcd2teq, \%trgd2teq, \%teqpairs ];
}

#memory efficient ready (ver 2.1).
sub getAllTeq4Dp( $$$ ) {
	my( $plextable, $srcd, $trgd ) = @_;
	my( $srcd2teq, $trgd2teq, $teqpairs ) = @{ $plextable };
	my( %teqs4dp ) = ();
	my( $intsrcd ) = $srcd2teq->{"DOC2INT"}->{$srcd};
	my( $inttrgd ) = $trgd2teq->{"DOC2INT"}->{$trgd};
	
	foreach my $itp ( keys( %{ $srcd2teq->{"INT2TEQ"}->{$intsrcd} } ) ) {
		if ( exists( $trgd2teq->{"INT2TEQ"}->{$inttrgd}->{$itp} ) ) {
			#Integer teq pair.
			$teqs4dp{$itp} = $teqpairs->{"INT2PROB"}->{$itp}->{"PROB"};
		}
	}
	
	#Integer teq pair.
	return \%teqs4dp;
}

#Input: the table.
#memory efficient ready (ver 2.1).
sub normalizePLexTable( $ ) {
	my( $tmref );
	
	timerStart( \$tmref );
	
	my( $plextable ) = $_[0];
	my( $srcd2teq, $trgd2teq, $teqpairs ) = @{ $plextable };
	
	foreach my $itp ( keys( %{ $teqpairs->{"INT2PROB"} } ) ) {
		my( $itpp, $dpnoinwhichitp ) = ( $teqpairs->{"INT2PROB"}->{$itp}->{"PROB"}, $teqpairs->{"INT2PROB"}->{$itp}->{"DPAIRS"} );
		my( $sum ) = $itpp * $dpnoinwhichitp;
		
		if ( defined( $sum ) && $sum > 0 ) {
			$teqpairs->{"INT2PROB"}->{$itp}->{"PROB"} /= $sum;
		} #end if sum is > 0
	} #end all pairs
	
	print( STDERR "emacc2::normalizePLexTable: elapsed time: " . timerStop( \$tmref ) . ".\n" );
}

#Input: the table and the best document alignment
#memory efficient ready (ver 2.1).
sub updatePLexTable( $$ ) {
	my( $tmref ) = 0;
	
	timerStart( \$tmref );
	
	my( $plextable, $align ) = @_;
	my( $srcd2teq, $trgd2teq, $teqpairs ) = @{ $plextable };
	my( %winteqmassprobs ) = ();
	my( $correctperc ) = int( $UPDATEPERCENT * scalar( @{ $align } ) );
	my( $corcnt ) = 0;

	#Through all the alignment
	foreach my $dpv ( @{ $align } ) {
		$corcnt++;
		
		last if ( $corcnt > $correctperc );
		
		my( $srcd, $trgd, $score ) = @{ $dpv };
		#Integer TEQs.
		my( $tpfordp ) = getAllTeq4Dp( $plextable, $srcd, $trgd );
		
		foreach my $itp ( keys( %{ $tpfordp } ) ) {
			my( $tpdpprob ) = $tpfordp->{$itp};
			
			next if ( $tpdpprob < $TEQPUPDATETHR );
			
			if ( ! exists( $winteqmassprobs{$itp} ) ) {
				$winteqmassprobs{$itp} = $tpdpprob;
			}
			else {
				$winteqmassprobs{$itp} += $tpdpprob;
			}
		}
	}
	
	#Set the updates
	foreach my $itp ( keys( %winteqmassprobs ) ) {
		$teqpairs->{"INT2PROB"}->{$itp}->{"PROB"} = $winteqmassprobs{$itp};
	} #end all pairs

	print( STDERR "emacc2::updatePLexTable: elapsed time: " . timerStop( \$tmref ) . ".\n" );
}

#Ok mem.
#memory efficient ready (ver 2.1).
sub probPLexAlign( $$$ ) {
	my( $den, $dro, $plextable ) = @_;
	my( $tpfordp ) = getAllTeq4Dp( $plextable, $den, $dro );
	my( $sum, $R ) = ( 0, 0 );
		
	foreach my $itp ( keys( %{ $tpfordp } ) ) {
		$sum += $tpfordp->{$itp};
		$R++;
	}
		
	#0 div guard.
	return $DEFSMALLP if ( $sum == 0 || $R == 0 );
	#Number of TEQ pairs in this document pair
	return $sum / $R;
}

####################################### End EM optimization #######################################################

sub readDocList1( $$ ) {
	my( $srcfile, $trgfile ) = @_;
	my( %docs ) = ( 
		$SRCL => { 
			"u" => {},
		},
		$TRGL => {
			"u" => {},
		},
		"ALIGNGS" => {
		}
	);
	
	open( SDL, "< $srcfile" ) or die( "emacc-nogs::readDocList1: cannot open file '$srcfile' !\n" );
	
	while ( my $line = <SDL> ) {
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		my( $doc ) = $line;
		
		if ( ! exists( $docs{$SRCL}->{"u"}->{$doc} ) ) {
			$docs{$SRCL}->{"u"}->{$doc} = 1;
		}
		else {
			warn( "emacc2::readDocList1: duplicate doc '$doc' !\n" );
		}
	}
	
	close( SDL );
	
	open( TDL, "< $trgfile" ) or die( "emacc2::readDocList1: cannot open file '$trgfile' !\n" );
	
	while ( my $line = <TDL> ) {
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;

		my( $doc ) = $line;
		
		if ( ! exists( $docs{$TRGL}->{"u"}->{$doc} ) ) {
			$docs{$TRGL}->{"u"}->{$doc} = 1;
		}
		else {
			warn( "emacc2::readDocList1: duplicate doc '$doc' !\n" );
		}
	}
	
	close( TDL );
	
	return %docs;
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

############################### Precompute Part that was emacc-precompute-p.pl #####################################
sub emaccPrecompP( $$ ) {
	my( %INDOCS ) = readDocLists2( $_[0] , $_[1] );
	#SRC
	my( @DOCSEN ) = @{ $INDOCS{"_${SRCL}_"} };
	#TRG
	my( @DOCSRO ) = @{ $INDOCS{"_${TRGL}_"} };
	my( %clusterinfo ) = readCluster( $CLUSTERFILE );
	my( $totalcpuno ) = scalar( keys( %clusterinfo ) );
	my( @cluster ) = keys( %clusterinfo );
	my( %checkoutfiles ) = ();
	my( @docpairbatches ) = distributeEvenly( \@DOCSEN, \@DOCSRO, $totalcpuno );

	for ( my $cnt = 0; $cnt < scalar( @docpairbatches ); $cnt++ ) {
		my( $dpbatchfile, $loadfactor ) = @{ $docpairbatches[$cnt] };
	
		#Here we run computing on the cluster ...
		my( $mach ) = shift( @cluster );
		my( $inoutfbname ) = 
			$loadfactor . "-" .
			$cnt . "-" .
			$mach . "-" .
			$SRCL . "-" .
			$TRGL . "-" .
			$LEXALSCORE . "-" .
			$PLEXTYPE;

		print( STDERR "emacc2::emaccPrecompP: starting worker '$inoutfbname' on machine '$mach' ...\n" );
		$checkoutfiles{$inoutfbname} = 1;
		runInParallel( $dpbatchfile, $mach, $clusterinfo{$mach}, $inoutfbname );
		emaccconf::portableRemoveFile( $dpbatchfile );
	}

	open( PCF, ">", $PCFILEOUT ) or die( "emacc2::emaccPrecompP: cannot open PC file '$PCFILEOUT' !\n" );
	open( TEQ, ">", $TEQFILEOUT ) or die( "emacc2::emaccPrecompP: cannot open F file '$TEQFILEOUT' !\n" );

	#Collect the results...
	#Here we must do a simple concatenation of the results...
	while ( scalar( keys( %checkoutfiles ) ) > 0 ) {
		foreach my $io ( keys( %checkoutfiles ) ) {
			my( $infile ) = File::Spec->catfile( $NFSMOUNTPOINT, $io . ".in" );
			my( $errfile ) = File::Spec->catfile( $NFSMOUNTPOINT, $io . ".err" );
			my( $readyfile ) = File::Spec->catfile( $NFSMOUNTPOINT, $io . ".ready" );
			my( $docalfile ) = File::Spec->catfile( $NFSMOUNTPOINT, $io . "-dpp.out" );
			my( $teqpfile ) = File::Spec->catfile( $NFSMOUNTPOINT, $io . "-teq.out" );
	
			if ( -f( $readyfile ) ) {
				open( PCFW, "< $docalfile" ) or die( "emacc2::emaccPrecompP: cannot open '$docalfile' file!\n" );
					while ( my $line = <PCFW> ) {
						print( PCF $line );
					}
				close( PCFW );

				open( TEQW, "< $teqpfile" ) or die( "emacc2::emaccPrecompP: cannot open '$teqpfile' file!\n" );
					while ( my $line = <TEQW> ) {
						print( TEQ $line );
					}
				close( TEQW );
		
				emaccconf::portableRemoveFile( $infile );
				emaccconf::portableRemoveFile( $readyfile );
				emaccconf::portableRemoveFile( $docalfile );
				emaccconf::portableRemoveFile( $teqpfile );
				emaccconf::portableRemoveFile( $errfile );
			
				delete( $checkoutfiles{$io} );
			}
		}
	}

	close( TEQ );
	close( PCF );
}

sub runInParallel( $$$$ ) {
	my( $docpairsbatchfile ) = $_[0];
	my( $mach, $machip ) = ( $_[1], $_[2] );
	my( $inoutbasename ) = $_[3];
	my( $infile ) = File::Spec->catfile( $NFSMOUNTPOINT, $inoutbasename . ".in" );

	#1. Write output for worker...
	open( IN, ">", $infile ) or die( "emacc2::runInParallel: cannot open file '" . $infile . "' !\n" );

	#Random code
	print( IN $inoutbasename . "\n" );

	#Print parameters from the command line:
	foreach my $p ( keys( %{ $emaccconfh } ) ) {
		print( IN "--param $p" . "=" . $emaccconfh->{$p} . "\n" );
	}
	
	#2. Open the batch file and write the pairs...
	open( BAT, "<", $docpairsbatchfile ) or die( "emacc2::runInParallel: cannot open file '" . $docpairsbatchfile . "' !\n" );
	binmode( BAT, ":utf8" );
	
	while( my $line = <BAT> ) {
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		my( $den, $dro, $dpsizekb ) = split( /\t+/, $line );
		
		print( IN $den . "\t" . $dro . "\n" );
	}
	
	close( BAT );
	close( IN );

	#Do work
	my( $thishostname ) = hostname();

	if ( $mach !~ /^${thishostname}/ ) {
		#Linux only. Sorry.
		system( "ssh rion\@$machip '.\/precompworker.pl $infile' &" );
	}
	else {
		warn( "emacc2::runInParallel: executing on localhost ...\n" );
		#If we are on this host, no ssh is needed.
		emaccconf::portableForkAndDetach( "perl precompworker.pl $infile" );
	}
}

#When you have a document alignment GS with <SRC DOC>\t<TRG DOC>\t<Corpus type: p, cs, cw>, use this function.
sub readDocList2( $ ) {
	my( %docs ) = ( "_${SRCL}_" => [], "_${TRGL}_" => [] );
	my( $lcnt ) = 0;
	
	open( DL, "< $_[0]" ) or die( "emacc2::readDocList2: cannot open file '$_[0]' !\n" );
	
	while ( my $line = <DL> ) {
		$lcnt++;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		#TAB between docs!
		my( $srcd, $trgd, $ctype ) = split( /\t+/, $line );
		
		if ( ! defined( $srcd ) || ! defined( $trgd ) || ! defined( $ctype ) || $srcd eq "" || $trgd eq "" || $ctype eq "" ) {
			warn( "emacc2::readDocList2: line $lcnt is bad in '$_[0]' !\n" );
			next;
		}
		
		push( @{ $docs{"_${SRCL}_"} }, $srcd );
		push( @{ $docs{"_${TRGL}_"} }, $trgd );
		
		if ( ! exists( $docs{$srcd} ) ) {
			$docs{$srcd} = [ $trgd, $ctype ];
		}
		else {
			warn( "emacc2::readDocList2: document '$srcd' is duplicated !\n" );
		}
	}
	
	close( DL );
	return %docs;
}

sub readDocLists2( $$ ) {
	my( %docs ) = ( "_${SRCL}_" => [], "_${TRGL}_" => [] );
	my( $lcnt ) = 0;
	
	open( DLS, "< $_[0]" ) or die( "emacc2::readDocLists2: cannot open file '$_[0]' !\n" );
	
	while ( my $line = <DLS> ) {
		$lcnt++;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		my( $docpath ) = $line;
		
		if ( ! defined( $docpath ) || $docpath eq "" ) {
			warn( "emacc2::readDocLists2: line $lcnt is bad in '$_[0]' !\n" );
			next;
		}
		
		push( @{ $docs{"_${SRCL}_"} }, $docpath );
	}
	
	close( DLS );

	open( DLT, "< $_[1]" ) or die( "emacc2::readDocLists2: cannot open file '$_[1]' !\n" );
	
	while ( my $line = <DLT> ) {
		$lcnt++;
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;
		
		my( $docpath ) = $line;
		
		if ( ! defined( $docpath ) || $docpath eq "" ) {
			warn( "emacc2::readDocLists2: line $lcnt is bad in '$_[0]' !\n" );
			next;
		}
		
		push( @{ $docs{"_${TRGL}_"} }, $docpath );
	}
	
	close( DLT );
	return %docs;
}

sub readCluster( $ ) {
	my( %cluster ) = ();

	open( CLST, "< $_[0]" ) or die( "emacc2::readCluster: cannot open file '$_[0]' !\n" );

	while ( my $line = <CLST> ) {
		$line =~ s/^\s+//;
		$line =~ s/\s+$//;

		next if ( $line =~ /^#/ );
		next if ( $line =~ /^$/ );

		my( $hostname, $ip, $cpuid ) = split( /\s+/, $line );

		$cluster{$hostname . "." . $cpuid} = $ip;
    }

	close( CLST );
	return %cluster;
}

sub distributeEvenly( $$$ ) {
	my( $docsen, $docsro, $totalcpu ) = @_;
	my( @dpbatches ) = ();
	my( $entotalsizekb ) = 0;
	my( $batchfnbase ) = File::Spec->catfile( $LOCALMOUNTPOINT, "$SRCL-$TRGL-batch" );
	
	foreach my $den ( @{ $docsen } ) {
		$entotalsizekb += ( ( -s $den ) / 1024 );
	}
	
	my( $crtcpu ) = 1;
	my( $crtbatchensizekb ) = 0;
	my( $crtbatchenrosizedp ) = 0;
	my( $batchfile ) = $batchfnbase . "-$crtcpu.dp";

	warn( "emacc2::distributeEvenly[$totalcpu CPUs]: distributing document pairs and writing to '$batchfile'...\n" );
	open( PAIRS, ">", $batchfile ) or die( "emacc-precomp-p::distributeEvenly: cannot open file '$batchfile' !\n" );
	binmode( PAIRS, ":utf8" );

	foreach my $den ( @{ $docsen } ) {
		my( $densizekb ) = ( -s $den ) / 1024;
		
		$crtbatchensizekb += $densizekb;

		if ( $crtbatchenrosizedp > 0 && $crtbatchensizekb > ( $entotalsizekb / $totalcpu ) && $crtcpu < $totalcpu ) {
			push( @dpbatches, [ $batchfile, $crtbatchenrosizedp ] );
			close( PAIRS );
			
			$crtbatchensizekb = $densizekb;
			$crtbatchenrosizedp = 0;
			$crtcpu++;
			$batchfile = $batchfnbase . "-$crtcpu.dp";
				
			warn( "emacc2::distributeEvenly[$totalcpu CPUs]: distributing document pairs and writing to '$batchfile'...\n" );
			open( PAIRS, ">", $batchfile ) or die( "emacc-precomp-p::distributeEvenly: cannot open file '$batchfile' !\n" );
			binmode( PAIRS, ":utf8" );
		}
		
		foreach my $dro ( @{ $docsro } ) {
			my( $psizekb ) = ( ( -s $den ) / 1024 ) * ( ( -s $dro ) / 1024 );
			
			print( PAIRS $den . "\t" . $dro . "\t" . $psizekb . "\n" );
			$crtbatchenrosizedp++;
		} #end all Romanian documents (aka target)
	} #end all English documents (aka source)
	
	push( @dpbatches, [ $batchfile, $crtbatchenrosizedp ] );
	close( PAIRS );
	
	return @dpbatches;
}

############################### End Precompute Part ################################################################

sub readCmdLineArguments( @ ) {
	my( @args ) = @_;
	my( %clconf ) = ();
	my( %allowedparams ) = (
		"EMLOOPS" => 1,
		"EMACCMODE" => 1,
		"MAXTARGETALIGNMENTS" => 1,
		"INIDISTRIB" => 1,
		"TEQPUPDATETHR" => 1,
		"LEXALSCORE" => 1,
		"PROBTYPE" => 1,
		"DICTINVERSE" => 1,
		"CLUSTERFILE" => 1
	);
	my( $input ) = 0;

	while ( scalar( @args ) > 0 ) {
		my( $opt ) = shift( @args );

		SWOPT: {
			$opt eq "--source" and do {
				$clconf{"SRCL"} = normalizeLang( shift( @args ) );
				last;
			};

			$opt eq "--target" and do {
				$clconf{"TRGL"} = normalizeLang( shift( @args ) );
				last;
			};

			$opt eq "--param" and do {
				my( $param, $value ) = split( /\s*=\s*/, shift( @args ) );

				$param = uc( $param );
				
				die( "emacc2::readCmdLineArguments: unknown parameter '$param' !\n" )
					if ( ! exists( $allowedparams{$param} ) );

				$clconf{$param} = $value;
				last;
			};

			$opt eq "--input" and do {
				$input++;
				$clconf{"INPUTFILE$input"} = shift( @args );
				last;
			};

			$opt eq "--output" and do {
				$clconf{"DALOUTFILE"} = shift( @args );
				last;
			};

			die( "emacc2::readCmdLineArguments: unknown option '$opt' !\n" );
        }
	}

	return \%clconf;
}

sub normalizeLang( $ ) {
	my( $lang ) = lc( $_[0] );
	my( %accuratlanguages ) = (
		#1
		"romanian" => "ro",
		"rum" => "ro",
		"ron" => "ro",
		"ro" => "ro",
		#2
		"english" => "en",
		"eng" => "en",
		"en" => "en",
		#3
		"estonian" => "et",
		"est" => "et",
		"et" => "et",
		#4
		"german" => "de",
		"ger" => "de",
		"deu" => "de",
		"de" => "de",
		#5
		"greek" => "el",
		"gre" => "el",
		"ell" => "el",
		"el" => "el",
		#6
		"croatian" => "hr",
		"hrv" => "hr",
		"hr" => "hr",
		#7
		"latvian" => "lv",
		"lav" => "lv",
		"lv" => "lv",
		#8
		"lithuanian" => "lt",
		"lit" => "lt",
		"lt" => "lt",
		#9
		"slovenian" => "sl",
		"slv" => "sl",
		"sl" => "sl"
	);

	return $accuratlanguages{$lang} if ( exists( $accuratlanguages{$lang} ) );
	die( "emacc2::normalizeLang: unknown language '$lang' !\n" );
}
