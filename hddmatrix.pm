package hddmatrix;

use strict;
use warnings;
use File::Spec;
use Fcntl qw( :seek );

#Classic matrix with m[1,2]
sub newInt;
#Matrix with keys like m['red', 'white']
sub newKey;
sub write;
sub read;
sub genPTypeValue;
sub DESTROY;

sub newInt {
	my( $classname, $m, $n, $packtype ) = @_;
	my( $this ) = {};
	
	$this->{"M"} = $m;
	$this->{"N"} = $n;
	$this->{"PTYPE"} = $packtype;
	
	my( $matrixfile ) =  File::Spec->catfile( File::Spec->tmpdir(), genRandomFileName() );

	open( my $mfhandle, "+>", $matrixfile ) or die( "hddmatrix::new: cannot open file '$matrixfile' !\n" );
	binmode( $mfhandle );
	
	$this->{"MFILE"} = $matrixfile;
	$this->{"HANDLE"} = $mfhandle;

	bless( $this, $classname );
	return $this;
}

#Pointers to key vectors for M and N dimensions.
sub newKey {
	my( $classname, $mkeys, $nkeys, $packtype ) = @_;
	my( $this ) = {};
	
	$this->{"M"} = scalar( @{ $mkeys } );
	$this->{"N"} = scalar( @{ $nkeys } );
	$this->{"MK2I"} = { "#" => 0 };
	$this->{"NK2I"} = { "#" => 0 };
	$this->{"PTYPE"} = $packtype;
	
	#Convert keys to integers M...
	foreach my $k ( @{ $mkeys } ) {
		next if ( $k eq "#" );
		
		if ( ! exists( $this->{"MK2I"}->{$k} ) ) {
			$this->{"MK2I"}->{$k} = $this->{"MK2I"}->{"#"};
			$this->{"MK2I"}->{"#"}++;
		}
		else {
			warn( "hddmatrix::newKey[M]: key '$k' is duplicated!...\n" );
		}
	}
	
	#Convert keys to integers M...
	foreach my $k ( @{ $nkeys } ) {
		next if ( $k eq "#" );
		
		if ( ! exists( $this->{"NK2I"}->{$k} ) ) {
			$this->{"NK2I"}->{$k} = $this->{"NK2I"}->{"#"};
			$this->{"NK2I"}->{"#"}++;
		}
		else {
			warn( "hddmatrix::newKey[N]: key '$k' is duplicated!...\n" );
		}
	}
	
	my( $matrixfile ) =  File::Spec->catfile( File::Spec->tmpdir(), genRandomFileName() );

	open( my $mfhandle, "+>", $matrixfile ) or die( "hddmatrix::new: cannot open file '$matrixfile' !\n" );
	binmode( $mfhandle );
	
	$this->{"MFILE"} = $matrixfile;
	$this->{"HANDLE"} = $mfhandle;

	bless( $this, $classname );
	return $this;
}

sub DESTROY {
	my( $this ) = shift();
	
	#Close and delete the file.
	close( $this->{"HANDLE"} );
	unlink( $this->{"MFILE"} );
}

sub getRowKeys {
	my( $this ) = shift();
	my( $rowkeyswithcount ) = $this->{"MK2I"};
	
	delete( $rowkeyswithcount->{"#"} );
	
	my( @rowkeys ) = keys( %{ $rowkeyswithcount } );
	
	return \@rowkeys;
}

sub getColumnKeys {
	my( $this ) = shift();
	my( $colkeyswithcount ) = $this->{"NK2I"};
	
	delete( $colkeyswithcount->{"#"} );
	
	my( @colkeys ) = keys( %{ $colkeyswithcount } );
	
	return \@colkeys;
}

sub addRowKey {
	my( $this ) = shift();
	my( $key ) = @_;
	
	return if ( $key eq "#" );
	
	if ( ! exists( $this->{"MK2I"}->{$key} ) ) {
		$this->{"MK2I"}->{$key} = $this->{"MK2I"}->{"#"};
		$this->{"MK2I"}->{"#"}++;
		$this->{"M"}++;
	}
	
}

sub addColumnKey {
	my( $this ) = shift();
	my( $key ) = @_;
	
	return if ( $key eq "#" );
	
	if ( ! exists( $this->{"NK2I"}->{$key} ) ) {
		$this->{"NK2I"}->{$key} = $this->{"NK2I"}->{"#"};
		$this->{"NK2I"}->{"#"}++;
		$this->{"N"}++;
	}
}

#0-based addressing.
sub write {
	my( $this ) = shift();
	my( $i, $j, $value ) = @_;
	
	#Check if i and j are integers.
	if ( $i !~ /^[0-9]+$/ || $j !~ /[0-9]+$/ ) {
		if ( exists( $this->{"MK2I"}->{$i} ) ) {
			$i = $this->{"MK2I"}->{$i};
		}
		else {
			warn( "hddmatrix::write: 'i' key '$i' is not mapped !\n" );
			return 0;
		}

		if ( exists( $this->{"NK2I"}->{$j} ) ) {
			$j = $this->{"NK2I"}->{$j};
		}
		else {
			warn( "hddmatrix::write: 'j' key '$j' is not mapped !\n" );
			return 0;
		}		
	}
	
	my( $pvalue ) = pack( $this->{"PTYPE"}, $value );
	my( $reclen ) = length( $pvalue );
	
	do {
		warn( "hddmatrix::write: 'i' value '$i' is larger than constructed 'M' (" . $this->{"M"}  . ") !\n" );
		return;
	} if ( $i >= $this->{"M"} );

	do {
		warn( "hddmatrix::write: 'j' value '$j' is larger than constructed 'N' (" . $this->{"N"}  . ") !\n" );
		return;
	} if ( $j >= $this->{"N"} );
	
	my( $offset ) = ( $i * $this->{"N"} + $j ) * $reclen;
	
	sysseek( $this->{"HANDLE"}, $offset, SEEK_SET );
	syswrite( $this->{"HANDLE"}, $pvalue );
	return 1;
}

#0-based addressing.
sub read {
	my( $this ) = shift();
	my( $i, $j ) = @_;

	#Check if i and j are integers.
	if ( $i !~ /^[0-9]+$/ || $j !~ /[0-9]+$/ ) {
		if ( exists( $this->{"MK2I"}->{$i} ) ) {
			$i = $this->{"MK2I"}->{$i};
		}
		else {
			warn( "hddmatrix::write: 'i' key '$i' is not mapped !\n" );
			return 0;
		}

		if ( exists( $this->{"NK2I"}->{$j} ) ) {
			$j = $this->{"NK2I"}->{$j};
		}
		else {
			warn( "hddmatrix::write: 'j' key '$j' is not mapped !\n" );
			return 0;
		}		
	}	

	#Random value to see the record length
	my( $pvaluetest ) = pack( $this->{"PTYPE"}, $this->genPTypeValue() );
	my( $reclen ) = length( $pvaluetest );
	
	do {
		warn( "hddmatrix::read: 'i' value is larger than constructed 'M' (" . $this->{"M"}  . ") !\n" );
		return;
	} if ( $i >= $this->{"M"} );

	do {
		warn( "hddmatrix::read: 'j' value is larger than constructed 'N' (" . $this->{"N"}  . ") !\n" );
		return;
	} if ( $j >= $this->{"N"} );
	
	my( $offset ) = ( $i * $this->{"N"} + $j ) * $reclen;
	
	sysseek( $this->{"HANDLE"}, $offset, SEEK_SET );
	
	my( $pvalue );
	
	sysread( $this->{"HANDLE"}, $pvalue, $reclen );
	
	return unpack( $this->{"PTYPE"}, $pvalue );
}

sub genPTypeValue {
	my( $this ) = shift();
	my( $type ) = $this->{"PTYPE"};
	
	SWTYPE: {
		( $type eq "f" || $type eq "d" ) and return 1.2e-3;
		( $type eq "i" || $type eq "l" ) and return -123;
		( $type eq "I" || $type eq "L" ) and return 123;
	}
	
	die( "hddmatrix::genPTypeValue: unknown type '$type' !\n" );
}

sub genRandomFileName {
	my( $this ) = shift();
	my( @letters ) = (
		"a", "b", "c", "d", "e", "f", "g", "h",
		"i", "j", "k", "l", "m", "n", "o", "p",
		"q", "r", "s", "t", "u", "v", "x", "z",
		"w", "y",
		"A", "B", "C", "D", "E", "F", "G", "H",
		"I", "J", "K", "L", "M", "N", "O", "P",
		"Q", "R", "S", "T", "U", "V", "X", "Z",
		"W", "Y"
	);
	my( @numbers ) = ( "0", "1", "2", "3", "4", "5", "6", "7", "8", "9" );
	my( $fname ) = "";

	for ( my $i = 1; $i <= 8; $i++ ) {
		my( $lettorno ) = int( rand( 2 ) );

		#Letter
		if ( $lettorno ) {
			my( $leti ) = int( rand( scalar( @letters ) ) );

			$fname .= $letters[$leti];
		}
		#Number
		else {
			my( $noi ) = int( rand( scalar( @numbers ) ) );
			$fname .= $numbers[$noi];
		}
	}

	return $fname . ".m";
}

1;
