#!/usr/bin/env perl

use strict;
use Getopt::Std;
use Storable qw(dclone);

$Getopt::Std::STANDARD_HELP_VERSION = 1;
our $VERSION = "1.0";

my $mode = 'S';
my $cMin = 4;

sub HELP_MESSAGE {
    print "\n
Description : Generates structures along the direct path between s and t
              that contain canonical stems from the exclusive base pairs
              of both structures.

Options     : -s string      : dot-bracket encoding of the start structure
              -t string      : dot-bracket encoding of the target structure
              -m char    [$mode] : generation mode (S)ingle canonical helix per structure
                               or (C)ombinations of canonical helices
              -l integer [$cMin] : minimal canonical helix length
              -h         : this help

Version     : $VERSION\n"
."";
};

# parsing
our %args;
getopt("hs:t:m:l:", \%args);

# help handling
if ( defined $args{h} ) { HELP_MESSAGE(); exit(0);}

# sanity check of dot-bracket structure 
# arg1 = structure string encoding
# arg2 = name of structures
sub checkStructure {
	my($p, $name) = @_; 
	# check dot-bracket encoding
	die ("$name structure is not in dot-bracket format") unless ($p =~ /^[\.\(\)]+$/);
	# check nested structure encoding
	my @pi = split //, $p;
	my @openBP = (); 
	for (my $i=0; $i<(scalar @pi); $i++) {
		if( $pi[$i] eq '(' ) { 
			# open bp
			push @openBP, $i; 
			$pi[$i] = $i;
		} elsif ( $pi[$i] eq ')' ) {
			# close bp 
			die("$name structure not nested") if (scalar @openBP == 0);
			my $lastOpen = pop @openBP;
			$pi[$lastOpen] = $i;
			$pi[$i] = $lastOpen; 
		} else {
			# unpaired
			$pi[$i] = -1;
		}
	}
	die("$name structure not nested") unless (scalar @openBP == 0);
	# return base pair index information
	return @pi;
};


# parse and check input
die ("start structure missing") if (! defined $args{s});
my $p0 = $args{s}; $p0 =~ s/^\s+|\s+$//g;
my @bp0 = checkStructure($p0, "start");
my $pn = $args{t}; $pn =~ s/^\s+|\s+$//g;
my @bpn = checkStructure($pn, "target");
die("start and target structure have different lengths") unless ( scalar @bp0 == scalar @bpn );
$mode = ($args{m}) if (defined $args{m});
die("mode has to be 'S'") unless ($mode =~ /^[SsCc]$/);
$cMin = ($args{l}) if (defined $args{l});
die("minimal canonical helix length has to be integer") unless ($cMin =~ /^\d+$/);
die("minimal canonical helix length has to be > 0") unless ($cMin > 0);

# identify exclusive base pairs
my $commonBpShift = -100;
for (my $i=0; $i < scalar @bp0; $i++) {
	if ($bp0[$i] >= 0 and $bp0[$i] == $bpn[$i]) {
		$bp0[$i] = $commonBpShift - $bp0[$i]; # mark as common base pairs from p0
		$bpn[$i] = $commonBpShift - $bpn[$i]; # mark as common base pairs from pn
	}
}

# get shared base pairs
my @bpShared = map {$_ < 0 ? $_ : -1} @bp0;

# generates dot-bracket string from bp index encoding
sub bp2p {
	my @bp = @_;
	my $db = "";
	for (my $i=0; $i< scalar @bp; $i++) {
		if ($bp[$i] == -1) { $db .= "."; }
		elsif ( $bp[$i] >= 0 ) {
			$db .= $bp[$i] > $i ? "(" : ")";
		} else { 
			$db .= $commonBpShift-$bp[$i] > $i ? "(" : ")";
		} 
	}
	return $db;
}

#print "start     = $p0\n";
#print "target    = $pn\n";
#print "shared    = ".(bp2p(@bpShared))."\n";
#print "start'    = ".(bp2p(@bp0))."\n";
#print "target'   = ".(bp2p(@bpn))."\n";
#print "start ex  = ".bp2p(map { $_ < 0 ? -1 : $_; } @bp0)."\n";
#print "target ex = ".bp2p(map { $_ < 0 ? -1 : $_; } @bpn)."\n";

# compute canonical stems (of minimal length) within exclusive base pairs
sub countCanonicalStackLength {
	my @bp = @_;
	my @csl = ();
	
	for (my $i=0; $i< scalar @bp; $i++) {
		if ($bp[$i] >= 0 and $bp[$i] > $i) {
			# find stacking length
			my $l = 1;
			while( $i+$l < $bp[$i] and $bp[$i+$l] == $bp[$i]-$l) {
				$csl[$i+$l] = 0;
				$csl[$bp[$i+$l]] = 0;
				$l++;
			}
			$csl[$i] = ($l >= $cMin) ? $l : 0;
			# set i to last stacked bp
			$i = $i+$l-1;
		} else {
			$csl[$i] = 0;
		}
	}
	
	return @csl;
}

# compile output structure ######################################

sub singleExclusiveStem {
	my @bp = @_;
	my @pList = ();
	
	my @csl = countCanonicalStackLength(@bp);
	#print join(",",@csl)."\n";
	
	for (my $i=0; $i<scalar @csl; $i++) {
		if ($csl[$i] > 0) {
			# copy shared base pairs into new array (deep copy needed)
			my $bpNew = dclone(\@bpShared);
			# add current canonical stem
			for (my $l=0; $l< $csl[$i]; $l++) {
				$bpNew->[$i+$l] = $bp[$i+$l];
				$bpNew->[ $bp[$i+$l] ] = $i+$l;
			}
			push @pList, bp2p(@$bpNew);
		}
	}
	
	return @pList;
};

# copied from List::PowerSet https://metacpan.org/source/NIKC/List-PowerSet-0.01/lib/List/PowerSet.pm
sub powerset {
  return [[]] if @_ == 0;
  my $first = shift;
  my $pow = &powerset;
  [ map { [$first, @$_ ], [ @$_] } @$pow ];
}

# check if helices are compatible
sub compatibleHelices {
	my @cs = @_;
	for my $i (@cs) {
		for my $j (@cs) {
			# check if from different structures
			if ($i != $j and $i * $j < 0) {
				# check if helices are overlapping
				my $l0 = ( $i > 0 ? $i : $j) -1;
				my $ln = abs( $i < 0 ? $i : $j) -1;
				if (   ($l0 <= $ln and $ln <= $bp0[$l0]) 
					or ($l0 <= $bpn[$ln] and $bpn[$ln] <= $bp0[$l0]) 
					or ($ln <= $l0 and $l0 <= $bpn[$ln]) 
					or ($ln <= $bp0[$l0] and $bp0[$l0] <= $bpn[$ln]) )
				{
					return 0;
				} 
			}
		}
	}
	return 1;
}

sub combinedExclusiveStems {
	my @pList = ();
	
	my @csl0 = countCanonicalStackLength(@bp0);
	my @csln = countCanonicalStackLength(@bpn);
	
	my @csIdx = ();
	for (my $i=0; $i<scalar @csl0; $i++) {
		if ($csl0[$i] > 0) { push @csIdx, $i+1; } # store 1-based index for p0
		if ($csln[$i] > 0) { push @csIdx, -$i-1; } # store negated 1-based index for pn
	}
	
	my $idxCombi = powerset( @csIdx );
	
	for my $subset (@$idxCombi) {
		next if (scalar @$subset == 0); # skip empty addition since already in list
		next if (compatibleHelices( @$subset ) == 0);
		# copy shared base pairs into new array (deep copy needed)
		my $bpNew = dclone(\@bpShared);
		# add helices
		for my $i (@$subset) {
			if ($i > 0) { # helix from p0
				for (my $l = 0; $l < $csl0[$i-1]; $l++) {
					$bpNew->[$i-1+$l] = $bp0[$i-1+$l];
					$bpNew->[ $bp0[$i-1+$l] ] = $i-1+$l;
				}
			} else { # helix from pn
				for (my $l = 0; $l < $csln[-$i-1]; $l++) {
					$bpNew->[-$i-1+$l] = $bpn[-$i-1+$l];
					$bpNew->[ $bpn[-$i-1+$l] ] = -$i-1+$l;
				}
			}
		}
		my $pNew = bp2p(@$bpNew);
		# store if not start/target structure since already in the list
		if($pNew ne $p0 and $pNew ne $pn) {
			push @pList, $pNew;
		}
	}
	
	
	
	return @pList;
}


# create list of final structures
my @pFinal = ();
push @pFinal, $p0, $pn, bp2p(@bpShared);

# single stem output mode (linear number of structures)
if ($mode =~ /^[Ss]$/ ) {
	
	push @pFinal, singleExclusiveStem( @bp0 );
	push @pFinal, singleExclusiveStem( @bpn );
	
} elsif ($mode =~ /^[Cc]$/ ) {
	
	push @pFinal, combinedExclusiveStems();
	
} 

# print final list of structures
print join "\n", @pFinal;



