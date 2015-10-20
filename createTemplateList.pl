#!/usr/bin/perl

use strict;
use warnings;
use Math::Libm qw(:all);
use List::Util qw(min max);

sub fisherYatesShuffle(@);
sub getAnatomyFilenames($@);
sub appendToTemplateFile($);

if ($#ARGV + 1 < 1) {
  die "You failed to provide the name of the directory in which the original anatomicals are located! Exiting.";
}

my $vbmDir = $ARGV[0];
my $templateFile = "$vbmDir/template_list";
unlink $templateFile;

opendir(VBMDIR, $vbmDir) || die "Can't open $vbmDir\n";

my @anatFiles = grep { -f "$vbmDir/$_" && !/^\./ && /.*[_.]anat\.nii(\.*)?/ } readdir(VBMDIR);
#my @anatFiles = grep { -f "$vbmDir/$_" && !/^\./ } readdir(VBMDIR);

print "Found the following anat files:\n";
print "@anatFiles" . "\n";

my %groupContents = ();

foreach my $anatFile  ( @anatFiles ) {
#  print "Group: $group\n";

  my @parts = split(/\./, $anatFile);
  #print "@parts" . "\n";
  my $group = $parts[0];
  my $subject = $parts[1];

  if (exists $groupContents{$group}) {
    push (@{ $groupContents{$group} }, $subject);
  }
  else {
    $groupContents{$group} = [ $subject ];
  }
  
} # end of foreach my $anatFile  ( @anatFiles )

my $min = 10000;
my $max = -1;

for my $group (keys %groupContents) {
  if (@{ $groupContents{$group} } < $min) {
    $min = @{ $groupContents{$group} } ;
  }
  if (@{ $groupContents{$group} } > $max) {
    $max = @{ $groupContents{$group} } ;
  }

  #print "$group group contains " . @{ $groupContents{$group} } . " subjects: @{ $groupContents{$group} }\n";
}

closedir(VBMDIR);

if ($min ne $max) {
  print "Unequal group sizes: creating template_list file\n";

  foreach my $key (sort keys %groupContents) {
    #print "$key\n";
    print "$key group has " . @{ $groupContents{$key} } . " subjects ";
    if ( @{ $groupContents{$key} }  == $min ) {
      print " *** SMALLEST GROUP ***\n";
    }
    elsif ( @{ $groupContents{$key} }  == $max ) {
      print " *** LARGEST GROUP ***\n";
    }
    else {
      print "\n";
    }
    if ( @{ $groupContents{$key} } eq $min) {
      #print "$key: @{$groupContents{$key}}\n\n";
      my $anatFiles = getAnatomyFilenames($key, $groupContents{$key});

      # the ls command (without args) appears to list things in
      # lexicographic order, the call to sort appears to mimic this
      my @subset = sort @{$anatFiles};
      appendToTemplateFile(\@subset);
    }
    else {
      #print "BEFORE: @{$groupContents{$key}}\n";
      # now randomly permute the contents of the array for the longer groups
      fisherYatesShuffle $groupContents{$key};
      #print "AFTER: @{$groupContents{$key}}\n\n";
      my $anatFiles = getAnatomyFilenames($key, $groupContents{$key});

      # the ls command (without args) appears to list things in
      # lexicographic order, the call to sort appears to mimic this
      my @subset = sort @{$anatFiles}[0..$min-1];
      appendToTemplateFile(\@subset);
    }
  }
}
else {
  print "Group sizes are equal, no need for a template file\n";
}

################################################################################

sub appendToTemplateFile ($) {
  my $anatFiles = shift;

  open TEMPLATE, ">>$templateFile";
  foreach my $file (@$anatFiles) {
    print TEMPLATE "$file\n";
  }
  close TEMPLATE;
}

sub getAnatomyFilenames($@) {
  my $group = shift;
  my $subjects = shift;
  my @anatFiles = ();
  
  for my $subject (@$subjects) {
    my $anatFile = join(".", $group, $subject, "anat.nii.gz");

    push @anatFiles, $anatFile;
    #print "@anatFiles\n";    
  }
  return \@anatFiles;
}

# fisher_yates_shuffle( \@array ) : generate a random permutation
# of @array in place
sub fisherYatesShuffle (@) {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}

