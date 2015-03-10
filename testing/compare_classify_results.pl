#!/usr/bin/env perl


############################# MNI Header #####################################
#@NAME       :  compare_classify_result.pl
#@DESCRIPTION:  check if the classify results is consistent
#@COPYRIGHT  :
#              Vladimir S. Fonov  2009
#              Montreal Neurological Institute, McGill University.
#              Permission to use, copy, modify, and distribute this
#              software and its documentation for any purpose and without
#              fee is hereby granted, provided that the above copyright
#              notice appear in all copies.  The author and McGill University
#              make no representations about the suitability of this
#              software for any purpose.  It is provided "as is" without
#              express or implied warranty.
###############################################################################


use strict;
use Getopt::Long;
use File::Basename;
use File::Temp qw/ tempdir /;

my $verbose=0;
my $clobber=1;
my $fake=0;
my $mask;
my $me = basename ($0);
my $method='ann';
my $priors;
my $tag_file= dirname($0)."/tags.tag";

GetOptions(
      'verbose'           => \$verbose,
      'clobber'           => \$clobber,
      'method=s'          => \$method,
      'tags=s'            => \$tag_file,
      'priors=s'          => \$priors
     );

my $Help = <<HELP;
  Usage: $me <input.mnc> <reference.mnc> 
    --verbose be verbose
    --method <ann,knn,min,bayes>
    --priors a,b,c
  Problems or comments should be sent to: vladimir.fonov\@gmail.com
HELP

die $Help if $#ARGV < 1;

my ($in,$ref)=@ARGV;

my $tol=0.00001;#every voxel counts

my $tmpdir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

my @args=("classify", '-tagfile',$tag_file,'-'.$method,'-nocache');

push @args,"-apriori", "-volume",$priors if $priors;

push @args,$in, "$tmpdir/cls_$method.mnc";

do_cmd(@args);

do_cmd('mincmath','-sub','-float',"$tmpdir/cls_$method.mnc" , $ref, "$tmpdir/diff.mnc",'-q');

my $mean=`mincstats -q -mean $ref`;
my $count=`mincstats -q -count $tmpdir/diff.mnc`;
my $sum2=`mincstats -q -sum2 $tmpdir/diff.mnc`;

chomp($mean);chomp($count);chomp($sum2);
print "mean=$mean count=$count sum2=$sum2\n" if $verbose;

my $rms=sqrt($sum2/$count);

my $rms_pc=$rms/$mean;

if($rms_pc>$tol) {
  die "relative RMS difference: $rms_pc larger then $tol\n";
} 

sub do_cmd { 
    print STDOUT "@_\n" if $verbose;
    if(!$fake){
      system(@_) == 0 or die "DIED: @_\n";
    }
}
