#!xPERLx -w
 
# ------------------------------ MNI Header ----------------------------------
#@NAME       : classify_clean
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: classify a stereotaxic T1-weighted MRI volume into Grey, 
#              White and CSF classes.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : Wed Feb 19, 1997, Louis Collins
#@MODIFIED   : Thu Dec 12, 1997, Alex Zijdenbos
#@MODIFIED   : Wed, Mar 21, 2002, Jason Lerch
#@VERSION    : $Id: classify_clean.pl,v 1.3 2004-01-23 19:32:12 jason Exp $
#-----------------------------------------------------------------------------

use MNI::Startup;
use MNI::Spawn;
use MNI::DataDir;
use MNI::FileUtilities qw(check_output_dirs check_files);
use MNI::PathUtilities qw(split_path);
use Getopt::Tabular;

use strict;

# ------------------------------ start here

# User-modifiable globals
my @Infiles;
my $Outfile;
my $CleanTags      = 0;
my $ClassifyOpt    = undef;
my $CleanedTags    = undef;
my $Mask           = undef;
my $MaskBinValue   = undef;
my $MaskTag        = 0;
my $MaskClassified = 0;
my $MaskCerebellum = 0;

my $TagDir         = MNI::DataDir::dir('classify');
my $BackgroundTags = 'ntags_1000_bg.tag';
my $StandardTags   = 'ntags_1000_prob_90_nobg.tag';
my $ICBM           = MNI::DataDir::dir('ICBM');
my $CerebellumMask = 'icbm_atlas_cerebellum_avg_53_mask_1mm.mnc';
my $MergeFile      = undef;
my $MergeLut       = undef;

# Other globals
my ($Usage, $Help);
my $Version = "0.2";
my $LongVersion = "version ${Version}: untested perl code. Beware!";

&Initialize;

my $Tags = $StandardTags;

if ($CleanTags) {
    my $tmp_classified = &initial_classify_volume(@Infiles);
    $Tags              = &customize_standard_tags($tmp_classified);
}

&compute_final_classification($Outfile, $Tags, @Infiles);

#-----------------------------------------------------------
sub customize_standard_tags {
  my ($classified) = @_ ;

  my ($indir, $filebase) = (split_path($classified))[0,1];

  my $cleanedtags = "${TmpDir}${filebase}_cleaned.tag";
  my $tags = defined($CleanedTags) ? $CleanedTags : "${TmpDir}${filebase}_custom.tag";

  my $cl = "${indir}/${filebase}";
  
  if (-e $tags && !$Clobber) {
      print "file $tags already exists\n";
  }
  else {
      unlink $cleanedtags if (-e $cleanedtags);

      # Figure out which fuzzy classifications were created
      my @fuzzy_files = glob("${cl}_*.mnc");
      
      # Build base cleantag command
      my($select) = "-mode 110 -threshold 0.7 -difference 0.3";
      my(@CleanTagCmd) = ('cleantag', '-oldtag', $StandardTags, '-newtag', $cleanedtags,
			  split(' ', $select), '-comment', $select);

      # Append fuzzy files and class IDs to cleantag command.
      # Skip class 0 (background)
      my($ffile);
      foreach $ffile (@fuzzy_files) {
	  if (($ffile =~ /\_(\d+)\.mnc/) && $1) {
	      push(@CleanTagCmd, ($ffile, $1));
	  }
      }

      # Run cleantag
      Spawn(\@CleanTagCmd);
      
      if ($MaskTag) {
	  $cleanedtags = &mask_tags($cleanedtags, $Mask, 1, 
				    "${TmpDir}masked_custom.tag");
      } 

      open (CL_TAG, $cleanedtags) ||
	  die "Error: can't open $cleanedtags\n$!\n";
      open (ST_TAG, $BackgroundTags) ||
	  die "Error: can't open $BackgroundTags\n$!\n";
      open (NW_TAG, ">$tags") ||
	  die "Error: can't open $tags\n$!\n";
      
      my $clean;
      while (defined($clean = <CL_TAG>)) {

	  print NW_TAG $clean;        # print line from cleanedtags
	  
	  if ($clean =~ /Points =/) {       # insert points from standard bkg tags
	      my $stand;
	      while (defined($stand = <ST_TAG>)) {
		  if ($stand =~ /^\s*[\d\.-]+\s+[\d\.-]+\s+[\d\.-]+/) {
		      $stand =~ s/;+$//;
		      print NW_TAG $stand;
		  }
	      }
	  }
      }
      close (CL_TAG);
      close (ST_TAG);
      close (NW_TAG);
  }
  
  return $tags;
}

#-----------------------------------------------------------
sub compute_final_classification {
    my @infiles    = @_;
    my $classified = shift @infiles;
    my $tags       = shift @infiles;
    my $filebase   = (split_path($infiles[0]))[1];

    # Just to make sure; this should never happen
    die "ERROR: output file $classified exists!\n" if (-e $classified);

    my $output = (defined $MergeFile) ? "${TmpDir}temp_classified_1.mnc" : $classified;
    
    Spawn(['classify', '-ann',  '-tagfile', $tags, @infiles, $output]);
#         "-fuzzy '011111' -fpath ${TmpDir}/ -fprefix ${filebase}_final ".

    if (defined $MergeFile) {
	my $mergefile = $MergeFile;
	if (defined $MergeLut) {
	    # Remap mergefile through LUT
	    my $temp = "${TmpDir}mergefile_remapped.mnc";
	    Spawn(['minclookup', '-byte', '-valid_range', 0, 255, '-discrete', 
		   '-lut_string', $MergeLut, $mergefile, $temp]);
	    $mergefile = $temp;
	}
	
	# Create BG mask of merge file. Take an interval around 0 to
	# avoid numerical errors in selecting 0.
	my $mergefile_bg = "${TmpDir}mergefile_bg.mnc";
	Spawn(['mincmath', '-byte', '-const2', -0.5, 0.5, '-segment', 
	       $mergefile, $mergefile_bg]);

	# Mask classified volume using BG of merge file
	my $temp_classified = "${TmpDir}temp_classified_2.mnc";
	Spawn(['mincmath', '-byte', '-range', 0, 255, '-mult', 
	       $output, $mergefile_bg, $temp_classified]);

	# Add merge file to classified volume
	Spawn(['mincmath', '-byte', '-range', 0, 255, '-add', 
	       $temp_classified, $mergefile, $classified]);
    }
}

#-----------------------------------------------------------
sub initial_classify_volume {
  my (@infiles) = @_ ;

  my $filebase = (split_path($infiles[0]))[1];
  
  my $classified = "${TmpDir}${filebase}_fuzzy.mnc";

  if (-e $classified) {
      print "file $classified already exists\n";
  }
  else {
      Spawn(['classify', '-min', '-tagfile', $StandardTags, '-fuzzy', 'all', 
	     '-fpath', $TmpDir, '-fprefix', "${filebase}_fuzzy", @Infiles,
	     $classified]);
  }
  
  return $classified;
}

#-----------------------------------------------------------
sub mask_tags {
    my($tagfile, $mask, $maskFG, $newtagfile) = @_;

    Spawn(['cleantag', '-oldtag', $tagfile, '-newtag', $newtagfile, '-mask', $mask,
	   '-maskbinvalue', $maskFG]);
    
    return $newtagfile;
}


#-----------------------------------------------------------
sub Initialize
{
   $Clobber     = 0;
   $Execute     = 1;
   $Verbose     = 1;

   &SetHelp;
   
   my @ArgInfo =
       (@DefaultArgs,
	["Specific options", "section"],
	["-version", "call", undef, \&print_version, 
	 "print version and quit"],
	["-clean_tags|-noclean_tags", "boolean", 1, \$CleanTags,
	 "clean tag file using mindist pre-classification [default is -noclean_tags]"],
	["-tagfile", "string", 1, \$StandardTags,
	 "specify the stereotaxic tag file [default: $TagDir$StandardTags]", "<tags.tag>"],
	["-cleanedtagfile", "string", 1, \$CleanedTags,
	 "save the cleaned tags in <tags.tag>", "<tags.tag>"],
	["-classify", "string", 1, \$ClassifyOpt,
	 "specify options to pass directly to classify (quoted, comma- or '::'- separated list)", "'<options>'"],
	
	["Mask options", "section"],
	["-mask", "string", 1, \$Mask,
	 "specify a mask volume", "<mask.mnc>"],
	["-maskbinvalue", "integer", 1, \$MaskBinValue,
	 "value of mask foreground [default: assume a binary mask]", "<int>"],
	["-mask_tag|-nomask_tag", "boolean", 1, \$MaskTag,
	 "apply mask to foreground tags prior to classification(s) [default: -nomask_tag]"],
	["-mask_classified|-nomask_classified", "boolean", 1, \$MaskClassified,
	 "apply mask to intermediate and final classified volumes [default: -nomask_classified]"],
	["-mask_cerebellum|-nomask_cerebellum", "boolean", 1, \$MaskCerebellum,
	 "remove cerebellar tags from standard set prior to classification"],
	["Merge options", "section"],
	["-merge", "string", 1, \$MergeFile,
	 "merge the argument file into the final classification. This implies first masking the classification with the background (0) of this file, and subsequently adding it in. Volume <merge.mnc> must be discrete (a label volume)", "<merge.mnc>"],
	["-merge_lut", "string", 1, \$MergeLut,
	 "a quoted, semicolon-separated lookup table, used to remap the values in <merge.mnc> prior to merging it with the final classification. The key/value pairs can be separated by spaces, commas, or '::'. E.g., \"-merge_lut '1 4; 2,5'\".", "<lut>"],
	);
   
   my (@argv) = @ARGV;
   GetOptions (\@ArgInfo, \@argv) || die "\n";
   
   if (@argv < 2)  {
       warn $Usage;
       die "Incorrect number of arguments\n";
   }
   
   if ($Execute) {
       check_output_dirs($TmpDir);
       if (!$ENV{'TMPDIR'}) {
	   $ENV{'TMPDIR'} = $TmpDir;
       }
   }
   
   MNI::Spawn::SetOptions (strict => 2);
   
   # Register required programs
   RegisterPrograms([qw(classify cleantag minclookup mincmath)]) || die;
   
   # They were found, so add options according to
   # $Debug, $Verbose, $Clobber flags
   AddDefaultArgs([qw(classify minclookup mincmath)], 
		  ($Verbose) ? '-verbose' : '-quiet');
   AddDefaultArgs([qw(classify minclookup mincmath)], '-clobber') if $Clobber;
   AddDefaultArgs('classify', '-debug 4') if $Debug;

   if (defined $ClassifyOpt) {
       $ClassifyOpt =~ s/,/ /g;
       $ClassifyOpt =~ s/::/ /g;

       AddDefaultArgs('classify', $ClassifyOpt);
   }

   @Infiles = @argv;
   $Outfile = pop(@Infiles);

   # Verify input files
   check_files(@Infiles) || die;
   
   # Verify output file
   if (-e $Outfile) {
       if (!$Clobber) {
	   die "$Outfile already exists (use -clobber to overwrite)\n";
       }
       else {
	   unlink $Outfile;
       }
   }

   # Verify tag files
   MNI::DataDir::check_data($TagDir, [$BackgroundTags, $StandardTags]);
   $BackgroundTags = "${TagDir}${BackgroundTags}";
   $StandardTags   = "${TagDir}${StandardTags}";
   
   # Verify cerebellar mask and mask standard tags
   if ($MaskCerebellum) {
       MNI::DataDir::check_data($ICBM, [$CerebellumMask]);
       $CerebellumMask = "${ICBM}${CerebellumMask}";

       $StandardTags = &mask_tags($StandardTags, $CerebellumMask, 0,
				  "${TmpDir}standard_nocereb.tag");
   }       
   
   if (defined($Mask)) {
       if (!$MaskTag && !$MaskClassified) {
	   die "Mask ${Mask} not used;\n  please specify -mask_tag and/or -mask_classify\n";
       }
       
       if (defined($MaskBinValue)) {
	   # Convert mask to binary volume
	   my $minval = $MaskBinValue - 0.5;
	   my $maxval = $MaskBinValue + 0.5;
	   my $newmask = "${TmpDir}mask.mnc";
	   Spawn(['mincmath', '-segment', '-const2', $minval, $maxval, $Mask, $newmask]);
	   $Mask = $newmask;
       }
       
       if ($MaskClassified) {
	   # Force classify to use the mask
	   AddDefaultArgs('classify', ['-mask', $Mask]);
       }
       
       if ($MaskTag) {
	   $StandardTags = 
	       &mask_tags($StandardTags, $Mask, 1, "${TmpDir}masked_standard.tag");
       }
   }
   elsif ($MaskTag || $MaskClassified) {
       die "Must use -mask with -mask_tag or -mask_classified\n";
   }

   # Verify merge file, if specified
   if (defined $MergeFile) {
       die "Merge file $MergeFile does not exist\n" if (! -e $MergeFile);
       if (defined $MergeLut) {
	   $MergeLut =~ s/,/ /g;
	   $MergeLut =~ s/::/ /g;
       }
   }
}

# ------------------------------ MNI Header ----------------------------------
#@NAME       : &SetHelp
#@INPUT      : none
#@OUTPUT     : none
#@RETURNS    : nothing
#@DESCRIPTION: Sets the $Help and $Usage globals, and registers them
#              with ParseArgs so that user gets useful error and help
#              messages.
#@METHOD     : 
#@GLOBALS    : $Help, $Usage
#@CALLS      : 
#@CREATED    : 95/08/25, Greg Ward (from code formerly in &ParseArgs)
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub SetHelp
{
    $Usage = <<USAGE;

Usage: $ProgramName [options] <in.mnc> [<in.mnc> ...] <classified.mnc>
       $ProgramName -help for details

USAGE

    $Help = <<HELP;

$ProgramName is used to classify stereotaxic MINC volumes.  It uses
the classify program (with mind-dist) and a set of standard sample
points to compute an initial volume classification.  This
classification is then optionally used to purge incorrect tag points
from the standard set, thus yielding a custom set of labels for the
particular subject.  The tag point set is then used by an ANN
classifier to classify the volume.

Input:
    <in.mnc>: input stereotaxic MRI volume(s)

Output: 
    <classified.mnc>: discrete classified MRI volume.

Steps involved:
   1. call classify -min_dist to compute initial classification
   2. correct standard tag file with min_dist classification
   3. call classify -ann to compute final classification
HELP

   Getopt::Tabular::SetHelp ($Help, $Usage);
}

sub print_version  {
  print "Program $ProgramName, built from:\n$LongVersion\n";
  exit;
}
