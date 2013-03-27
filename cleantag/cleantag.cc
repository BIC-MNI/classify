/*--------------------------------------------------------------------------
@COPYRIGHT  :
              Copyright 1997, Vasken Kollokian,
              McConnell Brain Imaging Centre,
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */
/* ----------------------------- MNI Header -----------------------------------
@NAME       : cleantag
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Nov 4, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
#include "config.h"

extern "C" {
#include <volume_io.h>
}
#include <ParseArgv.h>


/* function prototypes */

void parse_arguments(int argc, char* argv[]);
void load_volume(char *, VIO_Volume * );
void load_tag_file( char *);
void scan_and_clean_tags( char *clean_mode );
void write_tag_file(void);
int voxel_is_in_volume( VIO_Real vox1, VIO_Real vox2, VIO_Real vox3 );


/* global variables */

VIO_Status     status;
int        verbose = FALSE;
int        clobber = FALSE;
int        debug = 0;

VIO_Real       threshold = 0.9;      /* denote fuzzy acceptance level */
VIO_Real       diff_thresh = 0.1;    /* denote fuzzy proximity rejection level */

char       *tag_filename;                /* original tag filename */
char       *new_tag_filename;            /* cleaned tag filename */
char       *mask_filename = NULL;        /* mask volume */
int        mask_binvalue = 1;            /* value of mask foreground */
char       *comment;                     /* comment string */
char       *clean_mode = NULL;           /* string to denote cleaning mode*/

char       **fuzzy_label;
char       **input_filename;

VIO_Volume     *in_volume;
int        num_volumes;
VIO_Volume     mask_volume;
int        *volume_sizes;

VIO_Real       **tags;
char       **labels;
long       tagpoint;

VIO_Real       **new_tags;
char       **new_labels;

int        new_tagpoint = 0;
int        num_oldtags = 0;

int        fuzz_vol;                    /* index to  number of fuzzy volumes */
int        block_sizes[3] = { 1, 1, 1}; /* set the block size for vol. cache */

ArgvInfo argTable[] = {
  { NULL, ARGV_VERINFO, VERSION, NULL, NULL },
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Show progress"},
  
  {"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber,
     "Overwrite output file."},
  
  {"-debug", ARGV_INT, (char *) NULL, (char *) &debug,
     "Show debugging information"},
  
  {"-oldtag", ARGV_STRING, (char *) NULL, (char *) &tag_filename,
     "Specify the tag file to be cleaned"},      

  {"-newtag", ARGV_STRING, (char *) NULL, (char *) &new_tag_filename,
     "Specify the file name of the clean set of tag points"},      
   
  {"-mask", ARGV_STRING, (char *) NULL, (char *) &mask_filename,
     "Specify a mask to apply to the tag points"},      
   
  {"-maskbinvalue", ARGV_INT, (char *) NULL, (char *) &mask_binvalue,
     "Value of mask foreground"},

  {"-mode", ARGV_STRING, (char *) NULL, (char *) &clean_mode, 
     "Specify a quoted string 'xyz' to denote tag rejection mode.\n \
      Each of x,y,z being 0 or 1,  effects are additive.\n \
      x = 1, if tag class label <> voxel class (highest fuzzy voxel class),\n \
      y = 1, if fuzzy voxel class < threshold,\n \
      z = 1, if diff. between 2 highest fuzzy voxel classes < diff. threshold,\n"},

  {"-threshold", ARGV_FLOAT, (char *) 1, (char *) &threshold,
     "Set the fuzzy threshold rejection criterion (applies to 'y' under -mode)"},

  {"-difference", ARGV_FLOAT, (char *) NULL, (char *) &diff_thresh,
     "Set the difference threshold rejection criterion (applies to 'z' under -mode)"},

  {"-comment", ARGV_STRING, (char *) NULL, (char *) &comment, 
     "Specify a comment to be included in the cleaned tag file."},

  {NULL, ARGV_END, NULL, NULL, NULL}
};



int main(int argc, char *argv[])
{
  parse_arguments(argc, argv);

  /* set volume caching options to random */
  set_n_bytes_cache_threshold(0);
  set_default_max_bytes_in_cache(0);
  set_default_cache_block_sizes(block_sizes); 

  /* optimized for tag access */
  load_tag_file(tag_filename);

  ALLOC(volume_sizes, VIO_MAX_DIMENSIONS);

  if (mask_filename)
    load_volume( mask_filename, &mask_volume );

  if (num_volumes) {
    /* allocate memory for vol. pointer and start loading the volumes */
    ALLOC( in_volume, num_volumes );

    for_less ( fuzz_vol, 0, num_volumes ) 
      load_volume( input_filename[fuzz_vol], &in_volume[fuzz_vol] );

    get_volume_sizes(in_volume[0], volume_sizes);
  }
  else
    get_volume_sizes(mask_volume, volume_sizes);

  /* for each volume loaded, scan the the tag file and clean */
  scan_and_clean_tags( clean_mode );

  /* write the newly created tags to a file */
  write_tag_file();

  if (num_volumes)
    for_less (fuzz_vol, 0, num_volumes) 
      delete_volume(in_volume[fuzz_vol]);

  if (mask_filename)
    delete_volume(mask_volume);

  return (EXIT_SUCCESS);
} /* main */





/* ----------------------------- MNI Header -----------------------------------
@NAME       : parse_arguments
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 28, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void parse_arguments(int argc, char* argv[])
{

  int k;      /* volume filename index */
  int m = 2;  /* used to offset class argument */


  /* Call ParseArgv */
  if ( ParseArgv(&argc, argv, argTable, 0) ) {
    (void) fprintf(stderr, 
		   "\nUsage: %s [options] <fuzzy_class.mnc> <class_id> [<fuzzy_class.mnc> <class_id> ...]\n", 
		   argv[0]);
    (void) fprintf(stderr,   
		   "       %s [-help]\n\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if ( debug > 1) {

    fprintf(stdout, "the number of arguments left = %d\n", argc);

    for ( k = 0; k < argc; k++)
      fprintf(stdout, "argv[%d] = %s\n", k, argv[k]);
  }

  /* make sure you specify a new tag filename */
  if (!new_tag_filename) {

    printf("Please specify a new tag filename with -newtag switch\n");
    exit(EXIT_FAILURE);
  }
  
  /* make sure you specify a old tag filename */
  if (!file_exists(tag_filename)) {
    
    printf("Tag file %s does not exist.\n", tag_filename);
    exit(EXIT_FAILURE);
  }
  

  /* warn of clobbering if new tag file exists */
  if (new_tag_filename && !clobber && file_exists(new_tag_filename)) {
    
    printf("File %s exists.\n", new_tag_filename);
    printf("Use -clobber to overwrite.\n");
    exit(EXIT_FAILURE);
  }

  if ( !clean_mode && !mask_filename ) {
    fprintf( stdout, "Please specify -mode and/or -mask\n");
    exit(EXIT_FAILURE);
  }

  if (clean_mode) {
    /* make sure that clean_mode is only 3 characters long */
    if ( strlen( clean_mode ) != 3 )  {
      
      fprintf( stdout, "Cleaning mode should be 3 characters long. ex '011' \n");
      exit(EXIT_FAILURE);
    }
    
    for_less( k, 0, 3 ) {

      if ( debug >= 1 ) 
	fprintf( stdout, "clean_mode[%d] = %c\n", k, clean_mode[k]);
      
      if ( clean_mode[k] != '0' && clean_mode[k] != '1' ) {
	
	fprintf( stdout, "Invalid cleaning mode, please user 0 or 1.\n");
	exit(EXIT_FAILURE);
      }
    }
  }

  /* Check mask, if specified */
  if (mask_filename) {
    if (!file_exists(mask_filename)) {

      printf("Mask file %s does not exist.\n", mask_filename);
      exit(EXIT_FAILURE);
    }
  }
  
  /* make sure each volume has a class label */
  if ( ( argc - 1) % 2 != 0 ) {

    fprintf(stderr, "Please specify a class label.\n");   
    exit(EXIT_FAILURE); 
  }
     
  num_volumes = ( argc - 1 ) / 2;  /* count out progname */
 
  if ( debug > 1 ) 
    fprintf(stdout, "the number of volumes = %d\n", num_volumes);
 
  if ( threshold <= 0 ) {
    
    fprintf(stdout, "Please choose a positive threshold\n");
    exit(EXIT_FAILURE);
  }
  
  if (num_volumes) {
    ALLOC( input_filename, num_volumes );
    
    ALLOC( fuzzy_label, num_volumes );
    
    for_less ( k, 0, num_volumes) {
      
      if (!file_exists(argv[(k*m)+1])) {
	
	(void) fprintf(stderr, "filename `%s' not found. \n", argv[(k*m)+1]);
	exit(EXIT_FAILURE);
      }
      
      /* get the filename and class label */
      input_filename[k] = argv[(k*m)+1];
      
      fuzzy_label[k] = argv[(k*m)+2];
      
      if ( debug > 1 ) {
	
	fprintf( stdout, "input_filename[%d] = %s\n", k, input_filename[k]);
	fprintf( stdout, "fuzzy_label[%d] = %s\n", k, fuzzy_label[k]);
      }
      
    } /* for_less ( k...) */
  }
} /* parse_arguments() */



/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 28, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void load_volume(char *in_filename, VIO_Volume *volume)
{
  if (verbose)
    printf ("Processing volume %s...", in_filename);
  
  /* load the volume */
  status = input_volume( in_filename, 3, NULL, NC_UNSPECIFIED, 
			 FALSE, 0.0, 0.0,
			 TRUE, volume, (minc_input_options *) NULL ) ;
  
  if( status != VIO_OK )
    exit(EXIT_FAILURE);
  
  if (verbose)
    printf ("Done\n");

} /* load_volume(...) */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : scan_and_clean_tags
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Nov. 4, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void scan_and_clean_tags( char *clean_mode )
{
  VIO_Real wx, wy, wz;          /* world coordinates of tag */
  VIO_Real v1, v2, v3;          /* voxel ooordinates of tag */
  int  accept;              /* flag to accept or reject a tag point */

  int  hi_idx;                 /* highest index to fuzzy_class */
  int  hi2_idx;                /* second highest index to fuzzy_class */
  VIO_Real *fuzzy_class;           /* array to hold fuzzy values for classes */
  VIO_Real hi_fuzzy;               /* value of highest fuzzy class */
  VIO_Real hi2_fuzzy;              /* value of second highest fuzzy class */
  VIO_Real mask_value;             /* value of mask */
  int  lbl_idx;                /* index to fuzzy_class corresponding to tag label */

  int  *reject_count;           /* array to indicate number of rejected tags / class */
  int  *accept_count;           /* array to indicate number of accepted tags / class */

  VIO_Real mask_min = mask_binvalue - 0.5;
  VIO_Real mask_max = mask_binvalue + 0.5;

  /* reserve some memory for the fuzzy feature vector */
  if (num_volumes) {
    ALLOC( fuzzy_class, num_volumes );
    ALLOC( reject_count, num_volumes );
    ALLOC( accept_count, num_volumes );

    /* initialize reject_count, accept_count*/
    for_less( fuzz_vol, 0, num_volumes ) {
      
      reject_count[fuzz_vol] = 0;
      accept_count[fuzz_vol] = 0;
    }
  }

  if (debug > 1)
    fprintf(stdout, "Starting tag point loop, %d tag points\n", num_oldtags);

  /* repeat for each tag point */
  for_less( tagpoint, 0, num_oldtags) {

    /* initialize values */
    accept = TRUE;
    hi_idx = 0;
    hi2_idx = 0;
    hi_fuzzy = 0;
    hi2_fuzzy = 0;
    lbl_idx = 0;

    /* get the world coordinate of tag point */
    wx = tags[tagpoint][VIO_X];
    wy = tags[tagpoint][VIO_Y];
    wz = tags[tagpoint][VIO_Z];
   
    /* convert tag world to voxel */
    /* use first volume or mask   */
    convert_3D_world_to_voxel(num_volumes ? in_volume[0] : mask_volume, 
			      wx, wy, wz, &v1, &v2, &v3);

    /* if a mask is specified, check whether the tag point is in the mask */
    if (mask_filename) {
      mask_value = get_volume_real_value( mask_volume, VIO_ROUND(v1), VIO_ROUND(v2), VIO_ROUND(v3), 0, 0);
      
      if ((mask_value < mask_min) || (mask_value > mask_max)) {
	accept = FALSE;
	goto copy_tag;
      }
    }

    /* check to see if the tag point is in the volume before fetching the image
       value at that position - if not, then reject the point */
    if ( !voxel_is_in_volume( v1, v2, v3 ) ) {
	accept = FALSE;
	goto copy_tag;
    }

    if (!num_volumes)
      goto copy_tag;

    /* load fuzzy labels of the voxel into the fuzzy class vector */
    for_less( fuzz_vol, 0, num_volumes )
      fuzzy_class[fuzz_vol] = get_volume_real_value( in_volume[fuzz_vol], 
						    (int)v1, (int)v2, (int)v3, 0, 0);

    /* find the highest fuzzy index in the fuzzy_class[] array */
    for_less( fuzz_vol, 0, num_volumes ) {

      if ( hi_fuzzy < fuzzy_class[fuzz_vol] ) {

	 hi_fuzzy = fuzzy_class[fuzz_vol];
	 hi_idx = fuzz_vol;
       }

    }


    /* find the second highest fuzzy index in the fuzzy_class[] array */
    for_less( fuzz_vol, 0, num_volumes ) {

      if ( hi2_fuzzy < fuzzy_class[fuzz_vol] ) {

	if ( fuzz_vol != hi_idx ) {

	  hi2_fuzzy = fuzzy_class[fuzz_vol];
	  hi2_idx = fuzz_vol;
	}

      }

    } /* for_less */


    /* find the fuzzy_label index which corresponds to tag label */
    for_less( fuzz_vol, 0, num_volumes ) {

      if ( !strcmp (fuzzy_label[fuzz_vol], labels[tagpoint]) ) {
	
	lbl_idx = fuzz_vol;
	break;
      }

    } /* for_less */


    /* now start testing for the 'clean_mode string' */

    /* test for case of x = 1 */
    if ( clean_mode[0] == '1' ) {
      
      /* if tag label != highest class label, reject tag  */
      if ( strcmp( labels[tagpoint], fuzzy_label[hi_idx] ) ) {

	accept = FALSE;
	goto copy_tag;
      }

    }

    /* test for case of y = 1 */
    if ( clean_mode[1] == '1' ) {
      
      if ( fuzzy_class[lbl_idx] < threshold ) {

	accept = FALSE;
	goto copy_tag;
      }

    }

    /* test for case of z = 1 */
    if ( clean_mode[2] == '1' ) {
      
      if ( ( fuzzy_class[hi_idx] - fuzzy_class[hi2_idx] ) < diff_thresh ) {

	accept = FALSE;
	goto copy_tag;
      }

    }


  copy_tag:

    if ( accept ) {
	
      /* copy world coordinates */
      SET_ARRAY_SIZE( new_tags, new_tagpoint, new_tagpoint+1, 1000);
      ALLOC( new_tags[new_tagpoint], 3);
	
      new_tags[new_tagpoint][VIO_X] = wx;
      new_tags[new_tagpoint][VIO_Y] = wy;
      new_tags[new_tagpoint][VIO_Z] = wz;
      
      /* copy label */
      SET_ARRAY_SIZE( new_labels, new_tagpoint, new_tagpoint+1, 1000);
      ALLOC( new_labels[new_tagpoint], sizeof( labels[tagpoint] ) + 1);
      
      strcpy(new_labels[new_tagpoint], labels[tagpoint]);
      new_tagpoint++;
      
      if ( debug >= 4 ) {
	fprintf( stdout, "accepting tag %f %f %f -> %s\n",  
		wx, wy, wz, labels[tagpoint] );
      }

      if (num_volumes) 
	accept_count[lbl_idx]++;
      
    } /* if ( copy_tag )  */
    else {
      if ( debug >= 3 ) {
	
	fprintf( stdout, "rejecting tag %f %f %f -> %s\n",  
		wx, wy, wz, labels[tagpoint] );
	
      }

      if (num_volumes)
	reject_count[lbl_idx]++;
    }
    
  } /* for_less ( tagpoint...) */

  if ( verbose & num_volumes )
    for_less( fuzz_vol, 0, num_volumes ) 
      fprintf(stdout, "Rejected %d, accepted %d, class %s\n", 
	      reject_count[fuzz_vol], accept_count[fuzz_vol],
	      fuzzy_label[fuzz_vol]);
  
} /* scan_and_clean_tags(...) */

  
		
/* ----------------------------- MNI Header -----------------------------------
@NAME       : write_tag_file
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mar. 1, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void write_tag_file(void)
{


  if (verbose)
    printf("Writing tag file %s ...\n", new_tag_filename);
  
  if (output_tag_file(new_tag_filename, comment, 1, new_tagpoint, new_tags, 
		      NULL, NULL, NULL, NULL, new_labels) != VIO_OK)
    exit(EXIT_FAILURE);
  
}

    
/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_tag_file
@INPUT      : name of tag file
@OUTPUT     : 
@RETURNS    : number of tag points read
@DESCRIPTION: opens and loads a tag file
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 29, 1995 ( Vasco KOLLOKIAN)
@MODIFIED   : Oct 22, 1995 ( Vasco KOLLOKIAN)
---------------------------------------------------------------------------- */
void load_tag_file ( char *tag_filename )
{

  int num_vol, k;

  if (verbose)
    (void) fprintf(stdout, "Loading  tagfile %s\n", tag_filename);

  /* read the tag file */
  if ( input_tag_file(tag_filename, &num_vol, &num_oldtags,
		      &tags, NULL, NULL, NULL, NULL, &labels ) != VIO_OK ) {

    fprintf(stderr, "Error reading the tag file.\n");
    exit(EXIT_FAILURE);
  }

} /* load_tag_file (...) */



/* ----------------------------- MNI Header -----------------------------------
@NAME       : voxel_is_in_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: check to see if a voxel is in the volume (vol 0 is same as all)
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep 22, 1995 ( Vasco KOLLOKIAN)
@MODIFIED   : Mar 5, 1996 ( Vasco KOLLOKIAN)
---------------------------------------------------------------------------- */
int voxel_is_in_volume( VIO_Real vox1, VIO_Real vox2, VIO_Real vox3 )
{
 
  /* in_volume[0] is the volume against which the tags are verified */

  if ( vox1 < -0.5 || vox1 >= (VIO_Real) volume_sizes[0] - 0.5)
    return FALSE;
  
  else if ( vox2 < -0.5 || vox2 >= (VIO_Real) volume_sizes[1] - 0.5)
    return FALSE;
  
  else if ( vox3 < -0.5 || vox3 >= (VIO_Real) volume_sizes[2] - 0.5)
    return FALSE;

  else
    return TRUE;
}

