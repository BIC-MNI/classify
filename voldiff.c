/* ----------------------------- MNI Header -----------------------------------
@NAME       : voldiff
@INPUT      : argc, argv -command line arguments
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: compares bunch o volumes and produces a confusion matrix
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 8, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : $Log: voldiff.c,v $
@MODIFIED   : Revision 1.1  2011-01-06 19:18:58  jason
@MODIFIED   : added voldiff
@MODIFIED   :
@MODIFIED   : Revision 1.2  2001/07/08 21:03:39  crisco
@MODIFIED   : bug fix: under linux (egcs), the mask value test wasn't working if -O was used...; lowered the user_mask_value
@MODIFIED   :
@MODIFIED   : Revision 1.1  2001/07/08 20:58:19  crisco
@MODIFIED   : original (vasco\'s) version...
@MODIFIED   :
 * Revision 1.16  1996/03/24  03:37:59  vasco
 * Fixed the bug with confusion matrix alignment, to include 9 digit numbers.
 * reversed -ignore_outvoxel to -outvoxels_die, thus making ignore_outvoxel default.
 * Added set_default_max_bytes_in_cache(0), to enable cache of 1 block.
 *
 * Revision 1.15  1996/01/24  08:04:13  vasco
 * removed -exclude from the options list
 *
 * Revision 1.14  1996/01/24  07:48:57  vasco
 * Fixed a bug in calculating kappa_{numer,denom} by casting into (VIO_Real).
 *
 * Revision 1.13  1996/01/24  06:46:18  vasco
 * Removed sub_kappa and nozero_kappa, and all exclude related items.
 * Added modified version of per class kappa based on Bishop 1975.
 * Modified perviously displayed kappa to be renamed Ckappa for collapsed.
 * Added Total_Kappa, Brain_Kappa, Paren_Kappa as standard statistics.
 * Ultimately Brain_Kappa will be the most significant statistic to denote
 * similarity in healthy tissue classification.
 *
 * Revision 1.12  1995/12/04  18:41:16  vasco
 * added Total_kappa, and Xzero_kappa to report to help parsing.
 *
 * Revision 1.11  1995/09/24  22:22:23  vasco
 * added -mask option, to limit comparison to a given masked volume
 *
 * Revision 1.10  1995/09/24  18:18:59  vasco
 * I forgot to check it, thus the modifications I made
 *
 * Revision 1.9  1995/08/12  22:35:15  vasco
 * changed switchs, -matrix_size -> view_class_size
 * changed switch , -uppler_class -> c_matrix_size
 * separated compare_volume(i) to produce_confusion_matrix(i), and
 * calculate_similarity_measures()
 * fixed the alignment of vertical labels
 * removed -nozeroless, added -exclude <class>
 * changed total_voxels to c_matrix[class_total][class_total]
 * fixed excluded reporting scheme.
 * added kappa, ICC, AZms, and changed hit and miss into specificity
 * accuract, sensitivity, and error
 * started to work on weighted kappa
 *
 * Revision 1.8  1995/07/23  19:46:58  vasco
 * Also reports total hit and miss ratio without including class 0 in the
 * calculation.  And provided a -noexcluded switch to turn off
 * reporting of excluded calculation.
 *
 * Revision 1.7  1995/07/17  01:26:54  vasco
 * fixed a bug in calculating the percent hit ratio / class
 * total_voxels per class, should have been that of target and not source.
 * Also added %hit and %miss ratio of the whole volume.
 *
 * Revision 1.6  1995/06/15  06:33:14  vasco
 * give %hit and miss ratio per class, and combine stdout and file routines
 *
 * Revision 1.5  1995/06/07  05:20:25  vasco
 * a complete rewrite with parseargv,
 * loading of gold volume once, option of stdout and file
 * option of user specifying upper_class
 *
---------------------------------------------------------------------------- */

#include <ParseArgv.h>
#include <volume_io.h>

/* number of stats to report */
#define STATNUM 7

/* index of each stat */
#define CKAPPA        0   /* collapsed kappa */
#define AZSM          1
#define SENSITIVITY   2
#define ERROR         3
#define SPECIFICITY   4
#define ACCURACY      5
#define PCLASS_KAPPA  6   /* kappa per class */

#define DEFCMATSIZE   4  /* default confusion matrix size */
#define DEFVIEWSIZE   4  /* default confusion matrix view size */


/* Function declarations and prototyping */
void parse_arguments(int argc, char *argv[]);
void load_first_input_volume(void);
void load_input_volume(int volume_index);
void produce_confusion_matrix(int volume_index);
void produce_confusion_matrix_in_world_coor(int volume_index);
void set_confusion_matrix(void);
void calculate_similarity_measures(void);
void calculate_overall_measures(void);
void write_report(int volume_index);
void load_mask_volume(char *mask_file);
int voxel_is_in_volume( int index, VIO_Real vox1, VIO_Real vox2, VIO_Real vox3);
void  create_voxel_to_voxel_transform(VIO_Volume volume1, VIO_Volume volume2,
				      VIO_General_transform  *v1_to_v2 );

/* Global variables */  
VIO_Status     status;            /* status of loading and saving */

char       *pname;            /* the name of the invoked program */
int        verbose = FALSE;   /* inform user as things are processed */
int        clobber = FALSE;   /* overwrite existing report file */
int        exclude_class = 0; /* exclude class (VIO_X) in total hit and miss ratios */

int        debug = FALSE;     /* print debugging information */
int        stest = 0;         /* test stats on internal data */

VIO_Volume     *in_volume;              /* pointer to array of volumes */
VIO_Volume     mask_volume;             /* volume data to store the mask */

int        **in_volume_sizes;       /* 2D array to hold sizes */
int        num_volumes;            /* the number of volumes to compare */

char       **input_filename;        /* 1D char array to hold volume filenames*/
char       *output_filename;        /* confusion report filename */

int        view_class_size=DEFVIEWSIZE;   /* the size of the window to on CM */

VIO_Real       total_hit_ratio,
           excluded_total_hit_ratio, /* ratio without class zero */
           excluded_total_miss_ratio,
           total_miss_ratio,        /* total hit and miss ratio */
           total_kappa,             /* a measure of similarity ala Cohen */
           brain_kappa,             /* same measure for brain only */
           paren_kappa,             /* same measure for parenchyma (gry & wht) */
           sig_pclass_kappa,        /* same measure as above - Bishop */
           total_azsm;              /* a measure of similarity ala Zijdenbos */

int        c_matrix_size = DEFCMATSIZE;   /* the size of the confusion matrix */
int        class_total;        /* more discriptive index for marginals totals */
                               /* this is the same as c_matrix_size, init. later */
long       **c_matrix;              /* confusion matrix */
int        r,c;                     /* row and column counters */
VIO_Real       **c_table;               /* confusion table to report hit-miss/ class */ 

long       total_voxels_per_class, 
           total_hits,
           missed_voxels_per_class;

FILE       *matrep;                 /* FILE variable of output_filename */
char       *mask_filename;          /* filename of the mask volume */
VIO_Real       user_mask_value = 0.99;   /* default mask value to consider */

int        world_coor = FALSE;         /* flag for world coordinate traversal */
int        outvoxels_die = FALSE;      /* die if voxels fall outside volume 2 */

int        cache_set = FALSE;              /* set the caching option to false */
int        block_sizes[3] = {1,-1,-1};     /* default block size, (1 slice) */


ArgvInfo argTable[] = {

  {NULL, ARGV_HELP, NULL, NULL,
     " "},

  {"-debug", ARGV_CONSTANT, (char *) TRUE, (char *) &debug,
     "show debugging information."},

  {"-stest", ARGV_INT, (char *) NULL, (char *) &stest,
     "test stats on different internal data sets."},

  {"-clobber", ARGV_CONSTANT, (char *) TRUE, (char *) &clobber,
     "Overwrite confusion report file."},

  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Show progress"},

  {"-view_class_size", ARGV_INT, (char *) NULL, (char *) &view_class_size,
     "Specify the viewing size of the confusion matrix."},      

  {"-report", ARGV_STRING, (char *) NULL, (char *) &output_filename,
     "Write statistics to specified file"},      

  {"-c_matrix_size", ARGV_INT, (char *) NULL, (char *) &c_matrix_size,
     "Specify the actual size of the confusion matrix."},

  {"-mask", ARGV_STRING, (char *) NULL, (char *) &mask_filename, 
     "Specify volume which will be used as a mask"},

  {"-user_mask_value", ARGV_FLOAT, (char *) NULL, (char *) &user_mask_value,
     "Specify the mask value, in the mask volume to be used"},

  {"-world", ARGV_CONSTANT, (char *) TRUE, (char *) &world_coor,
     "Traverse the first volume in world coordinate space."},

  {"-outvoxels_die", ARGV_CONSTANT, (char *) TRUE, (char *) &outvoxels_die,
     "If voxels in first volume fall outside the second volume - abort."},

  {NULL, ARGV_HELP, NULL, NULL,
     "\nSpecify volume caching options.\n"  },

  {"-cached", ARGV_CONSTANT, (char *) TRUE, (char *) &cache_set, 
     "Specify that volume caching be used."},

  {"-block_sizes", ARGV_INT, (char *) 3, (char *) &block_sizes,
     "Set the size of the cache block for the volumes."},

  {NULL, ARGV_END, NULL, NULL, NULL}
};



/* Main program starts here */
int main(int argc, char *argv[])
{

  int vol_idx;

  parse_arguments(argc, argv);

  /* if you are testing measture of silumarity on dummy internal data */
  if ( stest  ) {

    set_confusion_matrix();

    calculate_similarity_measures();

    calculate_overall_measures();

    write_report(vol_idx);
	
    exit(EXIT_SUCCESS);

  }

  /* set specific caching options */
  if ( cache_set ) { 

    /* enable volume caching by setting threshold to 0 */ 
    set_n_bytes_cache_threshold(0);

    /* set max cache size to 0 to include only 1 block*/
    set_default_max_bytes_in_cache(0);

    /* randomize caching to access for 1 slice at a time - user selected */
    set_default_cache_block_sizes(block_sizes);

  }
  else {
    
    /* disable caching all together*/
    set_n_bytes_cache_threshold(-1);
  }

  load_first_input_volume();

  /* if a mask has been specified, load it in - invalid for world coordinates */
  if ( mask_filename && !world_coor) 
    load_mask_volume( mask_filename);

  for_less ( vol_idx, 1, num_volumes) {

    load_input_volume(vol_idx);

    /* if world coordinates have been specified, traverse in world */
    if ( world_coor) 
      produce_confusion_matrix_in_world_coor(vol_idx);
    else
      produce_confusion_matrix(vol_idx);

    calculate_similarity_measures();

    calculate_overall_measures();

    write_report(vol_idx);
	
    delete_volume(in_volume[vol_idx]);

  }

  exit(EXIT_SUCCESS);

} /* main */ 



/* ----------------------------- MNI Header -----------------------------------
@NAME       : parse_arguments
@INPUT      : argc, argv - command line arguments
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: parses command line arguments
@METHOD     : 
@GLOBALS    :
@CALLS      : 
@CREATED    : February 6, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void parse_arguments(int argc, char *argv[])
{

  int  i;

  pname = argv[0];

  /* Call ParseArgv - at least two volume has to be specified */
  if (ParseArgv(&argc, argv, argTable, 0) || ( argc < 2 ) ) { 

    if ( !stest  && ( argc < 3 ) ) {
    
      (void) fprintf(stderr, 
		 "\nUsage: %s <options> <infile1> <infile2> [infile2] ...\n", pname);
      (void) fprintf(stderr,   
		 "       %s [-help]\n\n", pname);
      exit(EXIT_FAILURE);
    }
  }
  
  if ( stest ) {

  }
  else {
  
    num_volumes = argc - 1;  /* count out progname */
  
    ALLOC( input_filename, num_volumes);    
  
    /* make sure that the typed filenames exist */
    for_less(i, 0, num_volumes ) {
      if (!file_exists(argv[i+1])) {
	(void) fprintf(stderr, "filename `%s' not found. \n", argv[i+1]);
	exit(EXIT_FAILURE);
      }
      /* get the filename of the volume */
      input_filename[i] = argv[i+1];
      if (debug)
	(void) fprintf(stderr, "input filename = `%s' \n", input_filename[i]);
    }
  
    /* if a report file is specified and exists, warn about clobbering */
    if (!clobber && output_filename && file_exists(output_filename)) {
      (void) fprintf(stderr,"File `%s' exists ! \n", output_filename);
      (void) fprintf(stderr,"Use -clobber to overwrite file.\n");
      exit(EXIT_FAILURE);
    }

    /* make sure view_class size is <= c_matrix_size */
    if ( view_class_size > c_matrix_size ) {
      fprintf(stderr, "view_class_size should be <= %d\n",  c_matrix_size);
      exit(EXIT_FAILURE);
    }

    /* make sure excluded class label is < c_matrix_size */
    if ( exclude_class >= c_matrix_size ) {
      fprintf(stderr, "excluded class label should be < %d\n", c_matrix_size);
      exit(EXIT_FAILURE);
    }


    /* open the report file */
    if ( output_filename) {
      matrep = fopen(output_filename, "w");
      if ( matrep == NULL) {
	fprintf(stderr, "Cannot open report file %s\n", output_filename);
	exit(EXIT_FAILURE);
      }
    }
  
    /* prevent use of mask while travering in world coordinate space,
       this is a temp measure not to deal with issues of world coor masks */
    if (mask_filename && world_coor ) {

      fprintf(stderr, "-mask and -world are mutually exclusive \n");
      exit(EXIT_FAILURE);
    }

    /* allocate area for filenames and derivatives */
    VIO_ALLOC2D( in_volume_sizes, num_volumes, VIO_MAX_DIMENSIONS); 
    ALLOC( in_volume, num_volumes);

  } /* if (test) */
   
  /* allocate area for confusion matrix */
  VIO_ALLOC2D( c_matrix, c_matrix_size+1, c_matrix_size+1); /* 1 extra for totals */
  
  /* confusion table has  STATNUM cols to report the following statistics, 
     kappa, azsm(Alex's), sensit., specif., accur. error, pclass_kappa */
  
  VIO_ALLOC2D( c_table, c_matrix_size, STATNUM); 
  

  /* to refer to class totals in the confusion matrix, it is more informative to
     refer to row-wise class total as c_matrix[<class>][class_total], than as
     c_matrix[<class>][c_matrix_size], and for columns as well,  functionally 
     they are the same */

  class_total = c_matrix_size ;



} /* parse_arguments */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_input_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : in_volume[j]  - the number of input volumes
              in_volume_sizes[MAX_VOLUMES][MAX_DIMESIONS] - sizes of volumes
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */
void load_input_volume(int volume_index)
{

  if (verbose) 
    fprintf (stderr, "Loading volume %s\n", input_filename[volume_index]);
    
  /* load the volume */
  status = input_volume(input_filename[volume_index], 
			3, 
			NULL,
			NC_BYTE, 
			FALSE, 
			0.0, 0.0,
			TRUE, 
			&in_volume[volume_index],
			(minc_input_options *) NULL ) ;
    
  if ( status != VIO_OK )
    exit(EXIT_FAILURE);
    

  /* get the volume sizes */
  get_volume_sizes(in_volume[volume_index], in_volume_sizes[volume_index]);
  
  /* check to see if volumes are of same size in each dim, if not in world*/
  if ( !world_coor ) {

    if (in_volume_sizes[volume_index][VIO_X] != in_volume_sizes[0][VIO_X]) {
      (void) fprintf(stderr,"Error - Volume size mismatch in X dimension\n");
      exit(EXIT_FAILURE);
    }
    
    if (in_volume_sizes[volume_index][VIO_Y] != in_volume_sizes[0][VIO_Y]) {
      (void) fprintf(stderr,"Error - Volume size mismatch in Y dimension\n");
      exit(EXIT_FAILURE);
    }
    
    if (in_volume_sizes[volume_index][VIO_Z] != in_volume_sizes[0][VIO_Z]) {
      (void) fprintf(stderr,"Error - Volume size mismatch in Z dimension\n");
      exit(EXIT_FAILURE);
    }

  } /*  if ( !world_coor ) */


} /* load_input_volume */
        

/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_first_input_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : load the first volume - the gold standard
@CALLS      : 
@CREATED    : June 6, 1995
@MODIFIED   : 
---------------------------------------------------------------------------- */
void load_first_input_volume(void)
{

    
  if (verbose) 
    fprintf (stderr, "Loading first volume %s\n", input_filename[0]);

  if ( debug && cache_set ) 
    fprintf (stderr, "Caching has been specified\n");

  /* load the volume */
  status = input_volume(input_filename[0], 
			3, 
			NULL,
			NC_BYTE, 
			FALSE, 
			0.0, 0.0,
			TRUE, 
			&in_volume[0],
			(minc_input_options *) NULL ) ;
    
  if ( status != VIO_OK )
    exit(EXIT_FAILURE);

  /* get the volume sizes */
  get_volume_sizes(in_volume[0], in_volume_sizes[0]);

} /* load_first_input_volume */
        


/* ----------------------------- MNI Header -----------------------------------
@NAME       : produce_confusion_matrix(int volume_index)
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: produces a confusion matrix for  a give volume (volume_index) 
              with the first volume.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Aug. 6, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void produce_confusion_matrix(int volume_index)
{
   
  int        v1, v2, v3;                    /* voxel indeces */
  VIO_Real       image_value1, image_value2;    /* values to be compared */
  VIO_Real       mask_value;                    /* value of the mask voxel */

  /* in the next two loops, and then in subsequent ones in this programm, 
     for_less (r, 0, c_matrix_size ) === for ( r = 0; r < c_matrix_size; r++)
     ie. if c_matrix_size = 8, r ranges in [0..7]
     don't get confused that matrix totals are being included also */

  /* initialize  the confusion matrix,  including marginal totals (+1) */ 
  for_less(r, 0, c_matrix_size+1)
    for_less(c, 0, c_matrix_size+1)
      c_matrix[r][c] = 0;

  /* initialize  the confusion table */
  for_less(r, 0, c_matrix_size) 
    for_less(c, 0, STATNUM) 
      c_table[r][c] = 0.0;  

  /* start traversing each of volume1, volume2 */
  if (verbose)
    fprintf(stderr, "Traversing volumes %s %s \n", 
	    input_filename[0], input_filename[volume_index]);

  for ( v1 = 0; v1 < in_volume_sizes[volume_index][0]; v1++) {
    for ( v2 = 0; v2 < in_volume_sizes[volume_index][1]; v2++) {
      for ( v3 = 0; v3 < in_volume_sizes[volume_index][2]; v3++) {

	
	/* if a mask is chosen, ( by having a non-NULL mask_filename )
	   and mask value is > user_mask_value, consider only that
	   voxel, in the confusion matrix. */

	if ( mask_filename ) {
	  
	  /* make sure the map is previously loaded, and has the same size */
	  GET_VALUE_3D( mask_value, mask_volume, v1, v2, v3);
	  
	  if (  mask_value  < user_mask_value ) {
	    
	    if( 0 && debug) 
	      fprintf( stdout, "mask_value= %f < %f\n", mask_value, user_mask_value);

	    goto NEXT_VOXEL;
	  }

	  /* this might be a future consideration spot for 
	     applying mask to only volume1, or volume2, separatly
	     
	     if (  mask_value  < user_mask_value ) {

	       switch ( option ) {
	     
	         case 1:
	     
	           goto NEXT_VOXEL;
	     
	         case 2:
	     
		   do something
	           break;
	     
	         case 3:
		 
		   do something
		   break;
	       }
	     }
	  */
	} /* if (mask_filename) */
	     
	GET_VALUE_3D(image_value1, in_volume[0], v1, v2, v3);
	GET_VALUE_3D(image_value2, in_volume[volume_index], v1, v2, v3);
	      
	if (image_value1 >= c_matrix_size) {
	  fprintf(stderr, "%s may not be a classified volume\n", input_filename[0]); 
	  fprintf(stderr, "Label value = %d, increase c_matrix_size\n",
		           (int)VIO_ROUND(image_value1) );
	  exit(EXIT_FAILURE);
	}

	if (image_value2 >= c_matrix_size) {
	  fprintf(stderr, "%s may not be a classified volume\n", 
       		           input_filename[volume_index]); 
	  fprintf(stderr, "Label value = %d, increase c_matrix_size\n",
		           (int)VIO_ROUND(image_value2) );

	  exit(EXIT_FAILURE);
	}
		
	c_matrix[VIO_ROUND(image_value1)][VIO_ROUND(image_value2)]++;
	
        NEXT_VOXEL: ;
      }
    }
  }
    
  if (verbose)
    fprintf(stderr, "Producing confusion Matrix...\n\n");
  
  /* calclate the total number of voxels for each class row-wise */
  for_less(r, 0, c_matrix_size)
    for_less(c, 0, c_matrix_size) 
      c_matrix[r][class_total] += c_matrix[r][c];

  /* calclate the total number of voxels for each class column-wise */
  for_less(c, 0, c_matrix_size)
    for_less(r, 0, c_matrix_size) 
      c_matrix[class_total][c] += c_matrix[r][c];


} 



/* ----------------------------- MNI Header -----------------------------------
@NAME       : produce_confusion_matrix_in_world(int volume_index)
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: produces a confusion matrix for  a give volume (volume_index) 
              with the first volume, while the first volume is being traversed
	      in world coordinate space. This is good to compare a high-res 
	      (0.5mm) phantom, with lower res phantoms, to assess the effect 
	      of partial volme.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mar 16, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void produce_confusion_matrix_in_world_coor(int volume_index)
{
   
  int        v1, v2, v3;                    /* voxel indeces */
  VIO_Real       image_value1, image_value2;    /* values to be compared */
  VIO_Real       wx, wy, wz;                    /* world coordinates  */
  VIO_Real       invol_v1, invol_v2, invol_v3;  /* voxel coordinates of in_volume[idx] */
  VIO_General_transform trans;                  /* trans from vol1 voxel to vol2 voxel */

  /* in the next two loops, and then in subsequent ones in this programm, 
     for_less (r, 0, c_matrix_size ) === for ( r = 0; r < c_matrix_size; r++)
     ie. if c_matrix_size = 8, r ranges in [0..7]
     don't get confused that matrix totals are being included also */

  /* initialize  the confusion matrix,  including marginal totals (+1) */ 
  for_less(r, 0, c_matrix_size+1)
    for_less(c, 0, c_matrix_size+1)
      c_matrix[r][c] = 0;

  /* initialize  the confusion table */
  for_less(r, 0, c_matrix_size) 
    for_less(c, 0, STATNUM) 
      c_table[r][c] = 0.0;  

  /* start traversing each of volume1, volume2 */
  if (verbose)
    fprintf(stderr, "Traversing volumes %s (in voxel) %s (in world)\n", 
	            input_filename[0], input_filename[volume_index]);

  /* compose the voxel to voxel transform */
  create_voxel_to_voxel_transform(in_volume[0], in_volume[volume_index], &trans);

  for ( v1 = 0; v1 < in_volume_sizes[0][0]; v1++) {
    for ( v2 = 0; v2 < in_volume_sizes[0][1]; v2++) {
      for ( v3 = 0; v3 < in_volume_sizes[0][2]; v3++) {

	
	GET_VALUE_3D(image_value1, in_volume[0], v1, v2, v3);


	/* voxel coor of in_volume[0] is checked against the world,
	   coordinate of the in_volume[voxel_idx] at the very same  location */

        /* the following method is a slow way of doing things, a better faster
	    method was suggested by David MacDonald, that concatanates transforms.
	    This is left here just for comparative purposed

	convert_3D_voxel_to_world(in_volume[0], 
				  (VIO_Real)v1,
				  (VIO_Real)v2,
				  (VIO_Real)v3,
				  &wx,
				  &wy,
				  &wz);

	convert_3D_world_to_voxel(in_volume[volume_index], 
				  wx,
				  wy,
				  wz,
				  &invol_v1,
				  &invol_v2,
				  &invol_v3);
      
        */

	general_transform_point( &trans, (VIO_Real) v1, (VIO_Real) v2, (VIO_Real) v3,
				 &invol_v1, &invol_v2, &invol_v3 );

	if ( voxel_is_in_volume( volume_index, invol_v1, invol_v2, invol_v3) ) {

	  image_value2 = get_volume_real_value(in_volume[volume_index], 
					       VIO_ROUND(invol_v1), 
					       VIO_ROUND(invol_v2), 
					       VIO_ROUND(invol_v3), 
					       0, 0);
          
	}
	else if ( outvoxels_die ) {

	  fprintf( stderr, "Voxels fall outside volume. Aborted... \n");
	  exit(EXIT_FAILURE);

	}

	else {

	  goto NEXT_VOXEL;

	}
	
	if (image_value1 >= c_matrix_size) {
	  fprintf(stderr, "%s may not be a classified volume\n", input_filename[0]); 
	  fprintf(stderr, "Label value = %d, increase c_matrix_size\n",
		           (int)VIO_ROUND(image_value1) );
	  exit(EXIT_FAILURE);
	}

	if (image_value2 >= c_matrix_size) {
	  fprintf(stderr, "%s may not be a classified volume\n", 
       		           input_filename[volume_index]); 
	  fprintf(stderr, "Label value = %d, increase c_matrix_size\n",
		           (int)VIO_ROUND(image_value2) );

	  exit(EXIT_FAILURE);
	}
		
	c_matrix[VIO_ROUND(image_value1)][VIO_ROUND(image_value2)]++;
	
        NEXT_VOXEL: ;
      }
    }
  }
    
  if (verbose)
    fprintf(stderr, "Producing confusion Matrix...\n\n");
  
  /* calclate the total number of voxels for each class row-wise */
  for_less(r, 0, c_matrix_size)
    for_less(c, 0, c_matrix_size) 
      c_matrix[r][class_total] += c_matrix[r][c];

  /* calclate the total number of voxels for each class column-wise */
  for_less(c, 0, c_matrix_size)
    for_less(r, 0, c_matrix_size) 
      c_matrix[class_total][c] += c_matrix[r][c];


} 




/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_confusion_matrix(void)
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: produces a dummy confusion matrix to test stats
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Aug. 15, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void set_confusion_matrix(void)
{
   

  /* in the next two loops, and then in subsequent ones in this programm, 
     for_less (r, 0, c_matrix_size ) === for ( r = 0; r < c_matrix_size; r++)
     ie. if c_matrix_size = 8, r ranges in [0..7]
     don't get confused that matrix totals are being included also */

  /* initialize  the confusion matrix,  including marginal totals (+1) */ 
  for_less(r, 0, c_matrix_size+1)
    for_less(c, 0, c_matrix_size+1)
      c_matrix[r][c] = 0;

  /* initialize  the confusion table */
  for_less(r, 0, c_matrix_size) 
    for_less(c, 0, STATNUM) 
      c_table[r][c] = 0.0;  

  /* start setting custom c_matrix values here */

  switch (stest) {

  case 1 : 

    /* set 1 - perfect agreement */
    c_matrix[0][0] = 6000;
    c_matrix[1][1] =  100;
    c_matrix[2][2] =  200;
    c_matrix[3][3] =  300; 
    break;

  case 2:

    /* set 2 */
    c_matrix[0][0] = 5970;
    c_matrix[0][1] =   30;
    c_matrix[1][1] =  100;
    c_matrix[2][2] =  200;
    c_matrix[3][3] =  300; 
    break;

  case 3:

    /* set 3 */
    c_matrix[0][0] = 6000;
    c_matrix[1][1] =  100;
    c_matrix[2][1] =   30;
    c_matrix[2][2] =  170;
    c_matrix[3][3] =  300; 
    break;

  case 4:

    /* set 4 */
    c_matrix[0][0] = 6000;
    c_matrix[1][0] =   30;
    c_matrix[1][1] =   70;
    c_matrix[2][2] =  200;
    c_matrix[3][3] =  300; 
    break;

  case 5:

    /* set 5 */
    c_matrix[0][0] = 6000;
    c_matrix[1][1] =   70;
    c_matrix[1][2] =   30;
    c_matrix[2][2] =  200;
    c_matrix[3][3] =  300; 
    break;

  case 6:

    /* set 6 */
    c_matrix[0][0] = 6000;
    c_matrix[1][1] =   70;
    c_matrix[1][3] =   30;
    c_matrix[2][2] =  200;
    c_matrix[3][3] =  300; 
    break;

  case 7:

    /* set 7 */
    c_matrix[0][0] = 6000;
    c_matrix[1][1] =   70;
    c_matrix[1][2] =   15;
    c_matrix[1][3] =   15;
    c_matrix[2][2] =  200;
    c_matrix[3][3] =  300; 
    break;

  case 8:

    /* set 8 */
    c_matrix[0][0] = 6000;
    c_matrix[1][1] =   70;
    c_matrix[1][2] =   15;
    c_matrix[2][1] =   15;
    c_matrix[2][2] =  200;
    c_matrix[3][3] =  300; 
    break;

  case 9:

    /* set 9 */
    c_matrix[0][0] = 6000;
    c_matrix[1][0] =   15;
    c_matrix[1][1] =   70;
    c_matrix[2][1] =   15;
    c_matrix[2][2] =  200;
    c_matrix[3][3] =  300; 
    break;

  case 10:

    /* set 10 */
    c_matrix[0][0] = 6000;
    c_matrix[0][1] =   15;
    c_matrix[1][0] =   15;
    c_matrix[1][1] =   70;
    c_matrix[2][2] =  200;
    c_matrix[3][3] =  300; 
    break;

  case 11:

    /* set 11 */
    c_matrix[0][0] = 6000;
    c_matrix[1][2] =   30;
    c_matrix[1][1] =   70;
    c_matrix[2][2] =  170;
    c_matrix[2][1] =   30;
    c_matrix[3][1] =   30; 
    c_matrix[3][3] =  270; 
    break;

  case 12:

    /* set 12 */
    c_matrix[0][0] = 6000;
    c_matrix[1][1] =   70;
    c_matrix[1][2] =   15;
    c_matrix[1][3] =   15;
    c_matrix[2][1] =   15;
    c_matrix[2][2] =  170;
    c_matrix[2][3] =   15;
    c_matrix[3][1] =   15;
    c_matrix[3][2] =   15;
    c_matrix[3][3] =  270; 
    break;

  case 13:

    /* set 13 */
    c_matrix[0][0] = 6000;
    c_matrix[0][1] =    0;
    c_matrix[0][2] =    0;
    c_matrix[0][3] =    0;
    c_matrix[1][0] =   30;
    c_matrix[1][1] =   70;
    c_matrix[2][0] =   30;
    c_matrix[2][2] =  170;
    c_matrix[3][0] =   30;
    c_matrix[3][3] =  270; 
    break;

  case 14:

    /* set 14 */
    c_matrix[0][0] = 5870;
    c_matrix[0][1] =   55;
    c_matrix[0][2] =   35;
    c_matrix[0][3] =   40;
    c_matrix[1][0] =   30;
    c_matrix[1][1] =   70;
    c_matrix[2][0] =   30;
    c_matrix[2][2] =  170;
    c_matrix[3][0] =   30;
    c_matrix[3][3] =  270; 
    break;

  case 15:

    /* set 15 */
    c_matrix[0][0] = 6000;
    c_matrix[0][1] =    0;
    c_matrix[0][2] =    0;
    c_matrix[0][3] =    0;
    c_matrix[1][0] =    0;
    c_matrix[1][1] =  100;
    c_matrix[2][0] =    0;
    c_matrix[2][2] =  200;
    c_matrix[3][0] =  300;
    c_matrix[3][3] =    0; 
    break;

  case 16:

    /* set 16 */
    c_matrix[0][0] =    0;
    c_matrix[0][1] = 2000;
    c_matrix[0][2] = 2000;
    c_matrix[0][3] = 2000;
    c_matrix[1][0] =    0;
    c_matrix[1][1] =  100;
    c_matrix[2][0] =    0;
    c_matrix[2][2] =  200;
    c_matrix[3][0] =    0;
    c_matrix[3][3] =  300; 
    break;

  case 17:

    /* set 17 */
    c_matrix[0][0] = 30000000;
    c_matrix[0][1] = 10;
    c_matrix[0][2] = 100;
    c_matrix[0][3] = 1000;
    c_matrix[1][1] =  10000;
    c_matrix[2][2] =  200000;
    c_matrix[3][3] =  30000000; 
    break;


  case 18:

    /* set 18 */
    c_matrix[0][1] = 2000;
    c_matrix[0][2] = 2000;
    c_matrix[0][3] = 2000;
    c_matrix[1][0] =  100;
    c_matrix[2][0] =  200;
    c_matrix[3][0] =  300; 
    break;



  default:
    fprintf( stderr, "non_existant set, sets available 1-18\n");
    exit(EXIT_FAILURE);

  }
 
  /* calclate the total number of voxels for each class row-wise */
  for_less(r, 0, c_matrix_size)
    for_less(c, 0, c_matrix_size) 
      c_matrix[r][class_total] += c_matrix[r][c];

  /* calclate the total number of voxels for each class column-wise */
  for_less(c, 0, c_matrix_size)
    for_less(r, 0, c_matrix_size) 
      c_matrix[class_total][c] += c_matrix[r][c];


} 



/* ----------------------------- MNI Header -----------------------------------
@NAME       : calcutate_similarity_measures(void)
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: compares a give volume (volume_index) with the first volume.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Aug 6, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void calculate_similarity_measures(void)
{
   
  /* much of the statistics in this routine is based on the following grid:

                               rater 2
	    	    class 0    class 1      total

	    class 0	a         b          a+b
   rater 1		   
            class 1     c         d          c+d

	    total      a+c       b+d         a+b+c+d
            
          +-
   ie.  + ab	    
        - cd
  
   and the entire confusion matrix is reduced to the above 2 elem matrix.
   where;


   a = c_matrix[class_idx][class_idx]               :  true positive
   b = c_matrix[class_idx][class_total]   - a       :  false negative
   c = c_matrix[class_total][class_idx]   - a       :  false positive
   d = c_matrix[class_total][class_total] - (a+b+c) :  true negative

   to calculate various measures of similarity the following formuale are used:

   kappa = (Po - Pc) / (1 - Pc) 

   Po = ( a + d ) / ( a + b + c + d )
   Pc = ( ( ( a + b ) * ( a + c ) / ( a + b + c + d ) ) +
          ( ( b + d ) * ( c + d ) / ( a + b + c + d ) ) ) / ( a + b + c + d )

   this is only for information - icc has been removed	  
   icc = ( ( 4*a*d ) - ( b*b + 2*b*c + c*c ) ) /
         ( ( 2*a + b + c ) * ( b + c + 2*d ) )

   azsm = ( 2 * a ) / ( (2 * a) + b + c ) (This is Alex's thesis measure)

   sensitivity = a / ( a + b )
   specificity = d / ( c + d )
   accuracy    = ( a + d ) / ( a + b + c + d)  same as Po        


   Also to calculate per class kappa without collapsing the matrix, the following
   formula was obtained from Bishop 1975 (p 397) :

   Ki = ( NXii - Xi+VIO_X+i ) / ( NXi+ - Xi+VIO_X+i )     kappa for each class i
   
   where N   = total number of voxels
         Xii = hits per class
         Xi+ = column total (row-wise)
         VIO_X+i = row total (column-wise)

   and   total kappa = Sig ( numerator ) / Sig (denominator )  */

  int class_idx;

  long hit_voxels_per_class;

  VIO_Real a, b, c, d, Po, Pc, ckappa, azsm;
  VIO_Real sensitivity, specificity, accuracy;

  VIO_Real pclass_kappa;        /* per class kappa ala Bishop 1975 */
  VIO_Real kappa_numer = 0.0;   /* used for calculating pclass_kappa */
  VIO_Real kappa_denom = 0.0;

  /* total kappa fraction components */
  VIO_Real sig_kappa_numer = 0.0; /* verification for total kappa */
  VIO_Real sig_kappa_denom = 0.0;

  /* brain kappa fraction components */
  VIO_Real sig_brain_kappa_numer = 0.0; 
  VIO_Real sig_brain_kappa_denom = 0.0;

  /* paren kappa fraction components */
  VIO_Real sig_paren_kappa_numer = 0.0; 
  VIO_Real sig_paren_kappa_denom = 0.0;

  total_hits   = 0;

  /* calculate the total number of elements in the matrix, by adding
     up class_total for all the classes.  store this value in the last
     diagonal element which is identical to (previously) total_voxels */
  for_less ( class_idx, 0, c_matrix_size ) {

    c_matrix[class_total][class_total] += c_matrix[class_idx][class_total];
  }


  /* calculate various similarity measures for each class */
  for_less ( class_idx, 0, c_matrix_size ) {


    /* get the definitions of a, b, c, and d */
    a = (VIO_Real) c_matrix[class_idx][class_idx];        
    b = (VIO_Real) c_matrix[class_idx][class_total]   - a;
    c = (VIO_Real) c_matrix[class_total][class_idx]   - a;
    d = (VIO_Real) c_matrix[class_total][class_total] - ( a + b + c );

    if (debug) {

      fprintf( stdout, "idx = %d\n", class_idx);
      fprintf( stdout, "a   = %f\n", a);
      fprintf( stdout, "b   = %f\n", b);
      fprintf( stdout, "c   = %f\n", c);
      fprintf( stdout, "d   = %f\n", d);
      fprintf( stdout, "tot = %f\n\n", (a+b+c+d));


    }


    /* calculate the hit and miss */
    hit_voxels_per_class = c_matrix[class_idx][class_idx];

    /* calculate total hits for an overall measure of the classifier hit ratio */
    total_hits += (int) a ;

    /* make sure the total voxels per class is that of master (target volume),
     rather than the volume being tested - the future of this statistics
     remains uncertain, however, I keep it here for comparative reasons,
     total_voxels_per_class = c_matrix[class_idx][class_total];
     missed_voxels_per_class = total_voxels_per_class - hit_voxels_per_class;


    /* if total voxels per class is not zero, calculate statistics
       and store results in the apropriate position in the c_table */
    
    if ( c_matrix[class_idx][class_total] != 0 ) {

      /* starts calculating kappa per class in this section */

      /* agreement due to observation */
      Po = ( a + d ) / ( a + b + c + d );

      /* agreement due to chance */
      Pc = (  ( (a+b) * (a+c) / (a+b+c+d) ) +
              ( (b+d) * (c+d) / (a+b+c+d) )   ) / ( a + b + c + d );

      /* calculate ckappa */
      ckappa = (Po - Pc) / (1 - Pc) ;

      /* store in c_table */
      c_table[class_idx][CKAPPA] = ckappa;

      /* starts calculating azsm per class in this section */
      azsm = ( 2*a ) / ( (2*a) + b + c );

      /* store in c_table */
      c_table[class_idx][AZSM] = azsm;

      /* calculate sens. spec. acc. and error in this section */
      sensitivity = a / ( a + b );
      specificity = d / ( c + d );
      accuracy    = Po;  /* this is same as Po */

      /* calculate per class kappa, Ki = ( NXii - Xi+VIO_X+i ) / ( NXi+ - Xi+VIO_X+i ) */


      kappa_numer = (  
		       ( (VIO_Real) c_matrix[class_total][class_total] * (VIO_Real) a ) -
		       ( (VIO_Real) c_matrix[class_idx][class_total] *  
			 (VIO_Real) c_matrix[class_total][class_idx] ) 
		     ); 


      kappa_denom =  (  
		       ( (VIO_Real) c_matrix[class_total][class_total] * 
		         (VIO_Real) c_matrix[class_idx][class_total] ) -
		       ( (VIO_Real) c_matrix[class_idx][class_total]  *
		         (VIO_Real) c_matrix[class_total][class_idx] ) 
		     );


      /* keep track of total kappa fraction components */
      if ( debug ) {

	fprintf( stdout, "kappa_numer = %f\n", kappa_numer);
	fprintf( stdout, "kappa_denom = %f\n\n", kappa_denom);
      }

      sig_kappa_numer += kappa_numer;
      sig_kappa_denom += kappa_denom;

      /* csf should be added to brain_kappa */
      if ( class_idx == 1 ) {

	sig_brain_kappa_numer += kappa_numer;
	sig_brain_kappa_denom += kappa_denom;
      }

      /* grey and white matter should be added to brain_kappa and paren_kappa */
      if ( class_idx == 2 || class_idx == 3 ) {

	sig_brain_kappa_numer += kappa_numer;
	sig_brain_kappa_denom += kappa_denom;
	sig_paren_kappa_numer += kappa_numer;
	sig_paren_kappa_denom += kappa_denom;
      }
	 

      if ( debug ) {

	fprintf( stdout, "sig_brain_kappa_numer = %f\n", sig_brain_kappa_numer);
	fprintf( stdout, "sig_brain_kappa_denom = %f\n", sig_brain_kappa_denom);
	fprintf( stdout, "sig_paren_kappa_numer = %f\n", sig_paren_kappa_numer);
	fprintf( stdout, "sig_paren_kappa_denom = %f\n\n", sig_paren_kappa_denom);

      }

      pclass_kappa = kappa_numer / kappa_denom;

      /* store results in c_table */
      c_table[class_idx][SENSITIVITY] = sensitivity;
      c_table[class_idx][ERROR] = 1.0 - (VIO_Real) c_table[class_idx][SENSITIVITY];
      c_table[class_idx][SPECIFICITY] = specificity;
      c_table[class_idx][ACCURACY] = accuracy;
      c_table[class_idx][PCLASS_KAPPA] = pclass_kappa;

     } /* if ( total_voxels_per_class != 0 ) */

  } /* for_less */    
    
  /* give an overall performance measure */
  total_hit_ratio  =  (VIO_Real) total_hits / 
                      (VIO_Real) c_matrix[class_total][class_total];
                      
  total_miss_ratio =  1.0 - total_hit_ratio;

  /* verification that sigma numer and denom result in total_kappa */
  sig_pclass_kappa = sig_kappa_numer / sig_kappa_denom;

  /* calculate brain_kappa and paren_kappa */
  brain_kappa = sig_brain_kappa_numer / sig_brain_kappa_denom;
  paren_kappa = sig_paren_kappa_numer / sig_paren_kappa_denom;

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : void calculate_overall_measures(void)
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION:  calculate total_kappa -a measure of similarity ala Jacob Cohen (1960) 
               total_kappa = ( Po - Pc ) / ( 1 - Pc )    where
	       Po = total hit ratio,   Sig ( diagonals ) / N 
	       Pc = agreement due to chance , Sig ( col_total * row_total / N) ) / N
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    :  Aug 6, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void  calculate_overall_measures(void)
{

  int class_idx;  

  VIO_Real Po = 0.0;
  VIO_Real Pc = 0.0;
  VIO_Real azsm_numerator = 0.0;           /* terms to denote num and denom of azsm */
  VIO_Real azsm_dinominator = 0.0;
  VIO_Real A, B, C, D;                   /* to calculate 2a+b+c, different than alex */
  int  i,j;                            /* counters */


  /* start calculating overall kappa here */

  /* this following is the same as total_hit_ratio */
  Po = (VIO_Real) total_hits / (VIO_Real) c_matrix[class_total][class_total] ;

  for_less ( class_idx, 0, c_matrix_size ) {

    Pc += ( (VIO_Real) c_matrix[class_idx][class_total]   *  /* col_total */
	    (VIO_Real) c_matrix[class_total][class_idx]   /  /* row_total */
	    (VIO_Real) c_matrix[class_total][class_total] ); /* N */
  }

  Pc = Pc /  (VIO_Real) c_matrix[class_total][class_total] ;  /* final division by N */

  if ( debug ) {

    fprintf ( stdout, "Pc = %f\n", Pc);

  }

  total_kappa = ( Po - Pc ) / ( 1.0 - Pc );


  /* start calculating overall AZSM here - this is derived by Alex Zijdenbos

     it is based on  azsm = 2a / ( 2a + b + c ), but for a matrix of more than
     two classes: for ex 5 class matrix will be as follows

     0 0 0 0 0       0 0 0 0 0
     0 1 1 1 1       0 1 0 0 0  
     0 1 1 1 1   +   0 0 1 0 0  
     0 1 1 1 1       0 0 0 1 0 
     0 1 1 1 1       0 0 0 0 1
   
    -----------divided by -------------------
      0 0 0 0 0       0 0 0 0 0      0 1 1 1 1 
      0 1 1 1 1       1 0 0 0 0      0 0 0 0 0  
   2x 0 1 1 1 1   +   1 0 0 0 0  +   0 0 0 0 0  
      0 1 1 1 1       1 0 0 0 0      0 0 0 0 0 
      0 1 1 1 1       1 0 0 0 0      0 0 0 0 0


      which is  ( Sig ij  + Sig ii ) / ( Sig ij + Sig i0 + Sig 0j) */


  /* first term of numerator, index starts at 1, to exclude class 0 */
  for_less ( i, 1, c_matrix_size )
    for_less ( j, 1, c_matrix_size )
      azsm_numerator += c_matrix[i][j];

  /* since first term of dinominator is twice that of numerator - calculate  */
  azsm_dinominator = 2 * azsm_numerator;

  /* continue adding second term */
  for_less ( i, 1, c_matrix_size )
    azsm_numerator += c_matrix[i][i];

  /* now start calculating dinominator */
  /* continue adding second term */
  for_less ( i, 1, c_matrix_size )
    azsm_dinominator += c_matrix[i][0]; /* 0 is index of class 0 */
 
  /* continue adding third term */
  for_less ( j, 1, c_matrix_size )
    azsm_dinominator += c_matrix[0][j];  /* 0 is index of class 0 */

  total_azsm = azsm_numerator / azsm_dinominator;

  D = c_matrix[0][0];
  C = c_matrix[0][class_total] - D; /* FN along first row */
  B = c_matrix[class_total][0] - D; /* FP along first col */
  A = c_matrix[class_total][class_total] - B - C - D; /* TP ie rest */

  if ( debug ) 
    fprintf( stderr, "A = %f\nB = %f\nC = %f\n", A, B, C );

}



/* ----------------------------- MNI Header -----------------------------------
@NAME       : void write_report(int volume_index)
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: prints the confusion matrix to a file or stdout
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 6, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void write_report(int volume_index)
{

  int i, flag = 1;

  if (output_filename)
    flag = 0;

  /* write the created matrix to a file or stdout */

  if ( !stest ) {
    /* write the filenames being compared */
    fprintf(flag?stdout:matrep, "Confusion matrix for the following volumes :\n\n");
    fprintf(flag?stdout:matrep, "Volume1 (row-left) = %s \n", input_filename[0]);
    fprintf(flag?stdout:matrep, "Volume2 (col-top)  = %s \n\n", input_filename[volume_index]);
  }

  fprintf(flag?stdout:matrep, "Confusion matrix of class distribution \n\n ");

  fprintf(flag?stdout:matrep, "class ");
  for_less(i, 0, view_class_size)
    fprintf(flag?stdout:matrep, "%9d ",i);
  fprintf(flag?stdout:matrep, "    Total\n\n");

  for_less(r, 0, view_class_size) {

    fprintf(flag?stdout:matrep, "%5d  ",r);
    for_less(c, 0, view_class_size) 
      fprintf(flag?stdout:matrep, "%9d ", c_matrix[r][c]);
    fprintf(flag?stdout:matrep, "%9d\n", c_matrix[r][class_total]);

  }

  fprintf(flag?stdout:matrep, "\nTotal  ",r);
  for_less(c, 0, view_class_size) 
    fprintf(flag?stdout:matrep, "%9d ", c_matrix[class_total][c]);

  fprintf(flag?stdout:matrep, "%9d\n\n", c_matrix[class_total][class_total] );  

  fprintf(flag?stdout:matrep, "Collapsed confusion matrix statistics discribing different similarity measures: \n\n ");
  fprintf(flag?stdout:matrep, "Class   Kappa   CKappa  AZSM    Sensit. Error   Specif. Accuracy\n\n");
  for_less(i, 0, view_class_size) 
    fprintf(flag?stdout:matrep, "%5d\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\n",
	    i, 
	    c_table[i][PCLASS_KAPPA], 
	    c_table[i][CKAPPA], 
	    c_table[i][AZSM], 
	    c_table[i][SENSITIVITY], 
	    c_table[i][ERROR], 
	    c_table[i][SPECIFICITY], 
	    c_table[i][ACCURACY]);
  
  /* last line to give an overall performance measure */
  fprintf(flag?stdout:matrep, "\nTotal\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%7.4f\n",
	  sig_pclass_kappa,
	  total_kappa, 
	  total_azsm, 
	  total_hit_ratio, 
	  total_miss_ratio);

  /* print non collapsed kappa statistics, per class, brain and parenchyma */
  fprintf(flag?stdout:matrep, "\nCombined Kappa statistics: \n\n");

  /* print various kappas on separate lines to make parsing easier */
  fprintf(flag?stdout:matrep, "Total_Kappa=%6.4f    ",  total_kappa);
  fprintf(flag?stdout:matrep, "Brain_Kappa=%6.4f    ",  brain_kappa);
  fprintf(flag?stdout:matrep, "Paren_Kappa=%6.4f\n\n\n",  paren_kappa);
}



/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_mask_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: loads the mask in form of minc file
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 24, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void load_mask_volume(char *mask_file)
{

  int   *mask_volume_sizes;

  if (verbose) 
    fprintf (stderr, "Processing mask volume %s\n", mask_file);

  /* load the volume */
  status = input_volume(mask_file,
			3, 
			NULL,
			NC_BYTE, 
			FALSE, 
			0.0, 0.0,
			TRUE, 
			&mask_volume,
			(minc_input_options *) NULL ) ;
    
  if ( status != VIO_OK )
    exit(EXIT_FAILURE);
    
  /* reserve memory for masked volume sizes */
  ALLOC( mask_volume_sizes, VIO_MAX_DIMENSIONS );  

  /* get the mask volume sizes */
  get_volume_sizes(mask_volume, mask_volume_sizes);
    
  if (in_volume_sizes[0][VIO_X] != mask_volume_sizes[VIO_X]) {
    (void) fprintf(stderr,"Error - Mask Volume size mismatch in X dimension\n");
    exit(EXIT_FAILURE);
  }

  if (in_volume_sizes[0][VIO_Y] != mask_volume_sizes[VIO_Y]) {
    (void) fprintf(stderr,"Error - Mask Volume size mismatch in Y dimension\n");
    exit(EXIT_FAILURE);
  }
         
  if (in_volume_sizes[0][VIO_Z] != mask_volume_sizes[VIO_Z]) {
    (void) fprintf(stderr,"Error - Mask Volume size mismatch in Z dimension\n");
    exit(EXIT_FAILURE);
  }
 
            
} /* load_mask_volume */
        
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
@MODIFIED   : Mar 16, 1996 ( Vasco KOLLOKIAN)
---------------------------------------------------------------------------- */
int voxel_is_in_volume( int index, VIO_Real vox1, VIO_Real vox2, VIO_Real vox3)
{
 
  /* see if voxel is in specific volume marked by volume index */

  if ( vox1 < -0.5 || vox1 >= (VIO_Real) in_volume_sizes[index][0] - 0.5) {

    return FALSE;
  }
  
  else if ( vox2 < -0.5 || vox2 >= (VIO_Real) in_volume_sizes[index][1] - 0.5) {
    
    return FALSE;
  }
  
  else if ( vox3 < -0.5 || vox3 >= (VIO_Real) in_volume_sizes[index][2] - 0.5) {
    
    return FALSE;
  }
  else
    return TRUE;
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_voxel_to_voxel_transform
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: creates transforms from vol1 to vol2 in one step
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mar 16, 1996 ( David MACDONALD )
@MODIFIED   : 
---------------------------------------------------------------------------- */
void  create_voxel_to_voxel_transform(VIO_Volume volume1, VIO_Volume volume2,
				      VIO_General_transform  *v1_to_v2 )
{
    VIO_General_transform  *v1_to_world, *v2_to_world, world_to_v2;

    v1_to_world = get_voxel_to_world_transform( volume1 );
    v2_to_world = get_voxel_to_world_transform( volume2 );
    create_inverse_general_transform( v2_to_world, &world_to_v2 );

    concat_general_transforms( v1_to_world, &world_to_v2, v1_to_v2 );
}

