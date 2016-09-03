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
---------------------------------------------------------------------------- 
$RCSfile: classify.cc,v $
$Revision: 1.8 $
$Author: claude $
$Date: 2011-05-27 20:47:01 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : classify
@INPUT      : argc, argv -command line arguments
@OUTPUT     : 
@RETURNS    : error status
@DESCRIPTION: calls a classifier routine, feeds input data, stores result
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 8, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : $Log: classify.cc,v $
@MODIFIED   : Revision 1.8  2011-05-27 20:47:01  claude
@MODIFIED   : mostly fixes for memory bugs
@MODIFIED   :
@MODIFIED   : Revision 1.7  2006-09-26 18:42:31  claude
@MODIFIED   : switched initialization of fuzzy volume after cache_set stuff, otherwise minc will not write fuzzy volume - strange
@MODIFIED   :
@MODIFIED   : Revision 1.6  2006/09/14 19:47:25  claude
@MODIFIED   : initialize fuzzy volume to zero
@MODIFIED   :
@MODIFIED   : Revision 1.4  2006/09/14 19:04:23  claude
@MODIFIED   : initialize fuzzy volume to zero
@MODIFIED   :
@MODIFIED   : Revision 1.3  2005/02/11 20:16:15  bert
@MODIFIED   : Minor changes, primarily for compilation issues
@MODIFIED   :
@MODIFIED   : Revision 1.2  2002/03/20 22:25:06  jason
@MODIFIED   : added a sadly forgotten configure.ac and took out one debugging line
@MODIFIED   :
@MODIFIED   : Revision 1.1.1.1  2002/03/20 22:16:34  jason
@MODIFIED   : first autoconfiscated version that compiles under linux gcc 3
@MODIFIED   :
@MODIFIED   : Revision 1.9  1999/01/14 19:32:50  alex
@MODIFIED   : Added a few printf statements
@MODIFIED   :
@MODIFIED   : Revision 1.8  1998/02/26 19:52:12  alex
@MODIFIED   : C++-ified things
@MODIFIED   :
@MODIFIED   : Revision 1.7  1997/12/16 18:33:27  alex
@MODIFIED   : Added -output_range option
@MODIFIED   :
@MODIFIED   : Revision 1.6  1997/12/14 16:02:15  alex
@MODIFIED   : Fixed small bug in '-fuzzy all' option
@MODIFIED   :
@MODIFIED   : Revision 1.5  1997/12/13 20:26:04  alex
@MODIFIED   : Added '-fuzzy all' option
@MODIFIED   :
@MODIFIED   : Revision 1.4  1997/09/11 14:26:05  alex
@MODIFIED   : Various changes directed at g++ compatibility
@MODIFIED   :
@MODIFIED   : Revision 1.3  1997/02/14 21:07:22  alex
@MODIFIED   : Modified ann training, added -clobber switch
@MODIFIED   :
@MODIFIED   : Revision 1.2  1997/02/11 05:48:31  alex
@MODIFIED   : Changed history from overwriting to appending to the history of the first input file
@MODIFIED   :
@MODIFIED   : Revision 1.1.1.1  1997/02/11 00:06:41  alex
@MODIFIED   : Sources for classify, copied from Vasken Kollokian
@MODIFIED   :
 * Revision 1.26  1996/08/25  05:20:49  vasco
 * Added support for Bayes classification.
 * Added support for apriori probabilities to be fed to any classifier
 * Fixed bug, so that fuzzy volume clobbering is checked before classifications starts
 * Made cache the default behavior by adding -nocache switch
 *
 * Revision 1.25  1996/03/14  06:30:02  vasco
 * eliminated some unused variables.
 *
 * Revision 1.24  1996/03/14  05:45:50  vasco
 * fixed a bug in cleanup_memory, and re-commented it out.
 *
 * Revision 1.23  1996/03/06  00:25:41  vasco
 * The is another major revamp. I added Fuzzy C Means, and I modularized fcm min and hcm
 *
 * Revision 1.22  1996/01/25  06:10:03  vasco
 * fixed max_class_index bug in -load_tr option. and commented out
 * cleanup_memory()
 *
 * Revision 1.21  1995/12/07  15:48:19  vasco
 * Added Hard C Means classification, by supplying load_tag_volume
 * function, also added supervised flag to distinguish between the two
 * type of classifiers, supervised and unsupervised.
 * Replaced volume and image min and max by max_class_label.
 *
 * Revision 1.20  1995/11/17  07:50:54  vasco
 * fixed a bug in -dump, and added -fpath, a fuzzy pathname for fuzzvols.
 *
 * Revision 1.19  1995/11/06  06:39:41  vasco
 * Added memory cleanup routines, and chaged some volio calls.
 *
 * Revision 1.18  1995/10/27  05:34:39  vasco
 * added volume caching options, and set class_names to NULL
 *
 * Revision 1.17  1995/10/26  06:28:36  vasco
 * added -fuzzy "01001" string to select which class is expected
 * to write a fuzzy classification.  This gives a finer control
 * over what class results oin fuzzy classification.  This is
 * quite useful in MS since they require the fuzziness of only
 * one 'lesion' class.
 *
 * Revision 1.16  1995/10/23  02:03:46  vasco
 * Prints "*" while classifying each slice.
 *
 * Revision 1.15  1995/10/09  00:45:11  vasco
 * added an option to user specified voxel_min and voxel_man, image_min, image_max
 * fuzzy_voxel_min, fuzzy_voxel_max, fuzzy_image_min, fuzzy_image_max.
 * fixed the -help crashing. added some newlines for more readability, and added punctuations.
 *
 * Revision 1.14  1995/10/03  03:03:05  vasco
 * added -dump_features flag, and changed some stderr to stdout.
 *
 * Revision 1.13  1995/09/24  17:31:26  vasco
 * Fixed a bug in counting the correct number of remaining tags. Once
 * the tags that fell outside the volume were ignored. Also fixed a bug
 * (rather volio modification) of accessing volume size by volume->size.
 *
 * Revision 1.12  1995/09/20  06:07:11  vasco
 * I added class_names array that houses the class names, in the order they
 * appear in the tag file.
 *
 * Revision 1.11  1995/09/15  20:46:51  vasco
 * Added -mask -user_mask_class -user_mask_value, so that masks could be supported
 *
 * Revision 1.10  1995/09/15  08:16:34  vasco
 * added -fuzzy option with all of its intricacies and dependencies.
 * these include, writing of different data types in the fuzzy volumes
 * availability of fuzziness of implemented classifiers
 * creation, initialization and writing of fuzzy volumes.
 * Also, isolated in set_classifier_functions(), all function initializations
 * moved data allocation just before load_training to avoid duplication of code.
 *
 * Revision 1.9  1995/09/14  03:24:37  vasco
 * removed test_classifier option, and converted include files to <> form.
 *
 * Revision 1.8  1995/08/28  05:11:48  vasco
 * first version of including c4.5
 *
 * Revision 1.7  1995/08/28  03:25:51  vasco
 * converted to g++, and added class_names
 *
 * Revision 1.6  1995/07/12  05:12:23  vasco
 * added -parameter switch to load custom classifier parameters
 *
 * Revision 1.5  1995/07/12  03:49:36  vasco
 * added command history saving feature.
 *
 * Revision 1.4  1995/07/12  03:21:32  vasco
 * added training tag inside volume checks
 *
 * Revision 1.3  1995/06/10  17:59:50  vasco
 * Radical changes
 *
 * Revision 1.2  1995/05/31  01:29:22  vasco
 * added some comments and potential location of initialization routine.
 *
 * Revision 1.1  1995/05/31  00:58:46  vasco
 * Initial revision
 *
---------------------------------------------------------------------------- */
#include "config.h"

#include <iostream>		// (bert) added to force -lCio?
using namespace std;		// (bert) added
extern "C" {
#include <volume_io.h>
#include "time_stamp.h"
#include <ParseArgv.h>		// (bert) moved inside extern "C"
}
#include <unistd.h>
#include "class_protos.h"	// (bert) use quotes
#include "classify.h"		// (bert) use quotes
#include "mindist/mindist.h"	// (bert) use quotes
#include "knn/knn.h"		// (bert) use quotes
#include "ann/ann.h"		// (bert) use quotes

//extern "C" {
//#include <c4_5.h>
//}
#include "hcm/hcm.h"		// (bert) use quotes
#include "fcm/fcm.h"		// (bert) use quotes
#include "bayes/bayes.h"	// (bert) use quotes

/* MAIN */
int main(int argc, char *argv[])
{

  parse_arguments(argc, argv);
  
  allocate_memory(); 

  /* set the generic classifier functions */
  set_classifier_functions( classifier );


  /* the three following verification are to be called after
     set_classifier_functions() that is why they do not reside in
     parse_argv */

  if ( apriori && !use_apriori ) {

    (void) fprintf(stderr,"Chosen classifier doesn't support -apriori.\n");
    exit(EXIT_FAILURE);

  }

  if ( supervised && trainvol_filename && !apriori) {

    (void) fprintf(stderr,"Supervised classifiers don't accept tag volumes without -apriori switch.\n");
    exit(EXIT_FAILURE);

  }

  if ( !supervised && tagfile_filename ) {

    (void) fprintf(stderr,"Unsupervised classifiers don't accept tag files.\n");
    exit(EXIT_FAILURE);

  }



  /* check for fuzziness, ie if a classifier supports a fuzzy option */
  if ( fuzzy && !fuzzy_available) {

    fprintf(stderr, "The specified classifier does not support a fuzzy option\n");
    exit(EXIT_FAILURE);
  }


  /* set specific caching options */
  if ( cache_set && supervised && !apriori ) { 

    /* enable volume caching by setting threshold to 0 */ 
    set_n_bytes_cache_threshold(0);

    /* set max cache size, originally 0 to include only 1 block*/
    set_default_max_bytes_in_cache(max_cache_size);

    /* randomize caching to access for 1 voxel at a time which is good
       for reading tag files throughout the volumes - no user choice */
    block_sizes[0] = 1;
    block_sizes[1] = 1;
    block_sizes[2] = 1;
    set_default_cache_block_sizes(block_sizes);

  }
  else if ( cache_set && ( apriori || !supervised)) {

    /* if unsupervised classifier used, or apriori , enable slice caching */

    /* enable volume caching by setting threshold to 0 */ 
    set_n_bytes_cache_threshold(0);

    /* set max cache size, originally 0 to include only 1 block*/
    set_default_max_bytes_in_cache(max_cache_size);

    block_sizes[0] = user_block_sizes[0]; 
    block_sizes[1] = user_block_sizes[1];
    block_sizes[2] = user_block_sizes[2];
  
    /* And set slice caching for volumes to be created */
    set_default_cache_block_sizes(block_sizes);

  }

  else {

    /* disable caching all together*/
    set_n_bytes_cache_threshold(-1);
  }


  load_input_volumes();


  /* train supervised classifiers here */
  if ( supervised && tagfile_filename) {

    load_tag_file( tagfile_filename );
    create_feature_matrix_from_tagfile();
    init_training(param_filename);

    /* load the apriori volumes (with the training volume mechanism) 
       since it is exactly as loading FCM training volumes */
    if ( apriori ) {

      load_train_volumes(trainvol_filename);

      /* make sure number of classes is the same number of apriori volumes */
      if ( num_classes != num_train_vols ) {
	
	fprintf( stderr, "Number of classes <> # of apriori volumes\n");
	exit(EXIT_FAILURE);
      }

      /* allocate some memory for apriori vector - for apriori classifier */
      ALLOC(apriori_vector, num_classes);  

    }

    if ( classifier != KNN )
      /* knn has to train in each sample, so it is placed on it own */
      train();

  }

  /* train unsupervised classifiers here */
  if ( !supervised && trainvol_filename) {

    load_train_volumes(trainvol_filename);
    init_training(param_filename);
    train();
  }

  /* load previously trained classifier here - if proper switch is supplied */
  if ( load_train_filename) {

    init_training(param_filename);
    load_training(load_train_filename);
    max_class_index = num_classes;
  }

  /* save trained classifier here - if proper switch is supplied */
  if ( save_train_filename )
    save_training(save_train_filename);


  /* if not training, create, initialize, load mask, classify and write */
  if ( !train_only ) {

    if ( cache_set )     /* set slice caching options */
      convert_features_to_slice_caching();

    if ( fuzzy ) {

      decide_fuzzy_volumes();

    }

    create_empty_classified_volume();

    if ( mask_filename ) 
      load_mask_volume( mask_filename);

    classify_volume();

    write_classified_volume();

  }

  if ( fuzzy || classifier == FCM)
    write_fuzzy_volumes();

  /* this is for testing memory leaks, and may not be necessary to call */
  cleanup_memory();

  return(EXIT_SUCCESS);

} /* main */ 



/* ----------------------------- MNI Header -----------------------------------
@NAME       : parse_arguments
@INPUT      : argc, argv - command line arguments
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: parses command line arguments
@METHOD     : 
@GLOBALS    :   pname          - program name
                verbose        - show progress
                clobber        - overwrite file
                debug          - show debug info
                sub_command    - subclassifer arguments
                type           - type to read
                signtype       - sign to read
                num_features    - number of input volumes
                output_filename- the name of the classified volume 
@CALLS      : 
@CREATED    : February 6, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void parse_arguments(int argc, char *argv[])
{

  int   i, dst_idx = 0, src_idx = 0, buff_idx = 0;

  /* buffer to hold fuzzy_training volume names, limit 10, kludge! */
  char  fuzzy_train_buffer[10][256]; 

  pname = argv[0];

  /* form the history string - before parseargv is called */
  history = time_stamp(argc, argv);


  /* Call ParseArgv */
  if ( !ParseArgv(&argc, argv, argTable, 0) ) {   /* switches supplied are VIO_OK */ 

    if ( ( train_only || dump_features ) && (argc < 2) ) { 
      /* at least one volume has to be specified : output volume is optional
	 if the -train_only, or -dump_features switches are specified */

      (void) fprintf(stderr, "-train, -dump: specify at least one input volume.\n");
      (void) fprintf(stderr, "       %s [-help]\n\n", pname);
      exit(EXIT_FAILURE);

    } 

    if ( !(train_only || dump_features) && (argc < 3) ) {
      /* at least two volumes have to be specified - source and destination  */
      
      (void) fprintf(stderr, 
	     "\nUsage: %s <options> <infile1> [infile2] ... <outfile>\n", pname);
      (void) fprintf(stderr, "       %s [-help]\n\n", pname);
      
      exit(EXIT_FAILURE);
    }
  
  } /* ParseArgv */


  /* determine the number of volumes - ie features given to classify */
  if ( !( train_only || dump_features) ) 
    num_features = argc - 2;  /* count out progname and outvol */
  else
    num_features = argc - 1;  /* -train, -dump: count out only progname - no output*/

  /* Set clobber flags */
  if (clobber_all)
    clobber = clobber_training = clobber_fuzzy = TRUE;

  /* this is a kludge, to prevent -help from core dumping I guess it
     is a shortcoming of parseargv or my ignorance, rather impatient
     to deal with this problem. I guess the way ParseArgv is structured */
  
  if ( num_features == 0 ) 
    exit(EXIT_SUCCESS);


  /* reserve array for input volume filenames */
  ALLOC(input_filename, num_features); 


  /* make sure that the typed filenames exist */
  for_less(i, 0, num_features ) {

    if (!file_exists(argv[i+1])) {
      (void) fprintf(stderr, "filename `%s' not found. \n", argv[i+1]);
      exit(EXIT_FAILURE);
    }

    /* get the filename of the volume */
    input_filename[i] = argv[i+1];
    if ( debug >= 1)
      (void) fprintf(stdout, "input_filename[%d] = `%s'\n", i, input_filename[i]);

  } /* for_less */
  
  if ( !(train_only || dump_features) ) 
    output_filename = argv[num_features + 1];


  if ( debug >= 1 )
    (void) fprintf(stdout, "output_filename = `%s' \n", output_filename);
  

  if (!clobber && !train_only && !dump_features && file_exists(output_filename)) {

    (void) fprintf(stderr,"File `%s' exists ! \n", output_filename);
    (void) fprintf(stderr,"Use -clob_out to overwrite output classified volume.\n");
    exit(EXIT_FAILURE);

  }

  /* if a mask is specified, make sure it exists */
  if (mask_filename && !file_exists(mask_filename)) {

    (void) fprintf(stderr,"File `%s' doesn't exists ! \n", mask_filename);
    exit(EXIT_FAILURE);

  }

  if ( train_only && save_train_filename == NULL ) {

    (void) fprintf(stderr,"-train_only : please specify -save_train <file>\n");
    exit(EXIT_FAILURE);

  }


  if ( train_only && fuzzy ) {

    (void) fprintf(stderr,"-train_only and -fuzzy are mutually exclusive.\n");
    exit(EXIT_FAILURE);
  }

  if ( dump_features && fuzzy ) {

    (void) fprintf(stderr,"-dump_features and -fuzzy are mutually exclusive.\n");
    exit(EXIT_FAILURE);
  }


  /* for load_train_filename, checking perhaps should be done by the
     actual routine that is doing the loading, in that case it could
     take default values in case a filename is not provided */

  if ( load_train_filename && !file_exists(load_train_filename)  ) {

    (void) fprintf(stderr,"File `%s' doesn't exist !\n ", load_train_filename);
    exit(EXIT_FAILURE);

  }

  if (!clobber_training && 
      save_train_filename && 
      file_exists(save_train_filename)) {

    (void) fprintf(stderr,"File `%s' exists !\n", save_train_filename);
    (void) fprintf(stderr,"Use -clob_tr to overwrite output training file.\n");
    exit(EXIT_FAILURE);

  }

  if (!tagfile_filename  &&  
      !load_train_filename &&
      !tagvolume_buffer && (argc > 2 ) ) {

    /* argc > 2 is used to supress the message only if -help is used */
    (void) fprintf(stderr,"Specify one of  -tagfile or -load_train.\n");   
    exit(EXIT_FAILURE);

  }
  
  if ( apriori && !tagvolume_buffer ) {

    (void) fprintf(stderr,"Specify apriori volumes with -volume\n");   
    exit(EXIT_FAILURE);

  }


  /* make sure the specified tag file exists */
  if (tagfile_filename != NULL && !file_exists(tagfile_filename) ) {

    (void) fprintf(stderr,"File `%s' doesn't exist !\n ", tagfile_filename);
    exit(EXIT_FAILURE);
  }


  if (tagfile_filename != NULL && 
      load_train_filename != NULL ) {
    
    (void) fprintf(stderr,"-tagfile,  -load_train are mutually exlusive.\n");
    exit(EXIT_FAILURE);
    
  }
  
  /* parse the argument supplied by the -volume switch, which is a
     comma separated collection of training volume
     filenames, num_train_vols is determined by the number of training
     volumes supplied, */
  
  if ( tagvolume_buffer) { 
    
    if ( debug > 4 )
      fprintf( stdout, "tagvolume_buffer = %s\n", tagvolume_buffer);
    
    while ( 1 )  {
      
      /* while buf is not a ',' or null, advance pointer and copy character */
      while ( tagvolume_buffer[src_idx] != ',' &&  
	      tagvolume_buffer[src_idx] != '\0') {
	
	fuzzy_train_buffer[buff_idx][dst_idx++] = tagvolume_buffer[src_idx++];

      }
      
      /* as soon as you encounter a ',' or null, end target string */
      fuzzy_train_buffer[buff_idx][dst_idx] = '\0';
      
      /* reset target string pointer */
      dst_idx = 0;

      /* advance to next storage slot */
      buff_idx++;
      
      /* increment the number of training volumes */
      num_train_vols++;


      /* if you reach the end of the string, break and carry on */
      if (tagvolume_buffer[src_idx] == '\0')
	break;
      else
	src_idx++;

    }

    /* allocate space for the filename array */
    ALLOC( trainvol_filename, num_train_vols);

    /* allocate space, and copy buffer content */
    for_less( i, 0, num_train_vols ) {

      ALLOC( trainvol_filename[i], strlen(fuzzy_train_buffer[i])+1 );
      strcpy( trainvol_filename[i], fuzzy_train_buffer[i] );

      if ( debug > 4 ) 
	fprintf( stdout, "%s\n", trainvol_filename[i]);

      if ( !file_exists(trainvol_filename[i]) ) {

	fprintf(stderr,"File `%s' doesn't exist !\n ", trainvol_filename[i]);
	exit(EXIT_FAILURE);
      }

    }

    /* allocate memory for training volume pointers */
    ALLOC( train_volume, num_train_vols); 

	     
  } /* if ( trainvol_buffer ) */


} /* parse_arguments */



/* ----------------------------- MNI Header -----------------------------------
@NAME       : allocate_memory
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: allocates some memory for data structures
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 10, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void allocate_memory(void) 
{

  /* allocate some memory for feature vector */
  ALLOC(feature_vector, num_features);  

  /* allocate memory for first volume sizes */
  ALLOC(first_volume_sizes, VIO_MAX_DIMENSIONS); 

  /* allocate memory for volume pointers */
  ALLOC( in_volume, num_features); 

}




/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_input_volumes
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : in_volume[j]  - the number of input volumes
@CALLS      : 
@CREATED    : 
@MODIFIED   : Feb 10, 1996 ( Vasco KOLLOKIAN) 
---------------------------------------------------------------------------- */
void load_input_volumes(void)
{

  int   j;

    
  /* read the volumes in one by one */
  for_less( j, 0, num_features ) { 
    
    if (verbose) 
      fprintf (stdout, "Loading volume %s\n", input_filename[j]);
    
    /* load the volume */
    status = input_volume(input_filename[j], 
                          3, 
                          NULL, /* File_order_dimension_names, */
                          type, 
                          sign, 
                          0.0, 0.0,
                          TRUE, 
                          &in_volume[j],
                          (minc_input_options *) NULL ) ;
    
    if ( status != VIO_OK )
      exit(EXIT_FAILURE);

    /* if loading the very first volume, get its sizes, number of dims and dim 
       names, dim starts and dim steps */
    if ( j == 0 ) {

      int k; /* local counter */

      get_volume_sizes(in_volume[0], first_volume_sizes);
      first_volume_num_dims = get_volume_n_dimensions(in_volume[0]);
      first_volume_dim_names = get_volume_dimension_names(in_volume[0]);

      if ( debug > 2 ) {
	fprintf(stdout, "Vol number of dims. = %d\n", first_volume_num_dims);

	fprintf(stdout, "Vol dimension names = ");
	for_less ( k, 0, first_volume_num_dims ) 
	  fprintf(stdout, "%s ", first_volume_dim_names[k]);
	fprintf(stdout, "\n");
	
      }

    }

    /* if you have more than one volume, check to see if volumes are
    of same size in each dim and make sure the dimension orders are
    also the same. this is done by volume_size_is_ok() function */

    if ( j >= 1 &&  !volume_size_is_ok( in_volume[j] )){

      (void) fprintf(stderr,"in volume %s\n", input_filename[j]);
      exit(EXIT_FAILURE);
    }

  } /* for_less j */ 
            
} /* load_input_volumes */
        

/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_mask_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: loads the mask in form of minc file
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 15, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void load_mask_volume(char *mask_file)
{



  if (verbose) 
    fprintf (stdout, "Loading mask volume %s\n", mask_file);
    


  /* load the volume */
  status = input_volume(mask_file,
			3, 
			NULL, /* File_order_dimension_names, */
			type, 
			sign, 
			0.0, 0.0,
			TRUE, 
			&mask_volume,
			(minc_input_options *) NULL ) ;
    
  if ( status != VIO_OK )
    exit(EXIT_FAILURE);
    

  if ( !volume_size_is_ok( mask_volume) ) {

    fprintf( stderr, "in mask volume %s\n", mask_file);
    exit(EXIT_FAILURE);
  }
 
            
} /* load_mask_volume */
        



/* ----------------------------- MNI Header -----------------------------------
@NAME       : create_empty_classified_volume
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIN: cre and initialized an empty classified volume
@METHOD     : 
@GLOBAL   : cssified_volume
@CALLS      : 
@CREATED    : February 6, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
-------------------------------------------------------------------------- */
void create_empty_classified_volume(void)
{
  VIO_Real minval = (output_range[0] == -MAXDOUBLE) ? 0 : output_range[0];
  VIO_Real maxval = (output_range[1] == -MAXDOUBLE) ? max_class_index : output_range[1];

  /* create the classification volume here */   
  if (verbose) { 
    write(2,"Creating output volume\n",23);
  }

  classified_volume = copy_volume_definition(in_volume[0],
                                             NC_BYTE,
                                             FALSE,
                                             minval, maxval);

  set_volume_voxel_range(classified_volume, minval, maxval);
  set_volume_real_range(classified_volume, minval, maxval); 

  if ( cache_set ) {
    set_cache_output_volume_parameters(classified_volume, 
				       output_filename, 
				       NC_BYTE, 
				       FALSE, 
				       minval, maxval,
				       input_filename[0], 
				       history, 
				       (minc_output_options *) NULL ) ;


  }

} /* create_empty_classified_volume */


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
@MODIFIED   : 
---------------------------------------------------------------------------- */
void load_tag_file ( char *tag_filename )
{
  int *structure_ids;
  int i;

  if (verbose)
    (void) fprintf(stdout, "Loading  tagfile %s\n", tagfile_filename);

  /* tag file should be opened here */
  if ( input_tag_file(tag_filename, &n_tag_volumes, &num_samples,
		      &tags, NULL, NULL, &structure_ids, NULL, &labels ) != VIO_OK ) {

    fprintf(stderr, "Error reading the tag file.\n");
    exit(EXIT_FAILURE);
  }

  /*
   * Check that the labels were specified as expected. If not, substitute
   * the structure_id if available. If neither field is present in a tag 
   * point, then the tag file is invalid and we cannot proceed.
   */
  for (i = 0; i < num_samples; i++) {
    if (labels[i] == NULL) {
      if (structure_ids[i] >= 0) {
        char label[128];
        sprintf(label, "%d", structure_ids[i]);
        labels[i] = create_string(label);
      }
      else {
        fprintf(stderr, "Error in tag file, neither structure id nor label found.\n");
        exit(EXIT_FAILURE);
      }
    }
  }

  if ( n_tag_volumes == 2 ) 
    printf("Tag file contains two volumes, using the first one.\n");

  FREE(structure_ids);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : load_train_volumes
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : Feb 11, 1996 ( Vasco KOLLOKIAN )
---------------------------------------------------------------------------- */
void load_train_volumes(char **trainvol_filename)
{

  int i;
    
    
  /* load the training volumes */
  for_less( i, 0, num_train_vols ) { 

    if (verbose) 
      fprintf (stdout, "Loading tag volume %s\n", trainvol_filename[i]);

    status = input_volume(trainvol_filename[i],
			  3, 
			  NULL, /* File_order_dimension_names, */
			  NC_BYTE, 
			  FALSE, 
			  0.0, 0.0,
			  TRUE, 
			  &train_volume[i],
			  (minc_input_options *) NULL ) ;
    

    if ( status != VIO_OK )
      exit(EXIT_FAILURE);


    if ( !volume_size_is_ok( train_volume[i] ) ) {
      
      fprintf( stderr, "in training volume %s\n", trainvol_filename[i]);
      exit(EXIT_FAILURE);
    }

  } /* for_less( i, 0, num_train_vols ) */


            
} /* load_train_volume */
        

/* ----------------------------- MNI Header -----------------------------------
@NAME       : void create_feature_matrix_from_tagfile
@INPUT      : 
@OUTPUT     : 
@RETURNS    : create a feature matrix and a class vector for each point
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 29, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void create_feature_matrix_from_tagfile(void) 
{

  int    i, j;              /* counters i, over samples, j over volumes */
  VIO_Real   wx, wy, wz;        /* world x y z coordinates */
  VIO_Real   v1, v2, v3;        /* voxel x y z coordinates */
  VIO_Real   value;             /* voxel value to go into feature matrix */   
  int    num_adj_samples;   /* var to hold the number of adjusted samples */
  VIO_Real   *class_counter;    /* int array to count the number of classes in samples */
  
  
  /* set adjusted samples to zero, it will also index into feature_matrix 
   and class column.*/
  num_adj_samples = 0;

  if (verbose)
    (void) fprintf(stdout, "Creating feature matrix from tagfile\n");

  /* allocate memory for feature matrix space, as big as your samples size */
  VIO_ALLOC2D( feature_matrix, num_samples, num_features); 

  /* allocate momory for the class column, also as big as your samples size*/
  ALLOC( class_column, num_samples );

  /* initialize the feature_matrix */
  for_less( j, 0, num_features) 
    for_less( i, 0, num_samples) 
      feature_matrix[i][j] = 0.0;   

  /* start traversing the tag file */
  for_less( i, 0, num_samples) {

    /* get the world coordinate from the tag file */
    wx = tags[i][0];
    wy = tags[i][1];
    wz = tags[i][2];
    
    /* convert world into voxel coordinates using first volume */
    convert_3D_world_to_voxel(in_volume[0], wx, wy, wz, &v1, &v2, &v3);
    
    /* NOTE: check to see if the training point (voxel coordinate)
     falls in the volume. Since all the volumes have the same sizes,
     take volume 0, (that is why 'voxel_is_in_volume( v1, v2, v3)'. If
     later a -world switch is supplied, each tag should be tested in
     each feature volume.  At this point, check only for volume 0 */
 
    /* if voxel is in the volume, get the feature volumes values */
    if ( voxel_is_in_volume( v1, v2, v3)) {

      for_less( j, 0, num_features) {

        GET_VALUE_3D(value, in_volume[j], VIO_ROUND(v1), VIO_ROUND(v2), VIO_ROUND(v3));
        feature_matrix[num_adj_samples][j] = value;
      }

      /* convert the label into integer */
      class_column[num_adj_samples] = atoi( labels[i] );

      if ( debug > 5 ) {
        fprintf(stdout, "tag %d (%.1f %.1f %.1f) : ", i, v1, v2, v3); 
        for_less ( j, 0, num_features ) 
          fprintf(stdout, "%f, ", feature_matrix[num_adj_samples][j]); 
        fprintf(stdout, "-> %d.\n", class_column[num_adj_samples]);
      }

      /* keep track of highest class label to set voxel and image max */
      if ( max_class_index <  class_column[num_adj_samples] )
        max_class_index = class_column[num_adj_samples] ;
	
      num_adj_samples++;  /* adjust new sample size */

    }
    else {
      fprintf(stderr, "tag %d is not in the volume - ignoring\n", i);

    }

  } /* for_less( i, 0, num_samples) */

  /* after checking and ignoring 'out of volume tags', re-adjust the num_samples */
  num_samples = num_adj_samples;

  /*  COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES */
  /*  COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES */
  /*  COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES */
  /*  COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES */


  /* reserve space for counting classes */
  ALLOC( class_counter, num_samples+1 ); 

  /* initialize class_counter vector */
  for_less( i, 0, num_samples) 
    class_counter[i] = 0;


  /* count for the number of classes */
  for_less( i, 0, num_samples) {

    if ( debug > 4 ) {
      
      /* dump the feature matrix and class_column here */
      fprintf(stdout, "feature[%d] : ", i); 
      for_less ( j, 0, num_features ) 
        fprintf(stdout, "%f, ", feature_matrix[i][j]); 
      fprintf(stdout, "-> %d.\n", class_column[i]);
    }

    /* increase the classes counter */
    class_counter[ class_column[i] ]++; 

  } /* for_less */


  for_less( i, 0, num_samples) {

    if ( class_counter[i] != 0) 
      num_classes++;
  }
  
  /*  COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES */
  /*  COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES */
  /*  COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES */
  /*  COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES COUNTING CLASSES */

  /* reserve a character pointer ( char *) for each class name */
  ALLOC( class_name, num_classes ); 

  /* reserve a array to keep count of each class  */
  ALLOC( class_count, num_classes ); 

  /* initialize the count of each class to zero, and class_name to NULL */
  for_less( i, 0, num_classes) {

    class_count[i] = 0;
    class_name[i] = NULL;
  }

  /* now that you know how many classes exist, get rid of class_counter */
  FREE(class_counter);

  /* NOTE: Now get the name of each class as it occures in
     class_column[i] vector and copy its name from labels[i] array,
     into the class_name array. At the same time, establish a count of
     each class, and put it in the array class_count.  The way this is
     done is as follows: take each sample from class_column[i] vector.
     If its name exists in the class_name array, increase its
     class_count by one, if not, take the first available spot in the
     class_name[i] array, and copy its name there from labels[i]
     array, and increase its class count by one.  */

  /* for each sample in class_column */
  for_less( i, 0, num_samples) {

    /* and for each class */
    for_less( j, 0, num_classes ) {

      /* see if there is an empty spot (NULL), OR it exists (strcmp) . */
      if (labels[i] == NULL) {
	fprintf(stderr, "Error in tag labels; error reading tag file (possibly incompatible tag file format)\n");
	fprintf(stderr, "class_name[%d]: %s labels[%d]: %s\n", 
		j, class_name[j], i, labels[i]);
	exit(EXIT_FAILURE);
      }
	
      if ( class_name[j] == NULL || !strcmp( class_name[j], labels[i] ) ) {

	/* increase the class_count of that class in that spot j */
	class_count[j]++;

	/* if the spot was empty, that it means that this is a new class name,
	   copy the name of the class in the array, by reserving memory 1st*/
	if ( class_name[j] == NULL ) {

	  /* reserve space for the class name */
	  ALLOC( class_name[j], strlen( labels[i] ) + 1 );
	  strcpy( class_name[j], labels[i] );
	}

	/* go to the next element in the class_column, don't bother with the
	   rest of the class cases, since you just increase a count */
	break;
      }
    } /* for j */
  } /* for i */

      
  if ( debug >= 1 ) {
    
    fprintf(stdout, "num_classes = %d\n", num_classes);
    fprintf(stdout, "num_features = %d\n", num_features);
    fprintf(stdout, "num_samples = %d\n", num_samples);

    for_less( j, 0, num_classes) {

      fprintf(stdout, "class_name[%d]  = %s\n", j, class_name[j]);
      fprintf(stdout, "class_count[%d] = %d\n", j, class_count[j]);
    }

  }

  /* dump the feature_matrix and class_column if requested */
  if ( dump_features ) {

    for_less ( i, 0, num_samples ) {

      for_less ( j, 0, num_features ) 
	fprintf(stdout, "%f, ", feature_matrix[i][j]); 
      fprintf(stdout, "%d.\n", class_column[i]);

    }

    /* if you are only dumping features, exit after completion */
    if ( dump_features ) 
      exit( EXIT_SUCCESS );


  }

}      


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
@MODIFIED   : 
---------------------------------------------------------------------------- */
int voxel_is_in_volume( VIO_Real vox1, VIO_Real vox2, VIO_Real vox3)
{
 
  if ( vox1 < -0.5 || vox1 >= (VIO_Real) first_volume_sizes[0] - 0.5) {

    return FALSE;
  }
  
  else if ( vox2 < -0.5 || vox2 >= (VIO_Real) first_volume_sizes[1] - 0.5) {
    
    return FALSE;
  }
  
  else if ( vox3 < -0.5 || vox3 >= (VIO_Real) first_volume_sizes[2] - 0.5) {
    
    return FALSE;
  }
  else
    return TRUE;
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : convert_features_to_slice_caching(void)
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: convert the feature volume into slice caching
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Oct 26, 1995 ( Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void convert_features_to_slice_caching(void) 
{

  int i;

  /* enable slice access per user selected block size def {1,-1,-1} slice*/
  block_sizes[0] = user_block_sizes[0]; 
  block_sizes[1] = user_block_sizes[1];
  block_sizes[2] = user_block_sizes[2];
  
  /* this is for the previously opened random access {1,1,1} feature volumes */
  for_less ( i, 0, num_features ) 
    set_volume_cache_block_sizes(in_volume[i], block_sizes);

  /* And set slice caching for volumes to be created */
  set_default_cache_block_sizes(block_sizes);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : write_classified_volume()
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: write out the classified volume into a minc file.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 6, 1995 ( Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void write_classified_volume(void)
{
  VIO_Real minval = (output_range[0] == -MAXDOUBLE) ? 0 : output_range[0];
  VIO_Real maxval = (output_range[1] == -MAXDOUBLE) ? max_class_index : output_range[1];

  /* write the classified volume to a file */
  if (verbose) 
    fprintf(stdout,"\nWriting classified volume %s to file ...\n", output_filename);
   

  status = output_modified_volume( output_filename, 
				   NC_BYTE, 
				   FALSE, 
				   minval, maxval,
				   classified_volume, 
				   input_filename[0],
				   history, 
				   (minc_output_options *) NULL ) ;

  /* now delete the classified volume so that cache is flushed - as per David */
  delete_volume( classified_volume );

  if ( status != VIO_OK )
     exit(EXIT_FAILURE);
} /* write_classified_volume */

/* ----------------------------- MNI Header -----------------------------------
@NAME       : decide_fuzzy_volumes
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIN    : initialize create_fuzzy_volume array to decide which
              fuzzy volume to initialzie and create
@METHOD     : 
@GLOBAL     :
@CALLS      : 
@CREATED    : Oct 26, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
-------------------------------------------------------------------------- */
void decide_fuzzy_volumes(void){

  int k, num_fuzzy_vols;
  int fuzzy_all = !strcmp(fuzzy, "all");

  /* reserve some space for fuzzy volume creation flag */
  ALLOC( create_fuzzy_volume, num_classes );

  /* for each class, set fuzzy volume creation flag to false */
  for_less( k, 0, num_classes ) 
    create_fuzzy_volume[k] = 0;

  /* if a fuzzy path is not defined, make it the current directory */
  if ( fuzzy_path == NULL ) {
    
    ALLOC( fuzzy_path, 2);
    sprintf( fuzzy_path, ".");
  }
  
  /* reserve some space for fuzzy filename prefix and give default name */
  if ( !fuzzy_prefix ) {

    char fuzzy_name[] = "fuzzy_class";
    ALLOC( fuzzy_prefix, strlen(fuzzy_name)+1);
    strcpy( fuzzy_prefix, fuzzy_name );

  }

  /* reserve memory for the fuzzy filename */
  /* the 8 accounts for '/' '_' '.' 'm' 'n' 'c' '\0' 'spare' */

  VIO_ALLOC2D( fuzzy_filename, num_classes, (  strlen(fuzzy_path)
					 + strlen(fuzzy_prefix) 
					 + strlen(class_name[0]) + 8) );

  if ( fuzzy_all || ( strlen(fuzzy) > num_classes ) )
    num_fuzzy_vols = num_classes;
  else
    num_fuzzy_vols = strlen(fuzzy);

  /* now for each fuzzy class specification that is not zero, set the flag on.
     make sure you go as far as the shorter of the strlen(fuzzy) or num_classes */
  for_less( k, 0, num_fuzzy_vols ) {

    if ( fuzzy_all || ( fuzzy[k] != '0' ) ) {

      create_fuzzy_volume[k] = 1;

      /* compose the fuzzy volume filename by adding '.mnc' to class_name */
      sprintf( fuzzy_filename[k], "%s/%s_%s.mnc", fuzzy_path, 
	                                          fuzzy_prefix, 
	                                          class_name[k]);

      /* check to see if the volume already exists */
      if ( file_exists(fuzzy_filename[k]) && !clobber_fuzzy ) {

	fprintf(stderr, "%s exists! use -clob_fuzzy to overwrite.\n",fuzzy_filename[k]);
	exit(EXIT_FAILURE);
      }

    }

  }

  if ( debug >= 1 ) {

    for_less ( k, 0, num_fuzzy_vols ) 
      fprintf( stdout, "\ncreate_fuzzy_volume[%d] = %d", k, create_fuzzy_volume[k]);
    
    fprintf( stdout, "\n");

  }



} /* decide_fuzzy_volumes */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : initialize_fuzzy_volumes
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIN    :      create and initialize empty fuzzy  volumes
@METHOD     : 
@GLOBAL     :
@CALLS      : 
@CREATED    : Sep 14, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
-------------------------------------------------------------------------- */
void initialize_fuzzy_volumes(void)
{

  int    k;
    
  ALLOC( fuzzy_volume, num_classes );
    
  /* create the fuzzy volume here */   
  for_less( k, 0, num_classes) {
  
    printf("create_fuzzy_volume[k] = %d\n", create_fuzzy_volume[k]);
  
    if ( create_fuzzy_volume[k] != 0 ) {
  
      if (verbose) { 
        fprintf(stdout, "creating fuzzy volume for class %s\n", class_name[k]);
      }

      fuzzy_volume[k] = copy_volume_definition(in_volume[0],
                                               fuzzy_type,
                                               FALSE,
                                               fuzzy_voxel_min, 
                                               fuzzy_voxel_max);
  
      set_volume_voxel_range(fuzzy_volume[k], fuzzy_voxel_min, fuzzy_voxel_max); 
      set_volume_real_range(fuzzy_volume[k], fuzzy_image_min, fuzzy_image_max);
 
      if ( cache_set ) {
        set_cache_output_volume_parameters(fuzzy_volume[k],
                                           fuzzy_filename[k],
                                           fuzzy_type, 
                                           FALSE, 
                                           fuzzy_voxel_min, 
                                           fuzzy_voxel_max,
                                           input_filename[0], 
                                           history, 
                                           (minc_output_options *) NULL ) ;

      } /* if ( cache_set ) */

      VIO_Real bg = 0.0;    
      for( v1_ptr = 0; v1_ptr < first_volume_sizes[0]; v1_ptr++ ) {
        for( v2_ptr = 0; v2_ptr < first_volume_sizes[1]; v2_ptr++ ) {
          for( v3_ptr = 0; v3_ptr < first_volume_sizes[2]; v3_ptr++ ) {
            set_volume_real_value(fuzzy_volume[k],
                                  v1_ptr, v2_ptr, v3_ptr,
                                  0, 0, bg );
          } // end for v3_ptr
        } // end for v2_ptr
      } // end for v1_ptr

    } /* if ( create_fuzzy ...) */

  } /* for_less */
	
} /* initialize_fuzzy_volumes */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : write_fuzzy_volumes()
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: write out the fuzzy class confidence levels of the classified volume
              in separate volumes for each class. VIO_Volume names are that of
	      class_name[num_classes].
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep 14, 1995 ( Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void write_fuzzy_volumes(void)
{

  int i;    /* counter */


  /* write the fuzzy volumes to separate files */

  for_less( i, 0, num_classes ) {

    if ( create_fuzzy_volume[i] != 0 ) {

      if (verbose) 
	fprintf(stdout,"Writing fuzzy volume %s to file ...\n", fuzzy_filename[i]);
    
      status = output_modified_volume(fuzzy_filename[i],
			     fuzzy_type, 
			     FALSE, 
			     fuzzy_voxel_min, 
			     fuzzy_voxel_max,
			     fuzzy_volume[i], 
			     input_filename[0],
			     history, 
			     (minc_output_options *) NULL ) ;

      /* now delete the fuzzy volumes so that cache is flushed - as per David */
      delete_volume( fuzzy_volume[i] );
 
      if ( status != VIO_OK )
	exit(EXIT_FAILURE);
      
    } /* if ( create_fuzzy_volume ) */

  } /* for_less (i...) */
  
} /* write_fuzzy_volumes */

/* ----------------------------- MNI Header -----------------------------------
@NAME       : set_classifier_functions
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: set the various classifier functions into a generic format functions
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Sep. 14, 1995 ( Vasco KOLLOKIAN )
@MODIFIED   : 
---------------------------------------------------------------------------- */
void set_classifier_functions( int classifier_index )
{


  switch ( classifier_index ) {

  case MD:

    supervised    = TRUE;
    fuzzy_available = TRUE;
    use_apriori   = FALSE;
    init_training = minimum_distance_init_training;
    load_training = minimum_distance_load_training;
    save_training = minimum_distance_save_training;
    train         = minimum_distance_train_samples;
    classify      = minimum_distance_classify_sample;
    break;

  case KNN:

    supervised    = TRUE;
    fuzzy_available = TRUE;
    use_apriori   = FALSE;
    init_training = knn_init_training;
    load_training = knn_load_training;
    save_training = knn_save_training;
    train         = knn_train_samples;
    classify      = knn_classify_sample;
    break;

  case ANN:

    supervised    = TRUE;
    fuzzy_available = TRUE;
    use_apriori   = FALSE;
    init_training = ann_init_training;
    load_training = ann_load_training;
    save_training = ann_save_training;
    train         = ann_train_samples;
    classify      = ann_classify_sample;
    break;

//   case C45:

//     supervised    = TRUE;
//     fuzzy_available = TRUE;
//     use_apriori   = FALSE;
//     init_training = c4_5_init_training;
//     load_training = c4_5_load_training;
//     save_training = c4_5_save_training; 
//     train         = c4_5_train_samples;
//     classify      = c4_5_classify_sample;
//     break;

  case HCM:

    supervised    = FALSE;
    fuzzy_available = FALSE;
    use_apriori   = FALSE;
    init_training = hcm_init_training;
    load_training = hcm_load_training;
    save_training = hcm_save_training; 
    train         = hcm_train_samples;
    classify      = hcm_classify_sample;
    break;

  case FCM:

    supervised    = FALSE;
    fuzzy_available = TRUE;
    use_apriori   = FALSE;
    init_training = fcm_init_training;
    load_training = fcm_load_training;
    save_training = fcm_save_training; 
    train         = fcm_train_samples;
    classify      = fcm_classify_sample;
    break;

  case BAY:

    supervised    = TRUE;
    fuzzy_available = TRUE;
    use_apriori   = TRUE;
    init_training = bayesian_init_training;
    load_training = bayesian_load_training;
    save_training = bayesian_save_training; 
    train         = bayesian_train_samples;
    classify      = bayesian_classify_sample;
    break;


  }

}



/* ----------------------------- MNI Header -----------------------------------
@NAME       : classify_volume()
@INPUT      : in_volume[j]
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: classifies the n-dimensional feature volumes
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 29, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void classify_volume(void) 
{

  int    j, k;                  /* counters for sample and volume, fuzzy vols */
  int    class_val;             /* class of the classified sample */
  VIO_Real   mask_value;            /* holds value of mask to consider */


  /* reserve space for class probabilities and their labels */
  if ( fuzzy || classifier == FCM ) { 

    ALLOC( class_probs, num_classes );
    ALLOC( class_labels, num_classes );
    initialize_fuzzy_volumes();

  } 

  if (verbose)
    (void) fprintf(stdout, "Classifying volume... \n");

  /* take the size of volume 0, since they should all be the same */
  for_less (v1_ptr, 0, first_volume_sizes[0]) {

    if ( verbose ) 
      write(2,"*",1);

    for_less (v2_ptr, 0, first_volume_sizes[1]) {

      if ( debug > 10 )
	write(2,"-",1);
      
      for_less (v3_ptr, 0, first_volume_sizes[2]) {

	if ( debug >= 20 )
	  write(2,"+",1);

	/* extract the feature vector from volumes */
	for_less(j, 0, num_features)  {


	  /*GET_VALUE_3D(feature_vector[j],in_volume[j],v1_ptr,v2_ptr,v3_ptr);*/
	  feature_vector[j] = get_volume_real_value(in_volume[j], 
						    v1_ptr, 
						    v2_ptr, 
						    v3_ptr,
						    0,0);


	  
	}
	
	/* and if apriori switch set, then extract apriori vector 
	   values from the probability maps, which in this case are
	   residing in train_volume[j] */
	if ( apriori ) {

	  if ( debug > 31 ) write(2,"%",1);

	  for_less(j, 0, num_classes)  
	    apriori_vector[j] = get_volume_real_value(train_volume[j], 
						      v1_ptr, 
						      v2_ptr, 
						      v3_ptr,
						      0,0);
	}

	/* if a mask is chosen, ( by having a non-NULL mask_filename )
	   and mask value is >= some specified value, classify only
	   that voxel, otherwise set to to user defines class value,
	   default being zero */

	if ( mask_filename ) {
	  /* make sure the map is previously loaded, and has the same size */
	  /*GET_VALUE_3D( mask_value, mask_volume, v1_ptr, v2_ptr, v3_ptr);*/
	   mask_value = get_volume_real_value(mask_volume, 
					      v1_ptr, 
					      v2_ptr, 
					      v3_ptr,
					      0,0);
  
	  if (  mask_value <  user_mask_value ) {
	
	    class_val = user_mask_class;
	    goto WRITE_VOXEL;
	  }
	}
	 
	/*********************  S T A R T  ***************************/

	/* knn has to train on each sample before classifying */
	if (classifier == KNN)
	  train();

	if ( fuzzy || classifier == FCM) {
	      
	  classify(&class_val, class_probs, class_labels);
	      
	  /* for each fuzzy_volume[k], get the voxel value, and write it */
	  for_less ( k, 0, num_classes ) {

	    if ( debug > 31 ) 
	      fprintf(stdout, "%f, ", class_probs[k] ) ; 


	    /* debug value of 100 is to test the volume cache - is temp */
	    if ( create_fuzzy_volume[k] != 0 && debug < 100) {    
	      /*SET_VOXEL_3D(fuzzy_volume[k],v1_ptr,v2_ptr,v3_ptr,voxel_prob );*/
	      set_volume_real_value(fuzzy_volume[k],
				    v1_ptr,
				    v2_ptr,
				    v3_ptr,
				    0,0,
				    class_probs[k] );

	    } /* if */

	  } /* for_less ( k, 0, num_classes )*/

	  if ( debug > 31 ) 
	    fprintf(stdout, "\n "); 

	} /* if (fuzzy || classifier == FCM) */
	
	else /* crisp */
	
	  classify(&class_val, 0, 0);

	/**********************  E N D  ******************************/
	
      WRITE_VOXEL:
	/*SET_VOXEL_3D(classified_volume,v1_ptr,v2_ptr,v3_ptr,class_val);*/
	set_volume_real_value(classified_volume,
				    v1_ptr,
				    v2_ptr,
				    v3_ptr,
				    0,0,
				    class_val );

      } /* for_less v3_ptr */

    } /* for_less v2_ptr */

  } /* for_less v1_ptr */
}



/* ----------------------------- MNI Header -----------------------------------
@NAME       : cleanup_memory
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: cleans memory to test for leaks in debuging 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Nov 5, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void cleanup_memory(void) 
{

  int loop_idx;

  /* get rid of input filenames array (filenames themselves are
     from argv and cannot be freed) */
  FREE( input_filename );

  /* get rid of class names */
  for_less( loop_idx, 0, num_classes ) 
    FREE( class_name[loop_idx] );
  FREE( class_name );

  FREE( feature_vector );
  FREE( first_volume_sizes );

  if ( supervised && tagfile_filename) {

    VIO_FREE2D( feature_matrix );
    FREE( class_column );
  }

  FREE( class_count );

  for_less ( loop_idx, 0, first_volume_num_dims ) {
    FREE( first_volume_dim_names[loop_idx] );
  }
  FREE( first_volume_dim_names );
  first_volume_num_dims = 0;

  /* now clean the feature volumes */
  for_less( loop_idx, 0, num_features ) {

      delete_volume( in_volume[loop_idx] );
  }

  FREE( in_volume ); 


  /* now clean the mask volume */
  if ( mask_filename ) {

    delete_volume( mask_volume );
  }

  /* now free tag points */
  free_tag_points(n_tag_volumes, num_samples,
		  tags, NULL, NULL, NULL, NULL, labels );

} /* cleanup_memory(void)  */

/* ----------------------------- MNI Header -----------------------------------
@NAME       : volume_size_is_ok
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: verifies that volume sizes are VIO_OK
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 10, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
int volume_size_is_ok( VIO_Volume loaded_volume) 
{


  int    *loaded_volume_sizes;
  int    loaded_volume_num_dims;
  VIO_STR *loaded_volume_dim_names;
  
  /* allocate memory for first volume sizes */
  ALLOC(loaded_volume_sizes, VIO_MAX_DIMENSIONS); 

  /* get dim size, nums, order */
  get_volume_sizes(loaded_volume, loaded_volume_sizes);
  loaded_volume_num_dims = get_volume_n_dimensions(loaded_volume);
  loaded_volume_dim_names = get_volume_dimension_names(loaded_volume);
  
  if ( debug > 2 ) {

    int k; /* local counter */      
    
    fprintf(stdout, "Vol number of dims. = %d\n", loaded_volume_num_dims);
    
    fprintf(stdout, "Vol dimension names = ");
    for_less ( k, 0, loaded_volume_num_dims ) 
      fprintf(stdout, "%s ", loaded_volume_dim_names[k]);
    fprintf(stdout, "\n");
    
  }

  /* all the volume dimensions should be the same as the first volume */

  /* check for number of dimensions mismatch */
  if (loaded_volume_num_dims != first_volume_num_dims ) {

    (void) fprintf(stderr,"Error - Number of dimensions mismatch ");
    return FALSE;    
  }

     
  /* check for volume size mismatches */
  if (loaded_volume_sizes[VIO_X] != first_volume_sizes[VIO_X]) {

    (void) fprintf(stderr,"Error - Volume size mismatch in X dimension ");
    return FALSE;
  }

  if (loaded_volume_sizes[VIO_Y] != first_volume_sizes[VIO_Y]) {

    (void) fprintf(stderr,"Error - Volume size mismatch in Y dimension ");
    return FALSE;
  }

  if (loaded_volume_sizes[VIO_Z] != first_volume_sizes[VIO_Z]) {

    (void) fprintf(stderr,"Error - Volume size mismatch in Z dimension ");
    return FALSE;
  }


  /* check for dimensions order mismatch */
  if ( strcmp(loaded_volume_dim_names[0], first_volume_dim_names[0]) ) {
      
    (void) fprintf(stderr,"Error - First dimension order mismatch ");
    return FALSE;
  }

  /* if there are more than 1 dimension - check dim order of second dim */
  if ( loaded_volume_num_dims > 1) 
    if ( strcmp(loaded_volume_dim_names[1], first_volume_dim_names[1]) ) {
      
      (void) fprintf(stderr,"Error - Second dimension order mismatch ");
      return FALSE;
    }

  /* if there are more than 2 dimensions - check dim order of third dim*/
  if ( loaded_volume_num_dims > 2) 
    if ( strcmp(loaded_volume_dim_names[2], first_volume_dim_names[2]) ) {
      
      (void) fprintf(stderr,"Error - Third dimension order mismatch ");
      return FALSE;
    }

  /* if there are more then 3 dimensions - warn and die */
  if ( loaded_volume_num_dims > 3) {

    (void) fprintf(stderr,"Support is limited to 3 spatial dimensions ");
    exit(EXIT_FAILURE);
  }

  /* free the reserved memory locations */
  delete_dimension_names( loaded_volume, loaded_volume_dim_names );
  FREE(loaded_volume_sizes);

  return TRUE;
            
} /* volume_size_is_ok */
        
/* kate: tab-width 8; indent-mode cstyle; indent-width 8; replace-tabs on; */
