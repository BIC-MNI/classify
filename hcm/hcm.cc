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
$RCSfile: hcm.cc,v $
$Revision: 1.1 $
$Author: jason $
$Date: 2002-03-20 22:16:35 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : hcm.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: hard c means unsupervised classifier
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Dec. 2, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : $Log: hcm.cc,v $
@MODIFIED   : Revision 1.1  2002-03-20 22:16:35  jason
@MODIFIED   : Initial revision
@MODIFIED   :
@MODIFIED   : Revision 1.1.1.1  1997/02/11 00:06:48  alex
@MODIFIED   : Sources for classify, copied from Vasken Kollokian
@MODIFIED   :
 * Revision 1.4  1996/08/29  03:53:38  vasco
 * I forgot
 *
 * Revision 1.3  1996/03/14  06:31:21  vasco
 * added static on some multiply defined variables.
 *
 * Revision 1.2  1996/03/09  18:24:31  vasco
 * removed the nested comment, - SGI's complained
 *
 * Revision 1.1  1995/12/07  16:01:04  vasco
 * Initial revision
 *
---------------------------------------------------------------------------- */
extern "C" {
#include <volume_io.h>
#include <limits.h>
}
#include "../unsuper_globals.h" 
#include <unistd.h>
#include <sys/types.h>

/* locally defined prototypes */
void hcm_count_classes(void);

/* locally defined global variables */
static VIO_Real   **mean_feature_matrix;     /* matrix to reflect mean features */
static int    *mean_feature_class_vector;/* vector to reflect mean feature classes */
static VIO_Real   *delta_vector;             /* eucledian vector of each sample */
VIO_Real   *fuzzy_hcm_vector;               /* fuzzy vector of each sample */
int    v1, v2, v3;                      /* voxel pointers in x, y, and z direction */
static VIO_Real   m;                        /* proxemity neighbourhood index */

/* ----------------------------- MNI Header -----------------------------------
@NAME       : hcm_init_training
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Dec 2, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void hcm_init_training(char *param_filename)
{


  FILE *hcm_par_file;

  /* if load_training is not used, then allocate space for struct */
  if ( !load_train_filename ) {

    /* count num_classes, and set their names */
    hcm_count_classes();
    
    /* reserve area for the delta_vector */
    ALLOC( delta_vector, num_classes );

    /* reserve area for the fuzzy_hcm_vector */
    ALLOC( fuzzy_hcm_vector, num_classes );

    /* reserve area for the mean feature matrix */
    VIO_ALLOC2D(mean_feature_matrix, num_classes, num_features);

    /* reserve area for the mean_feature_class_vector */
    ALLOC(mean_feature_class_vector, num_classes);

    /* reserve area for the number of class samples */
    ALLOC(class_count, num_classes);

  }

  /* check to see if the filename is there */
  if ( param_filename && !file_exists(param_filename)  ) {

    (void) fprintf(stderr,"File `%s' doesn't exist !\n ", param_filename);
    exit(EXIT_FAILURE);

  }

  if ( !param_filename ) {

    /* take the default here, in case a param file is not given */
    m = 2.0;

  }
  else {
    
    if (verbose) 
      fprintf(stdout, "Loading the parameter file %s\n", param_filename);      
    
    /* open the parameter file, and read the values  */
    hcm_par_file = fopen(param_filename, "r");
    
    if ( hcm_par_file == NULL) {

      fprintf(stderr, "Cannot open %s\n", param_filename);
      exit(EXIT_FAILURE);
    }
    
   /* scan for the neighbourhood proxemity number */
    fscanf( hcm_par_file, "m=%lf\n", &m);
     
    fclose(hcm_par_file);

    if ( m <= 1.0 ) {
 
     fprintf(stderr, "Invalid m value : %f\n", m);
      exit(EXIT_FAILURE);
    }
 
  } /* else */

  if ( debug > 2) {

    fprintf(stdout, "m=%lf\n", m);      
  }

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : hcm_train_samples
@INPUT      : 
@OUTPUT     : 
@RETURNS    : ?
@DESCRIPTION: takes a feature matrix and trains a classifier on it.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Dec 2, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void hcm_train_samples(void)
{
  
  int      i, j, k;              /* counters - samples, features, classes */  
  char     name_of_class[4];     /* temp space of char representation of 0-255 */
  VIO_Real     training_class;       /* traing voxel class */
  

  if (verbose)
    (void) fprintf(stdout, "Training on entire volume ... \n");


  /* initialize mean feature matrix, feature class vector & num of class samples */
  for_less( i, 0, num_classes) {

    class_count[i] = 0;
    mean_feature_class_vector[i] = INT_MAX;  /* big number to denote vacancy */

    for_less( j, 0, num_features) 
      mean_feature_matrix[i][j] = 0.0;

  }
   
  /* here the first_volume_sizes also denotes the size of the 
     feature volumes, and the tag volumes they should be identical */

  if ( num_train_vols > 1 ) {

    fprintf( stderr, "Hard C Means Classifier takes only one training volume\n");
    exit(EXIT_FAILURE);
  }

  /* repeat for the total number of samples - the entire volume */
  for_less( v1, 0, first_volume_sizes[0] ) {

    if ( verbose ) 
      write(2, "*", 1);

    for_less( v2, 0, first_volume_sizes[1] )
      for_less( v3, 0, first_volume_sizes[2] ) {
	
	/* get the training_class from the first train_volume */
	training_class =  get_volume_real_value(train_volume[0], 
						v1, 
						v2, 
						v3,
						0,0);
	
	if ( debug > 8 ) 
	  fprintf( stdout, "training voxel = %f\n", training_class);

	/* repeat for the total number of classes */   
	for_less( k, 0, num_classes ) { 
	  
	  if ( /* unoccupied spot */
	      mean_feature_class_vector[k] == INT_MAX ||      
	      /* already used */
	      mean_feature_class_vector[k] == (int) training_class ) { 
	    
	    /* extract a mean_feature_matrix row from volumes */
	    for_less(j, 0, num_features)  
	      mean_feature_matrix[k][j] += get_volume_real_value(in_volume[j],
								 v1, 
								 v2, 
								 v3,
								 0,0);
	    
	    class_count[k]++;
	    
	    /* if unoccupied, then assign class */
	    if ( mean_feature_class_vector[k] == INT_MAX ) 
	      mean_feature_class_vector[k] = (int) training_class;
	    
	    break;
	  }
	  
	} /* for k */
	
      } /* for_less v3 */
  } /* for_less v1 */

  if (verbose)
    (void) fprintf(stdout, "\nGenerating mean feature matrix, and class names.\n\n");

  for_less( i, 0, num_classes) {

    for_less( j, 0, num_features) 
      if (class_count[i] != 0)
	mean_feature_matrix[i][j] /=  (VIO_Real) class_count[i];

    /* set the class names - itoa the class number from mean_fcv */
    sprintf( name_of_class, "%d", mean_feature_class_vector[i]);

    /* reserve space for the class name */
    ALLOC( class_name[i], strlen( name_of_class ) );
    strcpy( class_name[i], name_of_class );

  }


  if (debug > 2 ) {

    fprintf( stdout, "Printing mean_feature_matrix ...\n");

    for_less( i, 0, num_classes) {

      fprintf( stdout, "class=%d, idx=%d, name=%s  ", i,
	                                              mean_feature_class_vector[i],
	                                              class_name[i]);
      for_less( j, 0, num_features) 
	fprintf( stdout, "%f ", mean_feature_matrix[i][j]);
      fprintf( stdout, "\n"); 

    }

    fprintf( stdout, "-----------\n");

  }

  /* get rid of the training volume, once finished */
  delete_volume(train_volume[0]);
  
} /* hcm_train_samples(void) */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : hcm_load_training
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Jun 3, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void hcm_load_training(char *load_train_filename)
{

  int i, j, k;
  VIO_Real feature_value;
  FILE *learn_file;


  /* open the input training  file */
  learn_file = fopen(load_train_filename, "r");

  if ( learn_file == NULL) {
    fprintf(stderr, "Cannot open %s\n", load_train_filename);
    exit(EXIT_FAILURE);
  }

  /* scan for the number of features */
  fscanf( learn_file, "num_of_features = %d\n", &num_features);

  /* scan for the max number of classes M\n*/
  fscanf( learn_file, "num_of_classes = %d\n", &num_classes);

  if (debug) {
    fprintf( stdout, "num_of_features = %d\n", num_features);
    fprintf( stdout, "num_of_classes = %d\n", num_classes);
  }


  /* reserve area for the delta_vector */
  ALLOC( delta_vector, num_classes );

  /* reserve area for the fuzzy_hcm_vector */
  ALLOC( fuzzy_hcm_vector, num_classes );

  /* reserve area for the mean feature matrix */
  VIO_ALLOC2D(mean_feature_matrix, num_classes, num_features);

  /* reserve area for the mean_feature_class_vector */
  ALLOC(mean_feature_class_vector, num_classes);

  /* reserve a character pointer ( char *) for each class name */
  ALLOC( class_name, num_classes ); 
  
  /* reserve area for the number of class samples */
  ALLOC(class_count, num_classes);

  /* load feature matrix from learn file */
  if (verbose)
    fprintf(stdout, "Loading the training file...\n");

  
  for_less( k, 0, num_classes) {
    
    fscanf(learn_file, "class = %d\n", &mean_feature_class_vector[k]);    

    /* allocate some space for each class_name[k], 5 for ex and
       simulate itoa, this could be replaced by a more elegant looking
       code, when I have time */
    ALLOC( class_name[k], 5);
    sprintf( class_name[k], "%d", mean_feature_class_vector[k]); 

    for_less( j, 0, num_features) {
      fscanf(learn_file, "%lf\n", &feature_value);
      mean_feature_matrix[k][j] = feature_value ;
    }

  }

  if (debug > 2 ) {

    fprintf( stdout, "Printing mean_feature_matrix ...\n");
    for_less( i, 0, num_classes) {
      for_less( j, 0, num_features) 
	fprintf( stdout, "%f ", mean_feature_matrix[i][j]);
      fprintf( stdout, "%d\n", mean_feature_class_vector[i]);      
    }

    fprintf( stdout, "-----\n");
  }

  fclose(learn_file);

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : hcm_count_classes_and_set_names
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: sets num_classes
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Dec. 3 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void hcm_count_classes(void)
{
  
  int      i, j = 0;       /* counter - classes */  
  VIO_Real     voxel_class;    /* class of voxel in training volume */
  int      bin_size = 5;   /* var to hold number of bins */
  long     *count_bin;     /* 1D array to hold the number of voxels / class */

  
  /* Allocate some memory for the count_bin array */
  ALLOC( count_bin, bin_size ) ;

  /* initialize the count bins */
  for_less ( i, 0, bin_size ) 
    count_bin[i] = 0;

  if (verbose)
    (void) fprintf(stdout, "Counting number of classes in training volume ... \n");

  /* repeat for the total number of samples - the entire volume */
  for_less( v1, 0, first_volume_sizes[0] ) 
    for_less( v2, 0, first_volume_sizes[1] )
      for_less( v3, 0, first_volume_sizes[2] ) {
	
	/* get the training_class from the first train volume */
	voxel_class =  get_volume_real_value(train_volume[0], 
					     v1, 
					     v2, 
					     v3,
					     0,0);
	
	/* this method won't handle cases where the number of 
	   classes is greater than 10, since the bin */


	/* if you exceed the bin size, increase apropriately */
	if ( (int) voxel_class > ( bin_size - 1 ) ) {
	  
	  int increment;

	  /* see how much we have to increase count_bin */
	  increment = (int) voxel_class - ( bin_size - 1 );

	  SET_ARRAY_SIZE( count_bin, bin_size, bin_size+increment, 100 );

	  /* now initialize the rest of the slots in count_bin */
	  for_less( i, bin_size, bin_size+increment ) 
	    count_bin[i] = 0;

	  bin_size += increment;

	}
	
	/* increase the bin count for that voxel */
	count_bin[ (int) voxel_class ]++;

      } /* for_less v3 */

  if ( debug >= 3 ) {

    for_less( i, 0, bin_size) 
      fprintf( stdout, "count_bin[%d] = %ld\n", i, count_bin[i] );

  }

  /* now count the number of classes, that have non-zero bin numbers */
  for_less( i, 0, bin_size)
    if ( count_bin[i] != 0) 
      num_classes++;

  /* reserve a character pointer ( char *) for each class name */
  ALLOC( class_name, num_classes ); 

  /* set the value of max_class_index  */
  for_less( i, 0, bin_size) {

    if ( count_bin[i] != 0) {

      if ( max_class_index <  i )
	max_class_index = i ;

    }

  }

  /* free up stuff not needed */
  FREE(count_bin);


} /* count_number_of_classes(...) */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : hcm_classify_sample
@INPUT      : 
@OUTPUT     : 
@RETURNS    : sample class
@DESCRIPTION: given a feature vector and its size, return a class
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 29, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void hcm_classify_sample(int *class_num, double *class_prob, int *class_labels)
{

  int           i, j, k;                       /* counters */
  VIO_Real          temp;                          /* temporary variable */
  VIO_Real          minimum;                       /* minimum distance class buffer */
  int           class_index;


  /* somehow load the learned classification from file and use it in this routine,
     but in order to do this there has to be an initialization routine that will do
     so and call this routine thereafter */

  /* initialize delta vector */
  for_less ( i, 0, num_classes ) 
    delta_vector[i] = 0.0;
  
  /* calculate Euclidian distance */
  for_less ( i, 0,  num_classes ) {

    for_less ( j, 0, num_features) {

      temp = feature_vector[j] - mean_feature_matrix[i][j]; /* X - Z */
      delta_vector[i] += temp * temp;                       /* Sig ( X - Z ) ^2 */
    }

  }

  //  minimum = DBL_MAX;  /* the largest double value the system supports */
  minimum = 10000000;
  class_index = 0;
    
  /* find the class with the MINIMUM DISTANCE */
  for_less ( k, 0,  num_classes ) 
    if ( delta_vector[k] < minimum ) {
      minimum = delta_vector[k];
      class_index = k;
    }
  
  if ( debug > 8 ) {
    fprintf( stdout, "index = %d, class = %d\n", class_index, 
	             mean_feature_class_vector[class_index]);
  }

  *class_num = mean_feature_class_vector[class_index];

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : hcm_save_training
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Jun 3, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void hcm_save_training(char *save_train_filename)
{

  int i, j;                 /* counters */
  FILE *train_file;         /* save filename */


  train_file = fopen(save_train_filename, "w");


  /* see if the file could be opened */
  if ( train_file == NULL) {

    printf("Cannot open %s\n", save_train_filename);
    exit(EXIT_FAILURE);

  }

  /* store in the file the number of features trained on */
  fprintf(train_file, "num_of_features = %d\n", num_features);

  /* store in the file the max number of classes trained on */
  fprintf(train_file, "num_of_classes = %d\n", num_classes);
  
  for_less ( i, 0, num_classes ) {
    fprintf(train_file, "class = %d\n", mean_feature_class_vector[i]); 
    for_less ( j, 0, num_features )
      fprintf(train_file, "%lf ", mean_feature_matrix[i][j]); 
    fprintf(train_file, "\n");
  }

  fclose(train_file);

}
