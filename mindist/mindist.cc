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
$RCSfile: mindist.cc,v $
$Revision: 1.1 $
$Author: jason $
$Date: 2002-03-20 22:16:35 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : mindist.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: minimum distance classifier
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 10, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : $Log: mindist.cc,v $
@MODIFIED   : Revision 1.1  2002-03-20 22:16:35  jason
@MODIFIED   : Initial revision
@MODIFIED   :
@MODIFIED   : Revision 1.1.1.1  1997/02/11 00:06:48  alex
@MODIFIED   : Sources for classify, copied from Vasken Kollokian
@MODIFIED   :
 * Revision 1.14  1996/08/29  03:55:22  vasco
 * added -static on some variables.
 *
 * Revision 1.13  1996/03/06  00:26:54  vasco
 * fixed the euclidian distance sqrt problem to avoid unneccesary calculations
 * Also modularized some of the routines.
 *
 * Revision 1.12  1996/01/25  06:10:27  vasco
 * fixed max_class_index bug in -load_tr option, by moving a copy
 * of alloc stuff to load_training, and making the one in init_training
 * conditional on num_classes.
 *
 * Revision 1.11  1995/12/07  16:04:28  vasco
 * Cleaned up the fuzzy section by moving the sqrt(delta...) there
 * Also Hand tested the algorithm, and its fuzzy versions as well.
 *
 * Revision 1.10  1995/12/05  03:49:24  vasco
 * replaced num_class_samples with class_count,
 * got rid of fuzzy_mathod, and decided to leave sqrt(delta_vector[i]) there.
 *
 * Revision 1.9  1995/11/14  20:35:43  vasco
 * added the new fuzziness from - nearest prototype algorithm
 *
 * Revision 1.8  1995/09/15  08:13:58  vasco
 * changes minimum = 320000 to DBL_MAX
 * moved all ALLOC calles into initialize_minimum_distance
 * replaced 9999 for feature occupancy, by INT_MAX
 * added a fuzzy version based on discussions with Dr. Shingal, (non-normalized version)
 *
 * Revision 1.7  1995/09/14  03:29:35  vasco
 * removed test_classifier
 *
 * Revision 1.6  1995/08/28  03:25:10  vasco
 * fixed a memory leak, in allocating feature_vector.
 *
 * Revision 1.3  1995/07/12  05:21:09  vasco
 * added header
 *
---------------------------------------------------------------------------- */
extern "C" {
#include <volume_io.h>
#include <limits.h>
}
#include "../class_globals.h" 


/* locally defined global variables */
static VIO_Real   **mean_feature_matrix;           /* matrix to reflect mean features */
static int    *mean_feature_class_vector;      /* vector to reflect mean feature classes */
static VIO_Real   *delta_vector;                   /* eucledian vector of each sample */
static VIO_Real   *fuzzy_min_vector;               /* fuzzy vector of each sample */
static VIO_Real   m;                               /* proxemity neighbourhood index */

/* ----------------------------- MNI Header -----------------------------------
@NAME       : minimum_distance_init_training
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
void minimum_distance_init_training(char *param_filename)
{


  FILE *min_par_file;

  /* if load_training is not used, then allocate space for struct */
  if ( !load_train_filename ) {

    /* reserve area for the delta_vector */
    ALLOC( delta_vector, num_classes );

    /* reserve area for the fuzzy_min_vector */
    ALLOC( fuzzy_min_vector, num_classes );

    /* reserve area for the mean feature matrix */
    VIO_ALLOC2D(mean_feature_matrix, num_classes, num_features);

    /* reserve area for the mean_feature_class_vector */
    ALLOC(mean_feature_class_vector, num_classes);

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
    min_par_file = fopen(param_filename, "r");
    
    if ( min_par_file == NULL) {

      fprintf(stderr, "Cannot open %s\n", param_filename);
      exit(EXIT_FAILURE);
    }
    
   /* scan for the neighbourhood proxemity number */
    fscanf( min_par_file, "m=%f\n", &m);
     
    fclose(min_par_file);

    if ( m <= 1.0 ) {
 
     fprintf(stderr, "Invalid m value : %f\n", m);
      exit(EXIT_FAILURE);
    }
 


  }

  if ( debug > 2) {

    fprintf(stdout, "m=%f\n", m);      
  }

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : minimum_distance_train_samples
@INPUT      : 
@OUTPUT     : 
@RETURNS    : ?
@DESCRIPTION: takes a feature matrix and trains a classifier on it.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 29, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void minimum_distance_train_samples(void)
{
  
  int      i, j, k;              /* counters - samples, features, classes */  


  if (verbose)
    (void) fprintf(stderr, "Training samples ... \n");


  /* initialize mean feature matrix, feature class vector & num of class samples */
  for_less( i, 0, num_classes) {

    mean_feature_class_vector[i] = INT_MAX;  /* big number to denote vacancy */

    for_less( j, 0, num_features) 
      mean_feature_matrix[i][j] = 0.0;

  }
   
  /* repeat for the total number of samples */
  for_less( i, 0, num_samples) {  

    /* repeat for the total number of classes */   
    for_less( k, 0, num_classes ) { 
    
      if ( mean_feature_class_vector[k] == INT_MAX ||      /* unoccupied spot */
	   mean_feature_class_vector[k] == class_column[i] ) { /* already has feat */

	for_less( j, 0, num_features ) 
	  mean_feature_matrix[k][j] += feature_matrix[i][j];

	/* if unoccupied, then assign class */
	if ( mean_feature_class_vector[k] == INT_MAX ) 
	  mean_feature_class_vector[k] = class_column[i];

	break;
      }

    } /* for k */

  } /* for i */


  if (verbose)
    (void) fprintf(stderr, "Generating mean feature matrix ...\n\n");

  for_less( i, 0, num_classes) 
    for_less( j, 0, num_features) 
      if (class_count[i] != 0)
	mean_feature_matrix[i][j] /=  (VIO_Real) class_count[i];

  if (debug > 2 ) {

    fprintf( stderr, "Printing mean_feature_matrix ...\n");

    for_less( i, 0, num_classes) {
      for_less( j, 0, num_features) 
	fprintf( stderr, "%f ", mean_feature_matrix[i][j]);
      fprintf( stderr, "%d\n", mean_feature_class_vector[i]);      
    }

    fprintf( stderr, "-----\n");

  }

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : minimum_distance_classify_sample
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
void minimum_distance_classify_sample(int *class_num, double *class_prob, int *class_labels)
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

  minimum =  10000000; //DBL_MAX;  /* the largest double value the system supports */
  class_index = 0;
    
  /* find the class with the MINIMUM DISTANCE */
  for_less ( k, 0,  num_classes ) 
    if ( delta_vector[k] < minimum ) {
      minimum = delta_vector[k];
      class_index = k;
    }
  

  *class_num = mean_feature_class_vector[class_index];

  /* THIS SECTION DEALS WITH FUZZY CLASSIFICATION ISSUES */
  /* THIS SECTION DEALS WITH FUZZY CLASSIFICATION ISSUES */
  /* THIS SECTION DEALS WITH FUZZY CLASSIFICATION ISSUES */

  if ( class_prob ) {

    VIO_Real sigma_fuzzy_min = 0.0;

    for_less( i, 0, num_classes) 
      if ( delta_vector[i] != 0 ) {

	/* here the sqrt(delta_vector[i]) is necessary to restore
	   the euclidian vector to || Xk - Zi || form, but in order
	   to avoid extra computations, the formula that could be
	   calculated as : 
	   
	   pow( sqrt(delta_vector[i]), (2/(m-1) ) )  is replaced by

	   pow( delta_vector[i], (1/(m-1) ) ) since

	   sqrt(x)^(2/(m-1)) === x^(1/(m-1))

	 */
	
	fuzzy_min_vector[i] = 1 / ( pow( delta_vector[i], (1/(m-1) ) ));
	sigma_fuzzy_min += fuzzy_min_vector[i];
      }
    
    for_less( i, 0, num_classes) {

      class_prob[i] = fuzzy_min_vector[i] / sigma_fuzzy_min;
      class_labels[i] = mean_feature_class_vector[i]; 
    }

  } /* if ( class_prob ) */

}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : minimum_distance_load_training
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
void minimum_distance_load_training(char *load_train_filename)
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
    fprintf( stderr, "num_of_features = %d\n", num_features);
    fprintf( stderr, "num_of_classes = %d\n", num_classes);
  }


  /* reserve area for the delta_vector */
  ALLOC( delta_vector, num_classes );

  /* reserve area for the fuzzy_min_vector */
  ALLOC( fuzzy_min_vector, num_classes );

  /* reserve area for the mean feature matrix */
  VIO_ALLOC2D(mean_feature_matrix, num_classes, num_features);

  /* reserve area for the mean_feature_class_vector */
  ALLOC(mean_feature_class_vector, num_classes);
 
 /* reserve a character pointer ( char *) for each class name */
  ALLOC( class_name, num_classes ); 
 
  /* load feature matrix from learn file */
  if (verbose)
    fprintf(stderr, "Loading the training file...\n");

  
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

    fprintf( stderr, "Printing mean_feature_matrix ...\n");
    for_less( i, 0, num_classes) {
      for_less( j, 0, num_features) 
	fprintf( stderr, "%f ", mean_feature_matrix[i][j]);
      fprintf( stderr, "%d\n", mean_feature_class_vector[i]);      
    }

    fprintf( stderr, "-----\n");
  }

  fclose(learn_file);

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : minimum_distance_save_training
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
void minimum_distance_save_training(char *save_train_filename)
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

