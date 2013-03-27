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
$RCSfile: fcm.cc,v $
$Revision: 1.1 $
$Author: jason $
$Date: 2002-03-20 22:16:34 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : fcm.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: fuzzy c means unsupervised classifier
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 12, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : $Log: fcm.cc,v $
@MODIFIED   : Revision 1.1  2002-03-20 22:16:34  jason
@MODIFIED   : Initial revision
@MODIFIED   :
@MODIFIED   : Revision 1.1.1.1  1997/02/11 00:06:48  alex
@MODIFIED   : Sources for classify, copied from Vasken Kollokian
@MODIFIED   :
 * Revision 1.3  1996/08/29  03:54:18  vasco
 * I forgot.
 *
 * Revision 1.2  1996/03/14  06:31:30  vasco
 * added static on some multiply defined variables.
 * removed some unused variables.
 *
 * Revision 1.1  1996/03/06  00:28:39  vasco
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
void fcm_count_classes_and_set_names(void);

/* locally defined global variables */
static VIO_Real   **mean_feature_matrix;  /* matrix to reflect mean features */
static int    *mean_feature_class_vector;/* vector to reflect mean feature classes */
static VIO_Real   *euclidian_vector;         /* eucledian vector of each sample */
VIO_Real   *fuzzy_fcm_vector;               /* fuzzy vector of each sample */
static int    v1, v2, v3;               /* voxel pointers in x, y, and z direction */
static VIO_Real   m;                        /* proxemity neighbourhood index */

/* ----------------------------- MNI Header -----------------------------------
@NAME       : fcm_init_training
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 12, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void fcm_init_training(char *param_filename)
{


  FILE *fcm_par_file;

  /* if load_training is not used, then allocate space for struct */
  if ( !load_train_filename ) {

    /* count num_classes, and set their names */
    fcm_count_classes_and_set_names();

    /* reserve area for the euclidianvector */
    ALLOC( euclidian_vector, num_classes );

    /* reserve area for the fuzzy_fcm_vector */
    ALLOC( fuzzy_fcm_vector, num_classes );

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
    fcm_par_file = fopen(param_filename, "r");
    
    if ( fcm_par_file == NULL) {

      fprintf(stderr, "Cannot open %s\n", param_filename);
      exit(EXIT_FAILURE);
    }
    
   /* scan for the neighbourhood proxemity number */
    fscanf( fcm_par_file, "m=%lf\n", &m);
     
    fclose(fcm_par_file);

    if ( m <= 1.0 ) {
 
     fprintf(stderr, "Invalid m value : %f\n", m);
      exit(EXIT_FAILURE);
    }
 
  } /* else */

  if ( debug > 2) {

    fprintf(stdout, "m=%f\n", m);      
  }

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : fcm_train_samples
@INPUT      : 
@OUTPUT     : 
@RETURNS    : ?
@DESCRIPTION: takes a feature matrix and trains a classifier on it.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 12, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void fcm_train_samples(void)
{
  
  int      i, j;              /* counters - samples, features, classes */  

  VIO_Real     u_ik, u_ikm, sig_u_ikm;


  /* initialize mean feature matrix, feature class vector & num of class samples */
  for_less( i, 0, num_classes) {

    for_less( j, 0, num_features) 
      mean_feature_matrix[i][j] = 0.0;

  }
   
 /* train on all input training volumes, one at a time, and one per class */
  for_less ( i, 0, num_classes ) {

    /* reset the u variable and the sig total */
    u_ik = 0.0;
    sig_u_ikm = 0.0;

    if (verbose)
      (void) fprintf(stdout, "Training on volume %s\n", trainvol_filename[i]);

    /* repeat for the total number of samples - the entire training volume */
    for_less( v1, 0, first_volume_sizes[0] ) {

      if ( verbose ) 
	write(2, "*", 1);

      for_less( v2, 0, first_volume_sizes[1] ) {

	if ( debug > 13 ) 
	  write(2, "#", 1);
	
	for_less( v3, 0, first_volume_sizes[2] ) {
	
	  /* get u_ik for class i, from the ith training volume */
	  u_ik =  get_volume_real_value(train_volume[i], 
					v1, 
					v2, 
					v3,
					0,0);
	
	  /* take it to the mth power  */
	  u_ikm = pow(u_ik, m);

	  /* calculate and store the total - Sig[ (uik)^m ] */
	  sig_u_ikm += u_ikm;

	  /* Sig [ (uik)^m * Xk ] for each and every feature */
	  for_less(j, 0, num_features)  
	    mean_feature_matrix[i][j] += u_ikm * get_volume_real_value(in_volume[j],
								       v1, 
								       v2, 
								       v3,
								       0,0);
	  
	    
	} /* for_less v3 */
      } /* for_less v2 */
    } /* for_less v1 */

    if (verbose)
      (void) fprintf(stdout, "\nGenerating class centroid for class %s\n",
		     class_name[i]);

    /* calculate for each Xk, Vi = (Sig[ (uik)^m *Xk ] ) / ( Sig[ (uik)^m ] ) */
    for_less(j, 0, num_features)  
      mean_feature_matrix[i][j] /= sig_u_ikm ;
    
    /* get rid of the training volume, once finished */
    delete_volume(train_volume[i]);

  } /* for_less ( i, 0, num_classes ) */

  if (debug > 2 ) {

    fprintf( stdout, "Printing class centroids  ...\n");

    for_less( i, 0, num_classes) {
      for_less( j, 0, num_features) 
	fprintf( stdout, "%f ", mean_feature_matrix[i][j]);
      fprintf( stdout, "%s\n", class_name[i]);      
    }

    fprintf( stdout, "-----------\n");

  }
  
} /* fcm_train_samples(void) */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : fcm_load_training
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 12, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void fcm_load_training(char *load_train_filename)
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


  /* reserve area for the euclidianvector */
  ALLOC( euclidian_vector, num_classes );

  /* reserve area for the fuzzy_fcm_vector */
  ALLOC( fuzzy_fcm_vector, num_classes );

  /* reserve area for the mean feature matrix */
  VIO_ALLOC2D(mean_feature_matrix, num_classes, num_features);

  /* reserve area for the mean_feature_class_vector */
  ALLOC(mean_feature_class_vector, num_classes);

  /* reserve a character pointer ( char *) for each class name */
  ALLOC( class_name, num_classes ); 
  
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
@NAME       : fcm_count_classes_and_set_names
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: sets num_classes, and class_names
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb. 16 1996 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void fcm_count_classes_and_set_names(void)
{
 
  int i, tmp_idx = 0, max_idx = 0, src_idx = 0, dst_idx = 0;

  /* make sure class names are provided to be matched with training
     volumes. They are parsed and checked to have a least the same
     count as the number of training volumes supplied. Parsing is
     done similar to the above method */
  
  if ( !classname_buffer) {
    
    fprintf( stderr, "Please specify class names in quotes, with -class_names \n");
    exit(EXIT_FAILURE);
  }

  else {

    VIO_ALLOC2D(class_name, num_train_vols, 10 );

    if ( debug > 4 )
      fprintf( stdout, "classname buffer = %s\n", classname_buffer);
    
    for_less( i, 0, num_train_vols ) {

      /* while buf is not a ',' or null, advance pointer and copy character */
      while ( classname_buffer[src_idx] != ',' &&
	      classname_buffer[src_idx] != '\0' )
	class_name[i][dst_idx++] = classname_buffer[src_idx++];

      /* as soon as you encounter a ',' or null, end target string */
      class_name[i][dst_idx] = '\0';
	
      /* keep track of the highest class index */
      tmp_idx = atoi(class_name[i]);

      if ( tmp_idx > max_idx ) 
	max_idx = tmp_idx;

      if ( debug > 4 )
	fprintf( stdout, "class_name[%d] = %s\n", i, class_name[i]);

      /* reset target string pointer */     
      dst_idx = 0;

      /* if you reach the end of the string, and the number of class names 
	 is less than the number of training volumes die, else advance */      
      if (classname_buffer[src_idx] == '\0' && i != ( num_train_vols - 1 )) {

	fprintf( stderr, "The number of class names and training vol mismatch\n");
	exit(EXIT_FAILURE);
      }
      else
	src_idx++;
      
    } /* for_less */
    
  } /* else */

  num_classes = num_train_vols;

  max_class_index = max_idx;

  /* reserve some space */
  ALLOC( create_fuzzy_volume, num_classes );
  
  /* for each class, set fuzzy volume creation flag to true, since fcm  */
  for_less( i, 0, num_classes ) 
    create_fuzzy_volume[i] = 1;

}



/* ----------------------------- MNI Header -----------------------------------
@NAME       : fcm_classify_sample
@INPUT      : 
@OUTPUT     : 
@RETURNS    : sample class
@DESCRIPTION: given a feature vector and its size, return a class
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 16,, 1996 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void fcm_classify_sample(int *class_num, double *class_prob, int *class_labels)
{

  int           i, j;                  /* counters */
  VIO_Real          temp;                  /* temporary variable */
  VIO_Real          max_prob = 0.0;

  /* initialize delta vector */
  for_less ( i, 0, num_classes ) 
    euclidian_vector[i] = 0.0;
  
  /* calculate Euclidian distance */
  for_less ( i, 0,  num_classes )
    for_less ( j, 0, num_features) {

      temp = feature_vector[j] - mean_feature_matrix[i][j]; /* X - Z */
      euclidian_vector[i] += temp * temp;                   /* Sig ( X - Z ) ^2 */
    }


  /* calculate denominator of Sig[ ( || Xk -Vi || / || Xk - Vj || )^(2(m-1))]
     which is || Xk - Vj ||, here euclidian_vector[i] =  || Xk - Vj || 
     for all classes Vj */
  for_less( i, 0, num_classes) {

    temp = 0.0;

    for_less( j, 0, num_classes) 
      /* Sig[ ( || Xk -Vi || / || Xk - Vj || ) */
      temp += euclidian_vector[i] / euclidian_vector[j];

    /* uik = { Sig[ ( || Xk -Vi || / || Xk - Vj || )^(2(m-1))] }^-1, Bezdek*/
    /* here also sqrt of euc_vec is unnecessary by changing power to (1/(m-1)) */
    class_prob[i] = 1 / pow( temp, (1/(m-1)));

    /* keep track of max probability value to assign the label to class_num */
    if ( class_prob[i] > max_prob ) {

      max_prob = class_prob[i];
      *class_num = i;
    }
    
    class_labels[i] = atoi(class_name[i]) ;
 
  } /* for_less (i,0,num_classes ) */


} /* fcm_classify_sample */ 


/* ----------------------------- MNI Header -----------------------------------
@NAME       : fcm_save_training
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
void fcm_save_training(char *save_train_filename)
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


