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
$RCSfile: knn.cc,v $
$Revision: 1.1 $
$Author: jason $
$Date: 2002-03-20 22:16:35 $
$State: Exp $
--------------------------------------------------------------------------*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : knn.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: k-nearest neigbour classifier
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : June 10, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : $Log: knn.cc,v $
@MODIFIED   : Revision 1.1  2002-03-20 22:16:35  jason
@MODIFIED   : Initial revision
@MODIFIED   :
@MODIFIED   : Revision 1.1.1.1  1997/02/11 00:06:48  alex
@MODIFIED   : Sources for classify, copied from Vasken Kollokian
@MODIFIED   :
 * Revision 1.9  1996/03/14  06:30:30  vasco
 * added static on some multiply defined variables.
 *
 * Revision 1.8  1996/03/09  18:21:56  vasco
 * added limits.h in the header area.
 *
 * Revision 1.7  1996/03/06  00:27:54  vasco
 * fixed the euclidian distance sqrt problem to avoid unneccesary calculations
 *
 * Revision 1.6  1995/12/07  15:54:33  vasco
 * Dramatic improvement and bug fix over previous version.
 * This implementation follows the Keller paper to the letter,
 * also implemented fuzzy-knn algorithm, and hand tested it.
 *
 * Revision 1.5  1995/09/14  03:31:08  vasco
 * removed test_classifier
 *
 * Revision 1.4  1995/08/28  03:25:40  vasco
 * *** empty log message ***
 *
 * Revision 1.3  1995/07/12  05:20:56  vasco
 * added header info
 *
---------------------------------------------------------------------------- */
extern "C" {
#include <volume_io.h>
#include <limits.h>
}
#include "../class_globals.h"


/* locally define structures */
struct record_str  {
  int   index;
  VIO_Real  value;
};

typedef struct record_str record;

/* locally define prototypes */
int  compare_ascending(const void *a, const void *b);
  

/* locally defined global variables */
int    knn;                  /* this is the k value of knn classifier */
record *euclidian_vector;    /* euclidian distance vector between unknown sample
				and each training sample */
int    *knn_class_vector;    /* class vector to find majority in knn */
VIO_Real   *knn_tie_vector;      /* class vector to find tie totals */
static VIO_Real   m;             /* proxemity neighbourhood index */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : knn_init_training
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Jun 9, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void knn_init_training(char *param_filename)
{

  FILE *knn_par_file;

  /* check to see if the filename is there */
  if ( param_filename && !file_exists(param_filename)  ) {

    (void) fprintf(stderr,"File `%s' doesn't exist !\n ", param_filename);
    exit(EXIT_FAILURE);

  }

  if ( !param_filename ) {
    /* take the default here, in case a param file is not given */
    knn = 5;
    m = 2.0;

  }
  else {
    
    if (verbose) 
      fprintf(stdout, "Loading the parameter file %s\n", param_filename);      
    
    /* open the parameter file, and read the values  */
    knn_par_file = fopen(param_filename, "r");
    
    if ( knn_par_file == NULL) {
      fprintf(stderr, "Cannot open %s\n", param_filename);
      exit(EXIT_FAILURE);
    }
    
    /* scan for the k nearest neighbour number */
    fscanf( knn_par_file, "knn=%d\n", &knn);
 
    /* scan for the neighbourhood proxemity number */
    fscanf( knn_par_file, "m=%f\n", &m);
     
    fclose(knn_par_file);

    if ( m <= 1.0 ) {
 
     fprintf(stderr, "Invalid m value : %f\n", m);
      exit(EXIT_FAILURE);
    }


  }

  if ( debug > 2) {

    fprintf(stdout, "knn = %d\n", knn);      
    fprintf(stdout, "m   = %f\n", m);      
  }

  /* for training */
  ALLOC(euclidian_vector, knn);   /* allocate memory for vector */

  /* for classification */
  ALLOC(knn_class_vector, num_classes);   /* allocate memory for knn class vector*/

  /* for calculating totals for resolving ties */
  ALLOC(knn_tie_vector, num_classes);     /* allocate memory for tie vector*/

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : knn_train_samples
@INPUT      : 
@OUTPUT     : 
@RETURNS    : ?
@DESCRIPTION: takes a feature matrix and determines its k nearest neihbours
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 29, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void knn_train_samples(void)
{
  int        i, j;                          /* counters */
  VIO_Real       temp, dist;                    /* temporary variable */

  /* calculate the euclidian distance between unknown sample and
     the training samples - and square function is not applied
     in the euclidian distance calculation, since it is a relative
     term and minimizes computation */

  /* if debugging is being done on knn table, clear it out
     to improve readability */
  if ( debug > 7 ) 
    for_less( j, 0, knn ) {
      euclidian_vector[j].value = 0.0;
      euclidian_vector[j].index = 255;
    }
  

  for_less( i, 0, num_samples ) {

  
    /* initialize the distance measure */
    dist = 0.0;

    /* calculate the dist to cluster centroid */
    for_less( j, 0, num_features ) {

      temp = feature_vector[j] - feature_matrix[i][j];    /* X - Z */
      dist += (temp * temp);                              /* Sig ( X - Z) ^ 2 */ 

    }


    if ( i < knn ) {

      /* if the first knn spots are free, fill them */
      euclidian_vector[i].value = dist;
      euclidian_vector[i].index = class_column[i];
    }
    else {
      
      /* sort the euclidian_vector of k - nearest neigbours */
      qsort(euclidian_vector, knn, sizeof(record), compare_ascending );  

      
      /* see if last (most distant) neighbour can be kicked out */
      if ( dist < euclidian_vector[knn-1].value ) {
	
	euclidian_vector[knn-1].value = dist;
	euclidian_vector[knn-1].index = class_column[i];
      }


    } /* else */

    if ( debug > 7 ) {

      fprintf(stdout, "Knn table for sample %d, Edist = %7.2f\n", i, dist);

      for_less( j, 0, knn ) {

	fprintf(stdout, "EV[%d].value = %7.2f\t\t", j, euclidian_vector[j].value);
	fprintf(stdout, "EV[%d].index = %d\n", j, euclidian_vector[j].index);
      }

      fprintf(stdout, "\n\n");      

    } /* if (debug > 7) */

  } /* for_less( i, 0, num_samples ) */
    
} /* knn_train_samples */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : knn_classify_sample
@INPUT      : 
@OUTPUT     : 
@RETURNS    : sample class
@DESCRIPTION: given a feature vector, return a class
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : May 29, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void knn_classify_sample(int *class_num, double *class_prob, int *class_labels)
{
  int                i, j, k;
  int                max_votes;           /* counter to keep max votes */
  VIO_Real               minimum;           /* counter to keep min tie total */
  int                tie = FALSE;

  
  /* initialize the knn_class_vector and tie vector */
  for_less( i, 0, num_classes ) {

    knn_class_vector[i] = 0;
    knn_tie_vector[i] = -1.0; /* denote a free spot */
  }

  /* take top knn elements of euclidian_vector and see what is the majority class */
  for_less( i, 0, knn )
    knn_class_vector[ euclidian_vector[i].index ]++;

  /* set max votes to zero to start from scratch */
  max_votes = 0;
  minimum =  10000000; //DBL_MAX;
  *class_num = 0;

  /* get the max number of votes */
  for_less( i, 0, num_classes )
    if ( max_votes <  knn_class_vector[i] ) {
      
      max_votes = knn_class_vector[i];
      *class_num = i;
    }

  if ( debug > 8 ) {

    for_less( i, 0, num_classes)
      fprintf(stdout, "knn_class_vector[%d] = %d\n", i, knn_class_vector[i]);

    fprintf(stdout, "max_vote = %d, index = %d\n", max_votes, *class_num);
  }

  /* check for tie votes in the highest votex indicated by max_votes, 
     if there are any, calculate tie totals */
  for_less( i, 0, num_classes ) {

    for_less( j, i+1, num_classes ) {

      /* max sure only ties in max_votes are considered */
      if ( knn_class_vector[i] == knn_class_vector[j] &&
	   knn_class_vector[i] == max_votes ) { 
	
	if ( debug > 8 ) 
	  fprintf(stdout, "tie occured at max_vote = %d\n", max_votes);
	
	/* a tie has occured */
	tie = TRUE;

	/* if tie total for first item is not calculated, then do so */
	if ( knn_tie_vector[i] < 0.0 ) {

	  /* reset the tie slot */
	  knn_tie_vector[i] = 0.0;
	  for_less( k, 0, knn ) 
	    if ( euclidian_vector[k].index == i )
	      knn_tie_vector[i] +=  euclidian_vector[k].value;
	}

	/* if tie total for second item is not calculated, then do so */
	if ( knn_tie_vector[j] < 0.0 ) {

	  /* reset the tie slot */
	  knn_tie_vector[j] = 0.0;
	  for_less( k, 0, knn ) 
	    if ( euclidian_vector[k].index == j )
	      knn_tie_vector[j] +=  euclidian_vector[k].value;
	}

      } /* ( knn_class_vector[i] == knn_class_vector[j] ) */

    } /* for_less(j...) */

  } /* for_less( i, 0, num_classes ) */

  if ( tie ) {
     
    for_less( i, 0, num_classes ) {

      if ( knn_tie_vector[i] <= minimum && knn_tie_vector[i] >= 0.0 ) { 

	/* the positive number check here denotes a calculated slot */

	/* last minimum found */
      	minimum = knn_tie_vector[i];
	*class_num = i;
      }

    }

    if ( debug > 8 ) {
      
      for_less( i, 0, num_classes)
	fprintf(stdout, "knn_tie_vector[%d] = %7.2f\n", i, knn_tie_vector[i]);
      
      fprintf(stdout, "minimum = %7.2f,  index = %d\n", minimum, *class_num);
    }

  } /* if (tie) */

  /*  FUZZY STARTS HERE, FUZZY STARTS HERE, FUZZY STARTS HERE, */
  /*  FUZZY STARTS HERE, FUZZY STARTS HERE, FUZZY STARTS HERE, */
  /*  FUZZY STARTS HERE, FUZZY STARTS HERE, FUZZY STARTS HERE, */
  /*  FUZZY STARTS HERE, FUZZY STARTS HERE, FUZZY STARTS HERE, */


  if ( class_prob ) {

    VIO_Real  sigma_denom = 0.0;     /* dinominator of sig expression */
    VIO_Real  sigma_numer = 0.0;     /* numerator of sig expression */

    for_less( i, 0, knn) {

      if ( euclidian_vector[i].value != 0 ) {
	
	/* here the sqrt(euclidian_vector[i].value) is necessary to restore
	   the euclidian vector to || Xk - Zi || form, but in order
	   to avoid extra computations, the formula that could be
	   calculated as : 
	   
	   pow( sqrt(euclidian_vector[i].value), (2/(m-1) ) )  is replaced by

	   pow( euclidian_vector[i].value, (1/(m-1) ) ) since

	   sqrt(x)^(2/(m-1)) === x^(1/(m-1))

	 */


	sigma_denom += 1 / ( pow( euclidian_vector[i].value, (1/(m-1) ) ));
      } /* if ( euclidian_vector[i].value != 0 ) */

    } /* for_less( i, 0, knn) */

    
    for_less( i, 0, num_classes) {

      class_prob[i] = 0.0;
      sigma_numer = 0.0;

      for_less( j, 0, knn) {
	
	if ( euclidian_vector[j].index == i && euclidian_vector[j].value != 0 ) {

	  /* restore it back to (|| Xk - Vi ||) form by the eariler method */
	  sigma_numer += 1 / ( pow( euclidian_vector[j].value, (1/(m-1) ) ));
	    
	} /* if ( euclidian_vecto ... ) */

      } /* for_less( j, 0, knn) */
 
      
      class_prob[i] = sigma_numer / sigma_denom;

      if ( class_labels ) 
	  class_labels[i] = euclidian_vector[i].index;

    } /* for_less( i, 0, num_classes) */

  } /* if ( class_prob ) */

} /* knn_classify_sample(...) */


/* ----------------------------- MNI Header -----------------------------------
@NAME       : knn_load_training
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Jun 9, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void knn_load_training(char *load_train_filename)
{

  fprintf(stderr, "k-nearest neigbout classifier has no -load_train option...\n");
  exit(EXIT_FAILURE);

}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : knn_init_training
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Jun 9, 1995 (Vasco KOLLOKIAN)
@MODIFIED   : 
---------------------------------------------------------------------------- */
void knn_save_training(char *save_train_filename)
{

  fprintf(stderr, "k-nearest neigbout classifier has no -save_train option...\n");
  exit(EXIT_FAILURE);

}


int  compare_ascending(const void *a, const void *b)
{

  VIO_Real left, right;

  left  = ((record *)a)->value;
  right = ((record *)b)->value;

  if ( left > right ) 
    return 1;
  else if ( left < right ) 
    return -1;
  else
    return 0;

}

