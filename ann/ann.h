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
$RCSfile: ann.h,v $
$Revision: 1.2 $
$Author: bert $
$Date: 2005-02-11 20:19:31 $
$State: Exp $
--------------------------------------------------------------------------*/
#ifndef _ANN_H
#define _ANN_H

#include <iostream>		/* (bert) changed from iostream */
using namespace std;		/* (bert) added */

void ann_init_training(char *filename);
void ann_load_training(char *filename);
void ann_save_training(char *filename);
void ann_train_samples(void);
void ann_classify_sample(int *class_num, double *class_prob = 0, int *class_labels = 0);
ostream& ann_save_train_param(ostream&);
istream& ann_load_train_param(istream&);
void     ann_error_monitor(unsigned, double);

#endif
