/*
 *  main.h
 *  Assignment One
 *
 *  Created by Derek Williams on 10-01-30.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#define FORWARD_SPACE	0
#define BACKWARD_SPACE	1
#define CENTRAL_SPACE	2
#define	LAX_FREDRICHS	3
#define LEAPFROG		4
#define EQUILIBRIUM	    5
#define LAX_WENDROFF    6
#define NUM_SCHEMES		7

#define SEPERATOR			"---------------------------"

#define FORWARD_SPACE_NAME	"Forward-Time Forward-Space"
#define BACKWARD_SPACE_NAME	"Forward-Time Backward-Space"
#define CENTRAL_SPACE_NAME	"Forward-Time Central-Space"
#define	LAX_FREDRICHS_NAME	"Lax Fredrichs"
#define LEAPFROG_NAME		"Leapfrog"
#define EQUILIBRIUM_NAME	"Equilibrium"
#define LAX_WENDROFF_NAME	"Lax Wendroff"

#define H_10			10
#define H_20			20
#define H_40			40
#define H_80			80

#define H_10_NAME		"1/10"
#define H_20_NAME		"1/20"
#define H_40_NAME		"1/40"
#define H_80_NAME		"1/80"

#define LAMBDA_8		100
#define LAMBDA_16		101

#define LAMBDA_8_NAME	"0.8"
#define LAMBDA_16_NAME	"1.6"

#define A_CONST			1000
#define A_SPACE			1001

#define A_CONST_NAME	"1"
#define A_SPACE_NAME	"1 + 1/4(3 - x)(1 + x)"

#define INIT_ONE		2000
#define INIT_TWO		2001
#define INIT_THREE		2002
#define INIT_FOUR		2003

#define INIT_ONE_NAME	"cos^2(PI x), |x| <= 0.5"
#define INIT_TWO_NAME	"1 - |x|, |x| <= 1"
#define INIT_THREE_NAME	"sin(2 PI x)"
#define INIT_FOUR_NAME	"sin(PI x)"

#define MAX_VM_N		3000

#define F_ZERO			4000
#define F_ONE			4001

#define F_ZERO_NAME		"F(t, x) = 0"
#define F_ONE_NAME		"F(t, x) = sin(PI (x-t))"

#define DEFAULT_SCHEME	BACKWARD_SPACE

#define BOUNDS_LEFT		0
#define BOUNDS_RIGHT	1

#define TIME_START		0
#define TIME_END		100