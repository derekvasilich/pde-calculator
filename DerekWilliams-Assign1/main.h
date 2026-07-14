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
#define NUM_SCHEMES		6

#define SEPERATOR			"---------------------------"
#define FORWARD_SPACE_NAME	"Forward-Time Forward-Space"
#define BACKWARD_SPACE_NAME	"Forward-Time Backward-Space"
#define CENTRAL_SPACE_NAME	"Forward-Time Central-Space"
#define	LAX_FREDRICHS_NAME	"Lax Fredrichs"
#define LEAPFROG_NAME		"Leapfrog"
#define EQUILIBRIUM_NAME	"Equilibrium"

#define H_10			10
#define H_20			20
#define H_40			40

#define H_10_NAME		"H = 1/10"
#define H_20_NAME		"H = 1/20"
#define H_40_NAME		"H = 1/40"

#define LAMBDA_8		100
#define LAMBDA_16		101

#define LAMBDA_8_NAME	"Lambda = 0.8"
#define LAMBDA_16_NAME	"Lambda = 1.6"

#define DEFAULT_SCHEME	BACKWARD_SPACE

#define BOUNDS_LEFT		-1
#define BOUNDS_RIGHT	3

#define TIME_START		0
#define TIME_END		2.4

#define OPENGL			1