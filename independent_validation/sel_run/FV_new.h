#ifndef _FV_new_h_
#define _FV_new_h_
//​
#include "TVector3.h"
//​
// Created by C Thorpe
// // Fiducial volume filter
// // Use to produce CC inclusive fiducial volume only samples
// ​
// // Fiducial volume definition
// ​
const double FVxmin =  10.0;
const double FVxmax = 246.0;
const double FVymin =-101.0; 
const double FVymax = 101.0; 
const double FVzmin =  10.0; 
const double FVzmax = 986.0; 
// ​

//for containment
//
const double FVxmin_1 = 2.0;            //12
const double FVxmax_1 = 256.35 - 2.0;  //244.35
const double FVymin_1 = -115.53 + 2.0;  //-80.53
const double FVymax_1 = 117.47 - 10.0;   //82
const double FVzmin_1 = 0.1 + 2.0;      //25.1
const double FVzmax_1 = 1036.9 - 2.0;   //951.9
//
//


 // Dead region to be cut
 const double deadzmin = 675.1;
 const double deadzmax = 775.1;
// ​
// ​
 bool inFV(TVector3 pos){
// ​
 if(pos.X() > FVxmax || pos.X() < FVxmin) return false;
// ​
 if(pos.Y() > FVymax || pos.Y() < FVymin) return false;
// ​
 if(pos.Z() > FVzmax || pos.Z() < FVzmin) return false;
// ​
 if(pos.Z() < deadzmax && pos.Z() > deadzmin) return false;
// ​
    return true;
       }

bool isContained(TVector3 endVec) {
	double padding = 0.0; // cm

	if (endVec.X() > (FVxmax_1 - padding) || endVec.X() < (FVxmin_1 + padding)) return false;
	if (endVec.Y() > (FVymax_1 - padding) || endVec.Y() < (FVymin_1 + padding)) return false;
	if (endVec.Z() > (FVzmax_1 - padding) || endVec.Z() < (FVzmin_1 + padding)) return false;

	return true;
}
//         ​
#endif
