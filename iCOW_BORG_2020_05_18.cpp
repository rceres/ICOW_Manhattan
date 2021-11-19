/* Copyright 2012-2014 The Pennsylvania State University
 *
 * This software was written by David Hadka and others.
 * 
 * The use, modification and distribution of this software is governed by the
 * The Pennsylvania State University Research and Educational Use License.
 * You should have received a copy of this license along with this program.
 * If not, contact <dmh309@psu.edu>.
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
//#include <chrono>
//#include <random>
#include<algorithm>
#include "borgms.h"
#include "borg.h"
//#include "borg.c"
#include "ICOW.h"
#include <string>
#include <cstring>



using namespace std;
int nvars;
int nobjs;
int ninfo;
int nconsts;


static uint GetUint()
{
  m_z = 36969 * (m_z & 65535) + (m_z >> 16);
  m_w = 18000 * (m_w & 65535) + (m_w >> 16);
  return (m_z << 16) + m_w;
}

double GetUniform()
{
  // adapted from code by John D. Cook 
  // 0 <= u < 2^32
  uint u = GetUint();
  // The magic number below is 1/(2^32 + 2).
  // The result is strictly between 0 and 1.
  return (u + 1.0) * 2.328306435454494e-10;
}

double GetGEV(double mu, double sigma, double xi)
{
  // adapted from code by John D. Cook 
  double C=GetUniform();
  return(
    -(pow(-(log(C)),-xi)*(-mu*xi*pow(-log(C),xi) + sigma*(pow(-log(C),xi)) - sigma))/xi
  );
  //return(-((-sigma-mu*xi*pow((-log(C)),xi)+sigma*pow((-log(C)),xi)*pow((-log(C)),(-xi))) /xi) );
}

void generateSurges(double *surges, double loc, double scale, double shape, double dLoc, double dScale, double dShape)
{double s;
  for (int i=0;i<maxStatesToEvaluate;i++)
  {
    for (int j=0;j<lengthSurgeSequences;j++)
    {
      //surges[i*lengthSurgeSequences+j]=GetGEV((loc+dLoc*j),(scale+dScale*j),(shape+dShape*j));
      s=GetGEV((loc+dLoc*j),(scale+dScale*j),(shape+dShape*j));
      /*if (s>(maxBaseSurge+j*surgeMaxRateIncrease/100)) 
      {
        s=maxBaseSurge+j*surgeMaxRateIncrease/100;
      }*/
      surges[i*lengthSurgeSequences+j]=s;
    }
  }
}


double CalculateDikeCost(double hd,double cd,double S,double W,double sd,double wdt,double ich){
  // hd height of dike
  // cd cost of dike per unit volume
  // S slope of ground
  // W width of dike
  // sd slope of the dike sides
  // wdt width of the top of the dike
  // ich initial cost height
  double result;
  
  
  double ch;  /* Cost height is the dike height plus the equivalent height for startup costs. */
double ch2; /* Cost height squared, used a lot, so calculate it once */
double vd;  /* volume of dike */
ch=hd+ich;
ch2=pow(ch,2);
vd=W*ch*(wdt+ch/pow(sd,2))+
  pow((
      -pow(ch,4)*pow((ch+1/sd),2)/pow(sd,2)-
        2*(pow(ch,5)*(ch+1/sd))/pow(S,4)-
        4*pow(ch,6)/(pow(sd,2)*pow(S,4))+
        4*pow(ch,4)*(2*dh*(ch+1/sd)-4*ch2+ch2/pow(sd,2))/(pow(sd,2)*pow(S,2))+
        2*pow(ch,3)*(ch+1/sd)/pow(S,2)
  ),1/2)/6+
    W*(ch2/pow(S,2));
// volume of front of dike + volume of straight part of two sides of dike + volume tetrahedron part of sides */
result=vd*cd;
return result;
}


double CalculateWithdrawalCost(double * cityChar)
  // vi value of initial infrastructure
  // hw amount of height to withdraw to
  // h total height change of city
  // cw percent of value required to withdraw
{
  double result;
  if (cityChar[wh]==0) result=0;
  else result=cityChar[tcvi]*cityChar[wh]/(CEC-cityChar[wh])*WithdrawelCostFactor;
  return(result);
}

double CalculateResiliencyCost1(double * cityChar){
  //dike base height is lower than resiliency height, there is an unprotected nonResiliant zone
  double fractionResilient=(Basement+cityChar[rh]/2) / BH;
  /*       double rcf = resistanceExponentialFactor *
  (cityChar[rp] +
  pow( (pow( cityChar[rp], RF2) +
  1) ,
  RF2) -
  1);*/
  double fcR = resistanceExponentialFactor*std::max(0.0,(cityChar[rp]-resistanceExponentialThreshold))/(1.0-cityChar[rp]) +
    cityChar[rp]*resistanceLinearFactor;
  return ( cityChar[vz1] * fcR * (fractionResilient));
}

//dike base height is lower than resiliency height, there is NOT an unprotected nonResiliant zone
// cases 2 and 6
double CalculateResiliencyCost2(double * cityChar){
  /*double rcf = resistanceExponentialFactor *
  (cityChar[rp] +
  pow( (pow( cityChar[rp], RF2) +
  1) ,
  RF2) -
  1);*/
  double fcR = resistanceExponentialFactor*std::max(0.0,(cityChar[rp]-resistanceExponentialThreshold))/(1.0-cityChar[rp]) +
    cityChar[rp]*resistanceLinearFactor;
  return( cityChar[vz1] * fcR * (Basement+cityChar[rh]-cityChar[dbh]/2) / BH);
}

double CalculateCostOfInfrastructureLostFromWithdrawal(double * cityChar)
{
  
  //        return(cityChar[tcvi]*cityChar[fw]*WithdrawelPercentLost);
  return(cityChar[tcvi]*cityChar[fw]*WithdrawelPercentLost);
}


double CalculateFinalValueOfInfrastructure(double vi,double vil)
  // vi intial value of all infrastructure
  // vil value of infrastructure leaving
{
  return(vi-vil);
}


double culateTotalCostAbatement(double cd,double cw,double cvlw,double cr)
  //   cd cost of dike
  //   cw cost of withdrawal
  //   cvlw cost of infrastructure lost due to withdrawal
  //   cr cost of resilancy
{
  double result;
  result=cd+cw+cvlw+cr;
  return(result);
}

void CharacterizeCity (double W,double B,double R,double P, double D, double * cityChar) {
  
  // check for base value and strategies above min heights
  if (W==baseValue||W<WminHeight) {cityChar[wh]=0.0;} else {cityChar[wh]=W;}
  if (R==baseValue || R<minHeight) 
  {
    cityChar[rh]=0.0;
    cityChar[rp]=0.0;
  }
  else
  {
    cityChar[rh]=R;
    cityChar[rp]=P;
  }
  if (D==baseValue || D<minHeight) {cityChar[dh]=0.0;} else {cityChar[dh]=D;}
  if (B==baseValue || B<minHeight) {cityChar[dbh]=0.0;} else {cityChar[dbh]=B;}
  
  // calculate the damage that results according to the resistance percent
  cityChar[dtr]=std::max(1-cityChar[rp],0.0);
  
  //check to see if the distance between top of withdrawal and dike base is too small
  if ((cityChar[dh]>=minHeight)&&(cityChar[dbh]<minHeight)&&(cityChar[rh]>=minHeight))
  {
    cityChar[dbh]=0;
    cityChar[rh]=0;
  }         
  int c=100; // which case?
  
  // (cityChar[dh]>0)
  if (cityChar[dh]>0) {
    
    // (cityChar[dh]>0) && (cityChar[dbh]>0)
    if (cityChar[dbh]>0) {
      
      // (cityChar[dh]>0) && (cityChar[dbh]>0) && (cityChar[rh]>0)
      if (cityChar[rh]>0) {
        
        // (cityChar[dh]>0) && (cityChar[dbh]>0) && (cityChar[rh]>0) && (cityChar[rh]<cityChar[dbh]) c=1
        if (cityChar[rh]<cityChar[dbh]) {
          c=1;
        }
        
        // (cityChar[dh]>0) && (cityChar[dbh]>0) && (cityChar[rh]>0) && (cityChar[rh]>=cityChar[dbh]) c=2
        else {
          c=2;
        }
        
      }
      // (cityChar[dh]>0) && (cityChar[dbh]>0) && (cityChar[rh]=0) c=3
      else {
        c=3;
      }
      
    }
    
    // (cityChar[dh]>0) && (cityChar[dbh]=0)
    else {
      c=4;             //there is a dike, there is no setback, there is or is not resilency
    }
  }
  
  // else (cityChar[dh]=0)
  else {
    if (cityChar[dbh]>0) {  // (cityChar[dh]=0) && (cityChar[dbh]>0)
      
      // if (cityChar[dh]=0) && (cityChar[dbh]>0) && (cityChar[rh]>0)
      if (cityChar[rh]>0) {
        
        // (cityChar[dh]=0) && (cityChar[dbh]>0) && (cityChar[rh]>0) && (cityChar[rh]<cityChar[dbh])
        if (cityChar[rh]<cityChar[dbh]) {
          c=5;  // no dike, but there is set back, and there is resiliency, resiliancy is lower than set back
        }
        // (cityChar[dh]=0) && (cityChar[dbh]>0) && (cityChar[rh]>0) && (cityChar[rh]>=cityChar[dbh])
        else {
          c=6;  // no dike, but there is set back, and there is resiliency, resiliancy is equal or higher than set back
        }
      }
      // (cityChar[dh]=0) && (cityChar[dbh]>0) && (cityChar[rh]=0)
      else {
        c=7;
      }
    }
    
    // (cityChar[dh]=0) && (cityChar[dbh]=0)
    else {
      
      // (cityChar[dh]=0) && (cityChar[dbh]=0) && (cityChar[rh]>0)
      if (cityChar[rh]>0) {
        c=8;
      }
      // (cityChar[dh]=0) && (cityChar[dbh]=0) && (cityChar[rh]=0)
      else  {
        c=9;
      }
    }
  }
  
  // implications of withdrawel are calculeted first since they will impact the rest of the calculations
  cityChar[tcvi]=TotalCityValueInitial;
  cityChar[wc]=CalculateWithdrawalCost(cityChar); //
  cityChar[fw]=cityChar[wh]/CEC; // Calculate Fraction Withdrawn
  cityChar[ilfw]=CalculateCostOfInfrastructureLostFromWithdrawal(cityChar);
  cityChar[tcvaw]=cityChar[tcvi]-cityChar[wc]; // calculate total city value after withdrawel
  cityChar[caseNum]=c;
  switch ( c ) {
  case 0:       // not valid
    break;
  case 1:        // (cityChar[dh]>0) && (cityChar[dbh]>0) && (cityChar[rh]>0) && (cityChar[rh]<cityChar[dbh])
    // city has all four zones
    
    cityChar[dc]=CalculateDikeCost(cityChar[dh],UnitCostPerVolumeDike,CitySlope,CityWidth,SlopeDike,WidthDikeTop,DikeStartingCostPoint);
    // and the dike is setback (cityChar[dbh]>0) and there is resiliency
    cityChar[vz1] = cityChar[tcvaw]*DikeUnprotectedValuationRatio*cityChar[rh]/(CEC-cityChar[wh]); // calculate value zone 1
    cityChar[vz2] = cityChar[tcvaw]*DikeUnprotectedValuationRatio*(cityChar[dbh]-cityChar[rh])/(CEC-cityChar[wh]); // calculate value zone 2
    cityChar[vz3] = cityChar[tcvaw]*ProtectedValueRatio*cityChar[dh]/(CEC-cityChar[wh]); // calculate value zone 3
    cityChar[vz4] = cityChar[tcvaw]*(CEC-cityChar[wh]-cityChar[dbh]-cityChar[dh])/(CEC-cityChar[wh]); // calculate value zone 4
    cityChar[fcv] = cityChar[vz1]+cityChar[vz2]+cityChar[vz3]+cityChar[vz4];
    cityChar[tz1] = cityChar[wh]+cityChar[rh];
    cityChar[tz2] = cityChar[wh]+cityChar[dbh];
    cityChar[tz3] = cityChar[wh]+cityChar[dbh]+cityChar[dh];
    cityChar[tz4] = CEC;
    cityChar[rc]  = CalculateResiliencyCost1(cityChar);
    cityChar[tic] = cityChar[wc]+cityChar[dc]+cityChar[rc];
    cityChar[tc]  = cityChar[tic]+cityChar[fcv]-cityChar[tcvi];
    break;
  case 2:         // (cityChar[dh]>0) && (cityChar[dbh]>0) && (cityChar[rh]>0) && (cityChar[rh]>=cityChar[dbh])
    cityChar[dc] = CalculateDikeCost(cityChar[dh],UnitCostPerVolumeDike,CitySlope,CityWidth,SlopeDike,WidthDikeTop,DikeStartingCostPoint);
    cityChar[vz1] = cityChar[tcvaw]*DikeUnprotectedValuationRatio*cityChar[dbh]/(CEC-cityChar[wh]); // calculate value zone 1
    cityChar[vz2] = 0; // there is no unprotected zone in front of the dike
    cityChar[vz3] = cityChar[tcvaw]*ProtectedValueRatio*cityChar[dh]/(CEC-cityChar[wh]); // calculate value zone 3
    cityChar[vz4] = cityChar[tcvaw]*(CEC-cityChar[wh]-cityChar[dbh]-cityChar[dh])/(CEC-cityChar[wh]); // calculate value zone 4
    cityChar[fcv] = cityChar[vz2]+cityChar[vz3]+cityChar[vz4];
    cityChar[tz1] = cityChar[wh]+cityChar[dbh];
    cityChar[tz2] = cityChar[wh]+cityChar[dbh];
    cityChar[tz3] = cityChar[wh]+cityChar[dbh]+cityChar[dh];
    cityChar[tz4] = CEC;
    cityChar[rc]  = CalculateResiliencyCost2(cityChar);
    cityChar[tic] = cityChar[wc]+cityChar[dc]+cityChar[rc];
    cityChar[tc]  = cityChar[tic]+cityChar[fcv]-cityChar[tcvi];
    break;
  case 3:        // (cityChar[dh]>0) && (cityChar[dbh]>0) && (cityChar[rh]=0)
    //the dike is not at the seawall and there is no resilancy
    cityChar[dc]  = CalculateDikeCost(cityChar[dh],UnitCostPerVolumeDike,CitySlope,CityWidth,SlopeDike,WidthDikeTop,DikeStartingCostPoint);
    cityChar[vz1] = 0; // not needed, we defined it this way
    cityChar[vz2] = cityChar[tcvaw]*DikeUnprotectedValuationRatio*cityChar[dbh]/(CEC-cityChar[wh]); // calculate value zone 2
    cityChar[vz3] = cityChar[tcvaw]*ProtectedValueRatio*cityChar[dh]/(CEC-cityChar[wh]); // calculate value zone 3
    cityChar[vz4] = cityChar[tcvaw]*(CEC-cityChar[wh]-cityChar[dbh]-cityChar[dh])/(CEC-cityChar[wh]); // calculate value zone 4
    cityChar[fcv] = cityChar[vz2]+cityChar[vz3]+cityChar[vz4];
    cityChar[tz1] = cityChar[wh];
    cityChar[tz2] = cityChar[wh]+cityChar[dbh];
    cityChar[tz3] = cityChar[wh]+cityChar[dbh]+cityChar[dh];
    cityChar[tz4] = CEC;
    cityChar[rc]  = 0; // there is no resiliency
    cityChar[tic] = cityChar[wc]+cityChar[dc]; // no resiliency cost
    cityChar[tc]  = cityChar[tic]+cityChar[fcv]-cityChar[tcvi];
    break;
  case 4:        // (cityChar[dh]>0) && (cityChar[dbh]=0)
    cityChar[dc]  = CalculateDikeCost(cityChar[dh],UnitCostPerVolumeDike,CitySlope,CityWidth,SlopeDike,WidthDikeTop,DikeStartingCostPoint);
    cityChar[vz1] = 0; // there is no protected zone in front of the dike
    cityChar[vz2] = 0; // there is no unprotected zone in front of the dike
    cityChar[vz3] = cityChar[tcvaw]*ProtectedValueRatio*cityChar[dh]/(CEC-cityChar[wh]); // calculate value zone 3
    cityChar[vz4] = cityChar[tcvaw]*(CEC-cityChar[wh]-cityChar[dh])/(CEC-cityChar[wh]); // calculate value zone 4
    cityChar[fcv] = cityChar[vz3]+cityChar[vz4];
    cityChar[tz1] = cityChar[wh];
    cityChar[tz2] = cityChar[wh];
    cityChar[tz3] = cityChar[wh]+cityChar[dh];
    cityChar[tz4] = CEC;
    cityChar[rc]  = 0; // width of the resiliency zone is zero
    cityChar[tic] = cityChar[wc]+cityChar[dc]; // no resiliency cost
    cityChar[tc]  = cityChar[tic]+cityChar[fcv]-cityChar[tcvi];
    break;
  case 5:        // (cityChar[dh]=0) && (cityChar[dbh]>0) && (cityChar[rh]>0) && (cityChar[rh]<cityChar[dbh])
    cityChar[dc]  = CalculateDikeCost(cityChar[dh],UnitCostPerVolumeDike,CitySlope,CityWidth,SlopeDike,WidthDikeTop,DikeStartingCostPoint);
    cityChar[vz1] = cityChar[tcvaw]*DikeUnprotectedValuationRatio*cityChar[rh]/(CEC-cityChar[wh]); // calculate value zone 1
    cityChar[vz2] = cityChar[tcvaw]*DikeUnprotectedValuationRatio*(cityChar[dbh]-cityChar[rh])/(CEC-cityChar[wh]); // calculate value zone 2
    cityChar[vz3] = 0; // dike height is 0
    cityChar[vz4] = cityChar[tcvaw]*(CEC-cityChar[wh]-cityChar[dbh])/(CEC-cityChar[wh]); // calculate value zone 4
    cityChar[fcv] = cityChar[vz1]+cityChar[vz2]+cityChar[vz4]; // , cityChar[dh]=0, so no zone 3
    cityChar[tz1] = cityChar[wh]+cityChar[rh];
    cityChar[tz2] = cityChar[wh]+cityChar[dbh];
    cityChar[tz3] = cityChar[tz2]; // dike height is 0
    cityChar[tz4] = CEC;
    cityChar[rc]  = CalculateResiliencyCost1(cityChar);
    cityChar[tic] = cityChar[wc]+cityChar[dc]+cityChar[rc];
    cityChar[tc]  = cityChar[tic]+cityChar[fcv]-cityChar[tcvi];
    break;
  case 6: // (cityChar[dh]=0) && (cityChar[dbh]>0) && (cityChar[rh]>0) && (cityChar[rh]>=cityChar[dbh])
    cityChar[dc]  = 0;
    cityChar[vz1] = cityChar[tcvaw]*DikeUnprotectedValuationRatio*cityChar[dbh]/(CEC-cityChar[wh]); // calculate value zone 1
    cityChar[vz2] = 0; // calculate value zone 2
    cityChar[vz3] = 0; // cityChar[dh]=0
    cityChar[vz4] = cityChar[tcvaw]*(CEC-cityChar[wh]-cityChar[dbh])/(CEC-cityChar[wh]); // calculate value zone 4
    cityChar[fcv] = cityChar[vz1]+cityChar[vz4];
    cityChar[tz1] = cityChar[wh]+cityChar[dbh];
    cityChar[tz2] = cityChar[tz1];
    cityChar[tz3] = cityChar[tz1];
    cityChar[tz4] = CEC;
    cityChar[rc]  = CalculateResiliencyCost2(cityChar);
    cityChar[tic] = cityChar[wc]+cityChar[dc]+cityChar[rc];
    cityChar[tc]  = cityChar[tic]+cityChar[fcv]-cityChar[tcvi];
    break;
  case 7: // (cityChar[dh]=0) && (cityChar[dbh]>0) && (cityChar[rh]=0)
    cityChar[dc]  = CalculateDikeCost(cityChar[dh],UnitCostPerVolumeDike,CitySlope,CityWidth,SlopeDike,WidthDikeTop,DikeStartingCostPoint);
    cityChar[vz1] = 0; // calculate value zone 1
    cityChar[vz2] = cityChar[tcvaw]*DikeUnprotectedValuationRatio*cityChar[dbh]/(CEC-cityChar[wh]); // calculate value zone 2
    cityChar[vz3] = 0; // calculate value zone 3
    cityChar[vz4] = cityChar[tcvaw]*(CEC-cityChar[wh]-cityChar[dbh])/(CEC-cityChar[wh]); // calculate value zone 4
    cityChar[fcv] = cityChar[vz2]+cityChar[vz4];
    cityChar[tz1] = cityChar[wh];
    cityChar[tz2] = cityChar[wh]+cityChar[dbh];
    cityChar[tz3] = cityChar[tz2];
    cityChar[tz4] = CEC;
    cityChar[rc]  = 0;
    cityChar[tic] = cityChar[wc]+cityChar[dc];
    cityChar[tc]  = cityChar[tic]+cityChar[fcv]-cityChar[tcvi];
    break;
  case 8:  // (cityChar[dh]=0) && (cityChar[dbh]=0) && cityChar[rh]>0
    cityChar[dc]=0;
    cityChar[vz1] = cityChar[tcvaw]*cityChar[rh]/(CEC-cityChar[wh]); // calculate value zone 1
    cityChar[vz2] = 0; // calculate value zone 2
    cityChar[vz3] = 0; // calculate value zone 3
    cityChar[vz4] = cityChar[tcvaw]*(CEC-cityChar[wh]-cityChar[rh])/(CEC-cityChar[wh]); // calculate value zone 4
    cityChar[fcv] = cityChar[vz1]+cityChar[vz4];
    cityChar[tz1] = cityChar[wh]+cityChar[rh];
    cityChar[tz2] = cityChar[tz1];
    cityChar[tz3] = cityChar[tz1];
    cityChar[tz4] = CEC;
    cityChar[rc]  = CalculateResiliencyCost1(cityChar);
    cityChar[tic] = cityChar[wc]+cityChar[rc];
    cityChar[tc]  = cityChar[tic]+cityChar[fcv]-cityChar[tcvi];
    break;
  case 9: // (cityChar[dh]=0) && (cityChar[dbh]=0) && cityChar[rh]=0
    cityChar[fcv] = cityChar[tcvaw];
    cityChar[dc]  = 0;
    cityChar[vz1] = 0; // calculate value zone 1
    cityChar[vz2] = 0; // calculate value zone 2
    cityChar[vz3] = 0; // calculate value zone 3
    cityChar[vz4] = cityChar[tcvaw]; // calculate value zone 4
    cityChar[fcv] = cityChar[vz4];
    cityChar[tz1] = cityChar[wh];
    cityChar[tz2] = cityChar[wh];
    cityChar[tz3] = cityChar[wh];
    cityChar[tz4] = CEC;
    cityChar[rc]  = 0;
    cityChar[tic] = cityChar[wc];
    cityChar[tc]  = cityChar[tcvi]-cityChar[fcv];
    break;
  }
  if (cityChar[tic]==0) {cityChar[tic]=1;}
}

double CalculateDamageResiliantUnprotectedZone1(double sl,double * cityChar)
  //dike base height is higher than resiliency height, there is an unprotected nonResiliant zone
{
  double washOver=sl-cityChar[wh];
  if (sl<=cityChar[tz1])
  {
    // surge is at or below resilient height
    return(
      cityChar[vz1]*damageFactor*cityChar[dtr]*
        washOver*(washOver/2+Basement)/
          (BH*cityChar[rh]));
  }
  else
  {
    //surge is higher than resilient height
    return(
      cityChar[vz1]*damageFactor*
        (cityChar[dtr]*(cityChar[rh]/2+Basement)+ // resilient zone is flooded to the top of the resilient height plus
        (washOver-cityChar[rh]))/  // resilient zone is flooded above the resilient height
        BH
    );
  }
}



double CalculateDamageResiliantUnprotectedZone2(double sl,double * cityChar)
  //used when cityChar[rh]>cityChar[dbh], there is no unprotected nonresiliant zone
{
  double washOver=sl-cityChar[wh];
  if (sl<cityChar[tz1])
  {
    // surge is below dike base height
    return(
      cityChar[vz1]*damageFactor*
        cityChar[dtr]*
        washOver*(washOver/2+Basement)/
          (BH*cityChar[dbh]));
  }
  else
  {
    //surge is higher than dike base height
    if (washOver<(cityChar[rh]))
      // surge is higher than the dike base, but not exceeding the cityChar[rh]
    {
      return(
        cityChar[vz1]*damageFactor*
          cityChar[dtr]*
          (Basement+washOver-cityChar[dbh]/2)/
            BH);
      //return(pow(cityChar[vz1]*cityChar[dtr]*cityChar[dbh]*cityChar[dbh]/2+sl-cityChar[tz1]/BH,damageratioexponent));
    }
    else
      // surge is higher than the dike base, and exceeds the cityChar[rh]
    {
      return(
        cityChar[vz1]*damageFactor*
          (cityChar[dtr]*
          (Basement+washOver-cityChar[dbh]/2)+
          washOver-cityChar[rh])/
            BH);
    }
  }
}


double CalculateDamageNonResiliantUnprotectedZone(double sl,double * cityChar)
{
  double washOver=sl-cityChar[tz1];
  if (sl<cityChar[tz2])
  {
    // surge is below dbh
    return(cityChar[vz2]*damageFactor*
           
           washOver*(washOver/2+Basement)/
             (BH*(cityChar[dbh]-cityChar[rh])));
  }
  else
  {
    // surge is higher than dbh
    return(
      cityChar[vz2]*damageFactor*
        (Basement+washOver+(cityChar[rh]-cityChar[dbh])/2)/
          BH);
  }
}

double CalculateDamageProtectedZone(double sl,double * cityChar)
  
{
  double washOver=sl-cityChar[tz2];
  if (sl>=cityChar[tz3]) return(cityChar[vz3]*FailedDikeDamageFactor*damageFactor*
      cityChar[dh]*(Basement+cityChar[dh]/2+
      (sl-cityChar[tz3]))/
        (BH*(cityChar[tz3]-cityChar[tz2])));
  else {
    double pf = std::max(pfBase, ((sl-cityChar[tz2])/cityChar[dh]-pfThreshold)/(1-pfThreshold));
    if((double)(rand() % 10000)/10000<(1.0-pf))return(
        cityChar[vz3]*intactDikeDamageFactor*damageFactor*
          washOver*(Basement+washOver/2)/(BH*(cityChar[tz3]-cityChar[tz2]))); // dike intact
    
    else return(cityChar[vz3]*FailedDikeDamageFactor*damageFactor*
                washOver*(Basement+washOver/2)/(BH*(cityChar[tz3]-cityChar[tz2]))); //dike failed
  }
}


double CalculateDamageAboveDikeProtectionZone(double sl,double * cityChar)
{
  double washOver=sl-cityChar[tz3];
  return(cityChar[vz4]*damageFactor*
         washOver*(washOver/2+Basement)/
           BH/(cityChar[tz4]-cityChar[tz3]));
}







//std::random_device rd;
using namespace std;
void CalculateDamageVector(double s,double * cityChar,double * damagevector)
  
{
  memset(damagevector, 0.0, sizeof(double)*8);
  
  int c=(int) cityChar[caseNum];
  double surge=0;
  if (s>Seawall) surge=s*runUpWave-Seawall;
  
  if(surge>cityChar[wh]) //otherwise there will be no damage
  {
    
    switch ( c ) {
    case 0:       // not valid
  {
    break;
  }
    case 1:       // (cityChar[dh]>0.0) && (cityChar[dbh]>0.0) && (cityChar[rh]>0.0) && (cityChar[rh]<cityChar[dbh])
  {
    damagevector[dvz1]=CalculateDamageResiliantUnprotectedZone1(surge,cityChar);
    damagevector[dvFE]=1;
    if (surge>cityChar[tz1])
    {
      damagevector[dvz2]=CalculateDamageNonResiliantUnprotectedZone(surge,cityChar);
      if (surge>cityChar[tz2])
      {
        damagevector[dvz3]=CalculateDamageProtectedZone(surge,cityChar);
        if (surge>cityChar[tz3]) damagevector[dvz4]=CalculateDamageAboveDikeProtectionZone(surge,cityChar);
      }
    }
    break;
  }
    case 2: // (cityChar[dh]>0.0) && (cityChar[dbh]>0.0) && (cityChar[rh]>0.0) && (cityChar[rh]>=cityChar[dbh])
  {
    damagevector[dvz1]=CalculateDamageResiliantUnprotectedZone2(surge,cityChar);
    damagevector[dvFE]=1;
    // damagevector[dvz2]=0; //there is no nonresilient unprotected zone
    if (surge>cityChar[tz2])
    {
      damagevector[dvz3]=CalculateDamageProtectedZone(surge,cityChar);
      if (surge>cityChar[tz3]) damagevector[dvz4]=CalculateDamageAboveDikeProtectionZone(surge,cityChar);
    }
    break;
  }
    case 3:       // (cityChar[dh]>0.0) && (cityChar[dbh]>0.0) && (cityChar[rh]=0.0)
  {
    // damagevector[dvz1]=0; //zone 1 does not exist
    damagevector[dvz2]=CalculateDamageNonResiliantUnprotectedZone(surge,cityChar);
    damagevector[dvFE]=1.0;
    if (surge>cityChar[tz2])
    {
      damagevector[dvz3]=CalculateDamageProtectedZone(surge,cityChar);
      if (surge>cityChar[tz3]) damagevector[dvz4]=CalculateDamageAboveDikeProtectionZone(surge,cityChar);
    }
    break;
  }
    case 4:       // (cityChar[dh]>0.0) && (cityChar[dbh]=0.0)
  {
    // damagevector[dvz1]=0;
    // damagevecotor[2=0;]
    if (surge>cityChar[tz2])
  {
    damagevector[dvz3]=CalculateDamageProtectedZone(surge,cityChar);
    if (surge>cityChar[tz3]) damagevector[dvz4]=CalculateDamageAboveDikeProtectionZone(surge,cityChar);
    if (damagevector[dvz3]>0) damagevector[dvFE]=1;
  }
    break;
  }
    case 5:       // (cityChar[dh]=0.0) && (cityChar[dbh]>0.0) && (cityChar[rh]>0.0) && (cityChar[rh]<cityChar[dbh])
  {
    damagevector[dvz1]=CalculateDamageResiliantUnprotectedZone1(surge,cityChar);
    damagevector[dvFE]=1;
    if (surge>cityChar[tz1])
    {
      damagevector[dvz2]=CalculateDamageNonResiliantUnprotectedZone(surge,cityChar);
      if (surge>cityChar[tz3])
      {
        // damagevector[dvz3]=0;
        damagevector[dvz4]=CalculateDamageAboveDikeProtectionZone(surge,cityChar);
      }
    }
    break;
  }
    case 6:       // (cityChar[dh]=0.0) && (cityChar[dbh]>0.0) && (cityChar[rh]>0.0) && (cityChar[rh]>=cityChar[dbh])
  {
    damagevector[dvz1]=CalculateDamageResiliantUnprotectedZone2(surge,cityChar);
    damagevector[dvFE]=1.0;
    // damagevector[dvz2]=0;
    if (surge>cityChar[tz3])
    {
      //damagevector[dvz3]=0;
      damagevector[dvz4]=CalculateDamageAboveDikeProtectionZone(surge,cityChar);
    }
    break;
  }
    case 7:       // (cityChar[dh]=0.0) && (cityChar[dbh]>0.0) && (cityChar[rh]=0.0)
  {
    damagevector[dvz2]=CalculateDamageNonResiliantUnprotectedZone(surge,cityChar);
    damagevector[dvFE]=1.0;
    if (surge>cityChar[tz3])
    {
      // damagevector[dvz3]=0;
      damagevector[dvz4]=CalculateDamageAboveDikeProtectionZone(surge,cityChar);
    }
    break;
  }
    case 8:       // (cityChar[dh]=0.0) && (cityChar[dbh]=0.0) && (cityChar[rh]>0.0)
  {
    damagevector[dvFE]=1.0;
    damagevector[dvz1]=CalculateDamageResiliantUnprotectedZone1(surge,cityChar);
    // damagevector[dvz2]=0;
    // damagevector[dvz3]=0
    if (surge>cityChar[tz1]) damagevector[dvz4]=CalculateDamageAboveDikeProtectionZone(surge,cityChar);
    break;
  }
    case 9:       // (cityChar[dh]=0.0) && (cityChar[dbh]=0.0) && (cityChar[rh]=0.0)
  {
    // damagevector[dvz1]=0;
    // damagevector[dvz2]=0;
    // damagevector[dvz3]=0;
    damagevector[dvFE]=1.0;
    damagevector[dvz4]=CalculateDamageAboveDikeProtectionZone(surge,cityChar);
    break;
  }
      
    }
    damagevector[dvt]=damagevector[dvz1]+damagevector[dvz2]+damagevector[dvz3]+damagevector[dvz4];
    if (damagevector[dvz3]>0) damagevector[dvBE]=1;
    if(damagevector[dvt]>threshold) {
      damagevector[dvTE]=1.0;
      damagevector[dvt]=damagevector[dvt]+
        pow(thresholdDamageFraction*(damagevector[dvt]-threshold),thresholdDamageExponent);
    }
    
    
    
  }
}



void iCOWmodule(double* levers, double* objs, double* info, double* consts)  
{
  double W=levers[0];
  double R=levers[1];
  double P=levers[2];
  double B=levers[3];
  double D=levers[4];
  
  double surge;
  
  int CountPositiveNPV = 0;
  // defended city
  double city[numCityChar];
  double damageVector[dvLength]={0,0,0,0,0,0,0,0};
  double damageVectorAccumulated[dvLength]={0,0,0,0,0,0,0,0};
  
  double seqDamage;
  
  // baseline city
  double city0[numCityChar];
  double damageVector0[dvLength]={0,0,0,0,0,0,0,0};
  double damageVectorAccumulated0[dvLength]={0,0,0,0,0,0,0,0};
  
  double seqDamage0;
  
  
  
  CharacterizeCity (0,0,0,P,0,city0);        
  CharacterizeCity (W,B,R,P,D,city);
  
  for (int s=0; s<statesToEvaluate; s++)                
  {
    seqDamage=0;
    seqDamage0=0;
    
    for (int y=0; y<years; y++)
    {
      memset(damageVector, 0, sizeof(damageVector));
      memset(damageVector0, 0, sizeof(damageVector0));
      surge=surges[s*lengthSurgeSequences+y];
      if (surge>Seawall){
        //undefended city
        CalculateDamageVector(surge,city0, damageVector0);
        damageVectorAccumulated0[dvFE]=damageVectorAccumulated0[dvFE]+damageVector0[dvFE];
        damageVectorAccumulated0[dvt]=damageVectorAccumulated0[dvt]+damageVector0[dvt];
        damageVectorAccumulated0[dvz1]=damageVectorAccumulated0[dvz1]+damageVector0[dvz1];
        seqDamage0=seqDamage0+damageVector0[dvt];
        damageVectorAccumulated0[dvBE]=damageVectorAccumulated0[dvBE]+damageVector0[dvBE];
        damageVectorAccumulated0[dvTE]=damageVectorAccumulated0[dvTE]+damageVector0[dvTE];
        
        //defended city
        CalculateDamageVector(surge,city, damageVector);
        damageVectorAccumulated[dvFE]=damageVectorAccumulated[dvFE]+damageVector[dvFE];
        damageVectorAccumulated[dvt]=damageVectorAccumulated[dvt]+damageVector[dvt];
        damageVectorAccumulated[dvz1]=damageVectorAccumulated[dvz1]+damageVector[dvz1];
        seqDamage=seqDamage+damageVector[dvt];
        damageVectorAccumulated[dvBE]=damageVectorAccumulated[dvBE]+damageVector[dvBE];
        damageVectorAccumulated[dvTE]=damageVectorAccumulated[dvTE]+damageVector[dvTE];
        
      }
    } //every year  
    if (seqDamage0==seqDamage)
    {
      seqDamage0=seqDamage0+1;
    }
    if (damageVector0[dvt]-damageVector[dvt]-city[tic]>0)
    {
      CountPositiveNPV=CountPositiveNPV+1;
    }
  } //every stateToEvaluate  
  //if (damageVectorAccumulated[dvt]==0) damageVectorAccumulated[dvt]=1.0;
  //if (city[dh]>1&&(0.1>city[dbh]-city[wh])) consts[0]=city[rh]; else consts[0]=0;
  if (damageVectorAccumulated[dvt]<damageVectorAccumulated0[dvt])
  {
    consts[0]=0.0;
  }
  else
  {
    consts[0]=(damageVectorAccumulated[dvt]-damageVectorAccumulated0[dvt])/statesToEvaluate/costUnit;
  }
  objs[0]=damageVectorAccumulated[dvt]/statesToEvaluate/costUnit;
  objs[1]=(city[dc]+city[rc]+city[wc])/costUnit;//total investment cost
  objs[2]=(double)(statesToEvaluate-CountPositiveNPV)/statesToEvaluate;//Freq Negative NPV
  objs[3]=-((damageVectorAccumulated0[dvt]-damageVectorAccumulated[dvt])/statesToEvaluate)/city[tic]; //negative ROI
  objs[4]=damageVectorAccumulated[dvTE]/statesToEvaluate/years; //Freq threashold events
  objs[5]=((damageVectorAccumulated[dvt]-damageVectorAccumulated0[dvt])/statesToEvaluate+city[tic])/costUnit; //neg npv
  objs[6]=city[dc]/costUnit;
  objs[7]=city[wc]/costUnit;
  objs[8]=city[rc]/costUnit;
  info[0]=(city[dc]+city[rc]+city[wc])/costUnit;//total investment cost
  info[1]=objs[0]/years;//total damage
  info[2]=-objs[3]; // average ROI
  info[3]=((damageVectorAccumulated0[dvt]-damageVectorAccumulated[dvt])/statesToEvaluate-city[tic])/costUnit; //average NPV
  info[4]=damageVectorAccumulated[dvTE]/statesToEvaluate/years; //freq of threashold events
  info[5]=(double)CountPositiveNPV/statesToEvaluate; //frequency positive NPV > 1
  info[6]=city[wc]/costUnit;
  info[7]=city[rc]/costUnit;
  info[8]=city[dc]/costUnit;
  info[9]=city[wh];
  info[10]=city[rh];
  info[11]=city[rp];
  info[12]=city[dbh];
  info[13]=city[dh];        
}




int main(int argc, char* argv[]) {
  double location = atof(argv[3]);
  double deltaLocation = atof(argv[4]); // in meters per year
  double scale =atof(argv[5]);  // in meters
  double deltaScale = atof(argv[6]); // in meters per year
  double shape = atof(argv[7]);
  double deltaShape = atof(argv[8]);
  double totTime = atof(argv[9]);
  statesToEvaluate = atoi(argv[10]);
  double damageE = atof(argv[11]);
  double dcE = atof(argv[12]);
  double freqNegNPVE = atof(argv[13]);
  double negBCRE = atof(argv[14]);
  double freqTEE = atof(argv[15]);
  double negNPVE = atof(argv[16]);
  double wcE = atof(argv[17]);
  double rcE = atof(argv[18]);
  srand(time(NULL));
  
  
  char runtime[256];
  //char fname[20]=*argv[1];
  double WRange[2] = {WHL,WHH};
  double RRange[2] = {RHL,RHH};  
  double PRange[2] = {PL,PH};
  double BRange[2] = {DBL,DBH};  
  double DRange[2] = {DHL,DHH};  
  
  
  nvars=5;
  nobjs=9;
  ninfo=14;
  nconsts=1;
  int rank;
  cout <<  "generate surges";
  generateSurges(surges, location, scale, shape, deltaLocation, deltaScale, deltaShape);
  cout <<  "surges generated";  
  
  // All master-slave runs need to call startup and set the runtime
  // limits.
  BORG_Algorithm_ms_startup(&argc, &argv);
  BORG_Algorithm_ms_max_time(totTime/10);
  
  // Define the problem.  Problems are defined the same way as the
  // serial example (see dtlz2_serial.c).
  BORG_Problem problem = BORG_Problem_create(nvars, nobjs, ninfo, nconsts, iCOWmodule);
  
  BORG_Problem_set_bounds(problem, 0, WRange[0], WRange[1]);  // W  
  BORG_Problem_set_bounds(problem, 1, RRange[0], RRange[1]);  // R
  BORG_Problem_set_bounds(problem, 2, PRange[0], PRange[1]); // P
  BORG_Problem_set_bounds(problem, 3, BRange[0], BRange[1]);  // B
  BORG_Problem_set_bounds(problem, 4, DRange[0], DRange[1]);  // D
  

  BORG_Problem_set_epsilon(problem, 0,  damageE); // .1 damages
  BORG_Problem_set_epsilon(problem, 1, dcE); // .1 total investment cost
  BORG_Problem_set_epsilon(problem, 2, freqNegNPVE); //f neg npv
  BORG_Problem_set_epsilon(problem, 3, negBCRE); //negBCR .005
  BORG_Problem_set_epsilon(problem, 4, freqTEE);//0.0002); //f TE
  BORG_Problem_set_epsilon(problem, 5, negNPVE); //neg npv 
  BORG_Problem_set_epsilon(problem, 6, wcE);//  withdrawal cost
  BORG_Problem_set_epsilon(problem, 7, rcE); //  resistance cost

 // BORG_Problem_set_epsilon(problem, 0, .05); // .1 damages
//  BORG_Problem_set_epsilon(problem, 1, .05); // .1 total investment cost
//  BORG_Problem_set_epsilon(problem, 2, 0.005); //f neg npv
//  BORG_Problem_set_epsilon(problem, 3, 0.02); //negBCR .005
//  BORG_Problem_set_epsilon(problem, 4, 0.0004);//0.0002); //f TE
//  BORG_Problem_set_epsilon(problem, 5, 0.01); //neg npv

 
  // Get the rank of this process.  The rank is used to ensure each
  // parallel process uses a different random seed.
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // When running experiments, we want to run the algorithm multiple
  // times and average the results.
  for (int i=0; i<10; i++) {
    // Save runtime dynamics to a file.  Only the master node
    // will write to this file.  Note how we create separate
    // files for each run.
    std::string runPath(argv[2]);
    std::string runFileString=runPath+"runtime%d.txt";
    const char *runFile = runFileString.c_str();
    sprintf(runtime, runFile, i);
    BORG_Algorithm_output_runtime(runtime);
    
    // Seed the random number generator.
    BORG_Random_seed(37*i*(rank+1));
    
    // Run the master-slave Borg MOEA on the problem.
    BORG_Archive result = BORG_Algorithm_ms_run(problem);
    
    // Only the master process will return a non-NULL result.
    // Print the Pareto optimal solutions to the screen.
    if (result != NULL) {
      FILE * oFile;
      oFile = fopen (argv[1],"w");
      ICOW_Archive_print(result, oFile);
      fclose (oFile);
      BORG_Archive_destroy(result);
    }
  }
  
  // Shutdown the parallel processes and exit.
  BORG_Algorithm_ms_shutdown();
  BORG_Problem_destroy(problem);

  return EXIT_SUCCESS;
}
