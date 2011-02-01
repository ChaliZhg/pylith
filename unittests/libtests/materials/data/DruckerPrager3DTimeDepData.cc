// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

// DO NOT EDIT THIS FILE
// This file was generated from python application druckerprager3dtimedep.

#include "DruckerPrager3DTimeDepData.hh"

const int pylith::materials::DruckerPrager3DTimeDepData::_dimension = 3;

const int pylith::materials::DruckerPrager3DTimeDepData::_numLocs = 2;

const int pylith::materials::DruckerPrager3DTimeDepData::_numProperties = 6;

const int pylith::materials::DruckerPrager3DTimeDepData::_numStateVars = 1;

const int pylith::materials::DruckerPrager3DTimeDepData::_numDBProperties = 6;

const int pylith::materials::DruckerPrager3DTimeDepData::_numDBStateVars = 6;

const int pylith::materials::DruckerPrager3DTimeDepData::_numPropsQuadPt = 6;

const int pylith::materials::DruckerPrager3DTimeDepData::_numVarsQuadPt = 6;

const double pylith::materials::DruckerPrager3DTimeDepData::_lengthScale =   1.00000000e+03;

const double pylith::materials::DruckerPrager3DTimeDepData::_timeScale =   1.00000000e+00;

const double pylith::materials::DruckerPrager3DTimeDepData::_pressureScale =   2.25000000e+10;

const double pylith::materials::DruckerPrager3DTimeDepData::_densityScale =   1.00000000e+03;

const double pylith::materials::DruckerPrager3DTimeDepData::_dtStableImplicit =   1.00000000e+10;

const int pylith::materials::DruckerPrager3DTimeDepData::_numPropertyValues[] = {
1,
1,
1,
1,
1,
1,
};

const int pylith::materials::DruckerPrager3DTimeDepData::_numStateVarValues[] = {
6,
};

const char* pylith::materials::DruckerPrager3DTimeDepData::_dbPropertyValues[] = {
"density",
"vs",
"vp",
"friction-angle",
"cohesion",
"dilatation-angle",
};

const char* pylith::materials::DruckerPrager3DTimeDepData::_dbStateVarValues[] = {
"plastic-strain-xx",
"plastic-strain-yy",
"plastic-strain-zz",
"plastic-strain-xy",
"plastic-strain-yz",
"plastic-strain-xz",
};

const double pylith::materials::DruckerPrager3DTimeDepData::_dbProperties[] = {
  2.50000000e+03,
  3.00000000e+03,
  5.19615242e+03,
  5.23598776e-01,
  3.00000000e+05,
  3.49065850e-01,
  2.00000000e+03,
  1.20000000e+03,
  2.07846097e+03,
  4.36332313e-01,
  1.00000000e+04,
  4.36332313e-01,
};

const double pylith::materials::DruckerPrager3DTimeDepData::_dbStateVars[] = {
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
  0.00000000e+00,
};

const double pylith::materials::DruckerPrager3DTimeDepData::_properties[] = {
  2.50000000e+03,
  2.25000000e+10,
  2.25000000e+10,
  2.30940108e-01,
  3.60000000e+05,
  1.48583084e-01,
  2.00000000e+03,
  2.88000000e+09,
  2.88000000e+09,
  1.89338478e-01,
  1.21811303e+04,
  1.89338478e-01,
};

const double pylith::materials::DruckerPrager3DTimeDepData::_stateVars[] = {
  4.10000000e-05,
  4.20000000e-05,
  4.30000000e-05,
  4.40000000e-05,
  4.50000000e-05,
  4.60000000e-05,
  1.10000000e-05,
  1.20000000e-05,
  1.30000000e-05,
  1.40000000e-05,
  1.50000000e-05,
  1.60000000e-05,
};

const double pylith::materials::DruckerPrager3DTimeDepData::_propertiesNondim[] = {
  2.50000000e+00,
  1.00000000e+00,
  1.00000000e+00,
  2.30940108e-01,
  1.60000000e-05,
  1.48583084e-01,
  2.00000000e+00,
  1.28000000e-01,
  1.28000000e-01,
  1.89338478e-01,
  5.41383567e-07,
  1.89338478e-01,
};

const double pylith::materials::DruckerPrager3DTimeDepData::_stateVarsNondim[] = {
  4.10000000e-05,
  4.20000000e-05,
  4.30000000e-05,
  4.40000000e-05,
  4.50000000e-05,
  4.60000000e-05,
  1.10000000e-05,
  1.20000000e-05,
  1.30000000e-05,
  1.40000000e-05,
  1.50000000e-05,
  1.60000000e-05,
};

const double pylith::materials::DruckerPrager3DTimeDepData::_density[] = {
  2.50000000e+03,
  2.00000000e+03,
};

const double pylith::materials::DruckerPrager3DTimeDepData::_strain[] = {
  1.10000000e-04,
  1.20000000e-04,
  1.30000000e-04,
  1.40000000e-04,
  1.50000000e-04,
  1.60000000e-04,
  4.10000000e-04,
  4.20000000e-04,
  4.30000000e-04,
  4.40000000e-04,
  4.50000000e-04,
  4.60000000e-04,
};

const double pylith::materials::DruckerPrager3DTimeDepData::_stress[] = {
  5.55066469e+05,
  5.54377286e+05,
  5.53688104e+05,
 -4.36889895e+03,
 -5.05808169e+03,
 -5.74726442e+03,
  9.73972396e+05,
  9.67659189e+05,
  9.61345981e+05,
 -2.85301193e+05,
 -2.93159178e+05,
 -3.01017164e+05,
};

const double pylith::materials::DruckerPrager3DTimeDepData::_elasticConsts[] = {
  2.69552119e+10,
  2.61171422e+10,
  2.52103070e+10,
 -1.14973070e+10,
 -1.33109780e+10,
 -1.51246488e+10,
  2.56403618e+10,
  2.47114216e+10,
  2.39200125e+10,
 -1.09057219e+10,
 -1.26260719e+10,
 -1.43464216e+10,
  2.42567462e+10,
  2.34432321e+10,
  2.25609525e+10,
 -1.03141367e+10,
 -1.19411658e+10,
 -1.35681945e+10,
 -8.77107957e+09,
 -8.47528712e+09,
 -8.17949468e+09,
  3.68143839e+09,
  4.34178877e+09,
  4.93337366e+09,
 -1.01546951e+10,
 -9.81224227e+09,
 -9.46978915e+09,
  4.34178907e+09,
  4.95792951e+09,
  5.71160063e+09,
 -1.15383109e+10,
 -1.11491971e+10,
 -1.07600836e+10,
  4.93337427e+09,
  5.71160093e+09,
  6.42106241e+09,
  2.65083753e+09,
  3.40097380e+09,
  3.37872031e+09,
 -2.01128150e+09,
 -2.06667779e+09,
 -2.12207396e+09,
  3.40097345e+09,
  2.60654697e+09,
  3.35689774e+09,
 -1.99182658e+09,
 -2.04668706e+09,
 -2.10154743e+09,
  3.37872049e+09,
  3.35689820e+09,
  2.56268622e+09,
 -1.97237171e+09,
 -2.02669634e+09,
 -2.08102085e+09,
 -1.00564075e+09,
 -9.95913200e+08,
 -9.86186147e+08,
  1.06799416e+08,
  9.03403590e+08,
  9.27618763e+08,
 -1.03333886e+09,
 -1.02334327e+09,
 -1.01334814e+09,
  9.03403910e+08,
  1.55897025e+08,
  9.53168172e+08,
 -1.06103672e+09,
 -1.05077369e+09,
 -1.04051051e+09,
  9.27618908e+08,
  9.53168172e+08,
  2.06328375e+08,
};

const double pylith::materials::DruckerPrager3DTimeDepData::_initialStress[] = {
  2.10000000e+04,
  2.20000000e+04,
  2.30000000e+04,
  2.40000000e+04,
  2.50000000e+04,
  2.60000000e+04,
  5.10000000e+04,
  5.20000000e+04,
  5.30000000e+04,
  5.40000000e+04,
  5.50000000e+04,
  5.60000000e+04,
};

const double pylith::materials::DruckerPrager3DTimeDepData::_initialStrain[] = {
  3.60000000e-05,
  3.50000000e-05,
  3.40000000e-05,
  3.30000000e-05,
  3.20000000e-05,
  3.10000000e-05,
  6.10000000e-05,
  6.20000000e-05,
  6.30000000e-05,
  6.60000000e-05,
  6.50000000e-05,
  6.40000000e-05,
};

const double pylith::materials::DruckerPrager3DTimeDepData::_stateVarsUpdated[] = {
  6.92302201e-05,
  8.02677575e-05,
  9.13052948e-05,
  1.07630420e-04,
  1.18667957e-04,
  1.29705495e-04,
  2.84142902e-04,
  2.94412556e-04,
  3.04682210e-04,
  4.32906457e-04,
  4.45444302e-04,
  4.57982146e-04,
};

pylith::materials::DruckerPrager3DTimeDepData::DruckerPrager3DTimeDepData(void)
{ // constructor
  dimension = _dimension;
  numLocs = _numLocs;
  numProperties = _numProperties;
  numStateVars = _numStateVars;
  numDBProperties = _numDBProperties;
  numDBStateVars = _numDBStateVars;
  numPropsQuadPt = _numPropsQuadPt;
  numVarsQuadPt = _numVarsQuadPt;
  lengthScale = _lengthScale;
  timeScale = _timeScale;
  pressureScale = _pressureScale;
  densityScale = _densityScale;
  dtStableImplicit = _dtStableImplicit;
  numPropertyValues = const_cast<int*>(_numPropertyValues);
  numStateVarValues = const_cast<int*>(_numStateVarValues);
  dbPropertyValues = const_cast<char**>(_dbPropertyValues);
  dbStateVarValues = const_cast<char**>(_dbStateVarValues);
  dbProperties = const_cast<double*>(_dbProperties);
  dbStateVars = const_cast<double*>(_dbStateVars);
  properties = const_cast<double*>(_properties);
  stateVars = const_cast<double*>(_stateVars);
  propertiesNondim = const_cast<double*>(_propertiesNondim);
  stateVarsNondim = const_cast<double*>(_stateVarsNondim);
  density = const_cast<double*>(_density);
  strain = const_cast<double*>(_strain);
  stress = const_cast<double*>(_stress);
  elasticConsts = const_cast<double*>(_elasticConsts);
  initialStress = const_cast<double*>(_initialStress);
  initialStrain = const_cast<double*>(_initialStrain);
  stateVarsUpdated = const_cast<double*>(_stateVarsUpdated);
} // constructor

pylith::materials::DruckerPrager3DTimeDepData::~DruckerPrager3DTimeDepData(void)
{}


// End of file
