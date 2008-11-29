// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "Constraint.hh" // implementation of object methods

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::Constraint::Constraint(void) :
  _normalizer(new spatialdata::units::Nondimensional),
  _useSolnIncr(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::feassemble::Constraint::~Constraint(void)
{ // destructor
  delete _normalizer; _normalizer = 0;
} // destructor

// ----------------------------------------------------------------------
// Set manager of scales used to nondimensionalize problem.
void
pylith::feassemble::Constraint::normalizer(const spatialdata::units::Nondimensional& dim)
{ // normalizer
  if (0 == _normalizer)
    _normalizer = new spatialdata::units::Nondimensional(dim);
  else
    *_normalizer = dim;
} // normalizer


// End of file 
