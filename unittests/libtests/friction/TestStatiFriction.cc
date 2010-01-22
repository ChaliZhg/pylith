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

#include "TestElasticIsotropic3D.hh" // Implementation of class methods

#include "data/ElasticIsotropic3DData.hh" // USES ElasticIsotropic3DData

#include "pylith/materials/ElasticIsotropic3D.hh" // USES ElasticIsotropic3D

#include "pylith/utils/constdefs.h" // USES MAXDOUBLE
#include "pylith/topology/Mesh.hh" // USES Mesh

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION( pylith::materials::TestElasticIsotropic3D );

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::materials::TestElasticIsotropic3D::setUp(void)
{ // setUp
  _material = new ElasticIsotropic3D();
  _matElastic = new ElasticIsotropic3D();
  _data = new ElasticIsotropic3DData();
  _dataElastic = new ElasticIsotropic3DData();
  setupNormalizer();
} // setUp

// ----------------------------------------------------------------------
// Test stableTimeStepImplicit().
void
pylith::materials::TestElasticIsotropic3D::testStableTimeStepImplicit(void)
{ // testStableTimeStepImplicit
  assert(0 != _matElastic);

  topology::Mesh mesh;

  const double dt = _matElastic->stableTimeStepImplicit(mesh);
  CPPUNIT_ASSERT_EQUAL(1.0, dt/pylith::PYLITH_MAXDOUBLE);
} // testStableTimeStepImplicit

// ----------------------------------------------------------------------
// Test hasProperty().
void
pylith::materials::TestElasticIsotropic3D::testHasProperty(void)
{ // testHasProperty
  ElasticIsotropic3D material;

  CPPUNIT_ASSERT(material.hasProperty("mu"));
  CPPUNIT_ASSERT(material.hasProperty("lambda"));
  CPPUNIT_ASSERT(material.hasProperty("density"));
  CPPUNIT_ASSERT(!material.hasProperty("aaa"));
} // testHasProperty

// ----------------------------------------------------------------------
// Test hasStateVar().
void
pylith::materials::TestElasticIsotropic3D::testHasStateVar(void)
{ // testHasStateVar
  ElasticIsotropic3D material;

  CPPUNIT_ASSERT(!material.hasStateVar("stress"));
} // testHasStateVar


// End of file 
