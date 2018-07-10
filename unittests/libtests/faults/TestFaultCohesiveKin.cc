// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TestFaultCohesiveKin.hh" // Implementation of class methods

#include "pylith/faults/FaultCohesiveKin.hh" // USES FaultCohesiveKin

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps::nondimensionalize()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/feassemble/AuxiliaryFactory.hh" // USES AuxiliaryFactory
#include "pylith/topology/VisitorMesh.hh" // USES VisitorMesh
#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii
#include "pylith/problems/SolutionFactory.hh" // USES SolutionFactory

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
// Setup testing data.
void
pylith::faults::TestFaultCohesiveKin::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _fault = new pylith::faults::FaultCohesiveKin();CPPUNIT_ASSERT(_fault);

    _data = NULL;
    _mesh = NULL;
    _solution = NULL;

    PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
// Tear down testing data.
void
pylith::faults::TestFaultCohesiveKin::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _fault;_fault = NULL;

    delete _data;_data = NULL;
    delete _mesh;_mesh = NULL;
    delete _solution;_solution = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ----------------------------------------------------------------------
// Test constructor.
void
pylith::faults::TestFaultCohesiveKin::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    FaultCohesiveKin* fault = new FaultCohesiveKin();CPPUNIT_ASSERT(fault);

    CPPUNIT_ASSERT_EQUAL(100, fault->id());

    const PylithReal tolerance = 1.0e-6;
    const PylithReal refDir1E[3] = { 0.0, 0.0, 1.0 };
    const PylithReal* refDir1 = fault->refDir1();
    for (int i = 0; i < 3; ++i) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(refDir1E[i], refDir1[i], tolerance);
    } // for

    const PylithReal refDir2E[3] = { 0.0, 1.0, 0.0 };
    const PylithReal* refDir2 = fault->refDir2();
    for (int i = 0; i < 3; ++i) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(refDir2E[i], refDir2[i], tolerance);
    } // for

    delete fault;fault = NULL;

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
/// Test accessors (dbTimeHistory, useInitial, useRate, useTimeHistory).
void
pylith::faults::TestFaultCohesiveKin::testAccessors(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_fault);

    // id()
    const PylithInt id = 543;
    _fault->id(id);
    CPPUNIT_ASSERT_EQUAL(id, _fault->id());

    // label()
    const std::string& faces = "fault_faces";
    _fault->label(faces.c_str());
    CPPUNIT_ASSERT_EQUAL(faces, std::string(_fault->label()));

    // edge()
    const std::string& edges = "buried_edges";
    _fault->edge(edges.c_str());
    CPPUNIT_ASSERT_EQUAL(edges, std::string(_fault->edge()));

    // refDir1(), refDir2
    const PylithReal tolerance = 1.0e-6;
    const PylithReal magDir = sqrt(1.1*1.1 + 0.1*0.1 + 0.2*0.2);
    const PylithReal refDir1E[3] = { 1.1/magDir, 0.1/magDir, 0.2/magDir };
    const PylithReal refDir2E[3] = { 0.1/magDir, 0.2/magDir, 1.1/magDir };
    _fault->refDir1(refDir1E);
    _fault->refDir2(refDir2E);

    const PylithReal* refDir1 = _fault->refDir1();
    const PylithReal* refDir2 = _fault->refDir2();
    for (int i = 0; i < 3; ++i) {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(refDir1E[i], refDir1[i], tolerance);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(refDir2E[i], refDir2[i], tolerance);
    } // for

    // eqsrcs()
    CPPUNIT_ASSERT_MESSAGE(":TODO: Implement test for eqsrcs().", false);

    PYLITH_METHOD_END;
} // testAccessors


// ----------------------------------------------------------------------
void
pylith::faults::TestFaultCohesiveKin::testAuxFieldDiscretization(void) {
    PYLITH_METHOD_BEGIN;

    const topology::FieldBase::Discretization infoDefault = pylith::topology::Field::Discretization(1, 1, true, pylith::topology::FieldBase::POLYNOMIAL_SPACE);
    const topology::FieldBase::Discretization infoA = pylith::topology::Field::Discretization(1, 2, false, pylith::topology::FieldBase::POLYNOMIAL_SPACE);
    const topology::FieldBase::Discretization infoB = pylith::topology::Field::Discretization(2, 2, true, pylith::topology::FieldBase::POINT_SPACE);

    CPPUNIT_ASSERT(_fault);
    _fault->auxSubfieldDiscretization("A", infoA.basisOrder, infoA.quadOrder, infoA.isBasisContinuous, infoA.feSpace);
    _fault->auxSubfieldDiscretization("B", infoB.basisOrder, infoB.quadOrder, infoB.isBasisContinuous, infoB.feSpace);

    CPPUNIT_ASSERT(_fault->_auxFactory());
    { // A
        const topology::FieldBase::Discretization& test = _fault->_auxFactory()->subfieldDiscretization("A");
        CPPUNIT_ASSERT_EQUAL(infoA.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoA.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoA.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoA.feSpace, test.feSpace);
    } // A

    { // B
        const topology::FieldBase::Discretization& test = _fault->_auxFactory()->subfieldDiscretization("B");
        CPPUNIT_ASSERT_EQUAL(infoB.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoB.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoB.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoB.feSpace, test.feSpace);
    } // B

    { // C (default)
        const topology::FieldBase::Discretization& test = _fault->_auxFactory()->subfieldDiscretization("C");
        CPPUNIT_ASSERT_EQUAL(infoDefault.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoDefault.feSpace, test.feSpace);
    } // C (default)

    { // default
        const topology::FieldBase::Discretization& test = _fault->_auxFactory()->subfieldDiscretization("default");
        CPPUNIT_ASSERT_EQUAL(infoDefault.basisOrder, test.basisOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.quadOrder, test.quadOrder);
        CPPUNIT_ASSERT_EQUAL(infoDefault.isBasisContinuous, test.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(infoDefault.feSpace, test.feSpace);
    } // default

    PYLITH_METHOD_END;
} // testAuxFieldDiscretization


// ----------------------------------------------------------------------
// Test auxFieldDB().
void
pylith::faults::TestFaultCohesiveKin::testAuxFieldDB(void) {
    PYLITH_METHOD_BEGIN;

    const std::string label = "test db";
    spatialdata::spatialdb::UserFunctionDB db;
    db.label(label.c_str());

    CPPUNIT_ASSERT(_fault);
    _fault->auxFieldDB(&db);

    CPPUNIT_ASSERT(_fault->_auxFactory());
    CPPUNIT_ASSERT(_fault->_auxFactory()->queryDB());

    CPPUNIT_ASSERT_EQUAL(label, std::string(_fault->_auxFactory()->queryDB()->label()));

    PYLITH_METHOD_END;
} // testAuxFieldDB


// ----------------------------------------------------------------------
// Test normalizer().
void
pylith::faults::TestFaultCohesiveKin::testNormalizer(void) {
    PYLITH_METHOD_BEGIN;

    spatialdata::units::Nondimensional normalizer;
    const double scale = 5.0;
    normalizer.lengthScale(scale);

    CPPUNIT_ASSERT(_fault);
    _fault->normalizer(normalizer);
    CPPUNIT_ASSERT_EQUAL(scale, _fault->_normalizer->lengthScale());

    PYLITH_METHOD_END;
} // testNormalizer


// ----------------------------------------------------------------------
// Test verifyConfiguration().
void
pylith::faults::TestFaultCohesiveKin::testVerifyConfiguration(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    _setupSolutionField();

    CPPUNIT_ASSERT(_fault);
    CPPUNIT_ASSERT(_solution);

    // Verify should pass.
    _fault->verifyConfiguration(*_solution);

    // Check for failure with displacement field or Lagrange multiplier field not in solution.
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad Implement testVerifyConfiguration() without displacement or Lagrange multiplier field.", false);
    //CPPUNIT_ASSERT_THROW(_fault->verifyConfiguration(displacement), std::runtime_error);
    //CPPUNIT_ASSERT_THROW(_fault->verifyConfiguration(lagrange), std::runtime_error);

    PYLITH_METHOD_END;
} // testVerifyConfiguration


// ----------------------------------------------------------------------
// Test initialize().
void
pylith::faults::TestFaultCohesiveKin::testInitialize(void) {
    PYLITH_METHOD_BEGIN;

    // Call initialize()
    _initialize(); // includes setting up auxField
    _setupSolutionField();

    CPPUNIT_ASSERT(_fault);
    CPPUNIT_ASSERT(_solution);

    _fault->initialize(*_solution);

    //_fault->auxField().view("AUXILIARY FIELD"); // DEBUGGING

    // Verify auxiliary field.
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_mesh);
    const pylith::topology::Field& auxField = _fault->auxField();
    CPPUNIT_ASSERT_EQUAL(std::string("auxiliary subfields"), std::string(auxField.label()));
    CPPUNIT_ASSERT_EQUAL(_mesh->dimension(), auxField.spaceDim());

    PylithReal norm = 0.0;
    PylithReal t = _data->t;
    const PetscDM dm = auxField.dmMesh();CPPUNIT_ASSERT(dm);
    pylith::topology::FieldQuery query(auxField);
    query.initializeWithDefaultQueryFns();
    CPPUNIT_ASSERT(_data->normalizer);
    query.openDB(_data->auxDB, _data->normalizer->lengthScale());
    PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, query.functions(), (void**)query.contextPtrs(), auxField.localVector(), &norm);CPPUNIT_ASSERT(!err);
    query.closeDB(_data->auxDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);

    PYLITH_METHOD_END;
} // testInitialize


// ----------------------------------------------------------------------
// Test preStep().
void
pylith::faults::TestFaultCohesiveKin::testPrestep(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    _setupSolutionField();

    CPPUNIT_ASSERT(_fault);
    CPPUNIT_ASSERT(_solution);
    _fault->initialize(*_solution);

    const pylith::topology::Field& auxField = _fault->auxField();
    CPPUNIT_ASSERT(auxField.hasSubfield("time_history_value"));

    pylith::topology::Field valueField(*_mesh);
    valueField.copySubfield(auxField, "time_history_value");
    CPPUNIT_ASSERT(valueField.sectionSize() > 0);
    CPPUNIT_ASSERT_EQUAL(std::string("time_history_value"), std::string(valueField.label()));

    PylithReal norm = 0.0;
    PylithReal t = _data->t;
    const PetscDM dm = auxField.dmMesh();CPPUNIT_ASSERT(dm);
    pylith::topology::FieldQuery query(valueField);
    query.initializeWithDefaultQueryFns();
    CPPUNIT_ASSERT(_data->normalizer);
    query.openDB(_data->auxDB, _data->normalizer->lengthScale());
    PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, query.functions(), (void**)query.contextPtrs(), valueField.localVector(), &norm);CPPUNIT_ASSERT(!err);
    query.closeDB(_data->auxDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);

    PYLITH_METHOD_END;
} // testPrestep


// ----------------------------------------------------------------------
// Test setSolution().
void
pylith::faults::TestFaultCohesiveKin::testComputeRHSResidual(void) {
    PYLITH_METHOD_BEGIN;

#if 1
    CPPUNIT_ASSERT_MESSAGE(":TODO: @brad not implemented.", false);
#else
    CPPUNIT_ASSERT(_fault);
    CPPUNIT_ASSERT(_solution);
    _fault->initialize(*_solution);

    // Initialize solution field.
    _solution->allocate();
    _solution->zeroLocal();
    _solution->createScatter(_solution->mesh(), "global");

    // Set solution field.
    CPPUNIT_ASSERT(_data);
    _fault->setSolution(*_solution, _data->t);

    // Verify number and DOF of constraints in solution field.
    int iConstraint = 0;
    PetscErrorCode err = 0;
    for (PetscInt v = vStart; v < vEnd; ++v) {
        PetscInt dof, cdof, fdof, fcdof;

        err = PetscSectionGetDof(fieldSection, v, &dof);PYLITH_CHECK_ERROR(err);
        err = PetscSectionGetConstraintDof(fieldSection, v, &cdof);PYLITH_CHECK_ERROR(err);
        err = PetscSectionGetFieldDof(fieldSection, v, 0, &fdof);PYLITH_CHECK_ERROR(err);
        err = PetscSectionGetFieldConstraintDof(fieldSection, v, 0, &fcdof);PYLITH_CHECK_ERROR(err);
        if (v != _data->constrainedPoints[iConstraint] + offset) {
            CPPUNIT_ASSERT_EQUAL(_data->numDOF, dof);
            CPPUNIT_ASSERT_EQUAL(0, cdof);
            CPPUNIT_ASSERT_EQUAL(_data->numDOF, fdof);
            CPPUNIT_ASSERT_EQUAL(0, fcdof);
        } else {
            CPPUNIT_ASSERT_EQUAL(_data->numDOF, dof);
            CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, cdof);
            CPPUNIT_ASSERT_EQUAL(_data->numDOF, fdof);
            CPPUNIT_ASSERT_EQUAL(_data->numFixedDOF, fcdof);
            ++iConstraint;
        } // if/else
    } // for

    // Verify values in solution field.
    PylithReal norm = 0.0;
    PylithReal t = _data->t;
    const PetscDM dmSoln = _solution->dmMesh();CPPUNIT_ASSERT(dmSoln);
    pylith::topology::FieldQuery* query = _db->_auxFieldsQuery;
    query->openDB(queryDB, _data->lengthScale);

    PetscErrorCode err = DMPlexComputeL2DiffLocal(dm, t, query->functions(), (void**)query->contextPtrs(), _solution->localVector(), &norm);CPPUNIT_ASSERT(!err);
    query->closeDB(queryDB);
    const PylithReal tolerance = 1.0e-6;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, norm, tolerance);
#endif

    PYLITH_METHOD_END;
} // testComputeRHSResidual


// ----------------------------------------------------------------------
// Test _auxFieldSetup().
void
pylith::faults::TestFaultCohesiveKin::testAuxFieldSetup(void) {
    PYLITH_METHOD_BEGIN;

    _initialize();
    _setupSolutionField();

    CPPUNIT_ASSERT(_fault);
    CPPUNIT_ASSERT(_solution);
    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);
    const PylithReal pressureScale = _data->normalizer->pressureScale();

    delete _fault->_faultMesh;_fault->_faultMesh = new pylith::topology::Mesh(_solution->mesh(), _data->faultLabel);
    CPPUNIT_ASSERT(_fault->_faultMesh);

    delete _fault->_auxField;_fault->_auxField = new pylith::topology::Field(*_fault->_faultMesh);CPPUNIT_ASSERT(_fault->_auxField);
    _fault->_auxFieldSetup();

    CPPUNIT_ASSERT(_mesh->coordsys());
    const size_t spaceDim = _mesh->coordsys()->spaceDim();
    const size_t numComponents = spaceDim;

    // Check discretizations
    int ifield = 0;
    { // slip}
        const char* label = "slip";
        CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(_data->auxSubfields[ifield]));
        const pylith::topology::Field::Discretization& discretization = _data->auxDiscretizations[ifield];

        const pylith::topology::Field::SubfieldInfo& info = _fault->_auxField->subfieldInfo(label);
        CPPUNIT_ASSERT_EQUAL(numComponents, info.description.numComponents);
        CPPUNIT_ASSERT_EQUAL(std::string(label), info.description.label);
        CPPUNIT_ASSERT_EQUAL(pylith::topology::Field::VECTOR, info.description.vectorFieldType);
        CPPUNIT_ASSERT_EQUAL(pressureScale, info.description.scale);
        CPPUNIT_ASSERT_EQUAL(discretization.basisOrder, info.fe.basisOrder);
        CPPUNIT_ASSERT_EQUAL(discretization.quadOrder, info.fe.quadOrder);
        CPPUNIT_ASSERT_EQUAL(discretization.isBasisContinuous, info.fe.isBasisContinuous);
        CPPUNIT_ASSERT_EQUAL(discretization.feSpace, info.fe.feSpace);
        ++ifield;
    } // slip

    PYLITH_METHOD_END;
} // testAuxFieldSetup


// ----------------------------------------------------------------------
void
pylith::faults::TestFaultCohesiveKin::_initialize(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_data);

    delete _mesh;_mesh = new pylith::topology::Mesh();CPPUNIT_ASSERT(_mesh);
    pylith::meshio::MeshIOAscii iohandler;
    CPPUNIT_ASSERT(_data->meshFilename);
    iohandler.filename(_data->meshFilename);
    iohandler.read(_mesh);
    _mesh->coordsys(_data->cs);
    CPPUNIT_ASSERT(_data->normalizer);
    pylith::topology::MeshOps::nondimensionalize(_mesh, *_data->normalizer);

    _fault->label(_data->faultLabel);
    _fault->auxFieldDB(_data->auxDB);
    _fault->normalizer(*_data->normalizer);
    for (int ifield = 0; ifield < _data->numAuxSubfields; ++ifield) {
        const pylith::topology::Field::Discretization& discretization = _data->auxDiscretizations[ifield];
        const char* name = _data->auxSubfields[ifield];
        _fault->auxSubfieldDiscretization(name, discretization.basisOrder, discretization.quadOrder, discretization.isBasisContinuous, discretization.feSpace);
    } // for

    //_fault->eqsrcs(_data->eqsrcs);

    PYLITH_METHOD_END;
} // _initialize


// ----------------------------------------------------------------------
void
pylith::faults::TestFaultCohesiveKin::_setupSolutionField(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    CPPUNIT_ASSERT(_data);
    CPPUNIT_ASSERT(_data->normalizer);

    delete _solution;_solution = new pylith::topology::Field(*_mesh);
    _solution->label("Solution (displacement, lagrange_multiplier_fault)");
    pylith::problems::SolutionFactory factory(*_solution, *_data->normalizer);
    factory.displacement(_data->solnDiscretizations[0]);
    factory.lagrangeMultiplierFault(_data->solnDiscretizations[1]);
    _solution->subfieldsSetup();

    PYLITH_METHOD_END;
} // setupSolutionField


// ----------------------------------------------------------------------
// Constructor
pylith::faults::TestFaultCohesiveKin_Data::TestFaultCohesiveKin_Data(void) :
    meshFilename(NULL),
    faultLabel(NULL),
    normalizer(new spatialdata::units::Nondimensional),
    numAuxSubfields(0),
    auxSubfields(NULL),
    auxDiscretizations(NULL),
    auxDB(new spatialdata::spatialdb::UserFunctionDB),
    t(0.0),
    solnNumSubfields(0),
    solnDiscretizations(NULL),
    solnDB(new spatialdata::spatialdb::UserFunctionDB){ // constructor
} // constructor


// ----------------------------------------------------------------------
// Destructor
pylith::faults::TestFaultCohesiveKin_Data::~TestFaultCohesiveKin_Data(void) {
    delete cs;cs = NULL;
    delete normalizer;normalizer = NULL;
    delete auxDB;auxDB = NULL;
    delete solnDB;solnDB = NULL;
} // destructor


// End of file
