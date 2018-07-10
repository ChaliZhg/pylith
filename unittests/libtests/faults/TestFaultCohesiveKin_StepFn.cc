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

#include "TestFaultCohesiveKin.hh" // Implementation of cases

#include "pylith/topology/Field.hh" // USES pylith::topology::Field::Discretization
#include "spatialdata/spatialdb/UserFunctionDB.hh" // USES UserFunctionDB
#include "spatialdata/geocoords/CSCart.hh" // USES CSCart
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

// ----------------------------------------------------------------------
namespace pylith {
    namespace faults {
        class TestFaultCohesiveKin_TriP1;
        class TestFaultCohesiveKin_TriP2;

        class TestFaultCohesiveKin_QuadQ1;
        class TestFaultCohesiveKin_QuadQ2;

    } // namespace faults
} // namespace pylith

class pylith::faults::TestFaultCohesiveKin_TriP1 : public TestFaultCohesiveKin {
    static const char* disp_units(void) {
        return "m";
    }

    static const char* pressure_units(void) {
        return "Pa";
    }

    /// Spatial database user functions for auxiliary subfields.

    static double slip_left_lateral(const double x,
                                    const double y) {
        return 2.0;
    } // slip_left_lateral

    static double slip_opening(const double x,
                               const double y) {
        return 0.5;
    } // slip_opening

    // Solution field at time t.

    static double disp_x(const double x,
                         const double y) {
        return 0.25;
    } // disp_x

    static double disp_y(const double x,
                         const double y) {
        return (x > 0) ? +1.0 : -1.0;
    } // disp_y

    static double lagrange_multiplier_fault_x(const double x,
                                              const double y) {
        return 0;
    } // lagrange_multiplier_fault_x

    static double lagrange_multiplier_fault_y(const double x,
                                              const double y) {
        return 0;
    } // lagrange_multiplier_fault_y

    CPPUNIT_TEST_SUB_SUITE(TestFaultCohesiveKin_TriP1, TestFaultCohesiveKin);
    CPPUNIT_TEST_SUITE_END();

protected:

    void setUp(void) {
        TestFaultCohesiveKin::setUp();
        _data = new TestFaultCohesiveKin_Data();CPPUNIT_ASSERT(_data);
        _data->meshFilename = "data/tri_small.mesh";
        _data->faultLabel = "fault";
        _data->cs = new spatialdata::geocoords::CSCart();CPPUNIT_ASSERT(_data->cs);
        _data->cs->setSpaceDim(2);
        _data->cs->initialize();

        CPPUNIT_ASSERT(_data->normalizer);
        _data->normalizer->lengthScale(1000.0);
        _data->normalizer->timeScale(10.0);
        _data->normalizer->pressureScale(0.1);
        _data->normalizer->densityScale(2.0);

        // eqsrcs

        _data->numAuxSubfields = 1;
        static const char* auxSubfields[1] = { "slip" };
        _data->auxSubfields = const_cast<const char**>(auxSubfields);
        static const pylith::topology::Field::Discretization auxDiscretizations[1] = {
            pylith::topology::Field::Discretization(0, 1), // initial_amplitude
        };
        _data->auxDiscretizations = const_cast<pylith::topology::Field::Discretization*>(auxDiscretizations);

        CPPUNIT_ASSERT(_data->auxDB);
        _data->auxDB->coordsys(*_data->cs);
        _data->auxDB->addValue("slip_opening", slip_opening, disp_units());
        _data->auxDB->addValue("slip_left_lateral", slip_left_lateral, disp_units());

        _data->t = 1.23;
        _data->solnNumSubfields = 2;
        static const pylith::topology::Field::Discretization solnDiscretizations[2] = {
            pylith::topology::Field::Discretization(0, 1), // displacement
            pylith::topology::Field::Discretization(0, 1), // lagrange_multiplier_fault
        };
        _data->solnDiscretizations = const_cast<pylith::topology::Field::Discretization*>(solnDiscretizations);

        CPPUNIT_ASSERT(_data->solnDB);
        _data->solnDB->coordsys(*_data->cs);
        _data->solnDB->addValue("displacement_x", disp_x, disp_units());
        _data->solnDB->addValue("displacement_y", disp_y, disp_units());
        _data->solnDB->addValue("lagrange_multiplier_fault_x", lagrange_multiplier_fault_x, pressure_units());
        _data->solnDB->addValue("lagrange_multiplier_fault_y", lagrange_multiplier_fault_y, pressure_units());

    } // setUp

}; // class TestFaultCohesiveKin_TriP1

CPPUNIT_TEST_SUITE_REGISTRATION(pylith::faults::TestFaultCohesiveKin_TriP1);

// End of file
