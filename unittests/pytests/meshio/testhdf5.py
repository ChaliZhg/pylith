#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/meshio/testhdf5.py

## @brief Python application for testing HDF5 related meshio code.

from pyre.applications.Script import Script

import unittest

class TestApp(Script):
  """
  Test application.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="testapp"):
    """
    Constructor.
    """
    Script.__init__(self, name)
    return


  def main(self):
    """
    Run the application.
    """
    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.options = [("malloc_dump", "true")]
    petsc.initialize()

    unittest.TextTestRunner(verbosity=2).run(self._suite())

    petsc.finalize()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _suite(self):
    """
    Setup the test suite.
    """

    suite = unittest.TestSuite()

    from TestDataWriterHDF5 import TestDataWriterHDF5Mesh
    suite.addTest(unittest.makeSuite(TestDataWriterHDF5Mesh))

    from TestDataWriterHDF5 import TestDataWriterHDF5SubMesh
    suite.addTest(unittest.makeSuite(TestDataWriterHDF5SubMesh))

    from TestDataWriterHDF5 import TestDataWriterHDF5SubSubMesh
    suite.addTest(unittest.makeSuite(TestDataWriterHDF5SubSubMesh))

    from TestDataWriterHDF5Ext import TestDataWriterHDF5ExtMesh
    suite.addTest(unittest.makeSuite(TestDataWriterHDF5ExtMesh))

    from TestDataWriterHDF5Ext import TestDataWriterHDF5ExtSubMesh
    suite.addTest(unittest.makeSuite(TestDataWriterHDF5ExtSubMesh))

    from TestDataWriterHDF5Ext import TestDataWriterHDF5ExtSubSubMesh
    suite.addTest(unittest.makeSuite(TestDataWriterHDF5ExtSubSubMesh))

    from TestXdmf import TestXdmf
    suite.addTest(unittest.makeSuite(TestXdmf))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
