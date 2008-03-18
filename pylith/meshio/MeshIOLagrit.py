#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file pyre/meshio/MeshIOLagrit.py
##
## @brief Python object for reading/writing finite-element mesh from
## LaGriT.
##
## Factory: mesh_io

from MeshIO import MeshIO

# MeshIOLagrit class
class MeshIOLagrit(MeshIO):
  """
  Python object for reading/writing finite-element mesh from LaGriT.

  Factory: mesh_io
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(MeshIO.Inventory):
    """
    Python object for managing MeshIOLagrit facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing MeshIOLagrit facilities and properties.
    ##
    ## \b Properties
    ## @li \b filename_gmv Name of mesh GMV file.
    ## @li \b filename_pset Name of mesh PSET file.
    ## @li \b flip_endian Flip endian type when reading/writing binary files.
    ## @li \b io_int32 PSET files use 64-bit integers.
    ## @li \b record_header_32bit Fortran record header is 32-bit.
    ##
    ## \b Facilities
    ## @li coordsys Coordinate system associated with mesh.

    import pyre.inventory

    filenameGmv = pyre.inventory.str("filename_gmv", default="mesh.gmv")
    filenameGmv.meta['tip'] = "Name of mesh GMV file."

    filenamePset = pyre.inventory.str("filename_pset", default="mesh.pset")
    filenamePset.meta['tip'] = "Name of mesh PSET file."

    flipEndian = pyre.inventory.bool("flip_endian", default=False)
    flipEndian.meta['tip'] = "Flip endian type when reading/writing binary " \
                             "files."

    ioInt32 = pyre.inventory.bool("io_int32", default=True)
    ioInt32.meta['tip'] = "PSTE files use 32-bit integers."

    isRecordHeader32Bit = pyre.inventory.bool("record_header_32bit",
                                              default=True)
    isRecordHeader32Bit.meta['tip'] = "Fortran record header is 32-bit."

    from spatialdata.geocoords.CSCart import CSCart
    coordsys = pyre.inventory.facility("coordsys", family="coordsys",
                                       factory=CSCart)
    coordsys.meta['tip'] = "Coordinate system associated with mesh."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="meshiolagrit"):
    """
    Constructor.
    """
    MeshIO.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    MeshIO._configure(self)
    self.filenameGmv = self.inventory.filenameGmv
    self.filenamePset = self.inventory.filenamePset
    self.coordsys = self.inventory.coordsys
    self.flipEndian = self.inventory.flipEndian
    self.ioInt32 = self.inventory.ioInt32
    self.isRecordHeader32Bit = self.inventory.isRecordHeader32Bit
    return


  def _sync(self):
    """
    Force synchronization between Python and C++.
    """
    if None == self.cppHandle:
      import pylith.meshio.meshio as bindings
      self.cppHandle = bindings.MeshIOLagrit()
    
    MeshIO._sync(self)
    self.cppHandle.filenameGmv = self.filenameGmv
    self.cppHandle.filenamePset = self.filenamePset
    self.cppHandle.flipEndian = self.flipEndian
    self.cppHandle.ioInt32 = self.ioInt32
    self.cppHandle.isRecordHeader32Bit = self.isRecordHeader32Bit
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def mesh_io():
  """
  Factory associated with MeshIOLagrit.
  """
  return MeshIOLagrit()


# End of file 
