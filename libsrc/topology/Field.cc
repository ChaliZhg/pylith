// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "Field.hh" // implementation of class methods

#include "Mesh.hh" // HASA Mesh
#include "pylith/utils/array.hh" // USES double_array

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::Field::Field(const Mesh& mesh) :
  _mesh(mesh)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::topology::Field::~Field(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get spatial dimension of domain.
int
pylith::topology::Field::spaceDim(void) const
{ // spaceDim
  const spatialdata::geocoords::CoordSys* cs = _mesh.coordsys();
  return (0 != cs) ? cs->spaceDim() : 0;
} // spaceDim

// ----------------------------------------------------------------------
// Create seive section.
void
pylith::topology::Field::newSection(void)
{ // newSection
  _section = new MeshRealSection(_mesh.comm(), _mesh.debug());  
} // newSection

// ----------------------------------------------------------------------
// Create sieve section and set chart and fiber dimesion.
void
pylith::topology::Field::newSection(
			  const ALE::Obj<SieveMesh::label_sequence>& points,
			  const int fiberDim)
{ // newSection
  if (fiberDim < 0) {
    std::ostringstream msg;
    msg
      << "Fiber dimension (" << fiberDim << ") for field '" << _name
      << "' must be nonnegative.";
    throw std::runtime_error(msg.str());
  } // if

  _section = new MeshRealSection(_mesh.comm(), _mesh.debug());

  const SieveMesh::point_type pointMin = 
    *std::min_element(points->begin(), points->end());
  const SieveMesh::point_type pointMax = 
    *std::max_element(points->begin(), points->end());
  _section->setChart(MeshRealSection::chart_type(pointMin, pointMax+1));
  _section->setFiberDimension(points, fiberDim);  
} // newSection

// ----------------------------------------------------------------------
// Create sieve section and set chart and fiber dimesion.
void
pylith::topology::Field::newSection(const DomainEnum domain,
				    const int fiberDim)
{ // newSection
  const ALE::Obj<SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  ALE::Obj<SieveMesh::label_sequence> points;
  if (VERTICES_FIELD == domain)
    points = sieveMesh->depthStratum(0);
  else if (CELLS_FIELD == domain)
    points = sieveMesh->heightStratum(1);
  else {
    std::cerr << "Unknown value for DomainEnum: " << domain << std::endl;
    assert(0);
  } // else

  newSection(points, fiberDim);
} // newSection

// ----------------------------------------------------------------------
// Create section given chart.
void
pylith::topology::Field::newSection(const MeshRealSection::chart_type& chart,
				    const int fiberDim)
{ // newSection
  if (_section.isNull())
    newSection();

  _section->setChart(chart);

  const MeshRealSection::chart_type::const_iterator chartEnd = chart.end();
  for (MeshRealSection::chart_type::const_iterator c_iter = chart.begin();
       c_iter != chartEnd;
       ++c_iter)
    _section->setFiberDimension(*c_iter, fiberDim);
  allocate();
} // newSection

// ----------------------------------------------------------------------
// Create section with same layout as another section.
void
pylith::topology::Field::newSection(const Field& src)
{ // newSection
  _vecFieldType = src._vecFieldType;

  const ALE::Obj<MeshRealSection>& srcSection = src.section();
  if (!srcSection.isNull() && _section.isNull())
    newSection();

  if (!_section.isNull()) {
    _section->setAtlas(srcSection->getAtlas());
    _section->allocateStorage();
    _section->setBC(srcSection->getBC());
  } // if
} // newSection

// ----------------------------------------------------------------------
// Clear variables associated with section.
void
pylith::topology::Field::clear(void)
{ // clear
  if (!_section.isNull())
    _section->clear();

  _scale = 1.0;
  _vecFieldType = OTHER;
  _dimensionsOkay = false;
} // clear

// ----------------------------------------------------------------------
// Allocate Sieve section.
void
pylith::topology::Field::allocate(void)
{ // allocate
  assert(!_section.isNull());

  const ALE::Obj<SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());
  sieveMesh->allocate(_section);
} // allocate

// ----------------------------------------------------------------------
// Zero section values.
void
pylith::topology::Field::zero(void)
{ // zero
  if (!_section.isNull())
    _section->zero();
} // zero

// ----------------------------------------------------------------------
// Complete section by assembling across processors.
void
pylith::topology::Field::complete(void)
{ // complete
  const ALE::Obj<SieveMesh>& sieveMesh = _mesh.sieveMesh();
  assert(!sieveMesh.isNull());

  if (!_section.isNull())
    ALE::Completion::completeSectionAdd(sieveMesh->getSendOverlap(),
					sieveMesh->getRecvOverlap(), 
					_section, _section);
} // complete

// ----------------------------------------------------------------------
// Copy field values and metadata.
void
pylith::topology::Field::copy(const Field& field)
{ // copy
  // Check compatibility of sections
  const int srcSize = (!field._section.isNull()) ? field._section->size() : 0;
  const int dstSize = (!_section.isNull()) ? _section->size() : 0;
  if (field.spaceDim() != spaceDim() ||
      field._vecFieldType != _vecFieldType ||
      field._scale != _scale ||
      srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot copy values from section '" << field._name 
	<< "' to section '" << _name << "'. Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    space dim: " << field.spaceDim() << "\n"
	<< "    vector field type: " << field._vecFieldType << "\n"
	<< "    scale: " << field._scale << "\n"
	<< "    size: " << srcSize
	<< "  Destination section:\n"
	<< "    space dim: " << spaceDim() << "\n"
	<< "    vector field type: " << _vecFieldType << "\n"
	<< "    scale: " << _scale << "\n"
	<< "    size: " << dstSize;
    throw std::runtime_error(msg.str());
  } // if
  assert( (_section.isNull() && field._section.isNull()) ||
	  (!_section.isNull() && !field._section.isNull()) );

  if (!_section.isNull()) {
    // Copy values from field
    const MeshRealSection::chart_type& chart = _section->getChart();
    const MeshRealSection::chart_type::const_iterator chartEnd = chart.end();

    for (MeshRealSection::chart_type::const_iterator c_iter = chart.begin();
	 c_iter != chartEnd;
	 ++c_iter) {
      assert(field._section->getFiberDimension(*c_iter) ==
	     _section->getFiberDimension(*c_iter));
      _section->updatePoint(*c_iter, field._section->restrictPoint(*c_iter));
    } // for
  } // if
} // copy

// ----------------------------------------------------------------------
// Add two fields, storing the result in one of the fields.
void
pylith::topology::Field::operator+=(const Field& field)
{ // operator+=
  // Check compatibility of sections
  const int srcSize = (!field._section.isNull()) ? field._section->size() : 0;
  const int dstSize = (!_section.isNull()) ? _section->size() : 0;
  if (field.spaceDim() != spaceDim() ||
      field._vecFieldType != _vecFieldType ||
      field._scale != _scale ||
      srcSize != dstSize) {
    std::ostringstream msg;

    msg << "Cannot add values from section '" << field._name 
	<< "' to section '" << _name << "'. Sections are incompatible.\n"
	<< "  Source section:\n"
	<< "    space dim: " << field.spaceDim() << "\n"
	<< "    vector field type: " << field._vecFieldType << "\n"
	<< "    scale: " << field._scale << "\n"
	<< "    size: " << srcSize
	<< "  Destination section:\n"
	<< "    space dim: " << spaceDim() << "\n"
	<< "    vector field type: " << _vecFieldType << "\n"
	<< "    scale: " << _scale << "\n"
	<< "    size: " << dstSize;
    throw std::runtime_error(msg.str());
  } // if
  assert( (_section.isNull() && field._section.isNull()) ||
	  (!_section.isNull() && !field._section.isNull()) );

  if (!_section.isNull()) {
    // Add values from field
    const MeshRealSection::chart_type& chart = _section->getChart();
    const MeshRealSection::chart_type::const_iterator chartEnd = chart.end();

    // Assume fiber dimension is uniform
    const int fiberDim = _section->getFiberDimension(*chart.begin());
    double_array values(fiberDim);

    for (MeshRealSection::chart_type::const_iterator c_iter = chart.begin();
	 c_iter != chartEnd;
	 ++c_iter) {
      assert(fiberDim == field._section->getFiberDimension(*c_iter));
      assert(fiberDim == _section->getFiberDimension(*c_iter));
      field._section->restrictPoint(*c_iter, &values[0], values.size());
      _section->updateAddPoint(*c_iter, &values[0]);
    } // for
  } // if
} // operator+=

// ----------------------------------------------------------------------
// Dimensionalize field.
void
pylith::topology::Field::dimensionalize(void)
{ // dimensionalize
  if (!_dimensionsOkay) {
    std::ostringstream msg;
    msg << "Cannot dimensionalize field '" << _name << "' because the flag "
	<< "has been set to keep field nondimensional.";
    throw std::runtime_error(msg.str());
  } // if

  if (!_section.isNull()) {
    const MeshRealSection::chart_type& chart = _section->getChart();
    const MeshRealSection::chart_type::const_iterator chartEnd = chart.end();

    // Assume fiber dimension is uniform
    const int fiberDim = _section->getFiberDimension(*chart.begin());
    double_array values(fiberDim);

    spatialdata::units::Nondimensional normalizer;

    for (MeshRealSection::chart_type::const_iterator c_iter = chart.begin();
	 c_iter != chartEnd;
	 ++c_iter) {
      assert(fiberDim == _section->getFiberDimension(*c_iter));
      
      _section->restrictPoint(*c_iter, &values[0], values.size());
      normalizer.dimensionalize(&values[0], values.size(), _scale);
      _section->updatePoint(*c_iter, &values[0]);
    } // for
  } // if
} // dimensionalize

// ----------------------------------------------------------------------
// Print field to standard out.
void
pylith::topology::Field::view(const char* label)
{ // view
  std::string vecFieldString;
  switch(_vecFieldType)
    { // switch
    case SCALAR:
      vecFieldString = "scalar";
      break;
    case VECTOR:
      vecFieldString = "vector";
      break;
    case TENSOR:
      vecFieldString = "tensor";
      break;
    case OTHER:
      vecFieldString = "other";
      break;
    case MULTI_SCALAR:
      vecFieldString = "multiple scalars";
      break;
    case MULTI_VECTOR:
      vecFieldString = "multiple vectors";
      break;
    case MULTI_TENSOR:
      vecFieldString = "multiple tensors";
      break;
    case MULTI_OTHER:
      vecFieldString = "multiple other values";
      break;
    default :
      std::cerr << "Unknown vector field value '" << _vecFieldType
		<< "'." << std::endl;
      assert(0);
    } // switch

  std::cout << "Viewing field '" << _name << "' "<< label << ".\n"
	    << "  vector field type: " << vecFieldString << "\n"
	    << "  scale: " << _scale << "\n"
	    << "  dimensionalize flag: " << _dimensionsOkay << std::endl;
  if (!_section.isNull())
    _section->view(label);
} // view


// End of file 
