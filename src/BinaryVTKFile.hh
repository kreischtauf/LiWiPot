/***************************************************************************
                           BinaryVTKFile.hh
                         -------------------
    begin                : Sun Dec 4 2011
    copyright            : (C) 2011 by Christof Kraus
    email                : christof.kraus-csrst@my.mail.de
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef BINARYVTKFILE_HH
#define BINARYVTKFILE_HH

#include <cstdint>
#include "defs.hh"
#include "utils.hh"

class BinaryVtkFile {
public:
    BinaryVtkFile();

    template<class FT,class C>
    void addGridExtent(const Field<FT,DIM,Mesh_t,C> & exampleField);

    template<class C>
    void addScalarField(const Field<Scalar_t,DIM,Mesh_t,C> & scalarField, const std::string & fieldName);
    template<class C>
    void addVectorField(const Field<Vector_t,DIM,Mesh_t,C> & vectorField, const std::string & fieldName);

    void writeFile(const std::string & baseFileName);

private:
    void addHeader();
    void addHeaderFileHeader();
    void addCoordinates(const Vector_t & dx);
    std::string addXCoordinates(const Vector_t & dx);
    std::string addYCoordinates(const Vector_t & dx);
    std::string addZCoordinates(const Vector_t & dx);
    void addGridExtent(const FieldLayout<DIM> & layout, const Vector_t & dx);
    void addScalarFieldHeader(const std::string & fieldName);
    void addVectorFieldHeader(const std::string & fieldName);
    void setLocalDomains(const FieldLayout<DIM> & layout);
    template<class C>
    void collectScalarData(const Field<Scalar_t,DIM,Mesh_t,C> & scalarField);
    template<class C>
    void collectVectorData(const Field<Vector_t,DIM,Mesh_t,C> & vectorField);
    std::string concatenateFieldNames() const;
    std::string concatenateVectorFieldNames() const;
    std::string concatenateScalarFieldNames() const;
    std::string concatenateFields() const;
    std::string concatenateFieldHeaders() const;
    std::string concatenatePieceExtents(const std::string & baseFileName) const;
    void writeHeaderFile(const std::string & baseFileName);

    const double _cutRelativeValue;
    bool _localDomainsSet;
    bool _gridExtentAdded;
    bool _oneFieldAdded;
    bool _fileWritten;

    std::string _header;
    std::string _headerFileHeader;
    std::string _globalGridExtent;
    std::string _localGridExtent;
    std::string _coordinates;
    std::vector<char> _data;

    std::map<std::string, std::string> _fields;
    std::map<std::string, std::string> _fieldHeaders;
    std::vector<NDIndex<DIM> > _localDomains;

    Vector_t _origin;
};

inline
BinaryVtkFile::BinaryVtkFile():
    _cutRelativeValue(1e-6),
    _localDomainsSet(false),
    _gridExtentAdded(false),
    _oneFieldAdded(false),
    _fileWritten(false),
    _header(""),
    _headerFileHeader(""),
    _globalGridExtent(""),
    _localGridExtent(""),
    _coordinates(""),
    _data(),
    _fields(),
    _fieldHeaders(),
    _localDomains()
{
    addHeader();
}

template<class FT,class C>
void BinaryVtkFile::addGridExtent(const Field<FT,DIM,Mesh_t,C> & exampleField)
{
    FieldLayout<DIM> & layout = exampleField.getLayout();
    Vector_t dx(exampleField.get_mesh().get_meshSpacing(0),
                exampleField.get_mesh().get_meshSpacing(1));
    addGridExtent(layout, dx);
}

inline
void BinaryVtkFile::addGridExtent(const FieldLayout<DIM> & layout,
                                  const Vector_t & dx)
{
    if (!_gridExtentAdded) {
        if (!_localDomainsSet) {
            setLocalDomains(layout);
        }
        NDIndex<DIM> gDom = layout.getDomain();
        NDIndex<DIM> localDomain = _localDomains[Ippl::myNode()];
        std::stringstream vtkPart;

        vtkPart << "\""
                << gDom[0].first() << " " << gDom[0].last() << " "
                << gDom[1].first() << " " << gDom[1].last() << " "
                << 1 << " " << 1
                << "\"";
        _globalGridExtent = vtkPart.str();
        vtkPart.str(std::string(""));
        vtkPart << "\""
                << localDomain[0].first() << " " << localDomain[0].last() << " "
                << localDomain[1].first() << " " << localDomain[1].last() << " "
                << 1 << " " << 1
                << "\"";
        _localGridExtent = vtkPart.str();

        _gridExtentAdded = true;

        addCoordinates(dx);
    }
}

inline
void BinaryVtkFile::setLocalDomains(const FieldLayout<DIM> & layout)
{
    Utils::getLocalDomains(layout, _localDomains);
    Utils::addGostCellToLocalDomains(layout, _localDomains);

    _localDomainsSet = true;
}

inline
void BinaryVtkFile::addScalarFieldHeader(const std::string & fieldName)
{
    std::stringstream vtkPart;

    vtkPart << indent_l0
            << "<PDataArray type=\"Float32\" Name=\"" << fieldName << "\" NumberOfComponents=\"1\"/>";

    _fieldHeaders.insert(std::pair<std::string, std::string>(fieldName, vtkPart.str()));
}

inline
void BinaryVtkFile::addVectorFieldHeader(const std::string & fieldName)
{
    std::stringstream vtkPart;

    vtkPart << indent_l0
            << "<PDataArray type=\"Float32\" Name=\"" << fieldName << "\" NumberOfComponents=\"3\"/>";

    _fieldHeaders.insert(std::pair<std::string, std::string>(fieldName, vtkPart.str()));
}

inline
void BinaryVtkFile::addHeader()
{
    std::stringstream vtkPart;
    short EndianTest_s;
    unsigned char *EndianTest = reinterpret_cast<unsigned char*>(&EndianTest_s);
    EndianTest[0] = 1;
    EndianTest[1] = 0;
    std::string endianness = (EndianTest_s == 1? "\"LittleEndian\"":"\"BigEndian\"");

    vtkPart << "<?xml version=\"1.0\"?>" << std::endl;
    vtkPart << indent_l0
            << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=" << endianness << ">\n";

    _header = vtkPart.str();

    addHeaderFileHeader();
}

inline
void BinaryVtkFile::addHeaderFileHeader()
{
    std::stringstream vtkPart;
    short EndianTest_s;
    unsigned char *EndianTest = reinterpret_cast<unsigned char*>(&EndianTest_s);
    EndianTest[0] = 1;
    EndianTest[1] = 0;
    std::string endianness = (EndianTest_s == 1? "\"LittleEndian\"":"\"BigEndian\"");

    vtkPart << "<?xml version=\"1.0\"?>" << std::endl;
    vtkPart << indent_l0
            << "<VTKFile type=\"PRectilinearGrid\" version=\"0.1\" byte_order=" << endianness << ">\n";

    _headerFileHeader = vtkPart.str();
}

inline
std::string BinaryVtkFile::addXCoordinates(const Vector_t & dx)
{
    NDIndex<DIM> localDomain = _localDomains[Ippl::myNode()];
    std::stringstream coordinatesStream;
    size_t offset = _data.size();
    const char* buffer;
    uint32_t datasize = sizeof(float) * (localDomain[0].last() - localDomain[0].first() + 1);

    coordinatesStream << "<DataArray type=\"Float32\" name=\"X_COORDINATES\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\"/>";

    _data.reserve(_data.size() + sizeof(uint32_t) + datasize);
    buffer = reinterpret_cast<const char*>(&datasize);
    _data.insert(_data.end(), buffer, buffer + sizeof(uint32_t));
    for (int i = localDomain[0].first(); i <= localDomain[0].last(); ++ i) {
        float tmp = i * dx[0] + _origin[0];
        buffer = reinterpret_cast<const char* >(&tmp);
        _data.insert(_data.end(), buffer, buffer + sizeof(float));
    }

    return coordinatesStream.str();
}

inline
std::string BinaryVtkFile::addYCoordinates(const Vector_t & dx)
{
    NDIndex<DIM> localDomain = _localDomains[Ippl::myNode()];
    std::stringstream coordinatesStream;
    size_t offset = _data.size();
    const char* buffer;
    uint32_t datasize = sizeof(float) * (localDomain[1].last() - localDomain[1].first() + 1);

    coordinatesStream << "<DataArray type=\"Float32\" name=\"Y_COORDINATES\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\"/>";

    _data.reserve(_data.size() + sizeof(uint32_t) + datasize);
    buffer = reinterpret_cast<const char* >(&datasize);
    _data.insert(_data.end(), buffer, buffer + sizeof(uint32_t));
    for (int j = localDomain[1].first(); j <= localDomain[1].last(); ++ j) {
        float tmp = j * dx[1] + _origin[1];
        buffer = reinterpret_cast<const char* >(&tmp);
        _data.insert(_data.end(), buffer, buffer + sizeof(float));
    }

    return coordinatesStream.str();
}

inline
std::string BinaryVtkFile::addZCoordinates(const Vector_t & dx)
{
    std::stringstream coordinatesStream;
    size_t offset = _data.size();
    uint32_t datasize = sizeof(float);
    const char* buffer;
    float tmp = 0.0;

    coordinatesStream << "<DataArray type=\"Float32\" name=\"Z_COORDINATES\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\"/>";

    _data.reserve(_data.size() + sizeof(uint32_t) + datasize);
    buffer = reinterpret_cast<const char* >(&datasize);
    _data.insert(_data.end(), buffer, buffer + sizeof(uint32_t));
    buffer = reinterpret_cast<const char* >(&tmp);
    _data.insert(_data.end(), buffer, buffer + sizeof(float));

    return coordinatesStream.str();
}

inline
std::string BinaryVtkFile::concatenateFieldNames() const
{
    std::stringstream fieldNames;

    fieldNames << "Scalars=" << concatenateScalarFieldNames()
               << " Vectors=" << concatenateVectorFieldNames();
    return fieldNames.str();
}

template<class C>
void BinaryVtkFile::addScalarField(const Field<Scalar_t,DIM,Mesh_t,C> & scalarField, const std::string & fieldName)
{
    _origin = scalarField.get_mesh().get_origin();
    if (!_gridExtentAdded) {
        addGridExtent(scalarField);
    }
    std::stringstream vtkPart;
    std::stringstream fieldNameAndComponents;
    size_t offset = _data.size();
    fieldNameAndComponents << fieldName << ",1";
    vtkPart << indent_l0 << "<DataArray type=\"Float32\" Name=\"" << fieldName << "\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << offset << "\"/>";
    collectScalarData(scalarField);


    _fields.insert(std::pair<std::string, std::string>(fieldNameAndComponents.str(), vtkPart.str()));
    addScalarFieldHeader(fieldName);
    _oneFieldAdded = true;
}

template<class C>
void BinaryVtkFile::addVectorField(const Field<Vector_t,DIM,Mesh_t,C> & vectorField, const std::string & fieldName)
{
    _origin = vectorField.get_mesh().get_origin();
    if (!_gridExtentAdded) {
        addGridExtent(vectorField);
    }
    std::stringstream vtkPart;
    std::stringstream fieldNameAndComponents;
    size_t offset = _data.size();
    fieldNameAndComponents << fieldName << "," << DIM;
    vtkPart << indent_l0 << "<DataArray type=\"Float32\" Name=\"" << fieldName << "\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"/>";
    collectVectorData(vectorField);

    _fields.insert(std::pair<std::string, std::string>(fieldNameAndComponents.str(), vtkPart.str()));
    addVectorFieldHeader(fieldName);
    _oneFieldAdded = true;
}

template<class C>
void BinaryVtkFile::collectScalarData(const Field<Scalar_t,DIM,Mesh_t,C> & scalarField)
{
    std::stringstream vtkPart;
    NDIndex<DIM> elem;
    const NDIndex<DIM> & localDomain = _localDomains[Ippl::myNode()];
    const char* buffer;
    uint32_t datasize = sizeof(float) *
        (localDomain[0].last() - localDomain[0].first() + 1) *
        (localDomain[1].last() - localDomain[1].first() + 1);

    _data.reserve(_data.size() + sizeof(uint32_t) + datasize);
    buffer = reinterpret_cast<const char* >(&datasize);
    _data.insert(_data.end(), buffer, buffer + sizeof(uint32_t));
    for (int j = localDomain[1].first(); j <= localDomain[1].last(); ++ j) {
        elem[1]=Index(j,j);
        for (int i = localDomain[0].first(); i <= localDomain[0].last(); ++ i) {
            elem[0]=Index(i,i);
            float tmp = scalarField.localElement(elem);
            buffer = reinterpret_cast<const char* >(&tmp);
            _data.insert(_data.end(), buffer, buffer + sizeof(float));
        }
    }
}

template<class C>
void BinaryVtkFile::collectVectorData(const Field<Vector_t,DIM,Mesh_t,C> & vectorField)
{
    std::stringstream vtkPart;
    NDIndex<DIM> elem;
    const NDIndex<DIM> & localDomain = _localDomains[Ippl::myNode()];
    const char* buffer;
    uint32_t datasize = sizeof(float) * 3 *
        (localDomain[0].last() - localDomain[0].first() + 1) *
        (localDomain[1].last() - localDomain[1].first() + 1);

    _data.reserve(_data.size() + sizeof(uint32_t) + datasize);
    buffer = reinterpret_cast<const char* >(&datasize);
    _data.insert(_data.end(), buffer, buffer + sizeof(uint32_t));
    for (int j = localDomain[1].first(); j <= localDomain[1].last(); ++ j) {
        elem[1]=Index(j,j);
        for (int i = localDomain[0].first(); i <= localDomain[0].last(); ++ i) {
            elem[0]=Index(i,i);
            Vector_t tmpVector = vectorField.localElement(elem);
            for (int l = 0; l < DIM; ++ l) {
                float tmp = tmpVector(l);
                buffer = reinterpret_cast<const char* >(&tmp);
                _data.insert(_data.end(), buffer, buffer + sizeof(float));
            }
            float tmp = 0.0;
            buffer = reinterpret_cast<const char* >(&tmp);
            _data.insert(_data.end(), buffer, buffer + sizeof(float));
        }
    }
}

#endif
