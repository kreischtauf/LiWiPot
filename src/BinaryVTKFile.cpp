/***************************************************************************
                          BinaryVTKFile.cpp
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

#include <fstream>
#include <boost/regex.hpp>

#include "BinaryVTKFile.hh"

extern std::ofstream dbg;

void BinaryVtkFile::writeFile(const std::string & baseFileName)
{
    const std::string nl("\n");
    const boost::regex pattern("\n");
    const std::string replaceFormat = nl + std::string(indent_l3);
    std::string pointData = concatenateFields();
    char buffer[50];
    sprintf(buffer, "%s_node%05d.vtr", baseFileName.c_str(), Ippl::myNode());
    std::string headerFileName = baseFileName + std::string(".pvtr");
    std::ofstream out(buffer);

    out << indent_l0 << _header
        << indent_l1 << "<RectilinearGrid WholeExtent=" << _localGridExtent << ">" << "\n"
        << indent_l2 << "<Piece Extent=" << _localGridExtent << ">\n"
        << indent_l3 << boost::regex_replace(pointData, pattern, replaceFormat) << "\n"
        << indent_l3 << boost::regex_replace(_coordinates, pattern, replaceFormat) << "\n"
        << indent_l2 << "</Piece>\n"
        << indent_l1 << "</RectilinearGrid>\n"
        << indent_l1 << "<AppendedData encoding=\"raw\">\n"
        << indent_l2 << "_";
    out.write(&(_data[0]), _data.size());
    out << nl
        << indent_l1 << "</AppendedData>\n"
        << indent_l0 << "</VTKFile>" << std::endl;
    out.close();

    writeHeaderFile(baseFileName);
}

void BinaryVtkFile::writeHeaderFile(const std::string & baseFileName)
{
    if (Ippl::myNode() == 0) {
        std::string fileName = baseFileName + std::string(".pvtr");
        std::ofstream out(fileName.c_str());

        const std::string nl("\n");
        const boost::regex pattern("\n");
        const std::string replaceFormat_l2 = nl + std::string(indent_l2);
        std::string fieldHeaders = concatenateFieldHeaders();
        std::string pieceExtents = concatenatePieceExtents(baseFileName);

        out << indent_l0 << _headerFileHeader
            << indent_l1 << "<PRectilinearGrid WholeExtent=" << _globalGridExtent << " GhostLevel=\"1\">\n"
            << indent_l2 << boost::regex_replace(fieldHeaders, pattern, replaceFormat_l2) << nl
            << indent_l2 << "<PCoordinates>\n"
            << indent_l3 << "<PDataArray type=\"Float32\" name=\"X_COORDINATES\" NumberOfComponents=\"1\"/>\n"
            << indent_l3 << "<PDataArray type=\"Float32\" name=\"Y_COORDINATES\" NumberOfComponents=\"1\"/>\n"
            << indent_l3 << "<PDataArray type=\"Float32\" name=\"Z_COORDINATES\" NumberOfComponents=\"1\"/>\n"
            << indent_l2 << "</PCoordinates>\n"
            << indent_l2 << boost::regex_replace(pieceExtents, pattern, replaceFormat_l2) << nl
            << indent_l1 << "</PRectilinearGrid>\n"
            << indent_l0 << "</VTKFile>" << std::endl;
    }
}

void BinaryVtkFile::addCoordinates(const Vector_t & dx)
{
    NDIndex<DIM> localDomain = _localDomains[Ippl::myNode()];
    const std::string nl("\n");
    boost::regex pattern("\n");
    const std::string replaceFormat = nl + std::string(indent_l1);
    std::stringstream vtkPart;

    std::string xCoordinates = addXCoordinates(dx);
    std::string yCoordinates = addYCoordinates(dx);
    std::string zCoordinates = addZCoordinates(dx);

    vtkPart << indent_l0 << "<Coordinates>\n"
            << indent_l1 << boost::regex_replace(xCoordinates, pattern, replaceFormat) << nl
            << indent_l1 << boost::regex_replace(yCoordinates, pattern, replaceFormat) << nl
            << indent_l1 << boost::regex_replace(zCoordinates, pattern, replaceFormat) << nl
            << indent_l0 << "</Coordinates>";

    _coordinates = vtkPart.str();
}

std::string BinaryVtkFile::concatenateScalarFieldNames() const
{
    std::stringstream scalarFieldNames, fieldNames;
    const boost::regex pattern("(.*),(\\d).*");
    boost::cmatch what;
    size_t counter = 0;

    for (std::map<std::string, std::string>::const_iterator fieldIterator = _fields.begin();
         fieldIterator != _fields.end();
         ++ fieldIterator) {

        std::string fieldNameAndComponents = (*fieldIterator).first;
        if (boost::regex_match(fieldNameAndComponents.c_str(), what, pattern)) {
            if (atoi(what[2].first) == 1) {
                scalarFieldNames << what[1] << ", ";
                ++ counter;
            }
        }
    }

    if (counter) {
        std::string tmpScalarFieldNames = scalarFieldNames.str();
        return std::string("\"") + tmpScalarFieldNames.substr(0, tmpScalarFieldNames.size() - 2) + std::string("\" ");
    } else {
        return "\"\"";
    }
}

std::string BinaryVtkFile::concatenateVectorFieldNames() const
{
    std::stringstream vectorFieldNames, fieldNames;
    const boost::regex pattern("(.*),(\\d).*");
    boost::cmatch what;
    size_t counter = 0;

    for (std::map<std::string, std::string>::const_iterator fieldIterator = _fields.begin();
         fieldIterator != _fields.end();
         ++ fieldIterator) {

        std::string fieldNameAndComponents = (*fieldIterator).first;
        if (boost::regex_match(fieldNameAndComponents.c_str(), what, pattern)) {
            if (atoi(what[2].first) == DIM) {
                vectorFieldNames << what[1] << ", ";
                ++ counter;
            }
        }
    }

    if (counter) {
        std::string tmpVectorFieldNames = vectorFieldNames.str();
        return std::string("\"") + tmpVectorFieldNames.substr(0, tmpVectorFieldNames.size() - 2) + std::string("\" ");
    } else {
        return "\"\"";
    }
}

std::string BinaryVtkFile::concatenateFields() const
{
    const std::string nl("\n");
    const boost::regex pattern("\n");
    const std::string replaceFormat = nl + std::string(indent_l3);
    std::string fieldNames = concatenateFieldNames();
    std::stringstream fieldOutput;

    fieldOutput << indent_l0 << "<PointData " << fieldNames << ">\n";
    for (std::map<std::string, std::string>::const_iterator fieldIterator = _fields.begin();
         fieldIterator != _fields.end();
         ++ fieldIterator) {
        std::string field = (*fieldIterator).second;
        std::string indentedField = boost::regex_replace(field, pattern, replaceFormat);
        fieldOutput << indent_l1 << indentedField << nl;
    }
    fieldOutput << indent_l0 << "</PointData>";

    return fieldOutput.str();
}

std::string BinaryVtkFile::concatenateFieldHeaders() const
{
    const std::string nl("\n");
    const boost::regex pattern("\n");
    const std::string replaceFormat = nl + std::string(indent_l1);
    std::stringstream fieldOutput;

    fieldOutput << indent_l0 << "<PPointData " << concatenateFieldNames() << ">\n";
    for (std::map<std::string, std::string>::const_iterator fieldIterator = _fieldHeaders.begin();
         fieldIterator != _fieldHeaders.end();
         ++ fieldIterator) {
        std::string fieldHeader = (*fieldIterator).second;
        fieldOutput << indent_l1 << boost::regex_replace(fieldHeader, pattern, replaceFormat) << nl;
    }
    fieldOutput << indent_l0 << "</PPointData>";

    return fieldOutput.str();
}

std::string BinaryVtkFile::concatenatePieceExtents(const std::string & pathFile) const
{
    char buffer[50];
    std::stringstream vtkout;

    size_t slashpos = pathFile.find_last_of("/");
    std::string baseFileName = pathFile.substr(slashpos+1);

    for (int l = 0; l < Ippl::getNodes(); ++ l) {
        NDIndex<DIM> domain = _localDomains[l];
        sprintf(buffer, "%s_node%05d.vtr", baseFileName.c_str(), l);
        std::string fileName(buffer);

        vtkout << indent_l0 << "<Piece Extent=\""
               << domain[0].first() << " " << domain[0].last() << " "
               << domain[1].first() << " " << domain[1].last() << " "
               << 1 << " " << 1 << "\" Source=\"" << fileName << "\"/>\n";
    }
    std::string returnValue = vtkout.str();
    return returnValue.substr(0, returnValue.size()-1);
}
