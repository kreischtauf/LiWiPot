/***************************************************************************
                           PlotBinaryVTK.hh
                         -------------------
    begin                : Tue Jun 23 2009
    copyright            : (C) 2009 by Christof Kraus
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

#ifndef PLOTBINARYVTKCOMMAND_HH
#define PLOTBINARYVTKCOMMAND_HH

#include <vector>
#include <fstream>

#include "defs.hh"

class PlotBinaryVTKCommand {
public:
    PlotBinaryVTKCommand(const VField_t & EFD,
                         const std::string & baseName = "c");

    void closeFile();

    void execute(const double & t);

private:

    void parseCollectionFile();

    const float _dx;
    const float _dy;
    const VField_t & _EFD;

    std::string _baseName;

    unsigned int _iteration;

    std::vector<NDIndex<DIM> > _lDoms;
    std::vector<double> _timesteps;
    IpplTimings::TimerRef _outputTimer;
};

inline
void PlotBinaryVTKCommand::closeFile()
{
    std::ofstream vtkout;
    short EndianTest_s;
    unsigned char *EndianTest = reinterpret_cast<unsigned char*>(&EndianTest_s);
    EndianTest[0] = 1;
    EndianTest[1] = 0;

    std::stringstream fname;
    fname << _baseName << "_collection.pvd";
    vtkout.open(fname.str().c_str());
    vtkout << "<?xml version=\"1.0\"?>\n"
           << indent_l0 << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"" << (EndianTest_s == 1? "LittleEndian": "BigEndian") << "\">\n"
           << indent_l1 << "<Collection>\n";
    for (unsigned int i = 0; i < _iteration; ++ i) {
        vtkout << indent_l2 << "<DataSet timestep=\"" << _timesteps[i] << "\" group=\"\" part=\"0\" file=\"" << _baseName << "_" << std::setw(4) << std::setfill('0') << i << ".pvtr\"/>\n";
    }
    vtkout << indent_l1 << "</Collection>\n"
           << indent_l0 << "</VTKFile>\n"
           << std::endl;
    vtkout.close();

}

#endif
