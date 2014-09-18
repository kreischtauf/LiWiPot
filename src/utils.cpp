/***************************************************************************
                              utils.cpp
                         -------------------
    begin                : Thu Jan 12 2012
    copyright            : (C) 2012 by Christof Kraus
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

#include "utils.hh"

extern std::ofstream dbg;
#define DBGOUT dbg << __FILE__ << ": " << __LINE__ << "\t"


Utils::KahanAccumulation::KahanAccumulation():
    total(0.0),
    correction(0.0)
{ }

Utils::KahanAccumulation & Utils::KahanAccumulation::sum(const double & value)
{
    long double y = value - correction;
    long double t = total + y;
    correction = (t - total) - y;
    total = t;
    return *this;
}
