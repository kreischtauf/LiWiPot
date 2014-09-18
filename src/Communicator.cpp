/***************************************************************************
                           Communicator.cpp
                         -------------------
    begin                : Sat Sep 21 2013
    copyright            : (C) 2013 by Christof Kraus
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

#include "Communicator.hh"
#include "utils.hh"

extern std::ofstream dbg;

#define DBGOUT dbg << "Communicator.cpp: " << __LINE__ << "\t"

std::vector<unsigned int> Communicator::_index2node;
std::vector<unsigned int> Communicator::_node2index;
std::vector<unsigned int> Communicator::_neighbours;
unsigned int Communicator::_nNodeCols = 0;
unsigned int Communicator::_nNodeRows = 0;
bool Communicator::_initialized = false;

void Communicator::initialize()
{
    if (_initialized) return;

    _index2node.resize(Ippl::getNodes());
    _node2index.resize(Ippl::getNodes());

    unsigned int depth = static_cast<unsigned int>(floor(log(1.0 * Ippl::getNodes()) / log(2.0) + 0.5));
    unsigned int numNodes = Utils::ipow(2, depth);
    if (numNodes != static_cast<unsigned int>(Ippl::getNodes()) && Ippl::myNode() == 0) {
        Ippl::exitAllNodes("\033[31;1mnumber of cores not a power of 2!!\033[0m", true);
    }
    unsigned int coldepth = static_cast<unsigned int>(floor(0.5 * (depth + 1.0)));
    unsigned int rowdepth = static_cast<unsigned int>(floor(0.5 * depth));
    _nNodeCols = Utils::ipow(2, coldepth);
    _nNodeRows = Utils::ipow(2, rowdepth);
    std::vector<std::pair<unsigned int, unsigned int> > n2mi(numNodes, std::make_pair(0,0));
    unsigned int nNodes = 1;
    for (unsigned int m = 1; m <= depth; ++ m) {
        for (unsigned int l = nNodes; l > 0; -- l) {
            if (m % 2 == 0) {
                coldepth = static_cast<unsigned int>(floor(0.5 * (m + 1.0)));
                unsigned int curcol = Utils::ipow(2, coldepth);
                n2mi[2*(l-1)].second = n2mi[l-1].second + (n2mi[l-1].second / curcol) * curcol;
                n2mi[2*(l-1)+1].second = n2mi[2*(l-1)].second + curcol;
            } else {
                n2mi[2*(l-1)].second = 2 * n2mi[l-1].second;
                n2mi[2*(l-1)+1].second = n2mi[2*(l-1)].second + 1;
            }
        }
        nNodes *= 2;
    }
    for (unsigned int i = 0; i < numNodes; ++ i) {
        n2mi[i].first = i;
        _node2index[i] = n2mi[i].second;
    }
    std::sort(n2mi.begin(), n2mi.end(),
              [](const std::pair<unsigned int, unsigned int> & a,
                 const std::pair<unsigned int, unsigned int> & b) { return a.second < b.second;});

    for (unsigned int i = 0; i < numNodes; ++ i) {
        _index2node[i] = n2mi[i].first;
    }

    unsigned int index = _node2index[Ippl::myNode()];
    unsigned int myCol = index % _nNodeCols;
    unsigned int myRow = index / _nNodeCols;
    unsigned int startCol = std::min(myCol, myCol - 1);
    unsigned int endCol   = std::min(_nNodeCols, myCol + 2);
    unsigned int startRow = std::min(myRow, myRow - 1);
    unsigned int endRow   = std::min(_nNodeRows, myRow + 2);

    for (unsigned int j = startRow; j < endRow; ++ j) {
        for (unsigned int i = startCol; i < endCol; ++ i) {
            if (i == myCol && j == myRow) continue;
            index = j * _nNodeCols + i;
            _neighbours.push_back(_index2node[index]);
        }
    }
    std::sort(_neighbours.begin(), _neighbours.end());

    _initialized = true;
}

std::vector<NDIndex<DIM> > Communicator::getAllLocalDomains(NDIndex<DIM> lDom)
{
    unsigned int myNode = Ippl::myNode();
    unsigned int numNodes = Ippl::getNodes();
    int *values = new int[2 * DIM * numNodes];
    for (unsigned int d = 0; d < DIM; ++ d) {
        values[2 * (myNode * DIM + d)    ] = lDom[d].first();
        values[2 * (myNode * DIM + d) + 1] = lDom[d].last();
    }

    MPI_Allgather(MPI_IN_PLACE, 2 * DIM, MPI_INT, values, 2 * DIM, MPI_INT, MPI_COMM_WORLD);

    std::vector<NDIndex<DIM> > ret(numNodes);
    for (unsigned int k = 0; k < numNodes; ++ k) {
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = values[2 * (k * DIM + d)    ];
            int upper = values[2 * (k * DIM + d) + 1];
            int sign = lower <= upper? 1: upper - lower;
            ret[k][d] = Index(lower, upper, sign);
        }
    }

    delete[] values;

    return ret;
}
