/***************************************************************************
                           Communicator.hh
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

#ifndef COMMUNICATOR_HH
#define COMMUNICATOR_HH

#include <numeric>

#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>

#include "defs.hh"
#include "utils.hh"
#include "FieldPatch.hh"

#ifdef OpenSSL_FOUND
#include <openssl/md5.h>
#else
#define MD5_DIGEST_LENGTH 0
#endif

#define P2PCOMM

class Communicator {
public:
    static std::vector<NDIndex<DIM> > getAllLocalDomains(NDIndex<DIM> lDom);

    template <class T>
    static void communicateFields(FieldPatch<T> & fp,
                                  const NDIndex<DIM> & dom,
                                  const bool  all2all = false);

    template <class T>
    static void communicateFields(FieldPatch<T> & fp,
                                  const std::vector<NDIndex<DIM> >  & localPDomains,
                                  const std::vector<NDIndex<DIM> >  & localFDomains,
                                  const bool all2all = false);

    template <class T>
    static void collectivelyCommunicateFields(FieldPatch<T> & fp,
                                              const NDIndex<DIM> & dom);

    template <class T>
    static void p2pCommunicateFields(FieldPatch<T> & fp,
                                     const std::vector<NDIndex<DIM> >  & localPDomains,
                                     const std::vector<NDIndex<DIM> >  & localFDomains,
                                     const bool all2all = false);

    static void initialize();
    static const std::vector<unsigned int> & getNeighbours();
    static const std::vector<unsigned int> & getNode2Index();
    static const std::vector<unsigned int> & getIndex2Node();
    static std::pair<unsigned int, unsigned int> getNodeMeshSize();

private:
    static std::vector<unsigned int> _index2node;
    static std::vector<unsigned int> _node2index;
    static std::vector<unsigned int> _neighbours;
    static unsigned int _nNodeCols;
    static unsigned int _nNodeRows;
    static bool _initialized;
};


inline
const std::vector<unsigned int> & Communicator::getNeighbours()
{
    initialize();
    return _neighbours;
}

inline
const std::vector<unsigned int> & Communicator::getNode2Index()
{
    initialize();
    return _node2index;
}

inline
const std::vector<unsigned int> & Communicator::getIndex2Node()
{
    initialize();
    return _index2node;
}

inline
std::pair<unsigned int, unsigned int> Communicator::getNodeMeshSize()
{
    initialize();
    return std::make_pair(_nNodeCols, _nNodeRows);
}

template <class T>
void Communicator::communicateFields(FieldPatch<T> & fp,
                                     const NDIndex<DIM> & dom,
                                     const bool all2all)
{
#ifdef P2PCOMM
    NDIndex<DIM> lPDom = fp.getDomain();

    auto localPDomains = getAllLocalDomains(lPDom);
    auto localFDomains = getAllLocalDomains(dom);

    p2pCommunicateFields(fp, localPDomains, localFDomains, all2all);
#else
    collectivelyCommunicateFields(fp, dom);
#endif
}

template <class T>
void Communicator::communicateFields(FieldPatch<T> & fp,
                                     const std::vector<NDIndex<DIM> >  & localPDomains,
                                     const std::vector<NDIndex<DIM> >  & localFDomains,
                                     const bool all2all)
{
    p2pCommunicateFields(fp, localPDomains, localFDomains, all2all);
}

template <class T>
void Communicator::collectivelyCommunicateFields(FieldPatch<T> & fp,
                                                 const NDIndex<DIM> & dom)
{
    boost::mpi::communicator world;

    std::vector<FieldPatch<T> > patches;
    boost::mpi::all_gather(world, fp, patches);

    fp.resize(dom);

    for (size_t i = 0; i < patches.size(); ++ i) {
        if (i == (size_t) Ippl::myNode() || patches[i].size() == 0) {
            patches[i].clear();
            continue;
        }
        fp.addWithoutResize(patches[i]);
        patches[i].clear();
    }
}

template <class T>
void Communicator::p2pCommunicateFields(FieldPatch<T> & fp,
                                        const std::vector<NDIndex<DIM> >  & localPDomains,
                                        const std::vector<NDIndex<DIM> >  & localFDomains,
                                        const bool all2all)
{
    boost::mpi::communicator world;
    initialize();

    const unsigned int myNode = Ippl::myNode();
    const unsigned int numNodes = Ippl::getNodes();

    std::vector<unsigned int> partners;
    if (all2all) {
        partners.resize(numNodes - 1);
        for (unsigned int k = 0; k < numNodes - 1; ++ k)
            partners[k] = k < myNode? k: k + 1;
    } else {
        partners.assign(_neighbours.begin(), _neighbours.end());
    }

    std::vector<FieldPatch<T> > mySplitPatches(partners.size());
    std::vector<FieldPatch<T> > othersSplitPatches(partners.size());

    int totalmsgsend = 0;
    std::vector<boost::mpi::request> requests;
    int tag = Ippl::Comm->next_tag(F_GUARD_CELLS_TAG, F_TAG_CYCLE);

    size_t totalmsgrecv = 0;
    for (unsigned int k: partners) {
        bool inside = true;
        NDIndex<DIM> dom;
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = std::max(localPDomains[k][d].first(), localFDomains[myNode][d].first());
            int upper = std::min(localPDomains[k][d].last(),  localFDomains[myNode][d].last());
            if (lower > upper) {
                inside = false;
                break;
            }
            dom[d] = Index(lower, upper);
        }
        if (!inside) continue;

        boost::mpi::request req = world.irecv(k, tag, &othersSplitPatches[totalmsgrecv], 1);
        ++ totalmsgrecv;
        requests.push_back(req);
    }

    for (unsigned int k: partners){
        NDIndex<DIM> dom;
        bool inside = true;
        for (unsigned int d = 0; d < DIM; ++ d) {
            int lower = std::max(localPDomains[myNode][d].first(), localFDomains[k][d].first());
            int upper = std::min(localPDomains[myNode][d].last(),  localFDomains[k][d].last());
            if (lower > upper) {
                inside = false;
                break;
            }
            dom[d] = Index(lower, upper);
        }
        if (!inside) continue;

        mySplitPatches[totalmsgsend].resize(dom);
        mySplitPatches[totalmsgsend].addWithoutResize(fp);

        // world.send(k, tag, &mySplitPatches[totalmsgsend], 1);
        boost::mpi::request req = world.isend(k, tag, &mySplitPatches[totalmsgsend], 1);
        requests.push_back(req);
        ++ totalmsgsend;
    }

    fp.resize(localFDomains[myNode]);

    boost::mpi::wait_all(&requests[0], &requests[0] + requests.size());

    for (size_t i = 0; i < totalmsgrecv; ++ i) {

        fp.addWithoutResize(othersSplitPatches[i]);
    }
}

class WrappedVector: public Vector_t {
public:
    WrappedVector():
        Vector_t()
    { }

    WrappedVector(const Vector_t &rhs):
        Vector_t(rhs)
    { }

    WrappedVector(const double& x00):
        Vector_t(x00)
    { }

    WrappedVector(const double& x00, const double& x01):
        Vector_t(x00, x01)
    { }

    WrappedVector(const double& x00, const double& x01, const double& x02):
        Vector_t(x00, x01, x02)
    { }

private:
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int version) const
    {
        for (unsigned int d = 0; d < DIM; ++ d) {
            double a = (*this)[d];
            ar & a;
        }
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int version)
    {
        for (unsigned int d = 0; d < DIM; ++ d) {
            double a;
            ar & a;
            (*this)[d] = a;
        }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

};

template <class T>
class Wrapper
{
public:
    typedef T W_t;
};

template<>
struct Wrapper<Vector_t>
{
    typedef WrappedVector W_t;
};

#endif
