#ifndef UTILS_HH
#define UTILS_HH

#include <cassert>
#include "defs.hh"

namespace Utils {

    void getLocalDomains(const FieldLayout<DIM> & FLayout,
                         std::vector<NDIndex<DIM> > & vectorLocalDomains);

    void getNodeIndices(const FieldLayout<DIM> & cellCenteredLayout,
                        std::vector<int> & vectorNodeIndices);

    void increaseLocalDomainsLast(const FieldLayout<DIM> & cellCenteredLayout,
                                  std::vector<NDIndex<DIM> > & vectorLocalDomains);

    void decreaseLocalDomainsLast(const FieldLayout<DIM> & cellCenteredLayout,
                                  std::vector<NDIndex<DIM> > & vectorLocalDomains);

    void increaseLocalDomainsFirst(const FieldLayout<DIM> & cellCenteredLayout,
                                   std::vector<NDIndex<DIM> > & vectorLocalDomains);

    void decreaseLocalDomainsFirst(const FieldLayout<DIM> & cellCenteredLayout,
                                   std::vector<NDIndex<DIM> > & vectorLocalDomains);

    void addGostCellToLocalDomains(const FieldLayout<DIM> & cellCenteredLayout,
                                   std::vector<NDIndex<DIM> > & vectorLocalDomains);

    template<class M, class T>
    class CompanionDeleter {
    public:
        CompanionDeleter(T* t);
        void operator()(M* m);
    private:
        T* companion_;
    };

    template<class T>
    T ipow(T base, unsigned int exp);

    float Q_rsqrt( float x );

    class KahanAccumulation
    {
    private:
        long double total;
        long double correction;

    public:
        KahanAccumulation();
        KahanAccumulation& sum(const double & value);

        KahanAccumulation& operator+(const double & value)
        {
            return sum(value);
        }

        KahanAccumulation& operator-(const double & value)
        {
            return sum(-value);
        }
    };

    template <typename T> int sgn(T val) {
        return (T(0) <= val) - (val < T(0));
    }
}

template<class M, class T>
Utils::CompanionDeleter<M,T>::CompanionDeleter(T* t):
companion_(t)
{ }

template<class M, class T>
void
Utils::CompanionDeleter<M,T>::operator()(M* master)
{
    delete master;
    delete companion_;
}

template<class T>
T Utils::ipow(T base, unsigned int exp)
{
    T result = 1;
    while (exp)
        {
            if (exp & 1)
                result *= base;
            exp >>= 1;
            base *= base;
        }

    return result;
}

inline
void Utils::getLocalDomains(const FieldLayout<DIM> & FLayout,
                            std::vector<NDIndex<DIM> > & vectorLocalDomains)
{
    vectorLocalDomains.resize(Ippl::getNodes());
    vectorLocalDomains[Ippl::myNode()] = FLayout.getLocalNDIndex();
    for (FieldLayout<DIM>::const_iterator_dv domainMapIterator = FLayout.begin_rdv();
         domainMapIterator != FLayout.end_rdv();
         ++ domainMapIterator) {
        RefCountedP<Vnode<DIM> > virtualNode = (*domainMapIterator).second;
        vectorLocalDomains[virtualNode->getNode()] = virtualNode->getDomain();
    }
}

inline
void Utils::getNodeIndices(const FieldLayout<DIM> & cellCenteredLayout,
                           std::vector<int> & vectorNodeIndices)
{
    vectorNodeIndices.resize(Ippl::getNodes());
    vectorNodeIndices[Ippl::myNode()] = Ippl::myNode();
    for (FieldLayout<DIM>::const_iterator_dv domainMapIterator = cellCenteredLayout.begin_rdv();
         domainMapIterator != cellCenteredLayout.end_rdv();
         ++ domainMapIterator)
        {
            RefCountedP<Vnode<DIM> > virtualNode = (*domainMapIterator).second;
            vectorNodeIndices[virtualNode->getNode()] = virtualNode->getNode();
        }
}

inline
void Utils::increaseLocalDomainsLast(const FieldLayout<DIM> & cellCenteredLayout,
                                     std::vector<NDIndex<DIM> > & vectorLocalDomains)
{
    const NDIndex<DIM> & globalDomain = cellCenteredLayout.getDomain();

    for (std::vector<NDIndex<DIM> >::iterator NDIndexIterator = vectorLocalDomains.begin();
         NDIndexIterator != vectorLocalDomains.end();
         ++ NDIndexIterator) {
        NDIndex<DIM> & domain = *NDIndexIterator;
        if (domain[0].last() == globalDomain[0].last()) {
            domain[0] = Index(domain[0].first(), globalDomain[0].last() + 1);
        }

        if (domain[1].last() == globalDomain[1].last()) {
            domain[1] = Index(domain[1].first(), globalDomain[1].last() + 1);
        }
    }
}

inline
void Utils::decreaseLocalDomainsLast(const FieldLayout<DIM> & cellCenteredLayout,
                                     std::vector<NDIndex<DIM> > & vectorLocalDomains)
{
    const NDIndex<DIM> & globalDomain = cellCenteredLayout.getDomain();

    for (std::vector<NDIndex<DIM> >::iterator NDIndexIterator = vectorLocalDomains.begin();
         NDIndexIterator != vectorLocalDomains.end();
         ++ NDIndexIterator) {
        NDIndex<DIM> & domain = *NDIndexIterator;
        if (domain[0].last() == globalDomain[0].last()) {
            domain[0] = Index(domain[0].first(), globalDomain[0].last() - 1);
        }

        if (domain[1].last() == globalDomain[1].last()) {
            domain[1] = Index(domain[1].first(), globalDomain[1].last() - 1);
        }
    }
}

inline
void Utils::increaseLocalDomainsFirst(const FieldLayout<DIM> & cellCenteredLayout,
                                      std::vector<NDIndex<DIM> > & vectorLocalDomains)
{
    const NDIndex<DIM> & globalDomain = cellCenteredLayout.getDomain();

    for (std::vector<NDIndex<DIM> >::iterator NDIndexIterator = vectorLocalDomains.begin();
         NDIndexIterator != vectorLocalDomains.end();
         ++ NDIndexIterator) {
        NDIndex<DIM> & domain = *NDIndexIterator;
        if (domain[0].first() == globalDomain[0].first()) {
            domain[0] = Index(globalDomain[0].first() - 1, domain[0].last());
        }

        if (domain[1].first() == globalDomain[1].first()) {
            domain[1] = Index(globalDomain[1].first() - 1, domain[1].last());
        }
    }
}

inline
void Utils::decreaseLocalDomainsFirst(const FieldLayout<DIM> & cellCenteredLayout,
                                      std::vector<NDIndex<DIM> > & vectorLocalDomains)
{
    const NDIndex<DIM> & globalDomain = cellCenteredLayout.getDomain();

    for (std::vector<NDIndex<DIM> >::iterator NDIndexIterator = vectorLocalDomains.begin();
         NDIndexIterator != vectorLocalDomains.end();
         ++ NDIndexIterator) {
        NDIndex<DIM> & domain = *NDIndexIterator;
        if (domain[0].first() == globalDomain[0].first()) {
            domain[0] = Index(globalDomain[0].first() + 1, domain[0].last());
        }

        if (domain[1].first() == globalDomain[1].first()) {
            domain[1] = Index(globalDomain[1].first() + 1, domain[1].last());
        }
    }
}


inline
void Utils::addGostCellToLocalDomains(const FieldLayout<DIM> & cellCenteredLayout,
                                      std::vector<NDIndex<DIM> > & vectorLocalDomains)
{
    const NDIndex<DIM> & globalDomain = cellCenteredLayout.getDomain();

    for (std::vector<NDIndex<DIM> >::iterator NDIndexIterator = vectorLocalDomains.begin();
         NDIndexIterator != vectorLocalDomains.end();
         ++ NDIndexIterator) {
        NDIndex<DIM> & domain = *NDIndexIterator;

        if (domain[0].first() != globalDomain[0].first()) {
            domain[0] = Index(domain[0].first() - 1, domain[0].last());
        }

        if (domain[0].last() != globalDomain[0].last()) {
            domain[0] = Index(domain[0].first(), domain[0].last() + 1);
        }

        if (domain[1].first() != globalDomain[1].first()) {
            domain[1] = Index(domain[1].first() - 1, domain[1].last());
        }

        if (domain[1].last() != globalDomain[1].last()) {
            domain[1] = Index(domain[1].first(), domain[1].last() + 1);
        }
    }
}

inline
float Utils::Q_rsqrt( float x )
{
    float xhalf = 0.5f*x;
    float xqhalf = 0.5010936f*x;
    uint32_t i;
    assert(sizeof(x) == sizeof(i));
    std::memcpy(&i, &x, sizeof(i));
    i = 0x5f3759df - (i>>1);
    std::memcpy(&x, &i, sizeof(i));
    x = x*(1.5f - xhalf*x*x);        // result is less than correct result for all x since 1/sqrt(x) is convex
    x = x*(1.5010936f - xqhalf*x*x); // optimized slope (steeper slope pushes result to higher values)
    return x;
}
#endif
