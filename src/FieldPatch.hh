/***************************************************************************
                            FieldPatch.hh
                         -------------------
    begin                : Thu Aug 23 2012
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

#ifndef FIELDPATCH_HH
#define FIELDPATCH_HH

#include "defs.hh"

#include <boost/mpi/operations.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>

#include <vector>

template<class T>
class BoundaryLayer;

template<class T>
class FieldPatch {
public:
    typedef T T_t;

    FieldPatch();
    FieldPatch(const FieldPatch & fp);
    FieldPatch(const NDIndex<DIM> & ldom);

    enum FieldType {EX, EY, HZ};

    void assign(const FieldPatch<T_t> & fp);
    void reset(const T_t & value = 0.0);
    void clear();

    T_t & getRelative(const size_t & i, const size_t & j);
    const T_t & getRelative(const size_t & i, const size_t & j) const;
    T_t & getAbsolute(const size_t & i, const size_t & j);
    const T_t & getAbsolute(const size_t & i, const size_t & j) const;
    T_t & operator()(const size_t & i, const size_t & j);
    const T_t & operator()(const size_t & i, const size_t & j) const;

    size_t size() const;
    void resize(NDIndex<DIM> dom);

    void add(const FieldPatch<T_t> & fp);
    void addWithoutResize(const FieldPatch<T_t> & fp);
    struct AddOp {
        FieldPatch operator()(FieldPatch<T_t> & fp1, const FieldPatch<T_t> & fp2){
            fp1.add(fp2);
            return fp1;
        }
    };

    FieldPatch<T_t> & operator=(const FieldPatch<T_t> & fp);
    FieldPatch<T_t> & operator=(const T_t & val);

    typename std::vector<T_t>::iterator begin();
    typename std::vector<T_t>::iterator end();

    void setOrigin(const Vector_t & orig);
    const Vector_t & getOrigin() const;

    void setSpacing(const Vector_t & dx);
    const Vector_t & getSpacing() const;

    const NDIndex<DIM> & getDomain() const;

    void setType(const FieldType & type);
    const FieldType & getType() const;

    void serialize(std::vector<char> & container) const;
    void deserialize(const std::vector<char> & container, size_t & startIndex);

private:
    friend class boost::serialization::access;
    friend class BoundaryLayer<T_t>;

    template<class Archive>
    void save(Archive& ar, const unsigned int version) const;

    template<class Archive>
    void load(Archive& ar, const unsigned int version);

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    NDIndex<DIM> _lDom;

    Vector_t _origin;
    Vector_t _spacing;

    FieldType _type;

    std::vector<T_t> _data;
};

template<class T>
FieldPatch<T>::FieldPatch():
    _data(1,0.0)
{
    _lDom[0] = Index(0, -1, -2);
    _lDom[1] = Index(0, -1, -2);
}

template<class T>
FieldPatch<T>::FieldPatch(const FieldPatch<T> & fp):
    _lDom(fp._lDom[0], fp._lDom[1]),
    _origin(fp._origin),
    _spacing(fp._spacing),
    _type(fp._type),
    _data(fp._data)
{ }

template<class T>
FieldPatch<T>::FieldPatch(const NDIndex<DIM> & ldom):
    _lDom(ldom)
{
    size_t s = std::max(size(), size_t(1));

    _data.resize(s, 0.0);
}

template<class T>
void FieldPatch<T>::assign(const FieldPatch<T> & fp)
{
    _lDom = fp._lDom;
    _origin = fp._origin;
    _spacing = fp._spacing;
    _data.assign(fp._data.begin(), fp._data.end());
}

template<class T>
void FieldPatch<T>::reset(const T & value)
{
    std::fill(_data.begin(), _data.end(), value);
}

template<class T>
void FieldPatch<T>::clear()
{
    _lDom[0] = Index(0,-1,-1);
    _lDom[1] = Index(0,-1,-1);
    _data.clear();
}

template<class T>
size_t FieldPatch<T>::size() const
{
    size_t s = _lDom[0].stride() > 0? _lDom[0].length(): 0;
    for (size_t i = 1; i < DIM; ++ i)
        s = _lDom[i].stride() > 0? s * _lDom[i].length(): 0;
    return s;
}

template<class T>
typename std::vector<T>::iterator FieldPatch<T>::begin()
{
    return _data.begin();
}

template<class T>
typename std::vector<T>::iterator FieldPatch<T>::end()
{
    return _data.end();
}

template<class T>
const NDIndex<DIM> & FieldPatch<T>::getDomain() const
{
    return _lDom;
}

template<class T>
void FieldPatch<T>::setType(const FieldType & type)
{
    _type = type;
}

template<class T>
const typename FieldPatch<T>::FieldType & FieldPatch<T>::getType() const
{
    return _type;
}

template<class T>
T & FieldPatch<T>::getRelative(const size_t & i, const size_t & j)
{
    size_t ii = j * _lDom[0].length() + i;
    // PAssert(j >= 0);
    // PAssert(i >= 0);
    // PAssert(ii < _data.size());
    return _data[ii];
}

template<class T>
const T & FieldPatch<T>::getRelative(const size_t & i, const size_t & j) const
{
    size_t ii = j * _lDom[0].length() + i;
    // PAssert(j >= 0);
    // PAssert(i >= 0);
    // PAssert(ii < _data.size());
    return _data[ii];
}

template<class T>
T & FieldPatch<T>::getAbsolute(const size_t & i, const size_t & j)
{
    return getRelative(i - _lDom[0].first(), j - _lDom[1].first());
}

template<class T>
const T & FieldPatch<T>::getAbsolute(const size_t & i, const size_t & j) const
{
    return getRelative(i - _lDom[0].first(), j - _lDom[1].first());
}

template<class T>
T & FieldPatch<T>::operator()(const size_t & i, const size_t & j)
{
    return getAbsolute(i,j);
}

template<class T>
const T & FieldPatch<T>::operator()(const size_t & i, const size_t & j) const
{
    return getAbsolute(i,j);
}

template<class T>
void FieldPatch<T>::setOrigin(const Vector_t & orig)
{
    _origin = orig;
}

template<class T>
const Vector_t & FieldPatch<T>::getOrigin() const
{
    return _origin;
}

template<class T>
void FieldPatch<T>::setSpacing(const Vector_t & dx)
{
    _spacing = dx;
}

template<class T>
const Vector_t & FieldPatch<T>::getSpacing() const
{
    return _spacing;
}

template<class T>
void FieldPatch<T>::serialize(std::vector<char> & container) const
{
    container.reserve(container.size() +
                      2 * DIM * sizeof(int) +
                      2 * DIM * sizeof(double) +
                      1 +
                      _data.size() * sizeof(T));
    const char* buffer;
    int ld[2 * DIM];
    double md[2 * DIM];
    for (unsigned int d = 0; d < DIM; ++ d) {
        ld[2 * d] = _lDom[d].first();
        ld[2 * d + 1] = _lDom[d].last();
        md[2 * d] = _origin(d);
        md[2 * d + 1] = _spacing(d);
    }
    buffer = reinterpret_cast<const char*>(ld);
    container.insert(container.end(), buffer, buffer + 2 * DIM * sizeof(int));

    buffer = reinterpret_cast<const char*>(md);
    container.insert(container.end(), buffer, buffer + 2 * DIM * sizeof(double));

    char fpt = _type;
    buffer = reinterpret_cast<const char*>(&fpt);
    container.insert(container.end(), buffer, buffer + 1);

    buffer = reinterpret_cast<const char*>(&(_data[0]));
    container.insert(container.end(), buffer, buffer + _data.size() * sizeof(T));
}

template<class T>
void FieldPatch<T>::deserialize(const std::vector<char> & container, size_t & startIndex)
{
    size_t it = startIndex;
    for (unsigned int d = 0; d < DIM; ++ d) {
        int lower = *reinterpret_cast<const int*>(&container[it]);
        it += sizeof(int);
        int upper = *reinterpret_cast<const int*>(&container[it]);
        it += sizeof(int);
        _lDom[d] = Index(lower, upper);
    }
    for (unsigned int d = 0; d < DIM; ++ d) {
        _origin(d) = *reinterpret_cast<const double*>(&container[it]);
        it += sizeof(double);
        _spacing(d) = *reinterpret_cast<const double*>(&container[it]);
        it += sizeof(double);
    }
    char fpt = container[it++];
    switch(fpt) {
    case 0:
        _type = EX;
        break;
    case 1:
        _type = EY;
        break;
    case 2:
        _type = HZ;
        break;
    }

    size_t sizeData = _lDom.size();
    _data.resize(sizeData);
    for (unsigned int l = 0; l < sizeData; ++ l) {
        _data[l] = *reinterpret_cast<const double*>(&container[it]);
        it += sizeof(double);
    }
    startIndex = it;
}

template<class T>
template<class Archive>
void FieldPatch<T>::save(Archive& ar, const unsigned int version) const
{
    const int & fi = _lDom[0].first();
    const int & li = _lDom[0].last();
    const int & fj = _lDom[1].first();
    const int & lj = _lDom[1].last();
    ar & fi;
    ar & li;
    ar & fj;
    ar & lj;
    for (unsigned int d = 0; d < DIM; ++ d) {
        double v = _origin(d);
        ar & v;
        v = _spacing(d);
        ar & v;
    }

    ar & _type;

    ar & _data;
}

template<class T>
template<class Archive>
void FieldPatch<T>::load(Archive& ar, const unsigned int version)
{
    int fi, li, fj, lj, s;
    ar & fi;
    ar & li;

    s = li >= fi? 1 : li - fi;
    _lDom[0] = Index(fi, li, s);

    ar & fj;
    ar & lj;

    s = lj >= fj? 1 : lj - fj;
    _lDom[1] = Index(fj, lj, s);

    for (unsigned int d = 0; d < DIM; ++ d) {
        double v;
        ar & v;
        _origin(d) = v;
        ar & v;
        _spacing(d) = v;
    }

    ar & _type;

    ar & _data;
}

// namespace boost {
//     namespace mpi {

//         //        template<>
//         template<class T>
//         struct is_commutative<FieldPatch<T>::AddOp, FieldPatch<T> >
//             : mpl::true_ { };

//     }
// }

template<class T>
void FieldPatch<T>::resize(NDIndex<DIM> dom)
{
    size_t s = 1;
    {
        for (size_t i = 0; i < DIM; ++ i) {
            if (dom[i].stride() <= 0) return;
            s *= dom[i].length();
        }
        if (s == 0) return;

        size_t i = 0;
        while (dom[i] == _lDom[i] && i < DIM) ++ i;
        if (i == DIM)  return;
    }

    std::vector<T> combinedData(s, 0.0);

    size_t ii = 0;
    size_t II = 0;

    int fi = std::max(_lDom[0].first(), dom[0].first());
    int li = std::min(_lDom[0].last(), dom[0].last());
    int fj = std::max(_lDom[1].first(), dom[1].first());
    int lj = std::min(_lDom[1].last(), dom[1].last());

    for (int j = fj; j <= lj; ++ j) {
        ii = (j - _lDom[1].first()) * _lDom[0].length() + fi - _lDom[0].first();
        II = (j - dom[1].first()) * dom[0].length() + fi - dom[0].first();

        for (int i = fi; i <= li; ++ i) {
            combinedData[II++] = _data[ii++];
        }
    }
    _data.assign(combinedData.begin(), combinedData.end());

    _lDom[0] = dom[0];
    _lDom[1] = dom[1];
}

template<class T>
void FieldPatch<T>::add(const FieldPatch<T> & fp)
{
    if (fp._lDom[0].stride() <= 0 ||
        fp._lDom[1].stride() <= 0) return;

    int fi = std::min(_lDom[0].first(), fp._lDom[0].first());
    int li = std::max(_lDom[0].last(), fp._lDom[0].last());
    int fj = std::min(_lDom[1].first(), fp._lDom[1].first());
    int lj = std::max(_lDom[1].last(), fp._lDom[1].last());

    size_t s = (li - fi + 1) * (lj - fj + 1);

    if (s > size()) {
        resize(NDIndex<DIM>(Index(fi,li), Index(fj,lj)));
    }

    for (int j = fp._lDom[1].first(); j <= fp._lDom[1].last(); ++ j) {
        size_t ii = (j - fp._lDom[1].first()) * fp._lDom[0].length();
        size_t II = (j -    _lDom[1].first()) *    _lDom[0].length() + fp._lDom[0].first() - _lDom[0].first();

        for (int i = fp._lDom[0].first(); i <= fp._lDom[0].last(); ++ i) {
            _data[II++] += fp._data[ii++];
        }
    }
    return;
}

template<class T>
void FieldPatch<T>::addWithoutResize(const FieldPatch<T> & fp)
{
    if (fp._lDom[0].stride() <= 0 ||
        fp._lDom[1].stride() <= 0) return;

    int fi = std::max(_lDom[0].first(), fp._lDom[0].first());
    int li = std::min(_lDom[0].last(), fp._lDom[0].last());
    int fj = std::max(_lDom[1].first(), fp._lDom[1].first());
    int lj = std::min(_lDom[1].last(), fp._lDom[1].last());

    for (int j = fj; j <= lj; ++ j) {
        size_t ii = (j - fp._lDom[1].first()) * fp._lDom[0].length() + fi - fp._lDom[0].first();
        size_t II = (j -    _lDom[1].first()) *    _lDom[0].length() + fi -    _lDom[0].first();

        for (int i = fi; i <= li; ++ i) {
            _data[II++] += fp._data[ii++];
        }
    }
    return;
}

template<class T>
FieldPatch<T> & FieldPatch<T>::operator=(const FieldPatch<T> & fp)
{
    _lDom = fp._lDom;
    _origin = fp._origin;
    _spacing = fp._spacing;
    _data.assign(fp._data.begin(), fp._data.end());

    return *this;
}

template<class T>
FieldPatch<T> & FieldPatch<T>::operator=(const T & val)
{
    for (int j = _lDom[1].first(); j <= _lDom[1].last(); ++ j) {
        for (int i = _lDom[0].first(); i <= _lDom[0].last(); ++ i) {
            (*this)(i,j) = val;
        }
    }

    return *this;
}
#endif
