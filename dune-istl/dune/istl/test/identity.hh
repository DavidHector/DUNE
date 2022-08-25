// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef IDENTITY_HH
#define IDENTITY_HH
#include <dune/istl/bcrsmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/scalarmatrixview.hh>

template<class B, class Alloc>
void setupIdentitySparsityPattern(Dune::BCRSMatrix<B, Alloc>& A, int N)
{
    typedef typename Dune::BCRSMatrix<B, Alloc> Matrix;
    A.setSize(N, N, N);
    A.setBuildMode(Matrix::row_wise);
    for (typename Dune::BCRSMatrix<B,Alloc>::CreateIterator i = A.createbegin(); i != A.createend(); ++i)
    {
        i.insert(i.index());
    }
}

template<class B, class Alloc>
void setupIdentity(Dune::BCRSMatrix<B,Alloc>& A, int N)
{
    typedef typename Dune::BCRSMatrix<B,Alloc>::field_type FieldType;
    setupIdentitySparsityPattern(A, N);
    for (typename Dune::BCRSMatrix<B,Alloc>::RowIterator i = A.begin(); i != A.end(); ++i) {
        i->operator[](i.index()) = 1.0;
    }
}
#endif