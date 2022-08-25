#ifndef Udune_ftworth_matrixutilities_HH
#define Udune_ftworth_matrixutilities_HH

#include <Eigen/Dense>

// check symmetry of a sparse matrix
template<typename MAT>
bool is_symmetric (const MAT& A, std::string s="", bool verbose=false, typename MAT::field_type tol=1e-13)
{
  unsigned count_unsym = 0;
  unsigned count_nonzeroes = 0;
  using FP = typename MAT::field_type;
  FP maxdiff = 0.0;
  
  for (auto rIt=A.begin(); rIt!=A.end(); ++rIt)
    {
      auto i = rIt.index();
      auto cIt = rIt->begin();
      auto cEndIt = rIt->end();
      for (; cIt!=cEndIt; ++cIt)
	{
	  count_nonzeroes += 1;
	  if (cIt.index()>=i)
	    {
	      auto j = cIt.index();
	      for (int compi=0; compi<(*cIt).N(); compi++)
		for (int compj=0; compj<(*cIt).M(); compj++)
		  {
		    auto mag = 0.5*(std::abs((*cIt)[compi][compj])+std::abs(A[j][i][compj][compi]));
		    if (mag<1e-12) continue; // do not consider zero entries
		    auto diff = std::abs((*cIt)[compi][compj]-A[j][i][compj][compi]);
		    auto reldiff = diff/mag;
		    maxdiff = std::max(maxdiff,reldiff);
		    if (reldiff > tol )
		      {
			count_unsym += 1;
			if (verbose)
			  std::cout << "!" << s << " : " << i << "," << j << "," << compi << "," << compj << " = " << (*cIt)[compi][compj] << " : " << A[j][i][compj][compi] << std::endl;
		      }
		  }
	    }
	}
    }
  if (verbose && count_unsym>0)
    std::cout << count_unsym << " entries not symmetric out of " << count_nonzeroes << " max difference =" << maxdiff << std::endl;
  return (count_unsym==0);
}

// check if two sparse matrices have the same entries
template<typename MAT>
bool is_equal (const MAT& A1, const MAT& A2, std::string s="", bool verbose=false, typename MAT::field_type tol=1e-13)
{
  if (A1.N() != A2.N())
    {
      if (verbose) std::cout << s << " row size is different" << std::endl;
      return false;
    }
  if (A1.M() != A2.M())
    {
      if (verbose) std::cout << s << " column size is different" << std::endl;
      return false;
    }

  unsigned count_nonzeroes = 0;
  using FP = typename MAT::field_type;
  FP maxdiff = 0.0;

  auto rIt1 = A1.begin();
  auto rIt2 = A2.begin();
  while (rIt1!=A1.end())
    {
      if (rIt2==A2.end())
	{
	  if (verbose) std::cout << s << " ups that should not happen" << std::endl;
	  return false;
	}
      auto i1 = rIt1.index();
      auto i2 = rIt2.index();
      if (i1!=i2)
	{
	  if (verbose) std::cout << s << " ups that should not happen" << std::endl;
	  return false;
	}

      // now compare rows
      auto cIt1 = rIt1->begin();
      auto cEndIt1 = rIt1->end();
      auto cIt2 = rIt2->begin();
      auto cEndIt2 = rIt2->end();
      while (cIt1!=cEndIt1)
	{
	  if (cIt2==cEndIt2)
	    {
	      if (verbose) std::cout << s << " rows have different number of entries" << std::endl;
	      return false;
	    }
	  auto j1 = cIt1.index();
	  auto j2 = cIt2.index();
	  if (j1!=j2)
	    {
	      if (verbose) std::cout << s << " entry with different column index found" << std::endl;
	      return false;
	    }

	  if ( (*cIt1).N()!=(*cIt2).N() || (*cIt1).M()!=(*cIt2).M() )
	    {
	      if (verbose) std::cout << s << " row or column size of block is different" << std::endl;
	      return false;
	    }
	    
	  // now compare this entry
	  for (int compi=0; compi<(*cIt1).N(); compi++)
	    for (int compj=0; compj<(*cIt1).M(); compj++)
	      {
		auto absdiff = std::abs( (*cIt1)[compi][compj] - (*cIt2)[compi][compj] );
		maxdiff = std::max(maxdiff,absdiff);
		if ( absdiff > tol && verbose)
		  std::cout << "!" << s << " : " << i1 << "," << j1 << "," << compi << "," << compj << " = " << (*cIt1)[compi][compj] << " : " << (*cIt2)[compi][compj] << std::endl;
	      }
	  
	  ++cIt1;
	  ++cIt2;
	}

      // next row
      ++rIt1;
      ++rIt2;
    }
  
  if (maxdiff>tol && verbose)
    std::cout << s << " matrices different. maxdiff=" << maxdiff << std::endl;
  
  return (maxdiff<=tol);
}

// check symmetry of a sparse matrix
template<typename MAT>
void make_symmetric (MAT& A)
{
  for (auto rIt=A.begin(); rIt!=A.end(); ++rIt)
    {
      auto i = rIt.index();
      auto cIt = rIt->begin();
      auto cEndIt = rIt->end();
      for (; cIt!=cEndIt; ++cIt)
	if (cIt.index()>i)
	  {
	    auto j = cIt.index();
	    for (int compi=0; compi<(*cIt).N(); compi++)
	      for (int compj=0; compj<(*cIt).M(); compj++)
		{
		  auto aij = (*cIt)[compi][compj];
		  auto aji = A[j][i][compj][compi];
		  (*cIt)[compi][compj] = A[j][i][compj][compi] = 0.5*(aij+aji);
		}
	  }
    }
}

// row equilibration
template<typename field_type> 
std::vector<field_type> row_equilibrate (Dune::DynamicMatrix<field_type>& A)
{
  std::vector<field_type> s(A.rows(),1.0);
  
  // row equilibration
  for (size_t i=0; i<A.rows(); i++)
    {
      field_type sum = 0.0;
      for (size_t j=0; j<A.cols(); j++)
	sum += std::abs(A[i][j]);
      if (sum>=1e-16)
	{
	  for (size_t j=0; j<A.cols(); j++)
	    A[i][j] /= sum;
	  s[i] = sum;
	}
    }
  return s;
}


// row equilibration
template<typename field_type> 
std::vector<field_type> row_equilibrate (Eigen::Matrix<field_type,Eigen::Dynamic,Eigen::Dynamic>& A)
{
  std::vector<field_type> s(A.rows(),1.0);
  
  // row equilibration
  for (size_t i=0; i<A.rows(); i++)
    {
      field_type sum = 0.0;
      for (size_t j=0; j<A.cols(); j++)
	sum += std::abs(A(i,j));
      if (sum>=1e-16)
	{
	  for (size_t j=0; j<A.cols(); j++)
	    A(i,j) /= sum;
	  s[i] = sum;
	}
    }
  return s;
}

// rank revealing LU decomposition using full pivoting
template<typename field_type, typename int_type> 
int_type rank_revealing_full_pivot_LU (Dune::DynamicMatrix<field_type>& A, std::vector<int_type>& rowpermut, std::vector<int_type>& colpermut, field_type threshold, field_type& min_pivot)
{
  if (colpermut.size()!=A.cols())
    DUNE_THROW(Dune::Exception, "rank_revealing_full_pivot_LU: col input size mismatch");
  if (rowpermut.size()!=A.rows())
    DUNE_THROW(Dune::Exception, "rank_revealing_full_pivot_LU: row input size mismatch");

  // initial norm for termination criterion
  auto initial_rowsumnorm = A.infinity_norm();
  
  // now determine rank of this matrix
  for (size_t j=0; j<colpermut.size(); j++) colpermut[j]=j;
  for (size_t j=0; j<rowpermut.size(); j++) rowpermut[j]=j;
  size_t matrix_rank=0;
  min_pivot = 1e100;
  for (size_t k=0; k<A.rows(); k++) // pivot row
    {
      // we treat also rectangular matrices
      if (k>=A.cols()) break;

      // compute rowsumnorm of rest matrix: this would be independent of pivoting
      field_type rowsumnorm = 0.0;
      for (size_t i=k; i<A.rows(); i++)
	{
	  field_type rowsum = 0.0;
	  for (size_t j=k; j<A.cols(); j++)
	    rowsum += std::abs(A[i][j]);
	  rowsumnorm = std::max(rowsumnorm,rowsum);
	}

      // termination criterion
      if (rowsumnorm<threshold*initial_rowsumnorm)
	break; // we cannot eliminate further

      // we can eliminate, increase rank by one
      matrix_rank++;

      // pivot search
      auto mag = std::abs(A[k][k]);
      size_t imax = k;
      size_t jmax = k;
      for (size_t i=k; i<A.rows(); i++)
      	for (size_t j=k; j<A.cols(); j++)
      	  if (std::abs(A[i][j])>mag)
      	    {
      	      mag = std::abs(A[i][j]);
      	      imax = i;
      	      jmax = j;
      	    }
      
      // interchange rows and columns if necessary
      if (imax!=k || jmax!=k)
	{
	  colpermut[k] = jmax;
	  rowpermut[k] = imax;
	  for (size_t i=0; i<A.rows(); i++)
	    std::swap(A[i][k],A[i][jmax]);
	  for (size_t j=0; j<A.cols(); j++)
	    std::swap(A[k][j],A[imax][j]);
	}

      min_pivot = std::min(std::abs(A[k][k]),min_pivot);
      
      // eliminate column
      for (size_t i=k+1; i<A.rows(); i++)
	{
	  field_type factor = A[i][k]/A[k][k];
	  A[i][k] = factor;
	  for (size_t j=k+1; j<A.cols(); j++)
	    A[i][j] -= factor*A[k][j];
	}
    }
  return matrix_rank;
}


// rank revealing Cholesky decomposition using full pivoting
// needs quadratic symmetric positive semidefinite matrix (really! it assumes nonnegativity of the diagonal!
// stops when all remaining diagonal elements are smaller than first_pivot*threshold
template<typename field_type, typename int_type> 
int_type rank_revealing_diagonal_pivot_Cholesky (Dune::DynamicMatrix<field_type>& A, std::vector<int_type>& permut, field_type threshold, field_type& min_pivot, int verbose=0)
{
  if (A.rows()!=A.rows())
    DUNE_THROW(Dune::Exception, "rank_revealing_diagonal_pivot_Cholesky: matrix not quadratic");
  if (permut.size()!=A.cols())
    DUNE_THROW(Dune::Exception, "rank_revealing_diagonal_pivot_Cholesky: permut array size does not match that of A");
  
  // initialize 
  for (size_t j=0; j<permut.size(); j++) permut[j]=j;
  size_t matrix_rank=0;
  field_type first_pivot;

  // the main loop
  for (size_t k=0; k<A.rows(); k++) // pivot row
    {
      // pivot search
      auto pivot = A[k][k];
      size_t imax = k;
      for (size_t i=k+1; i<A.rows(); i++)
	if (A[i][i]>pivot)
	  {
	    pivot = A[i][i];
	    imax = i;
	  }

      // termination criterion
      if (k==0)
	{
	  // first row/col
	  min_pivot = pivot;
	  if (pivot<=0.0)
	    {
	      if (verbose>0) std::cout << "cholesky: first pivot is not positive" << std::endl;
	      break; // we cannot eliminate further
	    }
	  first_pivot = pivot; // remember the first pivot, which is >0
	  if (verbose>0) std::cout << "first pivot is " << pivot << std::endl;
	}
      else
	{
	  // compare relative to the initial pivot
	  min_pivot = std::min(pivot,min_pivot);
	  if (pivot<first_pivot*threshold)
	    {
	      if (verbose>0) std::cout << " last good pivot was  " << A[k-1][k-1] << " remaining diagonal elements in step " << k << " are <=" << pivot << std::endl;
	      break;
	    }
	}

      // we can eliminate, increase rank by one
      matrix_rank++;
      
      // interchange rows and columns if necessary
      if (imax!=k)
	{
	  permut[k] = imax;
	  for (size_t i=0; i<A.rows(); i++)
	    std::swap(A[i][k],A[i][imax]);
	  for (size_t j=0; j<A.cols(); j++)
	    std::swap(A[k][j],A[imax][j]);
	}

      // eliminate column
      // this is not yet cholesky ... we could save half of the operations
      for (size_t i=k+1; i<A.rows(); i++)
	{
	  field_type factor = A[i][k]/A[k][k];
	  A[i][k] = factor;
	  for (size_t j=k+1; j<A.cols(); j++)
	    A[i][j] -= factor*A[k][j];
	}
    }
  return matrix_rank;
}

//! apply permutations to a right hand side vector
template<class P, class V>
void permute_forward (const P& p, V& b)
{
  if (b.size()!=p.size())
    DUNE_THROW(Dune::Exception,"permutation vector incompatible with b");
  
  for (std::size_t k=0; k<b.size()-1; ++k)
    if (p[k]!=k) std::swap(b[k],b[p[k]]);
}

//! apply permutations to a solution vector
template<class P, class V>
void permute_backward (const P& q, V& z)
{
  if (z.size()!=q.size())
    DUNE_THROW(Dune::Exception,"permutation vector incompatible with z");
  
  for (int k=z.size()-2; k>=0; --k)
    if (q[k]!=std::size_t(k)) std::swap(z[k],z[q[k]]);
}

// orthogonalize columns using modified Gram-Schmidt
// Q is an mxn dense matrix whose columns begin <= j < end (!) are to be orthonormalized
// When begin>0 then the columns 0 <= j < begin are assumed to be orthogonal and are not changed
// when end<n then all columns end <= j < n are not touched
// colpermut must be of size n; it is only modified in the range begin <= j < end
// the rank returned is the number of linearly independend columns in the range begin <= j < end
template<typename field_type, typename int_type> 
int_type rank_revealing_gram_schmidt (int begin, int end, Dune::DynamicMatrix<field_type>& Q,
				      std::vector<int_type>& colpermut, field_type threshold, field_type& min_norm)
{
  //std::cout << "QR " << Q.rows() << "x" << Q.cols() << std::endl;
  size_t rows=Q.rows();
  size_t cols=Q.cols();
    
  if (colpermut.size()!=cols)
    DUNE_THROW(Dune::Exception, "rank_revealing_gram_schmidt: input size mismatch");
  if (begin<0 || begin>end || end>cols)
    DUNE_THROW(Dune::Exception, "rank_revealing_gram_schmidt: begin or end out of bounds");

  // initialize rank
  int_type rank=0;

  // initialize permutation vector
  for (size_t j=begin; j<end; j++)
    colpermut[j] = j;
  
  // prepare norms
  std::vector<field_type> columnnorms2(cols,0.0);
  for (size_t i=0; i<rows; i++)
    for (size_t j=0; j<end; j++)
      columnnorms2[j] += Q[i][j]*Q[i][j];

  // startup phase !
  // subtract components of the first 0<=j<begin columns from columns begin<=j<end
  for (size_t k=0; k<begin; k++) // loop over all columns that are already orthogonal
    {
      // compute factors
      std::vector<field_type> alpha(cols,0.0);
      for (size_t i=0; i<rows; i++)
	{
	  auto Qik = Q[i][k];
	  for (size_t j=begin; j<end; j++)
	    alpha[j] += Qik*Q[i][j];
	}
      for (size_t j=begin; j<end; j++)
	{
	  alpha[j] /= columnnorms2[k];
	  columnnorms2[j] = 0.0; // will be updated below
	}

      // subtract component in direction k, forward looking; update norms
      for (size_t i=0; i<rows; i++)
	for (size_t j=begin; j<end; j++)
	  {
	    Q[i][j] -= alpha[j]*Q[i][k];
	    columnnorms2[j] += Q[i][j]*Q[i][j]; // recompute column norms as we go
	  }
    }

  // compute smallest norm of the columns in range(B)
  min_norm = 1e100;

  // modified Gram-Schmidt
  for (size_t k=begin; k<end; k++) // loop over all columns
    {	
      // find column with largest norm2
      size_t maxj=k;
      field_type maxnorm2=columnnorms2[k];
      for (size_t j=k+1; j<end; j++)
	if (columnnorms2[j]>maxnorm2)
	  {
	    maxj = j;
	    maxnorm2 = columnnorms2[j];
	  }
      // std::cout << "step k=" << k << " maxj=" << maxj << " norm2=" << maxnorm2 << std::endl;
      
      // interchange columns k and maxj if required
      if (maxj!=k)
	{
	  colpermut[k] = maxj;
	  std::swap(columnnorms2[k],columnnorms2[maxj]);
	  for (size_t i=0; i<rows; i++)
	    std::swap(Q[i][k],Q[i][maxj]);
	}

      // check if we have a zero column; then we stop
      //std::cout << "k=" << k << " mag=" << std::sqrt(columnnorms2[k]) << std::endl;
      if (std::sqrt(columnnorms2[k])<threshold)
	break;

      // we have column k with norm greater or equal threshold
      rank += 1;

      // compute factors
      std::vector<field_type> alpha(cols,0.0);
      for (size_t i=0; i<rows; i++)
	{
	  auto Qik = Q[i][k];
	  for (size_t j=k+1; j<end; j++)
	    alpha[j] += Qik*Q[i][j];
	}
      for (size_t j=k+1; j<end; j++)
	{
	  alpha[j] /= columnnorms2[k];
	  columnnorms2[j] = 0.0; // will be updated below
	}

      // subtract component in direction k, forward looking; update norms
      for (size_t i=0; i<rows; i++)
	for (size_t j=k+1; j<end; j++)
	  {
	    Q[i][j] -= alpha[j]*Q[i][k];
	    columnnorms2[j] += Q[i][j]*Q[i][j]; // recompute column norms as we go
	  }
	    	
      // normalize column k
      columnnorms2[k] = std::sqrt(columnnorms2[k]);
      min_norm = std::min(columnnorms2[k],min_norm);
      for (size_t i=0; i<rows; i++)
	Q[i][k] /= columnnorms2[k];
      columnnorms2[k] = 1.0;
    }

  return rank;
}


// orthogonalize columns using modified Gram-Schmidt
// this version uses column vectors stored in std::vector<V> where V is an istl block vector type
// mat is a matrix defining a scalar product. It is positive definite on the vectorspace spanned by precolumns and columns
// precolumns is a collection of linearly independent vectors already orthogonalized with respect to matrix
// colums are the columns to be orthogonalized w.r.t. precolums *and* itself using the scalar product given by matrix
// colpermut must be of same size as columns
// the rank returned is the number of linearly independend columns in columns
template<typename Mat, typename Vec, typename field_type, typename int_type> 
int_type rank_revealing_gram_schmidt (const Mat& Mpre, const Mat& M,
				      const std::vector<Vec> precolumns,
				      std::vector<Vec>& columns,
				      std::vector<int_type>& colpermut,
				      field_type threshold,
				      field_type& min_norm)
{
    
  if (colpermut.size()!=columns.size())
    DUNE_THROW(Dune::Exception, "rank_revealing_gram_schmidt: input size mismatch");

  // initialize rank
  int_type rank=0;

  size_t n; // length of the vectors
  if (columns.size()>0)
    n = columns[0].N();
  else
    return rank;  // input size is zero

  // initialize permutation vector
  for (size_t j=0; j<columns.size(); j++)
    colpermut[j] = j;

  // startup phase !
  // orthogonalize columns with respect to precolumns
  for (size_t k=0; k<precolumns.size(); k++) // loop over all columns that are already orthogonal
    {
      Vec Mv(n); // temporary vector
      Mpre.mv(precolumns[k],Mv); // use scalar product for pre phase
      auto norm2 = precolumns[k]*Mv;
      for (size_t j=0; j<columns.size(); j++)
	{
	  auto alpha = (columns[j]*Mv)/norm2;
	  columns[j].axpy(-alpha,precolumns[k]);
	}
    }

  // compute smallest norm of the columns in range(B)
  min_norm = 1e100;

  // orthogonalize columns using rank revealing modified Gram-Schmidt
  std::vector<field_type> norms2(columns.size(),0.0); // column norms for pivoting
  {
    // compute norms^2
    Vec Mv(n); // temporary vector
    for (size_t j=0; j<columns.size(); j++)
      {
	M.mv(columns[j],Mv);
	norms2[j] = columns[j]*Mv;
      }
  }
  for (size_t k=0; k<columns.size(); k++) // loop over all columns
    {
      // find column with largest norm2 (norms have been computed already)
      size_t maxj=k;
      field_type maxnorm2=norms2[k];
      for (size_t j=k+1; j<columns.size(); j++)
	if (norms2[j]>maxnorm2)
	  {
	    maxj = j;
	    maxnorm2 = norms2[j];
	  }
      // std::cout << "step k=" << k << " maxj=" << maxj << " norm2=" << maxnorm2 << std::endl;
      
      // interchange columns k and maxj if required
      if (maxj!=k)
	{
	  colpermut[k] = maxj;
	  std::swap(norms2[k],norms2[maxj]);
	  for (size_t i=0; i<n; i++)
	    std::swap(columns[k][i],columns[maxj][i]);
	}
      
      // check if we have a zero column; then we stop
      //std::cout << "k=" << k << " mag=" << std::sqrt(columnnorms2[k]) << std::endl;
      if (std::sqrt(norms2[k])<threshold)
	break;

      // we have column k with norm greater or equal threshold
      rank += 1;

      // orthogonalize columns j>k with column k
      Vec Mv(n); // temporary vector
      M.mv(columns[k],Mv);
      auto norm2 = columns[k]*Mv;
      for (size_t j=k+1; j<columns.size(); j++)
	{
	  auto alpha = (columns[j]*Mv)/norm2;
	  columns[j].axpy(-alpha,columns[k]);
	}
	    	
      // normalize column k
      columns[k] /= std::sqrt(norm2);

      // update norms for next round
      for (size_t j=k+1; j<columns.size(); j++)
	{
	  M.mv(columns[j],Mv);
	  norms2[j] = columns[j]*Mv;
	}
    }

  return rank;
}

#endif
