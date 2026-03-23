/// @file utils.hpp
/*-------------------------get a list of sparse matrices by their names----------------------*/
#include <unordered_map>

using namespace Eigen;
using namespace tmbutils;

template<class Type>
struct SparseMatrixList : vector<SparseMatrix<Type>> {
  
private:
  std::vector<std::string> names;
  std::unordered_map<std::string, size_t> nameIndex;
  
public:
  SparseMatrixList(SEXP listZs) {
    this->resize(Rf_length(listZs));
    // Extract names from R
    SEXP nm = PROTECT(Rf_getAttrib(listZs, R_NamesSymbol));
    if (nm != R_NilValue) {
      for (int i = 0; i < Rf_length(listZs); i++) {
        std::string matName = std::string(Rf_translateCharUTF8(STRING_ELT(nm, i)));
        names.push_back(matName);
        nameIndex[matName] = i;
      }
    }
    UNPROTECT(1);
    // Convert each R sparse matrix to Eigen sparse
    for (int idx = 0; idx < Rf_length(listZs); idx++) {
      SEXP spm = VECTOR_ELT(listZs, idx);
      if (!isValidSparseMatrix(spm))
        error("Matrix '%s' is not a sparse matrix", names[idx].c_str());
      (*this)[idx] = asSparseMatrix<Type>(spm);
    }
  }
  //getName=> search matrices by names;
  const SparseMatrix<Type>& getName(const std::string &name) const {
    auto it = nameIndex.find(name);
    if (it != nameIndex.end()) return (*this)[it->second];
    error(("Matrix name not found: " + name).c_str());
    return (*this)[0];
  }
  
};

// utils to check
namespace mmtmb {
  
  template<class Type>
  bool isEmpty(const vector<Type>& x) {
    return x.size() == 0;
  }

  template<class Type>  
  bool isEmpty(const array<Type>& x) {
    return x.size() == 0;
  }

  template<class Type>
  bool isEmpty(const matrix<Type>& X) {
    return (X.rows() == 0 || X.cols() == 0);
  }

  template<class Type>
  bool isEmpty(const Eigen::SparseMatrix<Type>& X) {
    return (X.rows() == 0 || X.cols() == 0);
  }
  
  // Catch-all for unknown types with error message
  template<class T>
  bool isEmpty(const T& x) {
    static_assert(sizeof(T) == 0, "mmtmb::isEmpty() called with unsupported type. Use vector, array, or matrix.");
    return false;
  }
  
  template<class Type>
  bool isSparseValid(const Eigen::SparseMatrix<Type>& X) {
    return (X.rows() > 1 || X.cols() > 1);
  }
  
  template<class Type>
  bool isNA(Type x){
    return R_IsNA(asDouble(x));
  }

  template<class Type>
  bool notFinite(Type x) {
    return (!R_FINITE(asDouble(x)));
  }

  template<class Type>
  void fillMatId(SparseMatrix<Type>& X) {
    X.setIdentity();
  }

  template<class Type>
  void fillMatId(matrix<Type>& X) {
    X.setIdentity();
  }

  // Error for any other type
  template<class T>
  void fillMatId(T& X) {
    static_assert(sizeof(T) == 0, "mmtmb::fillMatId() called with unknown-matrix type");
  }
  
}
