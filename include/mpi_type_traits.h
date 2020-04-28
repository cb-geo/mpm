#ifndef MPM_MPI_TYPE_TRAITS_H_
#define MPM_MPI_TYPE_TRAITS_H_

// MPI
#ifdef USE_MPI
#include "mpi.h"
namespace mpm {

// MPI Type traits
template <typename Ttype>
struct MPI_Type_Traits {
  static inline MPI_Datatype mpi_type(Ttype&& raw);
};

// Specialization of the MPI_Type_Traits for primitive data types types
#define PRIMITIVE_DATATYPES(Type, MPI_Type)                     \
  template <>                                                   \
  inline MPI_Datatype MPI_Type_Traits<Type>::mpi_type(Type&&) { \
    return MPI_Type;                                            \
  }

PRIMITIVE_DATATYPES(bool, MPI::BOOL);

PRIMITIVE_DATATYPES(char, MPI::CHAR);
PRIMITIVE_DATATYPES(int, MPI::INT);
PRIMITIVE_DATATYPES(long, MPI::LONG);
PRIMITIVE_DATATYPES(unsigned int, MPI::UNSIGNED);
PRIMITIVE_DATATYPES(unsigned long, MPI::UNSIGNED_LONG);
PRIMITIVE_DATATYPES(unsigned long long, MPI::UNSIGNED_LONG_LONG);

PRIMITIVE_DATATYPES(float, MPI::FLOAT);
PRIMITIVE_DATATYPES(double, MPI::DOUBLE);
PRIMITIVE_DATATYPES(long double, MPI::LONG_DOUBLE);

#undef PRIMITIVE_DATATYPES

// MPI Const type traits
template <typename Ttype>
struct MPI_Type_Traits<const Ttype> {
  static inline MPI_Datatype mpi_type(const Ttype& value) {
    return MPI_Type_Traits<Ttype>::mpi_type(Ttype());
  }
};

// MPI Vector type traits
template <typename Ttype>
struct MPI_Type_Traits<std::vector<Ttype>> {
  static inline MPI_Datatype mpi_type(std::vector<Ttype>&& vector_data) {
    return MPI_Type_Traits<Ttype>::mpi_type(Ttype());
  }
};

}  // namespace mpm

#endif  // MPI

#endif  // MPM_MPI_TYPE_TRAITS_H_
