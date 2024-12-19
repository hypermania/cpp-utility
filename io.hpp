#ifndef IO_HPP
#define IO_HPP
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>



template<typename Scalar>
void write_to_file(std::vector<Scalar> &vector, std::string filename){
  char *memblock = (char *)&vector[0];
  std::ofstream file(filename, std::ios::binary);
  if(file.is_open()){
    file.write(memblock, vector.size() * sizeof(Scalar));
  }
}

#ifndef NOT_USING_EIGEN

#include <Eigen/Dense>
template<typename Derived>
void write_to_file(const Eigen::PlainObjectBase<Derived> &obj, std::string filename){
  std::ofstream file(filename, std::ios::binary);
  if(file.is_open()){
    file.write((char *)obj.data(), obj.size() * sizeof(typename Eigen::DenseBase<Derived>::Scalar));
  }
}

template<typename Derived>
void write_to_filename_template(const Eigen::PlainObjectBase<Derived> &obj, const std::string format_string, const int idx)
{
  char filename[128];
  sprintf(filename, format_string.data(), idx);
  std::ofstream file(filename, std::ios::binary);
  if(file.is_open()){
    file.write((char *)obj.data(), obj.size() * sizeof(typename Eigen::DenseBase<Derived>::Scalar));
  }
}

#endif

#endif
