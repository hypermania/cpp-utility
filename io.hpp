#ifndef IO_HPP
#define IO_HPP
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

#include <Eigen/Dense>

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
