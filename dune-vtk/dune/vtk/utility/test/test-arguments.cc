#include <config.h>

#include <dune/vtk/utility/arguments.hh>
#include <dune/vtk/utility/errors.hh>

using namespace Dune::Vtk;

template <class... Args>
void test_found(Args const&... args)
{
  VTK_ASSERT(getArg<double>(args...,1.0) == 2.0);
  VTK_ASSERT(getArg<float>(args...,2.0f) == 3.0f);
  VTK_ASSERT(getArg<int>(args...,3) == 4);
  VTK_ASSERT(getArg<unsigned int>(args...,4u) == 5u);
}

template <class... Args>
void test_notfound(Args const&... args)
{
  VTK_ASSERT(getArg<double>(args...,1.0) == 1.0);
  VTK_ASSERT(getArg<float>(args...,2.0f) == 2.0f);
  VTK_ASSERT(getArg<int>(args...,3) == 3);
  VTK_ASSERT(getArg<unsigned int>(args...,4u) == 4u);
}

struct A{};

int main()
{
  test_found(2.0, 3.0f, 4, 5u);
  test_found(5u, 4, 3.0f, 2.0);

  test_notfound();
  test_notfound(4l);
  test_notfound(A{},4ul);
}
