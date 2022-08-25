# parts that all tests have in common
add_library(leastsquares leastsquares.cc)
target_include_directories(leastsquares PRIVATE ${CMAKE_SOURCE_DIR})

# solve test problems from Moré, Garbow and Hillstrom with L-BFGS
foreach(problem Bard Beale BiggsEXP6 Box3D BrownBadlyScaled BrownDennis
                ExtendedRosenbrock FreudensteinRoth Gaussian
                JennrichSampson KowalikOsborne Meyer Osborne1 Osborne2
                PowellBadlyScaled PowellSingular Rosenbrock Wood)

  string(TOLOWER ${problem} problem_lower)
  set(name ${problem})

  configure_file(test_functions.cc ${name}.cc)
  add_executable(${name} ${name}.cc)
  target_link_libraries(${name} PRIVATE leastsquares)
  target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/test)

  if (dune-common_FOUND)
    dune_add_test(TARGET ${name})
  else()
    add_test(${name} ${name})
  endif()
endforeach()

# solve Rosenbrock function with several different configurations
# requires Dune::ParameterTree, and therefore dune-common
if (dune-common_FOUND)
  foreach(conjugation fr prp+ hs dy hz none)
    foreach(quasi_newton true false)

      set(problem Rosenbrock)
      string(TOLOWER ${problem} problem_lower)

      set(verbosity 0)
      set(max_iter 9999)
      set(reset_iter 1000)
      set(use_scale true)
      set(storage_limit 10)
      set(linesearch hager-zhang)
      set(interpolation quadratic)

      set(extrapolation false)

      set(name ${problem}_${quasi_newton}_${extrapolation}_${conjugation})
      configure_file(test_configs.cc ${name}.cc)
      add_executable(${name} ${name}.cc)
      target_link_libraries(${name} PRIVATE leastsquares)
      target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/test)

      dune_add_test(TARGET ${name})

      # restrict N-GMRES tests to installations with Eigen3
      if (EIGEN3_FOUND)
        set(extrapolation true)

        set(name ${problem}_${quasi_newton}_${extrapolation}_${conjugation})
        configure_file(test_configs.cc ${name}.cc)
        add_executable(${name} ${name}.cc)
        target_link_libraries(${name} PRIVATE leastsquares)
        target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/test)

        dune_add_test(TARGET ${name})
      endif()

    endforeach()
  endforeach()
endif()

