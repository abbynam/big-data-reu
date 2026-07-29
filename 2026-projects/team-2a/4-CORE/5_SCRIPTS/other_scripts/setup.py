from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy

# ext_module = Extension("cc2csv",
#               ["cc2csv_cy.pyx"],
#               libraries=["m"],
#               extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp" ],
#               extra_link_args=['-fopenmp']
#               )
#
#
# setup(
#   ext_modules = cythonize("*.pyx")
# )

ext_modules=[
    # Extension("cc2csv_cy",
    #           ["cc2csv_cy.pyx"],
    #           libraries=["m"],
    #           extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp"],
    #           extra_link_args=['-fopenmp']
    #           ),

    Extension("utilities",
              ["utilities.pyx"],
              libraries=["m"],
              include_dirs=[numpy.get_include()],
              extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp"],
              extra_link_args=['-fopenmp']
              ),
    #
    # Extension("energy_matcher",
    #           ["energy_matcher.pyx"],
    #           libraries=["m"],
    #           extra_compile_args=["-O3", "-ffast-math", "-march=native", "-fopenmp"],
    #           extra_link_args=['-fopenmp']
    #           )
    Extension("energy_matcher",
              ["energy_matcher.pyx"],
              libraries=["m"],
              extra_compile_args=["-O3", "-ffast-math", "-march=native", "-fopenmp"],
              extra_link_args=['-fopenmp']
              )

]

setup(
  name = "thread_demo",
  cmdclass = {"build_ext": build_ext},
  ext_modules = ext_modules
)