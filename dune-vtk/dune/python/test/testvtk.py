import math
from dune.grid import structuredGrid, gridFunction
from dune.vtk import vtkUnstructuredGridWriter, vtkWriter, RangeTypes, FieldInfo
grid = structuredGrid([0,0],[1,1],[5,5])

def test1(writer):
    @gridFunction(grid)
    def f(x):
        return x[0]*x[1]
    writer.addPointData(f, name="f")
    del f

@gridFunction(grid)
def g(x):
    return [math.sin(2*math.pi*x[0]*x[1]), x[0]*x[1]]

writer1 = vtkUnstructuredGridWriter(grid)
test1(writer1)
writer1.write("test1")

writer2 = vtkUnstructuredGridWriter(grid)
writer2.addPointData(g, name="g")
writer2.write("test2")

@gridFunction(grid)
def g(x):
    return [math.sin(2*math.pi*x[0]*x[1]), x[0]*x[1]]*5

writer3 = vtkUnstructuredGridWriter(grid)
test1(writer3)
writer3.addPointData(g, name="g")
writer3.write("test3")

writer4 = vtkUnstructuredGridWriter(grid)
writer4.addPointData(g, name="g10",  components=(0,1))    # neither
writer4.addPointData(g, name="g012", components=[0,1,2])  # is vector
writer4.addPointData(g, name="g23",  components=[2,3])    # neiter
writer4.addPointData(g, name="g2",   components=[2])      # is scalar 
writer4.addPointData(g, name="g23V", range=RangeTypes.vector,
                                     components=[2,3])    # is vector
writer4.addPointData(g, name="g23S", range=RangeTypes.scalar,
                                     components=[2,3])    # is scalar
writer4.addPointData(g, info=FieldInfo(name="g2S",  range=RangeTypes.scalar),
                                     components=[2])      # is scalar
writer4.addPointData(g, info=FieldInfo(name="g2V",  range=RangeTypes.vector),
                                     components=[2])      # is vector
writer4.write("test4")


@gridFunction(grid)
def f(x):
    return math.sin(2*math.pi*x[0]*x[1]/(1+x[0]*x[0]))

"""
def vtkWriter(grid, name,
              version="unstructured",
              cellScalar=None, pointScalar=None,
              cellVector=None, pointVector=None,
              number=None, outputType=OutputType.appendedbase64,
              write=True,
              **kwargs      # additional ctor arguments for writer
"""
# Note: {"gV":g} has only one components in the data file although the range is auto
writer = vtkWriter( grid, "test",
                    pointData   = { "gV":g, ("g43V",(4,3)):g, ("g435V",(4,3,5)):g },
                    pointScalar = { "fS":f, "gS":g, ("g23S",(2,3)):g },
                  )
writer = vtkWriter( grid, "testLagrange",
                    version="lagrange", order=3,
                    pointData   = {"gV":g, ("g43V",(4,3)):g},
                    pointScalar = {("g435V",(4,3,5)):g},
                    cellData    = {("g23V",(2,3)):g, ("g234V",(2,3,4)):g},
                    cellScalar  = { "fS":f, "gS":g}
                  )
