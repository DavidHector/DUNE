import hashlib
from dune.generator.generator import SimpleGenerator

from ._vtk import *
from dune.grid import OutputType
from dune.grid.grid_generator import _writeVTKDispatcher as addDispatcher

generator = SimpleGenerator("VtkWriter", "Dune::Vtk")

def load(includes, typeName):
    includes += ["dune/python/vtk/writer.hh"]
    moduleName = "vtkwriter_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName)
    return module.VtkWriter

######################################################################

allCollectors = {
  "continuous":
    lambda grid:
        "ContinuousDataCollector<"+grid._typeName+">",
  "lagrange":
    lambda grid, **kwargs: 
        "Dune::Vtk::LagrangeDataCollector<"+grid._typeName+","+str(kwargs["order"])+">"
  }
allWriters = {
  "default":  [ lambda grid:
                       "Dune::VtkUnstructuredGridWriter<"+grid._typeName+">",
                ["dune/vtk/writers/vtkunstructuredgridwriter.hh"] ],
  "lagrange": [ lambda grid, **kwargs:
                       "Dune::VtkUnstructuredGridWriter<"+grid._typeName+","+\
                               allCollectors["lagrange"](grid,**kwargs)+">",
                ["dune/vtk/writers/vtkunstructuredgridwriter.hh",
                 "dune/vtk/datacollectors/lagrangedatacollector.hh"] ]
  }

def vtkWriter(grid, name,
              version="default",
              pointData=None, pointScalar=None,
              cellData=None,  cellScalar=None,
              number=None, write=True,
              format = FormatTypes.binary,
              datatype = DataTypes.Float32, headertype = DataTypes.UInt32,
              **kwargs      # additional ctor arguments for dataContainer or Writer 
             ):
    writer = load( allWriters[version][1],
                   allWriters[version][0](grid,**kwargs) )\
                   (grid, format, datatype, headertype, **kwargs)

    # method to extract the range dimension from a given grid function
    def rangeDimension(f):
        try: return len(f)
        except: pass
        try: return f.ufl_shape[0]
        except: pass
        try: return f.dimRange
        except: pass
        return None
    # reinterpret `functions` argument list as list of the form
    # return list with entries [f,{name,components}] or [f,{name}]
    def extract(functions, toScalar):
        ret = []
        if type(functions)==dict: # form name:function, or (name,components):function
            for k,v in functions.items():
                if type(k) == str:
                    ret += [ [v,{"name":k}] ]
                else: # name component pair assumed
                    ret += [ [v,{"name":k[0], "components":k[1]}] ]
        elif type(functions)==list: # each entry is a function with a 'name'
                                    # and possibly a 'components' attribute
            for v in functions:
                try:
                    ret += [ [v,{"name":v.name, "components":v.components}] ]
                except AttributeError:
                    ret += [ [v,{"name":v.name}] ]
        if toScalar: # convert scalar functions with multiply components
            scalarRet = []
            for r in ret:
                # can we extract range dimension from a grid function?
                dimR = rangeDimension(r[0])
                comp = r[1].get("components",range(dimR))
                if dimR is None or len(comp)==1:
                     scalarRet += [ [ r[0], {"name":r[1]["name"],
                                             "components":comp} ] ]
                else:
                    for i, c in enumerate(comp):
                        scalarRet += [ [ r[0], {"name":r[1]["name"]+"_"+str(i),
                                                "components":[c]} ] ]
            return scalarRet
        else:
            return ret
    # call the actual `add` method on the writer here other packages can
    # add dispatch method to convert `f` into a `Vtk::Function`, i.e., a
    # ufl expression.
    def addDefault(addMethod, f, **kwargs):
        try:
            addMethod(f,**kwargs)
            return
        except TypeError:
            pass
        for dispatch in addDispatcher:
            func = dispatch(grid,f)
            if func is not None:
                addMethod(func,**kwargs)
                return
        raise AttributeError("coundn't find a way to add function "+str(f)+
                             " of type "+str(type(f))+" to the VtkWriter")

    # add all the functions provided as arguments to the writer
    for f in extract(pointData, toScalar=False):
        addDefault( writer.addPointData, f[0], **f[1], range=RangeTypes.auto)
    for f in extract(pointScalar, toScalar=True):
        addDefault( writer.addPointData, f[0], **f[1], range=RangeTypes.scalar)
    for f in extract(cellData,False):
        addDefault( writer.addCellData,  f[0], **f[1], range=RangeTypes.vector)
    for f in extract(cellScalar,True):
        addDefault( writer.addCellData,  f[0], **f[1], range=RangeTypes.scalar)

    # the class is returned. It doesn't allow for adding additional data
    # since after the first `write` this doesn't make sense. The name is
    # fixed by the provided argument and in addition a consecutive number
    # is added for each execution of __call__.
    class _Writer:
        def __init__(self,writer,name,number):
            self.writer = writer
            self.name = name
            self.number = number
        def __call__(self):
            if self.number is None:
                self.writer.write(self.name)
            else:
                self.writer.write(self.name,self.number)
                self.number += 1
    _writer = _Writer(writer,name,number)
    if write:
        _writer()
    return _writer

#######################################################

def vtkUnstructuredGridWriter(grid):
    """create a VtkWriter for unstructured grids

    Args:
        grid:  grid view to use for the vtk file

    Returns:
        vtkWriter: the constructed writer
    """
    return load( allWriters["default"][1],
                 allWriters["default"][0](grid) )(grid)
