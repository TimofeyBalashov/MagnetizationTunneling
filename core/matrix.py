_modules = {'numpy': NumpyModule()
           ,'mpmath': MPmathModule()
           }
_default = _modules['numpy']

class MathModule(object):
    def matrix(a):
        raise NotImplementedError()

    def cmatrix(a):
        raise NotImplementedError()

    def diag(a, d=0):
        raise NotImplementedError()

    def zeros(*shape):
        raise NotImplementedError()

    def ones(*shape):
        raise NotImplementedError()
    
    def multiply(m1, m2):
        raise NotImplementedError()
    
    def eigh(m):
        raise NotImplementedError()
    
class NumpyModule(MathModule):
    def matrix(a):
        return np.matrix(a, dtype=double)
    
    def cmatrix(a):
        return np.matrix(a, dtype=complex)
    
    def diag(a, d=0):
        return np.matrix(np.diag(a, d))
    
    def zeros(*shape):
        return np.matrix(np.zeros(shape))
    
    def ones(*shape):
        return np.matrix(np.ones(shape))
    
    def multiply(m1, m2):
        return np.multiply(m1,m2)
    
    def eigh(m):
        raise np.linalg.eigh(m)   
    
class MPmathModule(MathModule):
    def set_precision(dps):
        mp.dps = dps
    
    def matrix(a):
        return mp.matrix(a)
    
    def cmatrix(a):
        return mp.matrix(a)*(1+0j)
    
    def diag(a, d=0):
        return mp.matrix(np.diag(a, d))
    
    def zeros(*shape):
        return mp.zeros(*shape)

    def ones(*shape):
        return mp.ones(*shape)

    def multiply(m1, m2):
        return mp.matrix(np.multiply(m1.tolist(),m2.tolist()))
    
    def eigh(m):
        raise mp.eigh(m)       
    
def _module(module)
    if module not in _modules:
        raise ValueError(
            "Module '{}' not supported. Supported modules are {}".format(
                modules, ', '.join(_modules)))
    else:
        return _modules[module]
    
def mathmodule(module=None):
    if module is None:
        return _default
    elif isinstance(module, MatrixModule):
        return module
    else:
        return _module(module)
    
def set_default(module):
    _default = _module(module)