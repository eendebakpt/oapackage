""" Orthogonal Array package

"""

import numpy as np
import oapackage
import oalib
import oapackage.Doptim


def autodoctest():
    """ Test the module using autodoc
    Example:
      >>> import oapackage
      >>> arrayclass=oapackage.arraydata_t(2, 40, 0, 7)
      >>> print(arrayclass)
      arrayclass: N 40, k 7, strength 0, s {2,2,2,2,2,2,2}, order 0

    """
    return


def test_numpy_interface(verbose=0):
    A = np.eye(3, 4).astype(int)
    A[0, :] = [10, 20, 30, 50]
    A[2, :] = [-3, -4, -5, -6]

    if verbose:
        print('makearraylink')
    al = oapackage.makearraylink(A)

    np.testing.assert_array_equal(np.array(al), A)

    if verbose:
        al.showarray()

    if verbose:
        print('direct')
    al = oapackage.array_link(A)
    if verbose:
        al.showarray()
    Ax = np.array(al)
    if verbose:
        print(A)
        print(Ax)

    if 0:
        # not possible right now...
        if verbose:
            print('direct float')
        A = np.eye(3).astype(float)
        al = oapackage.array_link(A)
        if verbose:
            al.showarray()


def test_nauty(verbose=0):
    if verbose:
        print('test_nauty: test reduction to normal form')
    al = oapackage.exampleArray(1, verbose)
    alr = al.randomperm()
    tr = oapackage.reduceOAnauty(alr)
    alx = tr.apply(alr)
    assert(alx == al)


def unittest(verbose=1):
    """ Perform some unit testing, return True if succesfull """
    if verbose:
        print('oapackage: unittest: oalib version %s' % oalib.version())
    al = oalib.array_link()
    ii = 0
    al = oalib.exampleArray(ii, 0)

    arrayclass = oalib.arraydata_t(2, 16, 2, 6)  # fine
    alist = [arrayclass.create_root()]
    r = oalib.extend_arraylist(alist, arrayclass)
    if len(r) != 3:
        raise Exception('extension generation for %s failed' % (arrayclass, ))

    if not isinstance(al.getarray(), np.ndarray):
        print(
            'oapackage: unittest: error: array interface not working properly')
    else:
        if not al[2, 0] == al.getarray()[2, 0]:
            print(
                'oapackage: unittest: error: array interface not working properly')

    arrayclass = oalib.arraylink2arraydata(al)

    if verbose >= 2:
        print('unittest: calculate efficiencies')
    Deff = al.Defficiency()
    aa = oalib.Aefficiencies(al)
    if verbose >= 2:
        print('## oapackage test: example array %d: Deff %.3f' % (ii, Deff))

    # DOP reduction
    if verbose >= 2:
        print('unittest: test delete-one-factor GWLP reduction')
    al = oalib.exampleArray(5, 1)
    al2 = al.randomperm()

    alr = al.reduceDOP()
    al2r = al2.reduceDOP()
    if not alr == al2r:
        print('error: DOP reduced arrays unequal!: %d' % (alr == al2r))
        print('alr')
        alr.showarraycompact()
        print('al2r')
        al2r.showarraycompact()
        return False

    at = oalib.reductionDOP(al)
    check = at.apply(al) == al.reduceDOP()
    if not check:
        print('error: DOP reduction transformation is invalid')

    # test graphtools
    if verbose >= 2:
        print('unittest: test graphtools')
    from oapackage.graphtools import oa2graph
    arrayclass = oalib.arraylink2arraydata(al)
    _ = oa2graph(al, arrayclass)

    test_numpy_interface()
    test_nauty()

    return True

#%%
if __name__ == "__main__":
    """ Test code """
    unittest()
