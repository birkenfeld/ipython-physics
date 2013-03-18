from physics import Q

def test_format():
    a = Q('1.45 m')
    s = "{}".format(a)
    assert s == '1.45 m'
    s = "{:.4e}".format(a)
    assert s == '1.4500E+0 m'
    s = "{.value}".format(a)
    assert s == '1.45'
    s = "%s" % a
    assert s == '1.45 m'

