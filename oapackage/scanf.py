import scanf as _scanfpip  # https://pypi.org/project/scanf/


def sscanf(inputString, formatString):
    """ Simple scanf function """
    return _scanfpip.scanf(formatString, inputString)
