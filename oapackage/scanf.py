import oapackage._scanf


def sscanf(inputString, formatString):
    """ Simple scanf function

    Args:
        inputString (str): string to be parsed
        formatString (str): specification
    Returns:
        list: list of parsed arguments
    """
    return oapackage._scanf.scanf(formatString, inputString)
