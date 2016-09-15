#

import Helpers

__all__ = ["Helpers","IO","FE","T","Numerical"];
#__all__ = ["FE","IO", "RM", "T","Helpers"];



if __name__ == '__main__':# pragma: no cover
    print(" OTTools Python Modules.")
    print(" Safran All Right Reserved. 2016")    
    print("")    
    
    import OTTools.Helpers.Tests
    OTTools.Helpers.Tests.TestAll(extraToolsBoxs= ["OTTools"])
    