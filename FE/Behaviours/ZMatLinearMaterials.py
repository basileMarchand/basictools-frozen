# -*- coding: utf-8 -*-

almanac = {}

almanac["Elastic"] = """
***behavior linear_elastic
 **elasticity isotropic
   young {YOUNG}
   poisson {POISSON}
***return"""

def GetMaterial(name,props):
    print(props)
    return almanac[name].format(**props)


def CheckIntegrity(GUI=False):
    props = {"YOUNG":4.e5,"POISSON":0.3}
    print(GetMaterial("Elastic",props))
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(True))# pragma: no cover
