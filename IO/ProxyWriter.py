# -*- coding: utf-8 -*-


class ProxyWriter(object):
    def __init__(self):
        super(ProxyWriter,self).__init__()
        self.writers = []

    def __getattribute__(self,name):
        writers = object.__getattribute__(self, "writers")
        def newfunc(*args, **kwargs):
            for w in writers:
                res = w.__getattribute__(name)(*args,**kwargs)
            return res
        return newfunc

def CheckIntegrity(GUI = False):
    return "ok"

if __name__ == '__main__':
    print(CheckIntegrity(GUI=True)) # pragma: no cover
