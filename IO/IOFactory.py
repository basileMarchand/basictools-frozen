# -*- coding: utf-8 -*-


from BasicTools.Helpers.Factory import Factory

def RegisterReaderClass(name, classtype, constructor=None, withError = True):
    return ReaderFactory.RegisterClass(name,classtype, constructor=constructor, withError = withError )

def CreateReader(name,ops=None):
    return ReaderFactory.Create(name,ops)

class ReaderFactory(Factory):
    _Catalog = {}
    def __init__(self):
        super(ReaderFactory,self).__init__()

def GetAvailableReaders():
    return list(ReaderFactory._Catalog.keys())

def InitAllReaders():
    try:
        import BasicTools.IO.InpReader as InpReader
    except:
        pass

    import BasicTools.IO.AscReader as AscReader
    import BasicTools.IO.GeofReader as GeofReader
    import BasicTools.IO.GeoReader as GeoReader
    import BasicTools.IO.GmshReader as GmshReader
    import BasicTools.IO.MeshReader as MeshReader
    import BasicTools.IO.GReader as GReader
    import BasicTools.IO.FemReader as FemReader
    from BasicTools.IO.StlReader import ReadStl
    from BasicTools.IO.XdmfReader import ReadXdmf
    from BasicTools.IO.PipeIO import PipeReader
    from BasicTools.IO.OdbReader import OdbReader
    from BasicTools.IO.UtReader import UtReader
    from BasicTools.IO.VtuReader import VtuReader



def InitAllWriters():

    from BasicTools.IO.GeofWriter import GeofWriter
    from BasicTools.IO.GmshWriter import GmshWriter
    from BasicTools.IO.MeshWriter import MeshWriter
    from BasicTools.IO.OdbWriter  import OdbWriter
    from BasicTools.IO.StlWriter  import StlWriter
    from BasicTools.IO.XdmfWriter import XdmfWriter
    from BasicTools.IO.PipeIO import PipeWriter
    from BasicTools.IO.CsvWriter import CsvWriter

def RegisterWriterClass(name, classtype, constructor=None, withError = True):
    WriterFactory.RegisterClass(name,classtype, constructor=constructor, withError = withError )

def CreateWriter(name,ops=None):
    return WriterFactory.Create("."+name.split(".")[-1],ops)

class WriterFactory(Factory):
    _Catalog = {}
    def __init__(self):
        super(WriterFactory,self).__init__()

def GetAvailableWriter():
    return list(WriterFactory._Catalog.keys())




def CheckIntegrity():
    from BasicTools.IO.IOFactory import WriterFactory
    from BasicTools.IO.IOFactory import GetAvailableReaders
    from BasicTools.IO.IOFactory import GetAvailableWriter
    ##
    InitAllReaders()
    print("Available Readers : ", GetAvailableReaders())

    InitAllWriters()
    print("Available Writers : ", GetAvailableWriter())

    return "ok"


if __name__ == '__main__':
    print(CheckIntegrity())# pragma: no cover