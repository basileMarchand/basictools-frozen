# -*- coding: utf-8 -*-



class FieldNotFound(ValueError):
     """Exception to treat Field Not found """
     
     def __init__(self, value):
         self.value = 'Field "' + value + '" not found, Sorry!!'
         
     def __str__(self):
         return repr(self.value) # pragma: no cover
         

def CheckIntegrity():
    
    FieldNotFound('toto');
    return 'OK'  
    
if __name__ == '__main__':
    print(CheckIntegrity()) # pragma: no cover 
      
