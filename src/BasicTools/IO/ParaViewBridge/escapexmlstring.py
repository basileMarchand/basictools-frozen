
import argparse

parser = argparse.ArgumentParser(description='Process python source to escape xml strings.')

parser.add_argument('-f', metavar='file', type=str, nargs='+',  help='files to escape')

parser.add_argument('-r', type=str, nargs=1,  help='to replace')

parser.add_argument('-o', type=str, nargs=1, help='output file')

args = parser.parse_args()
print(args)


def escapeForXmlAttribute(s):

    # http://www.w3.org/TR/2000/WD-xml-c14n-20000119.html#charescaping
    # In character data and attribute values, the character information items "<" and "&" are represented by "&lt;" and "&amp;" respectively.
    # In attribute values, the double-quote character information item (") is represented by "&quot;".
    # In attribute values, the character information items TAB (#x9), newline (#xA), and carriage-return (#xD) are represented by "&#x9;", "&#xA;", and "&#xD;" respectively.

    s = s.replace('&', '&amp;') # Must be done first!
    s = s.replace('<', '&lt;')
    s = s.replace('>', '&gt;')
    s = s.replace('"', '&quot;')
    s = s.replace('\r', '&#xD;')
    s = s.replace('\n', '&#xA;')
    s = s.replace('\t', '&#x9;')
    return s

def SearchAndReplace(data):

    source = data[0]

    for i in range(1,len(data)):
        filtered = escapeForXmlAttribute(data[i].rstrip())
        source = source.replace("KEY"+str(i),filtered)

    return source


if __name__ == '__main__':

    import sys

    source = open(args.r[0]).read()

    for i in range(1,len(args.f)+1):
      filtered = escapeForXmlAttribute(open(args.f[i-1]).read().rstrip())
      source = source.replace("KEY"+str(i),filtered)

    open(args.o[0],'w').write(source)



