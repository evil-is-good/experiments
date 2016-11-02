import inspect
import keyword

def patche(dimension, file_name):
    """Deleting all excess information from msh file

    :dimension: 2 od 3
    :file_name: file name
    :returns: no

    """
    if dimension == 2:
        type_e = '3 2 0'
    else:
        type_e = '5 2 0'

    fi = open(file_name + '.msh', 'r')
    fo = open(file_name + '_correct' + '.msh', 'w')

    for line in fi:
        if line.find("$Elements") == -1:
            fo.write(line)
        else:
            fo.write(line)
            break

    n_nodes = int(fi.readline())
    s = ''
    for line in fi:
        if line.find(type_e) == -1:
            n_nodes -= 1
        else:
            splt = line.split(' ')
            s += "1 " + ' '.join(splt[1:])
            break
    c = 2
    for line in fi:
        if line.find(type_e) != -1:
            splt = line.split(' ')
            s += str(c) + " " + ' '.join(splt[1:])
            c += 1
        else:
            s += line

    fo.write(str(n_nodes) + '\n')
    fo.write(s)

    fi.close()
    fo.close()
