import matplotlib.pyplot as plt


def splitout(handle):
    line = handle.readline().strip()
    if line.rstrip() != ('Name\tMotif'+'\t'+'Coordinate'+'\t'+'Strand'+'\t'+'Name'+'\t'+'Upstream'+'\t' +
                         'Upstream dist'+'\t'+'Downstream'+'\t'+'Downstream dist'):
        raise SyntaxError("Third line not recognised as OUT")
    # NAME=0
    # MOTIF=1
    # COORDINATE=2
    # STRAND=3
    # NAME=4
    # UPSTREAM=5
    # UPSTREAMD=6
    # DOWNSTREAM=7
    # DOWNSTREAMD=8
    for line in handle:
        parts1 = line.rstrip().split('\t')
        parts = parts1[6].split('/')  # coord is 0 and 1
        parts.append('wtf')
        dist = int(parts[0])
        # yield dist
        if parts[1] != 'wtf':
            yield int(parts[1])


handle = open('/home/eamag/Bio6A6T/Bio/Upstream.txt')
k = 0
mas = []
for i in splitout(handle):
    mas.append(int(i))
    # print (mas[k])
    k += 1
plt.hist(mas, 1000)
plt.show()
handle.close()
