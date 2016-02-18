# Figure out how far from gene every motif is placed
# open file and take distance of motif from exact gene
# then sum it, want to get some kind of plot
# its only for upstream part of gene, for downstream you need to make some changes



def splitout(handle):
    #
    # line = handle.readline().strip()
    # if line.rstrip() != ('Name\tMotif'+'\t'+'Coordinate'+'\t'+'Strand'+'\t'+'Name'+'\t'+'Upstream'+'\t' +
    #                      'Upstream dist'+'\t'+'Downstream'+'\t'+'Downstream dist'):
    #     raise SyntaxError("Third line not recognised as OUT")
    # NAME=0
    # MOTIF=1
    # COORDINATE=2
    # STRAND=3
    # NAME of gene, if motif in = 4
    # UPSTREAM=5
    # UPSTREAMD=6
    # DOWNSTREAM=7
    # DOWNSTREAMD=8
    for line in handle:
        parts = line.rstrip().split('\t')
        yield parts


file = open('C:\моё\универ\прогаем\Bio6A6T\Bio\AnalysisUpstream.txt')
for nParts in splitout(file):
    import matplotlib.pyplot as plt

    rangen = 300  # change here
    Gene = nParts[0]  # and here

    handle = open('C:\моё\универ\прогаем\Bio6A6T\Bio\motif_out.txt')
    result = [0]*rangen
    for parts in (splitout(handle)):
        up = parts[5].split('/')
        up.append('hypothetical protein')
        updist = parts[6].split('/')
        updist.append('100')
        if up[0] == Gene and int(updist[0]) < rangen:  # first - name, second - range of distance
            result[int(updist[0])] += 1
        if up[1] == Gene and int(updist[1]) < rangen:
            result[int(updist[1])] += 1
    handle.close()

    summ = 0
    for i in range(rangen):
        result[i] /= int(nParts[1])  # change here, full numbers of 'transcriptional regulator' (take from analysis)
        result[i] *= 100
        summ += result[i]
    print(str(summ) + '% is in first ' + str(rangen) + ' of ' + Gene)

    if summ > 25 and int(nParts[1]) > 50:
        plt.plot(result)
        # plt.savefig(str(summ) + '% is in first ' + str(rangen) + ' of ' + Gene + '.png', format='png')
        plt.show()
