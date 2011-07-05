import numpy as np


def replicates(aCGH):
    probes = aCGH['id']
    ratio = np.log2(aCGH['test_signal']) - np.log2(aCGH['reference_signal'])

    probe_positions = defaultdict(list)
    for i, probe in enumerate(probes):
        probe_positions[probe].append(i)

    for probe in probe_positions.keys():
        if len(probe_positions[probe]) == 1:
            del probe_positions[probe]

    print
    print 'There are %d replicated probes on the array' % len(probe_positions)

    values = list()
    for probe in probe_positions:
        spots = probe_positions[probe]
        values.append(ratio[spots].std()) #deviazione std dei replicati

    return values
    #out_replicated.append(np.array(values))
    #out_names.append(int(sample['Sample name'].split()[1]))


    #plt.figure()
    #plt.boxplot(out_replicated, notch=False, positions=range(len(out_names)))
    #plt.xticks(range(len(out_names)), out_names)
    #plt.show()
