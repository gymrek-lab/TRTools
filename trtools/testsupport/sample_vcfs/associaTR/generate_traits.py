#!/usr/bin/env python3

import pathlib
import numpy as np
import numpy.random

seed = 2

SCRIPT_DIR = pathlib.Path(__file__).parent.resolve()

with open(SCRIPT_DIR / 'samples.txt') as samples_file:
    samples = np.array([int(sample.strip()) for sample in samples_file.readlines() if 'IID' not in sample])
n_samples = len(samples)

n_traits_1 = 10
n_traits_2 = 5

traits = []
rng = np.random.default_rng(seed=seed)
for count, n_traits in enumerate((n_traits_1, n_traits_2)):
    traits.append(rng.random(size=(n_samples, n_traits)))
    np.save("traits_{}.npy".format(count), traits[-1])

all_traits = np.hstack(traits)
for trait_array, name in ((traits[0], 'single'), (all_traits, 'combined')):
    with open(name + '_traits_for_plink.tab', 'w') as out:
        out.write('IID\t' + '\t'.join('trait_' + str(num) for num in range(trait_array.shape[1])) + '\n')
        for row, sample in enumerate(samples):
            out.write(str(sample) + '\t' + '\t'.join('{:0.9}'.format(val) for val in trait_array[row, :]) + '\n')

# for sample merge tests
samples_40 = list(range(5, 45))
np.save('traits_0_40_samples.npy', np.hstack((samples.reshape(-1, 1), traits[0]))[50::-1, :][samples_40, :]) # samples 5-46 in reverse order
with open('samples_6_to_45.txt', 'w') as samples_40_file:
    samples_40_file.write('#IID\n')
    for sample in samples[samples_40]:
        samples_40_file.write(str(sample) + '\n')

samples_45 = [*range(21),21,23,25,27,29,*range(31, 50)]
np.save('traits_1_45_samples.npy', np.hstack((samples.reshape(-1, 1), traits[1]))[samples_45, :])  # all samples excluding 23, 25, 27, 29, 31 (base one)
with open('45_samples.txt', 'w') as samples_45_file:
    samples_45_file.write('#IID\n')
    for sample in samples[samples_45]:
        samples_45_file.write(str(sample) + '\n')

with open('35_samples.txt', 'w') as samples_35_file:
    samples_35_file.write('#IID\n')
    for idx, sample in enumerate(samples):
        if idx not in samples_45 or idx not in samples_40:
            continue
        samples_35_file.write(str(sample) + '\n')

# These regenerate test files for associaTR using a seed. There's currently a bug in the associaTR tests where sometimes a rounding error
# in the test comparisons causes the tests to fail, so rerunning this with a different seed my cause the tests to fail even though you in theory should be able to
# regenerate this data randomly and still have the tests pass. Another bug is that occasionally plink will not test a locus because VIF too high even though my code
# will which will cause a comparison to fail.
