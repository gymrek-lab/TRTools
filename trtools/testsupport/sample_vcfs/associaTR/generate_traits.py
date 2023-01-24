#!/usr/bin/env python3

import numpy as np
import numpy.random

with open('samples.txt') as samples_file:
    samples = samples_file.readlines()[0].strip().split()
n_samples = len(samples)

n_traits_1 = 10
n_traits_2 = 5

traits = []
for count, n_traits in enumerate((n_traits_1, n_traits_2)):
    traits.append(np.random.rand(n_samples, n_traits))
    np.save(f"traits_{count}.npy", traits[-1])

all_traits = np.hstack(traits)
for trait_array, name in ((traits[0], 'single'), (all_traits, 'combined')):
    with open(f'{name}_traits_for_plink.tab', 'w') as out:
        out.write('IID\t' + '\t'.join(f'trait_{num}' for num in range(trait_array.shape[1])) + '\n')
        for row, sample in enumerate(samples):
            out.write(f'{sample}\t' + '\t'.join(f'{val:0.9}' for val in trait_array[row, :]) + '\n')
            
        

