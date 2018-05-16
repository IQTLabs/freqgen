from freqgen import k_mer_frequencies

assert k_mer_frequencies("INQTEL", 1) == {'E': 0.16666666666666666,
                                          'I': 0.16666666666666666,
                                          'L': 0.16666666666666666,
                                          'N': 0.16666666666666666,
                                          'Q': 0.16666666666666666,
                                          'T': 0.16666666666666666}

assert k_mer_frequencies("GATGATGGC", 3) == {'ATG': 0.2857142857142857,
                                             'GAT': 0.2857142857142857,
                                             'GGC': 0.14285714285714285,
                                             'TGA': 0.14285714285714285,
                                             'TGG': 0.14285714285714285}

assert k_mer_frequencies("GATGATGGC", 2, include_missing="dna") == {'AA': 0,
                                                                    'AC': 0,
                                                                    'AG': 0,
                                                                    'AT': 0.25,
                                                                    'CA': 0,
                                                                    'CC': 0,
                                                                    'CG': 0,
                                                                    'CT': 0,
                                                                    'GA': 0.25,
                                                                    'GC': 0.125,
                                                                    'GG': 0.125,
                                                                    'GT': 0,
                                                                    'TA': 0,
                                                                    'TC': 0,
                                                                    'TG': 0.25,
                                                                    'TT': 0}
