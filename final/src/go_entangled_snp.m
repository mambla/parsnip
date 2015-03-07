function [ ytest ] = go_entangled_snp()

% search for the most entangled SNP, 9 SNPs before and after the missing one
ytest = go_generic(entangled_snp_algorithm(20));

end