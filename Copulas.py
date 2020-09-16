# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
from mixedvines.mixedvine import MixedVine, Marginal
from scipy.optimize import minimize
from mixedvines.copula import Copula, IndependenceCopula
from scipy.stats import rv_continuous, norm, gamma, poisson, binom, nbinom

class MixedVine_(MixedVine):
    class VineLayer_(MixedVine.VineLayer): 
        def fit(self, samples, is_continuous, is_adjusted=None, trunc_level=None):
            '''
            Fits the vine tree to the given samples.  This method is supposed
            to be called on the output layer and will recurse to its input
            layers.

            Parameters
            ----------
            samples : array_like
                n-by-d matrix of samples where n is the number of samples and d
                is the number of marginals.
            is_continuous : array_like
                List of boolean values of length d, where d is the number of
                marginals and element i is `True` if marginal i is continuous.
            manual_adjust : array_like
                list of distributions of length d, where d is the number of 
                marginals and element i is a scipy distribution if it was
                adjusted beforehand or None if the program has to ajust it.
            trunc_level : integer, optional
                Layer level to truncate the vine at.  Copulas in layers beyond
                are just independence copulas.  If the level is `None`, then
                the vine is not truncated.  (Default: `None`)

            Returns
            -------
            output_urvs : array_like
                The output uniform random variates of the layer.  Can be
                ignored if this is the output layer.
            '''
            if self.is_marginal_layer():
                output_urvs = np.zeros(samples.shape)
                if is_adjusted:
                    for i in range(samples.shape[1]):
                        self.marginals[i] = Marginal_.fit(samples[:, i],
                                                     is_continuous[i],
                                                     is_adjusted[i])
                        output_urvs[:, i] = self.marginals[i].cdf(samples[:, i])
                else:
                    for i in range(samples.shape[1]):
                        self.marginals[i] = Marginal_.fit(samples[:, i],
                                                     is_continuous[i])
                        output_urvs[:, i] = self.marginals[i].cdf(samples[:, i])
            else:
                input_urvs = self.input_layer.fit(samples, is_continuous,is_adjusted)
                truncate = trunc_level and samples.shape[1] \
                    - len(self.input_indices) > trunc_level - 1
                output_urvs = np.zeros((samples.shape[0],
                                        len(self.input_indices)))
                for i, i_ind in enumerate(self.input_indices):
                    if truncate:
                        self.copulas[i] = IndependenceCopula()
                    else:
                        self.copulas[i] = Copula.fit(input_urvs[:, i_ind])
                    output_urvs[:, i] \
                        = self.copulas[i].ccdf(input_urvs[:, i_ind])
            return output_urvs
    
    @staticmethod
    def _construct_c_vine(element_order):
        '''
        Constructs a c-vine tree without setting marginals or copulas.  The
        c-vine tree is constructed according to the input element order.  The
        index of the element with the most important dependencies should come
        first in the input argument.

        Parameters
        ----------
        element_order : array_like
            Permutation of all element indices.

        Returns
        -------
        root : VineLayer
            The root layer of the canonical vine tree.
        '''
        dim = len(element_order)
        marginals = np.empty(dim, dtype=Marginal)
        layer = MixedVine_.VineLayer_(marginals=marginals)
        identity_order = np.arange(dim - 1)
        for i in range(1, dim):
            input_indices = []
            if i == 1:
                order = element_order
            else:
                order = identity_order
            # For each successor layer, generate c-vine input indices
            for j in range(dim - i):
                input_indices.append(np.array([order[0], order[j+1]]))
            copulas = np.empty(len(input_indices), dtype=Copula)
            # Generate vine layer
            layer = MixedVine_.VineLayer_(input_layer=layer,
                                        input_indices=input_indices,
                                        copulas=copulas)
        root = layer
        return root
    
    @staticmethod
    def fit(samples, is_continuous, manual_adjust=None, trunc_level=None,
            do_refine=False, keep_order=True):
        '''
        Fits the mixed vine to the given samples.

        Parameters
        ----------
        samples : array_like
            n-by-d matrix of samples where n is the number of samples and d is
            the number of marginals.

        is_continuous : array_like
            List of boolean values of length d, where d is the number of
            marginals and element i is `True` if marginal i is continuous.

        trunc_level : integer, optional
            Layer level to truncate the vine at.  Copulas in layers beyond are
            just independence copulas.  If the level is `None`, then the vine
            is not truncated.  (Default: `None`)

        do_refine : boolean, optional
            If `True`, then all pair copula parameters are optimized jointly at
            the end.  (Default: `False`)

        keep_order : boolean, optional
            If `False`, then a heuristic is used to select the vine structure.
            (Default: `False`)

        Returns
        -------
        vine : MixedVine
            The mixed vine with parameters fitted to `samples`.
        '''

        dim = samples.shape[1]
        vine = MixedVine_(dim)
        if not keep_order:
            element_order = MixedVine_._heuristic_element_order(samples)
            vine.root = MixedVine_._construct_c_vine(element_order)
        vine.root.fit(samples, is_continuous, manual_adjust, trunc_level)
        if do_refine:
            # Refine copula parameters
            initial_point = vine.root.get_all_params()
            bnds = vine.root.get_all_bounds()

            def cost(params):
                '''
                Calculates the cost of a given set of copula parameters.
                '''
                vine.root.set_all_params(params.tolist())
                vals = vine.logpdf(samples)
                return -np.sum(vals)

            result = minimize(cost, initial_point, method='SLSQP',
                              bounds=bnds)
            vine.root.set_all_params(result.x.tolist())
        return vine
    
class Marginal_(Marginal):
    @staticmethod
    def fit(samples, is_continuous, is_final=False):
        '''
        Fits a distribution to the given samples.
        Parameters
        ----------
        samples : array_like
            Array of samples.
        is_continuous : bool
            If `True` then a continuous distribution is fitted.  Otherwise, a
            discrete distribution is fitted.
        is_final: A list with elements [True, distribution1], or [False]

        Returns
        -------
        best_marginal : Marginal
            The distribution fitted to `samples`.
        '''
        # Mean and variance
        if not is_final:
            #print("THIS IS NOT FINAL")
            mean = np.mean(samples)
            var = np.var(samples)
            # Set suitable distributions
            if is_continuous:
                if np.any(samples <= 0):
                    options = [norm]
                else:
                    options = [norm, gamma]
            else:
                if var > mean:
                    options = [poisson, binom, nbinom]
                else:
                    options = [poisson, binom]
            params = np.empty(len(options), dtype=object)
            marginals = np.empty(len(options), dtype=object)
            # Fit parameters and construct marginals
            for i, dist in enumerate(options):
                if dist == poisson:
                    params[i] = [mean]
                elif dist == binom:
                    param_n = np.max(samples)
                    param_p = np.sum(samples) / (param_n * len(samples))
                    params[i] = [param_n, param_p]
                elif dist == nbinom:
                    param_n = mean * mean / (var - mean)
                    param_p = mean / var
                    params[i] = [param_n, param_p]
                else:
                    params[i] = dist.fit(samples)
                rv_mixed = dist(*params[i])
                marginals[i] = Marginal(rv_mixed)
            # Calculate Akaike information criterion
            aic = np.zeros(len(options))
            for i, marginal in enumerate(marginals):
                aic[i] = 2 * len(params[i]) \
                     - 2 * np.sum(marginal.logpdf(samples))
            best_marginal = marginals[np.argmin(aic)]
        else:
            best_marginal = is_final
        return best_marginal


#### Test
#########
#########
values = [10, 20, 30]
probabilities = [0.2, 0.5, 0.3]
b=rv_discrete(values=(values, probabilities)).freeze()
dim = 3  # Dimension
vine = MixedVine(dim)
vine.set_marginal(0, exponweib(a=1, c=2.09, scale=10.895, loc=0)) #norm(0,1)
vine.set_marginal(1, poisson(5))
vine.set_marginal(2, b)
vine.set_copula(1, 0, GaussianCopula(0.5))
vine.set_copula(1, 1, FrankCopula(4))
vine.set_copula(2, 0, ClaytonCopula(5))
size=1000
samples = vine.rvs(size)
samples

order=vine._heuristic_element_order(samples)
samples[:,np.arange(samples.shape[1])]=samples[:,order]

is_continuous = [False, False, True]
is_final = [b,poisson(5),exponweib(a=1, c=2.09, scale=10.895, loc=0)]
vine2 = MixedVine_.fit(samples, is_continuous,is_final)

hj=r.shape[0]
count = [np.sum(r[:,11] == i)/hj for i in list(np.unique(r[:,11]))]