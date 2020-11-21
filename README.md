# BVEP
Bayesian Virtual Epileptic Patient (BVEP) to invert a individualized whole-brain model of epilepsy spread by PPLs using No-U-Turn Sam-
pler (NUTS) and Automatic Differentiation Variational Inference (ADVI).

Installation: 


For simulation using TVB:


https://www.thevirtualbrain.org/tvb/zwei

For inference using Stan:


https://mc-stan.org/



For inference using PyMC3:

https://docs.pymc.io/


### Ref:
M. Hashemi, A.N. Vattikonda, V. Sip, M. Guye, F. Bartolomei, M.M. Woodman, V.K. Jirsa,
The Bayesian Virtual Epileptic Patient: A probabilistic framework designed to infer the spatial map of epileptogenicity in a personalized large-scale brain model of epilepsy spread,
NeuroImage,
Volume 217,
2020,
116839,
ISSN 1053-8119,
https://doi.org/10.1016/j.neuroimage.2020.116839.
(http://www.sciencedirect.com/science/article/pii/S1053811920303268)


Abstract: Despite the importance and frequent use of Bayesian frameworks in brain network modeling for parameter inference and model prediction, the advanced sampling algorithms implemented in probabilistic programming languages to overcome the inference difficulties have received relatively little attention in this context. In this technical note, we propose a probabilistic framework, namely the Bayesian Virtual Epileptic Patient (BVEP), which relies on the fusion of structural data of individuals to infer the spatial map of epileptogenicity in a personalized large-scale brain model of epilepsy spread. To invert the individualized whole-brain model employed in this study, we use the recently developed algorithms known as No-U-Turn Sampler (NUTS) as well as Automatic Differentiation Variational Inference (ADVI). Our results indicate that NUTS and ADVI accurately estimate the degree of epileptogenicity of brain regions, therefore, the hypothetical brain areas responsible for the seizure initiation and propagation, while the convergence diagnostics and posterior behavior analysis validate the reliability of the estimations. Moreover, we illustrate the efficiency of the transformed non-centered parameters in comparison to centered form of parameterization. The Bayesian framework used in this work proposes an appropriate patient-specific strategy for estimating the epileptogenicity of the brain regions to improve outcome after epilepsy surgery.

Keywords: Bayesian inference; Personalized brain network model; Epileptic seizures; Epileptogenicity


### Fundings: 
the French National Research Agency (ANR) as part of the second “Investissements d’Avenir” program (ANR-17-RHUS-0004, EPINOV), European Union's Horizon 2020 research and innovation programme under grant agreement No. 785907 (SGA2), and No. 945539 (SGA3) Human Brain Project, and the SATT Sud-Est (827-SA-16-UAM).
