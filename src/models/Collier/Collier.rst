
``Collier``
===========

A model for lateral inhibition during cell differentiation motivated by the Notch-Delta pathway.
The model implementation follows the tutorial by Formosa-Jordan et al. [2014].

Model Fundamentals
------------------
The (extended) Collier model describes the per cell development of the protein
concentrations Notch (n), Delta (d), and NICD (Notch intra-cellular domain; r for repressor) using 3 coupled differential equations.
Transmembrane interactions of the proteins are included as next neighbor interactions using the mean protein concentrations of all next neighbors.

.. math::
    \dot{n_i} = \mu \left( \beta_n - k_t n_i \langle d_i \rangle - n_i \right) \\
    \dot{d_i} = \nu \left( \frac{\beta_d}{(1 + r_i)^h} - k_t d_i \langle n_i \rangle - d_i \right) \\
    \dot{r_i} = \frac{\beta (k_t d_i \langle d_i \rangle)^m}{1 + (k_t d_i \langle d_i \rangle)^m} - r_i

where the the quantities are unitless.
The time has been normalized to the decay rate of the repressor NICD.
These equations describe the intra-cellular inhibition of Delta by NICD, the trans-membrane cleavage of Notch and endocytosis of Delta causing the upregulation of the repressor, and the cis-inactivation of Notch and Delta.
The quantities in brackets denote the average over next neighbors.
For more details see Formosa-Jordan et al. [2014].


1.Formosa-Jordan, P. & Sprinzak, D. Modeling Notch Signaling: A Practical Tutorial. in Notch Signaling (eds. Bellen, H. J. & Yamamoto, S.) vol. 1187 285–310 (Springer New York, 2014).

2.Sprinzak, D. et al. Cis-interactions between Notch and Delta generate mutually exclusive signalling states. Nature 465, 86–90 (2010).

3.Sprinzak, D., Lakhanpal, A., LeBon, L., Garcia-Ojalvo, J. & Elowitz, M. B. Mutual Inactivation of Notch Receptors and Ligands Facilitates Developmental Patterning. PLoS Comput Biol 7, e1002069 (2011).
