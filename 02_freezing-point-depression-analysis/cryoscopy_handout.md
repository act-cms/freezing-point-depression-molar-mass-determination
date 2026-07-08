# FREEZING POINT DEPRESSION

> [!NOTE]
> Freezing point depression
> 
> Gibbs-Helmholtz equation
> 
> Chemical potential
> 
> $\Delta H_{fus}$, $\Delta S_{fus}$, $\Delta G_{fus}$
> 
> Data Filtering
> 
> Python Functions

In 1868, Alfred Nobel was awarded a US patent for inventing a way to safely transport and use the highly explosive compound nitroglycerin. Nitroglycerin shock-detonates with very little provocation, so his novel idea was to mix nitroglycerin with diatomaceous earth (an absorbent filler) making the explosive portable.[^1] His invention was known as dynamite. However, dynamite was still unstable at cold temperatures ($T_\text{fus}$ is 13.5$\degree$C /56.3$\degree$F for nitroglycerin). When the nitroglycerin crystallized, the friction between crystals was enough to provoke detonation. It was found that mixing another liquid explosive, ethylene glycol dinitrate ($T_\text{fus}$ is -23$\degree$C / -9$\degree$F) with the nitroglycerin produced a more stable product. A 50:50 mixture of nitroglycerine with ethylene glycol dinitrate results in a freezing temperature of -20$\degree$C (-4$\degree$F).[^2] This phenomenon is called *freezing point depression*.

A common, every-day application of this phenomenon is the "antifreeze" fluid for car engines. Car radiators are filled with a mixture of water and ethylene glycol as a coolant. A 50:50 solution of ethylene glycol and water has a freezing temperature of -37$\degree$C (-34.6$\degree$F) allowing for the coolant to circulate even in sub-zero temperatures.

<a id="Fig7-1-Interface"></a>

<p align="center">
<img src="media/Fig7-1.png" width="650">
</p>

**7.1.** The solid-liquid interface at thermal equilibrium

## Thermodynamic Background

The dynamic solid-liquid interface is shown in Figure [7.1](#Fig7-1-Interface). It is dynamic because, at thermal equilibrium and constant pressure, the rate of melting (solid$\rightarrow$liquid) is the same as the rate of freezing (liquid$\rightarrow$solid). The liquid-solid interface is in constant flux with the rates of melting and of freezing in balance. Adding a *solute* to the liquid phase decreases the rate of freezing while not affecting the rate of melting. Thus, with the addition of the solute to the liquid, more of the solid melts until thermal equilibrium may be reached again at a lower temperature.

From a qualitative thermodynamic point of view, we can understand freezing point depression by considering $\Delta G_{fus}$ ("fus" is fusion, solid$\rightarrow$liquid) and that addition of a solute to the liquid increases the entropy of the liquid. 

$$
\begin{align}
\Delta G_{fus} &=\Delta H_{fus} - T\Delta S_{fus} \\
\Delta S_{fus} &= S_{liquid} - S_{solid} \end{align}
$$

 When $S_{liquid}$ increases, so does $\Delta S_{fus}$ ($S_{solid}$ is unchanged by the addition of solute to the liquid). Thus, $\Delta G_{fus}$ becomes more negative, pushing the solid$\rightarrow$liquid equilibrium forward. It is important to note that the decrease in $\Delta G_{fus}$ is not because of a change in $\Delta H_{fus}$ (which would imply a change in the entropy of the surroundings). Rather, $\Delta G_{fus}$ decreases because $\Delta S_{fus}$, that is $\Delta S_\text{SYSTEM}$, increases. Thermodynamic analysis allows for a more quantitative application of the phenomenon of freezing point depression. Return again to the qualitative description of the equilibrium as when the rate of melting (solid$\rightarrow$liquid) is the same as the rate of freezing (liquid$\rightarrow$solid). Thermodynamically, this means that the chemical potential for the substance $A$ is equal in both phases, solid and liquid. 

$$
\mu_A(s) = \mu_A(l)
$$

 The solid will be considered a pure substance while the liquid is treated as an ideal solution with the solvent mole fraction, $\chi_A$.[^3] 

$$
\mu_{A^*(s)}=\mu_{A^*(l)} + RTln\chi_A
$$

 Since the chemical potential, $\mu$, is the molar Gibbs energy, $G_m$, we can rewrite Equation [7.4](#Eq7-4-ChemPotChi) as follows.

$$
\begin{align}
\mu_{A^*(s)} - \mu_{A^*(l)} &= -\Delta G_{fus} = RTln\chi_A \\
-\dfrac{\Delta G_{fus}}{RT} &= ln\chi_A \end{align}
$$

Recognizing and using Taylor series and Maclaurin series is an important take-away from 2nd semester Physical Chemistry. Let $\chi_B$ be the mole fraction of the *solute*, $B$. If $\chi_B \ll$1 , then we may make the following reduction: 

$$
\begin{eqnarray}
ln \chi_A = ln(1-\chi_B) = -\chi_B  \\
-\dfrac{\Delta G_{fus}}{RT} = -\chi_B \end{eqnarray}
$$

Now take the derivative of both sides with respect to T.

$$
\dfrac{1}{R}\dfrac{d}{dT}\bigg(\dfrac{\Delta G_{fus}}{T}\bigg) = \dfrac{d\chi_B}{dT}
$$

The left side of Equation [7.8](#Eq7-8-dGderiv) is recognizable as the Gibbs-Helmholtz equation.[^4]

$$
\begin{align}
-\dfrac{\Delta H_{fus}}{RT^2} &=\dfrac{d\chi_B}{dT} \\
-\dfrac{\Delta H_{fus}}{RT^2}dT&=d\chi_B \end{align}
$$

Now both sides are integrated, with T$^*$ being the melting point of the pure solvent, $A$.

$$
-\dfrac{\Delta H_{fus}}{R} \int_{T^*}^T \dfrac{dT}{T^2} =\int_0^{\chi_B} d\chi_B
$$

 
$$
\dfrac{\Delta H_{fus}}{R}
\bigg(\dfrac{1}{T}-\dfrac{1}{T^*}\bigg) = -\dfrac{\Delta H_{fus}}{R}
\bigg(\dfrac{1}{T^*}-\dfrac{1}{T}\bigg) =\chi_B
$$

Since $T^*$ and $T$ differ by a small amount relative to $T^*$, we can make the following reductions:

$$
\bigg(\dfrac{1}{T^*}-\dfrac{1}{T}\bigg) = \dfrac{1}{TT_{fus}^*}(T-T_{fus}^*)\approxeq\dfrac{\Delta T}{T_{fus}^{*2}}
$$


$$
-\dfrac{\Delta H_{fus}}{RT_{fus}^{*2}}\Delta T = \chi_B
$$

 

$$
\Delta T = -\dfrac{\chi_BRT_{fus}^{*2}}{\Delta H_{fus}} = -K_f\chi_B
$$

 $K_f$ is the freezing-point depression constant for the solvent. Before moving on, it is important to look at a *qualitative* conclusion from Equation [7.14](#Eq7-14-Reductions). $\Delta H_{fus}$, $R$, $T^*$, and $\chi_B$ are all *positive* numbers. Thus, $\Delta T$ must be negative which indicates that the melting/freezing temperature is *lowered* by the addition of a solute: *freezing point depression*. Lowering the melting/freezing temperature effectively pushes the solid$\rightarrow$liquid equilibrium forward. Here we have shown by thermodynamic analysis that which we reasoned above with Equations [7.2](#Eq7-2-DSfus) and [7.3](#Eq7-3-ChemicalPotential). Freezing point depression is a useful tool from determining the molar mass of an unknown solute that is soluble in a known solvent. Equation [7.15](#Eq7-15-DeltaT) must be manipulated to be in terms of the molality of $B$ in the solution, $b_B$. In the limit that $\chi_B\ll$1 and $\chi_A\approxeq$1, $b_B$ and $\chi_B$ are related by Equation [7.16](#Eq7-16-bB).

$$
b_B=\dfrac{n_B}{m_A} =\dfrac{\chi_B n_\text{TOT}}{\chi_A n_\text{TOT} M_A }=\dfrac{\chi_B}{M_A} 
$$

Using the relationship in Equation [7.16](#Eq7-16-bB), we can re-write Equation [7.15](#Eq7-15-DeltaT) in terms of the molality of the solute.

$$
\Delta T = -K_fM_Ab_B = -K_{f}^{'} b_B
$$

where $M_A$ is the molar mass of A. You will calculate $K_{f}^{'}$ for the solvent (cyclohexane) as part of the pre-lab assignment. Be careful with units, remembering that the units of molality are $\dfrac{\text{moles solute}}{\text{kg solvent}}.$

### Discussion Questions

1.  Comment on the qualitative aspects of the cooling curves including

    - Supercooling (if observed, provide a scientific explanation and cite your references)

    - the slope of the cooling curve (give a hypothesis)

2.  Explain the coding modules. Explain how the mathematical techniques used to select data points can justify the choices made. How would the points selected post-supercooling for that best fit line affect $T_{fus}$? Is there a chemical reason to do this? Or is this caused by \"overfitting\", i.e., we are cherry picking data to make the calculations work? How did we use rolling average and Savitzky-Golay filters to smooth the data? Did it make a difference what parameters were set?

3.  Based on the cooling curve, what can you say about the purity of the cyclohexane?

4.  How would impurities in your solvent affect your results?

5.  How does the concentration of the solution influence the accuracy of this method? Can you think of a way to minimize the concentration effects?

6.  If we were to consider the errors in the slopes and intercepts of the best fits in determining $T_\text{fus}$, how would you expect the predicted molar mass to be affected by the $\Delta_{95}$ range (shown as the plt.fill_between as we did in Lab [1](#Lab1-ExpData)) of the intersection? How sensitive is the value of $T_\text{fus}$ to the calculated molar mass?



[^1]: Alfred Nobel, Improved Explosive Compound, US Patent No 78,317, May 26, 1868.

[^2]: Johll, Matthew E. Investigating Chemistry. 3rd ed. (New York, NY: W. H. Freeman and Company, 2013),211.

[^3]: P W. Atkins, Physical Chemistry, 11th ed. (Oxford, United Kingdom: Oxford University Press, 2018), 161.

[^4]: Atkins, Physical Chemistry, 108.

[^5]: Figure 7.3 in David P. Shoemaker, Carl W. Garland, and Joseph W. Nibler, Experiments in Physical Chemistry, 5th ed. (New York: McGraw-Hill, 1989), 202.
