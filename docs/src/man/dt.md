# Darboux Transformations

Please see the paper to review the theory of Darboux Trasformations, and the [Darboux Transformation Examples](@ref) page for examples.

## Supported Equations

The most general equation supported by the Darboux transformation in `NonlinearSchrodinger.jl` is of the following form:

```math
i{\psi _x} + S[\psi (x,t)] - i\alpha H[\psi (x,t)] + \gamma P[\psi (x,t)] - i\delta Q[\psi (x,t)] = 0,
```

where
```math
\begin{aligned}
S[\psi (x,t)] &= \frac{1}{2}{\psi _{tt}} + {\left| \psi  \right|^2}\psi, \\
H[\psi (x,t)] &= {\psi _{ttt}} + 6{\left| \psi  \right|^2}{\psi _t}, \\
P[\psi (x,t)] &= {\psi _{tttt}} + 8{\left| \psi  \right|^2}{\psi _{tt}} + 6{\left| \psi  \right|^4}\psi + 4{\left| {{\psi _t}} \right|^2}\psi + 6{\psi _t}^2{\psi ^*} + 2{\psi ^2}\psi _{tt}^*, \\
Q[\psi (x,t)] &= {\psi _{ttttt}} + 10{\left| \psi  \right|^2}{\psi _{ttt}} + 30{\left| \psi  \right|^4}{\psi _t} + 10\psi {\psi _t}\psi _{tt}^* + 10\psi \psi _t^*{\psi _{tt}} + 20{\psi ^*}{\psi _t}{\psi _{tt}} + 10\psi _t^2\psi _t^*.
\end{aligned}
```

Special cases include the cubic nonlinear Schrodinger equation (``\alpha = \gamma = \delta = 0``), the Hirota equation (``\alpha \neq 0, \gamma = \delta = 0``) the Lakshmanan-Porsezian-Daniel (LPD) equation (``\gamma \neq 0, \alpha = \delta = 0``) and the Quintic nonlinear Schrodinger equation (``\delta \neq 0, \alpha = \gamma = 0``).

For this generalized NLSE, we support both the breather and soliton seeds. Additionally, for the cubic NLSE, we support the ``cn`` and ``dn`` seeds. Support for these seeds will be added at some point in the future for the generalized NLSE.