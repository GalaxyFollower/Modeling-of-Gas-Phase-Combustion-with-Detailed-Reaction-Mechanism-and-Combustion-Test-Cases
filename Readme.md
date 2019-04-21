# ABSTRACT

A mathematical model for Closed Homogeneous Constant Pressure Reactor (CHCPR) and
one-dimensional Premixed Reactive Flow (PRFM) is developed for the thesis. Detailed
gas phase reaction kinetics is coupled with thermodynamics and flow properties to develop
the reactor models. Global reaction mechanisms are used for simulating combustion in
complex flows. However, study of accurate combustion physics efficiently requires
detailed chemical kinetics in simplified flow system as it gives more accurate results for
studying ignition, explosion, reaction progression and other transient properties of flame.
Detailed reaction mechanism with three Arrhenius rate parameters and third-body
efficiencies are considered. Thermodynamic properties are calculated using the 7
coefficient NASA polynomials. CHCPR calculates the transient reaction properties like
temperature, species concentration, volumetric heat release, etc., during the progression of
reaction at constant pressure. PRFM solves both the reaction properties and flow properties
like pressure and velocity along the length of the combustor. Addition of source terms in
PRFM is considered to facilitate exchange of heat and mass from surrounding. Governing
equations for the reactor models are derived in the form of a system of nonlinear Ordinary
Differential Equations (ODEs), which are solved using ode15s for CHCPR and Backward
Differentiation Newton-Raphson Method with Finite Difference Jacobian for PRFM.
To incorporate complex flow characteristics like laminar and turbulent mixing and
transport properties in the combustion test cases, ANSYS Fluent is used. Preliminary
design of the gas combustor includes, design of geometry for temperature control at exit
for rated power. Calculation of the fuel mass flow rate for rated power in CHCPR and
diluent mass flow rate for temperature control is done in PRFM. Physics behind the simple
premixed flow, bluff body and cavity flame holding are studied in ANSYS Fluent.

# Closed Homogenous Constant Pressure Reactor Model
Closed homogenous constant pressure reactor is like a system of reactants confined in a
piston-cylinder arrangement that react at each and every location within the gas volume
at the same rate. Thus, there are no temperature or composition gradients within the
mixture, and a single temperature and set of species concentrations suffice to describe
the evolution of the system. For exothermic combustion reactions, both the temperature
and volume will increase with time, pressure held constant.

In CHCPR, the system’s mass is constant with no inflow and outflow of species from the
system. The thermodynamic properties are uniform and homogenous only varying with
time. The system is transient zero-dimensional with rate of conversion of reactants to
products is controlled by chemical reaction rates only. The conservation equations; mass
conservation, species conservation, energy conservation and equation of state, are used
to derive the equations for closed homogenous constant pressure reactor with
assumptions:

* Some reaction rate depends on pressure as well as temperature but only temperature dependence is considered in solver.
* No phase change during reaction.
* No surface reactions
* No heat transfers between the surface and surrounding

# Premixed Reactive Flow Model
General combustion devices are long and slender, as in a gas-turbine combustion
chamber, a cement kiln or a gasifier, in which air fuel mixture enter at one end and
products leave at the other. Such devices can be idealized as plug-flow thermochemical
reactors. In this reactor operating in steady state, properties such as velocity, temperature,
pressure, and mass fractions of species vary principally along the length of the reactor 
Therefore, for such reactors, equations of mass, momentum, and energy
can be simplified to their one-dimensional forms.

In the context of reactive flow, mixture mass conservation, species mass conservation,
momentum conservation and energy conservation was taken into account along with
ideal gas equation. Heat flow in and out of the reactor is modelled adding heat flux source
term in energy equation. Injection of any species or combination of species through any
portion of lateral surface of the reactor is incorporated as an additional source term in
RHS of conservation of mixture mass, species mass, momentum and energy.
The assumptions made for the reactor model are:

* Steady state.
* Fuel and air is premixed. The input mixture of fuel and air is fully mixed.
* Uniform properties in the direction of flow i.e. one dimensional flow. That means
the variation of velocity, temperature, pressure etc. in the direction perpendicular
to the flow is assumed negligible.
* Species transport due to thermal and mass diffusivities are neglected due to
domination of advective transport in the only one properties varying dimension.
* Ideal frictionless flow with ideal gas behavior.
* The rate of conversion of reactants to products is controlled by chemical
reaction rates only. Mixing process does not affect the conversion.
* No surface reactions.

