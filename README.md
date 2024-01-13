# N-Body Simulation of Celestial Bodies
## Multibody Jool Simulation
Earliest versions of this project involved simulating Jool and it's Moons from Kerbal Space Program. Each version attempted to improve readability and general form.
### Version Notes
- V1 established general form, calculation of initial conditions, and the core simulation loop
- V2 attempted to render the video alongside the simulation, each time a new timestep was simulated a new frame would be rendered
- V3 was the final version of this Joolian simulation before moving onto entire solar system simulation

[![Test](https://img.youtube.com/vi/R7MJvFjdgdw/0.jpg)](https://www.youtube.com/watch?v=R7MJvFjdgdw)

## Multibody Kerbol Simulation
Built upon Multibody Jool Simulation with emphasis on scalability and fixing form issues
### Version 1 Notes
- Moved Preallocation and Initilizations to beginning of document where possible
- Consistent Index, all arrays are same length with each index referring to same body
- Added Sphere Of Influence (SOI) to body property array
- Changed orbits array to match size of body_prop array
- added reference body to matrices in orbits array to indicate which body the orbital elements are in reference to (i.e. the parent body)
- Made various simulation constants seperate variables for ease of changes
- Cleaned up figure setup
- The trajectory and labels of moons (satelites) will not render until its distance from parent body is greater than SOI (i.e. it has left parent bodies SOI)
- Added progress tracker for longer simulations

### Version 2 Notes
- Changed simulation forloops to eliminate redundant calculations
- Greatly improved runtime

### Version 3 Notes


## Multibody Solar System Simulation
### Version 1 Notes
